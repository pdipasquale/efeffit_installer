c     XAFS Analysis Tool utilizing ifeffit
c     Wrritten by: Lucas Smale, Grant Schuster
c     Supervisor:  Chris Chantler
c     Started:     17 September 2003
c
c     This program relies on ifeffit version 1.2.3 or later and 
c     PGPLOT being installed.
c     This program sends commands to iffefit, gets the results back and
c     does a fitting procedure to minimise chisq difference in either
c     momentum or real space and with different k-weightings.
c
c     To execute this program compile using:
c     'make' with the Makefile configured according to the
c     ifeffit documentation
c     and, while using csh, add an alias to your .cshrc file to
c     make life easier:
c     alias gms_xafs home/your_user_name/path_to_program/gms_xafs
c     Then go to the directory where your FEFF8 and experimental data
c     is and type 'gms_xafs' to run.  
c ---------------------------------------------------------------
      implicit none

      integer i,j
      integer max_num_paths
      real*8 kmin,kmax,dk
      real*8 kweight
      
      real*8 rmin,rmax,dr
      real*8 e0
            
      real*8 scale
      real*8 amp_guess
      real*8 amp_target
      real*8 amp
      
      integer do_plots
      integer windowing
      integer full_model
      integer real_chi2
      character*128 arg1
      character*128 ifeffit_arg
      character*128 winstr,kwstr,realstr, filestr
      integer n_args

c     Tell the program where to find the ifeffit library file:
      include 'ifeffit.inc'
  
c     Default options
      do_plots = 0.00     ! No Plots
      windowing = 1       ! Do Widowing
      max_num_paths=19    ! 20 Paths
      full_model = 1      ! Use Full Model
      real_chi2=0.00      ! Use stansard analysis
      
c      amp_guess  = 0.85d0
c      amp_target = 0.85d0
c      scale = 100d0
c      amp = abs(amp_guess - amp_target) * scale
      
c     Initialise some parameters

c     e0=2.0008d4
      e0=7.1131d3
      kmin=3.6d0
      kmax=1.25d1
      dk=3.d-1
      
      rmin=0.0d0
      rmax=4.0d0
      dr=3.d-1
      
      kweight=0.0d0


     
      n_args = iargc()
      i = 1
      do while (i.le.n_args)
        call getarg(i,arg1)
        if (arg1.eq.'-kw') then
          i=i+1
          call getarg(i,arg1)
          read(arg1,*) kweight
          write(*,*)'kweight set'
        else if (arg1.eq.'-np') then
          i=i+1
          call getarg(i,arg1)
          read(arg1,*) max_num_paths
          write(*,*)'num paths set'
        else if (arg1.eq.'-plots') then
          do_plots = 1
          write(*,*)'will do plots'
        else if (arg1.eq.'-nw') then
          windowing = 0
          write(*,*)'nowindowing'
        else if (arg1.eq.'-rm') then
          full_model = 0
          write(*,*)'will load reduced model'
        else if (arg1.eq.'-frm') then
          full_model = -1
          write(*,*)'will load further reduced model'
        else if (arg1.eq.'-hc') then
          real_chi2 = 1
          write(*,*)'will use accurate chi2l'
        else if (arg1.eq.'-w') then
          i=i+1
          call getarg(i,arg1)
          read(arg1,*) kmin
          i=i+1
          call getarg(i,arg1)
          read(arg1,*) kmax
          i=i+1
          call getarg(i,arg1)
          read(arg1,*) dk
        else
          write(*,*)'Usage - lfs_xafs [-nw] [-kw n] [-w n n n]'
          return
        end if
        i=i+1
      end do


c          if(windowing.eq.0) then
c            write(winstr,'(A,I1,A,I2,A)')
c     $                   'w',int(kmin),'-',int(kmax),'_'
c          else
c            winstr = ''
c          endif

c          write(kwstr,'(A,I1,A)'),'k',int(kweight),'_'

c          if(real_chi2.eq.1) then
c            realstr = 'chi2'
c          else
c            realstr = ''
c          end if

c          write(filestr,'(A,A,A,A)')'fit_',winstr,kwstr,realstr
c          write(*,*) kwstr
c          write(*,*) filestr

c       return

c      amp_guess  = 0.85
c      amp_target = 0.85
       
c      scale = 100d0
      
c      write(ifeffit_arg,'(A,F20.9)')'amp_guess=',amp_guess
c      i = ifeffit(ifeffit_arg)
c      write(ifeffit_arg,'(A,F20.9)')'amp_target=',amp_target
c      i = ifeffit(ifeffit_arg)
c      write(ifeffit_arg,'(A,F20.9)')'scale=',scale
c      i = ifeffit(ifeffit_arg)
c      write(ifeffit_arg,'(A,F20.9)')'amp =',amp
c      i = ifeffit(ifeffit_arg)
      
c---------------------------------------------------------------------

c     Initialise ifeffit and
c     Set up for screen echo so that we can see what is being
c     reported in IFEFFIT:
      i = ifeffit(' ')
      i = ifeffit('reset')
      i = ifeffit('history(lfs_xafs_his.txt)')

      i = ifeffit('res1 = abs(min(theta,0))^16')
c ----- IFEFFIT MACROS:
c     Let's create a Fourier Filter Macro in IFEFFIT:

      i = ifeffit('macro do_kweight')
c     i = ifeffit('set $1.chik = $1.chi * $1.k^kweight')
      i = ifeffit('set $1.chik_w = $1.chik * $1.win')
      i = ifeffit('end macro')

      i = ifeffit('guess(e0_cor=-1.7354)')
      i = ifeffit('guess(alpha=0.005)')
      i = ifeffit('guess(ss2=0.00267)')
c     i = ifeffit('guess(s02=0.92)')
      
c  ------------ setting up restraint for S02 -----------------

c      i = ifeffit ('set scale = 2000')
c      i = ifeffit ('set s02a = 0.85')
       i = ifeffit('amp =scale*(s02-s02a)')
c ----------------------------------------------------      
      i = ifeffit ('set sa = 0.92')
      
      if (full_model.eq.1) then
        call load_full_model(max_num_paths)
      else if (full_model.eq.0) then
        call load_redu_model(max_num_paths) 
      else
        call load_fredu_model(max_num_paths) 
      end if

      i = ifeffit('read_data(file=data/chi_dat.dat,group=mudata)')

      if (windowing.eq.0) then
        i = ifeffit('kmax = ceil(mudata.k)+0.2')
        i = iffgetsca('kmax',kmax)
        i = ifeffit('kmin = floor(mudata.k)-0.2')
        i = iffgetsca('kmin',kmin)
        dk=0
      end if

c      i = ifeffit("$fit_space='k'")

c---------------------------------------------------------------------
c     set the E0 value
      e0=7.113d3

      write(ifeffit_arg,'(A,F20.9)')'kweight=',kweight
      i = ifeffit(ifeffit_arg)
      write(ifeffit_arg,'(A,F20.9)')'kmin   =',kmin
      i = ifeffit(ifeffit_arg)
      write(ifeffit_arg,'(A,F20.9)')'kmax   =',kmax
      i = ifeffit(ifeffit_arg)
      write(ifeffit_arg,'(A,F20.9)')'dk     =',dk
      i = ifeffit(ifeffit_arg)

      
c     r window setting

      write(ifeffit_arg,'(A,F20.9)')'rmin   =',rmin
      i = ifeffit(ifeffit_arg)
      write(ifeffit_arg,'(A,F20.9)')'rmax   =',rmax
      i = ifeffit(ifeffit_arg)
      write(ifeffit_arg,'(A,F20.9)')'dr     =',dr
      i = ifeffit(ifeffit_arg)
      
      i = ifeffit('show kweight')
      i = ifeffit('show kmin')
      i = ifeffit('show kmax')
      i = ifeffit('show dk')
      i = ifeffit('show rmin')
      i = ifeffit('show rmax')
      i = ifeffit('show dr')
      
c     i = ifeffit('show amp')

      i = ifeffit('do_bkg=true')
      if (real_chi2.ne.1) then
      
c-------------Do The Fit with feffit
        if (max_num_paths.lt.10) then
          write(ifeffit_arg,'(A,I1,A,A)')'feffit(1-',max_num_paths,
     &        ', k=mudata.k, chi=mudata.chi',
     &        ', group=fit)'
        else
          write(ifeffit_arg,'(A,I2,A,A)')'feffit(1-',max_num_paths,
     &        ', k=mudata.k, chi=mudata.chi',
     &        ', group=fit)'
        end if 
        i = ifeffit(ifeffit_arg)

c-------------Calculate the graph
        if (max_num_paths.lt.10) then
          write(ifeffit_arg,'(A,I1,A)')'ff2chi(1-',max_num_paths,
     &                ', group=feffit)'
        else
          write(ifeffit_arg,'(A,I2,A)')'ff2chi(1-',max_num_paths,
     &              ', group=feffit)'
        end if
        i = ifeffit(ifeffit_arg)
      

        if (max_num_paths.lt.10) then
          write(ifeffit_arg,'(A,I1,A,A)')'feffit2(1-',max_num_paths,
     &        ', k=mudata.k, chi=mudata.chi, dchi=mudata.dchi',
     &        ', group=fit,onlychi2=true)'
        else
          write(ifeffit_arg,'(A,I2,A,A)')'feffit2(1-',max_num_paths,
     &        ', k=mudata.k, chi=mudata.chi, dchi=mudata.dchi',
     &        ', group=fit,onlychi2=true)'
        end if
      i = ifeffit('set(tempkw = kweight)')
        i = ifeffit(ifeffit_arg)
      i = ifeffit('set(kweight = tempkw)')

      else
c-------------Do The Fit with feffit2
        if (max_num_paths.lt.10) then
          write(ifeffit_arg,'(A,I1,A,A)')'feffit2(1-',max_num_paths,
     &        ', k=mudata.k, chi=mudata.chi',
     &        ', dchi=mudata.dchi, group=fit)'
        else
          write(ifeffit_arg,'(A,I2,A,A)')'feffit2(1-',max_num_paths,
     &        ', k=mudata.k, chi=mudata.chi',
     &        ', dchi=mudata.dchi, group=fit)'
        end if
        i = ifeffit(ifeffit_arg)

c      endif
        i = ifeffit(ifeffit_arg)

       endif
      
      i = ifeffit('show @paths')      
      i = ifeffit('show @variables')
c     i = ifeffit('correl(@all,@all,print)')
      i = ifeffit('show chi_square, chi_reduced, r_factor')
      i = ifeffit('show epsilon_k, epsilon_r, n_idp, n_varys')
      i = ifeffit('show &fit_iteration')


      if (real_chi2.ne.1) then
        if (do_plots.eq.1) then
c         Plot the result of the paths

c         Let's set the plot xrange:
          i = ifeffit('xmin_k=0.95*(kmin-dk)')
          i = ifeffit('xmax_k=1.05*(kmax+dk)')

          i = ifeffit('do_kweight mudata')
          i = ifeffit('do_kweight fit')

          i = ifeffit('newplot(mudata.k, mudata.chik_w, color=blue,')
          i = ifeffit("title='k\\uk_weight\\d\\gx(k)(windowed) vs k',")
          i = ifeffit("ylabel='k\\uk_weight\\d\\gx(k)\',xmin=xmin_k,")
          i = ifeffit('xmax=xmax_k,')
          i = ifeffit("xlabel=' k (\\A\\u-1\\d) ', style=solid)")
          i = ifeffit('plot(fit.k,fit.chik_w, color=red, style=dashed)')

          i = ifeffit('plot(device="/cps",file="fit_allpaths.ps")')

c       Plot the theoretical fitted data
          i = ifeffit('newplot(mudata.k,mudata.chi,xmin=0,')
          i = ifeffit("title='\\gx(k) vs k',")
          i = ifeffit('xmin=0,xmax=xmax_k+2,ymin=-0.15,ymax=0.15,')
          i = ifeffit("ylabel='\\gx(k)',")
          i = ifeffit("xlabel=' k (\\A\\u-1\\d) ', style=solid)")
          i = ifeffit('plot(feffit.k, feffit.chi, style=dashed)')
          if (windowing.eq.1) then
            i = ifeffit('mudata.kwin2=0.13*mudata.win')
            i = ifeffit('plot(mudata.k, mudata.kwin2,
     &                    style=dotted, color=black)')
          end if
          i = ifeffit('plot(device="/cps",file="fit.ps")')
        end if
      else

        if (do_plots.eq.1) then
c         Plot the result of the paths

c         Let's set the plot xrange:
          i = ifeffit('xmin_k=0.95*(kmin-dk)')
          i = ifeffit('xmax_k=1.05*(kmax+dk)')


          i = ifeffit('mudata.win = fit.win')

          i = ifeffit('do_kweight mudata')
          i = ifeffit('do_kweight fit')

          i = ifeffit('newplot(mudata.k, mudata.chik_w, color=blue,')
          i =ifeffit("title='k\\uk_weight\\d\\gx(k)(windowed) vs k',")
          i = ifeffit("ylabel='k\\uk_weight\\d\\gx(k)',xmin=xmin_k,")
          i = ifeffit('xmax=xmax_k,')
          i = ifeffit("xlabel=' k (\\A\\u-1\\d) ', style=solid)")
          i = ifeffit('plot(fit.k,fit.chik_w, color=red, style=dashed)')


c       Plot the theoretical fitted data
          i = ifeffit('newplot(mudata.k,mudata.chi,xmin=0,')
          i = ifeffit("title='\\gx(k) vs k',")
          i = ifeffit('xmin=0,xmax=xmax_k+2,ymin=-0.15,ymax=0.15,')
          i = ifeffit("ylabel='\\gx(k)',")
          i = ifeffit("xlabel=' k (\\A\\u-1\\d) ', style=solid)")
          i = ifeffit('plot(fit.k, fit.chi, style=dashed)')

          if (windowing.eq.1) then
            i = ifeffit('mudata.kwin2=0.13*mudata.win')
            i = ifeffit('plot(mudata.k, mudata.kwin2,
     &                    style=dotted, color=black)')
          end if

          i = ifeffit('plot(device="/cps",file="fit.ps")')

        end if
      end if

c      if (windowing.eq.1) then
c        i = ifeffit('write_data(file=spec.dat, $title*, mudata.k,
c     &                          mudata.chi, fit.chi, mudata.kwin2)')
c      else
c        i = ifeffit('write_data(file=spec.dat, $title*, mudata.k,
c     &                          mudata.chi, fit.chi)')
c      end if

c     i = ifeffit('e0_final=e0+e0_cor')

      i = ifeffit('show kweight')

      i = ifeffit('log(file=variables.out)')
      i = ifeffit('&screen_echo=2')

      i = ifeffit('show kweight, kmin, kmax, dk')
      i = ifeffit('show rmin, rmax, dr')
      i = ifeffit('show chi_square, chi_reduced')
      i = ifeffit('show epsilon_k')
      i = ifeffit('show mychi_square, mychi_reduced')

      i = ifeffit('show @variables')
      i = ifeffit('log(close)')
      i = ifeffit('&screen_echo=0')


      i=ifeffit('history(off)')
      end


c--------------------------------------------------------
      subroutine load_full_model(num_paths)
      integer num_paths
      integer i,j
      character*128 ifeffit_arg
      
      if (num_paths.gt.99) then
        num_paths = 99
        write(*,*) 'Warning: Only 99 paths loaded'
      endif

c       i = ifeffit('ss2_norm_correction=0.00023')
c       i = ifeffit('guess(s02=0.92)')

      do j=1,num_paths
         if (j.lt.10) then
                   
            write(ifeffit_arg,'(A,I1,A,I1,A)')'path(index=',j,
     &                                   ',feff=feff/feff000',j,'.dat,'
            i = ifeffit(ifeffit_arg)
            write(ifeffit_arg,'(A,I1,A,I1,A)')'label=path ',j,
     &                                     ', e0=e0_cor)'
            i = ifeffit(ifeffit_arg)
            
            
            write(ifeffit_arg,'(A,I1,A)')'path(',j,
     &                         ', sigma2=max(ss2, 0.0001))'
            i = ifeffit(ifeffit_arg)       
            
            write(ifeffit_arg,'(A,I1,A)')'path(',j,
     &                              ', s02=sa)'
     
            i = ifeffit(ifeffit_arg)
              
            write(ifeffit_arg,'(A,I1,A,I1,A)')'path(',j,
     &                                   ', delr=reff*alpha)'
            i = ifeffit(ifeffit_arg)
         else

            write(ifeffit_arg,'(A,I2,A,I2,A)')'path(index=',j,
     &                                   ',feff=feff/feff00',j,'.dat,'
            i = ifeffit(ifeffit_arg)
             write(ifeffit_arg,'(A,I2,A,I2,A)')'label=path ',j,
     &                                     ', e0=e0_cor)'
            i = ifeffit(ifeffit_arg)
            
            write(ifeffit_arg,'(A,I2,A)')'path(',j,
     &                         ', sigma2=max(ss2, 0.0001))'
            i = ifeffit(ifeffit_arg)
            
            write(ifeffit_arg,'(A,I2,A)')'path(',j,
     &                            ', s02=sa)'
            i = ifeffit(ifeffit_arg)
             
            write(ifeffit_arg,'(A,I2,A,I2,A)')'path(',j,
     &                                   ', delr=reff*alpha)'
            i = ifeffit(ifeffit_arg)
         end if
      end do

      return
      end

c--------------------------------------------------------
      subroutine load_redu_model(num_paths)
      integer num_paths
      integer i,j
      character*128 ifeffit_arg
      
      if (num_paths.gt.99) then
        num_paths = 99
        write(*,*) 'Warning: Only 99 paths loaded'
      endif

c      i = ifeffit('s02_norm_correction=0.00023')

      do j=1,num_paths
         if (j.lt.10) then
c             write(ifeffit_arg,'(A,I1,A)')'guess(s02=0.85)'
c             write(ifeffit_arg,'(A,I1,A)')'guess(s02',j,'=0.85)'
c             i = ifeffit(ifeffit_arg)
            write(ifeffit_arg,'(A,I1,A,I1,A)')'path(index=',j,
     &                                   ',feff=feff/feff000',j,'.dat,'
            i = ifeffit(ifeffit_arg)
            write(ifeffit_arg,'(A,I1,A)')'label=path ',j,
     &                                   ', e0=e0_cor)'
            i = ifeffit(ifeffit_arg)
            write(ifeffit_arg,'(A,I1,A)')'path(',j,
     &                         ', sigma2=max(ss2, 0.0001))'
            i = ifeffit(ifeffit_arg)


             write(ifeffit_arg,'(A,I1,A)')'path(',j,
     &                              ', s02=sa)'
             i = ifeffit(ifeffit_arg)

            write(ifeffit_arg,'(A,I1,A,I1,A)')'path(',j,
     &                                   ', delr=alpha*reff)'
            i = ifeffit(ifeffit_arg)
         else
            write(ifeffit_arg,'(A,I2,A,I2,A)')'path(index=',j,
     &                                   ',feff=feff/feff00',j,'.dat,'
            i = ifeffit(ifeffit_arg)
            write(ifeffit_arg,'(A,I2,A)')'label=path ',j,
     &                                   ', e0=e0_cor)'
            i = ifeffit(ifeffit_arg)
            write(ifeffit_arg,'(A,I2,A)')'path(',j,
     &                         ', sigma2=max(ss2, 0.0001))'
            i = ifeffit(ifeffit_arg)
             write(ifeffit_arg,'(A,I2,A)')'path(',j,
     &                              ', s02=sa)'
             i = ifeffit(ifeffit_arg)

            write(ifeffit_arg,'(A,I2,A,I2,A)')'path(',j,
     &                                   ', delr=alpha*reff)'
            i = ifeffit(ifeffit_arg)
         end if
      end do

      return
      end

c--------------------------------------------------------
      subroutine load_fredu_model(num_paths)
      integer num_paths
      integer i,j
      character*128 ifeffit_arg
      
      if (num_paths.gt.99) then
        num_paths = 99
        write(*,*) 'Warning: Only 99 paths loaded'
      endif

c      i = ifeffit('s02_norm_correction=0.00023')

      do j=1,num_paths
         if (j.lt.10) then

            write(ifeffit_arg,'(A,I1,A,I1,A)')'path(index=',j,
     &                                   ',feff=feff/feff000',j,'.dat,'
            i = ifeffit(ifeffit_arg)
            write(ifeffit_arg,'(A,I1,A)')'label=path ',j,
     &                                   ', e0=e0_cor)'
            i = ifeffit(ifeffit_arg)
            write(ifeffit_arg,'(A,I1,A)')'path(',j,
     &                         ', sigma2=max(ss2, 0.0001))'
            i = ifeffit(ifeffit_arg)            
            
             write(ifeffit_arg,'(A,I1,A)')'path(',j,
     &                              ', s02=sa)'
             i = ifeffit(ifeffit_arg)
             
             write(ifeffit_arg,'(A,I1,A,I1,A)')'get_path(',j,
     &                                    ', prefix=path',j,')'
            i = ifeffit(ifeffit_arg)
            write(ifeffit_arg,'(A,I1,A,I1,A)')'path(',j,
     &                                   ', delr=alpha*reff)'
            i = ifeffit(ifeffit_arg)
         else
            write(ifeffit_arg,'(A,I2,A,I2,A)')'path(index=',j,
     &                                   ',feff=feff/feff00',j,'.dat,'
            i = ifeffit(ifeffit_arg)
            write(ifeffit_arg,'(A,I2,A)')'label=path ',j,
     &                                   ', e0=e0_cor)'
            i = ifeffit(ifeffit_arg)
            write(ifeffit_arg,'(A,I2,A)')'path(',j,
     &                         ', sigma2=max(ss2, 0.0001))'
            i = ifeffit(ifeffit_arg)            

             write(ifeffit_arg,'(A,I2,A)')'path(',j,
     &                              ', s02=sa)'
             i = ifeffit(ifeffit_arg)
            write(ifeffit_arg,'(A,I2,A,I2,A)')'get_path(',j,
     &                                    ', prefix=path',j,')'
            i = ifeffit(ifeffit_arg)
            write(ifeffit_arg,'(A,I2,A,I2,A)')'path(',j,
     &                                   ', delr=alpha*reff)'
            i = ifeffit(ifeffit_arg)
         end if
        end do

       return
      end