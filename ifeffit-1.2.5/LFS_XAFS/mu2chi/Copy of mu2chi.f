c     XAFS Analysis Tool utilizing ifeffit
c     Wrritten by: Grant Schuster
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     NEEDED FOR NON-INTERPOLATION
      integer m, n, p, q, i, j, first_energy,kpoints
      real*8 kmin,kmax,dk,nknots,rbkg
      real*8 e0,kweight
      real*8 pre2, pre1

      real*8 my_energy(100000),my_chik(100000),my_k(100000)
      real*8 my_pre(100000),my_bkg(100000)
      real*8 my_dchik(100000),my_dxmu(100000)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Interpolation stuff


c     Tell the program where to find the ifeffit library file:
     
      include 'ifeffit.inc'


c     Initialise ifeffit and
c     Set up for screen echo so that we can see what is being
c     reported in IFEFFIT:
      i = ifeffit(' ')
      i = ifeffit('reset')
      i = ifeffit('&screen_echo=0')
      i = ifeffit('history(mu2chi_his.txt)')
c     Initialise some parameters

    
      nknots=7
      rbkg=1.1d0
      e0=7.13331d3
      
      pre2=-30
      pre1=-150
      
      kmin=0.0d0
      kmax=1.50d1
      dk=1.d-1
      
      kweight=0.0d0
      
      
      i = iffputsca('i_nknots',nknots)
      i = iffputsca('i_rbkg',rbkg)


      i = iffputsca('pre2',pre2)
      i = iffputsca('pre1',pre1)


      i = iffputsca('i_kmin',kmin)
      i = ifeffit('set(kmin=i_kmin)')
      i = iffputsca('i_kmax',kmax)
      i = ifeffit('set(kmax=i_kmax)')
      i = iffputsca('i_dk',dk)
      i = ifeffit('set(dk=i_dk)')

      i = iffputsca('kweight',kweight)
      i = ifeffit('i_kweight=kweight')



c--------------------------------------------------------------------
c START
c--------------------------------------------------------------------
      i = ifeffit('read_data(data/mo.dat, group=mudata, type=xmu)')
       i = ifeffit('mudata.dxmu=mudata.3')
c      i = ifeffit('mudata.energy=mudata.1')
c      i = ifeffit('mudata.xmu=mudata.2')
c      i = ifeffit('mudata.dxmu=mudata.3')


c--------------------------------------------------------------------
c START PRE_EDGE CALC
c--------------------------------------------------------------------
c     sets mudata.pre to be the pre-edge subtracted xmu
      i = ifeffit('pre_edge(mudata.energy,mudata.xmu,
     &                      pre2=pre2,pre1=pre1,find_e0=true)')
c--------------------------------------------------------------------
c END PRE_EDGE CALC
c--------------------------------------------------------------------

     
     


c--------------------------------------------------------------------
c     Now I should spline the pre_edge_subtracted thingymo.
c     Let's rename the pre-edge subtracted data as gpre (Grant's pre)

      i = ifeffit ('mudata.gpre=mudata.pre')
     
      i = ifeffit('&e0_guess=e0')

c     sets mudata.bkg
      i = ifeffit('spline(mudata.energy, mudata.gpre,
     &     kweight=0, e0=&e0_guess, nknots=i_nknots,
     &        rbkg=i_rbkg)')
c--------------------------------------------------------------------

     
      write(*,*),'hithere50'
      m = iffgetarr('mudata.gpre',my_pre)
      write(*,*),'hithere51'
      p = iffgetarr('mudata.bkg',my_bkg)
      write(*,*),'hithere52'
      n = iffgetarr('mudata.energy',my_energy)
      write(*,*),'hithere53'
      q = iffgetarr('mudata.dxmu',my_dxmu)
      write(*,*),'hithere54'
      j=1
      do while (my_energy(j).LT.e0)
         j=j+1
      end do
      first_energy=j-1
      do while (j.LE.m)
         my_k(j-first_energy)=1.0d-10*SQRT(2*1.6d-19*9.11d-31*
     &        (my_energy(j)-e0))/1.055d-34

         my_chik(j-first_energy)=(my_pre(j)-my_bkg(j))/
     &        my_bkg(j)

         my_dchik(j-first_energy)=my_dxmu(j)/my_bkg(j)
         j=j+1
      end do
      kpoints=m-first_energy
      write(*,*),m,first_energy,kpoints
c     Put both my_chik and my_k back into IFEFFIT:
      i = iffputarr('mudata.gk',   kpoints,my_k)
      i = iffputarr('mudata.gchi', kpoints,my_chik)
      i = iffputarr('mudata.gdchi',kpoints,my_dchik)
c-------------------------------------------------------------------
     

      i = ifeffit('mudata.k    = mudata.gk')
      i = ifeffit('mudata.chi  = mudata.gchi')
      i = ifeffit('mudata.dchi = mudata.gdchi')
      i = ifeffit('mudata.e0 = ones(npts(mudata.k)) * e0')




      call print_echo

      i = ifeffit('write_data(file=data/chi_v_k.dat,
     &                        mudata.k,mudata.chi,mudata.dchi,
     &                        mudata.e0)')

      i = ifeffit('log(close)')
      i = ifeffit('&screen_echo=0')


      i=ifeffit('history(off)')

     
     
     
     
     
     
     
     
     

      end


c ------------------------------------------------------------------
      subroutine gms_interp(kpoints,my_k,my_chik,my_dchik,
     &                    kpoints_int,interp_k,interp_chik,interp_dchik)

      real*8 my_k(1000),my_chik(1000),my_dchik(1000)
      integer j, kpoints

      integer i,k
      integer i_int,kpoints_int,n_fit_pts
      integer ma, i_int_low,i_int_high

      parameter (ma=4,n_fit_pts=4)

      real*8 fitcoeffs(ma)       !cubic fit as ma=4
      real*8 u(n_fit_pts,n_fit_pts),v(n_fit_pts,n_fit_pts)
      real*8 w(n_fit_pts),intchisq
      real*8 cvm(ma,ma)
      real*8 wt,sx,sum_wt
      real*8 int_length
      real*8 a,b,c,d,da,db,dc,dd !coeffs of fit
 
      real*8 interp_k(5000),interp_chi(5000),interp_dchi(5000)
      real*8 interp_chitmp(n_fit_pts-1,5000)
      real*8 interp_dchitmp(n_fit_pts-1,5000)
      real*8 x_int(n_fit_pts),y_int(n_fit_pts),dy_int(n_fit_pts)

      external funcs


      int_length=5.0d-2
      do i=1,5000,1
         interp_k(i)=int_length*i
      end do
c     So interp_k=(0,0.05,0.10,0.15,...)


c     we start with my_chik and my_k
c     The standard deviations are my_dchik
c -----------------------------------------------------------------
c     do the interpolation over 5 points (n_int)
c     The interp data arrays are x_int,y_int and dy_int
c     Let's do the first interpolation for the first five points
c     Original data counter is i_unint
c     Interpolated data counter is i_int
      do i=1,(kpoints-n_fit_pts+1),1
         do j=1,n_fit_pts,1
            x_int(j)=my_k(j+i-1)
            y_int(j)=my_chik(j+i-1)
            dy_int(j)=my_dchik(j+i-1)
         end do
c     We now have the data arrays necessary to do a fit
      
         call svdfit(x_int,y_int,n_fit_pts,dy_int,fitcoeffs,ma,u,v,w,
     &        n_fit_pts,n_fit_pts,intchisq,funcs)

c     Now get the uncertainties in fitcoeffs:

         call svdvar(v,ma,n_fit_pts,w,cvm,ma)

c     the first ma diagonal elements of the cvm matrix
c     correspond to the variances of the fitted parameters

         a=fitcoeffs(1)
         b=fitcoeffs(2)
         c=fitcoeffs(3)
         d=fitcoeffs(4)
         da=sqrt(cvm(1,1))
         db=sqrt(cvm(2,2))
         dc=sqrt(cvm(3,3))
         dd=sqrt(cvm(4,4))

c      write(*,*),a,b,c,d,da,db,dc,dd
c     my_k(i_unint) was the maximum value of k used in the fit
c     The number of points we are going to interpolate up
c     to is i_int

         do j=1,(n_fit_pts-1),1
            i_int_low=aint(my_k(i+j-1)/int_length)+1
            if (i==1) then
               i_int_low=1
            end if
            i_int_high=aint(my_k(i+j)/int_length)
c            write(*,*)'i_high=',i_int_high,' i_low=',i_int_low
            do k=i_int_low,i_int_high,1
               interp_chitmp(j,k)=a*(k*int_length)**3
     &              +b*(k*int_length)**2+c*(k*int_length)+d
               interp_dchitmp(j,k)=da*(k*int_length)**3
     &              +db*(k*int_length)**2+dc*(k*int_length)+dd
            end do
         end do

c        Print them out
c         do i=1,i_int,1
c            write(*,*),interp_k(i),interp_chi(i),interp_dchi(i)
c         end do
      end do

c     Maximum k index
      kpoints_int=aint(my_k(kpoints)/int_length)
      i_int=aint(my_k(1)/int_length)+1
c      write(*,*)'i_int=',i_int

c     Average the data at each point
      do i=1,kpoints_int,1

         interp_chi(i)=0.d0
         interp_dchi(i)=0.d0
         sum_wt=0.d0

         do j=1,n_fit_pts-1,1
            if (interp_dchitmp(j,i).NE.0.d0) then
               wt=1.d0/(interp_dchitmp(j,i)**2)
            else 
               wt=0.d0
            end if
            sum_wt=sum_wt+wt
            interp_chi(i)=interp_chi(i)+wt*interp_chitmp(j,i)
         end do
         if (sum_wt.NE.0.d0) then
            interp_chi(i)=interp_chi(i)/sum_wt
         else
            interp_chi(i)=0.d0
         end if

         do j=1,n_fit_pts-1,1
            if (interp_dchitmp(j,i).NE.0.d0) then
               wt=1.d0/(interp_dchitmp(j,i)**2)
            else
               wt=0.d0
            end if   
            sx=(interp_chitmp(j,i)-interp_chi(i))**2
            interp_dchi(i)=interp_dchi(i)+sx*wt
         end do
         if (sum_wt.NE.0.d0) then
            interp_dchi(i)=sqrt(interp_dchi(i)/(sum_wt*(n_fit_pts-1)))
         else
            interp_dchi(i)=0.d0
         end if
      end do

      return
      end





c ------------------------------------------------------------------
      subroutine plot_mu_v_energy
      integer i
c     Plot the raw mu data
      i = ifeffit('newplot(mudata.energy,mudata.xmu,')
      i = ifeffit('style=points1, color=blue,')
      i = ifeffit('title=\' [\\gm/\\gr] vs E\',xlabel=\'E (eV)\',')
      i = ifeffit('ylabel=\' [\\gm/\\gr] (cm\\u2\\d/g) \',')
      i = ifeffit('xmin=19500,xmax=21000)')
      i = ifeffit('plot(device="/cps",file="mu_v_energy.ps")')
      return
      end


c ------------------------------------------------------------------
      subroutine plot_pre_edgea
      integer i
c     Plot the pre-edge background function:
      i = ifeffit('mudata.preedge=mudata.xmu-mudata.pre')
      i = ifeffit('newplot(mudata.energy,mudata.xmu,')
      i = ifeffit('style=solid, color=blue,')
      i = ifeffit('title=\' [\\gm/\\gr] vs E\',xlabel=\'E (eV)\',')
      i = ifeffit('ylabel=\' [\\gm/\\gr] (cm\\u2\\d/g) \')')
      i = ifeffit('plot(mudata.preedge,style=dashed,')
      i = ifeffit('color=red)')
      i = ifeffit('plot(device="/cps",file="pre_edgea.ps")')
      return
      end

c ------------------------------------------------------------------
      subroutine plot_ifeffit_autobk_error
      integer i
c     Plot the pre-edge function function:
      i = ifeffit('newplot(mudata.energy,mudata.pre,')
      i = ifeffit('style=solid, color=blue,')
      i = ifeffit('title=\' [\\gm/\\gr]\\dedge\\u vs E\',')
      i = ifeffit('xlabel=\'E (eV)\',')
      i = ifeffit('ylabel=\' [\\gm/\\gr]\\dedge\\u(cm\\u2\\d/g) \')')
      i = ifeffit('plot(device="/cps",file="ifeffit_autobk_error.ps")')
      return
      end

c ------------------------------------------------------------------
      subroutine plot_pre_suba
      integer i
c     Plot the pre-edge subtracted from the experimental data:
      i = ifeffit('newplot(mudata.energy,mudata.pre,')
      i = ifeffit('style=solid, color=blue,')
      i = ifeffit('title=\' [\\gm/\\gr]\\dedge\\u vs E\',')
      i = ifeffit('xlabel=\'E (eV)\',')
      i = ifeffit('ylabel=\' [\\gm/\\gr]\\dedge\\u (cm\\u2\\d/g) \')')
      i = ifeffit('plot(device="/cps",file="pre_suba.ps")')
      return
      end

c ------------------------------------------------------------------
      subroutine plot_spline_new
      integer i
      i = ifeffit('newplot(x=mudata.energy,y=mudata.pre)')
      i = ifeffit('plot(mudata.energy, mudata.bkg,')
      i = ifeffit('title=\' [\\gm/\\gr] vs E\',xlabel=\'E (eV)\',')
      i = ifeffit('ylabel=\' [\\gm/\\gr] (cm\\u2\\d/g) \',')
      i = ifeffit('xmin=(&e0_guess-200),')
      i = ifeffit('xmax=(&e0_guess+800),style=dashed)')
c      i = ifeffit('plot(mudata.gpre,style=dotted)')
      i = ifeffit('plot(device="/cps",file="spline_new.ps")')
      return
      end

c ------------------------------------------------------------------
      subroutine plot_my_chi
      integer i
c     Now plot both the chi from the spline and the gchi:
      i = ifeffit('newplot(mudata.k, mudata.chi, style=points20,')
      i = ifeffit('title=\' \\gx(k) vs k\',ylabel=\'\\gx(k)\',')
      i = ifeffit('xlabel=\' k (\\A\\u-1\\d) \')')
      i = ifeffit('plot(mudata.k,mudata.chi,style=solid)')
      i = ifeffit('plot(mudata.k,mudata.chi,style=points1)')
      i = ifeffit('plot(device="/cps",file="my_chi.ps")')
      return
      end



c ------------------------------------------------------------------
      subroutine plot_interp1
      integer i
      i = ifeffit('newplot(mudata.gk,mudata.gchi,style=points5,')
      i = ifeffit('title=\' \\gx(k) vs k\',')
      i = ifeffit('ylabel=\' \\gx(k)\')')
c      i = ifeffit('dy=myint.dchi)')
      i = ifeffit('plot(mudata.k, mudata.chi,style=points20)')
      i = ifeffit('plot(mudata.k, mudata.chi,style=solid,')
      i = ifeffit('color=red)')
c      i = ifeffit('plot(mydata.k,mudata.dchi,style=points1,color=black)')
      i = ifeffit('plot(device="/cps",file="myinterp1.ps")')
      return
      end


c ------------------------------------------------------------------
      subroutine plot_interp2
      integer i
      i=ifeffit('mudata.ddchi=10*mudata.dchi')

      i=ifeffit('newplot(mudata.k,mudata.chi,style=solid,
     & color=yellow,')
      i = ifeffit('title=\' \\gx(k) vs k\',')
      i = ifeffit('ylabel=\'\\gx(k)\')')

      i = ifeffit('plot(mudata.k,mudata.chi, style=points1,color=blue,')
      i=ifeffit('dy=mudata.ddchi)')
      
      i = ifeffit('plot(device="/cps",file="myinterp2.ps")')
      return
      end


c ------------------------------------------------------------------
      subroutine plot_interpe
      integer i
      i = ifeffit('newplot(mudata.k,mudata.dchi,')
      i = ifeffit('title=\' d\\gx(k) vs k\',')
      i = ifeffit('ylabel=\' d\\gx(k) \')')
      i = ifeffit('plot(device="/cps",file="myinterpe.ps")')
      return
      end

c ------------------------------------------------------------------
      subroutine plot_chi_vs_k_expt
      integer i
c     Plot the results of the background removal
      i = ifeffit('newplot(mudata.k, mudata.chi, style=solid,')
      i = ifeffit('title=\' \\gx(k) vs k\',ylabel=\'\\gx(k)\',')
      i = ifeffit('xlabel=\' k (\\A\\u-1\\d) \')')
      i = ifeffit('plot(device="/cps",file="chi_vs_k_expt.ps")')
      return
      end

c -------------------------------------------------------------------
      subroutine print_echo
      integer i,j, num_buffer
      double precision x_num_buffer
      character*128 x_text(64)

      i = iffgetsca('&echo_lines',x_num_buffer)
      num_buffer = min(64,int(x_num_buffer))
      write(*,*)'echo buffer=',num_buffer
      do 1010 j = 1, num_buffer
         i = iffgetecho(x_text(j))
         write(*,*)'echo line ', j ,' = ', x_text(j)(1:i)
 1010 continue

      return
      end

c --------------------------------------------------------------------
c     Based on Numerical Recipes in Fortran:
c     Solution by Use of Singular Value Decomposition
      subroutine svdfit(x,y,ndata,sig,a,ma,u,v,w,mp,np,chisq,funcs)
      integer ma,mp,ndata,np,NMAX,MMAX
      real*8 chisq,a(ma),u(mp,np),v(np,np),w(np)
      real*8 sig(ndata),x(ndata),y(ndata),TOL
      external funcs
      parameter (NMAX=1000,MMAX=50,TOL=1.e-5) !Max expected ndata and ma
      integer i,j
c     This subroutine uses svbksb,svdcmp
      real*8 sum,thresh,tmp,wmax,afunc(MMAX),b(NMAX)
c     Accumulate coefficients of the fitting matrix:
      do 1012 i=1,ndata     
         call funcs(x(i),afunc,ma)
         tmp=1./sig(i)
         do 1011 j=1,ma
            u(i,j)=afunc(j)*tmp
 1011      continue
         b(i)=y(i)*tmp
 1012   continue
      call svdcmp(u,ndata,ma,mp,np,w,v)
      wmax=0.
      do 1013 j=1,ma
         if(w(j).GT.wmax)wmax=w(j)
 1013   continue
      thresh=TOL*wmax
      do 1014 j=1,ma
         if(w(j).LT.thresh)w(j)=0
 1014   continue
      call svbksb(u,w,v,ndata,ma,mp,np,b,a)
      chisq=0.
      do 1016 i=1,ndata
         call funcs(x(i),afunc,ma)
         sum=0.
         do 1015 j=1,ma
            sum=sum+a(j)*afunc(j)
 1015      continue
         chisq=chisq+((y(i)-sum)/sig(i))**2
 1016   continue
      return
      END
c --------------------------------
c     The user defined funcs(x,afunc,ma) that returns the ma basis
c     functions evaluated at x=x in the array afunc.
c --------------------------------
      subroutine funcs(x,afunc,ma)
      integer ma
      real*8 x, afunc(ma)
      integer i

      do 1021 i=1,ma
         afunc(i)=x**(ma-i)
 1021   continue
      return
      END
c -----------------------------------
c     The Numerical Recipes for Fortran subroutine: svbksb
      subroutine svbksb(u,w,v,m,n,mp,np,b,x)
      integer m,mp,n,np,NMAX
      real*8 b(mp),u(mp,np),v(np,np),w(np),x(np)
      parameter (NMAX=500) !maximum anticipated n value
      integer i,j,jj
      real*8 s,tmp(NMAX)
      do 1032 j=1,n !calculate U^T.B
         s=0.
         if (w(j).NE.0.) then !nonzero result only if w_j is nonzero
            do 1031 i=1,m
               s=s+u(i,j)*b(i)
 1031       continue
            s=s/w(j) !this is the divide by w_j
         endif
         tmp(j)=s
 1032 continue
      do 1034 j=1,n !matrix multiply by V to get answer
         s=0.
         do 1033 jj=1,n
            s=s+v(j,jj)*tmp(jj)
 1033    continue
         x(j)=s
 1034 continue
      return
      END

c -----------------------------------
c     The Numerical Recipes for Fortran subroutine: svdcmp
      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=500) ! Maximum anticipated value of n
C     USES pythag
c      Given a matrix a(1:m,1:n),with physical dimensions 
c     mp by np this routine computes its
c      singular value decomposition,A = U · W · V T .
c     The matrix U replaces a on output.The
c     diagonal matrix of singular alues W is output as a
c     vector w(1:n).Thematrix V (not the
c      transpose V T )is output as v(1:n,1:n).
      INTEGER i,its,j,jj,k,l,nm
      REAL*8 anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),pythag
      g=0.0  ! Householder reduction to bidiagonal form.
      scale=0.0
      anorm=0.0
      do 25 i=1,n
         l=i+1
         rv1(i)=scale*g
         g=0.0
         s=0.0
         scale=0.0
         if(i.le.m)then
            do 11 k=i,m
               scale=scale+abs(a(k,i))
 11         continue
            if(scale.ne.0.0)then
               do 12 k=i,m
                  a(k,i)=a(k,i)/scale
                  s=s+a(k,i)*a(k,i)
 12            continue
               f=a(i,i)
               g=-sign(sqrt(s),f)
               h=f*g-s
               a(i,i)=f-g
               do 15 j=l,n
                  s=0.0
                  do 13 k=i,m
                     s=s+a(k,i)*a(k,j)
 13               continue
                  f=s/h
                  do 14 k=i,m
                     a(k,j)=a(k,j)+f*a(k,i)
 14               continue
 15            continue
               do 16 k=i,m
                  a(k,i)=scale*a(k,i)
 16            continue
            endif
         endif
         w(i)=scale *g
         g=0.0
         s=0.0
         scale=0.0
         if((i.le.m).and.(i.ne.n))then
            do 17 k=l,n
               scale=scale+abs(a(i,k))
 17         continue
            if(scale.ne.0.0)then
               do 18 k=l,n
                  a(i,k)=a(i,k)/scale
                  s=s+a(i,k)*a(i,k)
 18            continue
               f=a(i,l)
               g=-sign(sqrt(s),f)
               h=f*g-s
               a(i,l)=f-g
               do 19 k=l,n
                  rv1(k)=a(i,k)/h
 19            continue
               do 23 j=l,m
                  s=0.0
                  do 21 k=l,n
                     s=s+a(j,k)*a(i,k)
 21               continue
                  do 22 k=l,n
                     a(j,k)=a(j,k)+s*rv1(k)
 22               continue
 23            continue
               do 24 k=l,n
                  a(i,k)=scale*a(i,k)
 24            continue
            endif
         endif
         anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
 25   continue
      do 32 i=n,1,-1 ! Accumulation of right-hand transformations.
         if(i.lt.n)then
            if(g.ne.0.0)then
               do 26 j=l,n !Double division to avoid possible underflow.
                  v(j,i)=(a(i,j)/a(i,l))/g
 26            continue
               do 29 j=l,n
                  s=0.0
                  do 27 k=l,n
                     s=s+a(i,k)*v(k,j)
 27               continue
                  do 28 k=l,n
                     v(k,j)=v(k,j)+s*v(k,i)
 28               continue
 29            continue
            endif
            do 31 j=l,n
               v(i,j)=0.0
               v(j,i)=0.0
 31         continue
         endif
         v(i,i)=1.0
         g=rv1(i)
         l=i
 32   continue
      do 39 i=min(m,n),1,-1 !Accumulation of left-hand transformations.
         l=i+1
         g=w(i)
         do 33 j=l,n
            a(i,j)=0.0
 33      continue
         if(g.ne.0.0)then
            g=1.0/g
            do 36 j=l,n
               s=0.0
               do 34 k=l,m
                  s=s+a(k,i)*a(k,j)
 34            continue
               f=(s/a(i,i))*g
               do 35 k=i,m
                  a(k,j)=a(k,j)+f*a(k,i)
 35            continue
 36         continue
            do 37 j=i,m 
               a(j,i)=a(j,i)*g
 37         continue
         else
            do 38 j= i,m
               a(j,i)=0.0
 38         continue
         endif
         a(i,i)=a(i,i)+1.0
 39   continue
      do 49 k=n,1,-1 !Diagonalization of the bidiagonal form:
                     !Loop over singular values,and over allowed 
                     !iterations.
         do 48 its=1,30
            do 41 l=k,1,-1      !Test for splitting.
               nm=l-1           !Note that rv1(1) is always zero.
               if((abs(rv1(l))+anorm).eq.anorm) goto 2
               if((abs(w(nm))+anorm).eq.anorm) goto 1
 41         continue
 1          c=0.0               !Cancellation of rv1(l),if l > 1
            s=1.0
            do 43 i=l,k
               f=s*rv1(i)
               rv1(i)=c*rv1(i)
               if((abs(f)+anorm).eq.anorm) goto 2
               g=w(i)
               h=pythag(f,g)
               w(i)=h
               h=1.0/h
               c = (g*h)
               s=-(f*h)
               do 42 j=1,m
                  y=a(j,nm)
                  z=a(j,i)
                  a(j,nm)=(y*c)+(z*s)
                  a(j,i)=-(y*s)+(z*c)
 42            continue
 43         continue
 2          z=w(k)
            if(l.eq.k)then      !Convergence.
               if(z.lt.0.0)then !Singular alue is made nonnegative.
                  w(k)=-z
                  do 44 j=1,n
                     v(j,k)=-v(j,k)
 44               continue
               endif
               goto 3
            endif
            if(its.eq.30) pause !'no convergence in svdcmp'
            x=w(l)              !Shift from bottom 2-by-2 minor.
            nm=k-1
            y=w(nm)
            g=rv1(nm)
            h=rv1(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
            g=pythag(f,1.0d0)
            f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
            c=1.0               ! Next QR transformation:
            s=1.0
            do 47 j=l,nm
               i=j+1
               g=rv1(i)
               y=w(i)
               h=s*g
               g=c*g
               z=pythag(f,h)
               rv1(j)=z
               c=f/z
               s=h/z
               f= (x*c)+(g*s)
               g=-(x*s)+(g*c)
               h=y*s
               y=y*c
               do 45 jj=1,n
                  x=v(jj,j)
                  z=v(jj,i)
                  v(jj,j)= (x*c)+(z*s)
                  v(jj,i)=-(x*s)+(z*c)
 45            continue
               z=pythag(f,h)
               w(j)=z           ! Rotation can be arbitrary if z =0
               if(z.ne.0.0)then
                  z=1.0/z
                  c=f*z
                  s=h*z
               endif
               f= (c*g)+(s*y)
               x=-(s*g)+(c*y)
               do 46 jj=1,m
                  y=a(jj,j)
                  z=a(jj,i)
                  a(jj,j)= (y*c)+(z*s)
                  a(jj,i)=-(y*s)+(z*c)
 46            continue
 47         continue
            rv1(l)=0.0
            rv1(k)=f
            w(k)=x
 48      continue
 3       continue
 49   continue
      return
      END
c ----------------------------------
c     The Numerical Recipes for Fortran function: pythag
      function pythag(a,b)
      REAL*8 a,b,pythag
      REAL*8 absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb) then
         pythag=absa*sqrt(1.+(absb/absa)**2)
      else
         if(absb.eq.0.)then
            pythag=0.
         else
            pythag=absb*sqrt(1.+(absa/absb)**2)
         endif
      endif
      return
      END
c -----------------------------------
c     The Numerical Recipes for Fortran subroutine: svdvar
      subroutine svdvar(v,ma,np,w,cvm,ncvm)
      integer ma,ncvm,np,MMAX
      real*8 cvm(ncvm,ncvm),v(np,np),w(np)
      parameter (MMAX=20) !maximum number of fit parameters
      integer i,j,k
      real*8 sum,wti(MMAX)
      do 1051 i=1,ma
         wti(i)=0.
         if(w(i).NE.0.) wti(i)=1./(w(i)*w(i))
 1051 continue
      do 1054 i=1,ma !sum contributions to covariance matrix
         do 1053 j=1,i
            sum=0.
            do 1052 k=1,ma
               sum=sum+v(i,k)*v(j,k)*wti(k)
 1052       continue
            cvm(i,j)=sum
            cvm(j,i)=sum
 1053    continue
 1054 continue
      return
      END
