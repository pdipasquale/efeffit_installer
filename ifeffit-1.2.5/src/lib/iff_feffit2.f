       subroutine iff_feffit2(str)
c
c given a list of paths and data, fit sum of paths to data
c
c//////////////////////////////////////////////////////////////////////
c Copyright (c) 1997--2000 Matthew Newville, The University of Chicago
c Copyright (c) 1992--1996 Matthew Newville, University of Washington
c
c Permission to use and redistribute the source code or binary forms of
c this software and its documentation, with or without modification is
c hereby granted provided that the above notice of copyright, these
c terms of use, and the disclaimer of warranty below appear in the
c source code and documentation, and that none of the names of The
c University of Chicago, The University of Washington, or the authors
c appear in advertising or endorsement of works derived from this
c software without specific prior written permission from all parties.
c
c THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
c EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
c MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
c IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
c CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
c TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
c SOFTWARE OR THE USE OR OTHER DEALINGS IN THIS SOFTWARE.
c//////////////////////////////////////////////////////////////////////
c
       implicit none
       include 'consts.h'
       include 'keywrd.h'
       include 'arrays.h'
       include 'fft.h'
       include 'fefdat.h'
       include 'pthpar.h'
       include 'feffit.h'
       save

       double precision  toldef, epsmin, p01
       parameter (toldef = 1.d-7, epsmin = 1.d-8, p01 = 1.d-2)
       integer   lenwrk,  lenfvc
       parameter(lenwrk = 2*maxpts*mvarys + 20*mvarys + 2*maxpts)
       parameter(lenfvc = 2*maxpts)

       logical  do_pha, do_mag, do_re, isnum, dofit, do_bkg
       character*256 str*(*), pref
       character*256 name1(mdata), nam_chi(mdata), nam_x(mdata)
       character*256 nam_dchi(mdata)
       character*256 re_arr, dre_arr, xk_arr, winarr, winnam, fit_sp
       character*1024 list, cmdstr

       integer  idata, jchi, jdchi, ier, illen, jk,  ilen, jwin, n1, n2
       integer  istrln, k, i, jprint, nerstp, idx, ntmp, iup
       integer  nkpts,  np, u2ipth, lminfo, nfit1, j1, j0, j2
       integer  iwork(mvarys), ibadx(mvarys), ipos, jxk, ntchi, isp
       integer  ipath_tmp(max_pathindex)
       integer  get_array, imac_tmp, j, jqw
       logical  iftmac, isasca, isamac, uniq_dat
       double precision  xkmin, xkmax, dk1, dk2, xkw, getsca,q, chired
       double precision  xrmin, xrmax, toler, xolow, xohigh, tmp
       double precision  work(lenwrk), fvect(lenfvc), ftemp(lenfvc)
       double precision  fjac(lenfvc, mvarys), tmpchi(maxpts)
       double precision  alpha(mvarys, mvarys), varys(mvarys), small
       double precision  drin, qrange, bvalue,xspl(mtknot), qx
       double precision  xkdat(maxpts), xcdat(maxpts)
       double precision  xddat(maxpts), xwdat(maxpts)
       double precision  tmparr2(maxpts)
       double precision  wtmp, stmp, sumsqr, xidat, eps_k, eps_r
       double precision  qw_tmp(mqwfs), sum2chi(mdata)
       character*128     restr_tmp(max_restraint)
       integer  nrest_tmp
       integer  iff_eval, iff_eval_dp, iff_eval_in, nqw
       logical  onlychi2

       external iff_eval, iff_eval_dp, iff_eval_in, bvalue
       external istrln, u2ipth, isnum,  fitfun, fitfun2, getsca, sumsqr
       external isasca, isamac, get_array

       toler = zero
c
       fit_macro = undef
       fit_m_arg = ' '
       ifit_mac  = -1
       nerstp= 1
       j0    = 0
       jprint= 0
       xrmin = getsca('rmin')
       xrmax = getsca('rmax')
       xkmin = getsca('kmin')
       xkmax = getsca('kmax')
       xkw   = getsca('kweight')
       dk1   = getsca('dk')
       dk2   = getsca('dk')
       xidat = getsca('data_set')
       idata = int(xidat)
       if (idata.le.0) idata = 1
       nfdats= int(xidat)
       eps_k = zero
       eps_r = zero
       xk_arr= undef
       dre_arr = undef
       winarr= undef
       pref  = undef
       call gettxt('kwindow', winnam)
       call gettxt('fit_space',fit_sp)


c  interpret any and all keyword/value pairs for setting options
       call bkeys(str, mkeys, keys, values, nkeys)
       np     =  0
       illen  =  1
       list   =  ' '
       do_re  = .false.
       do_mag = .false.
       do_pha = .false.
       do_bkg = .false.
       onlychi2 = .false.

       lfs_kweight = getsca('kweight')

       nrest_tmp = 0
       
       do 100 i = 1, nkeys
          k = istrln( keys(i))
          if ((keys(i).eq.'prefix').or.(keys(i).eq.'group')) then 
             pref  = values(i)
             call smcase(pref,'a')
          elseif (keys(i).eq.'macro') then
             fit_macro = values(i)
          elseif ((keys(i).eq.'macro_arg')) then
             fit_m_arg = values(i)
          elseif (keys(i)(1:9).eq.'restraint') then
             nrest_tmp = nrest_tmp+1
             if (idata.le.0)  idata = 0
             if (nrest_tmp.gt.max_restraint) then
                call echo(" *** feffit2:: too many restraints")
                call echo("     restraint will be ignored:")
                call echo("      " // values(i))
                values(i) = ''
             endif
             if((values(i).ne.'').and.(values(i).ne.undef)) then
                restr_tmp(nrest_tmp) = values(i)
cc                print*, ' Feffit restraint : ', idata
cc                print*,  nrest_tmp, values(i)(:60)
             endif
          elseif ((keys(i).eq.'kmax')) then
             ier = iff_eval_dp(values(i), xkmax)
          elseif ((keys(i).eq.'kmin')) then
             ier = iff_eval_dp(values(i), xkmin)
          elseif ((keys(i).eq.'error_step')) then
             ier = iff_eval_in(values(i), nerstp)
          elseif ((keys(i).eq.'error_print')) then
             ier = iff_eval_in(values(i), jprint)
          elseif (keys(i).eq.'kwindow') then
             winnam = values(i)
          elseif (keys(i).eq.'altwindow') then
             winarr = values(i)
             call lower(winarr)
          elseif (keys(i).eq.'kweight') then
             ier = iff_eval_dp(values(i), lfs_kweight)
          elseif (keys(i).eq.'dk1') then
             ier = iff_eval_dp(values(i), dk1)
          elseif (keys(i).eq.'dk2') then
             ier = iff_eval_dp(values(i), dk2)
          elseif (keys(i).eq.'dk') then
             ier = iff_eval_dp(values(i), dk1)
             dk2 = dk1
          elseif (keys(i).eq.'chi') then
             re_arr = values(i)
             call lower(re_arr)
          elseif (keys(i).eq.'dchi') then
             dre_arr = values(i)
             call lower(dre_arr)
          elseif (keys(i).eq.'k') then
             xk_arr = values(i)
             call lower(xk_arr)
          elseif ((keys(i).eq.'epsilon_k')) then
             ier = iff_eval_dp(values(i), eps_k)
          elseif ((keys(i).eq.'epsilon_r')) then
             ier = iff_eval_dp(values(i), eps_r)
          elseif ((keys(i).eq.'fit_space')) then
             fit_sp = values(i)
          elseif ((keys(i).eq.'do_bkg')) then
             call str2lg(values(i), do_bkg, ier)
          elseif ((keys(i).eq.'toler')) then
             ier = iff_eval_dp(values(i), toler)
          elseif ((keys(i).eq.'onlychi2')) then
             onlychi2 = .true.
c
c path list:
          elseif (values(i).eq.undef) then
             call str2il(keys(i), max_pathindex, np, ipath_tmp,ier)
             if (ier .eq. 0) then
                jk    = istrln(keys(i))
                list  = list(1:illen)//keys(i)(1:jk)//','
                illen = illen+jk+1
             else
                call echo(' *** feffit2: error generating path list')
                call echo( keys(i)(1:k))
             end if
          else
             call echo(' *** feffit2: unknown key: '//keys(i)(1:k))
          end if
 100   continue 

c===================================================================
c  done reading input params
       idata  = max(1, min(mdata, idata))
       nfdats = max(1, min(mdata, nfdats))
       dofit  = (idata.ge.nfdats)
       bkgfit(idata) = do_bkg
cc       print*, ' feffit: ', idata, bkgfit(idata)

cc       print*, ' __ feffit __ : idata, nfdats = ', idata, nfdats
c===================================================================
c for this data set:
c    resolve chi(k) name, get index of chi data
c    set fft / fit params 
c    get eps_k and eps_r from iff_chieps (unless given explicitly!)
c    set path list

       name1(idata) = pref
       if (name1(idata) .eq. undef) name1(idata) = 'feffit'
       

       jchi  = get_array(re_arr, name1(idata), 0, xcdat)
       if (jchi.le.0) then
          call echo( ' feffit2: no chi(k) data array?')
          return
       end if
       nam_chi(idata) = re_arr
       call prenam(name1(idata),nam_chi(idata))

       jxk  = get_array(xk_arr, name1(idata), 0, xkdat)
       if (jxk.le.0) then
          call echo( ' feffit2: no k-data array?')
          return
       endif
       nam_x(idata) = xk_arr
       call prenam(name1(idata),nam_x(idata))

       jdchi  = get_array(dre_arr, name1(idata), 0, xddat)
       
       if (jxk.le.0) then
          call echo( ' feffit2: no dchi-data array?')
          return
       endif
       nam_dchi(idata) = dre_arr
       call prenam(name1(idata),nam_dchi(idata))


       ntchi = min(jxk,min(jchi,jdchi))

       lfs_nqfit = ntchi

c save this chi(k) data to the appropriate array for fitfun,
c accumulate sum-of-squares for later uniqueness check
       sum2chi(idata) = zero
       do 170 i = 1, ntchi
          lfs_chiq(i) = xcdat(i)
          lfs_dchiq(i) = xddat(i)
          lfs_kq(i) = xkdat(i)
          sum2chi(idata) = sum2chi(idata) + lfs_chiq(i)*lfs_chiq(i)
 170   continue 
       do 175 i = ntchi+1,maxpts
          lfs_chiq(i) = zero
          lfs_dchiq(i) = zero
          lfs_kq(i) = zero
 175   continue 
 

c support for refining bkg(k) with EXAFS parameters
       if (bkgfit(idata)) then
          rbkg_f(idata) = xrmin
          xrmin       = zero
       endif
c
c  check if fit_macro call will work
       if (fit_macro.ne.undef) then
          ifit_mac = 0
          iftmac = isamac(fit_macro,ifit_mac)
          if (.not.iftmac) then
             cmdstr = '  fitting macro not found: '//fit_macro
             call echo (cmdstr)
             fit_macro = undef
             ifit_mac  = 0
          endif
          if (ifit_mac.gt.0) then
             nmac_stop = nmacro
             imac_tmp = imac
          endif
       endif
c
c  set up fft/fit parameters, initialize window functions
       small      = rgrid * p01
       n1         = int( (xrmin + small) / rgrid )  + 1
       n2         = int( (xrmax + small) / rgrid )
       if (n2.le.n1) n2 = n1 + 1
       nrpts(idata)  = n2 - n1 + 1
       rmax(idata)   = rgrid * n2
       rmin(idata)   = rgrid * n1
       dq1(idata)    = dk1
       dq2(idata)    = dk2
       small         = qgrid * p01
       n1            = int( (xkmin + small) / qgrid )  + 1
       n2            = int( (xkmax + small) / qgrid )
       if (n2.le.n1) n2 = n1 + 1
       nqpts(idata)  = n2 - n1 + 1
       qmax(idata)   = qgrid * n2
       qmin(idata)   = qgrid * n1
       nqfit(idata)  = int((qmax(idata) + 5*dq2(idata))/qgrid)
       nqfit(idata)  = int(lfs_kq(lfs_nqfit)/qgrid)



       xinfo(idata)  = 2 * ((rmax(idata)-rmin(idata))*
     $      (qmax(idata)-qmin(idata))) /pi
       if (bkgfit(idata)) then
          drin        = pi / ( qmax(idata) - qmin(idata) )
          nbkg(idata) = 1 + 2 * int( rbkg_f(idata) / drin )
          nbkg(idata) = min(mtknot-korder-1,max(korder+1,nbkg(idata)))
c-bkg where to put knots
          do 220 i = 1, korder
             qknot(i,idata)      =  qmin(idata) - (korder-i) * qgrid
             qknot(nbkg(idata)+i,idata)= qmax(idata) + (i-1) * qgrid
 220      continue
          qrange = qmax(idata) - qmin(idata)
          ntmp   = nbkg(idata) - korder + 1
          do 240 i = korder+1, nbkg(idata)
             qknot(i, idata) = qmin(idata) + (i-korder)*qrange/ntmp
 240      continue

c-bkg add spline coefs to the variable list, initialize to zero
          do 260  i = 1, nbkg(idata)
             cmdstr = ' '
             write(cmdstr, 265)  'bkg',idata,'_',i,'  = 0.0'
             call iff_set('guess',cmdstr,1)
 260      continue
 265      format (a3,i2.2,a1,i2.2,a7)
          nvarys = nvarys + nbkg(idata)
       else
          do 280  i = 1, nbkg(idata)
 275         format (a3,i2.2,a1,i2.2)
             write(cmdstr, 275)  'bkg',idata,'_',i
             if (isasca(cmdstr)) then
                qx = getsca(cmdstr)
 285            format (a3,i2.2,a1,i2.2,a3,g15.9)
                write(cmdstr, 285)  'bkg',idata,'_',i,' = ',qx
                call iff_set('set', cmdstr,0)
             endif
 280      continue
       end if
c
c  kit in k-space? 
       ifft(idata)   = 0
cc       print*, ' Feffit ', fit_sp, ifft(idata), idata
c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
c
      write(*,*)'-----XKWCHECK0.6=',xkw,'-----'
c set k-weight, including multiple k-weights
      write(*,*)'-----NQWCHECK0=',nqw,'-----'
       if (nqw.le.0) nqw = 1
       nqwfs(idata) = nqw
      write(*,*)'-----XKWCHECK0.33=',xkw,'-----'
       do 295 i = 1, mqwfs
c          qwfs(i,idata) = qw_tmp(i)
          qwfs(i,idata) = lfs_kweight !qw_tmp(i) !<--CHANGED THIS!

 295   continue 
       qweigh(idata) = qwfs(1,idata)
      write(*,*)'-----XKWCHECK0.66=',xkw,'-----'
c loop over k-weights to set weight of fit
      write(*,*)'-----NQWCHECK1=',nqw,'-----'
      write(*,*)'-----IDATACHECK1=',idata,'-----'
       do 450 jqw = nqw, 1, -1
          xkw = qwfs(jqw,idata)
      write(*,*)'-----XKWCHECK0.7=',xkw,'-----'


c    determine eps_k/eps_r for this k-weight
c    call 'chieps ('chi_noise' command) to get eps_k / eps_r
          if ((eps_k .le. epsmin).and.(eps_r .le. epsmin)) then 
             ilen = istrln(nam_chi(idata))
 405         format('  chi=',a,',kmin=',f12.3,',kmax=',f12.3,
     $            ',kweight=',f12.3,',dk1=',f12.3,',dk2=',f12.3)
             write(cmdstr, 405) nam_chi(idata)(1:ilen),
     $            xkmin, xkmax,xkw,dk1,dk2
             if (nam_x(idata) .ne. undef) then
                j0 = istrln(cmdstr)
                j2 = istrln(nam_x(idata))
                cmdstr = cmdstr(:j0)//', k='//nam_x(idata)(:j2)
             endif
             call iff_chieps(cmdstr)
             sigdtk  = getsca('epsilon_k')
             sigdtr  = getsca('epsilon_r')
          elseif (eps_k .gt. epsmin) then 
             sigdtk  = eps_k
             wtmp    = 2 * xkw  + one
             sigdtr  = eps_k /  sqrt( 2 * pi * wtmp /
     $            (qgrid * (xkmax**wtmp - xkmin**wtmp )))
          elseif (eps_r .gt. epsmin) then 
             sigdtr  = eps_r
             wtmp    = 2 * xkw  + one
             sigdtk  = eps_r *  sqrt( 2 * pi * wtmp /
     $            (qgrid * (xkmax**wtmp - xkmin**wtmp )))
          end if
          call setsca('epsilon_k', sigdtk)
          call setsca('epsilon_r', sigdtr)

c set number of fit points, data-set weight, and
c a few other misc stuff based on which fit space is used
          rweigh(idata) = one

c k-space fit
          xolow      = qmin(idata)
          xohigh     = qmax(idata)
          nfit(idata)   = 2 * max (1, nqpts(idata))
          weight(jqw,idata) = sigdtk * sigdtk 


          weight(jqw,idata) =  sqrt(nfit(idata) * weight(jqw,idata)
     $         /xinfo(idata))
          if (weight(jqw,idata).le.zero) weight(jqw,idata) = one
cc          print*, 'FIT WEIGHT: ', jqw, idata , weight(jqw,idata)
 450   continue 
c
c set k-window 
       if (winarr.eq.undef)  then
          call window2(dk1,dk2,xkmin,xkmax,lfs_kq,maxpts,
     $         lfs_qwindo)
       else

          jwin = get_array(winarr, name1(idata), 0,xwdat)
          if (jwin.le.0) then
             call echo(' feffit2: no window array?')
             return
          end if

          ntchi = jwin
          do 480 i = 1, ntchi
             lfs_qwindo(i) =  xwdat(i)
 480      continue 
       end if

c convert list of path indices to iup array
       call str2il(list(1:illen), max_pathindex, np,
     $      ipath_tmp, ier)

       iup= 0
       do i = 1, np
          if (u2ipth(ipath_tmp(i)).ge.1) then
             iup=iup+1
             iulist(iup,idata) = ipath_tmp(i)
          endif
       enddo
c       print *, 'iulist ', idata, iup, np
c       do i = 1, iup
c          print*, iulist(i,idata)
c       enddo
cc       print*,' idata ', idata, nrest_tmp
       nrestraint(idata) = nrest_tmp
       do 570 i = 1, nrest_tmp
          restraint(i,idata) = restr_tmp(i)
 570   continue 
cc       print*,' idata / nrestraint ', idata, idx, nrest_tmp

c===================================================================
c if not doing fit (multiple data-set fit), bail out now:
       if (.not.dofit)  then
          call echo('  feffit2: not fitting until all data '//
     $         'sets defined')
          return
       end if


c everything below assumes we're really, really doing a fit
c
c read the needed feff arrays
       call fefinp
c


c   set up fit arrays
c   sync variables
c   check # of variables  and # fit points
c   call lmdif1 and fiterr
       write(cmdstr, '(1x,a,i3,a)')  ' feffit2 fitting ',
     $         nfdats, ' data sets'
       call echo(cmdstr)

       mfit = (lfs_nqfit + nrestraint(idata)) * 2

       mfit = min(mfit, lenfvc)

cc       print*, '  Feffit ', mfit, idx, nfit(idx), nqwfs(idx), nfdats,
cc     $      nrestraint(idx)
cc       if (nrestraint(idx) .ge.1) then
cc          print*, '  Feffit using restraints:',nrestraint(idx)
cc       endif
c===================================================================
c  set up arrays for non-linear fit with lmdif / fitfun
c===================================================================
c
       if (nfit(idata).le.4) then
          call echo('   feffit2: too few fit data points: ')
          if (ifft(idata).eq.1) then
             write(cmdstr, 625)  '  Fit R range = [',xolow,
     $            ' : ', xohigh, ']'
          else
             write(cmdstr, 625)  '  Fit k range = [',xolow,
     $            ' : ', xohigh, ']'
          endif
 625      format (a,f8.3,a,f8.3,a)
          call echo(cmdstr)
          return
       endif
       if (nvarys.gt.nfit(idata)) then
          call echo(' *** feffit2: fewer fit data points '//
     $         'than variables! fit may not work.')
       endif


c synchronize the math expressions, set up fit variables
       call iff_sync
       if (toler .le. zero) toler = toldef
       do 700 i =1, nvarys
          varys(i) = scalar(i)
 700   continue
       do 710 i =1, mfit
          fvect(i) = zero
 710   continue
       do 720 i =1, lenwrk
          work(i) = zero
 720   continue
       do 740 i =1, mvarys
          iwork(i)  = 0
 740   continue


c  the real fit
       if (nvarys.gt.0.and.onlychi2.eqv..false.)  then
          final = .false.
          itera = 0
          call setsca('&fit_iteration', zero)
          call chrdmp('  fitting ... ')

          call fitfun2(mfit, nvarys, varys, fvect, lminfo)
c          goto 881

          call lmdif1 (fitfun2, mfit, nvarys, varys, fvect,
     $         toler, lminfo, iwork, work, lenwrk)

          call lm_err(lminfo,toler)
cc          if ((lminfo.ne.0).and.(lminfo.ne.4).and.(lminfo.ne.5)) then
          if (lminfo.ne.0) then
             call chrdmp('  estimating uncertainties ... ')
             call fiterr(fitfun2, mfit, nvarys, lenfvc, mvarys, fvect,
     $            ftemp, fjac, alpha, jprint, nerstp, varys,
     $            delta, correl, ier, ibadx)
             if (ier.eq.0) then 
                call echo('done.')
             else 
                call echo('finished with warnings')
                call echo('  error bars not estimated. some '//
     $               'variable(s) may not')
                call echo('  affect the fit. Check these variables:')
                do 880 i = 1, nvarys
                   if (ibadx(i).gt.0) then
                      ilen = istrln(scanam(i))
                      call echo('        '//scanam(i)(:ilen))
                   end if
 880            continue
             endif
          else
             call echo(' no uncertainties estimated!')
          end if
       else
          if (onlychi2.eqv..true.)  then
            call echo(' feffit2:  Calculating mychi_squared.')
          else
            call echo(' feffit2:  no variables defined.')
          endif
       endif
c 
c save scalars to Program Variables
 881   continue
       final = .true.
       call synvar

       call fitfun2(mfit, nvarys, varys, fvect, lminfo)




c       do i = 1,lfs_nqfit
c         tmparr2(i) = lfs_qwindo(i)/lfs_dchiq(i)
c       end do
       chisqr = sumsqr(fvect, mfit)!/sumsqr(tmparr2,lfs_nqfit)
c       write(*,*) sumsqr(fvect, mfit),sumsqr(tmparr2,lfs_nqfit)
       
c       chisqr = chi2(fvect,mfit,lfs_dchiq,lfs_qwindo,lfs_nqfit)

c add total number of independent points
c  avoid double counting by inspecting chi(k) names and sum-of-(chi*chi):
c  sure, this can be fooled by adding 1.e-8 to chi(k=0), but....
       xnidp = 0
       do i = 1,lfs_nqfit
         xnidp = xnidp + lfs_qwindo(i)
       end do

c       xnidp = lfs_nqfit

       chired = chisqr / max(0.1d0, (xnidp - nvarys))
      write(*,*)'xnidp=',xnidp,', nvarys=',nvarys

       if (onlychi2.eqv..true.) then
         call setsca('mychi_square',  chisqr)
         call setsca('mychi_reduced', chired)
         write(*,*) 'onlychi2 = true'
         return
       else
         call setsca('chi_square',  chisqr)
         call setsca('chi_reduced', chired)
         write(*,*) 'onlychi2 = false'
       endif
       
       write(*,*) 'chisqr',chisqr
       write(*,*) 'chired',chired
       write(*,*) 'fvect[1:2]',fvect(22),fvect(55)
       
       call setsca('r_factor',    rfact_total)
       call setsca('n_idp', xnidp)
    
       call setsca('kweight', qwfs(1,idata))    
       call setsca('kmin',   xkmin)
       call setsca('kmax',   xkmax)
       call setsca('dk1',    dk1)
       call setsca('dk2',    dk2)
       call setsca('rmin',   xrmin)
       call setsca('rmax',   xrmax)

       tmp = nvarys * one
       call setsca('n_varys', tmp)
      
       do 960 i = 1, nvarys
          delta(i) = sqrt(abs(chired)) * delta(i)
          ilen = istrln(scanam(i))
          messg = 'delta_' // scanam(i)(1:ilen)
          call setsca(messg, delta(i))
 960  continue 

c
c  fit with a macro is done, so reset stopping point
c  and previous macro pointer 
      if (ifit_mac.gt.0) then
         nmac_stop = 0
         imac      = imac_tmp
      endif
c

c save arrays to Program Variables for each data set
c
      do 3000 idata = 1, nfdats
          call set_array('res', name1(idata), lfs_resq, ntchi,1)
          call set_array('k',   name1(idata), lfs_kq,   ntchi,1)
          call set_array('chi', name1(idata), lfs_thiq, ntchi,1)
          call set_array('win', name1(idata), lfs_qwindo, ntchi,1)
c 
c--background chi(k)
          if (bkgfit(idata)) then
             do 2310 isp = 1, nbkg(idata)
                write(cmdstr, 2325)  'bkg',idata,'_',isp
                xspl(isp) = getsca(cmdstr)
 2310        continue
 2325        format (a,i2.2,a,i2.2)
             do 2330 i = 1, ntchi
                qx = lfs_kq(i)
                tmparr(i) =  bvalue(qknot(1,idata), 
     $               xspl,nbkg(idata),korder,qx,0)
                lfs_thiq(i) = lfs_thiq(i) - tmparr(i)
 2330        continue 
             call set_array('kbkg', name1(idata), tmparr, ntchi, 1)
          endif 
 3000  continue 

c  guess what?   we're done!
       return
       end










c---------------------------------------------------------------------------
       subroutine window2(dx1,dx2,xmin,xmax,xa,mpts, wa)

       implicit none
       integer mpts, iw, i, istrln
       double precision wa(mpts),xa(mpts), halfpi, zero, one, half, eps
       double precision x, x1, x2, x3, x4, xmin,xmax, xgrid, dx1, dx2
       double precision del1, del2, del12, del22
       double precision bessi0, bki0, bkav, bkde, bkde2, bkx, bkxx, bkom
       external bessi0, istrln
       parameter (halfpi= 1.570796326795d0, eps= 1.4d-5)
       parameter ( zero=0.d0, one=1.d0, half= 0.5d0) 

       del1 = dx1
       del12= dx1 * half
       del2 = dx2
       del22= dx1 * half
C       del22= dx2 * half

       x1 = xmin - del12
       x2 = xmin + del12
       x3 = xmax - del22
       x4 = xmax + del22


       do 10 i=1,mpts
          x = xa(i)
          if ((x.ge.x1).and.(x.le.x2)) then
             wa(i) = sin(halfpi*(x-x1) / (x2-x1)) ** 2
          elseif ((x.ge.x3).and.(x.le.x4)) then
             wa(i) = cos(halfpi*(x-x3) / (x4-x3)) ** 2
          elseif ((x.lt.x3).and.(x.gt.x2)) then
             wa(i) = one
          else
             wa(i) = zero
          endif
 10    continue

       return
       end

c--------------------------------------------------------

c       function chi2(fvect,mfit, dchi,win,n)
c       implicit none
c       integer n,i,mfit
c       double precision chi2, norm
c       double precision dchi(n), win(n), fvect(mfit),tmp(n)
c       external sumsqrcc

c       do i = 1,n
c         tmp(i) = win(i)/dchi(i)
c         write(*,*) tmp(i)
c       end do

c       chi2 = sumsqr(fvect, mfit)/sumsqr(tmp,n)
c       return
c       end
