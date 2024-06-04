       subroutine fitfun2(mf, nx, xvar, fvec, iend)
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
c purpose:  feffit fitting function for lmdif1 minimization 
c
c arguments:
c      mf        number of points in ffit             [in]
c      nx        number of points in xv               [in]
c      xv        vector of variables                  [in]
c      fvec      output vector to be minimized        [out]
c      iend      information tag                      [out]
c Globals used as args:
c      final     is the fit done?                     [in]
c      nqfit     array of array sizes for regular k   [in]
c      lfs_nqfit sizes for real k array               [in]

c      thifit    vector of varibles to minamise       [out]

       implicit none
       include 'consts.h'
       include 'arrays.h'
       include 'keywrd.h'
       include 'fft.h'
       include 'fefdat.h'
       include 'pthpar.h'
       include 'feffit.h'
c local variables
       integer    nx, mf, iend, i,  id, isp, nx1, iter
       integer    jfit, nqdata, nfit1, nkpar, npaths, lfs_nqdata
       integer    iupath(mpaths), inxx, ixx, joff, jqw
       integer   ier, iff_eval_dp
       double precision xvar(nx), fvec(mf), xolow, xohigh, bvalue
       double precision getsca, xspl(mtknot), bx, qx, sum, tmpval
       double precision xsum, sumsqr, q, karr(maxpts)
       character  arg*2 
       integer    isarg, ilen, istrln, j, nmx, ipos, kweight, isumth
       external   bvalue, getsca,  istrln, iff_eval_dp, sumsqr
       save
c
       nkpar = 0
       id    = 1
       if (nx.ne.nvarys) iend = 1
       if (mf.ne.mfit)   iend = 2
       nx1 = nx - nbkg(id)
cc       print*, 'FITFUN N_vars:',nx,mf,nfdats
       do 20 i = 1, nx
          scalar(i) = xvar(i)
 20    continue
       call synvar
cc       print*, 'VARS ', xvar(1), xvar(2), xvar(3), xvar(4)
       if (final) rfact_total = zero
c
c  sum function to minimize over data sets
c   jfit is the counter (through all the data sets)
c   for the total number of fitting points

       jfit = 0

       nqdata = min(maxpts, max(2, nqfit(id)) + 10)
       lfs_nqdata = min(maxpts, max(2, lfs_nqfit) + 10)
 
       xolow  = qmin(id)
       xohigh = qmax(id)

       do 200 i = 1, lfs_nqfit+nrestraint(id)
          thifit(i)   = zero
 200   continue
c  re-initialize array for theoretical chi(k)
c  by assigning this to the background function
       do 300 i = 1, nqdata
          thiq(i, id) = zero
          thiqr(i, id) = zero
 300   continue

       do 350 i = 1, lfs_nqdata
          lfs_thiq(i) = zero
 350   continue
c===========================================
c
c sum over paths for theory chi for this data set
       npaths = 0
       do 400 i = 1, mpaths
          if (iulist(i,id).ne.0) then
             npaths = npaths + 1
             iupath(npaths) = iulist(i,id)
          end if
 400   continue 

       isumth = nqfit(id)+1
       call sum_paths(id, iupath, npaths, isumth,
     $      thiqr(1,id), thiq(1,id))


       do 410 i = 1, isumth
          karr(i) = (i-1) * qgrid
 410   continue


       ipos = 1
       do 420 i = 1, lfs_nqfit
          q = lfs_kq(i)

          call lintrp2(karr,thiq(1,id),
     $                 isumth, q,ipos,lfs_thiq(i))

 420   continue



c
c if refining background, include that here
       if ( bkgfit(id))  then
          do 450 isp = 1, nbkg(id)
             write(tmpstr, '(a,i2.2,a,i2.2)') 'bkg',id,'_',isp
             xspl(isp) = getsca(tmpstr)
cc                print*, 'FITFUN i ' , isp, xspl(isp)
 450      continue
          do 470 i = 1, lfs_nqfit
             qx = lfs_kq(i)
             lfs_thiq(i) = lfs_thiq(i) + 
     $           bvalue(qknot(1,id), xspl, nbkg(id),korder,qx,0)
 470      continue 
       end if

c subtract data from theory
       do 500 i = 1, lfs_nqfit
          lfs_resq(i) = (lfs_thiq(i)-lfs_chiq(i))
 500   continue 


c
c============= WRITE DATA TO FILE ===============
      open(unit=1,file=
     &    "./theory.dat")

      open(unit=2,file=
     &    "./experiment.dat")

      open(unit=3,file=
     &    "./residual.dat")

      open(unit=4,file=
     &    "./kwindow.dat")


      do 554 i=1, lfs_nqfit
         write(1,*)lfs_kq(i),	lfs_thiq(i)
         write(2,*)lfs_kq(i),	lfs_chiq(i), 	lfs_dchiq(i)
         write(3,*)lfs_kq(i), 	lfs_resq(i)
         write(4,*)lfs_kq(i),	lfs_qwindo(i)
 554   continue

      close(1)
      close(2)
      close(3)
      close(4)

c============================================================



       kweight = lfs_kweight

c   take fft of theory(+bkg) - data
c   note: imag part of theory chi(k) compares to real part of data!!!
c   see ff2chi as well
c       do 500 i = 1, lfs_nqdata
c          thifit(i) = lfs_thiq(i,id)-lfs_chiq(i,id)
c 500   continue


       do 600 i = 1, lfs_nqfit
      if(lfs_kq(i).ne.0) then
          thifit(i) = lfs_resq(i)*lfs_qwindo(i)*
     $                lfs_kq(i)**kweight
     $                /(lfs_dchiq(i)*lfs_kq(i)**kweight)
      else
      thifit(i)=0
      end if

 600   continue

c       call fitfft(thiq(1,id), maxpts, maxfft, wfftc, qgrid,
c     $      qwindo(1,id), qwfs(jqw,id), rwindo(1,id), one,
c     $      ifft(id),xolow,xohigh, nfit1, thifit)
c
c       if (nfit1.ne.nfit(id)) then
c          call echo(' fitfun2 fitfft failed internal test.')
c          iend = -10
c       end if


c  evaluate the contribution to fvec for this data set.  weight scales
c  chi-square properly to the number of independent points. this is
c  important for error analysis (if chi-square is to increase by one,
c  it  must be scaled correctly.), but only in the final pass, when
c  chi-square and r-factors will be calculated.
c             print*, 'fitfun ', jfit, jqw, id,
c     $            weight(jqw, id), nfit(id), thifit(1), thifit(2)
       do 700 i =  1, lfs_nqfit
          fvec(i) = thifit(i)
 700   continue

       jfit = lfs_nqfit

c construct r-factor
       if (final) then 

          do 750 i = 1, lfs_nqfit
      if(lfs_kq(i).ne.0) then
             thifit(i) = lfs_resq(i)*lfs_qwindo(i)*
     $                   lfs_kq(i)**kweight
     $                   /(lfs_dchiq(i)*lfs_kq(i)**kweight)
      else
      thifit(i)=0
      end if
 750      continue
 
 

c          call fitfft(chiq(1,id), maxpts, maxfft, wfftc, qgrid,
c     $         qwindo(1,id), qwfs(jqw,id), rwindo(1,id), one, 
c     $         ifft(id),xolow,xohigh, nfit1, chifit)
c          if (nfit1.ne.nfit(id)) then
c             call echo(' fitfun2 fitfft failed internal test.')
c             iend = -10
c          end if

          sum        = zero
          rfactr(id) = zero
          do 800 i = 1, nfit(id)
             sum        = sum        + chifit(i)*chifit(i)
             rfactr(id) = rfactr(id) + thifit(i)*thifit(i)
 800      continue
          if (sum.le.tiny) sum = tiny
          rfactr(id)  =  rfactr(id) / (sum * nqwfs(id))
          rfact_total = rfact_total + rfactr(id) 
cc                print*, ' final ', sum, rfactr(id), id, rfact_total
       end if
c
c  restraints
       if (nrestraint(id).ge.1) then
cc                print*, ' RESTRAINT ', nrestraint(id), id, jfit
          do 900 i =  1, nrestraint(id)
             if ( (restraint(i,id) .ne. undef).and.
     $            (restraint(i,id) .ne. '')) then 
                ier    = iff_eval_dp(restraint(i,id),tmpval)
                if (ier.eq.0) then
                   jfit   = jfit + 1
                   fvec(jfit) = tmpval
cc                         print*, ' t: ', tmpval, jfit, ier
                endif
             endif
 900      continue
       endif
cc          print*, ' end of ndata loop  ', id, nfit(id), jfit

       return
c end subroutine fitfun2
       end



c----------------------------------------------------------------
c  make sure ip is in range
c  find ip such that   x(ip) <= xin <= x(ip+1)
c  Assumes there is an i s.t. x(i) <= xin
       subroutine hunt2(x, npts, xin, ip)
       implicit none
       integer    npts, ip
       double precision   x(npts), xin
       do while (x(ip).lt.xin)
         if (ip.ge.npts) then
           ip = npts
           return
         end if
         ip=ip+1
       end do

       if (x(ip).ge.xin.and.ip.gt.1) ip=ip-1

       return
       end

c-------------------------------------------------------------------
       subroutine lintrp2(x, y, npts, xin, ip, yout)
c
c    linear interpolation for use in loops where xin increases 
c    steadily through the monotonically increasing array x. 
c  arguments:
c     x      array of ordinate values                   [in]
c     y      array of abscissa values                   [in]
c     npts   length of arrays x and y                   [in]
c     xin    value of x at which to interpolate         [in]
c     ip     index such that x(ip) <= xin <= x(ip+1)    [in/out]
c     y      interpolated abscissa at xin               [out]
c  note: this routine is called extremely often 
c        -- anything to improve efficiency should be done
       implicit none
       integer    npts, ip
       double precision   x(*), y(*), tiny, xin, yout
       parameter  (tiny = 1.d-9)
c  find ip such that   x(ip) <= xin < x(ip+1)
       call hunt2(x, npts, xin, ip)
       yout = y(ip)
       if ((x(ip+1)-x(ip)) .gt. tiny)  yout = yout +
     $     (y(ip+1)-y(ip)) * (xin-x(ip)) / (x(ip+1)-x(ip))
       return
c  end subroutine lintrp
       end
