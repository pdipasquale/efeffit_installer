      subroutine exjlnl (z, l, jl, nl)

c     purpose:  to calculate the spherical bessel functions jl and nl
c               for l = 0, 1, 2 or 3 using exact analytic expression
c
c     arguments:
c       z = argument of jl and nl
c       l = integer order of spherical bessel function
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c            Note that this nl = abramowitz yl.
c
c       analytic expressions from abramowitz 10.1.11 and 10.1.12
c       recurrence relation to get analytic j4,n4  eqns 10.1.19-22

      implicit double precision (a-h, o-z)

      complex*16 z, jl, nl

      complex*16 cosz, sinz

c     Exact formulae unstable for very small z, so use series
c     expansion there.  Limit of .3 chosen for 9 digit agreement.
      if (abs(z) .lt. 0.3)  then
         call bjnser (z, l, jl, nl, 0)
      else
c        use analytic formulae
         cosz = cos(z)
         sinz = sin(z)

         if (l .eq. 0)  then
            jl =  sinz / z
            nl = -cosz / z

         elseif (l .eq. 1)  then
            jl =  sinz/z**2 - cosz/z
            nl = -cosz/z**2 - sinz/z

         elseif (l .eq. 2)  then
            jl = ( 3/z**3 - 1/z)*sinz - 3*cosz/z**2
            nl = (-3/z**3 + 1/z)*cosz - 3*sinz/z**2

         elseif (l .eq. 3)  then
            jl = ( 15/z**4 - 6/z**2)*sinz + (-15/z**3 + 1/z)*cosz
            nl = (-15/z**4 + 6/z**2)*cosz + (-15/z**3 + 1/z)*sinz

         else
            call fstop(' at exjlnl, l out of range')

         endif
      endif

      return
      end
