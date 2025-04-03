      subroutine trap (x, y, n, sum)
      implicit double precision (a-h, o-z)

c     Trapeziodal integration of y(x), result in sum

      dimension x(n), y(n)

      sum = y(1) * (x(2) - x(1))
      do 10  i = 2, n-1
         sum = sum + y(i) * (x(i+1) - x(i-1))
   10 continue
      sum = sum + y(n) * (x(n) - x(n-1))
      sum = sum/2

      return
      end
