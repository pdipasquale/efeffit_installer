      double precision function dist (r0, r1)
c     find distance between cartesian points r0 and r1
      implicit double precision (a-h, o-z)
      dimension r0(3), r1(3)
      dist = 0
      do 10  i = 1, 3
         dist = dist + (r0(i) - r1(i))**2
   10 continue
      dist = sqrt (dist)
      return
      end
