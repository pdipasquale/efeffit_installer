       program diffkk 
c 
c  generate f' and f'' near x-ray resonances for an atom, including 
c  fine-structure due to solid-state effects (ie, xafs/dafs). 
c
c algorithm summary:
c 1 the brennan-cowan implementation of the cromer-libermann (cl) 
c   calculation is used as a starting set of a causal f' and f''. 
c   modifications were made to the bc code, mostly to make it easier 
c   to use, more closely f77 conforming, and smaller.  this data is
c   convolved with a lorenztian, typically with a width of a few ev.
c
c 2 an externally supplied file containing the xafs mu(e) is used to 
c   improve f''.  special support for xmu.dat files from feff is 
c   provided, or the external file can contain *measured* mu(e) for
c   the system of interest.  a simple matching procedure is done to 
c   make the supplied mu(e) match the f'' from cl.  
c
c 3 a differential kramers-kronig transform is used to convert the 
c   changes in f'' (ie f''_supplied - f''_cl) into the changes in  f' 
c   (ie f' - f'_cl).   the kk transform is done using a maclaurin 
c   series method, as suggested in the literature. 
c
c  the result is a causal pair of f' and f'' that reflect the 
c  presence of the atoms neighboring the central atom.
c
c  --  Further notes on the algorithms used and instructions for  --
c  --  program use are given in the program documentation.        --
c
c  copyright 1997,...,2003  matt newville
c  
c  acknowledgements: julie cross, chuck bouldin, john rehr, and
c            bruce ravel contributed to the design of this code.  
c            the best ideas were theirs.  all mistakes are mine. 
c  
c  1.21 using ifeffit atomic data for core-width if ewidth=0 on input
c  
       include "dkcom.f"
       include "../lib/ifeffit.inc"
       integer  i, ncol, mcol, ierr, ilen
       parameter (mcol = 6)
       double precision  df1(mpts), df2(mpts)
       double precision  o1(mpts), o2(mpts),  atz
       character*128     str
       integer  istrln, guess_iz
       external istrln, guess_iz

       versn = '1.3'
       ncol  = mcol
       ndoc  = mdoc
       npts  = mpts
c print version number
       i = ifeffit(" ")
       write(str,'(3a)') ' --  diffkk version ', versn, '--'
       ilen = istrln(str)
       call messag(str(1:ilen))
c read diffk.inp 
       call dkinp
c read mu(e) data, convert to f''(e) on an even energy grid
       write(str,'(2a)')  ' Reading experimental data: ',
     $      xmufil
       ilen = istrln(str)
       call messag(str(1:ilen))

       write(str,'(3a)')  'read_data(file=',
     $      xmufil(1:istrln(xmufil)),',group=dat,type=raw)'
       i = ifeffit(str)
       i = iffgetstr('column_label',str)
       if (str(1:13).eq.'--undefined--') then
          call messag( '  diffkk: fatal error.')
          goto 999
       endif
c
       
       write(str,'(a,i1)')  'set dat.energy = dat.',iencol
       if (iencol .ge.10) then
          write(str,'(a,i2)')  'set dat.energy = dat.',iencol          
       endif
       i = ifeffit(str)
       write(str,'(a,i1)')  'set dat.expdat = dat.',imucol
       if (imucol .ge.10) then
          write(str,'(a,i2)')  'set dat.expdat = dat.',imucol          
       endif
       i    = ifeffit(str)
       i    = iffgetarr("dat.energy", energy)
       npts = iffgetarr("dat.expdat", expdat)
       iatz = guess_iz(energy, expdat, npts, e0)
       write(str,'(a,g12.7)')  'set egrid = ', egrid 
       i = ifeffit(str)
       i = ifeffit('set elow  = floor(dat.energy)')
       i = ifeffit('set ehigh = ceil(dat.energy)')
       i = ifeffit('set d.energy = range(elow,ehigh,egrid)')
       i = ifeffit('set d.expdat = splint(dat.energy,'//
     $      'dat.expdat,d.energy)')
       i = iffgetsca('elow',  elow)
       i = iffgetsca('ehigh', ehigh)
       i = iffgetarr("d.energy", energy)
       if (f2tof1) then 
          i   = ifeffit('set d.expdat = d.expdat * '//
     $         'd.energy/(max(1,e0))')
          npts = iffgetarr("d.expdat",    expdat)
       else
          npts = iffgetarr("d.expdat",   expdat)
       endif
cc       print*,' npts =', npts
cc       i = ifeffit ('show @arrays')
c
c generate initial tables of f' f'' on the same energy grid
       atz= iatz * 1.d0
       i =  iffputsca("iz", atz)
       if (ewidth.gt.0.d0) then
          i =  iffputsca("wid", ewidth)
          i =  ifeffit("f1f2(energy=d.energy, z=iz,width=wid)")
       else
          i =  ifeffit("f1f2(energy=d.energy, z=iz)")
       endif
c broaden f1cl and f2cl
       i =  iffgetarr("d.f1", f1cl)
       i =  iffgetarr("d.f2", f2cl)
       i =  iffgetsca('core_width', ewidth)
       i =  ifeffit("rename d.f1 d.f1_cl")
       i =  ifeffit("rename d.f2 d.f2_cl")
cc       i =  ifeffit("show @arrays")
c
c align/shift/scale improved f''  to tabulated f''(df2 = expdat - f2cl)
       call chrdmp(' Matching experimental data to')
       if (f2tof1) then 
          call messag(' Cromer-Liberman f'''' ')
       else
          call messag(' Cromer-Liberman f'' ')
       end if
cc       print*,  ' f2tof1 ', f2tof1, npts
       call dkfit
       call dkfcn(npts,numvar,xvarys,df2,ierr)
c do kk transform of delta f''  -> delta f'
       call messag(' Doing difference Kramers-Kronig transform')
       if (f2tof1) then
          call kkmclr(npts, energy, df2, df1)
       else
          call kkmclf(npts, energy, df2, df1)
       endif
c add delta f' to initial f', delta f''  to initial f''
       do 100 i = 1, npts
          o1(i)  = df1(i) + f1cl(i)
          o2(i)  = df2(i) + f2cl(i)
 100   continue 
       i =  iffputarr("d.f1", npts, o1)
       i =  iffputarr("d.f2", npts, o2)
c play with docs (put user titles at top,
c      keep as many old doc lines as fit)
       i =  iffputstr('doc', title(1))

       if (active)  then
          call messag(' ')
          call messag('  Ready to write out data file:')
          call askstr('** output file name',outfil)
       end if

       write(str,'(3a)')  'write_data(file=',
     $      outfil(1:istrln(outfil)),
     $      ',d.energy,d.f1,d.f2,d.f1_cl,d.f2_cl,$doc)'

       i =  ifeffit(str)
       i =  ifeffit('save diffkk.sav')
       call messag(' writing summary to diffkk.log')
       call dklog
       call messag(' -- diffkk done -- ')
 999   continue 
       end
