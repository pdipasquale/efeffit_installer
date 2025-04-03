       subroutine getcom(jinit, line)
c
c  return next "real" command line from input file(s)
c    -  allows use of "include file" or "load file" for reading
c       from other files, and manages the set of include files
c    -  checks for and ignores comment lines and blank lines.
c    -  opens and closes all input files, including initial file.
c
c   jinit  initialization flag              [in]
c   line   next command line to parse   [in/out]
c
c notes:
c   1. to initialize, set jinit<0 and line= input_file_name.
c      if line=' ', commands will be read from standard input
c      (unit 5).
c   2. returned line will be sent through triml and untab.
c   3. uses routine iscomm to test if line is a comment line.
c   4. uses routine openfl to open files (which include automatic
c      assignment of next available unit number)
c   5. special returned values:
c        'getcom_end'  = done reading all inputs
c        'getcom_error'= an error has occurred. the calling routine
c                        should probably stop
c        'getcom_nofile'= on initialization, the file named by "line"
c                         could not be found
c matt newville march 1997
       implicit none
       integer mwords, ilen, i, jinit, mfil, nfil
       character*(*) line, stat*8
       parameter (mwords=2, mfil=10, stat = 'old')
       character*128  files(mfil), errmsg, words(mwords)
       integer   iunit(mfil), istrln, nwords, ierr, iex
       integer   iret, iread
       logical   iscomm
       external  istrln, iscomm, iread
       save      files, iunit, nfil
c
       if ((jinit.lt.0)) then
          jinit  = 1
          do 10 i = 1, mfil
             iunit(i) = 0
             files(i) = ' '
 10       continue
          nfil     = 1
          files(1) = line
          call triml(files(1))
          if (files(1) .eq. ' ') then
             iunit(1) = 5
          else
             call openfl(iunit(1), files(1), stat, iex, ierr)
             if (iex.lt.0) then
                line = 'getcom_nofile'
                return
             elseif (ierr.ne.0) then
                line = 'getcom_error'
                return
             end if
          end if
       end if
c  read next line from current input file
 100   continue
       line   = ' '
       iret = iread(iunit(nfil), line)
       if (iret.eq.-1) goto  500
       if (iret.eq.-2) goto 1000
       if (iret.eq. 0) goto  100
c
c  check if command line is 'include filename'.
c  if so, open that file, and put it in the files stack

       call triml(line)
       if (iscomm(line)) go to 100
       nwords = mwords
       words(2) = ' '
       call bwords(line, nwords, words)
       call lower(words(1))
       if (((words(1) .eq. 'include').or.(words(1) .eq. 'load'))
     $      .and. (nwords .gt. 1)) then
          nfil = nfil + 1
          if (nfil .gt. mfil) go to 2000
          call getfln(words(2), files(nfil), ierr)
          if (ierr. ne. 0) go to 2400
c  test for recursion:
          do 400 i = 1, nfil - 1
             if (files(nfil) .eq. files(i)) go to 3000
 400      continue
          call openfl(iunit(nfil), files(nfil), stat, iex, ierr)
          if (iex .lt. 0) go to 2600
          if (ierr.lt. 0) go to 2800
          go to 100
       end if
       return
c
c  end-of-file for command line file: drop nfil by 1,
c  return to get another command line
 500   continue
       if (iunit(nfil) .ne. 5) close(iunit(nfil))
       iunit(nfil) = 0
       files(nfil) = ' '
       nfil = nfil - 1
       if (nfil.gt.0) go to 100
       line = 'getcom_end'
       return
c   error messages
 1000  continue
       call echo(' # getcom error: general read error')
       go to 4500
 2000  continue
       call echo(' # getcom error: too many nested "include"s')
       write(errmsg, '(1x,a,i3)') ' # current limit is ', mfil
       ilen  = istrln(errmsg)
       call echo(errmsg(1:ilen))
       go to 4500
 2400  continue
       call echo(' # getcom error: cannot determine "include" file')
       go to 4500
 2600  continue
       call echo(' # getcom error: cannot find "include"d file')
       go to 4500
 2800  continue
       call echo(' # getcom error: cannot open "include"d file')
       go to 4500
 3000  continue
       call echo(' # getcom error: recursive "include" of file')
       go to 4500
 4500  continue
       errmsg = ' # reading file: '//files(nfil)
       if (files(nfil) .eq. ' ')
     $      errmsg = ' # reading from standard input'
       ilen   = istrln(errmsg)
       call echo(errmsg(1:ilen) )
       line = 'getcom_error'
       return
c end subroutine getcom
       end
