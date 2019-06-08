      subroutine Metreader

!      subroutine that reads meteorological data files
!      The data must have three header lines

      use Module1
      use precis_param

      implicit none

      integer      :: ilayer,lengthz
      character*1  :: answer
      character*3  :: metfileformat                                     !equals either 'old' or 'uwy'
      character*80 :: lineinput, line1, line2, line3, line4, line5, line6, line7
      character*10 :: zval
      real(kind=ip)         :: DwptAirLayer(2000), garbage1, garbage2, psat
      integer               :: ios

      metfileformat = 'old'           !assume default file type=old type
!     Read data from file

      write(6,*) 'Opening ',Metfile
      open(unit=20,file=Metfile,status='old')
      !open(unit=20,file=Metfile,status='old',err=1000)

      ilayer = 1

      !If the sounding was copied from the text output of the NOAA READY web site;
      !http://ready.arl.noaa.gov/READYcmet.php
      !The first three lines will be header lines as follows:
      !   PRESS  HGT(MSL)  TEMP  DEWPT  WNDDIR  WNDSPEED
      !   HPA  M  C  C  DEG  M/S
      !   E Estimated surface height
      !
      !skip the first four lines

      !If the sounding came from the Univ. Wyoming radiosonde database it should have the format:
      !04018 BIKF Keflavikurflugvollur Observations at 12Z 15 Apr 2010
      !-----------------------------------------------------------------------------
      !   PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV
      !    hPa     m      C      C      %    g/kg    deg   knot     K      K      K 
      !-----------------------------------------------------------------------------

      !In that case we skip the first five lines

      read(20,'(a80)',err=400,end=500) line1
      read(20,'(a80)') line2
      read(20,'(a80)') line3
      read(20,'(a80)') line4

      !See whether we're reading the old-fashioned metfiles or the standard radiosonde
      !output from the Univ. Wyoming
      if (line2(1:10).eq.'----------') then
         !Looks like a Univ. Wyoming radiosonde with the table header starting on line 2
         metfileformat = 'uwy'                 !assume we're reading from a Univ. Wyoming radiosonde format
         write(6,*) 'metfile format=Univ. Wyoming radiosonde'
         write(6,*) line1
         write(9,*) 'metfile format=Univ. Wyoming radiosonde'
         write(9,*) line1
         read(20,'(a80)') line5
         if (line5(1:10).ne.'----------') go to 600  !Doesn't follow the pattern of table header
       else if (line3(1:10).eq.'----------') then
         !Looks like a Univ. Wyoming radiosonde with the table header starting on line 3
         metfileformat = 'uwy'                 !assume we're reading from a Univ. Wyoming radiosonde format
         write(6,*) 'metfile format=Univ. Wyoming radiosonde'
         write(6,*) line1
         write(6,*) line2
         write(9,*) 'metfile format=Univ. Wyoming radiosonde'
         write(9,*) line1
         write(9,*) line2
         read(20,'(a80)') line5
         read(20,'(a80)') line6
         if (line6(1:10).ne.'----------') go to 600  !Doesn't follow the pattern of table header
       else if (line4(1:10).eq.'----------') then
         !Looks like a Univ. Wyoming radiosonde with the table header starting on line 4
         metfileformat = 'uwy'                 !assume we're reading from a Univ. Wyoming radiosonde format
         write(6,*) 'metfile format=Univ. Wyoming radiosonde'
         write(6,*) line1
         write(6,*) line2
         write(9,*) 'metfile format=Univ. Wyoming radiosonde'
         write(9,*) line1
         write(9,*) line2
         read(20,'(a80)') line5
         read(20,'(a80)') line6
         read(20,'(a80)') line7
         if (line7(1:10).ne.'----------') go to 600  !Doesn't follow the pattern of table header
       else
         write(6,*) 'metfile format=old-style plumeria'
      end if

      write(6,1)
      write(9,1)
1     format('layer     p (hPa)  z a.s.l. (m)       T (C) dew point (C)        rh    wind dir     wnd spd',/, &
             '                                                                     deg CW of N        m/s')

       !Read lines from sounding.
       !After the height variable, there may be an "E" for estimated height.
       !So, we may need to use one of two read statements

100    continue
       !Assume we're reading from the old-fashoned metfiles.
       if (metfileformat.eq.'old') then
          read(20,'(a80)',end=999) lineinput
          read(lineinput,*,err=400,end=999) pAirLayer(ilayer), zval, TairLayer(ilayer), &
                     DwptAirLayer(ilayer), WinddirLayer(ilayer), WindspeedLayer(ilayer)
          lengthz = index(zval,'E')
          if (DwptAirLayer(ilayer).lt.-270.0)  then
              DwptAirLayer(ilayer)=-270.0
              write(6,*) 'Dewpoint was <-270 C.  Changed to -270.'
          end if

          if (lengthz.eq.0) then            !if there's not an "E" in zval, read z directly
             read(zval,*) ZairLayer(ilayer)
           else                             !if there's an "E", strip it out
             read(zval(1:lengthz-1),*) ZairLayer(ilayer)
          end if
         else
          read(20,'(a80)',end=999) lineinput
          if (lineinput(1:19).eq.'Station information') go to 999
          read(lineinput,*,err=400,end=999) pAirLayer(ilayer), ZairLayer(ilayer), &
                     TairLayer(ilayer), DwptAirLayer(ilayer), garbage1, garbage2, &
                     WinddirLayer(ilayer), WindspeedLayer(ilayer)
           WindSpeedLayer(ilayer) = WindSpeedLayer(ilayer)*0.51444 !convert from knots to m/s
       end if
       TairLayer(ilayer) = TAirLayer(ilayer)+273.15          !convert to Kelvin
       pAirLayer(ilayer) = pAirLayer(ilayer)*100.            !convert to Pascals
       HumidAirLayer(ilayer) = psat(DwptAirLayer(ilayer) + 273.15)/psat(TairLayer(ilayer))
       write(6,2) ilayer, pAirLayer(ilayer)/100.,ZairLayer(ilayer),TairLayer(ilayer), &
                          DwptAirLayer(ilayer)+273.15,HumidAirLayer(ilayer), &
                          WinddirLayer(ilayer), WindspeedLayer(ilayer)
       write(9,2) ilayer, pAirLayer(ilayer)/100.,ZairLayer(ilayer),TairLayer(ilayer), &
                          DwptAirLayer(ilayer)+273.15,HumidAirLayer(ilayer), &
                          WinddirLayer(ilayer), WindspeedLayer(ilayer)
2      format(i5,7f12.4)
       ilayer = ilayer + 1
       goto 100

       !once we're at the bottom of the file
999    iatmlayers = ilayer-1

       return

       !Error traps
400    write(6,*) 'error reading from meteorological input file.  Program stopped'
       write(9,*) 'error reading from meteorological input file.  Program stopped'
       close(9)
       stop
500    write(6,*) 'error reading from meteorological input file.  End of file reached'
       write(9,*) 'error reading from meteorological input file.  End of file reached'
       close(9)
       stop

600    write(6,*) 'Error:  The header lines of the input file should either follow this format:'
       write(6,*) 'PRESS  HGT(MSL)  TEMP  DEWPT  WNDDIR  WNDSPEED'
       write(6,*) 'HPA  M  C  C  DEG  M/S'
       write(6,*) 'E Estimated surface height'
       write(6,*) 
       write(6,*) '. . . . or this format:'
       write(6,*) '04018 BIKF Keflavikurflugvollur Observations at 12Z 15 Apr 2010'
       write(6,*) '-----------------------------------------------------------------------------'
       write(6,*) '   PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV'
       write(6,*) '    hPa     m      C      C      %    g/kg    deg   knot     K      K      K '
       write(6,*) '-----------------------------------------------------------------------------'
       write(6,*)
       write(6,*) 'The first lines of your header file are:'
       write(6,*) line1
       write(6,*) line2
       write(6,*) line3
       write(6,*) line4
       write(6,*) line5
       write(6,*) 'Program stopped'
       stop
  
1000   write(6,*) 'Error: Meteorological input file does not exist.  Program stopped'
       write(9,*) 'Error: Meteorological input file does not exist.  Program stopped'
       close(9)
       stop

       end subroutine Metreader

