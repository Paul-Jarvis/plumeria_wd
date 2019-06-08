       subroutine read_input

!      subroutine that reads from the input file

       use Module1
       use Module2
       use precis_param

       implicit none
       character  :: infile*60, yesno*3
       real(kind=ip)       :: diam
       character(len=80)   :: inputline
       character(len=3)    :: gastype                     !='air' if magmatic gas is air
       integer             :: ioerr, iargc, nargs         !number of program arguments
       real(kind=ip)       :: AirTemp
       real(kind=ip)       :: wind, windnow

       !SET DEFAULT VALUES
       
       ReadMetFile = .false. 
       !do i=1,50                                            !calculate mass of each particle size
       !   diam     = real(20*i)
       !   pmass(i) = rho_m*pi*diam**3/6.
       !end do

        !OPEN LOG FILE, WRITE CURRENT DATE & TIME
        open(unit=9,file='plumeria.lst')        !open log file
        call date_and_time(date,time2,zone,values)  !get current date & time
        !write(6,1) values(1), values(2), values(3), values(5), values(6)  !write date & time
        write(9,1) values(1), values(2), values(3), values(5), values(6)
1       format(4x,'Plumeria run ',i4,'.',i2.2,'.',i2.2,i4,':',i2.2,/)

!       TEST READ COMMAND LINE ARGUMENTS
        nargs = iargc()
        if (nargs.eq.1) then                !If there's 1 command-line argument,
                                            ! it's the input file name.  Read it.
           call getarg(1,infile)
         else
            write(6,*)'Enter name of ESP input file:$'
           read(5,*) infile
       end if

       write(6,*)  'opening input file ',infile
       write(9,*)  'opening input file ',infile
       write(6,*)
       write(9,*)
       open(unit=10,file=infile,status='old',err=2000)

       read(10,2)
2      format(2/)                     !skip some lines
       read(10,*)

       !READ OUTPUT FILE NAME, STDOUT OPTIONS
       read(10,'(a80)') OutFileName
       open(unit=11,file=OutFileName)
       write(6,*) 'output file name=',OutFileName
       write(9,*) 'output file name=',OutFileName


       !DETERMINE WHETHER TO READ METFILE
       read(10,2)                       !skip some lines
       read(10,*)
       write(6,*) 'reading whether to read from met. file'
       write(9,*) 'reading whether to read from met. file'
       read(10,'(a3)',err=1000) yesno
       !IF WE'RE SUPPOSED TO READ FROM THE MET. FILE . . .
       if (yesno.eq.'yes') then
           ReadMetFile = .true.
           write(6,3)
           write(9,3)
3          format('Reading atmospheric data from file',/,a50)
           read (10,'(a80)') Metfile       !read metfile filename
           call Metreader                  !call subroutine that read the metfile
           read(10,7)                      !skip over the lines in the input file with the simplified atmosphere.
7          format(8/)
           read(10,*)
        !IF NOT . . .
        else if (yesno.eq.'no') then
          iatmlayers = 1
          write(6,*)  'No external met. file input.    Reading atmospheric parameters'
          write(9,*)  'No external met. file input.    Reading atmospheric parameters'
          read(10,2)
          read(10,*)
          read(10,*,err=1005) T_air(1)         !atm. temperature at vent elevation
          T_air(1) = T_air(1)+273.15_ip        !convert from C to K
          write(6,9) T_air(1)
          write(9,9) T_air(1)
          read(10,*,err=1000) Humidity
          Humidity = Humidity/100._ip          !convert from percent to fraction
          write(6,4) Humidity
          write(9,4) Humidity
          read(10,*,err=1001) dTdz             !lapse rate in troposphere (K/m)
          write(6,10) dTdz
          write(9,10) dTdz
          read(10,*,err=1002) z_trop          !tropopause elevation (m asl)
          write(6,11) z_trop
          write(9,11) z_trop
          read(10,*,err=1003) H_trop           !tropopause thickness (m)
          write(6,12) H_trop
          write(9,12) H_trop
          read(10,*,err=1004) dTdz_strat       !lapse rate in stratosphere (K/m)
          write(6,13) dTdz_strat
          write(9,13) dTdz_strat
          read(10,'(a80)') inputline
             read(inputline,*,err=1200) windconst, windslope, winddir   !see if windslope is included
             go to 1201
1200         read(inputline,*,err=1015) windconst, winddir     !If not, just read wind speed (m/s), direction (deg. E of N)
             windslope = 0.0
1201      write(6,14) windconst, windslope, winddir
          write(9,14) windconst, windslope, winddir
          T_trop = T_air(1) + dTdz * z_trop_above_vent    !temperature at tropopause
9         format('                Air temperature at vent (K)=',f8.1)
4         format('                                Humidity, %=',f8.1)
10        format('                  dTdz in troposphere (K/m)=',f8.4)
11        format('               Tropopause elevation (m asl)=',f8.0)
12        format('                    Tropopause thicknes (m)=',f8.0)
13        format('                 dTdz in stratosphere (K/m)=',f8.4)
14        format('               wind speed at sea level(m/s)=',f8.2,/, &
                 'increase in wind speed with elevation (1/s)=',f8.4,/, &
                 '                wind direction (deg E of N)=',f8.2)
        else
          write(6,*) 'Error: on line 10 of the input file you must state either "yes" or "no"'
          write(6,*) 'regarding whether to read a meteorological input file.'
          write(6,*) 'You gave: ',yesno
          write(6,*) 'Program stopped.'
          stop 1
       end if

       !READ VENT PROPERTIES
       write(6,*)   'Reading vent properties'
       write(9,*)   'Reading vent properties'
       read(10,2)                           !skip some lines
       read(10,*,err=1016) vent_elevation   !vent elevation, m
       write(6,15) vent_elevation
       write(9,15) vent_elevation
       if (ReadMetFile) then 
           T_air(1)=AirTemp(0.)
           write(6,16) T_air(1)
           write(9,16) T_air(1)
          else
            z_trop_above_vent = z_trop-vent_elevation      !calculate height of trop above the vent
            T_asl  = T_air(1)-dTdz*vent_elevation !calculate air temp. at sea level
       end if
       read(10,*,err=1006) diam             !vent diameter, m
       write(6,17) diam
       write(9,17) diam
       r(1) = diam/2.                       !vent radius
       read(10,*,err=1007) u(1)             !exit velocity (m/s)
       write(6,18) u(1)
       write(9,18) u(1)
       read(10,*,err=1008) mw               !mass fraction added water
       write(6,19) mw
       write(9,19) mw
15      format('         vent elevation (m asl)=',f8.0)
16      format('    Air temperature at vent (K)=',f8.2)
17      format('              vent diameter (m)=',f8.1)
18      format('            exit velocity (m/s)=',f8.1)
19      format('      mass fraction added water=',f8.3)

       !READ MAGMA PROPERTIES
       write(6,*)    'Reading magma properties'
       write(9,*)    'Reading magma properties'
       read(10,2)                           !skip some lines
       read(10,*,err=1009) T_mag            !magma temperature (Celsius)
       T_mag = T_mag+273.15                 !convert from C to K
       read(10,'(a80)')   inputline
           read(inputline,*,iostat=ioerr) n_0air, gastype !try and see if the air is gas
           if ((ioerr.eq.0).and.(gastype.eq.'air')) then
                n_0 = 0.
             else
2001            read(inputline,*,err=1010) n_0              !mass fraction gas
                n_0air  = 0.001
                gastype = 'H2O'
           end if
       read(10,*,err=1011) Cp_m             !magma specific heat, J/kg K
       read(10,*,err=1012) rho_m            !magma density, kg/m3
       write(6,6)   T_mag, n_0, n_0air, Cp_m, rho_m
       write(9,6)   T_mag, n_0, n_0air, Cp_m, rho_m
6      format('       magma temperature (K)=',f6.1,/, &
              '       mass fraction H2O gas=',f7.4,/, &
              '           mass fraction air=',f7.4,/, &
              'magma specific heat (J/kg K)=',f6.1,/, &
              '       magma density (kg/m3)=',f6.1)

!       !READ GRAIN SIZE DISTRIBUTION
!       write(9,6) 
!       write(6,7) 
!7      format(/,'Grain sizes:',/,'size (um)    mass fraction')
!       read(10,2)
!       do i=1,50
!          read(10,*,err=1015) garbage1, garbage2, psize(1,i)
!          if (psize(1,i).gt.0.) then
!             write(9,8) garbage2, psize(1,i)
!             write(6,8) garbage2, psize(1,i)
!8            format(i8,f16.3)
!          end if
!       end do

!       !MAKE SURE MASS FRACTIONS ADD TO 1
!       write(6,*)  'Checking to make sure mass fractions add to 1'
!       write(9,*)  'Checking to make sure mass fractions add to 1'
!       if (abs(sum(psize(1,1:50))-1.0).gt.0.01) then
!          write(6,9)
!          write(9,9)
!9         format('Sum of mass fractions of grain sizes is not equal',/,&
!                  'to 1.  Program stopped.')        
!        else
!          write(6,*)  'Good.  They do.'
!          write(9,*)  'Good.  They do.'
!       end if
       
       close(10)

       write(6,*) 'finished with read_input'
       return

1000   write(6,*) 'Error reading Humidity.  Program stopped.'
       write(9,*) 'Error reading Humidity.  Program stopped.'
       close(9)
       stop
1001   write(6,*) 'Error reading dTdz.  Program stopped.'
       write(9,*) 'Error reading dTdz.  Program stopped.'
       close(9)
       stop
1002   write(6,*) 'Error reading z_trop.  Program stopped.'
       write(9,*) 'Error reading z_trop.  Program stopped.'
       close(9)
       stop
1003   write(6,*) 'Error reading H_trop.  Program stopped.'
       write(9,*) 'Error reading H_trop.  Program stopped.'
       close(9)
       stop
1004   write(6,*) 'Error reading dTdz_strat.  Program stopped.'
       write(9,*) 'Error reading dTdz_strat.  Program stopped.'
       close(9)
       stop
1005   write(6,*) 'Error reading T_air(1).  Program stopped.'
       write(9,*) 'Error reading T_air(1).  Program stopped.'
       close(9)
       stop
1015   write(6,*) 'Error reading wind speed, slope, or direction.  Program stopped.'
       write(9,*) 'Error reading wind speed, slope, or direction.  Program stopped.'
       close(9)
       stop
1006   write(6,*) 'Error reading vent diameter.  Program stopped.'
       write(9,*) 'Error reading vent diameter.  Program stopped.'
       close(9)
       stop
1007   write(6,*) 'Error reading exit velocity.  Program stopped.'
       write(9,*) 'Error reading exit velocity.  Program stopped.'
       close(9)
       stop
1008   write(6,*) &
         'Error reading mass fraction added water.  Program stopped.'
       write(9,*) &
         'Error reading mass fraction added water.  Program stopped.'
       close(9)
       stop
1009   write(6,*) 'Error reading magma temperature.  Program stopped.'
       write(9,*) 'Error reading magma temperature.  Program stopped.'
       close(9)
       stop
1010   write(6,*) 'Error reading mass fraction gas.  Program stopped.'
       write(9,*) 'Error reading mass fraction gas.  Program stopped.'
       close(9)
       stop
1011   write(6,*) 'Error reading magma specific heat.  Program stopped.'
       write(9,*) 'Error reading magma specific heat.  Program stopped.'
       close(9)
       stop
1012   write(6,*) 'Error reading magma density.  Program stopped.'
       write(9,*) 'Error reading magma density.  Program stopped.'
       close(9)
       stop
1016   write(6,*) 'Error reading vent elevation. Program stopped.'
       write(9,*) 'Error reading vent elevation. Program stopped.'
       close(9)
       stop
2000   write(6,*) 'Error: specified input file does not exist.  Program stopped'
       write(9,*) 'Error: specified input file does not exist.  Program stopped'
       close(9)
       stop

       end subroutine read_input
