        program main
        
        !This is the Fortran version of Plumeria.  Developed by Larry G. Mastin, U.S. Geological
        !Survey (lgmastin@usgs.gov).
        
        use Module1
        use Module2
        use precis_param
        
        implicit none
        real(kind=ip) snow, du, hdid_min, umin, vdot, sstep
        real(kind=ip) yscale(13)
        real(kind=ip) h_v, h_a, h_m, h_l, airtemp, AirPres, pathmax, &
                       s_zmax, uzchgmin, uzlast, zmax, &
                       z_nb, r_nb, u_nb, rho_nb, Vdot_nb, mdot_solids
        real(kind=ip) gperm3_l(30000), gperm3_i(30000)  
        real(kind=ip) zeta                      !z function
        real(kind=ip) vx, vy, wind              !functions of z
        character(len=1) answer
        logical          StillRising
        integer          imax, step_nb
        real(kind=ip)    nextheight
        logical          stopflag

        external cashkarpqs

        external findT
        external derivs

        !SET DEFAULT VALUES
        z_nb       = 0.0                         !neutral buoyancy elevation
        u_nb       = 0.0                         !upward velocity at nbl
        rho_nb     = 0.0                         !density at nbl
        r_nb       = 0.0                         !radius at nbl
        Vdot_nb    = 0.0                         !volume flux at nbl
        step_nb    = 1
        p_asl      = 1.013e05                    !pressure at sea level (used only if iatmlayers=1)
        kgmole_air = 8.314    / R_air            !specific gas constant for air
        kgmole_w   = 0.0180152                   !molar weight of water
        R_w        = 8.314    / 0.0180152        !R for water
        iatmlayers = 1                           !number of layers in atmosphere
        imax       = 0                           !index at top of plume
        pathmax    = 2.0_ip                      !value of s/s_zmax at which calculations stop
        uzchgmin   = 1.0e-07                     !minimum change in uz between successive steps
        hdid_min   = 0.000001_ip                    !execution stops when hdid<hdid_min*r(1)
        umin       = 0.01_ip                      !velocity below which the simulation stops
        s_zmax     = 99999._ip                   !value of s at zmax
        zmax       = 0.0_ip                      !maximum plume height
        stopatthetop = .true.                    !=.true if we want to stop at the max. elevation
        StillRising  = .true.                    !.true if the plume has not yet reached its peak height
        nextheight = 200._ip
        stopflag   = .false.

        call read_input                          !read input values
        
        !TROPOSPHERE PROPERTIES
        p(1)   = AirPres(vent_elevation)        !air pressure at the vent
        rho_air(1) = p(1) / (R_air * T_air(1))  !air density at the vent

        !MAGMA PROPERTIES
        gamma = 1.25                             !cp/cv for magmatic gas
        
        !VENT PROPERTIES
        z(1) = 0.                                !z is height above the vent, NOT above sea level
        s(1) = 0.                                !distance along plume axis
        x(1) = 0.
        y(1) = 0.
        Time(1) = 0.                             !time (s) after leaving the vent
        
        !WATER PROPERTIES
        rho_w = 1000.             !density of water
        T_water = 273.15          !Temperature of external water mixed with magma at beginning
        rho_ice = 900.            !ice density
        T_ColdWater = 266.65      !top of temperature range at which liquid water & ice coexist
        T_ice = 258.15            !bottom of temperature range which liquid water & ice coexist
        
        !SET INTEGRATION PARAMETERS
        maxsteps = 30000          !Program will stop if this number of steps is exceeded.
        alpha = 0.09              !entrainment coefficient
        gamma_w = 0.5             !crossflow entrainmnent coefficient
        n_exp   = 1.0             !Devenish exponent to entrainment equation

        !Trap possible errors
        If (u(1).lt.5.) Then 
                write(6,*) &
                   'This model can!t handle such low exit velocities.'
                write(6,*) 'Please enter a higher value.'
                stop
          Else if (dtdz.gt.0.) Then 
                write(6,*) 'Please enter a NEGATIVE thermal lapse rate.'
                write(6,*)  &
                   '(a positive rate implies an intrease in T with z)'
                stop
          else if ((iatmlayers.eq.1).and.(humidity.gt.1.)) Then 
                write(6,*)  &
                   'Please enter a relative humidity between 0 and 100'
                stop
        end if
        
        !SET INITIAL CONDITIONS
        istep = 1
        sstep = 1.
        if (r(1).lt.10.) then
          sstep = r(1)/10.
        end if
        hnext = sstep
        m_gas = (1. - mw) * n_0        !mass fraction magmatic gas
        m_w(1) = (1.-n_0air) * (mw + m_gas)  !mass fraction total water in mixture
        m_m(1) = (1.-n_0air) - m_w(1)
        m_a(1) = n_0air                 !mass fraction air (added to keep equations from blowing up)
        ux(1)  = 0.
        uy(1)  = 0.
        uz(1)  = u(1)
        Time(1) = 0.
        phi(1) = pi/2.                 !initial plume angle from horizontal
        idensflag = 0                  !set this flag to zero at the beginning.  It is set to 1 once we cross
        uzlast    = 0.0_ip

                
        !Calculate new mixture temperature
        hmix = m_m(1) * h_m(T_mag, p(1)) + m_gas * h_v(T_mag) + m_a(1) * h_a(T_mag) + mw * h_l(T_water)
        Call FindT_init(m_m(1), m_a(1), m_w(1), hmix, T_mix(1), p(1), m_v(1), m_l(1), m_i(1)) 
        rho_mix(1) = 1 / (m_l(1) / rho_w + m_m(1) / rho_m + m_i(1) / rho_ice + &
             (m_a(1) * R_air + m_v(1) * R_w) * T_mix(1) / p(1))
        mdot = rho_mix(1) * u(1) * 3.14159 * r(1) ** 2 !mass flux
        mdot_solids = mdot* (1.0-m_w(1))                !mass fraction of solids
        gperm3_l(1) = 1000.*m_l(1)*rho_mix(1)          !g/m3 liquid water
        gperm3_i(1) = 1000.*m_i(1)*rho_mix(1)          !g/m3 ice

        !calculate mixture sound speed
        c_mix = sqrt(gamma*p(1)/rho_mix(1))
        
        !Make sure there's enough momentum to get the plume to the first dz step
        du = -((rho_mix(1)-rho_air(1)/rho_mix(1))*g/u(1) + &
              (2.*alpha*u(1)/r(1))*(rho_air(1)/rho_mix(1)))*sstep
        If ((u(1)+du)<0.) Then
                write(6,*) 'Insufficient momentum to lift the plume.  Program stopped'
                stop
        End If
     
        !WRITE DATE, TIME, AND INPUT VALUES TO STDOUT AND TO THE OUTPUT FILE.
        write(6,1) values(1), values(2), values(3), values(5), values(6)  !write date & time
        write(11,1) values(1), values(2), values(3), values(5), values(6)                  !write date & time
1       format(/,'Output from Plumeria 2, fortran version, run date ',i4,'.',i2.2,'.',i2.2,i4,':',i2.2,//,'Input values:')
        If (iatmlayers.eq.1) Then
         write(6,2)  z_trop_above_vent/1000., H_trop/1000., dTdz*1000., humidity*100., &
                    T_air(1)-273.15, p(1)/101300., windconst, windslope, winddir
         write(11,2) z_trop_above_vent/1000., H_trop/1000., dTdz*1000., humidity*100., &
                    T_air(1)-273.15, p(1)/101300., windconst, windslope, winddir
2        format('        Tropopause elevation above vent (km):',f8.1,/, &
                '                   Height of tropopause (km):',f8.1,/, &
                '                   Thermal lapse rate (C/km):',f8.2,/, &
                '                       Relative humidity, % :',f8.1,/, &
                '                 Air temperature at vent (C):',f8.1,/, &
                '          Air pressure  at vent, atmospheres:',f8.3,/, &
                '               wind speed at sea level (m/s):',f8.3,/, &
                '   change in wind speed with elevation (1/s):',f8.3,/, &
                '                wind direction (deg. E of N):',f8.1)
        Else
         write(6,3) Metfile
         write(11,3) Metfile
3        format('      Using multilayered atmosphere from file ',a80)
      End If
      write(6,4)   r(1)*2., vent_elevation, u(1), T_mag-273.15, n_0, n_0air, Cp_m, &
                            rho_m, rho_mix(1), mw, mdot, mdot_solids, c_mix
      write(11,4)  r(1)*2., vent_elevation, u(1), T_mag-273.15, n_0, n_0air, Cp_m, &
                            rho_m, rho_mix(1), mw, mdot,  mdot_solids, c_mix
4        format('                           Vent diameter (m):',f8.1,/, &
                '                          Vent elevation (m):',f8.1,/, &
                '                      Initial velocity (m/s):',f8.1,/, &
                '                       Magma temperature (C):',f8.0,/, &
                '                         Weight fraction gas:',f8.3,/, &
                '                         Weight fraction air:',f8.3,/, &
                '               Magma specific heat, (J/kg K):',f8.1,/, &
                '                      Magma density, (kg/m3):',f8.1,/, &
                '                     Mixture density (kg/m3):',f8.3,/, &
                '                   Mass fraction water added:',f8.3,/, &
                '                           Mass flux, (kg/s):',e11.3,/, &
                '                 Mass flux of solids, (kg/s):',e11.3,/, &
                '                   Mixture sound speed (m/s):',f8.1)

!        write(6,*) 'crossflow length scale=',mdot*(rho_air(1)-rho_mix(1))/(rho_mix(1)**2*windconst**3)
!        stop
        write(6,5)
        write(11,5)
5       format('**********************************************************************************************************************************************************************************************',/, &
               '   i        s        z        x        y       m_m       m_a       m_w       m_l       m_i      v       u          r   T_mix   T_air   rho_mix   rho_air    Time     p_air     water       ice',/, &
               '       meters   meters   meters   meters                                                      m/s     m/s          m  Kelvin  Kelvin     kg/m3     kg/m3 seconds        Pa      g/m3      g/m3')

        !BEGIN LOOP
        do while (istep<maxsteps)
        !do while ((u(istep).gt.umin).and.(istep<maxsteps))

            yarr(1)=pi*r(istep)**2*u(istep)*rho_mix(istep)      !column mass flux
            yarr(2)=pi*r(istep)**2*u(istep)*ux(istep)*rho_mix(istep)   !column momentum
            yarr(3)=pi*r(istep)**2*u(istep)*uy(istep)*rho_mix(istep)   !column momentum
            yarr(4)=pi*r(istep)**2*u(istep)*uz(istep)*rho_mix(istep)   !column momentum
            yarr(5)=yarr(1)*(u(istep)**2/2+g*z(istep)+hmix)     !column total energy
            yarr(6)=p(istep)                                    !column pressure
            yarr(7)=m_m(istep)                                 !mass fraction magma
            yarr(8)=m_a(istep)                                 !mass fraction air
            yarr(9)=m_w(istep)                                 !mass fraction w. vapor
            yarr(10)=Time(istep)                                !time
            yarr(11)=x(istep)                                   !meters north
            yarr(12)=y(istep)                                   !meters east
            yarr(13)=z(istep)                                   !z
                
            yscale(1) = 100.*yarr(1)!Values by which rkqs normalizes errors to check for accuracy
            yscale(2) = 100.*(yarr(2)+yarr(3)+yarr(4))
            yscale(3) = 100.*(yarr(2)+yarr(3)+yarr(4))
            yscale(4) = 100.*(yarr(2)+yarr(3)+yarr(4))
            yscale(5) = 100.*yarr(5)               !column energy
            yscale(6) = yarr(6)               !p
            yscale(7) = 1.0_ip               !m_m
            yscale(8) = 1.0_ip               !m_a
            yscale(9) = 1.0_ip               !m_w
            yscale(10) = 100._ip              !time
            yscale(11)= 1000.                 !x
            yscale(12)= 1000.                 !y
            yscale(13)= 1000.                 !z
                
            snow = s(istep)
                
            !WRITE VALUES AT EACH TIME STEP
            if ((istep.eq.1).or.(mod(istep,1).eq.0)) then
              write(6,6) istep, s(istep), z(istep), x(istep), y(istep), m_m(istep), m_a(istep), m_w(istep), &
                             m_l(istep), m_i(istep), wind(z(istep)), u(istep), &
                             r(istep), T_mix(istep), AirTemp(z(istep)), rho_mix(istep), rho_air(istep), &
                             Time(istep), p(istep), gperm3_l(istep), gperm3_i(istep)
              write(11,6) istep, s(istep), z(istep), x(istep), y(istep), m_m(istep), m_a(istep), m_w(istep), &
                             m_l(istep), m_i(istep), wind(z(istep)), u(istep), &
                             r(istep), T_mix(istep), AirTemp(z(istep)), rho_mix(istep), rho_air(istep), &
                             Time(istep), p(istep), gperm3_l(istep), gperm3_i(istep)
6             format(i4,4f9.1,5e10.3,f8.2,f9.2,f10.1,2f8.1,2f10.3,f8.1,f10.1,2e10.3)
            end if
 
            !write(6,*) 'Continue?'
            !read(5,'(a1)') answer
            !if (answer.eq.'n') stop

            !if (stopflag) then
            !   stopflag = .false.
            !   !write(6,*) 'y(5)/y(1)=',yarr(5)/yarr(1)
            !   !write(6,*) 'hmix=',hmix,', u=',u(istep),', z=',z(istep),', g=',g
            !   !write(6,*) 'dTdz=',dTdz,', R_air=',R_air
            !   !write(6,*) 'T_air=',T_air(istep),', rho_air=',rho_air(istep)
            !   !write(6,*) 'p(istep)=',p(istep),', AirPres(znow)=',AirPres(z(istep))
            !   write(6,*) 'z=',z(istep),', zeta=',zeta(z(istep))*180./pi,', wind=',wind(z(istep))
            !   write(6,*) 'vx=',vx(z(istep)),', vy=',vy(z(istep))
            !   write(6,*) 'Continue (y/n)?'
            !   read(5,'(a1)') answer
            !   if (answer.eq.'n') stop
            !end if

            !ADJUST STEP SIZE BASED ON RESULTS OF LAST STEP
            if ((snow+hnext).gt.nextheight) then
                  sstep = nextheight-snow
                  nextheight = nextheight + 1000._ip
                  stopflag   = .true.
               else
                 sstep = hnext
            end if

            call derivs(snow,yarr,dydx)
            Call cashkarpqs(yarr, dydx, snow, sstep, yscale, hdid, hnext, derivs)

            !sstep = hnext
            istep = istep + 1
            s(istep) = snow
            uzlast   = uz(istep-1)
                
            !ASSIGN PLUME PROPERTIES
            p(istep)        = yarr(6)                                                !air pressure, Pa, at z
            m_m(istep)      = yarr(7)
            m_a(istep)      = yarr(8)                                                !mass fraction air
            m_w(istep)      = yarr(9)                                                !mass fraction total water
            Time(istep)     = yarr(10)                                                !time
            x(istep)        = yarr(11)                                               !x
            y(istep)        = yarr(12)                                               !y
            z(istep)        = yarr(13)                                               !z
            m_m(istep)      = m_m(1)*mdot/yarr(1)
            !p(istep)        = AirPres(z(istep))
            
            !FIND TEMPERATURE, MASS FRACTIONS OF AQUEOUS PHASES AT THIS STEP BASED ON NEW ENTHALPY
            if (z(istep).gt.zmax)  zmax=z(istep)
            ux(istep)       = yarr(2)/yarr(1)
            uy(istep)       = yarr(3)/yarr(1)
            uz(istep)       = yarr(4)/yarr(1)
            if (isnan(uz(istep))) then
                write(6,*) 'istep=',istep,'z=',z(istep)
                write(6,*) 'yarr(4)=',yarr(4),', yarr(1)=',yarr(1)
                stop
            end if
            u(istep)        = dsqrt(ux(istep)**2+uy(istep)**2+uz(istep)**2)
            hmix            = (yarr(5) / yarr(1)) - u(istep) ** 2 / 2 - g * z(istep) !enthalpy
            Call FindT(m_m(istep), m_a(istep), m_w(istep), hmix, T_mix(istep), p(istep), &
                       m_v(istep), m_l(istep), m_i(istep))
            rho_mix(istep) = 1. / (m_m(istep) / rho_m + m_l(istep) / rho_w + m_i(istep) / &
                       rho_ice + (m_a(istep) * R_air + m_v(istep) * R_w) * T_mix(istep) / p(istep))
            r(istep)  = sqrt(yarr(1) / (pi * u(istep) * rho_mix(istep)))
            
            !ATMOSPHERIC PROPERTIES
            T_air(istep)    = AirTemp(z(istep))
            rho_air(istep)  = p(istep) / (R_air * T_air(istep))
            gperm3_l(istep) = 1000.*m_l(istep)*rho_mix(istep)                     !g/m3 liquid water
            gperm3_i(istep) = 1000.*m_i(istep)*rho_mix(istep)                     !g/m3 ice

            !FIND NEUTRAL BUOYANCY ELEVATION
            if (stillrising.and.((rho_mix(istep-1).lt.rho_air(istep-1)).and. &
                  ((rho_mix(istep)  .ge.rho_air(istep))))) then
                 z_nb    = z(istep)
                 rho_nb  = rho_mix(istep)
                 r_nb    = r(istep)
                 u_nb    = u(istep)
                 Vdot_nb = pi*r_nb**2 * u_nb
                 step_nb = istep
            end if

            !SEE IF THE PLUME IS STILL RISING
            if (StillRising.and.(uz(istep).lt.0.0_ip)) then
                StillRising = .false.
                imax        = istep
                s_zmax      = s(istep)
            end if
                
            if ((stopatthetop).and.(uz(istep).lt.0.0_ip)) then
               write(6,10)
               write(11,10)
               write(6,*) 'stopped because uz<0'
               write(11,*) 'stopped because uz<0'
               exit
            end if
            !if (abs(x(istep)).gt.600._ip) then
            !   write(6,10)
            !   write(11,10)
            !   write(6,*) 'stopped because x>600m'
            !   write(11,*) 'stopped because x>600m'
            !   exit
            !end if
            if ((s(istep)/s_zmax).gt.pathmax) then
               write(6,10)
               write(11,10)
               write(6,*) 'stopped because s/z>',pathmax
               write(11,*) 'stopped because s/z>',pathmax
               exit
            end if
            !if ((iatmlayers.eq.1).and.(windconst.eq.0.0).and.(uz(istep).lt.0.1)) then
            !   write(6,10)
            !   write(11,10)
            !   write(6,*) 'stopped because uz(istep)<0.1'
            !   write(11,*) 'stopped because uz(istep)<0.1'
            !   exit
            !end if
            !if (abs(uz(istep)-uzlast)/hdid.lt.uzchgmin) then
            !   write(6,10)
            !   write(11,10)
            !   write(6,*) 'stopped because (uz-uzlast)/hdid<',uzchgmin
            !   write(11,*) 'stopped because (uz-uzlast)/hdid<',uzchgmin
            !   exit
            !end if
            !if ((abs(uz(istep)).lt.0.1_ip).and.(abs(u(istep)-wind(z(istep))).lt.0.1_ip)) then
            !   write(6,10)
            !   write(11,10)
            !   write(6,*) 'stopped because uz<0.1 or (u-wind)<0.1'
            !   write(11,*) 'stopped because uz<0.1 or (u-wind)<0.1'
            !   exit
            !end if
            !if (hdid.lt.(hdid_min*r(1))) then
            !   write(6,10)
            !   write(11,10)
            !   write(6,*) 'stopped because step size < ',hdid_min,'* vent radius'
            !   write(11,*) 'stopped because step size < ',hdid_min,'* vent radius'
            !   exit
            !end if

        end do 

10          format('*********************************************************************************************************************************************************************************************')

        If (istep.eq.maxsteps) Then
           write(6,10)
           write(11,10)
           write(6,7) maxsteps
           write(11,7) maxsteps
7          format('maximum number of integration steps reached (',i5,').')
        End If

        if (u(istep).le.umin) then
            write(6,10)
            write(11,10)
            write(6,8)
            write(11,8)
8           format('u(istep)<umin')
        end if
        
        !calculate height from Sparks curve
        vdot = mdot*(1.-mw)/rho_m
        write(6,9)  zmax/1000., z_nb/1000., step_nb, r_nb, u_nb, rho_nb, Vdot_nb, mdot, &
                    mdot_solids, 1.67*vdot**0.259, 2.0*vdot**0.241
        write(11,9) zmax/1000., z_nb/1000., step_nb, r_nb, u_nb, rho_nb, Vdot_nb, mdot, &
                    mdot_solids,  1.67*vdot**0.259, 2.0*vdot**0.241
9       format('                                   Maximum height =',f7.3,' km',/, &
               '                          Neutral buoyancy height =',f7.2,' km at step ',i4,/, &
               '                                         r at nbl =',f10.2,' m',/, &
               '                                         u at nbl =',f10.2,' m/s',/, &
               '                                       rho at nbl =',f6.4,' kg/m3',/, &
               '                                      Vdot at nbl =',e12.4,' m3/s',/, &
               '                                             mdot =',e12.4', kg/s',/, &
               '                                      mdot_solids =',e12.4', kg/s',/, &
               '              height calculated from Sparks curve =',f7.3,/, &
               'height calculated from Mastin et al. (2009, eq. 1)=',f7.3,/, &
               'Successful completion')

        nsteps(numruns) = istep
        
        close(9)
        close(11)
        end  program main
        
