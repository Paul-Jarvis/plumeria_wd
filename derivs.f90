        subroutine derivs(x, yarr, dydx)
                
        use Module2
        use precis_param

        implicit none
        real(kind=ip) x, y, yarr(13), dydx(13)
        real(kind=ip) dMxds, dMyds, dMzds, dQds, h_airnow, &
             h_mixnow, HumidityNow,  &
             Massflux, ma_airnow,  m_anow, m_inow, m_lnow, m_mnow, &
             m_wnow, Mx, My, Mz, mv_airnow, &
             m_vnow, MomentumFlux_now, phinow, pnow, &
             rnow, rho_airnow, rho_mixnow, sigma, T_airnow, &
             T_mixnow, TimeNow, unow, ux, uy, uz, &
             v_perpnow, v_s1, v_s2, v_s3, &             
             w_airnow, xnow, ynow, znow
        real(kind=ip) AirPres, AirTemp, AirHumid, psat, h_a, h_v, &
                       vx, vy, wind, zeta     !functions
        real(kind=ip) R_w
        character(len=1) answer
        real(kind=ip)   dQds2
                
        R_w   = 8.314_ip/0.0180152_ip         !specific gas constant for water (J/kg K)
                
        !GET VARIABLES FROM YARR
        Massflux = yarr(1)
        Mx       = yarr(2)
        My       = yarr(3)
        Mz       = yarr(4)
        pnow     = yarr(6)
        m_mnow   = yarr(7)
        m_anow   = yarr(8)
        m_wnow   = yarr(9)
        TimeNow  = yarr(10)
        xnow     = yarr(11)
        ynow     = yarr(12)
        znow     = yarr(13)

        !CALCULATE ADDITIONAL VARIABLES
        ux       = yarr(2)/yarr(1)
        uy       = yarr(3)/yarr(1)
        uz       = yarr(4)/yarr(1)
        phinow   = datan2(uz,dsqrt(ux**2+uy**2))
        if ((ux.ne.0.).or.(uy.ne.0.)) then
            sigma = datan2(uy,ux)
          else
            sigma = pi/2.0_ip-zeta(znow)
        end if
        unow     = dsqrt(ux**2 + uy**2 + uz**2)
        h_mixnow = (yarr(5) / Massflux) - unow ** 2 / 2.0_ip - g * znow
        !when u-->0, MassFlux-->0, and h_mixnow and unow can get very unstable.
        !the lines below stabilize these values until u increases enough to stabilize
        !the calculations
        if (h_mixnow.lt.0.0_ip) then
           uz=0.001_ip
           unow = sqrt(ux**2+uy**2+uz**2)
           h_mixnow = (yarr(5) / Massflux) - unow ** 2 / 2.0_ip - g * znow
        end if
        MomentumFlux_now = MassFlux*unow
                
        !CALCULATE COMPONENT PROPERTIES
        T_airnow = AirTemp(znow)
        HumidityNow = AirHumid(znow)
        rho_airnow = pnow / (R_air * T_airnow)             !air density, kg/m3
        w_airnow = (kgmole_w/kgmole_air)*(psat(T_airnow)* &
                   HumidityNow/pnow)
        ma_airnow = 1. / (1. + w_airnow)                       !mass fraction dry air in ambient air
        mv_airnow = w_airnow / (1. + w_airnow)                 !mass fraction water vapor in ambient air
        h_airnow = ma_airnow * h_a(T_airnow) + mv_airnow * h_v(T_airnow)

        Call FindT(m_mnow, m_anow, m_wnow, h_mixnow, T_mixnow, pnow, &
                   m_vnow, m_lnow, m_inow) !find temperature
        rho_mixnow = 1 / (m_lnow / rho_w + m_mnow / rho_m + m_inow / &
                      rho_ice + &
                     (m_anow * R_air + m_vnow * R_w) * T_mixnow / pnow)
        rnow = sqrt(yarr(1) / (pi * unow * rho_mixnow)) !present radius

        !CALCULATE COPONENTS OF WIND PARALLEL AND PERPENDICULAR TO THE PLUME
        v_s1  = (vx(znow)*dcos(sigma)+vy(znow)*dsin(sigma))*dcos(phinow)
        v_s2 =  (-vx(znow)*dsin(sigma)+vy(znow)*dcos(sigma)) 
        v_s3  = (vx(znow)*dcos(sigma)+vy(znow)*dsin(sigma))*dsin(phinow)
        v_perpnow = dsqrt(v_s2**2 + v_s3**2)  !max wind component perpendicular to plume

        !CALCULATE NEW VALUES OF DYDX
        If ((idensflag.eq.0).and.(rho_mixnow>rho_airnow)) Then
                dydx(1) = 2. * pi * rnow * &
                              ((alpha   * abs(unow-v_s1))**n_exp + &
                               (gamma_w * v_perpnow )**n_exp) **(1./n_exp) * &
                              sqrt(rho_airnow * rho_mixnow)      !Used in Woods 1993 & later
                !write(6,*) 'In derivs, idensflag=0'
                !write(6,*) 'znow=',x
                !write(6,*) 'rho_airnow=',rho_airnow,', rho_mixnow=',rho_mixnow
                !write(6,*) 'alpha=',alpha,', unow=',unow,', v_s1=',v_s1
                !write(6,*) 'n_exp=',n_exp,', v_perpnow=',v_perpnow
                !write(6,*) 'dydx(1)=',dydx(1)
                !write(6,*) 'continue (y/n)?'
                !read(5,'(a1)') answer
                !if (answer.eq.'n') stop
          Else
                idensflag = 1
                dydx(1) = 2. * rho_airnow * pi * rnow * &
                              ((alpha   * abs(unow-v_s1))**n_exp + &
                               (gamma_w * v_perpnow)**n_exp ) ** (1./n_exp)
        End If

        dQds = dydx(1)
      
        !CALCULATE MOMENT COMPONENTS IN X, Y, Z DIRECTION
        dydx(2) = dQds*vx(znow)                                      !dMx/ds
        dydx(3) = dQds*vy(znow)                                      !dMy/ds
        !dydx(4) = pi*rnow**2*g*(rho_airnow-rho_mixnow)*dsin(phinow) !dMz/ds
        dydx(4) = pi*rnow**2*g*(rho_airnow-rho_mixnow)               !dMz/ds
        dMxds   = dydx(2)
        dMyds   = dydx(3)
        dMzds   = dydx(4)
        dydx(5) = dQds * (h_airnow + wind(znow)**2/2. + g * znow)       !energy gradient
        dydx(6) = -rho_airnow * g * sin(phinow)                         !pressure gradient
        dydx(7) = -(m_mnow / Massflux) * dQds                           !dm_m/ds
        dydx(8) = (ma_airnow - m_anow) * dQds / Massflux                !dm_a/ds
        dydx(9) = (mv_airnow - m_wnow) * dQds / Massflux                !dm_w/ds
        dydx(10) = 1. / unow                                             !dt/ds
        dydx(11) = ux/unow                                               !dx/ds
        dydx(12) = uy/unow                                               !dy/ds
        dydx(13) = uz/unow                                               !dz/ds

        !if (x.gt.74634.0_ip) then
        !   write(6,*) 'Massflux=',Massflux,', unow=',unow,', znow=',znow
        !   write(6,*) 'yarr(2)=',yarr(2),', yarr(3)=',yarr(3),', yarr(4)=',yarr(4)
        !   write(6,*) 'yarr(5)=',yarr(5),', yarr(5)/Massflux=',yarr(5)/MassFlux
        !   write(6,*) 'uz=',uz
        !    write(6,*) 'dydx(1)=',dydx(1),', dydx(1)/Massflux=',dydx(1)/MassFlux
        !    write(6,*) 'dydx(5)=',dydx(5)
        !    write(6,*) 'rho_airnow=',rho_airnow,', rho_mixnow=',rho_mixnow
        !    write(6,*) 'rnow=',rnow,', unow=',unow
        !    write(6,*) 'unow=',unow,', ux=',ux,', uy=',uy,', uz=',uz
        !    write(6,*) 'v_s1=',v_s1,', v_perpnow=',v_perpnow
        !    write(6,*) 'm_lnow=',m_lnow,', m_mnow=',m_mnow,', m_inow=',m_inow
        !    write(6,*) 'm_anow=',m_anow,', m_vnow=',m_vnow,', T_mixnow=',T_mixnow
        !    write(6,*) 'pnow=',pnow
        !    write(6,*) 'h_mixnow=',h_mixnow
        !    write(6,*) 'continue?'
        !    read(5,'(a1)') answer
        !    if (answer.eq.'n') stop
        !end if
        end subroutine derivs
