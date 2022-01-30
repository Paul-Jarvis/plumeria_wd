!******************************************************************************         
        function psat(t)
                
                !function that returns the partial pressure of water (in Pascals) at saturation, given
                !a temperature (t) in Kelvin.  Taken from  Ghiorso.
                use Module2
                use precis_param
                use precis_param
                
                implicit none
                real(kind=ip) t
                real(kind=ip) A(8)
                real(kind=ip) psat, wsq, v, ff, e0, MassFrac, w, Zhaar
                integer i
                
                data A/-7.8889166,2.5514255,-6.716169,33.239495,-105.38479,174.35319,-148.39348,48.631602/
                
                e0 = 611. !partial pressure of water vapor at 273.15 K
                
                If (t.lt.T_ice) Then
                        psat = exp(log(e0) + 22.49 - 6142 / t)
                        return
                  ElseIf (t.lt.T_ColdWater) Then 
                        !in the temperature range (T_ice<t<=T_coldwater) where liquid & ice coexist,
                        !the partial pressure is weighted between the partial pressure of ice and of liquid.
                        MassFrac = (t - T_ice) / (T_ColdWater - T_ice)
                        psat = (1. - MassFrac) * exp(log(e0) + 22.49 - 6142 / t) + MassFrac * &
                               (100000. * exp(6.3573118 - 8858.843 / t + 607.56335 / (t ** 0.6)))
                        return
          ElseIf (t.le.314.) Then 
                        psat = 100000. * exp(6.3573118 - 8858.843 / t + 607.56335 / (t ** 0.6))
                        return
                Else
                        v = t / 647.25
                        w = abs(1. - v * 0.9997)
                        wsq = sqrt(w)
                        ff = 0.
                        do i = 1,8
                                Zhaar = i
                                ff = ff + A(i) * w ** ((Zhaar + 1.) / 2.)
                        end do
                        psat = 100000. * 1.001 * 220.93 * exp(ff / v)
                End If
        End Function psat

!****************************************************************************** 
        function AirTemp(znow)
        
        !Function that calculates air temperature in a layered atmosphere, 
        !given the elevation znow above the vent
        use Module1
        use Module2
        use precis_param
                
        implicit none
        integer ilayer 
        real(kind=ip) znow, AirTemp
        !character(len=1)  :: answer
                
!        write(6,*) 'in AirTemp: iatmlayers=',iatmlayers

        !For a 1-layer atmosphere
        If (iatmlayers==1) Then
           If ((znow+vent_elevation).lt.z_trop) Then            !If we!re still in the troposphere
              AirTemp = T_air(1) + dTdz * znow
           ElseIf ((znow+vent_elevation).lt.(z_trop+H_Trop)) Then  !in 10-km isothermal layer at tropopause
              AirTemp = T_air(1)+dTdz*(z_trop-vent_elevation)
           Else !in stratosphere
              AirTemp = T_air(1)+dTdz*(z_trop-vent_elevation)+ &
                       dTdz_strat*((znow+vent_elevation)-(z_trop+H_Trop))
           End If
           !For a multilayer atmosphere
         else
           If (((znow+vent_elevation).gt.ZairLayer(1)).and. &
               ((znow+vent_elevation).lt.ZairLayer(iatmlayers))) Then
               do ilayer = 2,iatmlayers
                  If ((znow + vent_elevation).lt.ZairLayer(ilayer)) Then
                     AirTemp = TairLayer(ilayer-1)+((znow+vent_elevation) - &
                             ZairLayer(ilayer-1))* &
                             (TairLayer(ilayer)-TairLayer(ilayer-1)) / &
                             (ZairLayer(ilayer) - ZairLayer(ilayer - 1))
                      exit
                   End If
                end do
           ElseIf ((znow + vent_elevation)<=ZairLayer(1)) Then  !If we!re below the lowest sounding
              AirTemp = TairLayer(1) + ((znow + vent_elevation) - ZairLayer(1)) * &
                        (TairLayer(2) - TairLayer(1)) / &
                        (ZairLayer(2) - ZairLayer(1))
           Else                                                            !If we!re above the highest sounding
              AirTemp = TairLayer(iatmlayers) + ((znow + vent_elevation) - &
                        ZairLayer(iatmlayers - 1)) * &
                        (TairLayer(iatmlayers) - TairLayer(iatmlayers - 1)) / &
                        (ZairLayer(iatmlayers) - ZairLayer(iatmlayers - 1))
           End If
        End If
        return
        End Function AirTemp
        
!****************************************************************************** 

        function AirPres(VentElevation)
        !This function is called only to find air pressure at the vent
        use Module1
        use Module2
        use precis_param
                
        implicit none
        real(kind=ip) VentElevation
        integer ilayer
        real(kind=ip) airpres, Tnow, p0, LapseRate, T0, z0
        real(kind=ip) AirTemp
                
        Tnow = AirTemp(0.0_ip) !gives air temperature at the vent elevation
      
        !FOR A SIMPLE ATMOSPHERE
        if (iatmlayers.eq.1) then
           p0 = p_asl
           T0 = T_asl        
           z0 = 0.
           LapseRate = dTdz
        !FOR A MULTI-LAYER ATMMOSPHERE
         else
            If ((VentElevation>=ZairLayer(1)).and. &
                (VentElevation<ZairLayer(iatmlayers))) Then
               do ilayer = 2,iatmlayers
                  If (VentElevation<ZairLayer(ilayer)) Then
                     LapseRate = (TairLayer(ilayer) - TairLayer(ilayer - 1)) / &
                                 (ZairLayer(ilayer) - ZairLayer(ilayer - 1))
                     p0 = pAirLayer(ilayer - 1)
                     z0 = ZairLayer(ilayer - 1)
                     T0 = AirTemp(ZairLayer(ilayer - 1) - VentElevation)
                     exit
                  End If
               end do
              ElseIf (VentElevation<ZairLayer(1)) Then  !If we!re below the lowest sounding
                    LapseRate = (TairLayer(2) - TairLayer(1)) / (ZairLayer(2) - ZairLayer(1))
                    p0 = pAirLayer(1)
                    z0 = ZairLayer(1)
                    T0 = AirTemp(ZairLayer(1) - VentElevation)
              Else !If we!re above the highest sounding
                    LapseRate = (TairLayer(iatmlayers) - TairLayer(iatmlayers - 1)) / &
                                (ZairLayer(iatmlayers) - ZairLayer(iatmlayers - 1))
                    p0 = pAirLayer(iatmlayers)
                    z0 = ZairLayer(iatmlayers)
                    T0 = AirTemp(ZairLayer(iatmlayers) - VentElevation)
            End If
        end if
                
        if (abs(LapseRate).gt.1.e-08) Then
               AirPres = p0 * (T0 / Tnow) ** (g / (R_air * LapseRate))
          Else
               AirPres = p0 * exp(-g * (VentElevation - z0) / (R_air * Tnow))
        End If
        End Function AirPres
        
!****************************************************************************** 
        
        function AirHumid(znow)
        !Function that calculates humidity in a layered atmosphere, 
        !given the elevation z above the vent
        use Module1
        use Module2
        use precis_param
                
        implicit none
        integer ilayer
        real(kind=ip) znow
        real(kind=ip) AirHumid
                
        !For a 1-layer atmosphere
        If (iatmlayers.eq.1) Then
                AirHumid = humidity
         !For a multilayer atmosphere
          Else
             If (((znow + vent_elevation)>ZairLayer(1)).and. &
                 ((znow + vent_elevation)<ZairLayer(iatmlayers))) Then
                  do ilayer = 2,iatmlayers
                    If ((znow + vent_elevation)<ZairLayer(ilayer)) Then
                       AirHumid = HumidAirLayer(ilayer - 1) + &
                             ((znow + vent_elevation) - ZairLayer(ilayer - 1)) * &
                             (HumidAirLayer(ilayer) - HumidAirLayer(ilayer - 1)) / &
                             (ZairLayer(ilayer) - ZairLayer(ilayer - 1))
                       exit
                    End If
                  end do
             ElseIf ((znow + vent_elevation)<=ZairLayer(1)) Then  !If we!re below the lowest sounding
               AirHumid = HumidAirLayer(1) + ((znow + vent_elevation) - ZairLayer(1)) * &
                           (HumidAirLayer(2) - HumidAirLayer(1)) / (ZairLayer(2) - ZairLayer(1))
             Else                                            !If we!re above the highest sounding
               AirHumid = HumidAirLayer(iatmlayers) + &
                     ((znow + vent_elevation) - ZairLayer(iatmlayers - 1)) * &
                     (HumidAirLayer(iatmlayers) - HumidAirLayer(iatmlayers - 1)) / &
                     (ZairLayer(iatmlayers) - ZairLayer(iatmlayers - 1))
             End If
        End If
        End Function AirHumid

!******************************************************************************

       function zeta(znow)
       !Function that calculates the wind direction parallel to the plume
       use Module2
       use Module1
       use precis_param

       implicit none
       integer        :: ilayer
       real(kind=ip)  :: speednow, xwind1, xwind2, xwindnow, ywind1, ywind2, ywindnow, zeta, znow
       real(kind=ip)  :: DEG2RAD, RAD2DEG

       DEG2RAD = pi/180.0_ip
       RAD2DEG = 180.0_ip/pi
       ilayer = 0
      
       if (iatmlayers.eq.1) then
           !For a single-layered atmosphere
           !(we add 180 degrees because to convert it from the direction FROM WHICH the wind is blowing
           !the the direction TO WHICH the wind is blowing.  Zeta is measured in degrees CW from north.
           zeta = (winddir+180.0_ip)*DEG2RAD
         else
           !For a multi-layered atmosphere
           If (((znow + vent_elevation)>ZairLayer(1)).and. &
              ((znow + vent_elevation)<ZairLayer(iatmlayers))) Then
               do ilayer = 2,iatmlayers
                 If ((znow + vent_elevation)<ZairLayer(ilayer)) Then
                    xwind1 = WindspeedLayer(ilayer-1)*sin(DEG2RAD*(WinddirLayer(ilayer-1)+180.))
                    xwind2 = WindspeedLayer(ilayer)  *sin(DEG2RAD*(WinddirLayer(ilayer)+180.))
                    ywind1 = WindspeedLayer(ilayer-1)*cos(DEG2RAD*(WinddirLayer(ilayer-1)+180.))
                    ywind2 = WindspeedLayer(ilayer)  *cos(DEG2RAD*(WinddirLayer(ilayer)+180.))
                    xwindnow = xwind1 + ((znow + vent_elevation) - ZairLayer(ilayer-1)) * &
                               (xwind2 - xwind1) / (ZairLayer(ilayer) - ZairLayer(ilayer-1))
                    ywindnow = ywind1 + ((znow + vent_elevation) - ZairLayer(ilayer-1)) * &
                               (ywind2 - ywind1) / (ZairLayer(ilayer) - ZairLayer(ilayer-1))
                    speednow = sqrt(xwindnow**2 + ywindnow**2)
                    zeta = atan2(xwindnow,ywindnow)
                    exit
                 End If
               end do
             ElseIf ((znow + vent_elevation)<=ZairLayer(1)) Then  !If we!re below the lowest sounding
               zeta = DEG2RAD*(Winddirlayer(1)+180.0_ip)
             Else                                            !If we!re above the highest sounding
               zeta = DEG2RAD*(WinddirLayer(iatmlayers) + 180.0_ip)
             End If
       end if

       end function zeta

!******************************************************************************

       function wind(znow)
       !Function that calculates the wind speed at elevation znow
       use Module2
       use Module1
       use precis_param

       implicit none
       integer        :: ilayer
       real(kind=ip)  :: wind, xwind1, xwind2, xwindnow, ywind1, ywind2, ywindnow, znow
       real(kind=ip)  :: DEG2RAD
       character(len=1) :: answer

       ilayer = 0
       DEG2RAD = 3.14159/180.

       if (iatmlayers.eq.1) then
           wind = windconst + windslope*znow
         else
           !For a multi-layered atmosphere
           If (((znow + vent_elevation)>ZairLayer(1)).and. &
              ((znow + vent_elevation)<ZairLayer(iatmlayers))) Then
               do ilayer = 2,iatmlayers
                 If ((znow + vent_elevation)<ZairLayer(ilayer)) Then
                    xwind1 = WindspeedLayer(ilayer-1)*sin(DEG2RAD*(WinddirLayer(ilayer-1)+180.))
                    xwind2 = WindspeedLayer(ilayer)  *sin(DEG2RAD*(WinddirLayer(ilayer)+180.))
                    ywind1 = WindspeedLayer(ilayer-1)*cos(DEG2RAD*(WinddirLayer(ilayer-1)+180.))
                    ywind2 = WindspeedLayer(ilayer)  *cos(DEG2RAD*(WinddirLayer(ilayer)+180.))
                    xwindnow = xwind1 + ((znow + vent_elevation) - ZairLayer(ilayer-1)) * &
                               (xwind2 - xwind1) / (ZairLayer(ilayer) - ZairLayer(ilayer-1))
                    ywindnow = ywind1 + ((znow + vent_elevation) - ZairLayer(ilayer-1)) * &
                               (ywind2 - ywind1) / (ZairLayer(ilayer) - ZairLayer(ilayer-1))
                    wind = sqrt(xwindnow**2 + ywindnow**2)
                    exit
                 End If
               end do
             ElseIf ((znow + vent_elevation)<=ZairLayer(1)) Then  !If we!re below the lowest sounding
               wind = WindspeedLayer(1)
             Else                                            !If we!re above the highest sounding
               wind = WindspeedLayer(iatmlayers)
             End If
       end if

       end function wind

!******************************************************************************

       function vx(znow)
       !Function that calculates the x component of the wind vector

       use Module1
       use Module2       
       use precis_param

       implicit none
       real(kind=ip)  :: vx, wind, zeta, znow

       vx = wind(znow)*dcos(pi/2.0_ip-zeta(znow))

       end function vx
       
       !******************************************************************************

       function vy(znow)
       !Function that calculates the y component of the wind vector
       
       use Module1
       use Module2       
       use precis_param

       implicit none
       real(kind=ip)  :: vy, wind, zeta, znow

       vy = wind(znow)*dsin(pi/2.0_ip-zeta(znow))

       end function vy
