      subroutine FindT(m_m,m_a,m_w,h_mix,Tmix,pnow,m_v,m_l,m_i)

      !subroutine that finds the mixture temperature given a mixture enthalpy and composition.
      use module2
      use precis_param

      implicit none
      real(kind=ip) m_m,m_a,m_w,h_mix,Tmix,pnow,m_v,m_l,m_i                !arguments
      real(kind=ip) Cp_mix,Cpa_avg,Cpwi_avg,Cpwl_avg,Cpwv_avg, &
           H_ColdWater,H_freezing, &
           H_sat,hdif,Hmixhi,hmixlast,Hmixlo,hmixnow, &
           m_vColdWater,m_vfreezing, &
           TmixHi,TmixLst,TmixLo,T_sat,w_ColdWater,w_freezing, &
           w_s,x_a,x_w
      real(kind=ip) h_a,h_i,h_l,h_m,h_v,psat,tsat                           !functions
      integer i
      character*1  answer

      !SET DEFAULT VALUES
      kgmole_air = 8.314 / 286.98         !molar weight of air
      kgmole_w = 0.0180152                            !molar weight of water

      Cpa_avg = 1060.                 !average value of cp_air between 270 and 1000 K
      Cpwv_avg = 2150.                !approximate average specific heat of water vapor between 100 c and 900 C
      Cpwl_avg = 4190.                !Approximate average water specific heat
      Cpwi_avg = 1850.                !average cp for ice
 
      x_w = (m_w / kgmole_w) / ((m_w / kgmole_w) + (m_a / kgmole_air))  !mole fraction water vapor at saturation
      x_a = 1. - x_w                                                    !mole fraction air at saturation
      T_sat = Tsat(x_w * pnow)                                          !Saturation temperature
      H_sat = m_m*h_m(T_sat, pnow) + m_a*h_a(T_sat) + m_w*h_v(T_sat)    !Enthalpy at saturation, assumine all water is vapor

      !FIND H_freezing
      if (pnow>psat(T_ice)) then                      !if the ambient pressure exceeds the boiling pressure at T_ice
         w_freezing = (kgmole_w/kgmole_air) * &
                      (psat(T_ice)/(pnow - psat(T_ice)))
         m_vfreezing = m_a * w_freezing
         if (m_w > m_vfreezing) then                   !if there!s enough water vapor to saturate the plume
             H_freezing = m_m * h_m(T_ice, pnow) + &
                          m_a * h_a(T_ice) + m_vfreezing*h_v(T_ice) + &
                     (m_w - m_vfreezing) * h_i(T_ice)
           else                                          !if air is non water-saturated at T_ice
             m_vfreezing = m_w
             H_freezing = m_m * h_m(T_ice, pnow) + m_a * h_a(T_ice) + &
                          m_vfreezing * h_v(T_ice)
         end if
       else                                
         !if pnow is less than the boiling pressure at T_ice
         m_vfreezing = m_w
         H_freezing = m_m * h_m(T_ice, pnow) + m_a * h_a(T_ice) + &
                      m_vfreezing * h_v(T_ice)
      end if
     
      !FIND H_coldwater
      if (pnow>psat(T_ColdWater)) then   
        !if pnow exceeds the boiling pressure at T_ColdWater (i.e. pnow > psat)
         w_ColdWater = (kgmole_w / kgmole_air) * (psat(T_ColdWater) / &
                       (pnow - psat(T_ColdWater)))
         m_vColdWater = m_a * w_ColdWater
         if (m_w > m_vColdWater) then          !if there!s enough water vapor to saturate the plume 
            !Enthalpy at top of mixed water-ice temperature range
             H_ColdWater = m_m * h_m(T_ColdWater, pnow) + &
                           m_a * h_a(T_ColdWater) + &
                           m_vColdWater * h_v(T_ColdWater) + &
                     (m_w - m_vColdWater) * h_l(T_ColdWater)
           else           
             !if there!s not enough water vapor to saturate the plume
             m_vColdWater = m_w
             H_ColdWater = m_m * h_m(T_ColdWater, pnow) + &
                           m_a * h_a(T_ColdWater) + &
                           m_vColdWater * h_v(T_ColdWater)
          end if
       else                      
         !if we're above the boiling point at this pressure
         m_vColdWater = m_w
         H_ColdWater = m_m * h_m(T_ColdWater, pnow) +  &
                       m_a * h_a(T_ColdWater) + &
                       m_vColdWater * h_v(T_ColdWater)
      end if

      !The following criteria are met only as a result of inaccuracies 
      !in calculating T_sat, when it's close
      !to 273.15 K.  if T_sat<273.15 K and h_mix > H_freezing 
      !and h_mix<h_sat, the program tries to calculate
      !the enthalpy of liquid water at T<273.15, which causes it to blow up.
      !if (h_mix < H_sat) And (h_mix > H_freezing) And (T_sat < 273.15) then
      !    T_sat = 273.15
      !end if

      if ((h_mix > H_sat).and.(h_mix > H_ColdWater)) then  
        !if we're above water saturation, and above freezing
         m_v = m_w
         m_l = 0.
         m_i = 0.
         Cp_mix = m_m * Cp_m + m_a * Cpa_avg + m_v * Cpwv_avg
         !Estimate temperature, based on average specific heats.
         if (H_sat > H_ColdWater) then
                 Tmix = T_sat + (h_mix - H_sat) / Cp_mix
           else
                 Tmix = T_ColdWater + (h_mix - H_ColdWater) / Cp_mix
         end if
         hmixnow = m_m * h_m(Tmix, pnow) + &
                   m_v * h_v(Tmix) + m_a * h_a(Tmix)
         do while (Abs(hmixnow - h_mix) / h_mix > 0.001)      
                 !Iterate on final solution
                 Tmix = Tmix + (h_mix - hmixnow) / Cp_mix
                 hmixnow = m_m * h_m(Tmix, pnow) + &
                           m_v * h_v(Tmix) + &
                           m_a * h_a(Tmix)
         end do
        else if ((h_mix > H_ColdWater).and.(T_sat > T_ColdWater)) then
          !if we're within the saturated regime but above the freezing regime
          i = 1
          m_i = 0.
          !Take a first stab at temperature
          Tmix = T_ColdWater + (T_sat - T_ColdWater) * &
                 (h_mix - H_ColdWater) / (H_sat - H_ColdWater)
          w_s = (kgmole_w / kgmole_air) *  &
                (psat(Tmix) / (pnow - psat(Tmix)))
          if (m_w >= m_a * w_s) then        !if we're water saturated
                   m_v = m_a * w_s
                   m_l = m_w - m_v
             else
                   m_v = m_w
                   m_l = 0.
          end if
          Hmixhi = H_sat
          TmixHi = T_sat
          Hmixlo = H_ColdWater
          TmixLo = T_ColdWater
          if (Hmixhi < Hmixlo) then
                  write(6,*) 'Problem with enthalpy calculations.'
                  write(6,*) 'Calculations stopped'
                  stop
          end if
          hmixnow = m_m * h_m(Tmix, pnow) + &
                    m_v * h_v(Tmix) + &
                    m_l * h_l(Tmix) + &
                    m_a * h_a(Tmix)
          hmixlast = 0.0
          Do while (((Abs(hmixnow - h_mix) / h_mix > 0.001)).and. &
                    ((Abs(hmixnow - hmixlast)/hmixnow)>1.0e-07)) !Iterate on final solution
                  hmixlast = hmixnow
                  if (hmixnow > h_mix) then
                          Hmixhi = hmixnow
                          TmixHi = Tmix
                    else
                          Hmixlo = hmixnow
                          TmixLo = Tmix
                  end if
                  if (Hmixhi < Hmixlo) then
                        write(6,*) 'Problem with enthalpy '
                        write(6,*) 'calculations. Calculations stopped'
                        exit
                  end if
                  Tmix = TmixLo + (TmixHi - TmixLo) * (h_mix - Hmixlo) &
                         / (Hmixhi - Hmixlo)
                  if (Tmix < T_ColdWater) then
                        write(6,*) 'Problem with enthalpy '
                        write(6,*) 'calculations. Calculations stopped'
                        exit
                  end if
                  w_s = (kgmole_w / kgmole_air) * (psat(Tmix) / &
                        (pnow - psat(Tmix)))        !recalculate m_v
                  if (m_w >= m_a * w_s) then
                          m_v = m_a * w_s
                          m_l = m_w - m_v
                    else
                          m_v = m_w
                          m_l = 0.
                  end if
                  hmixnow = m_m * h_m(Tmix, pnow) + &
                            m_v * h_v(Tmix) + &
                            m_l * h_l(Tmix) + &
                            m_a * h_a(Tmix)
                  i = i + 1
                  !if i > 100 then
                  !    Stop
                  !end if
          end do
       else if (h_mix > H_freezing) then      !if we're at T_ice < tnow < T_ColdWater
            i = 1
            !Take a first stab at temperature
            Tmix = T_ice + (T_ColdWater-T_ice) * (h_mix-H_freezing) / &
                   (H_ColdWater - H_freezing)
            w_s = (kgmole_w/kgmole_air) * (psat(Tmix)/(pnow-psat(Tmix)))
            if (w_s.lt.0.0_ip) w_s=m_w/m_a    !w_s<0 only very high, perhaps at z>40 km
            if (m_w >= m_a * w_s) then        !if we're water saturated
                     m_v = m_a * w_s
                     m_l = (m_w-m_v)*(Tmix-T_ice)/(T_ColdWater - T_ice)
                     if (m_l.gt.1.0_ip) then
                        write(6,*) 'in FindT.  m_l>1.  m_l=',m_l
                        write(6,*) 'pnow=',pnow,', psat(Tmix)=',psat(Tmix)
                        write(6,*) 'm_a=',m_a,', w_s=',w_s
                        write(6,*) 'm_w=',m_w,', m_v=',m_v
                        write(6,*) 'Tmix=',Tmix,', T_ice=',T_ice
                        write(6,*) 'T_Coldwater=',T_ColdWater
                        stop
                     end if
                     m_i = m_w - m_v - m_l
               else
                     m_v = m_w
                     m_l = 0.
                     m_i = 0.
            end if
            Hmixhi = H_ColdWater
            TmixHi = T_ColdWater
            Hmixlo = H_freezing
            TmixLo = T_ice
            if (Hmixhi < Hmixlo) then
                    write(6,*) 'Problem with enthalpy calculations.'
                    write(6,*) 'Calculations stopped'
                    stop
            end if
            hmixnow = m_m * h_m(Tmix, pnow) + &
                      m_v * h_v(Tmix) + &
                      m_l * h_l(Tmix) + &
                      m_i * h_i(Tmix) + &
                      m_a * h_a(Tmix)
            hdif = abs(hmixnow-h_mix)/h_mix
            do while ((Abs(hmixnow - h_mix) / h_mix > 0.001).and.&
                      ((TmixHi - TmixLo) > 0.25).and. &
                      abs(TmixLst - Tmix) > 0.1)                       !Iterate on final solution
               if (hmixnow > h_mix) then
                       Hmixhi = hmixnow
                       TmixHi = Tmix
                 else
                       Hmixlo = hmixnow
                       TmixLo = Tmix
               end if
               if (Hmixhi < Hmixlo) then
                       write(6,*) 'Problem with enthalpy calculations.'
                       write(6,*) ' Calculations stopped'
                       stop
               end if
               TmixLst = Tmix
               Tmix = TmixLo + (TmixHi-TmixLo)*(h_mix-Hmixlo)/&
                               (Hmixhi-Hmixlo)
               if (Tmix < T_ice) then
                       write(6,*) 'Problem with enthalpy calculations.'
                       write(6,*) ' Calculations stopped'
                       exit
               end if
               w_s = (kgmole_w/kgmole_air)*(psat(Tmix)/ &
                                           (pnow-psat(Tmix)))        !recalculate m_v
               if (m_w >= m_a * w_s) then
                       m_v = m_a * w_s
                       m_l = (m_w - m_v) * (Tmix-T_ice) / &
                                           (T_ColdWater-T_ice)
                       m_i = m_w - m_v - m_l
                 else
                       m_v = m_w
                       m_l = 0.
                       m_i = 0.
               end if
               hmixnow = m_m * h_m(Tmix, pnow) + &
                         m_v * h_v(Tmix) + &
                         m_l * h_l(Tmix) + &
                         m_i * h_i(Tmix) + m_a * h_a(Tmix)
               i = i + 1
               hdif = abs(hmixnow-h_mix)/h_mix
            end do
          else                                    !if we're below T_ice
            i = 1
            m_l = 0.
            Hmixhi = H_freezing
            Hmixlo = m_m * h_m(100., pnow) + &
                     m_w * h_i(100.) + &
                     m_a * h_a(100.)
            TmixHi = T_ice
            TmixLo = 100.
            Tmix = 100. + (TmixHi - TmixLo) * (h_mix - Hmixlo) / &
                                               (Hmixhi - Hmixlo) !Take a first stab at temperature
            w_s = (kgmole_w / kgmole_air) * (psat(Tmix) / &
                                            (pnow - psat(Tmix)))
            if (m_w >= m_a * w_s) then
                     m_v = m_a * w_s
                     m_i = m_w - m_v
               else
                     m_v = m_w
                     m_i = 0.
            end if
            hmixnow = m_m * h_m(Tmix, pnow) + &
                      m_v * h_v(Tmix) + &
                      m_i * h_i(Tmix) + &
                      m_a * h_a(Tmix)
            Do while ((abs(hmixnow - h_mix) / h_mix > 0.001).and. &
                      ((TmixHi - TmixLo) > 0.25)) !Iterate on final solution
              if (hmixnow > h_mix) then
                      TmixHi = Tmix
                      Hmixhi = hmixnow
                else
                      TmixLo = Tmix
                      Hmixlo = hmixnow
              end if
              Tmix = TmixLo + (TmixHi - TmixLo) * (h_mix - Hmixlo) / &
                                                  (Hmixhi - Hmixlo)
              w_s = (kgmole_w / kgmole_air)*(psat(Tmix) / &
                                            (pnow - psat(Tmix)))        !recalculate m_v
              if (m_w >= m_a * w_s) then
                      m_v = m_a * w_s
                      m_i = m_w - m_v
                else
                      m_v = m_w
                      m_i = 0.
              end if
              hmixnow = m_m * h_m(Tmix, pnow) + &
                        m_v * h_v(Tmix) + &
                        m_i * h_i(Tmix) + &
                        m_a * h_a(Tmix)
              if (i > 20) stop
            end do
        end if
                            
        end subroutine FindT
