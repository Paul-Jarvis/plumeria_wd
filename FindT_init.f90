      subroutine FindT_init(m_m,m_a,m_w,h_mix,Tmix,pnow,m_v,m_l,m_i)

      !subroutine that finds the mixture temperature at the vent given a mixture enthalpy and composition.
      !use Module1
      use Module2
      use precis_param

      implicit none 
      real(kind=8), intent(out) :: Tmix
      real(kind=8)  :: m_m,m_a,m_w,h_mix,pnow,m_v,m_l,m_i                                  !arguments
      real(kind=8)  :: cp_a,cp_l,h_a,h_l,h_m,h_v,psat,tsat                                      !functions
      real(kind=8)  :: Cp_mix,Cpa_avg,Cpwi_avg,Cpwl_avg,Cpwv_avg,H_boiling,H_ColdWater, &       !other real variables
                 H_toboil,hmixnow,m_vfreezing,Tboil,w_freezing

      Cpa_avg = 1060.                 !average value of cp_air between 270 and 1000 K
      Cpwv_avg = 2150.                !approximate average specific heat of water vapor between 100 c and 900 C
      Cpwl_avg = 4190.                !Approximate average water specific heat
      Cpwi_avg = 1850.                !average cp for ice

      Tboil = Tsat(pnow)
      H_toboil = m_m * h_m(Tboil, pnow) + m_a * h_a(Tboil) + m_w * h_v(Tboil)         !Enthalpy at boiling, assumine all water is vapor
      H_boiling = m_m * h_m(Tboil, pnow) + m_a * h_a(Tboil) + m_w * h_l(Tboil)        !Enthalpy at boiling, assuming all water is liquid


      w_freezing = (kgmole_w / kgmole_air) * (psat(273.15) / (pnow - psat(273.15)))
      m_vfreezing = m_a * (w_freezing / (1. + w_freezing))
      If (m_vfreezing.gt.m_w) then
           m_vfreezing = m_w
      End If

      !Enthalpy at freezing, assuming all water is liquid
      H_ColdWater = m_m * h_m(273.15, pnow) + m_a * h_a(273.15) + &
                    m_vfreezing * h_v(273.15) + (m_w - m_vfreezing) * h_l(273.15)      

      print *, h_m(273.15, pnow), h_v(273.15), h_l(273.15)
      If (h_mix.gt.H_toboil) then                  !If we're above the boiling regime
         print *, 'Above boiling regime'
          m_v = m_w
          m_l = 0.
          m_i = 0.
          Cp_mix = m_m * Cp_m + m_a * Cpa_avg + m_v * Cpwv_avg
          !Estimate temperature, based on average specific heats.
          Tmix = Tboil + (h_mix - H_toboil) / Cp_mix
          hmixnow = m_m * h_m(Tmix, pnow) + m_v * h_v(Tmix) + m_a * h_a(Tmix)
          Do while ((abs(hmixnow - h_mix)/h_mix).ge.0.001)       !Iterate on final solution
            Tmix = Tmix + (h_mix - hmixnow) / Cp_mix
            hmixnow = m_m * h_m(Tmix, pnow) + m_v * h_v(Tmix) + m_a * h_a(Tmix)
          end do
       else if (h_mix.gt.H_boiling) then       !If we're within the boiling regime
         print *, 'Within boiling regime'
          Tmix = Tboil
          m_v = m_w * (h_mix - H_boiling) / (H_toboil - H_boiling)
          m_l = m_w - m_v
          m_i = 0.
       else if (h_mix.gt.H_ColdWater) then           !If we're below the boiling regime but above the freezing regime
         print *, 'Below boiling regime and above freezing'
          m_l = m_w
          m_i = 0.
          m_v = 0.
          Cp_mix = m_m * Cp_m + m_a * Cpa_avg + m_w * Cpwl_avg
          Tmix = 273.15 + (h_mix - H_ColdWater) / Cp_mix                       !Take a first stab at temperature
          hmixnow = m_m * h_m(Tmix, pnow) + m_l * h_l(Tmix) + m_a * h_a(Tmix)
          Do while ((abs(hmixnow-h_mix)/h_mix).ge.0.001)       !Iterate on final solution
            Cp_mix = m_m * Cp_m + m_a * cp_a(Tmix) + m_l * cp_l(Tmix)
            Tmix = Tmix + (h_mix - hmixnow) / Cp_mix
            hmixnow = m_m * h_m(Tmix, pnow) + m_v * h_v(Tmix) + m_l * h_l(Tmix) + m_a * h_a(Tmix)
         end do
      else
         print *, 'Does not enter any clause'
      End If

    end subroutine findT_init
