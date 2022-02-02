!##############################################################################
      module precis_param

      !module that states precision of various parameters

      integer, parameter :: ip = selected_real_kind(10,50)

      end module precis_param

!##############################################################################

!##############################################################################
        module Module1
        
        !The following variables are used in this program.
        
        !Real Variables:
        !alpha      = entrainment coefficient
        !c_mix      = sound speed at the vent
        !Cp_m       = specific heat of magma
        !dTdz       = thermal lapse rate in atmosphere (T per meter)
        !dTdz_strat = thermal lapse rate above tropopause (T per meter)
        !gamma_w    = crossflow entrainment coefficient
        !gamma      = cp/cv for magmatic gas
        !H_avail    = magmatic heat available to boil water
        !H_trop     = tropopause thickness, meters
        !HeightPred = height predicted from empirical relation (Sparks et al., 1997, p. 118)
        !hmix       = mixture enthalpy
        !humidity   = relative humdity (fraction of saturation value)
        !kgmole_air = molar weight of air (kg/mole
        !kgmole_w   = molar weight of water (kg/mole)
        !m_gas      = mass fraction magmatic gas
        !mdot       = mass flux, kg/s
        !mw         = mass fraction water added to erupting mixture
        !g          = gravitational acceleration
        !n_0        = gas mass fraction at the vent
        !n_0air     = gas mass fraction at the vent, assuming all gas is air.
        !p_asl      = atm. pressure at sea level (used only if iatmlayers=1)
        !p_trop     = air pressure at tropopause
        !pi         = pi
        !Quality    = fraction of water composed of vapor
        !rho_0      = bulk density of mixture at the vent
        !rho_m      = magma density (kg/m3)
        !rho_w      = density of water (1,000 kg/m3)
        !rho_pumice = pyroclast density
        !rho_trop   = air density at tropopause
        !r_0        = initial vent radius
        !T_asl      = air temperature (K) at sea level (used only if iatmlayers=1)
        !T_water    = ambient water temperature
        !T_mag      = magma temperature (K)
        !T_trop     = temperature at tropopause
        !u_0        = initial velocity at vent
        !windconst  = wind speed (constant)
        !winddir    = wind direction (degrees E of N)
        !windslope  = rate of increase in wind speed with elevation (m/s per m)
        !z_trop_above_vent  = elevation of tropopause (m above vent)
        !zstep      = vertical step of integration
        
        use precis_param

        real(kind=ip) t, mdot, hmix, H_Trop, &
             dTdz_strat, c_mix, dTdz, enthalpy, H_avail, HeightPred, &
             humidity, ml, p_asl, Quality, T_asl, vent_elevation

        !stopatthetop = .true. if we want to stop execution at the maximum elevation
        
        !Parameters used in atmospheric layers
        integer, parameter    :: nvars = 15
        integer               :: iatmlayers
        logical               :: ReadmetFile                       
        logical               :: stopatthetop
        real(kind=ip), dimension(2000) :: TairLayer, Zairlayer, HumidAirLayer, pAirLayer,  &
                                          WinddirLayer, WindspeedLayer
        
        real(kind=ip)         :: z_trop, u_0, T_trop, &
                                 R_w, p_trop, n, m_gas, mw, n_0, n_0air, rho_0, rho_trop, &
                                 T_mag, T_water, z_trop_above_vent, zstep
        real(kind=ip) hnext, yarr(13)
        real(kind=ip) hdid
        real(kind=ip) dydx(13)
        real(kind=ip) yout(13)
        
        !variables used for integration
        !yarr(1)         = column mass flux
        !yarr(2)         = column momentum flux in x direction
        !yarr(3)         = column momentum flux in y direction
        !yarr(4)         = column momentum flux in z direction
        !yarr(5)         = column total energy
        !yarr(6)         = pressure
        !yarr(7)         = m_m
        !yarr(8)         = m_a
        !yarr(9)         = m_w
        !yarr(10)        = time
        !yarr(11)        = x
        !yarr(12)        = y
        !yarr(13)        = z
        
        !Real 1-D variable arrays that vary with elevation:
        !m_a()       = mass fraction dry air in column
        !m_i()         = mass fraction ice in water column
        !m_l()         = mass fraction liquid water in column
        !m_m()         = mass fraction magma in column
        !m_v()         = mass fraction water vapor in column
        !m_w()         = total water in column
        !n()          = gas mass fraction in column
        !p()          = air pressure (Pascals)
        !rho_air()    = air density (kg/m3)
        !r()          = plume radius (m)
        !rho_mix()    = bulk density (kg/m3)
        !T_air()      = air temperature (K)
        !T_mix()      = temperature of column (K)
        !phi()        = angle of plume from horizontal
        !u()          = velocity (m/s)
        !z()          = vertical position (m, upward positive)
        
        real(kind=ip), dimension(30000) :: Cp_mix, m_a, m_i, m_l, m_m, m_w, m_v, p, phi, &
                             r, rho_air, rho_mix, s, T_air, T_mix, Time, u, &
                             ux, uy, uz, x, y, z
        
        !Integer variables:
        !idensflag  = flag set to zero at beginning of each run, then 
                      !set to 1 once we cross into the convective thrust region.
        !istep      = step number in model
        !nsteps     = total number of steps in model
        !ierr       = 0 if upward velocity is positive, 1 if <=0
        !maxsteps   = maximum number of integration steps
        
        integer maxsteps, istep, i, ierr, j, numruns
        integer nsteps(5)
        
        !File containing meteorology data
        character*80  :: Metfile

        !output file
        character*80  ::  OutFileName
        integer       :: values(8)                        !used by "date_and_time" function
        character     :: date*8, zone*5, time2*10         !used by "date_and_time" function

     end module Module1

!******************************************************************************

   module Module2
   
   use precis_param
   !module containing some parameters that are used by otherwise private functions
   real(kind=ip) alpha, Cp_m, gamma, gamma_w, kgmole_air, n_exp, kgmole_w, rho_m, rho_w, &
                 rho_ice, T_coldwater, T_ice, windconst, winddir, windslope
   real(kind=ip), parameter       ::  g=9.81_ip
   real(kind=ip), parameter       ::  pi=3.1415926535897932384626433_ip
   real(kind=ip), parameter       ::  R_air = 286.98

   integer idensflag
   !real(kind=ip) psize(30000,50)                                      !particle sizes
   !real(kind=ip) pmass(50)                                            !particle mass
   
   end module Module2
