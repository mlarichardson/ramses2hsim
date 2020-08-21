MODULE ramses_info
  !-----------------------------------------------------------------------
  ! Module to read and contain the information contained in
  ! info_xxxxx.txt and info_rt_xxxxx.txt (if it exists).
  ! J. Rosdahl, Jan 2011, modified from code by R. Teyssier.
  !-----------------------------------------------------------------------

  implicit none

   character(LEN=5)::nchar
 ! All units are cgs
  integer::ncpu,ndim,npart
  integer::nlevelmax,levelmin
  integer::ngridmax,nstep_coarse
  real::t,aexp,scale,h0,tot_sfr=0.d0,tot_sfr2=0.d0
  real(kind=8) gamma
  real::omega_m,omega_l,omega_k,omega_b
  real(kind=8)::unit_l,unit_d,unit_t,unit_m
  real(kind=8)::unit_v, unit_P, unit_T2, unit_nH
  real::t_myr
  real,parameter::pi=3.14159265359e0
  real,parameter::kb= 1.3806200e-16
  real,parameter::mH= 1.6600000e-24
  real,parameter::c = 2.9979250e+10
  real,parameter::e_lya = 1.634e-11
  real,parameter::yr2s = 3.155693d7
  real,parameter::Mpc=3.08e24
  real,parameter::rhoc=1.88d-29
  real(kind=8)::SFRtoHA=1.26d41
  real(kind=8) :: hubble
  real(kind=8),parameter::eV_to_erg=1.6022d-12  ! eV to erg conv. constant
  real(kind=8),parameter,dimension(3)::ion_egy  &
       = (/13.60d0, 24.59d0, 54.42d0/)*eV_to_erg

  character(LEN=80)::ordering
  real(kind=8),dimension(:),allocatable::bound_key

  ! SF Parameters
! Nut 10pc
  real::n_star=400.e0,eps_star=-0.01,del_star=1.0
  real::t_star=0.316228,t0,d0,nISM0,nCOM

! Nut 50pc
!  real::n_star=25.e0,eps_star=-0.01,del_star=1.0
!  real::t_star=1.264912,t0,d0,nISM0,nCOM

  ! RT variables:
  integer:: nHvar=8, nRTvar=0, nIons=3, nGroups=0, iIons=6 !7 if metals
  integer::iPrad,iPT
  real:: X_frac=0.76, Y_frac=0.24
  real::unit_Np=1., unit_Fp=1., rt_c_frac=1.
  real(kind=8),parameter::c_cgs =2.9979250d+10
  real(kind=8)::rt_c_cgs = 2.9979250d+10
  real(kind=8)::rt_c=0d0,rt_kIR=1d1

  ! Doppler variables
  real, parameter :: max_deltav=3.e0

  !real::n_star=0.1, T2_star=1.d4, g_star=1.6        ! Polytropic variables
  !real::n_star=10d0, T2_star=100d0, g_star=2d0        ! Polytropic variables
  real::T2_star=0d0, g_star=2d0        ! Polytropic variables

  ! Emissivity variables
  character(len=20)::emType
  integer::emVar
  logical::volMult=.false., noPol = .false., do_inverse=.false.
  logical::div_nH2=.false. ! Divide result by number density squared
  logical::calc_TK=.false. ! Result includes temperature calculation
  logical::RT_OTSA=.true.
  real::defUnit=1.,defSubtr=0.

  logical::nml_ok,cooling,haardt_madau,metal,isothermal,bondi,neq_chem   &
       ,sink_agn,sinkprops,drag,delayed_cooling,sf_hopkins               &
       ,self_shielding,smbh,agn,sn2_real_delay, log_mfb
  real(kind=8)::m_star=-1                                                &
       & ,jeans_ncells,eta_sn,yield,rbubble,f_ek,ndebris,f_w,mass_gmc    &
       & ,J21,a_spec,z_ave,z_reion,n_sink,ind_rsink,rsink_max            &
       & ,units_density,units_time,units_length,ir_feedback,ir_eff       &
       & ,merge_stars,larson_lifetime,flux_accretion,t2thres_sf          &
       & ,sn_delay_myr,sn_dt2_min,mseed,x_floor,sigmav_max,r_gal         &
       & ,T2max_AGN,t_delay,t_delcool,sig_delcool,msink_max,M_SNII
  integer::mechanical_feedback
  character(LEN=80)::star_maker, star_imf
  namelist/physics_params/cooling,haardt_madau,metal,isothermal,bondi     &
       & ,m_star,t_star,n_star,T2_star,g_star,del_star,eps_star           &
       & ,jeans_ncells,eta_sn,yield,rbubble,f_ek,ndebris,f_w,mass_gmc     &
       & ,J21,a_spec,z_ave,z_reion,n_sink,ind_rsink                       &
       & ,delayed_cooling,self_shielding,smbh,agn,rsink_max,msink_max     &
       & ,units_density,units_time,units_length,neq_chem,ir_feedback      &
       & ,ir_eff,merge_stars                                              &
       & ,larson_lifetime,flux_accretion,t2thres_sf,sn_delay_myr          &
       & ,sn_dt2_min,mseed, x_floor, sigmav_max,r_gal,T2max_AGN           &
       & ,t_delay,t_delcool,sig_delcool                                   &
       & ,sink_agn,sinkprops,bondi,drag,delayed_cooling,sf_hopkins        &
       & ,star_maker, mechanical_feedback, star_imf, sn2_real_delay       &
       & ,log_mfb,M_SNII

  logical:: rt_smooth, SED_isEgy, rt_isoPress, rt_isIR                    &
       & ,rt_isOpt, rt_output_coolstats, rt_star, rt_isIRtrap=.false.     &
       & ,isH2,rt_AGN                                                     &
       & ,is_kIR_T,rt_metal_cooling, rt_vc, rt_T_rad, rt_isIRTrapPress    &
       & ,rt_cool_eq, extra_output
  real(kind=8)::rt_esc_frac,rt_courant_factor   &
       & ,rt_c_fraction,rt_pressBoost,rt_UV_nhss, rt_uvsrc_nhmax &
       & ,aexp_nocosmo
  integer::SEDprops_update, rt_nregion, rt_nsource, nUVgroups
  integer,dimension(1,50)::rt_reg_group,rt_src_group
  real(kind=8),dimension(1,50)::rt_reg_x_center,rt_reg_y_center           &
       ,rt_reg_z_center, rt_reg_length_x, rt_reg_length_y                 &
       ,rt_reg_length_z, rt_exp_region, rt_n_region, rt_u_region          &
       ,rt_v_region, rt_w_region, rt_xion_region                          &
       ,rt_src_x_center,rt_src_y_center                                   &
       ,rt_src_z_center, rt_src_length_x, rt_src_length_y                 &
       ,rt_src_length_z, rt_exp_source, rt_n_source, rt_u_source          &
       ,rt_v_source, rt_w_source, rt_wind_source,X,Y,rt_kscir,cs_kn

  character(LEN=80)::rt_flux_scheme,SED_dir,UV_file,hll_evals_file
  character(LEN=80),dimension(50)::rt_region_type,rt_source_type
  integer,dimension(50)::rt_movie_vars
  namelist/rt_params/rt_smooth, SED_isEgy, rt_isoPress, rt_isIR           &
       & ,rt_output_coolstats, rt_star, rt_otsa, rt_isIRtrap              &
       & ,SEDprops_update                                                 &
       & ,rt_esc_frac, rt_courant_factor                                  &
       & ,rt_c_fraction,rt_pressBoost                                     &
       & ,rt_flux_scheme, SED_dir, UV_file,rt_UV_nhss                     &
       & ,rt_nregion,rt_nsource,rt_region_type,rt_source_type             &
       & ,rt_reg_x_center, rt_reg_y_center                                &
       & ,rt_reg_z_center, rt_reg_length_x, rt_reg_length_y               &
       & ,rt_reg_length_z, rt_exp_region, rt_n_region, rt_u_region        &
       & ,rt_v_region, rt_w_region, rt_xion_region, rt_reg_group          &
       & ,rt_src_x_center,rt_src_y_center                                 &
       & ,rt_src_z_center, rt_src_length_x, rt_src_length_y, rt_src_group &
       & ,rt_src_length_z, rt_exp_source, rt_n_source, rt_u_source        &
       & ,rt_v_source, rt_w_source, rt_wind_source,X,Y,isH2,rt_kscir      &
       & ,cs_kn, rt_AGN, is_kIR_T, rt_metal_cooling, nUVgroups            &
       & ,hll_evals_file, rt_vc, rt_T_rad, rt_isIRtrap, rt_uvsrc_nhmax    &
       & ,aexp_nocosmo, rt_cool_eq, extra_output, rt_movie_vars

  real(kind=8),dimension(50,50)::group_csn,group_cse
  real(kind=8),dimension(50)::kappaAbs,kappaSc=1d1,groupL0,groupL1,group_egy
  namelist/rt_groups/kappaAbs, kappaSc, groupL0, groupL1, group_egy       &
       ,group_csn,group_cse

  ! Cooling
  integer :: cool_format=1
  integer, parameter :: in=161,iT=101
  double precision, dimension(in) :: cool_n
  double precision, dimension(iT) :: cool_T
  double precision, dimension(in,iT,6) :: cool_spec
  real(KIND=8), allocatable, dimension(:,:) :: spec

CONTAINS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE read_info(nsnap)

! Read info files
!-------------------------------------------------------------------------
  character*50::filename,dirname
  character(len=128)::nsnap,ifmt,dfmt,varname
!  character*5::nchar
  integer::iargc,ipos
  logical::ok
  integer::ilun, i, impi
  real::value
!-------------------------------------------------------------------------
  ilun=10
  dirname=trim(nsnap)//'/'
  ipos=INDEX(dirname,'output_')
  nchar=dirname(ipos+7:ipos+13)

  ! Try to read namelist (for eta_sn)
  INQUIRE(file='nml.nml',exist=nml_ok)
  if(.not. nml_ok)then
     write(*,*)'Could not open namelist'
  else
     open(1,file='nml.nml')
     rewind(1)
     read(1,NML=physics_params,END=101)
     rewind(1)
     read(1,NML=rt_params,END=101)
     rewind(1)
     read(1,NML=rt_groups,END=101)
101 continue                                   ! No harm if no rt namelist
     close(1)
  end if

  print*,'read the namelist, and eta_sn = ',eta_sn

  ! Read info file
  ! Open HYDRO file and read gamma
  filename=TRIM(dirname)//'/hydro_'//TRIM(nchar)//'.out00001'
  inquire(file=filename, exist=ok)
  if (.not. ok) then
     write(*,*)'File '//TRIM(filename)//' not found'
     gamma=1.4
     !STOP
  else
     open(unit=ilun,file=filename,status='old',form='unformatted')
     read(ilun)
     read(ilun) nhvar
     read(ilun)
     read(ilun)
     read(ilun)
     read(ilun) gamma
     close(ilun)
  endif

  filename=TRIM(dirname)//'info_'//TRIM(nchar)//'.txt'
  inquire(file=filename, exist=ok)
  if (.not. ok) then
     write(*,*)'File '//TRIM(filename)//' not found'
     STOP
  else
     open(unit=ilun,file=filename,status='old',form='formatted')
     read(ilun,'(13x,I11)')ncpu
     read(ilun,'(13x,I11)')ndim
     read(ilun,'(13x,I11)')levelmin
     read(ilun,'(13x,I11)')nlevelmax
     read(ilun,'(13x,I11)')ngridmax
     read(ilun,'(13x,I11)')nstep_coarse
     read(ilun,*)
     read(ilun,'(13x,E23.15)')scale
     read(ilun,'(13x,E23.15)')t
     read(ilun,'(13x,E23.15)')aexp
     read(ilun,'(13x,E23.15)')h0
     read(ilun,'(13x,E23.15)')omega_m
     read(ilun,'(13x,E23.15)')omega_l
     read(ilun,'(13x,E23.15)')omega_k
     read(ilun,'(13x,E23.15)')omega_b
     read(ilun,'(13x,E23.15)')unit_l
     read(ilun,'(13x,E23.15)')unit_d
     read(ilun,'(13x,E23.15)')unit_t
     read(ilun,*)
     read(ilun,'(14x,A80)')ordering
     unit_m=unit_d*unit_l**3
     hubble = h0*1d5/Mpc*sqrt(omega_m/aexp**2 + omega_l)*unit_t
     
     !read(10,'(A14,A80)')tmp,ordering
     write(*,'(" ordering type=",A20)'),TRIM(ordering)
     read(10,*)
     if(TRIM(ordering).eq.'hilbert')then
        allocate(bound_key(0:ncpu))
        do impi=1,ncpu
           read(10,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
        end do
     endif
     close(ilun)
  endif

  ! Try to read rt info file
  filename=TRIM(dirname)//'info_rt_'//TRIM(nchar)//'.txt'
  inquire(file=filename, exist=ok)
  if (.not. ok) then
     write(*,*)'File '//TRIM(filename)//' not found'
  else
     open(unit=ilun,file=filename,status='old',form='formatted')
     call read_int( ilun, 'nRTvar', nRTvar)
     call read_int( ilun, 'nIons', nIons)
     call read_int( ilun, 'nGroups', nGroups)
     call read_int( ilun, 'iIons', iIons)
     call read_real(ilun, 'X_fraction', X_frac)
     call read_real(ilun, 'Y_fraction', Y_frac)
     call read_real(ilun, 'unit_np', unit_np)
     call read_real(ilun, 'unit_pf', unit_fp)
     call read_real(ilun, 'rt_c_frac', rt_c_frac)
     call read_real(ilun, 'n_star', n_star)
     call read_real(ilun, 'T2_star', T2_star)
     call read_real(ilun, 'g_star', g_star)
     call read_groups(ilun)
     close(ilun)
  endif
  !nHvar=iIons+nIons-1

!  if(eps_star>0) &
!     t_star=0.1635449*(n_star/0.1)**(-0.5)/eps_star

  unit_v  = unit_l/unit_t
  unit_P  = unit_d*unit_v**2
  unit_T2 = mH/kB*unit_v**2
  unit_nH = X_frac*unit_d/mH
  !unit_nH = unit_d/mH
  t_myr   = t*unit_t/31.55693d12

  rt_c_cgs = c_cgs*rt_c_frac
  rt_c = rt_c_cgs/unit_v

  iPrad=ndim+2

  iPT = ndim+2 ! Thermal pressure variable index
  if(rt_isIRTrap) iPT=iPT+1 ! one up if trapping.
  if(rt_isIRTrap) print*,'Trapping is activated!!!'

  !iPT=iPT+1 ! one up if trapping.

  rt_kIR=kappaSc(1)

END SUBROUTINE read_info

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE read_real(lun, param_name, value)

! Try to read a parameter from lun
!-------------------------------------------------------------------------
  integer::lun
  character(*)::param_name
  character(128)::line,tmp
  real::value
!-------------------------------------------------------------------------
  rewind(unit=lun)
  do
     read(lun, '(A128)', end=222) line
     if(index(line,trim(param_name)) .eq. 1) then
        read(line,'(A13,E23.15)') tmp, value
        return
     endif
  end do
222 return                        ! eof reached, didn't find the parameter

END SUBROUTINE read_real

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE read_groups(lun)

! Try to read photon group properties from lun
!-------------------------------------------------------------------------
  integer::lun,i
  character(128)::line,tmp
  real,dimension(:),allocatable::group_tmp
!-------------------------------------------------------------------------
  !allocate(group_egy(nGroups))
  allocate(group_tmp(nIons))
  rewind(unit=lun)
  i=0 ! Read group_egy
  do
     read(lun, '(A128)', end=220) line
     if(index(line,trim('egy')) .eq. 3) then
        i=i+1
        read(line,'(A19,100F12.3)') tmp, group_egy(i)
     endif
  end do
220 continue
  group_egy=group_egy*eV_to_erg
  i=0 ! Read group_csn
  rewind(unit=lun)
  do
     read(lun, '(A128)', end=221) line
     if(index(line,trim('csn')) .eq. 3) then
        i=i+1
        read(line,'(A19,3(1pe12.3))') tmp, group_tmp
        group_csn(i,1:nIons)=group_tmp
     endif
  end do
221 continue
  i=0 ! Read group_cse
  rewind(unit=lun)
  do
     read(lun, '(A128)', end=222) line
     if(index(line,trim('cse')) .eq. 3) then
        i=i+1
        read(line,'(A19,100(1pe12.3))') tmp, group_tmp
        group_cse(i,1:nIons) = group_tmp
     endif
  end do
222 continue
  deallocate(group_tmp)
  return                        ! eof reached, didn't find the parameter

END SUBROUTINE read_groups

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE read_int(lun, param_name, value)

! Try to read a parameter from lun
!-------------------------------------------------------------------------
  integer::lun
  character(*)::param_name
  character(128)::line,tmp
  integer::value
!-------------------------------------------------------------------------
  rewind(unit=lun)
  do
     read(lun, '(A128)', end=223) line
     if(index(line,trim(param_name)) .eq. 1) then
        read(line,'(A13,I30)') tmp, value
        return
     endif
  end do
223 return                        ! eof reached, didn't find the parameter

END SUBROUTINE read_int

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION getTemperature(rho, p, xHII, xHeII, xHeIII)

! Calculate Temperature in cell
! All in arguments are in user units and the result is in Kelvin.
!-------------------------------------------------------------------------
  real(kind=8),intent(in)::rho, p, xHII, xHeII, xHeIII
  real(kind=8)::getTemperature
  real(kind=8),parameter::Tmin=1d-2
!-------------------------------------------------------------------------
  getTemperature = p/rho * unit_T2 &
                 - T2_star*(rho*unit_nH/n_star)**(g_star-1d0) !New polytr.
  getTemperature = MAX(getTemperature, Tmin)
  getTemperature = getTemperature * getMu(xHII, xHeII, xHeIII)
END FUNCTION getTemperature

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION getTmu(rho, p, subPol)

! Calculate T/mu in cell
! All in arguments are in user units and the result is in Kelvin.
!-------------------------------------------------------------------------
  real(kind=8),intent(in)::rho, p
  logical,intent(in)::subPol
  real(kind=8)::getTmu
!-------------------------------------------------------------------------
  getTmu = p/rho * unit_T2
  if(subPol) getTmu = getTmu - T2_star*(rho*unit_nH/n_star)**(g_star-1d0)
END FUNCTION getTmu

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION getMu(xHII, xHeII, xHeIII)

! Calculate average mass per particle in units of mH
!-------------------------------------------------------------------------
  real(kind=8),intent(in)::xHII, xHeII, xHeIII
  real(kind=8)::getMu
!-------------------------------------------------------------------------
  getMu = 1./( X_frac*(1.d0+xHII) + 0.25*Y_frac*(1.d0+xHeII+2.d0*xHeIII) )
END FUNCTION getMu

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION isPolytrope(rho, p)

! Calculate whether cells are polytropic (1.) or not (0.)
! All in arguments are in user units
!-------------------------------------------------------------------------
  real(kind=8),intent(in)::rho, p
  real(kind=8)::isPolytrope
  real(kind=8)::T2min,T2
!-------------------------------------------------------------------------
  T2min = T2_star*(rho*unit_nh/n_star)**(g_star-1.) + 10.d0
  T2 = p/rho*unit_T2
  isPolytrope=0.
  if(T2 .lt. T2min) isPolytrope=1.
END FUNCTION isPolytrope

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION get_AlphaA_HII(TK)

! Returns case A rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)
! TK           => Temperature [K]
!------------------------------------------------------------------------
  real(kind=8),intent(in)::TK
  real(kind=8)::get_AlphaA_HII, lambda
!------------------------------------------------------------------------
  lambda = 315614./TK
  get_AlphaA_HII =  1.269d-13 * lambda**1.503 &
                     / ( ( 1.d0+(lambda/0.522)**0.47 )**1.923 )
END FUNCTION get_AlphaA_HII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION get_AlphaB_HII(TK)

! Returns case B rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)
! TK           => Temperature [K]
!------------------------------------------------------------------------
  real(kind=8),intent(in)::TK
  real(kind=8)::get_AlphaB_HII, lambda
!------------------------------------------------------------------------
  lambda = 315614./TK
  get_AlphaB_HII =                                                       &
       2.753d-14 * lambda**1.5 / ( (1.d0+(lambda/2.74)**0.407)**2.242 )
END FUNCTION get_AlphaB_HII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION get_collExrate_HI(TK)

! Gives Collisional Excitation rate coefficient for the 1->2
! transition in hydrogen atom energies, in [cm3 s-1], according to
! Callaway, Unnikrishnan & Oza 1987.
! TK      => Temperatue in Kelvin
!-------------------------------------------------------------------------
  real(kind=8),intent(in)::TK
  real(kind=8)::get_collExrate_HI
!-------------------------------------------------------------------------
  get_collExrate_HI =                                                   &
       2.41d-6/sqrt(TK) * (TK/1.d4)**0.22 * exp(-1.63d-11/(kb*TK))
END FUNCTION get_collExrate_HI

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION getRecLyaEmi(TK, rho, p, xHII, xHeII, xHeIII)

! Calculate recombinative Lyman alpha emissivity of cells.  All in
! arguments are in user units, except TK which is in K, and the result is
! in erg per sec per cc.
!-------------------------------------------------------------------------
  real(kind=8),intent(in)::TK, rho, p, xHII, xHeII, xHeIII
  real(kind=8)::getRecLyaEmi
  real(kind=8),parameter::f=0.68d0
  real(kind=8)::nH, nHe, nHII, n_e
!-------------------------------------------------------------------------
  if (noPol) then
     if(isPolytrope(rho, p) .gt. 0.d0) then
        getRecLyaEmi = 0.d0
        return
     endif
  endif

  nH   = rho*unit_nH
  nHe  = 0.25*nH*Y_frac/X_frac
  nHII = xHII*nH
  n_e  = nHII+nHe*(xHeII+2.d0*xHeIII)

  getRecLyaEmi = f * nHII * n_e * get_AlphaB_HII(TK) *  e_lya
END FUNCTION getRecLyaEmi

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION getCollLyaEmi(TK, rho, p, xHII, xHeII, xHeIII)

! Calculate collisional Lyman alpha emissivity of cells.  All in
! arguments are in user units, except TK which is in K, and the result is
! in erg per sec per cc.
!-------------------------------------------------------------------------
  real(kind=8),intent(in)::TK, rho, p, xHII, xHeII, xHeIII
  real(kind=8)::getCollLyaEmi
  real(kind=8)::nH, nHe, nHII, n_e
!-------------------------------------------------------------------------
  if (noPol) then
     if(isPolytrope(rho, p) .gt. 0.d0) then
        getCollLyaEmi = 0.d0
        return
     endif
  endif

  nH   = rho*unit_nH
  nHe  = 0.25d0*nH*Y_frac/X_frac
  nHII = xHII*nH
  n_e  = nHII+nHe*(xHeII+2.d0*xHeIII)

  getCollLyaEmi = (nH-nHII) * n_e * get_collExrate_HI(TK) * e_lya
END FUNCTION getCollLyaEmi

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION getcPyy(cE, Fx, Fy)

! Returns yy element of the Eddington tensor (in code units)
!-------------------------------------------------------------------------
  real(kind=8),intent(in)::cE, Fx, Fy
  real(kind=8)::getcPyy
  real(kind=8)::chi, fred, Fmag, nbold
!-------------------------------------------------------------------------
  Fmag = sqrt( Fx**2 + Fy**2 )
  fred = Fmag/cE
  chi = (3d0 + 4d0*fred**2)/(5d0+2d0*sqrt(4d0-3d0*fred**2))
  nbold = Fy / Fmag
  getcPyy = (1d0-chi)/2d0 + (3d0*chi-1d0)/2d0 * nbold**2
  getcPyy = getcPyy * cE
END FUNCTION getcPyy

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE getPHrate(var,ind,emi)

! Calculate photoheating rate in cells for all photon groups.
! The result (emi) is in erg per sec per cc.
!-------------------------------------------------------------------------
  real(kind=8)::var(:,:,:), emi(:)
  integer::ind,iGroup,indexGroup,iIon,sz
  real(kind=8),allocatable::nH(:), nHe(:), nN(:,:)
  real(kind=8)::homPHrates(3)= (/4.719772E-24, 7.311075E-24, 8.531654E-26/)
!-------------------------------------------------------------------------
  sz = size(emi)
  allocate(nH(sz))
  allocate(nHe(sz))
  allocate(nN(3,sz))
  nH    = var(:,ind,1)*unit_nH
  nN(1,:) = nH * max(0.,1.d0-var(:,ind,iIons))        !  HI number density
  nHe   = 0.25d0 * nH * Y_frac/X_frac
  nN(2,:) = nHe * max(0.,1.d0-var(:,ind,iIons+1)-var(:,ind,iIons+2)) ! HeI
  nN(3,:) = nHe * var(:,ind,iIons+1)                  !HeII number density

  emi=0.d0
  if(nGroups .eq. 0) then
     do iIon=1,3
        emi = emi + nN(iIon,:) *   homPHrates(iIon)
     enddo
  else
     do iGroup=1,nGroups
        indexGroup=nHvar+1+(iGroup-1)*(1+ndim)
        do iIon=1,3
           emi = emi                                                     &
                + nN(iIon,:) * var(:,ind,indexGroup) * unit_Fp           &
                *(group_cse(iGroup,iIon)*group_egy(iGroup)               &
                -group_csn(iGroup,iIon)*ion_egy(iIon))
        enddo
     enddo
  endif
  deallocate(nH,nHe,nN)
END SUBROUTINE getPHrate

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE getPIrate(var,ind,emi,iIon)

! Calculate photoionization rate in cells for all photon groups.
! The result (emi) is in s-1
!-------------------------------------------------------------------------
  real(kind=8)::var(:,:,:), emi(:)
  integer::ind,iSpecies,iGroup,indexGroup,iIon
  real(kind=8)::homPIrates(3)= (/4.719772E-24, 7.311075E-24, 8.531654E-26/)
!-------------------------------------------------------------------------
  emi=0.d0
  if(nGroups .eq. 0) then
     !do iIon=1,3
        emi = emi + homPIrates(iIon)
     !enddo
  else
     do iGroup=1,nGroups
        indexGroup=nHvar+1+(iGroup-1)*(1+ndim)
        !do iIon=1,3
           emi = emi                                                     &
               ! + var(:,ind,indexGroup) * unit_Np * rt_c_cgs            &
                + var(:,ind,indexGroup) * unit_Fp  & ! Fix bc output is cN
                * group_csn(iGroup,iIon)
        !enddo
     enddo
  endif
END SUBROUTINE getPIrate

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE calcSFR(var,ind,sfr)

! Calculate star formation rate density, in Msun/yr/kpc^3
! The result is returned in sfr
!-------------------------------------------------------------------------
  real(kind=8)::var(:,:,:), sfr(:),t0_yr
  integer::ind,sz
  real(kind=8),allocatable::nH(:),Tmu(:)
!-------------------------------------------------------------------------
  sz = size(sfr)
  allocate(nH(sz))
  allocate(Tmu(sz))
  nH = var(:,ind,1) * unit_nH
  t0_yr=t_star*1d9
  sfr = nH / (t0_yr*sqrt(n_star/nH)) ! atoms/yr/cm3  Msun/yr/cm3
  sfr = sfr /unit_nH*unit_d/1.98892d33 * (3.08d21)**3 !->Msun/yr/kpc3
  Tmu = getTmu( var(:,ind,1), var(:,ind,iPT), .true. )
  where (Tmu .ge. 3d3) sfr=0d0
  where (nH .lt. n_star) sfr=0d0
  deallocate(nH,Tmu)
END SUBROUTINE calcSFR

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION getSFR(rho,P)

! Calculate star formation rate density, in Msun/yr/cm^3
!-------------------------------------------------------------------------
  real(kind=8),intent(in)::rho,P
  real(kind=8)::getSFR
  real(kind=8)::nH,Tmu,t0_yr
!-------------------------------------------------------------------------
  nH = rho * unit_nH
  if (nH .lt. n_star) getSFR=0d0
  Tmu = getTmu( rho, P, .true. )
  if (Tmu .ge. 3d3) getSFR=0d0
  t0_yr=t_star*1d9
  getSFR = nH / (t0_yr*sqrt(n_star/nH)) ! atoms/yr/cm3
  getSFR = getSFR /unit_nH*unit_d/1.98892d33 !->Msun/yr/cm3

END FUNCTION getSFR

!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!SUBROUTINE calc_CIE_xHII(var,ind)
!
!! Replace hydrogen ionization states with CIE ionization states
!! vars     => Cell variables (rho, v, u, w, etc)
!!------------------------------------------------------------------------
!  real(kind=8)::var(:,:,:)
!  integer::ind,i
!  real(kind=8),dimension(nIons)::xion
!  real(kind=8)::nH, T2, mu
!  real(kind=8),dimension(nIons)::phI_rates      ! Photoionization rates [s-1]
!  real(kind=8),dimension(6)::nSpec              !          Species abundances
!!------------------------------------------------------------------------
!  pHI_rates = 0.d0
!  do i=1,size(var,1)
!     nH=var(i,ind,1)*unit_nH
!     if(nh .ge. 1.d-2) then
!        T2=var(i,ind,iPT)/var(i,ind,1) * unit_T2
!        call cmp_Equilibrium_Abundances(T2, nH, pHI_rates, mu, nSpec)
!        var(i,ind,iIons) = nSpec(3)/(nSpec(2)+nSpec(3))
!     endif
!  end do
!END SUBROUTINE calc_CIE_xHII

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE getAccRad_dust(var, ind, emi, dx)

! Calculate local acceleration of gas due to radiation pressure.
! We assume acc_rad = Sum(F_group*kappa_group*Z/c) + 1/3 Erad/(dx*rho).
!-------------------------------------------------------------------------
  real(kind=8)::var(:,:,:), emi(:), dx
  integer::ind,i,iGroup!,sz
!-------------------------------------------------------------------------
  !sz = size(emi)

  do i=1,nGroups ! Radiation acceleration [cm s-2]
     iGroup=nHvar+1+(i-1)*(1+ndim)
     emi = emi + var(:,ind,iGroup) * unit_Fp / 3d10                      &
               * max(kappaAbs(i),kappaSc(i))* group_egy(i)               &
               * var(:,ind,iPT+1) / 0.02 ! Metallicity (solar)
  enddo

END SUBROUTINE getAccRad_dust

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE getAccRad_trapped(var, ind, emi, dx)

! Calculate local acceleration of gas due to trapped radiation pressure.
! We assume acc_rad = 1/3 Erad/(dx*rho).
!-------------------------------------------------------------------------
  real(kind=8)::var(:,:,:), emi(:), dx
  integer::ind
!-------------------------------------------------------------------------
  if(.not. rt_isIRTrap) return ! No trapping
  emi = emi + var(:,ind,ndim+2) / var(:,ind,1) / dx / scale &
       * unit_P / unit_d / unit_l

END SUBROUTINE getAccRad_trapped

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE getAccRad_direct(var,ind,emi)

! Calculate local acceleration of gas due to direct radiation pressure.
! We assume acc_rad = Sum(F_group/c*csn_group*xion/mp).
!-------------------------------------------------------------------------
  real(kind=8)::var(:,:,:), emi(:)
  integer::ind,iGroup,indGroup,iIon,sz
  real(kind=8),allocatable::xn(:,:)
!-------------------------------------------------------------------------
  sz = size(emi)
  allocate(xn(3,sz))
  ! Calculate neutral fractions
  xn(1,:) = max(0.,1.d0-var(:,ind,iIons))
  xn(2,:) = max(0.,1.d0-var(:,ind,iIons+1)-var(:,ind,iIons+2))
  xn(3,:) = var(:,ind,iIons+1)

  do iGroup=1,nGroups
     indGroup=nHvar+1+(iGroup-1)*(1+ndim)
     do iIon=1,3
        emi = emi                                                        &
             + var(:,ind,indGroup) * unit_Fp * group_egy(iGroup) /3d10   &
             * group_csn(iGroup,iIon) * xn(iIon,:) / 1.6726d-24
     enddo
  enddo
  deallocate(xn)
END SUBROUTINE getAccRad_direct

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE emissivity(var, emi, ind, dx,ngrida, temp_in)

! Calculate emissivity of a variable in cells.
! All the cells should have the same dimensions, dx
!-------------------------------------------------------------------------
  real(kind=8),intent(out),optional::temp_in(:)
  real(kind=8)::var(:,:,:), emi(:), hnu, ion_frac, dx
  real(kind=8),allocatable::tmp(:),tmp2(:), TK(:)
  integer::ind, i, ngrida
!-------------------------------------------------------------------------
  !calc_TK=.true. !DEBUG joki
  if(calc_TK) then
     allocate(TK(size(emi)))
     TK = getTemperature( var(:,ind,1), var(:,ind,iPT)                   &
          ,var(:,ind,iIons), var(:,ind,iIons+1), var(:,ind,iIons+2) )
     ! TK is now temperature in Kelvin
  endif
  !BEGIN DEBUG: Trying out a maximum Fluorescence Lya emissivity
  !call calc_CIE_xHII(var,ind)
  !var(:,ind,iIons)=1. ! H is completely ionized
  !END DEBUG

  if (present(temp_in)) then
    temp_in=getTmu( var(:,ind,1), var(:,ind,iPT), .true. )
  endif

  SELECT CASE (emType)
     CASE ('oneVar') ! One variable times a unit
        emi = var(:,ind,emVar)*defUnit - defSubtr
     CASE ('A') ! One variable times a unit - CATCH NEGATIVE VALUES ONLY
        emi = var(:,ind,emVar)*defUnit
        emi=min(emi,0d0)
        emi=-emi
     CASE ('b') ! Speed
        !emi = (var(:,ind,2)+4554021.4/unit_v)**2
        !if(ndim .gt. 1) emi = emi + (var(:,ind,3)+3532903.5/unit_v)**2
        !if(ndim .gt. 2) emi = emi + (var(:,ind,4)+9367328.1/unit_v)**2
        emi = var(:,ind,2)**2
        if(ndim .gt. 1) emi = emi + var(:,ind,3)**2
        if(ndim .gt. 2) emi = emi + var(:,ind,4)**2
        emi = sqrt(emi)*unit_v
     CASE ('c') ! Tmu
        emi=getTmu( var(:,ind,1), var(:,ind,iPT), .true. )
     CASE ('Tmunt') ! 'Nonthermal' temperature
        emi=getTmu( var(:,ind,1), var(:,ind,iPrad), .false. )
     CASE ('d') ! T
        emi = TK
     CASE ('e') ! Ekin
        emi = var(:,ind,2)**2
        if(ndim .gt. 1) emi = emi + var(:,ind,3)**2
        if(ndim .gt. 2) emi = emi + var(:,ind,4)**2
        emi = 0.5d0 * emi*unit_v**2 * var(:,ind,1)*unit_d
     CASE ('f') ! Mach
        emi = var(:,ind,2)**2
        if(ndim .gt. 1) emi = emi + var(:,ind,3)**2
        if(ndim .gt. 2) emi = emi + var(:,ind,4)**2
        emi=sqrt(emi)/sqrt(gamma*var(:,ind,iPT)/var(:,ind,1))
     CASE ('g') ! xHI
        emi = 1.d0-var(:,ind,emVar)
     CASE ('h') ! xHeI
        emi = 1.d0-var(:,ind,iIons+1)-var(:,ind,iIons+2)
     CASE ('L') ! Lya
        emi = getCollLyaEmi(TK, var(:,ind,1), var(:,ind,iPT)             &
              ,var(:,ind,iIons), var(:,ind,iIons+1), var(:,ind,iIons+2)) &
            + getRecLyaEmi(TK, var(:,ind,1), var(:,ind,iPT)              &
              ,var(:,ind,iIons), var(:,ind,iIons+1), var(:,ind,iIons+2))
     CASE ('C') ! CLya
        emi = getCollLyaEmi(TK, var(:,ind,1), var(:,ind,iPT)             &
                ,var(:,ind,iIons), var(:,ind,iIons+1), var(:,ind,iIons+2))
     CASE ('R') ! RLya
        emi = getRecLyaEmi(TK, var(:,ind,1), var(:,ind,iPT)              &
                ,var(:,ind,iIons), var(:,ind,iIons+1), var(:,ind,iIons+2))
     CASE ('i') ! Polytrope
        emi = isPolytrope(var(:,ind,1), var(:,ind,iPT))
     CASE ('j') ! Entropy inversed
        emi=var(:,ind,1)**(5./3.)/var(:,ind,iPT)
     CASE ('k') ! Ion density
        emi = var(:,ind,1)*(1.d0-var(:,ind,emVar))*defUnit
     CASE ('K') ! Ion density
        emi = var(:,ind,1)*(1.d0-var(:,ind,emVar)-var(:,ind,emVar-1))*defUnit
     CASE ('l') ! Grav luminosity
        ! Calculate Lya cooling rate
        emi = getCollLyaEmi(TK, var(:,ind,1), var(:,ind,iPT)             &
              ,var(:,ind,iIons), var(:,ind,iIons+1), var(:,ind,iIons+2)) &
            + getRecLyaEmi(TK, var(:,ind,1), var(:,ind,iPT)              &
              ,var(:,ind,iIons), var(:,ind,iIons+1), var(:,ind,iIons+2))
        allocate(tmp(size(emi))) ; call getPHrate(var, ind, tmp)
        emi = emi - tmp                       ! Subtract photoheating rate
        call getPIrate(var, ind, tmp,1)
        tmp = tmp * 0.68d0                                               &
                    * var(:,ind,1)*(1.d0-var(:,ind,iIons))*unit_nH * e_lya
        emi = emi - tmp             ! Subtract photoionzation contribution
        deallocate(tmp)
        emi=max(0.d0,emi)                         ! Should not be negative
     CASE ('m') ! UV -> lya luminosity
        ! Calculate Lya cooling rate
        emi = getCollLyaEmi(TK, var(:,ind,1), var(:,ind,iPT)             &
              ,var(:,ind,iIons), var(:,ind,iIons+1), var(:,ind,iIons+2)) &
            + getRecLyaEmi(TK, var(:,ind,1), var(:,ind,iPT)              &
              ,var(:,ind,iIons), var(:,ind,iIons+1), var(:,ind,iIons+2))
        allocate(tmp(size(emi)))  ; call getPHrate(var, ind, tmp)
        allocate(tmp2(size(emi))) ; call getPIrate(var, ind, tmp2,1)
        tmp2= tmp2 * 0.68d0                                              &
                    * var(:,ind,1)*(1.d0-var(:,ind,iIons))*unit_nH * e_lya
        emi=min(emi,tmp+tmp2)              ! Should not be more than Lyemi
        deallocate(tmp,tmp2)
     CASE ('n') ! Photoheating
        call getPHrate(var, ind, emi)
     CASE ('o') ! Photoionization rate
        call getPIrate(var, ind, emi,1)
     CASE ('p') ! Ion density
        emi = var(:,ind,1)*(var(:,ind,emVar))*defUnit
     CASE ('q') ! Gas flux magnitude [g cm-2]
        emi = var(:,ind,2)**2                          ! x velocity
        if(ndim .gt. 1) emi = emi + var(:,ind,3)**2    ! y velocity
        if(ndim .gt. 2) emi = emi + var(:,ind,4)**2    ! z velocity
        emi = sqrt(emi)*var(:,ind,1)*unit_v*unit_d
     CASE ('r') ! Gas flux component (x,y,z) [g cm-2]
        emi = var(:,ind,1)*(var(:,ind,emVar))*unit_v*unit_d

     CASE ('acc_rad') ! Acceleration of gas due to radiation
        emi=0d0
        call getAccRad_dust(var, ind, emi, dx)
        call getAccRad_trapped(var, ind, emi,dx)
        call getAccRad_direct(var, ind, emi)

     CASE ('Fpreduced') ! Reduced photon flux
        emi = var(:,ind,emvar+1)**2
        if(ndim .gt. 1) emi = emi + var(:,ind,emvar+2)**2
        if(ndim .gt. 2) emi = emi + var(:,ind,emvar+3)**2
        emi=sqrt(emi)
        emi=emi/var(:,ind,emvar)
     CASE ('Fpyreduced') ! Reduced photon flux in y direction
        emi = var(:,ind,emvar+2)/var(:,ind,emvar)
     CASE ('FeIR') ! Streaming + trapped IR photon flux
        emi = var(:,ind,nHvar+1)*unit_Fp*group_egy(1)   ! Streaming
        emi = emi + var(:,ind,iPrad) * unit_P*c_cgs*3d0 ! Trapped
     CASE ('t') ! SFR density
        call calcSFR(var, ind, emi)
        emi=emi*defUnit
     CASE ('tauIR') ! tauIR
        emi = var(:,ind,1)*unit_d*rt_kIR*var(:,ind,iPT+1)/0.02           &
             * dx*scale*unit_l
     CASE ('StauIR') ! StauIR
        emi = var(:,ind,1)*unit_d*rt_kIR*var(:,ind,iPT+1)/0.02           &
             * scale*unit_l

     CASE ('TrKrumh')
        emi = var(:,ind,nHvar+1)*unit_Fp*group_egy(1)
        ! Add trapped photon density?
        if(rt_isIRtrap) emi=emi +var(:,ind,ndim+2)                       &
                                *unit_d*unit_v**2*c_cgs*3d0
        Emi=(Emi/2.99d10/7.56d-15)**0.25 / 82d0

     CASE ('TgKrumh')
        emi = var(:,ind,iPT)/var(:,ind,1)*unit_T2*2.33 / 82d0

     CASE ('TratioKrumh') ! TrKrumh/TgKrumh
        emi = var(:,ind,nHvar+1)*unit_Fp*group_egy(1) ! Tr
        ! Add trapped photon density?
        if(rt_isIRtrap) emi=emi+var(:,ind,ndim+2)                        &
                               *unit_d*unit_v**2*c_cgs*3d0
        Emi=(Emi/2.99d10/7.56d-15)**0.25              ! Tr
        emi = emi/ ( var(:,ind,iPT)/var(:,ind,1)*unit_T2*2.33 ) !Tr/Tg

     CASE ('fKrumh') ! f_rad / f_g
        emi = var(:,ind,iPT)/var(:,ind,1)*unit_T2*2.33  ! Tg
        emi = 0.0316 * (emi/1d1)**2                     ! Kappa
        emi = emi * var(:,ind,nHvar+3) * unit_Fp / 3d10 ! acc_rad
        emi = emi / 1.46d-6                             ! acc_rad / acc_grav
     CASE ('tauKrumh') ! tau_cell
        emi = var(:,ind,iPT)/var(:,ind,1)*unit_T2*2.33  ! Tg
        emi = 0.0316 * (emi/1d1)**2                     ! Kappa [cm2/g]
        emi = emi * var(:,ind,1) * dx * scale          &! tau_cell
            * unit_d * unit_l

     CASE ('Trad')
        emi = var(:,ind,nHvar+1)*unit_Fp*group_egy(1)
        ! Add trapped photon density?
        if(rt_isIRtrap) emi=emi+var(:,ind,ndim+2)*unit_d*unit_v**2*c_cgs*3d0
        Emi=(Emi/2.99d10/7.56d-15)**0.25

     CASE ('cPyy') ! Pyy/E (Eddington pressure tensor element)
        emi = getcPyy(var(:,ind,nHvar+1), var(:,ind,nHvar+2) &
             , var(:,ind,nHvar+3)) / var(:,ind,nHvar+1)

     CASE ('mass')
        emi = var(:,ind,1) *unit_d               &     ! rho in g/cm3
             * (dx * scale * unit_l)**ndim       &     ! Cell volume in cm3
             / 2d33                                    ! Units of Msun
            ! / 1.98892d33                              ! Units of Msun

     CASE ('Zmass')
        emi = var(:,ind,1)*var(:,ind,iPT+1) * unit_d & ! rho * Z in g/cm3
             * (dx * scale * unit_l)**ndim           & ! Cell volume
             / 2d33                                    ! Units of Msun
            ! / 1.98892d33                              ! Units of Msun

     CASE ('sfr_Ha') ! Straight integral of star formation density, as calculated in NUT.
                  ! Corresponds to density of cell / free-fall time of cell * efficiency
                  ! Converted to Halpha from Canonical IMF fit in Pflamm-Altenburg 2007
                  ! The reason we don't use the IGIMF is that we are resolving the SF regions. 
        emi = 0.d0
        where ( (var(:,ind,1)*unit_d/1.66d-24*0.76) .ge. d0 ) 
          emi = var(:,ind,1)*unit_d*(dx*unit_l)**3/(t0*sqrt(d0/(var(:,ind,1)*unit_d/1.66d-24*0.76)))/1.989d33
        elsewhere
          emi = 0.d0
        endwhere

        tot_sfr = tot_sfr + sum(emi)
        where (emi .lt. 1.5d-5)
          emi = 0.d0
        endwhere
        tot_sfr2 = tot_sfr2 + sum(emi)


        where (emi .gt. 1.5d-5) 
          emi = SFRtoHA*emi
!         Standard IGIMF: emi = 1.d41/(1.27d-14/emi**5.7 + 0.18d0/emi**1.8 + 1.2d-4/emi**2.86 + 2.5d0/emi) 
                             ! to get a line density that is integrated back in amr2Intcube
!         Used to have 3.3e41 from Pflamm-Altenburg, now the default is the simple Kennicut value for Salpeter: 1.26e41
          emi = emi / (dx*unit_l/3.08e21)**2 ! This gives the surface brightness in units erg/s/kpc^2. 
        elsewhere
          emi = 0.d0
        end where

     CASE ('Ost_Ha') ! Osterbrock-based Ha estimate, following Valentin Perret's dissertation:
        emi = 0.d0

        if (present(temp_in)) then
          call get_spec(in,iT,cool_n,cool_T,cool_spec,ngrida,var(:,ind,1)*unit_d/1.66d-24*0.76,temp_in,spec)
          emi = 2.59e-13*(temp_in* ( var(:,ind,1)*unit_d*6.02e23 / &
                              (spec(:,1)+spec(:,2)+spec(:,3)+spec(:,4)+spec(:,5)+spec(:,6)) )/1.e4)**(-0.845)
        else
          emi=getTmu( var(:,ind,1), var(:,ind,iPT), .true. )
          call get_spec(in,iT,cool_n,cool_T,cool_spec,ngrida,var(:,ind,1)*unit_d/1.66d-24*0.76,emi,spec)
          emi = 2.59e-13*(emi*( var(:,ind,1)*unit_d*6.02e23 / &
                              (spec(:,1)+spec(:,2)+spec(:,3)+spec(:,4)+spec(:,5)+spec(:,6)) )/1.e4)**(-0.845)
                               ! assumption 1: here mu is 1: but fully neutral is 1.22, fully ionized is 0.588
              ! here roughly assuming electron density is roughly 20% hydrogen density. This is a function of T,rho --> see cool files. 
        endif

        hnu = 3.02679e-12

        emi = (dx*unit_l)**3*( spec(:,1) )**2 * emi * hnu
        tot_sfr2 = tot_sfr2 + sum(emi)/SFRtoHA
        tot_sfr = tot_sfr2
        emi = emi / (dx*unit_l/3.08e21)**2 ! This gives the surface brightness in units erg/s/kpc^2.
     CASE DEFAULT
        write(*,*) 'emType not found!!!'
        STOP
  END SELECT

  if(do_inverse) emi=1.d0/emi
  if(volMult) emi = emi * (dx * scale * unit_l)**ndim / 1.d20
  ! Need that division by 1.d20 because the output is single precision
  if(div_nH2) emi = emi / (var(:,ind,1)*unit_nH)**2

  !where (TK .ge. 1d5) emi=0.d0 !DEBUG JOKI
  !where (TK .ge. 1d5) var(:,ind,1)=0.d0 !DEBUG JOKI
  if(calc_TK .and. allocated(TK)) deallocate(TK)

END SUBROUTINE emissivity

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE prepareEmissivity(varName)

!-------------------------------------------------------------------------
  character(*)::varName
  integer::i
!-------------------------------------------------------------------------
  volMult=.false.
  div_nH2=.false.
  calc_TK=.false.
  do_inverse=.false.
  if(index(varName, 'ref') .eq. 1) then !---------------------------------
     emType='' ; emVar=1 ; defUnit=1.
  else if(index(varName, 'ramvar') .eq. 1) then !-------------------------
     read( varName(7:), '(i10)' ) i
     emType='oneVar' ; emVar=i ; defUnit=1d0
  else if(varname .eq. 'nH')       then !---------------------------------
     emType='oneVar' ; emVar=1 ; defUnit=unit_nH
  else if(varname .eq. 'SnH')       then !--------------------------------
     emType='oneVar' ; emVar=1 ; defUnit=unit_nH*scale*unit_l
  else if(varname .eq. 'negnH')       then !------------------------------
     emType='A' ; emVar=1 ; defUnit=unit_nH

  else if(index(varName, 'rhoKrumh') .eq. 1) then !-----------------------
     ! Krumholz/Davis experiment
     emType='oneVar' ; emVar=1 ; defUnit=1d-5
  else if(varName .eq. 'TrKrumh') then !----------------------------------
     emType='TrKrumh'
  else if(varName .eq. 'TgKrumh') then !----------------------------------
     emType='TgKrumh'
  else if(varName .eq. 'TratioKrumh') then !------------------------------
     emType='TratioKrumh'
  else if(varName .eq. 'fKrumh') then !-----------------------------------
     emType='fKrumh'
  else if(varName .eq. 'tauKrumh') then !---------------------------------
     emType='tauKrumh'

  else if(index(varName, 'rho') .eq. 1) then !----------------------------
     emType='oneVar' ; emVar=1 ; defUnit=unit_d
  else if(index(varName, 'Srho') .eq. 1) then !---------------------------
     emType='oneVar' ; emVar=1 ; defUnit=unit_d*scale*unit_l
  else if(index(varName, 'speed') .eq. 1) then !--------------------------
     emType='b' ; defUnit= unit_v
  else if(index(varName, 'vx') .eq. 1) then !-----------------------------
     emType='oneVar' ; emVar=2 ; defUnit=unit_v !; defSubtr=-4554021.4
  else if(index(varName, 'vy') .eq. 1) then !-----------------------------
     emType='oneVar' ; emVar=3 ; defUnit=unit_v !; defSubtr=-3532903.5
  else if(index(varName, 'vz') .eq. 1) then !-----------------------------
     emType='oneVar' ; emVar=4 ; defUnit=unit_v !; defSubtr=-9367328.1
  else if(varName .eq. 'P') then !----------------------------------------
     emType='oneVar' ; emVar=iPT ; defUnit=unit_P
  else if(varName .eq. 'Tmu') then !--------------------------------------
     emType='c' ; defUnit=unit_T2
  else if(varName .eq. 'T' .or. varName .eq. 'TInv') then !---------------
     emType='d' ; defUnit=unit_T2; calc_TK=.true.
     if(index(varName, 'Inv') .ne. 0) do_inverse=.true.
  else if(varName .eq. 'Trad') then !-------------------------------------
     emType='Trad'
  else if(varName .eq. 'Etherm') then !-----------------------------------
     emType='oneVar' ; emVar=iPT ; volMult=.true.
     defUnit = unit_d * unit_l**2 / unit_t**2 / (gamma-1.)
  else if(varName .eq. 'Ekin') then !-------------------------------------
     emType='e' ; defUnit=unit_T2 ; volMult=.true.
  else if(varName .eq. 'Z') then !----------------------------------------
     emType='oneVar' ; emVar=iPT+1
  else if(varName .eq. 'Mach') then !-------------------------------------
     emType='f' ; defUnit=unit_T2
  else if(varName .eq. 'xH2') then !--------------------------------------
     emType='h' ; emVar=iIons
  else if(varName .eq. 'xHI') then !--------------------------------------
     emType='g' ; emVar=iIons ; if(isH2) emType='oneVar'
  else if(varName .eq. 'xHII') then !-------------------------------------
     emType='oneVar' ; emVar=iIons ; if(isH2) emVar=iIons+1
  else if(varName .eq. 'xHeI') then !-------------------------------------
     emType='h' ; emVar=iIons+1 ; if(isH2) emVar=iIons+2
  else if(varName .eq. 'xHeII') then !------------------------------------
     emType='oneVar' ; emVar=iIons+1 ; if(isH2) emVar=iIons+2
  else if(varName .eq. 'xHeIII') then !-----------------------------------
     emType='oneVar' ; emVar=iIons+2 ; if(isH2) emVar=iIons+2
  else if(index(varName, 'Lya') .eq. 1) then !----------------------------
     ! Can be 'LyaEmi', 'LyaEmi_noPol', 'LyaLum', 'LyaLum_noPol'
     ! Can also contain 'coll' or 'rec'
     emType='L'
     calc_TK=.true.
     if(index(varName, 'coll')  .ne. 0) emType='C'
     if(index(varName, 'rec')   .ne. 0) emType='R'
     if(index(varName, 'grav')  .ne. 0) emType='l'
     if(index(varName, 'UV')    .ne. 0) emType='m'
     if(index(varName, 'noPol') .ne. 0) noPol=.true.
     if(index(varName, 'Lum')   .ne. 0) volMult=.true.
  else if(index(varName, 'Np')  .eq. 1) then !----------------------------
     read( varName(3:4), '(i1)' ) i
     emType='oneVar' ; emVar = nHvar+1+(i-1)*(ndim+1)
     defUnit = unit_Np/rt_c_cgs
  else if(index(varName, 'Fpx') .eq. 1) then !----------------------------
     read( varName(4:5), '(i1)' ) i
     emType='oneVar' ; emVar = nHvar+2+(i-1)*(ndim+1) ; defUnit = unit_Fp
  else if(varName .eq. 'Fpreduced' .or. varName .eq. 'Fpyreduced') then !-
     i=1 ! First photon group
     emType=varName ; emVar = nHvar+1+(i-1)*(ndim+1) ; defUnit = 1d0
  else if(index(varName, 'FpyPre')  .eq. 1) then !------------------------
     emType='oneVar' ; emVar = iIons+3 ; defUnit = unit_Fp
  else if(index(varName, 'Fpy')  .eq. 1) then !---------------------------
     read( varName(4:5), '(i1)' ) i
     emType='oneVar' ; emVar = nHvar+3+(i-1)*(ndim+1) ; defUnit = unit_Fp
  else if(index(varName, 'Fpz') .eq. 1) then !----------------------------
     read( varName(4:5), '(i1)' ) i
     emType='oneVar' ; emVar = nHvar+4+(i-1)*(ndim+1) ; defUnit = unit_Fp
  else if(index(varName, 'Fey')  .eq. 1) then !---------------------------
     read( varName(4:5), '(i1)' ) i
     emType='oneVar' ; emVar = nHvar+3+(i-1)*(ndim+1)
     defUnit = unit_Fp*group_egy(i)
  else if(index(varName, 'FeIR')  .eq. 1) then !--------------------------
     emType='FeIR'
  else if(index(varName, 'FeStream') .eq. 1) then !-----------------------
     emType='oneVar' ; emVar = nHvar+1 ; defUnit = unit_Fp*group_egy(1)
  else if(index(varName, 'FeTrap') .eq. 1) then !-------------------------
     emType='oneVar' ; emVar = iPrad ; defUnit = unit_P*c_cgs*3d0
  else if(index(varName, 'Fe')  .eq. 1) then !----------------------------
     read( varName(3:4), '(i1)' ) i
     emType='oneVar' ; emVar = nHvar+1+(i-1)*(ndim+1)
     defUnit = unit_Fp*group_egy(i)
  else if(index(varName, 'Fp') .eq. 1) then !-----------------------------
     read( varName(3:4), '(i1)' ) i
     emType='oneVar' ; emVar = nHvar+1+(i-1)*(ndim+1) ; defUnit = unit_Fp
     ! Temporary, to reproduce old maps:
     ! defUnit = unit_Fp * unit_v/rt_c_cgs
     ! Temporary, for old outputs:
     !defUnit = unit_Np * rt_c_cgs
     if(index(varName, 'Inv') .ne. 0) do_inverse=.true.
  else if(index(varName, 'Polytrope') .eq. 1) then !----------------------
     emType='i'
  else if(index(varName, 'Einv') .eq. 1) then !---------------------------
     emType='j'
  else if(varName .eq. 'SnH2') then !-------------------------------------
     emType='K' ; emVar=iIons+1 ; defUnit=unit_nH*scale*unit_l
  else if(varName .eq. 'nHI') then !--------------------------------------
     emType='k' ; emVar=iIons ; defUnit=unit_nH
     if(isH2) emType='p'
  else if(varName .eq. 'SnHI') then !-------------------------------------
     emType='k' ; emVar=iIons ; defUnit=unit_nH*scale*unit_l
     if(isH2) emType='p'
  else if(varName .eq. 'SnHII') then !-------------------------------------
     emType='p' ; emVar=iIons ; defUnit=unit_nH*scale*unit_l
     if(isH2) emVar=iIons+1
  else if(index(varName, 'PHrate') .eq. 1) then !-------------------------
     emType='n'
  else if(index(varName, 'Gamma') .eq. 1) then !--------------------------
     if(index(varName, 'Inv') .ne. 0) do_inverse=.true.
     emType='o'
  else if(varName .eq. 'nHII') then !-------------------------------------
     emType='p' ; emVar=iIons ; defUnit=unit_nH
  else if(varName .eq. 'Mom') then !--------------------------------------
     emType='q'
  else if(varName .eq. 'Momx') then !-------------------------------------
     emType='r' ; emVar=2
  else if(varName .eq. 'Momy') then !-------------------------------------
     emType='r' ; emVar=3
  else if(varName .eq. 'Momz') then !-------------------------------------
     emType='r' ; emVar=4
  else if(varName .eq. 'Pnt') then !--------------------------------------
     emType='oneVar'; emVar=iPrad ; defUnit=unit_P
  else if(varName .eq. 'FIRtrap') then !----------------------------------
     !!!emType='oneVar'; emVar=iIRtrap ; defUnit=unit_Fp
  else if(varName .eq. 'dCool') then !------------------------------------
     emType='oneVar'; emVar=7 ; defUnit=1d0
  else if(varName .eq. 'dSFR') then  !------------------------------------
     emType='t'; defUnit=1d0
  else if(varName .eq. 'sSFR') then  !------------------------------------
     emType='t'; defUnit=scale*unit_l/3.08d21
  else if(varName .eq. 'tauIR') then  !-----------------------------------
     emType='tauIR'; defUnit=1
  else if(varName .eq. 'StauIR') then  !----------------------------------
     emType='StauIR'; defUnit=1
  else if(varName .eq. 'cPyy') then  !-------------------------------------
     emType=varName; defUnit=unit_Fp*group_egy(1)
  else
     !write(*,*)'Invalid emissivity type'
     !write(*,*),'varName=',varName
     !STOP
     emtype=varName ! the emissivity function can perhaps deal with it
  endif

  if(index(varName, 'volMult') .ne. 0) volMult = .true.
  if(index(varName, 'divnH2')  .ne. 0) div_nH2=.true.

END SUBROUTINE prepareEmissivity

!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif


end subroutine title

subroutine read_cool(file_format,outdir,nout,cool_n,cool_T,cool_spec)

  implicit none

  integer :: file_format !1 = Nut_Old, 2=HAGN, Nut_New
  integer :: n1, n2,i

  logical :: debug=.false.

  real(kind=8),dimension(161) :: cool_n
  real(kind=8),dimension(101) :: cool_T
  real(kind=8),dimension(161,101,6) :: cool_spec

  character(len=128) :: cool_file,outdir
  character(len=5) :: nout

  cool_file=trim(outdir)//'/cooling_'//nout//'.out'

  open(42,file=cool_file,status='old',form='unformatted')

  read(42) n1,n2
  if ( (n1 .ne. 161) .or. (n2 .ne. 101) ) then
    write(6,*) "Error expect n_dens = 161, n_Temp = 101, got: ", n1, n2
    stop
  endif

  read(42) cool_n
  read(42) cool_T
  if (file_format .eq. 1) then
    read(42) ! Teq
  endif

  read(42) ! cool
  read(42) ! heat
  if (file_format .eq. 2) then
    read(42) ! cool_com
    read(42) ! heat_com
  endif
  read(42) ! metal

  read(42) ! cool_prime
  read(42) ! heat_prime
  if (file_format .eq. 2) then
    read(42) ! cool_com_prime
    read(42) ! heat_com_prime
  endif
  read(42) ! metal_prime

  read(42) ! mu

  read(42) cool_spec
  close(42)
  cool_spec = 10.d0**cool_spec
  if (debug) then
    write(6,*) "cool_n, cool_T, cool_spec"
    write(6,*) cool_n(1),cool_n(161)
    write(6,*) cool_T(1),cool_T(101)
    write(6,*) (cool_spec(1,1,i), cool_spec(161,101,i),i=1,6)
  endif
  ! species = {n_E, n_HI, n_HII, n_HeI, n_HeII, n_HeIII}

  write(6,*) "Cooling file successfully read in."

end subroutine read_cool


subroutine get_spec(in,iT,cool_n,cool_T2,cool_spec,ngrid,n_H,T2,spec)
  implicit none

  integer :: in, iT, ngrid
  double precision, dimension(in) :: cool_n
  double precision, dimension(iT) :: cool_T2
  double precision, dimension(in,iT,6) :: cool_spec
  double precision, dimension(ngrid)  :: n_H, T2
  double precision, dimension(ngrid,6) :: spec

  double precision :: dn, dT, logn, logT
  integer :: igrid, i, in_0, in_1, iT_0, iT_1

  logical :: debug=.false.

  dn = cool_n(2) - cool_n(1)
  dT = cool_T2(2) - cool_T2(1)

  if (debug) then
    write(6,*) "in, iT, cool_n, cool_T2, dn, dT, cool_spec, ngrid, n_H, T2"
    write(6,*) in
    write(6,*) iT
    write(6,*) cool_n(1),cool_n(161)
    write(6,*) cool_T2(1),cool_T2(101)
    write(6,*) dn
    write(6,*) dT
    write(6,*) (cool_spec(1,1,i), cool_spec(161,101,i),i=1,6)
    write(6,*) ngrid
    write(6,*) n_H(1), n_H(ngrid)
    write(6,*) T2(1), T2(ngrid)

  endif

  do igrid=1,ngrid
    logn = log10(n_H(igrid))
    logT = log10(T2(igrid))

    in_0 = int((logn - cool_n(1))/dn) + 1
    iT_0 = int((logT - cool_T2(1))/dT) + 1


    i=MAX(1,in_0-5)
    do while (cool_n(i) .lt. logn .and. i .lt. in)
      i=i+1
    enddo
    in_0 = i-1
    in_1 = i

    i=MAX(1,iT_0-5)
    do while (cool_T2(i) .lt. logT .and. i .lt. iT)
      i=i+1
    enddo
    iT_0 = i-1
    iT_1 = i

    spec(igrid,:) = 0.d0
    spec(igrid,:) = ( cool_spec(in_0,iT_0,:)*(cool_n(in_1)-logn)/dn + &
              cool_spec(in_1,iT_0,:)*(logn-cool_n(in_0))/dn ) * &
              (cool_T2(iT_1)-logT)/dT + &
            ( cool_spec(in_0,iT_1,:)*(cool_n(in_1)-logn)/dn + &
              cool_spec(in_1,iT_1,:)*(logn-cool_n(in_0))/dn ) * &
              (logT-cool_T2(iT_0))/dT
  enddo

end subroutine get_spec

End MODULE ramses_info
