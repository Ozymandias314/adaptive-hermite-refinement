MODULE constants


  !Nuno Loureiro
  !Constants and run specifications .
  !Last modified 07/26/06

  !08/02/2013 AVK : Changing this to read in input parameters
  ! Removed the parameter attribute from everywhere
  ! Have retained some arbitrary default values from the constants.f90 I was using
  ! It might be a good idea to change these to something more reasonable
  ! Added a read_parameters subroutine which reads in the input, and defines nkx, nky etc.
 

  implicit none
  save

  !********** CONSTANT DECLARATION *************
!BOX PARAMETERS:
  logical :: three_d=.true.
  real, parameter :: PI=3.1415926535
  real :: Lx=1.0	   !X box size in units of 2*pi
  real :: Ly=1.0	   !Y box size in units of 2*pi
  real :: Lz=1.0       !Z box size in units of 2*pi
  integer :: NLx=16  !number of grid points in X
  integer :: NLy=16  !number of grid points in Y
  integer :: NLz=16  !number of grid points in Z
  integer :: mm=2
  integer :: nn=2
  logical :: g_inc=.false.   !true included g's in the calculation
  integer :: ngtot = 0!number of Hermite moments; set to zero to run without g's
  integer :: gmin
  integer :: dim_vec=2   !dimension of the solution vector. 
                       !E.g., if g_inc=.false. dim_vec=2 (ne and Apar)
  integer :: NPE=1  ! number of processors
  integer :: NPEz=1
!time parameters:
  integer :: tmax=1000	!maximum number of iterations	
  real::init_aa0_fac=0.1  !multiplier on the SI operator   
  
  real :: CFL_frac=0.25
  real :: epsilon=1.e-10     !sets the tolerance for the pth iteration at each time step
!  real, parameter :: g_epsilon=1.e-3
  integer :: pmax=1  !maximum number of iterations per time step 
!FLR:      
  real :: rhoi = 1.e-8
  real :: rhos=1.e-8
  real :: de=1.e-8
  real :: rhos_de
!MHD:
  real :: small_rhoi = 1.e-6!if rhoi<small_rhoi, code evolves the RMHD eqs.
!ANJOR:
  logical :: anjor=.true.
  real :: notanj = 1.0
  real :: zcharge = 1.0
  real :: tite = 1.0
  real :: betai = 1.0
  integer :: plusorminus = 1
  integer :: inkx=1, inky=1, inkz = 1

  real :: lambda = 1.0
  real:: sigma = 1.0
  real :: rhos_diag
!Background electron temp gradient:
  real :: rhoe_LTe=0.0
  real :: kc0=0.0
!DIFFUSION:
  real :: nu_ei=0.00	  !collisions
  real :: res=0.0  !nu_ei*de**2
  real :: niu=0.0  !viscosity
  real :: hyper_coef=1.0 !0.25   !set to zero for no hyper-diffusion
  real :: hyper_coef_g=0.0 !0.25   !set to zero for no hyper-diffusion on the g's
  real :: hyperm_coef = 2.0  !0.25
  logical :: hyper_fixed=.false.  ! Fixes the hyper coeffs, if false, scales with d
  integer :: hyper_order=3  !this means it's k^(2*hyper_order)
  integer :: hyper_order_g=3  !this means it's k^(2*hyper_order_g) in the g eqs
  integer :: hyper_morder=3  !this means it's k^(2*hyper_order)
  real :: hyper_eta=1.0 !3.0e-8  
  real :: hyper_nu=1.0 !3.0e-8  
  real :: hyper_nu_g = 0.0  !0.25
  logical :: hyper_colls_fixed=.false.
  real :: hyper_colls=0.0

  !********** EQUILIBRIUM VALUES **************
  real :: A0=0.0	
  real :: PHI0=0.0   !not used
  real :: Leq=0.3
!  real :: perturb_A=0.0  !1.e-5*A0	   	 
!  real :: perturb_PHI=0.0  !1.0e-8  !1.e-5*PHI0
!  real :: g_perturb=0.0
  real :: perturb_amp=0.0
  character(len=4):: equilib_type="none" !AVK: Changed default equil to none
  character(len=4):: perturb_type="none"
  !*********************************
  !**** TURBULENT SOURCE PARAMETERS ******
  logical::turb=.true.
  !logical :: gturb=.false.
  real :: kfp1=1
  real :: kfp2=2
  real :: kfz1=1
  real :: kfz2=1
  real :: feps=0.0
  !***************************************
 
  !*** data saving stuff: ***************
  integer :: save_energyfiles=50      !saving frequency for energy files and cuts
  integer :: save_checkpoints=50      !saving frequency for fields/gfields checkpoints
  integer :: save_datafiles=50		!saving frequency for data files
  integer :: restart=0         !if 0 start from scratch; anything else start from file
  real :: rs_time=0.0  	  !if above is zero, so is this one
  logical,target,save :: use_hac_checkpoint = .true. !> MMH: if .true. and WANT_HAC then use HAC

  character(len=100) :: oldrun
  character(len=100) :: ppfile !AVK: to write out savetimes for postprocessing
  character(len=100) :: oldppfile !AVK: to read times from oldrun
  
  !integer :: nly_old=64
  character(len=100)::Apar_data="_Apar_data.dat"
  character(len=100)::gm_ttrace="_gm_timetrace.dat"
  character(len=100)::invariants="_invariants.dat"
  character(len=100)::hypercoeffs="_hypercoeffs.dat"
  character(len=100)::outflow="_outflow.dat"
  character(len=100)::resol="_resol.dat"
  character(len=100)::timestep="_timestep.dat"
  character(len=100)::grates="_grates.dat"
  !character(len=100)::energy_inj="_energy.dat"
  character(len=5)::cutx="cutx_"
  character(len=5)::cuty="cuty_"
  character(len=5)::cutz="cutz_"
  character(len=8)::energyk="energyk_"
  character(len=9)::energykp="energykp_"
  character(len=6)::zonal="zonal_" 
  character(len=9)::energy_gm="energygm_"

  character(len=100)::PATH=""
  character(len=100)::apar_in
  character(len=100)::gfields_in
  !logical::file_old=.true.     !if restart file is in old (formatted) format, set to .true.
 !this option is not working yet, set to .true.
  !integer,parameter::len1=len(apar_in)
  !integer,parameter::len11=len(gfields_in)
  integer,parameter::len3=len(PATH)
  !character(len=len1+len3)::restart_file1
  !character(len=len11+len3)::restart_file2
  character(len=100)::restart_file1
  character(len=100)::restart_file2
  character(len=8)::fieldsfile="_fields_"
  character(len=9)::gfieldsfile="_gfields_"
  character(len=7)::ktimefile="_error_"
  character(len=100)::kfieldsfile="kfieldsfile"

  integer :: NKx_parmax

  !***********Linear_Test !LMM
  logical::linear
  !********** !LMM
  real::amp_ratio
  logical::twomodes
  integer::kx1
  integer::ky1
  integer::kz1
  integer::kx2
  integer::ky2
  integer::kz2
  !************DEFINITIONS -- DO NOT CHANGE!!!!!!!!!!!!!!!!*******************************

  integer :: nkx, nky, nkz
  integer :: x_loc, y_loc, z_loc
  integer :: k_max, kperpmax
  integer :: nkx_par, nly_par, nlypar_old, nlz_par
  real :: dx, dy, dz
  real :: etaz = 0., etaz_g = 0. ! AVK: adding z viscosity term
  real:: scale
  
!  INTEGER, PARAMETER :: NKx=(NLx-1)/3 + 1
!  INTEGER, PARAMETER :: NKy=2*((NLy-1)/3)+1
!  INTEGER, PARAMETER :: x_loc=Nlx/2+1
!  INTEGER, PARAMETER :: y_loc=Nly/2+1
!  INTEGER, PARAMETER :: z_loc=Nlz/2+1
!  INTEGER, PARAMETER :: k_max=nkx-nkx/2   !maximum k after filtering
!  INTEGER, PARAMETER :: NKx_par = (NKx-1)/NPE+1
!  INTEGER, PARAMETER :: NLy_par = (NLy-1)/NPE+1	
!  INTEGER, PARAMETER :: NLypar_old = (NLy_old-1)/NPE+1	
!  INTEGER, PARAMETER :: nlz_par = (nlz-1)/NPEz+1
!  real,parameter::dx=lx/(nkx*1.)
!  real,parameter::dy=ly/(nky*1.)
!  real,parameter::dz=lz/(nlz*1.)
!  REAL, PARAMETER :: scale=1./(NLx*NLy) ! scale factor for FFTs

 !**** ANTENNA DRIVE ********************
  real :: amplitude=0.0
  real :: omega0=0.9
  integer :: kpar0=1
  integer :: kperp0=1  !this is ky
  integer :: j1, j2
  real :: facpm
  !**************************************

  ! ZINTEGRATION METHOD ****************
  character(len=7)::zmethod='cormack'  !alternative is 'highord'

 !******************************************************************************************
contains


!*********************************************************
  subroutine set_anjor(anjor, rhos_de, gmin)
    implicit none
    logical::anjor
    integer::gmin
    real::rhos_de

    if(anjor) then
       rhos_de=1.0/sqrt(2.0)  !the ratio of rho_s/d_e
       rhos_diag = 1.0 !rhos used in diagnostics
       gmin=0
       notanj=0.0
       sigma = 1.0 + tite/zcharge + 1./betai + sqrt((1. + tite/zcharge)**2 + &
            1./betai**2)
       lambda = - tite/zcharge + 1./betai + &
            plusorminus*sqrt((1.+tite/zcharge)**2 + 1./betai**2)
       
       if(plusorminus == 1) then 
          facpm = (1./((zcharge/tite)*(betai/2.) +1. + betai)) * &
               (1.+ (1./sigma)*(1.+zcharge/tite))
       else
          facpm = (1./((zcharge/tite)*(betai/2.) +1. + betai)) * &
               (1. + 2.*tite/(sigma*zcharge*betai))
       endif
    else
       gmin=2
       rhos_de=rhos/de
       rhos_diag = rhos !AVK: added rhos for diagnostics
       notanj=1.0
       !NFL: lambda needs to be given a value even if anjor=false because it appears in funcgm
       !AVK: Changed this to lambda = 1.0 by default in constants.f90 (look above)
       !lambda=1.0
    end if
  end subroutine set_anjor


!*********************************************************
  subroutine read_parameters(ipfile)

    implicit none

    character(len=100), intent(in) :: ipfile

    integer::ierr

    namelist /box_parameters/ lx, ly, lz, nlx, nly, nlz, &
         & g_inc, ngtot, npe, npez, three_d, mm, nn
    
    namelist /time_parameters/ tmax, init_aa0_fac, cfl_frac, epsilon, &
         & pmax
    
    namelist /flr/ rhoi, rhos, de
    
    namelist /mhd/ small_rhoi
    
    namelist /anjor_param/ anjor, inkx, inky, inkz, tite, zcharge, betai
    
    namelist /elec_temp_grad/ rhoe_LTe
    
    namelist /diffusion/ nu_ei, res, niu, hyper_coef, hyper_coef_g, &
         hyperm_coef, hyper_fixed, hyper_order, hyper_order_g, hyper_morder,&
         hyper_eta, hyper_nu, hyper_nu_g, etaz, etaz_g, hyper_colls_fixed, &
         hyper_colls
    
    namelist /equil/ a0, phi0, leq, perturb_amp, equilib_type, &
         perturb_type
    
    namelist /turb_src/ turb, kfp1, kfp2, kfz1, kfz2, feps
    
    namelist /data_sav/ save_energyfiles, save_checkpoints, save_datafiles, restart, rs_time, PATH, oldrun, use_hac_checkpoint !> MMH added use_hac_checkpoint
    !namelist /data_sav/ save_energyfiles, save_checkpoints, save_datafiles, restart, rs_time, PATH, oldrun
    
    namelist /antenna/ amplitude, omega0, kpar0, kperp0

    namelist /zintegration/ zmethod

    namelist /lineartest/ linear 

    namelist /seeding/ amp_ratio, twomodes, kx1, ky1, kz1, kx2, ky2, kz2  
    
    open(unit=10, file=trim(ipfile), status='old')
    read (10, nml=box_parameters,iostat=ierr)
    if(ierr /= 0) then
       write(*,*) "Reading box_parameters failed"
       write(*,*) "lx, ly, lz, nlx, nly, nlz, g_inc, ngtot, npe, npez, three_d =",&
            lx, ly, lz, nlx, nly, nlz, g_inc, ngtot, npe, npez, three_d
    end if
    read(10, nml=time_parameters,iostat=ierr)
    if(ierr /= 0) then
       write(*,*) "Reading time_parameters failed"
       write(*,*) "tmax, init_aa0_fac, cfl_frac, epsilon, pmax =",&
            tmax, init_aa0_fac, cfl_frac, epsilon, pmax
    end if
    read(10, nml=flr,iostat=ierr)
    if(ierr /= 0) then
       write(*,*) "Reading flr failed"
       write(*,*) "rhoi, rhos, de =", rhoi, rhos, de
    end if
    read(10, nml=mhd,iostat=ierr)
    if(ierr /= 0) then
       write(*,*) "Reading mhd failed"
       write(*,*) "small_rhoi =", small_rhoi
    end if
    read(10, nml=anjor_param,iostat=ierr)
    if(ierr /= 0) then
       write(*,*) "Reading anjor_param failed"
       write(*,*) "anjor, inkx, inky, inkz, tite, zcharge, betai =", anjor, inkx, inky, inkz, tite, zcharge, betai
    end if
    read(10, nml=elec_temp_grad,iostat=ierr)
    if(ierr /= 0) then
       write(*,*) "Reading elec_temp_grad failed"
       write(*,*) "rhoe_lte =", rhoe_lte
    end if
    read(10, nml=diffusion,iostat=ierr)
    if(ierr /= 0) then
       write(*,*) "Reading diffusion failed"
       write(*,*) "nu_ei, res, niu =", nu_ei, res, niu
       write(*,*) "hyper_coef, hyper_coef_g, hyperm_coef =", hyper_coef, hyper_coef_g, hyperm_coef
       write(*,*) "hyper_fixed, hyper_order, hyper_order_g =", hyper_fixed, hyper_order, hyper_order_g
       write(*,*) "hyper_morder, hyper_eta, hyper_nu, hyper_nu_g = ", hyper_morder, hyper_eta, hyper_nu, hyper_nu_g
    end if
    read(10, nml=equil,iostat=ierr)
    if(ierr /= 0) then
       write(*,*) "Reading equil failed"
       write(*,*) "a0, phi0, leq, perturb_amp, equilib_type, perturb_type =",&
            a0, phi0, leq, perturb_amp, equilib_type, perturb_type
    end if
    read(10, nml=turb_src,iostat=ierr)
    if(ierr /= 0) then
       write(*,*) "Reading turb_src failed"
       write(*,*) "turb, kfp1, kfp2, kfz1, kfz2, feps =",&
            turb, kfp1, kfp2, kfz1, kfz2, feps
    end if
    read(10, nml=data_sav,iostat=ierr)
    if(ierr /= 0) then
       write(*,*) "Reading data_sav failed"
       write(*,*) "save_energyfiles, save_checkpoints, save_datafiles, restart, rs_time, path, oldrun =",&
            save_energyfiles, save_checkpoints, save_datafiles, restart, rs_time, path, oldrun
    end if
    read(10, nml=antenna,iostat=ierr)
    if(ierr /= 0) then
       write(*,*) "Reading antenna failed"
       write(*,*) "amplitude, omega0, kpar0, kperp0 =",&
            amplitude, omega0, kpar0, kperp0
    end if
    read(10, nml=zintegration,iostat=ierr)
    if(ierr /= 0) then
       write(*,*) "Reading zintegration failed"
       write(*,*) "zmethod =", zmethod
    end if
    read(10, nml=lineartest, iostat=ierr)
    if(ierr /= 0) then
       write(*,*) "Reading lineartest failed"
       write(*,*) "lineartest =", linear
    end if
    read(10, nml=seeding, iostat=ierr)
    if(ierr /= 0) then
       write(*,*) "Reading seeding failed"
       write(*,*) "amp_ratio, twomodes, kx1, ky1, kz1, kx2, ky2, kz2 =", amp_ratio, twomodes, kx1, ky1, kz1, kx2, ky2, kz2
    end if
    close(10)
    
    lx = 2.*pi*lx
    ly = 2.*pi*ly
    lz = 2.*pi*lz

!AVK: dim_vec should always be 2 + (ngtot-gmin+1) if g_inc is true
   if(g_inc) then
     dim_vec = ngtot-gmin+3
   else
     dim_vec = 2
   end if

    
    !************DEFINITIONS -- DO NOT CHANGE!!!!!!!!!!!!!!!!*******************************

!    nkx=(nlx-1)/3 + 1
    nkx=nlx/2+1
!    nky=2*((nly-1)/3)+1
!    nky=2*((nly-1)/2)+1
    nky=nly
    nkz=2*((nlz-1)/3)+1
    x_loc=nlx/2+1
    y_loc=nly/2+1
    z_loc=nlz/2+1
    k_max=nkx-nkx/2   !maximum k after filtering
    kperpmax=sqrt((nkx*1.0)**2+(nky/2.)**2)
    nkx_par = (nkx-1)/npe+1
    nly_par = (nly-1)/npe+1	
    !nlypar_old = (nly_old-1)/npe+1	
    nlz_par = (nlz-1)/npez+1
    dx=lx/(nkx*1.)
    dy=ly/(nky*1.)
    dz=lz/(nlz*1.)
!    etaz = etaz*3./(nlz*1.) !AVK z viscosity will act on kz .ge. nlz/3 == dealiasing?
!    etaz_g = etaz_g*3./(nlz*1.) !AVK z viscosity will act on kz .ge. nlz/3 == dealiasing?
    etaz = etaz*dz**2
    etaz_g = etaz_g*dz**2
    scale=1./(nlx*nly) ! scale factor for FFTs

    j1 = kperp0*ly/lx+1
    j2 = nky-j1+2

  end subroutine read_parameters

END module constants

