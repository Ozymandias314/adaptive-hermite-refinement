program REGK
!*************************************************************************
!VIRIATO, 3D VERSION
!Original serial code by Nuno Loureiro.
!Parallelization by Bill Dorland and N. Loureiro.
!
!Reduced GK (KREHM model) + slow mode eqs (see tome) + RMHD eqs
!
!
!Last changes:
!NFL, 15/08/2012: moved the z-step to a subroutine (z_step)
!NFL, 05/08/2012: added background temperature gradient terms --- still checking
!
!
!08/02/2013 AVK
! Added input file capability. If the input file is runname.in, the run command is viriato
! runname
! Made all arrays allocatable which are now allocated after reading in run parameters
!
! forchk -l fck_report -define perf,gasca2d -r8 -cpp -ff REGK.F90 mpi_mod.F90 mp_mpi_r8.f90 constants.f90 fft_work_fftw.f90 grid.f90 redistribute_mpi.f90 transforms.F90 diag.F90 forcing.f90 functions.f90 initialize.F90 brackets.F90 stepping.F90 fluxes.f90 aux.F90 $FCKDIR/share/forcheck/MPI.flb
!
! Jan-Mar 2016: Michele Martone (MMH):
! Introduced optional Parallel ! I/O using the HLST ADIOS-Checkpoint (HAC) module.
! Its preprocessor conditional is WITH_HAC, file hlst_adios_checkpoint.F90
!******************************************************************

  !TTR
  !. FPP macro to simplify RZG Perflib calls (performance-counter)
# include "perfmacro.h"

  use constants,  only: set_anjor, read_parameters, nlx, nly, nly_par, nlz_par, nkx, nkx_par, nky, &
       &                g_inc, gmin, ngtot,  Lx, Ly, Lz, inkx, inky, inkz, x_loc, pi, &
       &                Apar_data, gm_ttrace, invariants, hypercoeffs,outflow, resol, timestep, &
       &                ppfile, restart, oldppfile, oldrun, rs_time, apar_in, gfields_in, &
       &                fieldsfile, gfieldsfile, anjor, PATH, &
       &                rhos, rhos_de, rhoi, small_rhoi, rhoe_LTe, CFL_frac, &
       &                perturb_amp, npe, npez, dx, dy, dz, de, restart_file1, restart_file2, &
       &                hyper_coef, hyper_coef_g, hyper_fixed, hyper_order, hyper_order_g, &
       &                hyper_morder, hyper_colls_fixed, hyper_colls, init_aa0_fac, &
       &                hyper_nu, hyper_nu_g, hyper_eta, hyper_coef, hyperm_coef, tmax, turb, feps, &
       &                kfp1, kfp2, kfz1, kfz2, kperpmax, niu, pmax, epsilon, three_d, zmethod, &
       &                save_datafiles, save_energyfiles, save_checkpoints, nu_ei, j1, &
       &                nlz, use_hac_checkpoint
  use mp
  use transforms, only: init_transforms, FFT2d_direct, FFT2d_inv
  use grid,       only: kperp, xx, yy, zz, gama0, init_grid
  use diag,       only: proc0_reads_to_all, proc0_reads_to_all_nogs, flows, b_field, predatasave, &
       &                datasavetest, kdatasave3d, test_of_test, cuts, gm_spectrum, k_energy, k_energy_pert,injected, &
       &                diagnostics, gamma_ky, zsteploss, gm_timetrace, zonalflows,  &
       &                Convol, Convol2, dataloadtest, zderivgk2
  use forcing,    only: force, uniran
  use Initialize, only: init_X_point, init_perturb, equilibrium
  use Functions,  only: exp_nu, exp_eta, exp_ng, &
       &                nu_func1, nu_func2, nu_func3, nu_func4, nu_func5, nu_func6
  use Brackets,   only: funcne_i, funcAkpar_i, funcg2, funcgm, func_lastg, bracket_3
  use Stepping,   only: z_step, z_step2
  use Fluxes,     only: rightflux, leftflux, pass_z_array_left_flux, pass_z_array_right_flux
  use Aux,        only: dtnext, increase_fac, hyper_diff, z_diffusion, omegakaw, &
       &                PHI_POT, SEMI_IMP_OP, resolution_check, double_resolution, &
       &                comp_eigen, inv, calc_leftover, closure
#if WITH_HAC
  use diag,       only: dataloadtest, hv_tstime, hv_lbl, hv_vv !< for hac
  use mpi,        only: MPI_COMM_WORLD !< for hac
  use hlst_adios_checkpoint, only: hac_init,hac_stats,hac_exit,HAC_QUIET,HAC_VERBOSE !< for hac
#endif


  implicit none

  !******* VARIABLE DECLARATION **********

  integer :: i, j, jj, k, t, p, z, jglob, y_loc_loc, proc_X, pos_isl, p_layer, n, ng, ip
!  real:: elapsed_time
  
  complex, allocatable, dimension(:,:,:) :: nek,neKnew,guess,nek_star,dummy,dummy2,dzgk2,braakpargk2 !the last element of this allocation is for the zderiv !helicity
  complex, allocatable, dimension(:,:,:) :: PHIK,phik_star, AKpar,PHIKnew, SI_oper,rel_error
  complex, allocatable, dimension(:,:,:) :: AKparnew, Akpar_star!, dadtK
  complex, allocatable, dimension(:,:,:) :: ueKpar, ueKparnew,uekpar_star,braakparuekpar
  complex, allocatable, dimension(:,:,:) :: AKpar_eq_double_prime,bkxx,bkyy
  complex, allocatable, dimension(:,:,:) :: AKpar_eq, ueKpar_eq, neK_perturb,akpar_perturb
  complex, allocatable, dimension(:,:,:) :: Fne_old,FApar_old,Fnepred,FAparpred,Fapar_star
  !AVK: not being used anywhere
  !complex, allocatable, dimension(:,:,:) :: Fg2_old,Fg2_pred
  complex, allocatable, dimension(:,:,:,:) :: gk, fgm_old, Fgm_pred, gk_star, gknew
  complex, allocatable, dimension(:,:,:) :: f_lastg_old, f_lastg, f_firstg_old

  real, allocatable, dimension(:,:,:) :: ne, ne_buff, vex, vey, vrhosx, vrhosy, temp, dummy_real, dummy_real2, ne_real
  real, allocatable, dimension(:,:,:) :: phi, Apar, Apar_buff, real_error
  real, allocatable, dimension(:,:,:) :: Epar, uepar, bx, by, bperp, uxavg, uyavg, ux_temp, uy_temp
  real, allocatable, dimension(:,:,:) :: Apar_eq_double_prime, Apar_eq, uepar_eq, phi_eq
  real, allocatable, dimension(:,:,:) :: apar_perturb, phi_perturb
  real, allocatable, dimension(:,:,:) :: Dxphi, Dyphi, value_phix, value_phiy, Dxphi_star, Dyphi_star
  real, allocatable, dimension(:,:,:) :: Dxgm, Dygm
  real, allocatable, dimension(:,:,:) :: Dxne, Dyne, value_nex, value_ney, Dxne_star, Dyne_star
  real, allocatable, dimension(:,:,:) :: Dxapar, Dyapar
  real, allocatable, dimension(:,:,:) :: Dxapar_star, Dyapar_star
  real, allocatable, dimension(:,:,:) :: Dxuepar, Dyuepar
  real, allocatable, dimension(:,:,:) :: Dxuepar_star, Dyuepar_star
  real, allocatable, dimension(:,:,:,:) :: Dxg, Dyg, Dxg_star, Dyg_star, value_gx, value_gy, g, gdummy
  !AVK: not being used anywhere
  !real, allocatable, dimension(:) :: phizero,phizero_star, phizero_new
  real, allocatable, dimension(:,:,:) :: dumq1, dumq2, dumq3, dumq4
  real, allocatable,dimension(:,:,:,:) :: dumqg
  real:: vxmax, vymax, CFL_flow, savetime
  real:: width,W2,L_sheet,LCS2 !,L_sheet_alex,deltaprime
  real:: outflow_max, Bup_max, Bup_jemella, delta, delta_x,delta_y,a_X, u_X, x !,alex_width
  real::bxmax,bymax,bperp_max
  real::omega_kaw,rel_err,RELATIVE_ERROR,old_error
  real::res2,niu2,nu_g,hyper_nuei  !hyper-resisitvity and viscosity
  real::zdiffcoeff
  real::low,energy,old_energy,denergy_dt,av_denergy_dt,total,lost
  real::aa0,dti,x1,x2,x3,ek_u,ek_b, ek_e, hk_a, hk_t, ik_t, ik_a, ek_n, hk_aa, hk_ta, ik_ta, ik_aa,hk_1,hk_2,chk, sumleftover !,new_dti !,N_aa0
  real,allocatable,dimension(:)::xg
  real::epar_avg,vrms,vrms2,inj,energyzloss,energyzloss1,energyzloss2,ux,uy,totalgk,totalgknew
  integer :: howlong, p_iter, p_count, number
  character(len=100)::r1,r2,r3,r4
  integer:: keepatit=0
  logical::repeat,noinc,divergent  !,g_continue
  logical,save::first=.true.
  logical::flagout=.false.
  real:: kfz,fac  !, g_kfz
  real:: time1, time2  ! to see how long the code takes inside main loop
  logical::debugging=.false. !LMM - Flag for debugging

  !AVK: added to read in input file
  character(len=100) :: runname, ipfile
  character(len=20) :: rs_time_str
  character(len=10):: time, date
#if WITH_HAC
  !> MMH: Variables to get HAC-related statistics
  logical,save :: wvhac = .true. !< Want Verbose HLST-ADIOS-Checkpoint. Turn this on/off at your option.
  integer,save :: rdi, rd=1 !< Repeat datasavetest and loop index
  real(kind=8) :: ogwb, ogrb !< Overall global written bytes
  real(kind=8) :: cstep_et, save_et, lsave_et !< Computing Step/Save/Loop Save elapsed time (in s)
#endif

!**************************************************************************************
!
!		MAIN SECTION OF THE CODE STARTS HERE
!
!**************************************************************************************

  !. Perflib initialization
  Myperfinit

  !. Perflib time measurement
  Myperfon('Main')

  !. Code initialization
  Myperfon('CodeInit')
  ! AVK: read ipfile and allocate arrays > 
  call getarg(1, runname)
  ipfile = trim(runname)//".in"
  call read_parameters(ipfile)
  ! Allocate arrays
  allocate(dzgk2(nky, nkx_par, nlz_par)) !zderivative
  allocate(braakpargk2(nky,nkx_par,nlz_par)) !helicity
  allocate(nek(nky, nkx_par, nlz_par))
  allocate(ne_buff(nlx, nly_par, nlz_par))
  allocate(neknew(nky, nkx_par, nlz_par))
  allocate(guess(nky, nkx_par, nlz_par))
  allocate(nek_star(nky, nkx_par, nlz_par))
  allocate(dummy(nky, nkx_par, nlz_par))
  allocate(dummy2(nky, nkx_par, nlz_par))
  allocate(phik(nky, nkx_par, nlz_par))
  allocate(phik_star(nky, nkx_par, nlz_par))
  allocate(akpar(nky, nkx_par, nlz_par))
  allocate(phiknew(nky, nkx_par, nlz_par))
  allocate(si_oper(nky, nkx_par, nlz_par))
  allocate(rel_error(nky, nkx_par, nlz_par))
  allocate(akparnew(nky, nkx_par, nlz_par))
  !allocate(dadtk(nky, nkx_par, nlz_par))
  allocate(akpar_star(nky, nkx_par, nlz_par))
  allocate(uekpar(nky, nkx_par, nlz_par))
  allocate(uekparnew(nky, nkx_par, nlz_par))
  allocate(uekpar_star(nky, nkx_par, nlz_par))
  allocate(braakparuekpar(nky, nkx_par, nlz_par))
  allocate(akpar_eq_double_prime(nky, nkx_par, nlz_par))
  allocate(bkxx(nky, nkx_par, nlz_par))
  allocate(bkyy(nky, nkx_par, nlz_par))
  allocate(akpar_eq(nky, nkx_par, nlz_par))
  allocate(uekpar_eq(nky, nkx_par, nlz_par))
  allocate(nek_perturb(nky, nkx_par, nlz_par))
  allocate(akpar_perturb(nky, nkx_par, nlz_par))
  allocate(fne_old(nky, nkx_par, nlz_par))
  allocate(fapar_old(nky, nkx_par, nlz_par))
  allocate(fnepred(nky, nkx_par, nlz_par))
  allocate(faparpred(nky, nkx_par, nlz_par))
  allocate(fapar_star(nky, nkx_par, nlz_par))
  !allocate(fg2_old(nky, nkx_par, nlz_par))
  !allocate(fg2_pred(nky, nkx_par, nlz_par))
  
  allocate(                  ne(nlx, nly_par, nlz_par))
  allocate(                 vex(nlx, nly_par, nlz_par))
  allocate(                 vey(nlx, nly_par, nlz_par))
  allocate(              vrhosx(nlx, nly_par, nlz_par))
  allocate(              vrhosy(nlx, nly_par, nlz_par))
  allocate(                temp(nlx, nly_par, nlz_par))
  allocate(          dummy_real(nlx, nly_par, nlz_par))
  allocate(         dummy_real2(nlx, nly_par, nlz_par))
  allocate(             ne_real(nlx, nly_par, nlz_par))
  allocate(                 phi(nlx, nly_par, nlz_par))
  allocate(                apar(nlx, nly_par, nlz_par))
  allocate(           apar_buff(nlx, nly_par, nlz_par))
  allocate(          real_error(nlx, nly_par, nlz_par))
  allocate(                epar(nlx, nly_par, nlz_par))
  allocate(               uepar(nlx, nly_par, nlz_par))
  allocate(                  bx(nlx, nly_par, nlz_par))
  allocate(                  by(nlx, nly_par, nlz_par))
  allocate(               bperp(nlx, nly_par, nlz_par))
  allocate(               uxavg(nlx, nly_par, nlz_par))
  allocate(               uyavg(nlx, nly_par, nlz_par))
  allocate(             ux_temp(nlx, nly_par, nlz_par))
  allocate(             uy_temp(nlx, nly_par, nlz_par))
  allocate(apar_eq_double_prime(nlx, nly_par, nlz_par))
  allocate(             apar_eq(nlx, nly_par, nlz_par))
  allocate(            uepar_eq(nlx, nly_par, nlz_par))
  allocate(              phi_eq(nlx, nly_par, nlz_par))
  allocate(        apar_perturb(nlx, nly_par, nlz_par))
  allocate(         phi_perturb(nlx, nly_par, nlz_par))
  allocate(               dxphi(nlx, nly_par, nlz_par))
  allocate(               dyphi(nlx, nly_par, nlz_par))
  allocate(          value_phix(nlx, nly_par, nlz_par))
  allocate(          value_phiy(nlx, nly_par, nlz_par))
  allocate(          dxphi_star(nlx, nly_par, nlz_par))
  allocate(          dyphi_star(nlx, nly_par, nlz_par))
  allocate(                dxgm(nlx, nly_par, nlz_par))
  allocate(                dygm(nlx, nly_par, nlz_par))
                      
  allocate(                Dxne(nlx, nly_par, nlz_par))
  allocate(                Dyne(nlx, nly_par, nlz_par))
  allocate(           value_nex(nlx, nly_par, nlz_par))
  allocate(           value_ney(nlx, nly_par, nlz_par))
  allocate(           Dxne_star(nlx, nly_par, nlz_par))
  allocate(           Dyne_star(nlx, nly_par, nlz_par))
  allocate(              Dxapar(nlx, nly_par, nlz_par))
  allocate(              Dyapar(nlx, nly_par, nlz_par))
  allocate(         Dxapar_star(nlx, nly_par, nlz_par))
  allocate(         Dyapar_star(nlx, nly_par, nlz_par))
  allocate(             Dxuepar(nlx, nly_par, nlz_par))
  allocate(             Dyuepar(nlx, nly_par, nlz_par))
  allocate(        Dxuepar_star(nlx, nly_par, nlz_par))
  allocate(        Dyuepar_star(nlx, nly_par, nlz_par))
  
  !allocate(phizero(nlz))
  !allocate(phizero_star(nlz))
  !allocate(phizero_new(nlz))
  
  ! Create diagnostic filenames
  Apar_data = trim(runname)//trim(Apar_data)
  gm_ttrace = trim(runname)//trim(gm_ttrace)
  invariants = trim(runname)//trim(invariants)
  hypercoeffs = trim(runname)//trim(hypercoeffs)
  outflow = trim(runname)//trim(outflow)
  resol = trim(runname)//trim(resol)
  timestep = trim(runname)//trim(timestep)
  ppfile = trim(runname)//trim(".times") !AVK: To write out savetimes
  !energy_inj = trim(runname)//trim(energy_inj)
  rs_time_str=''! MMH
  
  if(restart/=0) then
     if(debugging) then
     print*, 'restart started'
     endif
     !AVK: automatically reads from oldrun.times to find restart time
     oldppfile = trim(oldrun)//trim(".times")
     open(unit = 1728, file=trim(oldppfile))
     i = 0
     do while( i==0)
        read(1728, *, iostat=i) rs_time
     end do
     close(1728)

     if(rs_time .lt. 10) then
        write(rs_time_str, '(f5.3)') rs_time
     end if
     if(rs_time .ge. 10 .and. rs_time .lt. 100) then
        write(rs_time_str, '(f6.3)') rs_time
     end if
     if(rs_time .ge. 100 .and. rs_time .lt. 1000) then
        write(rs_time_str, '(f7.3)') rs_time
     end if
     if(rs_time .ge. 1000 .and. rs_time .lt. 10000) then
        write(rs_time_str, '(f8.3)') rs_time
     end if
     if(rs_time .ge. 10000 .and. rs_time .lt. 100000) then
        write(rs_time_str, '(f9.3)') rs_time
     end if
     
     apar_in = trim(oldrun)//trim(fieldsfile)//trim(rs_time_str)//".dat"
     gfields_in = trim(oldrun)//trim(gfieldsfile)//trim(rs_time_str)//".dat"
     if(debugging) then
     print*, 'gfields_in'
     endif
  end if
  
  ! <AVK
  call set_anjor(anjor,rhos_de,gmin)
  open (unit = 1729, file=trim(ppfile)) !AVK: Savetimes will be written to ppfile
   open (unit=1408, File=trim(hypercoeffs)) !LMM
        write(1408, '(a12,a12,a12,a12)') 'res2', &
             &'niu2', 'nu_g', 'hyper_nuei'

  width=0.0
  w2=0.0
  vrms=0.0
  vrms2=0.0
  
  low=0.92
  call init_mp
  call init_grid
  call init_transforms
  call init_X_point(proc_X, y_loc_loc)
  
  savetime=rs_time
  
  if(g_inc) then
     allocate (gk          (nky, nkx_par, nlz_par, gmin:ngtot))
     allocate (gk_star     (nky, nkx_par, nlz_par, gmin:ngtot))
     allocate (gknew       (nky, nkx_par, nlz_par, gmin:ngtot))
     allocate (Fgm_old     (nky, nkx_par, nlz_par, gmin:ngtot))
     allocate (Fgm_pred    (nky, nkx_par, nlz_par, gmin:ngtot))
     allocate (f_lastg     (nky, nkx_par, nlz_par))
     allocate (f_lastg_old (nky, nkx_par, nlz_par))
     allocate (f_firstg_old(nky, nkx_par, nlz_par))
     !     allocate (array_zparbc(nky,nkx_par,1:ngtot-gmin+5))
     allocate (       g(nlx, nly_par, nlz_par, gmin:ngtot))
     allocate (  gdummy(nlx, nly_par, nlz_par, gmin:ngtot))
     allocate (value_gx(nlx, nly_par, nlz_par, gmin:ngtot))
     allocate (value_gy(nlx, nly_par, nlz_par, gmin:ngtot))
     allocate (     Dxg(nlx, nly_par, nlz_par, gmin:ngtot))
     allocate (     Dyg(nlx, nly_par, nlz_par, gmin:ngtot))
     allocate (Dxg_star(nlx, nly_par, nlz_par, gmin:ngtot))
     allocate (Dyg_star(nlx, nly_par, nlz_par, gmin:ngtot))
  else
     allocate (gk(1,1,1,1))
     allocate (gk_star(1,1,1,1))
     allocate (gknew(1,1,1,1))
     allocate (Fgm_old(1,1,1,1))
     allocate (F_lastg_old(1,1,1))
     allocate (F_lastg(1,1,1))
     allocate (Fgm_pred(1,1,1,1))
     !     allocate (array_zparbc(nky,nkx_par,4))
     
     allocate (g(1,1,1,1))
     allocate (gdummy(1,1,1,1))
     allocate (value_gx(1,1,1,1))
     allocate (value_gy(1,1,1,1))
     allocate (Dxg(1,1,1,1))
     allocate (Dyg(1,1,1,1))
     allocate (Dxg_star(1,1,1,1))
     allocate (Dyg_star(1,1,1,1))
  end if
  gk(:,:,:,:)=0.0
  gk_star(:,:,:,:)=0.0
  gknew(:,:,:,:)=0.0
  Fgm_old(:,:,:,:)=0.0
  F_lastg_old=0.0
  f_lastg=0.0
  Fgm_pred(:,:,:,:)=0.0
  
  g(:, :, :, :) = 0.0
  value_gx(:, :, :, :) = 0.0
  value_gy(:, :, :, :) = 0.0
  Dxg(:, :, :, :) = 0.0
  Dyg(:, :, :, :) = 0.0
  Dxg_star(:, :, :, :) = 0.0
  Dyg_star(:, :, :, :) = 0.0
  !  end if
  dummy_real = 0.0
  gdummy = 0.0
  if(debugging) then
  print*, 'allocation completed'
  endif

!
#if WITH_HAC
  if_hacinit: if( use_hac_checkpoint ) then
    !> MMH: Benchmark mode.
    !! Assuming production will use large tmax this turns on repeated checkpointing.
    if ( tmax .lt. 5  ) then
      if ( proc0 ) write(*,'(a,i0,a)') hv_lbl//'Given the small tmax (',tmax,') will run in benchmark mode.'
      if ( proc0 ) write(*,'(a,i0,a)') hv_lbl//'Will repeat checkpointing ',rd,' times.'
      wvhac=.true.
      rd=5
    end if
  
    if ( wvhac ) then
      if( proc0 ) write(*,'(a)') hv_lbl//'Checkpoint/restart will be more verbose.'
      hv_vv = HAC_VERBOSE
    end if
    call hac_init (MPI_COMM_WORLD, hv_vv )
  end if if_hacinit
#else
  if_hacwhac: if( use_hac_checkpoint ) then
    if(proc0) write(*,'(a)') "ERROR: User requested use_hac_checkpoint=.true. without having activated HAC !"
    stop
  end if if_hacwhac
#endif
!
!  open (unit=2020, File='debugregi')
 if (debugging) then
 print*, 'checkpoint a'
 endif
  if (proc0) then
     !   call timer(ti)
     open (unit=343, file= 'exit')
     write (343,'(i1)') keepatit
     close (unit=343)
  end if
  !  open(unit=344,file = 'exit')
  if (restart==0) Then
     if (proc0) then
        !       open (unit=2020, File='debugregi')
        open (Unit=33, File= trim(timestep))
        open (Unit=666, File=trim(outflow)) 
        !write(666,'(a5,a15,a15,a15,a15,a15,a15,a15,a15)') 'time', 'Bup_max','Bup_jem',&
        !     & 'outflow_max','delta','L_sheet','niu2'

        write(666,'(a15,a15,a15,a15,a15,a15,a15,a15)') 'niu2', 'res2', 'nu_g', 'hyper_nuei' !Z.Liu 7/10/2020
        !        open (unit=667, file=energy_inj)
        open (Unit=221,File=trim(resol))
        open (unit=99,file=trim(gm_ttrace))
!        open (unit = 1729, file=trim(ppfile)) !AVK: Savetimes will be written to ppfile
     end if
     if (iproc==proc_X)then
        open (Unit=88, File= trim(invariants)) 
        !        write(88,'(a5,a15,a15,a15,a15,a15)') 'time', 'asum', 'nesum', 'A2sum','neAparsum',&
        !             & 'energy'
        write(88,'(a15,a15,a15,a15,a18,a17,a15,a15,a15,a15,a15,a18,a18,a18,a18,a18,a15,a15)') 'savetime',&
             &'b_energytot','de_energytot','ne_energytot','phine_energytot',&
             &'sumgsquaretot','dwdt','inj','ohmic_disstot','visc_disstot','gm_disstot','hyper_eta_disstot',&
             &'hyper_nu_disstot','hyper_gm_disstot', 'hyper_gm_kdisstot','energyzloss','sumleftover' ! Luis                                                         
        open (Unit=77, File= trim(Apar_data))
        write(77,'(a15,a15,a15,a15,a15,a15)') 'time','rec_flux',&
             'uepar_X','width_theo','width_real','vrms','vrms2'
     end if
     !********** INITIAL CONDITIONS *************
     if (debugging) then
     print*, 'checkpointb'
     endif
     ! if equilib_type = 'none' set in input file, all eqs are zero.
     call equilibrium(Apar_eq, Akpar_eq, uepar_eq, uekpar_eq, Apar_eq_double_prime, &
          &           Akpar_eq_double_prime, phi_eq)

     apar_perturb=0.0
     phi_perturb=0.0
     ! if perturb_type = 'none' set in input file, all perturbs are zero.
     call init_perturb(apar_perturb)

     !TTR
#    ifdef gasca2d

     do k = 1, nlz_par
        do j = 1, nly_par
           do i = 1, nlx
              Apar(i, j, k) = Apar_eq(i, j, k) + Apar_perturb(i, j, k)
              phi(i, j, k) = phi_eq(i, j, k) + phi_perturb(i, j, k)
              epar(i, j, k) = 0.0
           end do
        end do

        call FFT2d_direct (Apar(:, :, k), AKpar(:, :, k))
        call FFT2d_direct ( phi(:, :, k),  phiK(:, :, k))

!        if(g_inc) then
        if(anjor) then
           g = 0.0
           gk = 0.0
           do j = 1,nly_par
              do i = 1, nlx
                 g(i, j, k, gmin) = perturb_amp &
                      *cos( inkx*xx(i) * 2.*pi/lx + inky*yy(j) * 2.*pi/ly + inkz*zz(k) * 2.0*pi/Lz )
              end do
           end do
           !BUG to be fixed: ng is not set here
           !deprec. call FFT1d_direct (g(:, :, k, ng), gk(:, :, k, ng))
           call FFT2d_direct (g(:, :, k, ng), gk(:, :, k, ng))
           write(*,*)'BUG: anjor!!!', anjor
        end if
        if (debugging) then
        print *, 'checkpoint c'
        endif
        if (rhoi .lt. small_rhoi) then    !basically if rhoi=0
           do i = 1, nkx_par
              do j = 1, nky
                 neK_perturb(j, i, k) = -kperp(j, i)**2 * phiK(j, i, k)
              end do
           end do
        else
           do i = 1, nkx_par
              do j = 1, nky 
                 neK_perturb(j, i, k) = 2./rhoi**2 * ( gama0(kperp(j,i)**2 * rhoi**2/2.) -1. )*&
                      phiK(j, i, k)
              end do
           end do
        endif

        do i = 1, nkx_par
           do j = 1, nky
              neK(j, i, k) = neK_perturb(j, i, k)
              uekpar(j, i, k) = -kperp(j, i)**2 * Akpar(j, i, k) !uekpar = - \nabla_{\perp}^2 A_{\parallel}  Z.Liu
           end do
        end do

        call FFT2d_inv (ueKpar(:, :, k), uepar(:, :, k))
        call FFT2d_inv (   neK(:, :, k),    ne(:, :, k))     

     end do

#    elif defined(gasca3d)
     if (debugging) then
     print*, 'checkpoint d'
     endif
     do k = 1, nlz_par
        do j = 1, nly_par
           do i = 1, nlx
              Apar(i, j, k) = Apar_eq(i, j, k) + Apar_perturb(i, j, k)
              phi(i, j, k) = phi_eq(i, j, k) + phi_perturb(i, j, k)
              epar(i, j, k) = 0.0
           end do
        end do
     end do

     call FFT2d_direct (Apar(:, :, :), AKpar(:, :, :))
     call FFT2d_direct ( phi(:, :, :),  phiK(:, :, :))
     call FFT2d_direct (Apar(:, :, :), AKpar(:, :, :)) !why repeat?  Z.Liu


!     if(g_inc) then
     if(anjor) then
        !. NOT yet implemented: first fix the BUG in the original code above (gasca2d)
     end if

     if (rhoi .lt. small_rhoi) then    !basically if rhoi=0
        do k = 1, nlz_par
           do i = 1, nkx_par
              do j = 1, nky
                 neK_perturb(j, i, k) = -kperp(j, i)**2 * phiK(j, i, k)
              end do
           end do
        end do
     else
        do k = 1, nlz_par
           do i = 1, nkx_par
              do j = 1, nky 
                 neK_perturb(j, i, k) = 2./rhoi**2 * ( gama0(kperp(j,i)**2 * rhoi**2/2.) -1. )*&
                      phiK(j, i, k)
              end do
           end do
        end do
     endif
        
     do k = 1, nlz_par
        do i = 1, nkx_par
           do j = 1, nky
              neK(j, i, k) = neK_perturb(j, i, k)
              uekpar(j, i, k) = -kperp(j, i)**2 * Akpar(j, i, k)
           end do
        end do
     end do

     call FFT2d_inv (ueKpar(:, :, :), uepar(:, :, :))
     call FFT2d_inv (   neK(:, :, :),    ne(:, :, :))
     call FFT2d_direct (  ne(:, :, :),   nek(:, :, :)) ! ???  Z.Liu

#    endif

     !save initial data:
     !call cuts(1,savetime,proc_X,y_loc_loc,uepar-uepar_eq, apar-apar_eq, phi,vey,runname)
     
  else ! restart=1
     if (iproc==proc_X) then
        open (Unit=88, File= trim(invariants)) 
        open (Unit=77, File= trim(Apar_data)) !,STATUS='UNKNOWN',ACCESS='APPEND') 	    
     end if
     if (proc0) then
        open (Unit=33, File= trim(timestep)) !,STATUS='UNKNOWN',ACCESS='APPEND') 	 
        open (Unit=666, File=trim(outflow)) !,STATUS='UNKNOWN',ACCESS='APPEND')
!       open (Unit=667, File=energy_inj) !,STATUS='UNKNOWN',ACCESS='APPEND')
        open (Unit=221, File=trim(resol)) !,STATUS='UNKNOWN',ACCESS='APPEND')
        open (unit=99, file=trim(gm_ttrace))
!***************************************************************************************
        !read quantities from data file
        if(debugging) then
        print*, 'about to load data'
        endif
        open (unit=2020, File=trim(trim(runname)//'restart_files'))
        write(2020,*) trim(PATH),trim(apar_in)
        write(2020,*) trim(PATH),trim(gfields_in)
        if (debugging) then
        print*, 'finished loading data'
        endif
        close(unit=2020)
        open (unit=2021, File=trim(trim(runname)//'restart_files'))
        read(2021,'(a100)') r1
        read(2021,'(a100)') r3
        close(unit=2021)
        r2=adjustl(r1)
        restart_file1=trim(r2)
        r4=adjustl(r3)
        restart_file2=trim(r4)

!CHANGES MADE ON 09/21/06 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     end if
     !call broadcast(restart_file1)
     !call broadcast(restart_file2)
     !allocate (dumq1(nlx, nlypar_old,nlz_par))
     !allocate (dumq2(nlx, nlypar_old,nlz_par))
     !allocate (dumq3(nlx, nlypar_old,nlz_par))
     call dataloadtest(oldrun, apar,Apar_buff,ne,ne_buff,g,gdummy,epar,rs_time_str)

!     if (file_old) then
!        open (Unit=13, File=restart_file1,STATUS='OLD')!,form='unformatted') 
!        do j=1,iproc*nlypar_old
!           do k=1,nlz
!              do i=1, nlx
!                 READ (13,*) x1,x2
!              end do
!           end do

!           do k=iproc*nlypar_old+1, (iproc+1)*nlypar_old
!              do i=1, nlx
!                 READ (13,*) dumq1(i,k-iproc*nlypar_old), dumq2(i,k-iproc*nlypar_old) 
!           end do
!        end do
!        close(unit=13)
!
!        if(g_inc) then
!           open (Unit=14, File=restart_file2,STATUS='OLD')!,form='unformatted') 
!           !           do ng=2,ngtot
!           do k=1,iproc*nlypar_old
!              do i=1, nlx
!                 READ (14,*) xg(:)
!              end do
!           end do
!           !       end do
!
!           !           do ng=2,ngtot
!           do k=iproc*nlypar_old+1, (iproc+1)*nlypar_old
!              do i=1, nlx
!                 READ (14,*) dumqg(i,k-iproc*nlypar_old,:)
!              end do
!           end do
!           !           end do
!           close(unit=14)
!        end if
!     else
!        !THIS IS NOT WORKING YET!!!!!*************************************!
!
!        open (Unit=13, File=restart_file1,STATUS='OLD',position='rewind',form='unformatted') 
!        do k=1,iproc*nlypar_old
!           read(13) x1
!           read(13) x2
!        end do
!        do k=iproc*nlypar_old+1, (iproc+1)*nlypar_old
!           read(13) dumq1
!           read(13) dumq2
!        end do
!        close(unit=13)
!     end if
!     !*****************************************
!     !double the resolution:
!     if (nly==nly_old) then
!        Apar=dumq1
!        ne=dumq2
!        do ng=2,ngtot
!           do j=1,NLy_par     
!              do i=1,NLx
!                 g(i,j,ng)=dumqg(i,j,ng)
!              end do
!           end do
!        end do
!     else
!        call double_resolution(dumq1,Apar)
!        call double_resolution(dumq2,ne)
!        do ng=2,ngtot
!           call double_resolution(dumqg(:,:,ng),g(:,:,ng))
!        end do
     !     end if

     !TTR
     if (debugging) then
     print*, 'doing the fourier transform'
     endif 
#    ifdef gasca2d
     do k = 1, nlz_par
        call FFT2d_direct (Apar(:, :, k), AKpar(:, :, k))
        call FFT2d_direct (  ne(:, :, k),   nek(:, :, k))
        do i = 1, nkx_par
           do j = 1, nky
              uekpar(j, i, k) = -kperp(j, i)**2 * Akpar(j, i, k)
           end do
        end do
     end do
#    elif defined(gasca3d)
     call FFT2d_direct (Apar(:, :, :), AKpar(:, :, :))
     call FFT2d_direct (  ne(:, :, :),   nek(:, :, :))
     do k = 1, nlz_par
        do i=1, nkx_par
           do j=1, nky
              uekpar(j, i, k) = -kperp(j, i)**2 * Akpar(j, i, k)
           end do
        end do
     end do
#    endif

     if(g_inc) then
        !TTR
#       ifdef gasca2d
        do ng = gmin, ngtot
           do k = 1, nlz_par
              call FFT2d_direct (g(:, :, k, ng), gk(:, :, k, ng))
           end do
        end do
#       elif defined(gasca4d)
        call FFT2d_direct (g(:, :, :, :), gk(:, :, :, :))
#       elif defined(gasca3d)
        do ng = gmin, ngtot
           call FFT2d_direct (g(:, :, :, ng), gk(:, :, :, ng))
        end do
#       endif
     end if
     
     if(debugging) then
     print*, 'almost done'
     endif
     !TTR
     !     call FFT2d_inv (ueKpar, uepar)
     call PHI_POT(neK, phiK) ! compute phik by nek in aux.F90, why cannot be directly computed in loops like uekpar above Z.LIU
     !     call FFT2d_inv(phiK,phi)
     !     !********** INITIAL CONDITIONS *************
     call equilibrium(Apar_eq, Akpar_eq, uepar_eq, uekpar_eq, Apar_eq_double_prime, &
          &           Akpar_eq_double_prime, phi_eq) !??? why here Z.Liu
!     if (equilib_type /= 'tear') then
!        Apar_eq=0.0
!        Akpar_eq=0.0
!        uepar_eq=0.0
!        uekpar_eq=0.0
!     end if


     !     Epar=0.0
     !     vex=0.0
     !     vey=0.0
     !     vrhosx=0.0
     !     vrhosy=0.0
     !!     call cuts(1,savetime,proc_X,y_loc_loc,uepar-uepar_eq)
     !!     call cuts(1,savetime,proc_X,y_loc_loc,apar)
     !!     call predatasave(savetime,howlong)
     !!     call datasavetest(runname, savetime,howlong,apar,ne,epar)
  end if
  if(debugging) then
  print*, 'done with restarting'
  endif 
  !*******restarting done********!

  !TTR
# ifdef gasca2d
  call Convol(phik,   Dxphi,   Dyphi)
  call Convol(nek,    Dxne,    Dyne)
  call Convol(akpar,  Dxapar,  Dyapar)
  call Convol(uekpar, Dxuepar, Dyuepar)
# elif defined(gasca3d)
  call Convol2(phik,   Dxphi,   Dyphi)
  call Convol2(nek,    Dxne,    Dyne)
  call Convol2(akpar,  Dxapar,  Dyapar)
  call Convol2(uekpar, Dxuepar, Dyuepar)
# endif
  if(g_inc) then
#    ifdef gasca2d
     do ng = gmin, ngtot
        call Convol(gk(:, :, :, ng), Dxg(:, :, :, ng), Dyg(:, :, :, ng))
     end do
#    elif defined(gasca4d)
     call Convol2(gk(:, :, :, :), Dxg(:, :, :, :), Dyg(:, :, :, :))
#    elif defined(gasca3d)
     do ng = gmin, ngtot
        call Convol2(gk(:, :, :, ng), Dxg(:, :, :, ng), Dyg(:, :, :, ng))
     end do
#    endif
  end if

 !**********calculate time interval and hyper terms***************

  call flows(dxphi,dyphi,dxne,dyne,vex,vey,vrhosx,vrhosy) ! calculate the exb flow and density gradient (rhos) flow
  vxmax=max(maxval(abs(vex)),maxval(abs(vrhosx)))
  vymax=max(maxval(abs(vey)),maxval(abs(vrhosy)))
  call max_allreduce (vxmax)
  call max_allreduce (vymax)

  call b_field(dxapar,dyapar,bx,by) ! calculate b fields
  bxmax=maxval(abs(bx))
  bymax=maxval(abs(by))
  call max_allreduce(bxmax)
  call max_allreduce(bymax)
  bperp=sqrt(bx**2+by**2)
  bperp_max=maxval(bperp)
  call max_allreduce(bperp_max)
  
  !  call tot_energy(vex,vey,bx,by,energy)
  
  !     call hyper_diff(bperp_max,omega_kaw,niu2,res2)
  call omegakaw(bperp_max,omega_kaw)

  if(g_inc) then
     CFL_flow=min(dx/vxmax,dy/vymax,2./omega_kaw, &
          (1./rhos_de)*1./sqrt(ngtot*1.0)*min(dx/bxmax,dy/bymax), &
          (1./rhos_de)*dz/sqrt(ngtot*1.0), &
          dz*(1.0+max(nky,nkx_par*Npe)*de**2)**0.5/(max(nky,nkx_par*Npe)*rhos),&
          dz/(rhos_de*sqrt(ngtot + 1.0)),1./rhos_de/rhoe_Lte/nky)
  else
! IGNORE this, for now, since hyper_nuei depends on "dti";
!      zdiffcoeff = rhos_de**2*(ngtot+1)/((ngtot+1)*nu_ei+(ngtot+1)**(2*hyper_morder)*hyper_nuei)
!      CFL_flow=min(dx/vxmax,dy/vymax,2./omega_kaw, &
!           (1./rhos_de)*1./sqrt(ngtot*1.0)*min(dx/bxmax,dy/bymax), &
!           dz*(1.0+max(nky,nkx_par*Npe)*de**2)**0.5/(max(nky,nkx_par*Npe)*rhos),&
!           dz/(rhos_de*sqrt(ngtot + 1.0)),sqrt(2.0)/rhos_de/rhoe_Lte/nky,&
!           dz**2/(zdiffcoeff), 1/(zdiffcoeff*kperpmax**2*(max(bxmax,bymax))**2)) ! Two new terms due to Hermite closure, LUIS 25/5/2015

     !   CFL_flow=min(dx/vxmax,dy/vymax,2./omega_kaw)
     CFL_flow=min(dx/vxmax,dy/vymax,dy/bymax,dx/bxmax,dz,2./omega_kaw, &
          dz*(1.0+max(nky,nkx_par*Npe)*de**2)**0.5/(max(nky,nkx_par*Npe)*rhos))
     !     sqrt(2.0)/rhos_de/rhoe_Lte/nky)   ! for the case of d_e=0.1
  endif
  
  !CFL_flow=min(dx/vxmax,dy/vymax,dx/Bxmax,dy/Bymax)  
  dti=CFL_frac*CFL_flow

  if(proc0) then
      write(666,50) savetime, CFL_frac*dx/vxmax, CFL_frac*dy/vymax, &
      CFL_frac*dx/bxmax, CFL_frac*dy/bymax, CFL_frac*2./omega_kaw, &
      CFL_frac*(1./rhos_de)*1./sqrt(ngtot*1.0)*min(dy/bymax,dx/bxmax), &
      CFL_frac*(1./rhos_de)*dz/sqrt(ngtot*1.0), & 
      CFL_frac*dz*(1.0+max(nky,nkx_par*Npe)*de**2)**0.5/(max(nky,nkx_par*Npe)*rhos), &
      CFL_frac*dz/(rhos_de*sqrt(ngtot + 1.0)),&
      CFL_frac*min(dx/vxmax,dy/vymax,2./omega_kaw, &
      (1./rhos_de)*1./sqrt(ngtot*1.0)*min(dx/bxmax,dy/bymax), &
      (1./rhos_de)*dz/sqrt(ngtot*1.0), dz/(rhos_de*sqrt(ngtot + 1.0)), &
      dz*(1+max(nky,nkx_par*Npe)*de**2)**0.5/(max(nky,nkx_par*Npe)*rhos)),dti!p_count*1.0   !p_count*1.0 
50         format(18g16.8)  !,f6.1)
  endif

  if (hyper_fixed) then 
     nu_g=hyper_nu_g
     niu2=hyper_nu
     res2=hyper_eta
  else
!     nu_g=hyper_coef_g/dti/(NLx/2.)**(2*hyper_order_g)
     nu_g=hyper_coef_g/dti/(kperpmax)**(2*hyper_order_g)
!     niu2=hyper_coef/dti/(NLx/2.)**(2*hyper_order)
     niu2=hyper_coef/dti/(kperpmax)**(2*hyper_order)
!     if ((nkx_par*npe)**2*de**2>1) then
     if (kperpmax**2*de**2>1) then
!        res2=hyper_coef/dti/(NLx/2.)**(2*hyper_order-2)*de**2
        res2=hyper_coef/dti/(kperpmax)**(2*hyper_order-2)*de**2
     else
!       cd ..
        res2=hyper_coef/dti/(NLx/2.)**(2*hyper_order)
        res2=hyper_coef/dti/(kperpmax)**(2*hyper_order)
     end if
  end if
  !Hyper-collisions:
  if (hyper_colls_fixed) then
     hyper_nuei=hyper_colls
  else
     hyper_nuei=hyperm_coef/dti/(Ngtot+1)**(2*hyper_morder) 
     ! Luis - as per PRL13, eq. 4, hypercollisions;                                                    
  end if
  !**********End of calculating time interval and hyper terms***************

  !initializations:
  RELATIVE_ERROR=0.0
  !dti=0.001
  aa0=init_aa0_fac !*.5
  total=0.0
  uxavg=0.0
  uyavg=0.0
  inj=0.0
  energyzloss=0.0
  sumleftover=0.0
  energyzloss1=0.0
  energyzloss2=0.0
!  print*, iproc, maxval(abs(Akpar))
  call diagnostics(phiK, nek, uekpar,uekpar_eq,Akpar,&
       Akpar_eq,gk,sumleftover,savetime,res2,niu2,nu_g,hyper_nuei,dti,proc_X,inj,energyzloss)

!  if(proc0) then !Z.Liu 7/10/2020
!     write(666,50) niu2, res2, nu_g, hyper_nuei
!50   format(10g16.8)
!  endif
  
  if(debugging) then
  print*, 'may God help us'
  endif
!  stop
  !******* THIS ENDS THE INITIAL CONDITIONS SECTION **************
  ! AVK debugging
923 continue
  if(proc0) then
     call cpu_time(time1)
     write(*,*) "Initial conditions done, t=", time1
  endif
  Myperfoff ! CodeInit
  !***************** TIME LOOP *******************************

  p=0
  z=0
  repeat=.false.
  noinc=.false.
  divergent=.false.

  Myperfon('TimeLoop')
  time_loop: do t=0, tmax
     if (iproc .eq. 0 .and. mod(t,1) .eq. 0) &
          write(*,'(a,i6)') 'time iteration: ', t
!     if(proc0) then
!        call date_and_time(date, time)
!        write(*,*) "Time, timestep, clock = ", savetime, t, time
!     end if
     p=t-z

!     if(three_D .eqv. .true.) then!
!
!        !Z-STEP (strang split)                                                                         
!        !collect energy information 
!        call zsteploss(phiK, nek, uekpar,Akpar,gk,dti/2.,proc_X,energyzloss)
!        !advance eqs:
!     if(zmethod=='cormack') then
!        call z_step(phik,nek,uekpar,akpar,gk,dti/2.)
!     else
!        call z_step2(phik,nek,uekpar,akpar,gk,dti/2.)
!     end if
!        call zsteploss(phiK, nek, uekpar,Akpar,gk,dti/2.,proc_X,energyzloss1)
!     end if
     
     !call force(only called once per timestep)

!     dummy(:,:,:)=0.0  !NFL: unnecessary
     if(debugging) then
     print*, 'about to start with turb :)'
     endif
     if(turb) then
!     if(turb .and. t==0) then
        !choosing kz
        kfz=nint(kfz1+1.d0*uniran()*(kfz2-kfz1))
        ! random sign
        kfz=kfz*nint(sign(1d0,(2d0*uniran()-1d0)))
        !NFL there's a bug in the forcing routine, line 136,137 if j=1
        call force(kfz,dti,kfp1,kfp2,feps,dummy)
     else
        dummy=0.0
     end if
     if(debugging) then
     print*, 'done with turb'
     endif
     !AVK: adding forcing for g
!     if(g_inc .and. gturb) then
!        !fac=Lz/2./pi
!        g_kfz=nint(g_kfz1*fac-.5d0+uniran()*(fac*g_kfz2-fac*g_kfz1+1d0))/fac
!        ! random sign
!        g_kfz=g_kfz*nint(sign(1d0,(2d0*uniran()-1d0)))
!        call force(g_kfz,dti,g_kfp1,g_kfp2,g_feps,dummy2)
!        !AVK: This is for gplus calculation, gminus will have a different scaling
!        dummy2 = dummy2 * ((1./sigma)*(1.+zcharge/tite) +1.) *(1./lambda)
!     else
!        dummy2=0.0
!     end if


     if (repeat) then
        repeat=.false.
     else
        if (divergent) then
           divergent=.false.
        else 
           p_count=0
           if (first) then
              !TTR
#             ifdef gasca2d
              call Convol(phiK,  Dxphi,    Dyphi)
              call Convol(nek,   Dxne,     Dyne)
              call Convol(AKpar, DxApar,   DyApar)
              call Convol(uekpar, Dxuepar, Dyuepar)
#             elif defined(gasca3d)
              call Convol2(phiK,  Dxphi,    Dyphi)
              call Convol2(nek,   Dxne,     Dyne)
              call Convol2(AKpar, DxApar,   DyApar)
              call Convol2(uekpar, Dxuepar, Dyuepar)
#             endif
              if (g_inc) then
                 Myperfon('ngTLoop')
#                ifdef gasca2d
                 do ng = gmin, ngtot
                    call Convol( gk(:, :, :, ng), Dxg(:, :, :, ng), Dyg(:, :, :, ng) )
                 end do
                 call Convol2( gk(:, :, :, :), Dxg(:, :, :, :), Dyg(:, :, :, :) )
#                elif defined(gasca3d)
                 do ng = gmin, ngtot
                    call Convol2( gk(:, :, :, ng), Dxg(:, :, :, ng), Dyg(:, :, :, ng) )
                 end do
#                endif
                 Myperfoff ! ngTLoop
              end if
              first=.false.
              !         else
              !            DxApar(:,:)=Dxaparnew(:,:)
              !            dyapar(:,:)=dyaparnew(:,:)
              !            dxuepar(:,:)=dxueparnew(:,:)
              !            dyuepar(:,:)=dyueparnew(:,:)
           end if
           if(debugging) then 
           print*, 'starting to do the first computation'
           endif

           call funcne_i(Dxphi,Dyphi,Dxne,Dyne,DxApar,DyApar,Dxuepar,Dyuepar,Fne_old, &
                braakparuekpar)
           if(g_inc) then
              call funcAKpar_i(DxApar,DyApar,Dxphi,Dyphi,Dxne,Dyne,Dxuepar,Dyuepar, &
                   Dxg(:,:,:,gmin),Dyg(:,:,:,gmin), akpar, FApar_old, savetime)
           else
              call funcAKpar_i(DxApar,DyApar,Dxphi,Dyphi,Dxne,Dyne,Dxuepar,Dyuepar, &
                   dummy_real,dummy_real, akpar, FApar_old, savetime)
           end if

           if(g_inc) then
              call funcg2(Dxg(:, :, :, gmin), Dyg(:, :, :, gmin), Dxphi, Dyphi, &
                   &      Dxapar, Dyapar, Dxg(:, :, :, gmin+1), Dyg(:, :, :, gmin+1), &
                   &      braakparuekpar, phik, fgm_old(:, :, :, gmin))

              Myperfon('ngTLoop-')
              !TTR
#             ifdef gasca4d
              call funcgm(Dxg(:, :, :, :), Dyg(:, :, :, :), &
                   &      Dxg(:, :, :, :), Dyg(:, :, :, :), &
                   &      Dxphi(:, :, :), Dyphi(:, :, :), Dxapar(:, :, :), Dyapar(:, :, :), Akpar(:, :, :), &
                   &      fgm_old(:, :, :, :), savetime)
#             else
              do ng = gmin+1, ngtot-1
                 !               if (ng<ngtot) then
                 call funcgm(ng, Dxg(:, :, :, ng-1), Dyg(:, :, :, ng-1), &
                      &          Dxg(:, :, :, ng),   Dyg(:, :, :, ng), &
                      &          Dxg(:, :, :, ng+1), Dyg(:, :, :, ng+1), &
                      &          Dxphi, Dyphi, Dxapar, Dyapar, akpar, &
                      &          fgm_old(:, :, :, ng), savetime)
                 !               else
                 !                  call funcgm(ng,Dxg(:,:,ng-1),Dyg(:,:,ng-1),Dxg(:,:,ng),Dyg(:,:,ng), &
                 !                       & dummy_real,dummy_real,Dxphi,Dyphi,Dxapar,Dyapar,fgm_old(:,:,ng))
                 !               end if
              end do
#             endif
              Myperfoff ! ngTLoop-

              !>
              !NFL, 06/12/2013:
              call func_lastg(hyper_nuei, niu2, Dxg(:, :, :, ngtot-1), Dyg(:, :, :, ngtot-1), &
                   &                            Dxg(:, :, :, ngtot),   Dyg(:, :, :, ngtot), &
                   &          Dxphi, Dyphi, Dxapar, Dyapar, f_lastg_old)
              !<
           end if

        end if
     end if

     call SEMI_IMP_OP(dti,bperp_max,aa0,SI_oper) ! ??? in aux.F90

     !***** Predictor step calculation begins *******
     do i=1,nkx_par
        do j=1,nky
           guess(j,i,:)=Akpar(j,i,:)
           
           nek_star(j,i,:)=exp_nu(j,i,niu2,dti)*nek(j,i,:)+ &
                dti/2.*(1.0+exp_nu(j,i,niu2,dti))*Fne_old(j,i,:)
           !add the forcing (if it's on -- see constants):
           nek_star(j,i,:)=nek_star(j,i,:)+ &
                dti*(1.0-0.5*(niu*kperp(j,i)**2+niu2*kperp(j,i)**(2*hyper_order))*dti)*dummy(j,i,:)  
           
           Akpar_star(j,i,:)=exp_eta(j,i,res2,dti)*Akpar(j,i,:)+ &
                dti/2.*(1.0+exp_eta(j,i,res2,dti))*Fapar_old(j,i,:)+ &
                (1.0-exp_eta(j,i,res2,dti))*Akpar_eq(j,i,:)!+ &
             !   dti*(1.0-0.5*(res2*kperp(j,i)**2+res2*kperp(j,i)**(2*hyper_order))*dti)*dummy(j,i,:)
           
           uekpar_star(j,i,:)=-kperp(j,i)**2*Akpar_star(j,i,:)
        end do
     end do
     if(debugging) then
     print*, 'moving forward'
     endif
     !AVK: gmin is done differently because no collisions on gmin
     if (g_inc) then
        do i = 1, nkx_par
           do j = 1, nky
              gk_star(j,i,:,gmin)= &
                   exp_nu(j, i, nu_g, dti) * gk(j, i, :, gmin)+ &
                   dti/2.*(1.0+exp_nu(j,i,nu_g,dti))*Fgm_old(j,i,:,gmin)
                !   nu_func1(j,i,nu_g,dti)*Fgm_old(j,i,:,gmin)
           end do
        end do
	  !AVK: Removed white noise forcing from g1
        do ng = gmin+1, ngtot-1
           do i=1,nkx_par
              do j=1,nky
                 gk_star(j,i,:,ng)=exp_ng(ng,hyper_nuei,dti)*exp_nu(j,i,nu_g,dti)*gk(j,i,:,ng)+ &
                      dti/2.*(1.0+exp_ng(ng,hyper_nuei,dti)*exp_nu(j,i,nu_g,dti))*Fgm_old(j,i,:,ng)
                      !nu_func2(j,i,nu_g,ng,hyper_nuei,dti)*fgm_old(j,i,:,ng)
              end do
           end do
        end do
        !This is where the closure should be called if not func_lastg        
!        call closure(gk_star(:,:,ngtot-1)) !---- the arguments depend on the order of the closure
        do i=1,nkx_par
           do j=1,nky
              gk_star(j,i,:,ngtot)=exp_ng(ngtot,hyper_nuei,dti)*exp_nu(j,i,nu_g,dti)*gk(j,i,:,ngtot)+ &
                   dti/2.*(1.0+exp_ng(ngtot,hyper_nuei,dti)*exp_nu(j,i,nu_g,dti))*f_lastg_old(j,i,:)
                   !nu_func2(j,i,nu_g,ngtot,hyper_nuei,dti)*f_lastg_old(j,i,:)
           end do
        end do
   !***** Predictor step calculation ends *******

        Myperfon('ngTLoop')
        !TTR
#       ifdef gasca2d
        do ng = gmin, ngtot
           call Convol(gk_star(:, :, :, ng), Dxg_star(:, :, :, ng), Dyg_star(:, :, :, ng))
        end do
#       elif defined(gasca4d)
        call Convol2(gk_star(:, :, :, :), Dxg_star(:, :, :, :), Dyg_star(:, :, :, :))
#       elif defined(gasca3d)
        do ng = gmin, ngtot
           call Convol2(gk_star(:, :, :, ng), Dxg_star(:, :, :, ng), Dyg_star(:, :, :, ng))
        end do
#       endif
        Myperfoff ! ngTLoop
     end if
     if(debugging) then
     print*, 'let us see how far we can get'
     endif
!	 !AVK: Is gmin done differently because no collisions on gmin?
!     if (g_inc) then
!        gk_star(:,:,:,gmin)=exp_nu(j,i,nu_g,dti)*gk(:,:,:,gmin)+ &
!             dti/2.*(1.0+exp_nu(j,i,nu_g,dti))*Fgm_old(:,:,:,gmin)
!        call CONVOL(gk_star(:,:,:,gmin),dxg_star(:,:,:,gmin),dyg_star(:,:,:,gmin))
!		! AVK: add forcing for g1 here?
!        do ng=gmin+1,Ngtot-1
!           gk_star(:,:,:,ng)=exp_ng(ng,hyper_nuei,dti)*exp_nu(j,i,nu_g,dti)*gk(:,:,:,ng)+ &
!                dti/2.*(1.0+exp_ng(ng,hyper_nuei,dti)*exp_nu(j,i,nu_g,dti))*Fgm_old(:,:,:,ng)
!           call CONVOL(gk_star(:,:,:,ng),dxg_star(:,:,:,ng),dyg_star(:,:,:,ng))
!        end do
!        !M closure:
!        
!        !      gk_star(:,:,ngtot)=  4.0*gk_star(:,:,ngtot-1) - 6.0*gk_star(:,:,ngtot-2) &
!        !           + 4.0*gk_star(:,:,ngtot-3) - gk_star(:,:,ngtot-4)
!        !      gk_star(:,:,ngtot) = 2.0*gk_star(:,:,ngtot-1) - gk_star(:,:,ngtot-2)
!        gk_star(:,:,:,ngtot) = gk_star(:,:,:,ngtot-1)
!        !      gk_star(:,:,ngtot) = 0.0
!        
!        call CONVOL(gk_star(:,:,:,ngtot),dxg_star(:,:,:,ngtot),dyg_star(:,:,:,ngtot))
!     end if
     
     call PHI_POT(nek_star,phik_star)

     !TTR
#    ifdef gasca2d
     call Convol(phik_star,   Dxphi_star,   Dyphi_star)
     call Convol(nek_star,    Dxne_star,    Dyne_star)
     call Convol(AKpar_star,  DxApar_star,  DyApar_star)
     call Convol(uekpar_star, Dxuepar_star, Dyuepar_star)
#    elif defined(gasca3d)
     call Convol2(phik_star,   Dxphi_star,   Dyphi_star)
     call Convol2(nek_star,    Dxne_star,    Dyne_star)
     call Convol2(AKpar_star,  DxApar_star,  DyApar_star)
     call Convol2(uekpar_star, Dxuepar_star, Dyuepar_star)
#    endif
   
     if(g_inc) then
        call funcAKpar_i(DxApar_star,DyApar_star,Dxphi_star,Dyphi_star,&
             & Dxne_star,Dyne_star,Dxuepar_star,Dyuepar_star, &
             & Dxg_star(:,:,:,gmin),Dyg_star(:,:,:,gmin), akpar_star, FApar_star, savetime)
     else
        call funcAKpar_i(DxApar_star,DyApar_star,Dxphi_star,Dyphi_star,&
             & Dxne_star,Dyne_star,Dxuepar_star,Dyuepar_star, &
             & dummy_real,dummy_real, akpar_star, FApar_star, savetime)
     end if
     
     !***** P_Loop begins, this is to achieve desired convergence in one time step!!*****
     Myperfon('PLoop')
     P_LOOP:   do p_iter = 0, pmax
        
        p_count=p_count+1
        rel_err=0.0 
        rel_error=0.0
        old_error=0.0
        
        if(p_iter == 0) then
           FAparpred = FApar_star 
           value_nex = Dxne_star 
           value_ney = Dyne_star 
           value_phix = Dxphi_star 
           value_phiy = Dyphi_star 
           value_gx = Dxg_star
           value_gy = Dyg_star

        else
           
           old_error=RELATIVE_ERROR
           if (g_inc) then
              call funcAKpar_i(DxApar,DyApar,Dxphi,Dyphi,Dxne,Dyne,&
                   Dxuepar,Dyuepar,Dxg(:,:,:,gmin),Dyg(:,:,:,gmin),akparnew,FAparpred,savetime)
           else
              call funcAKpar_i(DxApar,DyApar,Dxphi,Dyphi,Dxne,Dyne,&
                   Dxuepar,Dyuepar,dummy_real,dummy_real,akparnew, FAparpred, savetime)
           end if

           value_nex = Dxne
           value_ney = Dyne 
           value_phix = Dxphi 
           value_phiy = Dyphi 
           value_gx = Dxg
           value_gy = Dyg

        end if
        

        RELATIVE_ERROR=0.0
        do k=1,nlz_par
           do i=1,NKx_par
              do j=1,NKy
                 AKparnew(j,i,k)=1./(1.+SI_oper(j,i,k)/4.)*&
                      (exp_eta(j,i,res2,dti)*Akpar(j,i,k)+&
                      dti/2.*exp_eta(j,i,res2,dti)*Fapar_old(j,i,k)+&
                      dti/2.*Faparpred(j,i,k)+&
                      (1.0-exp_eta(j,i,res2,dti))*Akpar_eq(j,i,k)+&
                      SI_oper(j,i,k)/4.*guess(j,i,k))

                 akparnew(j,i,k)=akparnew(j,i,k)! + &
                  ! dti*(1.-0.5*(niu*kperp(j,i)**2+niu2*kperp(j,i)**(2*hyper_order))*dti)*dummy(j,i,k)

                 ueKparnew(j,i,k)=-kperp(j,i)**2*AKparnew(j,i,k)

                 rel_err=rel_err+(abs(Akparnew(j,i,k)-Akpar(j,i,k)))**2 ! rel_err is the change of Akpar in one time step, accumulated                
              end do
           end do
        end do

        call sum_allreduce(rel_err)

        !TTR
#       ifdef gasca2d
        call Convol(AKparnew,  DxApar,  DyApar)
        call Convol(uekparnew, Dxuepar, Dyuepar)   
#       elif defined(gasca3d)
        call Convol2(AKparnew,  DxApar,  DyApar)
        call Convol2(uekparnew, Dxuepar, Dyuepar)   
#       endif
        call funcne_i(value_phix, value_phiy, value_nex, value_ney, DxApar, DyApar, & 
             Dxuepar, Dyuepar, Fnepred, braakparuekpar)

        do i=1,NKx_par
           do j=1,NKy
              neknew(j,i,:)=exp_nu(j,i,niu2,dti)*nek(j,i,:)+&
                   dti/2.*exp_nu(j,i,niu2,dti)*Fne_old(j,i,:)+&
                   dti/2.*Fnepred(j,i,:)
              
              neknew(j,i,:)=neknew(j,i,:)+ &
                   dti*(1.-0.5*(niu*kperp(j,i)**2+niu2*kperp(j,i)**(2*hyper_order))*dti)*dummy(j,i,:)
              
              !CALCULATE THE ERROR FOR THIS P-ITERATION
              rel_error(j,i,:)=abs(SI_oper(j,i,:)/4.*(Akparnew(j,i,:)-guess(j,i,:)))/ &
                   sqrt(rel_err/(Nkx*Nky)) !NOTE that guess:=aparnew at the end of one loop
           end do
        end do
        
        call PHI_POT(neKnew, phiKnew)

        !TTR
#       ifdef gasca2d
        call Convol(phiKnew, Dxphi, Dyphi)
        call Convol(neknew,  Dxne,  Dyne)
#       elif defined(gasca3d)
        call Convol2(phiKnew, Dxphi, Dyphi)
        call Convol2(neknew,  Dxne,  Dyne)
#       endif

        !      totalgk=0.0
        !      totalgknew=0.0

        if(g_inc) then
           call funcg2(value_gx(:, :, :, gmin), value_gy(:, :, :, gmin), &
                &      Dxphi, Dyphi, Dxapar, Dyapar, &
                &      value_gx(:, :, :, gmin+1), value_gy(:, :, :, gmin+1), &
                &      braakparuekpar, phiknew, Fgm_pred(:,:,:,gmin))

           do i = 1, nkx_par
              do j = 1, nky
                 gknew(j, i, :, gmin) = exp_nu(j,i,nu_g,dti)*gk(j,i,:,gmin)+ &
                      dti/2.*exp_nu(j,i,nu_g,dti)*Fgm_old(j,i,:,gmin)+dti/2.*Fgm_pred(j,i,:,gmin)
                      !nu_func3(j,i,nu_g,dti)*fgm_old(j,i,:,gmin) +&
                      !nu_func4(j,i,nu_g,dti)*fgm_pred(j,i,:,gmin)

                 !               totalgknew=totalgknew+abs(gknew(j,i,2))**2
              end do
           end do

           !         call sum_allreduce(totalgknew)

           !TTR
#          ifdef gasca2d
           call Convol ( gknew(:, :, :, gmin), Dxg(:, :, :, gmin), Dyg(:, :, :, gmin) )
#          elif defined(gasca3d)
           call Convol2( gknew(:, :, :, gmin), Dxg(:, :, :, gmin), Dyg(:, :, :, gmin) )
#          endif

           Myperfon('ngPLoop-')
           !. here the gasca4d can NOT be used because Dxgm and Dygm depend on previous ng-value calculation
           do ng = gmin+1, ngtot-1
              !totalgk=totalgknew
              
              Dxgm(:, :, :) = Dxg(:, :, :, ng-1) ! need the p+1 gm-1 values, which are updated within this loop
              Dygm(:, :, :) = Dyg(:, :, :, ng-1)

              !            if (ng<ngtot) then
              call funcgm(ng,     Dxgm(:, :, :),           Dygm(:, :, :), &
                   &          value_gx(:, :, :, ng),   value_gy(:, :, :, ng), &
                   &          value_gx(:, :, :, ng+1), value_gy(:, :, :, ng+1), &
                   &          Dxphi, Dyphi, DxApar, DyApar, Akparnew, &
                   &          fgm_pred(:, :, :, ng), savetime)
              !            else
              !               call funcgm(ng,Dxgm,Dygm,value_gx(:,:,ng),value_gy(:,:,ng), &
              !                    & dummy_real, dummy_real,Dxphi,Dyphi,Dxapar,Dyapar,fgm_pred(:,:,ng))
              !            end if

              do i=1, nkx_par
                 do j=1,nky
                    !gknew(j,i,ng)=exp(-(ng*nu_ei+ng**4*hyperm_coef*nu_ei/Ngtot**3)*dti)*gk(j,i,ng)+ &
                    !dti/2.*exp(-(ng*nu_ei+ng**4*hyperm_coef*nu_ei/Ngtot**3)*dti)*fgm_old(j,i,ng)+ &
                    !dti/2.*fgm_pred(j,i,ng)
                    
                    gknew(j,i,:,ng) = exp_ng(ng,hyper_nuei,dti)*exp_nu(j,i,nu_g,dti)*gk(j,i,:,ng)+ &
                         dti/2.*exp_ng(ng,hyper_nuei,dti)*exp_nu(j,i,nu_g,dti)*fgm_old(j,i,:,ng)+ &
                         dti/2.*fgm_pred(j,i,:,ng)
                         !nu_func5(j,i,nu_g,ng,hyper_nuei,dti)*fgm_old(j,i,:,ng)+&
                         !nu_func6(j,i,nu_g,ng,hyper_nuei,dti)*fgm_pred(j,i,:,ng)
                    
                    !                  totalgknew=totalgknew+abs(gknew(j,i,ng))**2
                 end do
              end do
              
              !           call sum_allreduce(totalgknew)
              

              !. calculation for the next Dxgm and Dygm's, which prevents using gasca4d here
              !TTR
#             ifdef gasca2d
              call Convol( gknew(:, :, :, ng), Dxg(:, :, :, ng), Dyg(:, :, :, ng) )
#             elif defined(gasca3d)
              call Convol2( gknew(:, :, :, ng), Dxg(:, :, :, ng), Dyg(:, :, :, ng) )
#             endif
           end do
           Myperfoff ! ngPLoop-


           call func_lastg(hyper_nuei, niu2, Dxg(:,:,:,ngtot-1), Dyg(:,:,:,ngtot-1),&
                value_gx(:,:,:,ngtot), value_gy(:,:,:,ngtot), Dxphi, Dyphi, Dxapar, Dyapar, f_lastg)

           do i=1, nkx_par
              do j=1, nky
                 gknew(j,i,:,ngtot) = exp_ng(ngtot,hyper_nuei,dti)*exp_nu(j,i,niu2,dti)*gk(j,i,:,ngtot)+ &
                      dti/2.*exp_ng(ngtot,hyper_nuei,dti)*exp_nu(j,i,niu2,dti)*f_lastg_old(j,i,:)+ &
                      dti/2.*f_lastg(j,i,:)
                      !nu_func5(j,i,nu_g,ng,hyper_nuei,dti)*f_lastg_old(j,i,:)+&
                      !nu_func6(j,i,nu_g,ng,hyper_nuei,dti)*f_lastg(j,i,:)
              end do
           end do

!         call closure(gknew(:,:,ngtot-2),gknew(:,:,ngtot-1),ngtot-1,hyper_nuei,gknew(:,:,ngtot))

           !TTR
#          ifdef gasca2d
           call Convol( gknew(:, :, :, ngtot), dxg(:, :, :, ngtot), dyg(:, :, :, ngtot) )
#          elif defined(gasca3d)
           call Convol2( gknew(:, :, :, ngtot), dxg(:, :, :, ngtot), dyg(:, :, :, ngtot) )
#          endif
        end if

        RELATIVE_ERROR=maxval(abs(rel_error))
        call max_allreduce(RELATIVE_ERROR)
        
        if (.not. turb) then
           if (p_iter.GE. 1 .AND. relative_error/old_error .GE. 1.0) then
              dti=low*dti
              z=z+1
              !      repeat=.true.
              divergent=.true.
              exit p_loop
           end if
           
           if (RELATIVE_ERROR .LE. epsilon) then
              exit P_LOOP
           end if
           
           if ((RELATIVE_ERROR .GT. epsilon) .AND. (p_iter==pmax)) then   !REDUCE TIMESTEP AND REDO P_LOOP
              dti=low*dti
              z=z+1
              repeat=.true.
              exit p_loop
           end if
        end if

        do i=1, NKx_par
           do j=1, NKy
              !         neKpred(j,i)=neKnew(j,i)
              !         ueKparpred(j,i)=ueKparnew(j,i)
              !         AKparpred(j,i)=AKparnew(j,i)
              !         phiKpred(j,i)=phiKnew(j,i)
              guess(j,i,:)=Akparnew(j,i,:)
           end do
        end do

      
   end do P_LOOP
   Myperfoff ! PLoop


     if (divergent) cycle time_loop
     
     if (repeat) then
        noinc=.true.
        cycle time_loop
     end if
     
     !********* UPDATE THE VARIABLES ***********
     
     Myperfon('UpdateVa')
     neK=neKnew
     AKpar=AKparnew
     phiK=phiKnew
     do k=1,nlz_par
        do i=1,nkx_par
           do j=1,nky
              ueKpar(j,i,k)=-kperp(j,i)**2*Akpar(j,i,k)
           end do
        end do
     end do
     if (g_inc) gk=gknew
     
     if(three_D .eqv. .true.) then
        
        !NOW DO THE Z-STEP;
        !collect energy information
        call zsteploss(phiK, nek, uekpar,Akpar,gk,dti,proc_X,energyzloss)
        if(zmethod=='cormack') then
           call z_step(phik,nek,uekpar,akpar,gk,dti)
        else
!           call z_step2(phik,nek,uekpar,akpar,gk,dti/2.)
        call z_step2(phik,nek,uekpar,akpar,gk,dti)
        end if

!        call zsteploss(phiK, nek, uekpar,Akpar,gk,dti/2.,proc_X,energyzloss2)
        call zsteploss(phiK, nek, uekpar,Akpar,gk,dti,proc_X,energyzloss2)
        !NFL, 20/03/2013:
        energyzloss=energyzloss1+energyzloss2 ! why this? so that we remember that for Strang 
        !there would be 2 variables (2 calls) needed... LUIS

        !NFL, 25/11/2014: call d^2/dz^2 closure term 
        if(g_inc) then
           !zdiffcoeff = rhos_de**2*(ngtot+1)/((ngtot+1)*nu_ei+(ngtot+1)**(2*hyper_morder)*hyper_nuei)
           zdiffcoeff = 0.0 !Z.Liu 10/15/2020
           call z_diffusion(zdiffcoeff,dti,gk(:,:,:,ngtot),gk(:,:,:,ngtot))
        end if
     end if
     
     savetime=savetime+dti
     Myperfoff ! UpdateVa

     !*********     

     !DIAGNOSTICS:
     
     !Now check if dti is smaller than timestep required by the flows,
     !and determine the hyper-diffusion coefficients (if needed)

     Myperfon('Diags')
     
     call flows(dxphi,dyphi,dxne,dyne,vex,vey,vrhosx,vrhosy)
     vxmax=max(maxval(abs(vex)),maxval(abs(vrhosx)))
     vymax=max(maxval(abs(vey)),maxval(abs(vrhosy)))
     call max_allreduce (vxmax)
     call max_allreduce (vymax)
   
     call b_field(dxapar,dyapar,bx,by)
     bxmax=maxval(abs(bx))
     bymax=maxval(abs(by))
     call max_allreduce(bxmax)
     call max_allreduce(bymax)
     bperp=sqrt(bx**2+by**2)
     bperp_max=maxval(bperp)
     call max_allreduce(bperp_max)
     call omegakaw(bperp_max,omega_kaw)

     !*******************************************************
     !calculate urms:
     !First average the total flow over time to get the mean flow 
     !as a function of (x,y):
   
     uxavg =(uxavg +vex )/(p+1.)
     uyavg =(uyavg +vey )/(p+1.)

!*******************************************************


     !  old_energy=energy
     !  call tot_energy(vex,vey,bx,by,energy)
     
     !  denergy_dt=(energy-old_energy)/dti
     !  total=total+denergy_dt
     !  number=number+1
     !  av_denergy_dt = total/savetime

     if (turb) then
        if(g_inc) then
            CFL_flow=min(dx/vxmax,dy/vymax,2./omega_kaw, &
                  (1./rhos_de)*1./sqrt(ngtot*1.0)*min(dx/bxmax,dy/bymax), &
                  (1./rhos_de)*dz/sqrt(ngtot*1.0), &
                  dz*(1.0+max(nky,nkx_par*Npe)*de**2)**0.5/(max(nky,nkx_par*Npe)*rhos),&
                  dz/(rhos_de*sqrt(ngtot + 1.0)),1./rhos_de/rhoe_Lte/nky)
           !zdiffcoeff = rhos_de**2*(ngtot+1)/((ngtot+1)*nu_ei+(ngtot+1)**(2*hyper_morder)*hyper_nuei)
           !CFL_flow=min(dx/vxmax,dy/vymax,2./omega_kaw, &
           !     (1./rhos_de)*1./sqrt(ngtot*1.0)*min(dx/bxmax,dy/bymax), &
           !     dz*(1.0+max(nky,nkx_par*Npe)*de**2)**0.5/(max(nky,nkx_par*Npe)*rhos),&
           !     dz/(rhos_de*sqrt(ngtot + 1.0)),sqrt(2.0)/rhos_de/rhoe_Lte/nky,&
           !     dz**2/(zdiffcoeff), 1/(zdiffcoeff*kperpmax**2*(max(bxmax,bymax))**2)) ! Two new terms due to Hermite closure, LUIS 25/5/2015
         
        else
           !         CFL_flow=min(dx/vxmax,dy/vymax,2./omega_kaw)

     CFL_flow=min(dx/vxmax,dy/vymax,dy/bymax,dx/bxmax,dz,2./omega_kaw, &
          dz*(1.0+max(nky,nkx_par*Npe)*de**2)**0.5/(max(nky,nkx_par*Npe)*rhos))
!           CFL_flow=min(dx/vxmax,dy/vymax,dy/bymax,dx/bxmax,dz,&
 !               dz*(1+max(nky,nkx_par*Npe)*de**2)**0.5/(max(nky,nkx_par*Npe)*rhos),&
 !               sqrt(2.0)/rhos_de/rhoe_Lte/nky)
        endif
        !      CFL_flow=min(dx/vxmax,dy/vymax,dx/Bxmax,dy/Bymax)  
     else
        CFL_flow=min(dx/vxmax,dy/vymax,dz,&
             dz*(1+max(nky,nkx_par*Npe)*de**2)**0.5/(max(nky,nkx_par*Npe)*rhos))  !,1./(ky(nky/2+1)*bymax))
     end if
   
     x=CFL_frac*CFL_flow ! CFL_frac = 0.25 in constants 
     !*******
   
     If(mod(p,save_datafiles)==0)THEN 
!        if(proc0) then
!           call date_and_time(date, time)
!           write(*,*) "Time, timestep, clock = ", savetime, t, time
!        end if
        !*******************************************************	
        !calculate urms:						
        !First average the total flow over time to get the mean flow	
        !as a function of (x,y): 
      
        !    uxavg(:,:)=(uxavg(:,:)+vex(:,:))/(p+1.)
        !    uyavg(:,:)=(uyavg(:,:)+vey(:,:))/(p+1.)
      
        !Now, subtract the mean flow from the total flow:
        ux_temp = vex - uxavg 
        uy_temp = vey - uyavg 
      
        ux=0.0
        uy=0.0
        do k=1,nlz_par
           do j=1,nly_par
              do i=1,nlx
                 ux=ux+ux_temp(i,j,k)**2
                 uy=uy+uy_temp(i,j,k)**2
              end do
           end do
        end do
        call sum_reduce(ux,proc_X)
        call sum_reduce(uy,proc_X)
        !call zderivzgk2(gk,dzgk2) !zderivative
        dzgk2(:,:,:) = 0.0
        call bracket_3(dxapar,dyapar,dxg(:,:,:,gmin),dyg(:,:,:,gmin),braakpargk2)
        call k_energy(1,savetime,p,akpar,phik,nek,gk,dzgk2,braakpargk2,ek_u,ek_b,ek_e,hk_a,hk_t,ik_t,ik_a,ek_n,hk_aa, hk_ta, ik_ta, ik_aa,hk_1,hk_2,chk,runname)      
        !call k_energy_pert(1,savetime,p,akpar-akpar_eq,phik,nek,gk,dzgk2,braakpargk2,ek_u,ek_b,ek_e,hk_a,hk_t,ik_t,ik_a,ek_n,hk_aa, hk_ta, ik_ta, ik_aa,hk_1,hk_2,chk,runname)   ! Deprecated 2/28/23 AAV

!        call zonalflows(652,savetime,p,phik(1,:,1),akpar(1,:,1)-akpar_eq(1,:,1),runname)
!        if (iproc==proc_X) then
!           vrms=sqrt(vex(x_loc,y_loc_loc)**2+vey(x_loc,y_loc_loc)**2)
!           vrms2=sqrt(1./(nlx*nly)*(ux+uy))
!        end if

        !TTR
#       ifdef gasca2d
        do k = 1, nlz_par
           call FFT2d_inv ( Akpar(:, :, k),        Apar(:, :, k))
           call FFT2d_inv (uekpar(:, :, k),       uepar(:, :, k))
           call FFT2d_inv (  phiK(:, :, k),         phi(:, :, k))
           call FFT2d_inv ( dummy(:, :, k), dummy_real2(:, :, k))
           call FFT2d_inv (   nek(:, :, k),    ne_real(:, :, k)) ! Z.Liu 6/29/2020
        end do
#       elif defined(gasca3d)
        call FFT2d_inv ( Akpar(:, :, :),        Apar(:, :, :))
        call FFT2d_inv (uekpar(:, :, :),       uepar(:, :, :))
        call FFT2d_inv (  phiK(:, :, :),         phi(:, :, :))
        call FFT2d_inv ( dummy(:, :, :), dummy_real2(:, :, :))
        call FFT2d_inv (   nek(:, :, :),    ne_real(:, :, :)) !Z.Liu 6/29/2020
        call FFT2d_direct (Apar(:, :, :), AKpar(:, :, :))

#       endif

 !       call layer_xy_widths(uekpar, uepar, y_loc_loc, proc_X, delta_x, delta_y)
 !       if(iproc==proc_X) then 
 !          call X_point(Apar,uepar-uepar_eq,y_loc_loc,a_X,u_X)
 !          call layer_width(uepar-uepar_eq,y_loc_loc,delta,p_layer)
 !          !        call alex_sheet_width(uepar-uepar_eq,y_loc_loc,alex_width)
 !          !           call send(delta,0)
 !          !           call broadcast(delta,proc_X)
 !          call send(a_X,0)
 !          !        call send(alex_width,0)
 !       end if
 !       call broadcast(delta,proc_X)
 !       call sheet_length(uepar-uepar_eq,L_sheet)
 !       call sheet_length2(vey,LCS2)
 !       !      call epar_field (dadtk,epar)
 !       !      call average_epar(epar,0.075,0.4,epar_avg)
 !       !    call alex_sheet_length(uepar-uepar_eq,L_sheet_alex)
 !       !    call epar_field(dadtK,AKpar,phiK,Epar,gradphi,dAdt)
 !       !        call b_field(Akpar,bx,by)
 !       call SP(vey,by,p_layer,proc_X,y_loc_loc,outflow_max,Bup_max,Bup_jemella)
 !       if (iproc==proc_X)then
 !          call send(Bup_max,0)
 !          call send(Bup_jemella,0)
 !          width=4*sqrt(abs((Apar(x_loc,y_loc_loc,nlz/2)-Apar_eq(x_loc,y_loc_loc,nlz/2))&
 !               /Apar_eq_double_prime(x_loc,y_loc_loc,nlz/2)))
 !       end if
        !call injected(dummy_real2,phi,inj)
        call injected(dummy_real2,phi,ne_real,inj) !Z.Liu 6/29/2020
        !call injected(dummy_real2,apar,inj) !LMM
        call calc_leftover(dxapar,dyapar,dxg(:,:,:,ngtot),dyg(:,:,:,ngtot),gk(:,:,:,ngtot-1),sumleftover)
        !call diagnostics(phiK, nek, uekpar,uekpar_eq,Akpar,& !Z.Liu 10/15/2020
        !     Akpar_eq,gk,sumleftover,savetime,res2,niu2,nu_g,hyper_nuei,dti,proc_X,inj,energyzloss)  
        !>NFL, 04/11/2013
        !call gamma_ky(savetime,nek,apar-apar_eq,runname)
        !<


        !     call broadcast(width)
        !        call barrier
        if(proc0) then
!           call receive(a_X,proc_X)
!           call receive(Bup_max,proc_X)
!           call receive(Bup_jemella,proc_X)
!           call isl_width(Apar(:,:,nlz/2),proc_X,y_loc_loc,a_X,W2,pos_isl)
!           call send(W2,proc_X) 
!           write(666,50) savetime, Bup_max, Bup_jemella, outflow_max,&
!                & delta,delta_x,L_sheet,delta_y,LCS2,niu2
!50         format(10g16.8)
!           write(667,51) savetime,&
!                & ek_u,ek_b,energy,denergy_dt,av_denergy_dt,lost,epar_avg
!           write(667,51) savetime,&
!                & ek_u,ek_b,epar_avg
           
!51         format(8g16.8)
           write(33,53) savetime, CFL_frac*dx/vxmax, CFL_frac*dy/vymax, &
                CFL_frac*dx/bxmax, CFL_frac*dy/bymax, CFL_frac*2./omega_kaw, &
                CFL_frac*(1./rhos_de)*1./sqrt(ngtot*1.0)*min(dy/bymax,dx/bxmax), &
                CFL_frac*(1./rhos_de)*dz/sqrt(ngtot*1.0), & 
                CFL_frac*dz*(1+max(nky,nkx_par*Npe)*de**2)**0.5/(max(nky,nkx_par*Npe)*rhos), &
                CFL_frac*dz/(rhos_de*sqrt(ngtot + 1.0)),&
                CFL_frac*min(dx/vxmax,dy/vymax,2./omega_kaw, &
                (1./rhos_de)*1./sqrt(ngtot*1.0)*min(dx/bxmax,dy/bymax), &
                (1./rhos_de)*dz/sqrt(ngtot*1.0), &
                dz*(1+max(nky,nkx_par*Npe)*de**2)**0.5/(max(nky,nkx_par*Npe)*rhos),&
                dz/(rhos_de*sqrt(ngtot + 1.0)),dz**2/(zdiffcoeff),1/(zdiffcoeff*kperpmax**2*(max(bxmax,bymax))**2)),dti, &
                RELATIVE_ERROR, abs(akpar(j1,2,1)), &
                CFL_frac*dz**2/(zdiffcoeff), &
                CFL_frac/(zdiffcoeff*kperpmax**2*(max(bxmax,bymax))**2)!p_count*1.0   !p_count*1.0 
53         format(18g16.8)  !,f6.1)

        end if  
        if (iproc==proc_X) then
!           call receive(W2,0)
           write(77,49) savetime,Apar(x_loc,y_loc_loc,1)-apar_eq(x_loc&
                ,y_loc_loc,1), uepar(x_loc,y_loc_loc,1),width,W2,vrms,vrms2 
           !49         format(1p35e13.5,1p35e13.5,1p35e13.5,1p35e13.5,1p35e13.5&
           !                &,1p35e13.5) !,1p35e13.5)
49         format(7g16.8)
        end if

        if(g_inc) call gm_timetrace(savetime,gk(:,:,:,:) )

        if(proc0) then
           open(unit=344, file='exit')
           read(344,'(i1)') keepatit
           close (unit=344)
        end if 
        call broadcast (keepatit)
        if (keepatit /= 0) exit time_loop
     end if !end if of mod(p,save_datafiles)==0 
     if(p .GT. 1  .AND. mod(p,save_energyfiles)==0)THEN

        call predatasave(savetime,howlong)
        !TTR
#       ifdef gasca2d
        do k = 1, nlz_par
           call FFT2d_inv (Akpar(:, :, k), Apar(:, :, k))
           call FFT2d_inv (  neK(:, :, k),   ne(:, :, k)) 
        end do
#       elif defined(gasca3d)
        call FFT2d_inv (Akpar(:, :, :), Apar(:, :, :))
        call FFT2d_inv (  neK(:, :, :),   ne(:, :, :))
        call FFT2d_direct (  ne(:, :, :),   nek(:, :, :))
        call FFT2d_direct (Apar(:, :, :), AKpar(:, :, :))
 
#       endif   
        if(g_inc) then
           call gm_spectrum(1, savetime, gk, runname)
           Myperfon('ngTLoop')
#          ifdef gasca2d
           do ng = gmin, ngtot
              do k = 1, nlz_par
                 call FFT2d_inv (gk(:, :, k, ng), g(:, :, k, ng))
              end do
           end do
#          elif defined(gasca4d)
           call FFT2d_inv (gk(:, :, :, :), g(:, :, :, :))
#          elif defined(gasca3d)
           do ng = gmin, ngtot
              call FFT2d_inv (gk(:, :, :, ng), g(:, :, :, ng))
           end do
#          endif
           Myperfoff ! ngTLoop
           !"call datasavetest(runname, savetime,howlong,apar,ne,epar,g)
         !else
           !call datasavetest(runname, savetime,howlong,apar,ne,epar,gdummy)
        end if

        call cuts(1,savetime, proc_X, y_loc_loc, uepar-uepar_eq, apar-apar_eq, phi, vey, runname)	
     end if
     Myperfoff ! Diags
     if(p .gt. 1  .and. mod(p,save_checkpoints)==0)then 
#if WITH_HAC
     datasavebenchloop: do rdi = 1, rd
#endif
  
        call predatasave(savetime, howlong)

        !TTR
#       ifdef gasca2d
        do k = 1, nlz_par
           call FFT2d_inv (Akpar(:, :, k), Apar(:, :, k))
           call FFT2d_inv (  neK(:, :, k),   ne(:, :, k))
        end do
#       elif defined(gasca3d)
        call FFT2d_inv (Akpar(:, :, :), Apar(:, :, :))
        call FFT2d_inv (  neK(:, :, :),   ne(:, :, :))
        call FFT2d_direct (  ne(:, :, :),   nek(:, :, :))
        call FFT2d_direct (Apar(:, :, :), AKpar(:, :, :))

#       endif
   
        if(g_inc) then
#          ifdef gasca2d
           do ng = gmin, ngtot
              do k = 1, nlz_par
                 call FFT2d_inv (gk(:, :, k, ng), g(:, :, k, ng))
              end do
           end do
#          elif defined(gasca4d)
           call FFT2d_inv (gk(:, :, :, :), g(:, :, :, :))
#          elif defined(gasca3d)
           do ng = gmin, ngtot
              call FFT2d_inv (gk(:, :, :, ng), g(:, :, :, ng))
           end do
#          endif
           call datasavetest(runname, savetime,howlong,apar,ne,epar,g)
           call FFT2d_direct (  ne(:, :, :),   nek(:, :, :))
           call FFT2d_direct (Apar(:, :, :), AKpar(:, :, :))
           !call kdatasave3d(runname, savetime, howlong, akpar, nek, phik)
           call init_transforms !LMM
           if (proc0) then
                 write(1408, '(4E11.4)') res2, &
                   &niu2, nu_g, hyper_nuei
           endif
        else
           call datasavetest(runname, savetime,howlong,apar,ne,epar,gdummy)
           call FFT2d_direct (  ne(:, :, :),   nek(:, :, :))
           call FFT2d_direct (Apar(:, :, :), AKpar(:, :, :))
           !call kdatasave3d(runname, savetime, howlong, akpar, nek, phik)
           if (proc0) then
                 write(1408, '(4E11.4)') res2, &
                   &niu2, nu_g, hyper_nuei
           endif
           call init_transforms !LMM
        end if
#if WITH_HAC
     end do datasavebenchloop
#endif
     end if
 
     !CALCULATE NEW TIMESTEP
     if (turb) then
        dti=x
     else
        call dtnext(RELATIVE_ERROR,x,noinc,dti)
     end if
!     call hyper_diff (bperp_max,omega_kaw,niu2,res2)
     if (hyper_fixed) then 
        nu_g=hyper_nu_g
        niu2=hyper_nu
        res2=hyper_eta
     else
!        nu_g=hyper_coef_g/dti/(NLx/2.)**(2*hyper_order_g)
        nu_g=hyper_coef_g/dti/(kperpmax)**(2*hyper_order_g)
!        niu2=hyper_coef/dti/(NLx/2.)**(2*hyper_order)
        niu2=hyper_coef/dti/(kperpmax)**(2*hyper_order)
!        if ((nkx_par*npe)**2*de**2>1) then
        if (kperpmax**2*de**2>1) then
!           res2=hyper_coef/dti/(NLx/2.)**(2*hyper_order-2)*de**2
           res2=hyper_coef/dti/(kperpmax)**(2*hyper_order-2)*de**2
        else
!        res2=hyper_coef/dti/(NLx/2.)**(2*hyper_order)
           res2=hyper_coef/dti/(kperpmax)**(2*hyper_order)
        end if
     end if

     if (hyper_colls_fixed) then
        hyper_nuei=hyper_colls
     else
        hyper_nuei=hyperm_coef/dti/(Ngtot+1)**(2*hyper_morder)
     end if

     If(mod(p,save_datafiles)==0)THEN
         call diagnostics(phiK, nek, uekpar,uekpar_eq,Akpar,&
             Akpar_eq,gk,sumleftover,savetime,res2,niu2,nu_g,hyper_nuei,dti,proc_X,inj,energyzloss)
     end if
  
  end do time_loop		!time loop
#if WITH_HAC
  lsave_et=hv_tstime
#endif
  Myperfoff ! TimeLoop

  if(proc0) then   ! to see how long the code takes    
     call cpu_time(time2)
     write(*,*) 'Total time = ', time2-time1
  endif

  if (flagout) then ! don't write more output; not needed, for timings

     !call predatasave(savetime,howlong)
     !TTR
#    ifdef gasca2d
     do k = 1, nlz_par
        call FFT2d_inv (Akpar(:, :, k), Apar(:, :, k))
        call FFT2d_inv (  neK(:, :, k),   ne(:, :, k)) 
     end do
#    elif defined(gasca3d)
     call FFT2d_inv (Akpar(:, :, :), Apar(:, :, :))
     call FFT2d_inv (  neK(:, :, :),   ne(:, :, :)) 
     call FFT2d_direct (  ne(:, :, :),   nek(:, :, :))
     call FFT2d_direct (Apar(:, :, :), AKpar(:, :, :))
#    endif
     if(g_inc) then
        call gm_timetrace(savetime,gk(:,:,:,:))
        call gm_spectrum(1,savetime,gk,runname)
#       ifdef gasca2d
        do ng = gmin, ngtot
           do k = 1, nlz_par
              call FFT2d_inv (gk(:, :, k, ng), g(:, :, k, ng))
           end do
        end do
#       elif defined(gasca4d)
        call FFT2d_inv (gk(:, :, :, :), g(:, :, :, :))
#       elif defined(gasca3d)
        do ng = gmin, ngtot
           call FFT2d_inv (gk(:, :, :, ng), g(:, :, :, ng))
        end do
#       endif
        call datasavetest(runname, savetime, howlong, apar, ne, epar, g)
        call FFT2d_direct (  ne(:, :, :),   nek(:, :, :))
        call FFT2d_direct (Apar(:, :, :), AKpar(:, :, :))
        !call kdatasave3d(runname, savetime, howlong, akpar, nek, phik)
        call init_transforms !LMM
        if (proc0) then
                 write(1408, '(4E11.4)') res2, &
                   &niu2, nu_g, hyper_nuei
           endif
     else
        call datasavetest(runname, savetime, howlong, apar, ne, epar, gdummy)
        call FFT2d_direct (  ne(:, :, :),   nek(:, :, :))
        call FFT2d_direct (Apar(:, :, :), AKpar(:, :, :))
        !call kdatasave3d(runname, savetime, howlong, akpar, nek, phik)
        call init_transforms !LMM 
        if (proc0) then
         write(1408, '(4E11.4)') res2, &
                   &niu2, nu_g, hyper_nuei
           endif
     end if
     call cuts(1,savetime,proc_X,y_loc_loc,uepar-uepar_eq, apar-apar_eq, phi, vey, runname)
     call k_energy(1,savetime,p,akpar,phik,nek,gk,dzgk2,braakpargk2,ek_u,ek_b,ek_e,hk_a,hk_t,ik_t,ik_a,ek_n,hk_aa, hk_ta, ik_ta, ik_aa,hk_1,hk_2,chk,runname)
!     call zonalflows(652,savetime,p,phik(1,:,1),akpar(1,:,1)-akpar_eq(1,:,1),runname)

  endif
 
  !  close (unit=344)
  if (proc0) then
     !     close (unit=2020)
     close (unit=33)
     close (Unit=666)
     close (unit=221)
     close (unit=99)
     close (unit=1729)
     close (unit=1408)
!     close (Unit=44) 
!     close (Unit=46) 
  end if
  if (iproc==proc_X) then
     close(unit=77)
     close(unit=88)
  end if

  !. Perflib de-initialization
  Myperfoff ! main
  do ip = 0, 0!nproc-1
     !call barrier
     if (iproc .eq. ip) then
        write(*,'(A,I2)')'Perfout PE', iproc
        Myperfout('Main')
     end if
  end do

999  continue
#if WITH_HAC
  if_hacexit: if( use_hac_checkpoint ) then
    if( hv_vv /= HAC_QUIET ) then
      call hac_stats(ogwb_=ogwb,ogrb_=ogrb)
      save_et=(hv_tstime)/(rd*save_checkpoints)
      cstep_et=((time2-time1-lsave_et)/(tmax+1))
      if( proc0 ) &
      &write(*,'(3(a,i6),5(a,i3),8(a,es9.2))') hv_lbl//&
        &" nproc: ", (nproc),&
        &" csteps_n: ", (tmax+1),&
        &" saves_n: ", (rd*save_checkpoints), &
        &" nlx: ", nlx, &
        &" nly: ", nly, &
        &" nlz: ", nlz, &
        &" gmin: ", gmin, &
        &" ngtot: ", ngtot, &
        & &
        &" saved_bytes: ", ogwb,&
        &" loaded_bytes: ", ogrb, &
        &" tot_save_et: ", (hv_tstime),&
        &" init_et: ", (time1),&
        &" cstep_et: ", cstep_et, &
        &" save_et: ", save_et,&
        &" save_to_cstep_et_ratio: ", save_et/MAX(cstep_et,TINY(cstep_et)),&
        &" loop_et: ", (time2-time1)
    end if
    CALL hac_exit ()
  end if if_hacexit
#endif
  call finish_mp
  stop

end program REGK


!********** end of main section of code**************************



