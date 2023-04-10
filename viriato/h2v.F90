!
!  This is a "HAC to VIRIATO" conversion example program.
!
! Given a run name (oldrun), will read it and the associated checkpoint file.
! It can be run from either the original number of npe*npez MPI processes or 1.
! In both cases it can process the data and save it in the same format.
!
! Invocation: ./h2v <oldrun> [dowhat=dodefault]
!
! Jan-Mar 2016: Michele MARTONE (MMH):
!
      PROGRAM H2V
#if WITH_HAC
        use mpi,       only: MPI_Abort, MPI_COMM_WORLD
        use mp,        only: barrier, proc0, iproc, nproc, init_mp, finish_mp, sum_reduce, min_reduce, max_reduce
        use diag,      only: dataloadtest_hac, datasavetest_hac, hv_lbl, hv_vv !< for hac
        use constants, only: nlz, read_parameters, gmin, ngtot, nlx, nly, nly_par, nlz_par
        use constants, only: set_anjor, anjor, rhos_de, use_hac_checkpoint, oldrun, npe, npez, g_inc
        use hlst_adios_checkpoint, only: hac_init,hac_stats,hac_exit,HAC_QUIET,HAC_VERBOSE !< for hac
        use grid,       only: init_grid
!
        implicit none
        character(len=100) :: ipfile, newrun
        character(len=100) :: dowhat
        character(len=*), parameter :: dodefault='imcwzrcw'
        real, allocatable, dimension(:,:,:) :: Apar, Apar_buff
        real, allocatable, dimension(:,:,:) :: ne, ne_buff
        real, allocatable, dimension(:,:,:) :: Epar
        real, allocatable, dimension(:,:,:,:) :: g, gdummy
        real :: savetime
        real :: asv, gsv, miv, mav ! absolute sum, sum, minimal, maximal  value
        integer :: howlong
        integer,save :: dwi, rep, repi, rproc
        integer*4,save :: ferr=-1, ierr=0
!
10      format (10(f10.0))
#if WITH_HAC
        call getarg(1, oldrun)
        call getarg(2, dowhat)
        ipfile = trim(oldrun)//".in"
!
        if ( len(trim(oldrun)) == 0 ) then
          if ( iproc == 0 ) write(*,'(a)') &
           hv_lbl//'No oldrun string specified, terminating !'
          goto 9999    
        end if
!
        call read_parameters(ipfile)
        ! TODO: we need a test case scenario here.
        hv_vv = HAC_VERBOSE
        call set_anjor(anjor,rhos_de,gmin) ! gmin
        call init_mp
        call barrier
        if ( iproc == 0 ) then
          write(*,'(a,i0,"/",i0,", ",i0,"x",i0)') &
          & hv_lbl//'Running converter on ',iproc,nproc,npe,npez
          write(*,'(a,3(a," ",i0," "))') &
          & hv_lbl//'Dims: ','nlx:',nlx,'nly:',nly,'nlz:',nlz
        end if
        call barrier
        if ( nproc /= 1 .and. nproc /= npe*npez ) then
          if ( iproc == 0 ) write(*,'(a,i0," /= ",i0,": ABORTING!")') &
          & hv_lbl//'Running on non-matching proc count: ',npe*npez,nproc
          CALL MPI_Abort(MPI_COMM_WORLD, ferr, ierr)
        end if
        call barrier
        if ( nproc == 1 ) then
          allocate(                apar(nlx, nly, nlz))
          allocate(           apar_buff(nlx, nly, nlz))
          allocate(                  ne(nlx, nly, nlz))
          allocate (ne_buff(nlx,nly,nlz))
          allocate (       g(nlx, nly, nlz, gmin:ngtot))
          allocate (  gdummy(nlx, nly, nlz, gmin:ngtot))
          allocate(                epar(nlx, nly, nlz))
        else
          allocate(                apar(nlx, nly_par, nlz_par))
          allocate(           apar_buff(nlx, nly_par, nlz_par))
          allocate(                  ne(nlx, nly_par, nlz_par))
          allocate (ne_buff(nlx,nly_par,nlz_par))
          allocate (       g(nlx, nly_par, nlz_par, gmin:ngtot))
          allocate (  gdummy(nlx, nly_par, nlz_par, gmin:ngtot))
          allocate(                epar(nlx, nly_par, nlz_par))
        end if
!
        call init_grid
        apar=0.0
        apar_buff=0.0
        ne=0.0
        ne_buff=0.0
        g=0.0
        gdummy=0.0
        epar=0.0
        call barrier
        if ( nproc == 1 ) then
          write(*,'(a,3(i0," "),i0,":",i0)') &
          & hv_lbl//'Read the entire array:',nlx,nly,nlz,gmin,ngtot
        else
          if (proc0) then
            write(*,'(a,i0,"x",i0)') hv_lbl//&
            & 'Assuming running with same task count as original run: ',npe,npez
            write(*,'(a,3(i0," "),i0,":",i0)') &
            & hv_lbl//'Read the local slice:',nlx,nly_par,nlz_par,gmin,ngtot
          end if
        end if
!
        call barrier
        call hac_init (MPI_COMM_WORLD, hv_vv )
!
        if ( len(trim(dowhat)) == 0 ) then
          dowhat=dodefault
          if ( iproc == 0 ) write(*,'(a)') &
           & hv_lbl//'No command string specified, so setting to "'//trim(dowhat)//'"'
        else
          if ( iproc == 0 ) write(*,'(a)') &
           & hv_lbl//'User specified command string is "'//trim(dowhat)//'"'
        end if
!
        rep=0
        do dwi=1, len(trim(dowhat))
        select case(dowhat(dwi:dwi))
        case('0','1','2','3','4','5','6','7','8','9')
         rep=10*rep+(iachar(dowhat(dwi:dwi))-iachar('0'))
         cycle
        end select
!
        if(rep == 0) rep = 1
!
        if ( iproc == 0 ) write(*,'(a)') &
        & hv_lbl//'Entering I/O demo command cycle.'
        do repi=1, rep
        if ( iproc == 0 .and. rep .gt. 1 ) &
         & write(*,'(a,i0,a,i0)') &
         & hv_lbl//'Executing command: Repeat:',repi,'/',rep
        select case(dowhat(dwi:dwi))
        case('l', 'r')
         if ( iproc == 0 ) write(*,'(a)') &
         & hv_lbl//'Executing command: Read.'
         call dataloadtest_hac(oldrun, apar,Apar_buff,ne,ne_buff,g,gdummy,epar)
        case('d')
         if ( iproc == 0 ) write(*,'(a)') &
         & hv_lbl//'Executing command: Double.'
         g=g*2
         ne=ne*2
         apar=apar*2
         epar=epar*2
        case('z')
         if ( iproc == 0 ) write(*,'(a)') &
         & hv_lbl//'Executing command: Zero.'
         g=0
         ne=0
         apar=0
         epar=0
        case('i')
         if ( iproc == 0 ) write(*,'(a)') &
         & hv_lbl//'Executing command: Init.'
         CALL hv_lai4(   g,iproc)
         CALL hv_lai3(  ne,iproc)
         CALL hv_lai3(apar,iproc)
         CALL hv_lai3(epar,iproc)
        case('p')
         if ( iproc == 0 ) write(*,'(a)') &
         & hv_lbl//'Executing command: Print.'
         if ( PRODUCT(SHAPE(g))*nproc < 500 ) then
          print *, "g array on proc ",iproc," :"
          print 10, g
         endif
        case('m')
         if ( iproc == 0 ) write(*,'(a)') &
         & hv_lbl//'Executing command: Modify.'
         if (proc0)write(*,'(a)')hv_lbl//'Sample processing of data ...'
         g=g+1
        case('w', 's')
         if ( iproc == 0 ) write(*,'(a)') &
         & hv_lbl//'Executing command: Write.'
         if (proc0) write(*,'(a)')hv_lbl//'Saving preprocessed data ...'
         if (proc0) write(*,'(a)') &
           & hv_lbl//'NOTE: overwriting old checkpoint ...'
           newrun = trim(oldrun)
           call datasavetest_hac(newrun , savetime, howlong, apar, ne, epar, g)
        case('c')
         if(g_inc) then
          asv=sum(abs(g))+sum(abs(apar))+sum(abs(ne))
          gsv=sum(g)+sum(apar)+sum(ne)
          miv=min(minval(g),minval(apar),minval(ne))
          mav=max(maxval(g),maxval(apar),maxval(ne))
         else
          asv=sum(abs(epar))+sum(abs(apar))+sum(abs(ne))
          gsv=sum(epar)+sum(apar)+sum(ne)
          miv=min(minval(epar),minval(apar),minval(ne))
          mav=max(maxval(epar),maxval(apar),maxval(ne))
         endif
         call sum_reduce(asv, rproc)
         call sum_reduce(gsv, rproc)
         call min_reduce(miv, rproc)
         call max_reduce(mav, rproc)
         if (proc0) write(*,'(4(a,e20.10))') hv_lbl//&
          &'Global  sum: ',gsv, &
          &'  sum of abs: ',asv, &
          &'  min: ',miv,'  max: ',mav
        case default
         if (proc0) write(*,'(a)') &
          &hv_lbl//'Command char '//dowhat(dwi:dwi)//&
          &' unknown: skipping it.'
        end select
        end do
        rep=0
        end do
        if ( iproc == 0 ) write(*,'(a)') &
        & hv_lbl//'Leaving I/O demo command cycle.'
!
        CALL hac_exit ()
#else
        write(*,*) "No HAC: no converter functionality possible."
        stop -1
#endif
! WITH_HAC
        call finish_mp
9999    stop
      CONTAINS
       SUBROUTINE hv_lai3(AS,erank)
        ! Local Array Initialize
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: erank !
        REAL,DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT) :: AS !< Array Slice
        INTEGER :: LB(3), UB(3) !< Lower/Upper Bounds
        INTEGER,PARAMETER :: MF = 100
        INTEGER :: i1, i2, i3
        LB=LBOUND(AS)
        UB=UBOUND(AS)
        AS =        10**8 * (erank + 1)
        FORALL(i1=LB(1):UB(1),i2=LB(2):UB(2),i3=LB(3):UB(3))&
                        &AS(i1,i2,i3) = i1+MF*i2+i3*MF**2
       END SUBROUTINE hv_lai3
       SUBROUTINE hv_lai4(AS,erank)
        ! Local Array Initialize
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: erank !
        REAL,DIMENSION(:,:,:,:),ALLOCATABLE,INTENT(INOUT) :: AS !< Array Slice
        INTEGER :: LB(4), UB(4) !< Lower/Upper Bounds
        INTEGER,PARAMETER :: MF = 100
        INTEGER :: i1, i2, i3, i4
        LB=LBOUND(AS)
        UB=UBOUND(AS)
        AS =        10**8 * (erank + 1)
        FORALL(i1=LB(1):UB(1),i2=LB(2):UB(2),i3=LB(3):UB(3),i4=LB(4):UB(4))&
                        &AS(i1,i2,i3,i4) = i1+MF*i2+i3*MF**2+MF**3*i4
       END SUBROUTINE hv_lai4
#else
        write(*,*) "No HAC (HLST-ADIOS-Checkpoint) enabled."
#endif
      END PROGRAM H2V
!
