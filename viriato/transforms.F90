module transforms

  !TTR
  !. FPP macro to simplify RZG PERFLIB calls (performance-counter)
# include "perfmacro.h"

  use constants
  use fft_work, only: fft_type, POINTER_KIND
  use redistribute, only: redist_type

  implicit none

  private

  !TTR
  type (redist_type), public, save :: r2k3, r2k4
  type (fft_type), public :: fft_x2k, fft_k2x, fft_y2k, fft_k2y!, fft_z2k

  integer(Kind=POINTER_KIND) :: iplan_x2k, iplan_k2x, iplan_y2k, iplan_k2y!, iplan_z2k

  public :: init_transforms, oneDfourninv

  !deprecated: replaced by FFT2d_[*]2. Can be REMOVED!
  public :: FFT1d_direct, FFT1d_inv
  public :: FFT2d_direct, FFT2d_inv

  interface FFT2d_direct
     module procedure FFT2d_direct2
     module procedure FFT2d_direct3
     module procedure FFT2d_direct4
  end interface

  interface FFT2d_inv
     module procedure FFT2d_inv2
     module procedure FFT2d_inv3
     module procedure FFT2d_inv4
  end interface

contains
  

!*********************************************************  
  subroutine init_transforms
  
    implicit none

    !TTR
    call init_redistribute3
    call init_redistribute4

  end subroutine init_transforms


!TTR
!> this subroutine completely replaces the original "init_redistribute"
!> it can be used together with the original FFT2D (2D) method "FFT1d_[*]", 
!> as well as with the new 2D and 3D versions of it "FFT2d_[*]2" and "FFT2d_[*]3"
!*********************************************************  
  subroutine init_redistribute3

    use grid
    use mp
    use redistribute, only: index_list_type, init_redist, delete_list

    implicit none

    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension (3) :: from_low, to_low, from_high, to_high
    logical :: initialized = .false.
    integer :: i, j, ip, n, k

  !  if (initialized) return !LMM
  !  initialized = .true.

    ! count number of elements to be redistributed to/from each processor

    nn_to = 0
    nn_from = 0
    do k = 1, nlz
       do j = 1, NLy
          do i = 1, Nkx
             if (idx_local(r_variable, j, k)) &
                  nn_from(proc_id(k_variable, i, k)) = nn_from(proc_id(k_variable, i, k)) + 1
             if (idx_local(k_variable, i, k)) &
                  nn_to(proc_id(r_variable, j, k)) = nn_to(proc_id(r_variable, j, k)) + 1
          end do
       end do
    end do

    ! allocate lists to hold local indices of elements distributed to/from other processors

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first (nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third (nn_from(ip))) !.3
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first (nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
          allocate (to_list(ip)%third (nn_to(ip))) !.3
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    !.repeated calculation of nn_to & nn_from (?)

    nn_to = 0
    nn_from = 0
    do k = 1, nlz
       do j = 1, NLy
          do i = 1, Nkx   ! this is dealiasing in x
             if (idx_local(r_variable, j, k)) then
                ip = proc_id(k_variable, i, k)
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n)  = i
                from_list(ip)%second(n) = 1+mod(j-1, NLy_par)
                from_list(ip)%third(n)  = 1+mod(k-1, nlz_par)
             end if
             if (idx_local(k_variable, i, k)) then
                ip = proc_id(r_variable, j, k)
                n = nn_to(ip) + 1
                nn_to(ip) = n
                ! local indices of elements to be sent in the plane (kx_par, y)
                to_list(ip)%first(n)  = j
                to_list(ip)%second(n) = 1+mod(i-1, NKx_par)
                to_list(ip)%third(n)  = 1+mod(k-1, nlz_par)
             end if
          end do
       end do
    end do

    from_low(1) = 1
    from_low(2) = 1
    from_low(3) = 1 !.3

    to_low(1) = 1
    to_low(2) = 1
    to_low(3) = 1 !.3

    to_high(1) = NLy
    to_high(2) = NKx_par
    to_high(3) = nlz_par !.3

    from_high(1) = NLx/2+1
    from_high(2) = NLy_par
    from_high(3) = nlz_par !.3

    call init_redist (r2k3, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

  end subroutine init_redistribute3


!TTR
!> this subroutine initialises the original FFT method but carrying 
!> the z- and ng-dimensions (4D)
!> To be used together with the original FFT2D (4D) method "FFT2d_[*]4 only"
!*********************************************************  
  subroutine init_redistribute4

    use grid
    use mp
    use redistribute, only: index_list_type, init_redist, delete_list

    implicit none

    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    !.4 integer, dimension (2) :: from_low, to_low, from_high, to_high
    integer, dimension (4) :: from_low, to_low, from_high, to_high
    logical :: initialized = .false.
    integer :: i, j, ip, n, k, ng

   ! if (initialized) return
   ! initialized = .true.

    ! count number of elements to be redistributed to/from each processor

    nn_to = 0
    nn_from = 0

    do ng = gmin, ngtot
    do k = 1, nlz
       do j = 1, NLy
          do i = 1, Nkx
             if (idx_local(r_variable, j, k)) &
                  nn_from(proc_id(k_variable, i, k)) = nn_from(proc_id(k_variable, i, k)) + 1
             if (idx_local(k_variable, i, k)) &
                  nn_to(proc_id(r_variable, j, k)) = nn_to(proc_id(r_variable, j, k)) + 1
          end do
       end do
    end do
    end do

    ! allocate lists to hold local indices of elements distributed to/from other processors

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_list(ip)%first (nn_from(ip)))
          allocate (from_list(ip)%second(nn_from(ip)))
          allocate (from_list(ip)%third (nn_from(ip))) !.3
          allocate (from_list(ip)%fourth(nn_from(ip))) !.4
       end if
       if (nn_to(ip) > 0) then
          allocate (to_list(ip)%first (nn_to(ip)))
          allocate (to_list(ip)%second(nn_to(ip)))
          allocate (to_list(ip)%third (nn_to(ip))) !.3
          allocate (to_list(ip)%fourth(nn_to(ip))) !.4
       end if
    end do

    ! get local indices of elements distributed to/from other processors
    !.repeated calculation of nn_to & nn_from (?)

    nn_to = 0
    nn_from = 0
    do ng = gmin, ngtot
    do k = 1, nlz
       do j = 1, NLy
          do i = 1, Nkx   ! this is dealiasing in x
             if (idx_local(r_variable, j, k)) then
                ip = proc_id(k_variable, i, k)
                n = nn_from(ip) + 1
                nn_from(ip) = n
                from_list(ip)%first(n)  = i
                from_list(ip)%second(n) = 1+mod(j-1, NLy_par)
                from_list(ip)%third(n)  = 1+mod(k-1, nlz_par) !.3
                from_list(ip)%fourth(n) = ng!+1                !.4
             end if
             if (idx_local(k_variable, i, k)) then
                ip = proc_id(r_variable, j, k)
                n = nn_to(ip) + 1
                nn_to(ip) = n
                ! local indices of elements to be sent in the plane (kx_par, y)
                to_list(ip)%first(n)  = j
                to_list(ip)%second(n) = 1+mod(i-1, NKx_par)
                to_list(ip)%third(n)  = 1+mod(k-1, nlz_par) !.3
                to_list(ip)%fourth(n) = ng!+1                !.4
             end if
          end do
       end do
    end do
    end do

    from_low(1) = 1
    from_low(2) = 1
    from_low(3) = 1    !.3
    from_low(4) = gmin !.4

    to_low(1) = 1
    to_low(2) = 1
    to_low(3) = 1    !.3
    to_low(4) = gmin !.4

    to_high(1) = NLy
    to_high(2) = NKx_par
    to_high(3) = nlz_par !.3
    to_high(4) = ngtot   !.4

    from_high(1) = NLx/2+1
    from_high(2) = NLy_par
    from_high(3) = nlz_par !.3
    from_high(4) = ngtot   !.4

    call init_redist (r2k4, 'c', to_low, to_high, to_list, &
         from_low, from_high, from_list)

    call delete_list (to_list)
    call delete_list (from_list)

  end subroutine init_redistribute4


!deprecated: replaced by FFT2d_direct2. This subroutine can be REMOVED!
!> 2D FFT forward (real, 2D)
!*********************************************************  
  subroutine FFT1d_direct(array, arrayk)

    !... transform to k space

    use redistribute, only: gather
    use grid, only: ky!, kx,kperp

    implicit none
!   
! example case: "array" has dimensions
!    NLx = 1024    <- nlx
!    NLy = 128     <- nly
!
! and "arrayk" has dimensions
!    NKy = 85      <- nky
!    NKx = 342     <- nkx
!
! After x transform, resulting array ("array_temp") is of size
!     kx = 513     <- NLx/2+1
!     ly = 128     <- nly
!
! which is then dealiased down to 
!    Nkx = 342     <- nkx
!     ly = 128     <- nly
!
! at the same time that the array is transposed into 
! arrayk, which has dimensions
!
!     ly = 128
!    Nkx = 342
! 
! after y transform, resulting array is of size 
!     ky = 128
!    Nkx = 342
!
! which is then dealiased down to 
!    Nky = 85 
!    Nkx = 342
!
! in the "map" routine
!
    real,    intent(in),  dimension(:, :) :: array   ! size is nlx * nly_par
    complex, intent(out), dimension(:, :) :: arrayk  ! size is nky * nkx_par

    complex, allocatable, dimension(:, :) :: array_temp, ak
    integer :: i, j!, iglobal, ip

#   ifdef perf
    logical, save :: first = .true.
#   endif

    Myperfon('2DFT_F1')

    allocate( array_temp(nlx/2+1, nly_par) )
    allocate( ak(nly, nkx_par) )

    array_temp = cmplx(0.0, 0.0)
    ak = cmplx(0.0, 0.0)

    !... FFT in x-direction (forward)
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Ff1_x')
    end if
#   endif
    call Nfourndirect_X(array, array_temp)
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Ff1_x
    end if
#   endif

    do i = 1, nkx
       array_temp(i,:)=array_temp(i,:)*exp(-36.0*((i*1.0-1.0)/((nkx-1)*1.0))**36) !Hou-Li expression
    end do
    !
    !... data transposition xy
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Ff1_gath')
    end if
#   endif
    call gather(r2k3, array_temp, ak)
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Ff1_gath
    end if
#   endif

!this is dealising. The y-direction delaising is done in the 
!map subroutine.
   
!    do i=1, nkx_par
!       iglobal=i+mod(iproc,NPE)*nkx_par
!       if (iglobal <= nkx)then
!          ak(:,i)=ak(:,i)
!       else
!          ak(:,i)=0.0
!       end if
!    end do
 !   end if

    !
    !... FFT in y-direction (forward)
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Ff1_y')
    end if
#   endif
    call Nfourndirect_Y(ak)

!    call map(ak, arrayk)

    !NFL 29/05/2013: Hou-Li filter:  
    !I could instead do rho(kx)rho(ky) (where rho is the filter function).

    do j = 1, nky
       arrayk(j, :) = ak(j, :) * exp( -36.0*(abs(ky(j))/(nky/2.*lx/ly))**36 )
    end do
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Ff1_y
    end if
    if (first) first = .false.
#   endif

!    do i=1, nkx_par 
!       do j=1,nky
!          arrayk(j,i)=ak(j,i)*exp(-36.0*&
!               (kperp(j,i)/sqrt((nky/2.*lx/ly)**2+(nlx/2.)**2))**(36.0))
!       end do
!    end do

    deallocate (ak, array_temp)

    Myperfoff ! 2DFT_F1

  end subroutine FFT1d_direct


!> 2D FFT forward (real, 2D)
!*********************************************************  
  subroutine FFT2d_direct2(array, arrayk)

    !... transform to k space

    use redistribute, only: gather
    use grid, only: ky

    implicit none
!   
! example case: "array" has dimensions
!    NLx = 1024    <- nlx
!    NLy = 128     <- nly
!
! and "arrayk" has dimensions
!    NKy = 85      <- nky
!    NKx = 342     <- nkx
!
! After x transform, resulting array ("array_temp") is of size
!     kx = 513     <- NLx/2+1
!     ly = 128     <- nly
!
! which is then dealiased down to 
!    Nkx = 342     <- nkx
!     ly = 128     <- nly
!
! at the same time that the array is transposed into 
! arrayk, which has dimensions
!
!     ly = 128
!    Nkx = 342
! 
! after y transform, resulting array is of size 
!     ky = 128
!    Nkx = 342
!
! which is then dealiased down to 
!    Nky = 85 
!    Nkx = 342
!
! in the "map" routine
!
    real,    intent(in),  dimension(:, :) :: array   ! size is nlx * nly_par
    complex, intent(out), dimension(:, :) :: arrayk  ! size is nky * nkx_par

    complex, allocatable, dimension(:, :) :: array_temp, ak
    integer :: i, j

#   ifdef perf
    logical, save :: first = .true.
#   endif

    Myperfon('2DFT_Fo2')

    allocate( array_temp(nlx/2+1, nly_par) )
    allocate( ak        (nly,     nkx_par) )

    array_temp = cmplx(0.0, 0.0)
    ak = cmplx(0.0, 0.0)

    !... FFT in x-direction (forward)
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Ffo2_x')
    end if
#   endif
    call Nfourndirect_X(array, array_temp)
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Ffo2_x
    end if
#   endif

    do i = 1, nkx
       array_temp(i, :) = array_temp(i, :) &
            &           * exp(-36.0*((i*1.0-1.0)/((nkx-1)*1.0))**36) !Hou-Li expression
    end do
    !
    !... data transposition xy
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Ffo2_gath')
    end if
#   endif
    call gather(r2k3, array_temp, ak)
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Ffo2_gath
    end if
#   endif

!this is dealising. The y-direction delaising is done in the 
!map subroutine.
   
    !do i=1, nkx_par
    !   iglobal=i+mod(iproc,NPE)*nkx_par
    !   if (iglobal <= nkx)then
    !      ak(:,i)=ak(:,i)
    !   else
    !      ak(:,i)=0.0
    !   end if
    !end do
    !   end if

    !
    !... FFT in y-direction (forward)
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Ffo2_y')
    end if
#   endif
    call Nfourndirect_Y(ak)

!    call map(ak, arrayk)

    !NFL 29/05/2013: Hou-Li filter:  
    !I could instead do rho(kx)rho(ky) (where rho is the filter function).

    do j = 1, nky
       arrayk(j, :) = ak(j, :) * exp( -36.0*(abs(ky(j))/(nky/2.*lx/ly))**36 )
    end do
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Ffo2_y
    end if
    if (first) first = .false.
#   endif

    !do i=1, nkx_par 
    !   do j=1,nky
    !      arrayk(j,i)=ak(j,i)*exp(-36.0*&
    !           (kperp(j,i)/sqrt((nky/2.*lx/ly)**2+(nlx/2.)**2))**(36.0))
    !   end do
    !end do

    deallocate (ak, array_temp)

    Myperfoff ! 2DFT_Fo2

  end subroutine FFT2d_direct2


!TTR
!> 2D FFT forward (real, 3D)
!*********************************************************  
  subroutine FFT2d_direct3(array, arrayk)

    !... transform to k space

    use redistribute, only: gather
    use grid, only: ky

    implicit none
!   
! example case: "array" has dimensions
!    NLx = 1024    <- nlx
!    NLy = 128     <- nly
!
! and "arrayk" has dimensions
!    NKy = 85      <- nky
!    NKx = 342     <- nkx
!
! After x transform, resulting array ("array_temp") is of size
!     kx = 513     <- NLx/2+1
!     ly = 128     <- nly
!
! which is then dealiased down to 
!    Nkx = 342     <- nkx
!     ly = 128     <- nly
!
! at the same time that the array is transposed into 
! arrayk, which has dimensions
!
!     ly = 128
!    Nkx = 342
! 
! after y transform, resulting array is of size 
!     ky = 128
!    Nkx = 342
!
! which is then dealiased down to 
!    Nky = 85 
!    Nkx = 342
!
! in the "map" routine
!
    !. size is nlx * nly_par * nlz_par
    real,    intent(in),  dimension(:, :, :) :: array
    !. size is nky * nkx_par * nlz_par
    complex, intent(out), dimension(:, :, :) :: arrayk

    complex, allocatable, dimension(:, :, :) :: array_temp, ak

    integer :: i, j, k

#   ifdef perf
    logical, save :: first = .true.
#   endif
    Myperfon('2DFT_Fo3')

    allocate( array_temp(nkx, nly_par, nlz_par) )
    allocate( ak        (nly, nkx_par, nlz_par) )

    array_temp = cmplx(0.0, 0.0)
    ak = cmplx(0.0, 0.0)

    !... FFT in x-direction (forward)
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Ffo3_x')
    end if
#   endif
    do k = 1, nlz_par
       call Nfourndirect_X(array(:, :, k), array_temp(:, :, k))
    end do
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Ffo3_x
    end if
#   endif

    do k = 1, nlz_par
       do i = 1, nkx
          array_temp(i, :, k) = array_temp(i, :, k) &
               &              * exp(-36.0*((i*1.0-1.0)/((nkx-1)*1.0))**36) !Hou-Li expression
       end do
    end do
    !
    !... data transposition xy
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Ffo3_gath')
    end if
#   endif
    call gather(r2k3, array_temp, ak)
    !
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Ffo3_gath
    end if
#   endif

!this is dealising. The y-direction delaising is done in the 
!map subroutine.
   
    !do i=1, nkx_par
    !   iglobal=i+mod(iproc,NPE)*nkx_par
    !   if (iglobal <= nkx)then
    !      ak(:,i)=ak(:,i)
    !   else
    !      ak(:,i)=0.0
    !   end if
    !end do
    !   end if

    !
    !... FFT in y-direction (forward)
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Ffo3_y')
    end if
#   endif
    do k = 1, nlz_par
       call Nfourndirect_Y(ak(:, :, k))

!       call map(ak(:, :, k), arrayk(:, :, k))

       !NFL 29/05/2013: Hou-Li filter:  
       !I could instead do rho(kx)rho(ky) (where rho is the filter function).

       do j = 1, nky
          arrayk(j, :, k) = ak(j, :, k) * exp( -36.0*(abs(ky(j))/(nky/2.*lx/ly))**36 )
       end do
    end do
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Ffo3_y
    end if
    if (first) first = .false.
#   endif

    !do k = 1, nlz_par
    !do i=1, nkx_par 
    !   do j=1,nky
    !      arrayk(j,i, k)=ak(j,i, k)*exp(-36.0*&
    !           (kperp(j,i)/sqrt((nky/2.*lx/ly)**2+(nlx/2.)**2))**(36.0))
    !   end do
    !end do
    !end do

    deallocate (ak, array_temp)

    Myperfoff ! 2DFT_Fo3

  end subroutine FFT2d_direct3


!TTR
!> 2D FFT forward (real, 4D)
!*********************************************************  
  subroutine FFT2d_direct4(array, arrayk)

    !... transform to k space

    use redistribute, only: gather
    use grid, only: ky

    implicit none
!   
! example case: "array" has dimensions
!    NLx = 1024    <- nlx
!    NLy = 128     <- nly
!
! and "arrayk" has dimensions
!    NKy = 85      <- nky
!    NKx = 342     <- nkx
!
! After x transform, resulting array ("array_temp") is of size
!     kx = 513     <- NLx/2+1
!     ly = 128     <- nly
!
! which is then dealiased down to 
!    Nkx = 342     <- nkx
!     ly = 128     <- nly
!
! at the same time that the array is transposed into 
! arrayk, which has dimensions
!
!     ly = 128
!    Nkx = 342
! 
! after y transform, resulting array is of size 
!     ky = 128
!    Nkx = 342
!
! which is then dealiased down to 
!    Nky = 85 
!    Nkx = 342
!
! in the "map" routine
!
    !. size is nlx * nly_par * nlz_par * (ngtot-gmin+1)
    real,    intent(in),  dimension(:, :, :, gmin:) :: array
    !. size is nky * nkx_par * nlz_par * (ngtot-gmin+1)
    complex, intent(out), dimension(:, :, :, gmin:) :: arrayk

    complex, allocatable, dimension(:, :, :, :) :: array_temp, ak

    integer :: i, j, k, ng

#   ifdef perf
    logical, save :: first = .true.
#   endif
    Myperfon('2DFT_Fo4')

    allocate( array_temp(nkx, nly_par, nlz_par, gmin:ngtot) )
    allocate( ak        (nky, nkx_par, nlz_par, gmin:ngtot) )

    array_temp = cmplx(0.0, 0.0)
    ak = cmplx(0.0, 0.0)

    !... FFT in x-direction (forward)
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Ffo4_x')
    end if
#   endif
    do ng = gmin, ngtot
       do k = 1, nlz_par
          call Nfourndirect_X(array(:, :, k, ng), array_temp(:, :, k, ng))
       end do
    end do
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Ffo4_x
    end if
#   endif

    do ng = gmin, ngtot
       do k = 1, nlz_par
          do i = 1, nkx
             array_temp(i, :, k, ng) = array_temp(i, :, k, ng) &
                  &              * exp(-36.0*((i*1.0-1.0)/((nkx-1)*1.0))**36) !Hou-Li expression
          end do
       end do
    end do
    !
    !... data transposition xy
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Ffo4_gath')
    end if
#   endif
    call gather(r2k4, array_temp, ak)
    !
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Ffo4_gath
    end if
#   endif

!this is dealising. The y-direction delaising is done in the 
!map subroutine.
   
    !do ng = gmin, ngtot
    !   do k = 1, nlz_par
    !do i=1, nkx_par
    !   iglobal=i+mod(iproc,NPE)*nkx_par
    !   if (iglobal <= nkx)then
    !      ak(:,i, k, ng)=ak(:,i, k, ng)
    !   else
    !      ak(:,i, k, ng)=0.0
    !   end if
    !end do
    !   end do
    !end do
    !   end if

    !
    !... FFT in y-direction (forward)
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Ffo4_y')
    end if
#   endif
    do ng = gmin, ngtot
       do k = 1, nlz_par
          call Nfourndirect_Y(ak(:, :, k, ng))

!         call map(ak(:, :, k, ng), arrayk(:, :, k, ng))

          !NFL 29/05/2013: Hou-Li filter:  
          !I could instead do rho(kx)rho(ky) (where rho is the filter function).

          do j = 1, nky
             arrayk(j, :, k, ng) = ak(j, :, k, ng) * exp( -36.0*(abs(ky(j))/(nky/2.*lx/ly))**36 )
          end do
       end do
    end do
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Ffo4_y
    end if
    if (first) first = .false.
#   endif

    !do ng = gmin, ngtot
    !do k = 1, nlz_par
    !do i=1, nkx_par 
    !   do j=1,nky
    !      arrayk(j,i, k, ng)=ak(j,i, k, ng)*exp(-36.0*&
    !           (kperp(j,i)/sqrt((nky/2.*lx/ly)**2+(nlx/2.)**2))**(36.0))
    !   end do
    !end do
    !end do
    !end do

    deallocate (ak, array_temp)

    Myperfoff ! 2DFT_Fo4

  end subroutine FFT2d_direct4


!deprecated: replaced by FFT2d_inv2. This subroutine can be REMOVED!
!> 2D FFT backward (complex, 2D)
!*********************************************************  
  subroutine FFT1d_inv(arrayk, array)
    
    !... transform to real space

    use redistribute, only: scatter

    implicit none

    complex, intent(in),  dimension(:, :) :: arrayk
    real,    intent(out), dimension(:, :) :: array

    complex, allocatable, dimension(:, :) :: array_temp, ak

#   ifdef perf
    logical, save :: first = .true.
#   endif

    Myperfon('2DFT_B1')

    allocate (ak(nly, nkx_par))
!    call unmap(arrayk, ak)
    allocate (array_temp(nlx/2+1, nly_par))

    !... Initialize arrays
    ak = arrayk
    array_temp = cmplx(0.0, 0.0)

#   ifdef perf
    if (.not. first) then
       Myperfon('F1b_y')
    end if
#   endif
    call Nfourninv_Y(ak)
#   ifdef perf
    if (.not. first) then
       Myperfoff ! F1b_y
    end if
#   endif

    !. data transposition yx
#   ifdef perf
    if (.not. first) then
       Myperfon('F1b_scat')
    end if
#   endif
    call scatter (r2k3, ak, array_temp)
#   ifdef perf
    if (.not. first) then
       Myperfoff ! F1b_scat
    end if
#   endif

!    do i=nkx+1,nlx/2+1
!       array_temp(i,:)=0.
!    end do

!if (iproc==0)then
!   do j=1,nly_par 
!      do i=1, nlx/2+1 
!         print*, i, j, real(array_temp(i,j))
!      end do
!   end do
!end if

#   ifdef perf
    if (.not. first) then
       Myperfon('F1b_x')
    end if
#   endif
    call Nfourninv_X(array_temp, array)
#   ifdef perf
    if (.not. first) then
       Myperfoff ! F1b_x
    end if
    if (first) first = .false.
#   endif

    deallocate (ak, array_temp)

    Myperfoff ! 2DFT_B1

  end subroutine FFT1d_inv


!TTR
!> 2D FFT backward (complex, 2D)
!*********************************************************  
  subroutine FFT2d_inv2(arrayk, array)
    
    !... transform to real space

    use redistribute, only: scatter

    implicit none

    complex, intent(in),  dimension(:, :) :: arrayk
    real,    intent(out), dimension(:, :) :: array

    complex, allocatable, dimension(:, :) :: array_temp, ak

#   ifdef perf
    logical, save :: first = .true.
#   endif

    Myperfon('2DFT_Bo2')

    allocate (ak(nly, nkx_par))
!    call unmap(arrayk, ak)
    allocate (array_temp(nlx/2+1, nly_par))

    !... Initialize arrays
    ak(:, :) = arrayk(:, :)

#   ifdef perf
    if (.not. first) then
       Myperfon('Fbo2_y')
    end if
#   endif
    call Nfourninv_Y(ak)
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Fbo2_y
    end if
#   endif

    !. data transposition yx
#   ifdef perf
    if (.not. first) then
       Myperfon('Fbo2_scat')
    end if
#   endif
    call scatter (r2k3, ak, array_temp)
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Fbo2_scat
    end if
#   endif

#   ifdef perf
    if (.not. first) then
       Myperfon('Fbo2_x')
    end if
#   endif
    call Nfourninv_X(array_temp, array)
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Fbo2_x
    end if
    if (first) first = .false.
#   endif

    deallocate (array_temp)
    deallocate (ak)

    Myperfoff ! 2DFT_Bo2

  end subroutine FFT2d_inv2


!TTR
!> 2D FFT backward (complex, 3D)
!*********************************************************  
  subroutine FFT2d_inv3(arrayk, array)

    !... transform to real space

    use redistribute, only: scatter

    implicit none

    complex, intent(in),  dimension(:, :, :) :: arrayk
    real,    intent(out), dimension(:, :, :) :: array

    complex, allocatable, dimension(:, :, :) :: array_temp, ak

    integer :: k

#   ifdef perf
    logical, save :: first = .true.
#   endif

    !. Viriato: dim3 = nlz_par  &  dim4 = ngtot+1

    Myperfon('2DFT_Bo3')

    allocate (ak        (nly, nkx_par, nlz_par))
    allocate (array_temp(nkx, nly_par, nlz_par))

    !... Initialize arrays
    ak(:, :, :) = arrayk(:, :, :)

    !
    !... FFT in y-direction (backward)
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Fbo3_y')
    end if
#   endif
    do k = 1, nlz_par
       call Nfourninv_Y(ak(:, :, k))
    end do
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Fbo3_y
    end if
#   endif
    !
    !... data transposition yx
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Fbo3_scat')
    end if
#   endif
    !
    call scatter (r2k3, ak, array_temp)
    !
#      ifdef perf
    if (.not. first) then
       Myperfoff ! Fbo3_scat
    end if
#   endif
    !
    !... FFT in x-direction (backward)
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Fbo3_x')
    end if
#   endif
    do k = 1, nlz_par
       call Nfourninv_X(array_temp(:, :, k), array(:, :, k))
    end do
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Fbo3_x
    end if
    if (first) first = .false.
#   endif

    deallocate (ak, array_temp)

    Myperfoff ! 2DFT_Bo3

  end subroutine FFT2d_inv3


!TTR
!> 2D FFT backward (complex, 4D)
!*********************************************************  
  subroutine FFT2d_inv4(arrayk, array)

    !... transform to real space

    use redistribute, only: scatter

    implicit none

    complex, intent(in),  dimension(:, :, :, gmin:) :: arrayk
    real,    intent(out), dimension(:, :, :, gmin:) :: array

    complex, allocatable, dimension(:, :, :, :) :: array_temp, ak

    integer :: k, ng

#   ifdef perf
    logical, save :: first = .true.
#   endif

    Myperfon('2DFT_Bo4')

    allocate (ak        (nly, nkx_par, nlz_par, gmin:ngtot))
    allocate (array_temp(nkx, nly_par, nlz_par, gmin:ngtot))

    !... Initialize arrays
    ak(:, :, :, :) = arrayk(:, :, :, :)

    !
    !... FFT in y-direction (backward)
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Fbo4_y')
    end if
#   endif
    do ng = gmin, ngtot
       do k = 1, nlz_par
          call Nfourninv_Y(ak(:, :, k, ng))
       end do
    end do
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Fbo4_y
    end if
#   endif
    !
    !... data transposition yx
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Fbo4_scat')
    end if
#   endif
    !
    call scatter (r2k4, ak, array_temp)
    !
#      ifdef perf
    if (.not. first) then
       Myperfoff ! Fbo4_scat
    end if
#   endif
    !
    !... FFT in x-direction (backward)
    !
#   ifdef perf
    if (.not. first) then
       Myperfon('Fbo4_x')
    end if
#   endif
    do ng = gmin, ngtot
       do k = 1, nlz_par
          call Nfourninv_X(array_temp(:, :, k, ng), array(:, :, k, ng))
       end do
    end do
#   ifdef perf
    if (.not. first) then
       Myperfoff ! Fbo4_x
    end if
    if (first) first = .false.
#   endif

    deallocate (ak, array_temp)

    Myperfoff ! 2DFT_Bo4

  end subroutine FFT2d_inv4


!*********************************************************  
  subroutine Nfourndirect_X(array, ak)	   
    
    use fft_work
    implicit none

    real,    intent(in),  dimension(:,:) :: array
    complex, intent(out), dimension(:,:) :: ak

    logical, save :: first = .true.
    real, save :: scale
    integer :: imax, jmax

    imax=size(array, 1)
    jmax=size(array, 2)

    if (first) then
       call init_rcfftw (fft_x2k, -1, imax, iplan_x2k)
       scale = 1./sqrt(real(imax))
       first = .false.
    end if

    call rfftwnd_f77_real_to_complex (iplan_x2k, jmax, array, 1, imax, aK, 1, imax/2+1)

    ak = ak*scale

  end subroutine Nfourndirect_X


!*********************************************************  
  subroutine Nfourndirect_Y(ak)	   
    
    use fft_work

    implicit none
    
    complex, intent(inout), dimension(:,:):: ak
    
    logical, save :: first = .true.
    complex :: dummy
    real, save :: scale
    integer :: imax, jmax
    
    imax=size(ak, 1)
    jmax=size(ak, 2)

    if (first) then
       iplan_y2k = 10  ! do fft in place
       call init_ccfftw (fft_y2k, -1, imax, iplan_y2k)
       scale = 1./sqrt(real(imax))
       first = .false.
    end if

    call fftwnd_f77 (iplan_y2k, jmax, aK, 1, imax, dummy, 0, 0)  

    ak = ak*scale

  end subroutine Nfourndirect_Y
  

!*********************************************************  
  subroutine Nfourninv_Y(ak)
    
    use fft_work
    implicit none
    
    complex, intent(inout), dimension(:,:) :: ak

    logical, save :: first = .true.
    complex :: dummy
    real, save :: scale
    integer :: i, j!, k1, k2
 
    i=size(ak, 1)
    j=size(ak, 2)
    
    if (first) then
       iplan_k2y = 10 ! do fft in place
       call init_ccfftw (fft_k2y, 1, i, iplan_k2y)
       scale = 1./sqrt(real(i))
       first = .false.
    end if
    
    call fftwnd_f77 (iplan_k2y, j, ak, 1, i, dummy, 0, 0)
    
    ak = ak*scale

  end subroutine Nfourninv_Y
  

!*********************************************************  
  subroutine Nfourninv_X (akx, array)

    use fft_work
    implicit none
    
    complex, intent(in), dimension(:,:):: akx
    real, intent(inout), dimension(:,:):: array
    
    logical, save :: first = .true.
    real, save :: scale
    integer :: i, j

    i=size(array, 1)
    j=size(array, 2)
   
    if (first) then
       call init_crfftw (fft_k2x, 1, i, iplan_k2x)
       scale = 1./sqrt(real(i))
       first = .false.
    end if
    
    call rfftwnd_f77_complex_to_real (iplan_k2x, j, aKx, 1, i/2+1, array, 1, i)
    
    array=array*scale
    
  end subroutine Nfourninv_X
  

!*********************************************************  
  subroutine map (ak, arrayk)

    use mp

    implicit none

    complex, dimension (:,:) :: ak, arrayk
    integer :: i, j

    arrayK = 0.    
    do i=1,NKx_par
       do j=1,Nky/2+1
          arrayK(j,i) = ak(j,i) !*sqrt(scale)
       end do
    end do
    do i= 1, NKx_par
       do j=NLy-Nky/2+1, NLy
          arrayK(j-NLy+nky,i) = ak(j,i)!*sqrt(scale)
       end do
    end do

    
  end subroutine map
  

!*********************************************************  
  subroutine unmap(arrayk,ak)

    implicit none
    complex, dimension(:,:) :: arrayk, ak
    integer::i,j
    
    ak=0.0
    do i=1,NKx_par
       do j=1,NKy/2+1
          ak(j,i)=arrayK(j,i)
       end do
    end do
    do i= 1, NKx_par
       do j=NLy-NKy/2+1, NLy
          ak(j,i)=arrayK(j-NLy+NKy,i)
       end do
    end do
    
  end subroutine unmap
  

!*********************************************************  
  subroutine oneDfourndirect(array, arrayk,dim)	   

    !!!use constants
    use fft_work

    implicit none

    integer::dim
    real, dimension(dim):: array
    complex, dimension(dim/2+1):: arrayk

    type (fft_type), save :: fft_x2k
    integer(Kind=POINTER_KIND), save :: iplan_x2k
    logical, save :: first = .true.

    if (first) then
       call init_rcfftw (fft_x2k, -1, dim, iplan_x2k)
       first = .false.
    end if

    call rfftwnd_f77_one_real_to_complex (iplan_x2k, array, arrayK)

    arrayK = arrayk*sqrt(1./dim)

  end subroutine oneDfourndirect


!*********************************************************  
  subroutine oneDfourninv(arrayK, array,dim)

    !!!use constants
    use fft_work
    implicit none

    integer::dim
    complex, dimension(dim/2+1)::arrayk
    real, DIMENSION(dim):: array

    type (fft_type), save :: fft_k2x
    integer(Kind=POINTER_KIND), save :: iplan_k2x
    logical, save :: first = .true.

    if (first) then
       call init_crfftw (fft_k2x, 1, dim, iplan_k2x)
       first = .false.
    end if

    call rfftwnd_f77_one_complex_to_real (iplan_k2x, arrayK, array)

    array = array*sqrt(1./dim)

  end subroutine oneDfourninv


!*********************************************************  
! AVK: Z fourier transform subroutine
  subroutine four_direct_z(array, arrayk)

    use constants, only : nlz, nkz
    use fft_work

    implicit none
   
    !integer :: dim
    complex, dimension(nlz) :: temp
    complex, dimension(nlz), intent(in) :: array
    complex, dimension(nkz), intent(out) :: arrayk

    type (fft_type), save :: fft_z2k
    integer(Kind=POINTER_KIND), save :: iplan_z2k
    logical, save :: first = .true.

    if (first) then
       call init_ccfftw (fft_z2k, -1, nlz, iplan_z2k)
       first = .false.
    end if

    !call fftw_f77_one (iplan_z2k, array, temp)  
    !call fftwnd_f77_one (iplan_z2k, array, arrayk)  
    call fftwnd_f77_one (iplan_z2k, array, temp)  

    call mask_z (temp, arrayk)

    arrayk = arrayk/(1.*nlz)

  end subroutine four_direct_z


!*********************************************************  
  subroutine mask_z (array_nlz, array_nkz)

    use constants, only : nlz, nkz
    implicit none

    complex, dimension(:) :: array_nlz, array_nkz
    integer :: iz

    do iz = 1, nkz/2 + 1
       array_nkz(iz) = array_nlz(iz)
    end do
    do iz = nlz-nkz/2 + 1, nlz
       array_nkz(iz - nlz + nky) = array_nlz(iz)
    end do

  end subroutine mask_z

end module transforms
