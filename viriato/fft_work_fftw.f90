module fft_work

  implicit none

!  integer, private, parameter :: one_star_eight=1
  integer*8, private, parameter :: one_star_eight=1
  integer, public, parameter :: POINTER_KIND=Kind(one_star_eight)

  public :: fft_type
  public :: init_ccfftw, init_crfftw, init_rcfftw, init_z
  
  interface init_crfftw
     module procedure init_crfftw_1d
     module procedure init_crfftw_2d
  end interface

  interface init_rcfftw
     module procedure init_rcfftw_1d
     module procedure init_rcfftw_2d
  end interface

  interface init_ccfftw
     module procedure init_ccfftw_1d
     module procedure init_ccfftw_2d
  end interface

  type :: fft_type
     integer :: n, is, type
  end type fft_type

  integer, parameter :: fftw_estimate   =  0
  integer, parameter :: fftw_measure    =  1
  integer, parameter :: fftw_in_place   =  8
  integer, parameter :: fftw_use_wisdom = 16

  private

contains


!*********************************************************  
  subroutine init_z (fft, is, n, iplan)

    type (fft_type), intent (in out) :: fft
    integer, intent (in) :: is, n
    integer(Kind=POINTER_KIND) :: iplan
    integer :: j
    
    fft%n = n
    fft%is = is
    fft%type = 1
    
    j = fftw_measure + fftw_use_wisdom
    call fftwnd_f77_create_plan(iplan,1,n,is,j)

  end subroutine init_z


!*********************************************************  
  subroutine init_ccfftw_1d (fft, is, n, iplan)

    type (fft_type), intent (in out) :: fft
    integer, intent (in) :: is, n
    integer(Kind=POINTER_KIND) :: iplan
    integer :: j, one
    
    fft%n = n
    fft%is = is
    fft%type = 1
    
    one = 1
    if(iplan.eq.int(10,Kind=POINTER_KIND)) then
       j = fftw_in_place + fftw_measure + fftw_use_wisdom
    else
       j = fftw_measure + fftw_use_wisdom
    end if
       
    call fftwnd_f77_create_plan(iplan,one,n,is,j)

  end subroutine init_ccfftw_1d


!*********************************************************  
  subroutine init_ccfftw_2d (fft, is, nx, ny, iplan)

    type (fft_type), intent (in out) :: fft
    integer, intent (in) :: is, nx, ny
    integer(Kind=POINTER_KIND) :: iplan
    integer :: j!!!, two
    
    fft%n = nx
    fft%is = is
    fft%type = 1
    
    !!!two = 2
    if(iplan.eq.int(10,Kind=POINTER_KIND)) then
       j = fftw_in_place + fftw_measure + fftw_use_wisdom
    else
       j = fftw_measure + fftw_use_wisdom
    end if
       
    call fftw2d_f77_create_plan(iplan, nx, ny, is, j)

  end subroutine init_ccfftw_2d


!*********************************************************  
  subroutine init_rcfftw_1d (fft, is, n, iplan)

    type (fft_type), intent (in out) :: fft
    integer, intent (in) :: is, n
    integer(Kind=POINTER_KIND) :: iplan
    integer :: j

    fft%n = n
    fft%is = is

    fft%type = 0

    j = fftw_measure + fftw_use_wisdom
    call rfftwnd_f77_create_plan(iplan,1,N,is,j)

  end subroutine init_rcfftw_1d
  

!*********************************************************  
  subroutine init_crfftw_1d (fft, is, n, iplan)

    type (fft_type), intent (in out) :: fft
    integer, intent (in) :: is, n
    integer(Kind=POINTER_KIND) :: iplan
    integer :: j

    fft%n = n
    fft%is = is

    fft%type = 0

    j = fftw_measure + fftw_use_wisdom
    call rfftwnd_f77_create_plan(iplan,1,N,is,j)

  end subroutine init_crfftw_1d


!*********************************************************  
  subroutine init_rcfftw_2d (fft, is, m, n, iplan)

    type (fft_type), intent (in out) :: fft
    integer, intent (in) :: is, m, n
    integer(Kind=POINTER_KIND) :: iplan
    integer :: j

    j = fftw_measure + fftw_use_wisdom
    call rfftw2d_f77_create_plan(iplan,m,n,is,j)

  end subroutine init_rcfftw_2d
  

!*********************************************************  
  subroutine init_crfftw_2d (fft, is, m, n, iplan)

    type (fft_type), intent (in out) :: fft
    integer, intent (in) :: is, m, n
    integer(Kind=POINTER_KIND) :: iplan
    integer :: j

    j = fftw_measure + fftw_use_wisdom
    call rfftw2d_f77_create_plan(iplan,m,n,is,j)

  end subroutine init_crfftw_2d

end module fft_work
