module mp
!
! Easier Fortran90 interface to the MPI Message Passing Library.
!
!     (c) Copyright 1991 to 1998 by Michael A. Beer, William D. Dorland, 
!     P. B. Snyder, Q. P. Liu, and Gregory W. Hammett. ALL RIGHTS RESERVED.
!     
! Note: mp_mpi_r8.f90 is a version of mp_mpi.f90 to use when compiling 
! with -r8 (where the default real type is taken to be 8 bytes).  Just 
! replaced all occurances of MPI_REAL with MPI_DOUBLE_PRECISION and 
! MPI_COMPLEX with MPI_DOUBLE_COMPLEX.
!
! Added mpi_mod.F90 to allow for compile-time argument checks by replacing 
! "include 'mpif.h'" with "use mpi". For backward compatibility, the 
! former is still kept, and can be invoked by using "FPPFLAGS_helios = -Dmpi1" 
! in the Makefile

  implicit none
  private

  public :: init_mp, finish_mp
  public :: broadcast, sum_reduce, sum_allreduce
  public :: max_reduce, max_allreduce
  public :: min_reduce, min_allreduce
  public :: nproc, iproc, proc0
  public :: send, receive, zsend, zreceive, isend, ireceive
  public :: barrier

  integer :: nproc, iproc
  logical :: proc0

  interface broadcast
     module procedure broadcast_integer 
     module procedure broadcast_integer_array 

     module procedure broadcast_real    
     module procedure broadcast_real_array    
     module procedure broadcast_real_2d_array

     module procedure broadcast_complex 
     module procedure broadcast_complex_array

     module procedure broadcast_logical 
     module procedure broadcast_logical_array 

     module procedure bcastfrom_integer 
     module procedure bcastfrom_integer_array 

     module procedure bcastfrom_real    
     module procedure bcastfrom_real_array    

     module procedure bcastfrom_complex 
     module procedure bcastfrom_complex_array 

     module procedure bcastfrom_logical 
     module procedure bcastfrom_logical_array 

     module procedure broadcast_character
     module procedure bcastfrom_character
  end interface

  interface sum_reduce
     module procedure sum_reduce_integer
     module procedure sum_reduce_integer_array

     module procedure sum_reduce_real
     module procedure sum_reduce_real_array

     module procedure sum_reduce_complex
     module procedure sum_reduce_complex_array
  end interface

  interface sum_allreduce
     module procedure sum_allreduce_integer
     module procedure sum_allreduce_integer_array

     module procedure sum_allreduce_real
     module procedure sum_allreduce_real_array

     module procedure sum_allreduce_complex
     module procedure sum_allreduce_complex_array
  end interface

  interface max_reduce
     module procedure max_reduce_integer
     module procedure max_reduce_integer_array

     module procedure max_reduce_real
     module procedure max_reduce_real_array
  end interface

  interface max_allreduce
     module procedure max_allreduce_integer
     module procedure max_allreduce_integer_array

     module procedure max_allreduce_real
     module procedure max_allreduce_real_array
  end interface

  interface min_reduce
     module procedure min_reduce_integer
     module procedure min_reduce_integer_array

     module procedure min_reduce_real
     module procedure min_reduce_real_array
  end interface

  interface min_allreduce
     module procedure min_allreduce_integer
     module procedure min_allreduce_integer_array

     module procedure min_allreduce_real
     module procedure min_allreduce_real_array
  end interface

  interface send
     module procedure send_integer
     module procedure send_integer_array

     module procedure send_real
     module procedure send_real_array

     module procedure send_complex
     module procedure send_complex_array

     module procedure send_logical
     module procedure send_logical_array
  end interface

  interface receive
     module procedure receive_integer
     module procedure receive_integer_array

     module procedure receive_real
     module procedure receive_real_array

     module procedure receive_complex
     module procedure receive_complex_array

     module procedure receive_logical
     module procedure receive_logical_array
  end interface

contains

!*********************************************************
  subroutine init_mp
    use mpi_mod
    implicit none
    integer :: ierror, rank

    call mpi_init (ierror)
    call mpi_comm_size (mpi_comm_world, nproc, ierror)
    call mpi_comm_rank (mpi_comm_world, iproc, ierror)
    proc0 = iproc == 0
  end subroutine init_mp

!*********************************************************
  subroutine finish_mp
    use mpi_mod
    implicit none
    integer :: ierror

    call mpi_finalize (ierror)
  end subroutine finish_mp

! ************** broadcasts *****************************

!*********************************************************
  subroutine broadcast_character (char)
    use mpi_mod
    implicit none
    character(*), intent (in out) :: char
    integer :: ierror
    call mpi_bcast (char, len(char), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_character

!*********************************************************
  subroutine broadcast_integer (i)
    use mpi_mod
    implicit none
    integer, intent (in out) :: i
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_integer

!*********************************************************
  subroutine broadcast_integer_array (i)
    use mpi_mod
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_integer_array

!*********************************************************
  subroutine broadcast_real (x)
    use mpi_mod
    implicit none
    real, intent (in out) :: x
    integer :: ierror
    call mpi_bcast (x, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_real

!*********************************************************
  subroutine broadcast_real_array (x)
    use mpi_mod
    implicit none
    real, dimension (:), intent (in out) :: x
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_real_array

!*********************************************************
  subroutine broadcast_real_2d_array (x)
    use mpi_mod
    implicit none
    real,dimension(:,:), intent (in out) :: x
    integer :: ierror
    call mpi_bcast(x, size(x), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_real_2d_array

!*********************************************************
  subroutine broadcast_complex (z)
    use mpi_mod
    implicit none
    complex, intent (in out) :: z
    integer :: ierror
    call mpi_bcast (z, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_complex

!*********************************************************
  subroutine broadcast_complex_array (z)
    use mpi_mod
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_complex_array

!*********************************************************
  subroutine broadcast_logical (f)
    use mpi_mod
    implicit none
    logical, intent (in out) :: f
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_logical

!*********************************************************
  subroutine broadcast_logical_array (f)
    use mpi_mod
    implicit none
    logical, dimension (:), intent (in out) :: f
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_logical_array

!*********************************************************
  subroutine bcastfrom_logical (f, src)
    use mpi_mod
    implicit none
    logical, intent (in out) :: f
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_logical

!*********************************************************
  subroutine bcastfrom_logical_array (f, src)
    use mpi_mod
    implicit none
    logical, dimension (:), intent (in out) :: f
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_logical_array

!*********************************************************
  subroutine bcastfrom_character (c, src)
    use mpi_mod
    implicit none
    character(*), intent (in out) :: c
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (c, len(c), MPI_CHARACTER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_character

!*********************************************************
  subroutine bcastfrom_integer (i, src)
    use mpi_mod
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_integer

!*********************************************************
  subroutine bcastfrom_integer_array (i, src)
    use mpi_mod
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_integer_array

!*********************************************************
  subroutine bcastfrom_real (x, src)
    use mpi_mod
    implicit none
    real, intent (in out) :: x
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (x, 1, MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_real

!*********************************************************
  subroutine bcastfrom_real_array (x, src)
    use mpi_mod
    implicit none
    real, dimension (:), intent (in out) :: x
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_real_array

!*********************************************************
  subroutine bcastfrom_complex (z, src)
    use mpi_mod
    implicit none
    complex, intent (in out) :: z
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (z, 1, MPI_DOUBLE_COMPLEX, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_complex

!*********************************************************
  subroutine bcastfrom_complex_array (z, src)
    use mpi_mod
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_complex_array

! ************** reductions ***********************

!*********************************************************
  subroutine sum_reduce_integer (i, dest)
    use mpi_mod
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_integer

!*********************************************************
  subroutine sum_reduce_integer_array (i, dest)
    use mpi_mod
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_integer_array

!*********************************************************
  subroutine sum_reduce_real (a, dest)
    use mpi_mod
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_real

!*********************************************************
  subroutine sum_reduce_real_array (a, dest)
    use mpi_mod
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_real_array

!*********************************************************
  subroutine sum_reduce_complex (z, dest)
    use mpi_mod
    implicit none
    complex, intent (in out) :: z
    integer, intent (in) :: dest
    complex :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_complex

!*********************************************************
  subroutine sum_reduce_complex_array (z, dest)
    use mpi_mod
    implicit none
    complex, dimension (:), intent (in out) :: z
    integer, intent (in) :: dest
    complex, dimension (size(z)) :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, size(z), MPI_DOUBLE_COMPLEX, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_complex_array

!*********************************************************
  subroutine sum_allreduce_integer (i)
    use mpi_mod
    implicit none
    integer, intent (in out) :: i
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_integer

!*********************************************************
  subroutine sum_allreduce_integer_array (i)
    use mpi_mod
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_integer_array

!*********************************************************
  subroutine sum_allreduce_real (a)
    use mpi_mod
    implicit none
    real, intent (in out) :: a
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_real

!*********************************************************
  subroutine sum_allreduce_real_array (a)
    use mpi_mod
    implicit none
    real, dimension (:), intent (in out) :: a
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_real_array

!*********************************************************
  subroutine sum_allreduce_complex (z)
    use mpi_mod
    implicit none
    complex, intent (in out) :: z
    complex :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_complex

!*********************************************************
  subroutine sum_allreduce_complex_array (z)
    use mpi_mod
    implicit none
    complex, dimension (:), intent (in out) :: z
    complex, dimension (size(z)) :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, size(z), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_complex_array

!*********************************************************
  subroutine max_reduce_integer (i, dest)
    use mpi_mod
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_integer

!*********************************************************
  subroutine max_reduce_integer_array (i, dest)
    use mpi_mod
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_integer_array

!*********************************************************
  subroutine max_reduce_real (a, dest)
    use mpi_mod
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_real

!*********************************************************
  subroutine max_reduce_real_array (a, dest)
    use mpi_mod
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_real_array

!*********************************************************
  subroutine max_allreduce_integer (i)
    use mpi_mod
    implicit none
    integer, intent (in out) :: i
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_integer

!*********************************************************
  subroutine max_allreduce_integer_array (i)
    use mpi_mod
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_integer_array

!*********************************************************
  subroutine max_allreduce_real (a)
    use mpi_mod
    implicit none
    real, intent (in out) :: a
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_real

!*********************************************************
  subroutine max_allreduce_real_array (a)
    use mpi_mod
    implicit none
    real, dimension (:), intent (in out) :: a
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_real_array

!*********************************************************
  subroutine min_reduce_integer (i, dest)
    use mpi_mod
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_integer

!*********************************************************
  subroutine min_reduce_integer_array (i, dest)
    use mpi_mod
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_integer_array

!*********************************************************
  subroutine min_reduce_real (a, dest)
    use mpi_mod
    implicit none
    real, intent (in out) :: a
    integer, intent (in) :: dest
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_real

!*********************************************************
  subroutine min_reduce_real_array (a, dest)
    use mpi_mod
    implicit none
    real, dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_real_array

!*********************************************************
  subroutine min_allreduce_integer (i)
    use mpi_mod
    implicit none
    integer, intent (in out) :: i
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_integer

!*********************************************************
  subroutine min_allreduce_integer_array (i)
    use mpi_mod
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_integer_array

!*********************************************************
  subroutine min_allreduce_real (a)
    use mpi_mod
    implicit none
    real, intent (in out) :: a
    real :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_real

!*********************************************************
  subroutine min_allreduce_real_array (a)
    use mpi_mod
    implicit none
    real, dimension (:), intent (in out) :: a
    real, dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_real_array

! ********************* barrier **********************

!*********************************************************
  subroutine barrier
    use mpi_mod
    implicit none
    integer :: ierror
    call mpi_barrier (MPI_COMM_WORLD, ierror)
  end subroutine barrier

! ********************* sends **********************

!*********************************************************
  subroutine send_integer (i, dest, tag)
    use mpi_mod
    implicit none
    integer, intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, 1, MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_integer

!*********************************************************
  subroutine send_integer_array (i, dest, tag)
    use mpi_mod
    implicit none
    integer, dimension (:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, size(i), MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_integer_array

!*********************************************************
  subroutine send_real (a, dest, tag)
    use mpi_mod
    implicit none
    real, intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, 1, MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_real

!*********************************************************
  subroutine send_real_array (a, dest, tag)
    use mpi_mod
    implicit none
    real, dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, size(a), MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_real_array

!*********************************************************
  subroutine send_complex (z, dest, tag)
    use mpi_mod
    implicit none
    complex, intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, 1, MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_complex

!*********************************************************
  subroutine isend (z, dest, request, tag)
    use mpi_mod
    implicit none
    complex, intent (in) :: z
    integer, intent (in) :: dest
    INTEGER, INTENT (OUT) :: request
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    !ttr call mpi_isend (z, 1, MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
    call mpi_isend (z, 1, MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, request, ierror)
  end subroutine isend

!*********************************************************
  subroutine send_complex_array (z, dest, tag)
    use mpi_mod
    implicit none
    complex, dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, size(z), MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_complex_array

!*********************************************************
  subroutine zsend (z, dest, tag)
    use mpi_mod
    implicit none
    complex, dimension (:,:,:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer :: dim1, dim2, dim3, k, i
    tagp = 0
    if (present(tag)) tagp = tag
    
    dim1 = size(z,1)
    dim2 = size(z,2)
    dim3 = size(z,3)
    
    ! AVK: instead of sending the whole array, we send it in chunks to avoid buffer
    !overflow
    do k = 1,dim3
       do i = 1,dim2
          call mpi_send (z(:,i,k), dim1, MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
       end do
    end do
  end subroutine zsend

!*********************************************************
  subroutine send_logical (f, dest, tag)
    use mpi_mod
    implicit none
    logical, intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, 1, MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_logical

!*********************************************************
  subroutine send_logical_array (f, dest, tag)
    use mpi_mod
    implicit none
    logical, dimension (:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, size(f), MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_logical_array

!*********************************************************
  subroutine send_character (s, dest, tag)
    use mpi_mod
    implicit none
    character(*), intent (in) :: s
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send &
         (s, len(s), MPI_CHARACTER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_character

! ********************* receives  **********************

!*********************************************************
  subroutine receive_integer (i, src, tag)
    use mpi_mod
    implicit none
    integer, intent (out) :: i
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, 1, MPI_INTEGER, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_integer

!*********************************************************
  subroutine receive_integer_array (i, src, tag)
    use mpi_mod
    implicit none
    integer, dimension (:), intent (out) :: i
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, size(i), MPI_INTEGER, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_integer_array

!*********************************************************
  subroutine receive_real (a, src, tag)
    use mpi_mod
    implicit none
    real, intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, 1, MPI_DOUBLE_PRECISION, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_real

!*********************************************************
  subroutine receive_real_array (a, src, tag)
    use mpi_mod
    implicit none
    real, dimension (:), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, size(a), MPI_DOUBLE_PRECISION, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_real_array

!*********************************************************
  subroutine receive_complex (z, src, tag)
    use mpi_mod
    implicit none
    complex, intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, 1, MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_complex

!*********************************************************
  subroutine ireceive(z, src, request, tag)
    use mpi_mod
    implicit none
    complex, intent (out) :: z
    integer, intent (in) :: src
    INTEGER, INTENT (OUT) :: request
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (z, 1, MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        request, ierror)
  end subroutine ireceive

!*********************************************************
  subroutine receive_complex_array (z, src, tag)
    use mpi_mod
    implicit none
    complex, dimension (:), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, size(z), MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_complex_array

!*********************************************************
  subroutine zreceive (z, src, tag)
    use mpi_mod
    implicit none
    complex, dimension (:,:,:), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
	integer :: dim1, dim2, dim3, k, l
    tagp = 0
    if (present(tag)) tagp = tag


	dim1 = size(z,1)
	dim2 = size(z,2)
	dim3 = size(z,3)

	! AVK: instead of sending the whole array, we send it in chunks to avoid buffer
	!overflow
	do k =1,dim3
		do l = 1,dim2
			call mpi_recv (z(:,l,k), dim1, MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
	end do 
	end do
  end subroutine zreceive

!*********************************************************
  subroutine receive_logical (f, src, tag)
    use mpi_mod
    implicit none
    logical, intent (out) :: f
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, 1, MPI_LOGICAL, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_logical

!*********************************************************
  subroutine receive_logical_array (f, src, tag)
    use mpi_mod
    implicit none
    logical, dimension (:), intent (out) :: f
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, size(f), MPI_LOGICAL, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_logical_array

end module mp
