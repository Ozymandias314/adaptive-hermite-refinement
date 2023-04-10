module grid
!
! has information about the real and k-space grids, including
! how this information is shared among processors.
!
  interface idx_local
     module procedure idx_local_r, idx_local_k
  end interface

  interface proc_id
     module procedure proc_id_r, proc_id_k
  end interface

  type :: r_layout_type
     integer :: iproc
     integer :: NLx, NLy, NLy_par
     integer :: llim_world, ulim_world, llim_proc, ulim_proc
  end type r_layout_type

  type :: k_layout_type
     integer :: iproc
     integer :: NKx, NKy, NKx_par
     integer :: llim_world, ulim_world, llim_proc, ulim_proc
  end type k_layout_type

  type (r_layout_type) :: r_variable
  type (k_layout_type) :: k_variable

contains


!*********************************************************  
  subroutine init_grid
    
    use mp, only: iproc
    use constants
    implicit none
    !logical, save :: initialized = .false.

    r_variable%iproc = iproc
    r_variable%nlx = nlx
    r_variable%nly = nly
    r_variable%nly_par = nly_par
    r_variable%llim_world = 1 
    r_variable%ulim_world = NLy
    r_variable%llim_proc = 1   
    r_variable%ulim_proc = NLy_par

    k_variable%iproc = iproc
    k_variable%nkx = nkx
    k_variable%nky = nky
    k_variable%nkx_par = nkx_par
    k_variable%llim_world = 1 
    k_variable%ulim_world = NKx
    k_variable%llim_proc = 1   
    k_variable%ulim_proc = NKx_par

  end subroutine init_grid


!*********************************************************  
  real function kperp (j, i)
    implicit none
    integer :: j,i
    kperp=sqrt(ky(j)**2+kx(i)**2)
  end function kperp


!TTR (not used: meant for the case with the XOR transpose not implemented in Viriato)
!*********************************************************  
  real function kperp_new (j, i)
    implicit none
    integer :: j,i
    kperp_new=sqrt(ky(j)**2+kx_new(i)**2)
  end function kperp_new


!*********************************************************  
  real function xx(i)
    ! returns x value corresponding to i index
    ! in the real space layout.
    use constants
    implicit none
    integer::i
    xx=lx*(i-1.-nlx/2)/nlx

  end function xx


!*********************************************************  
  real function yy(j)
    ! returns y value corresponding to a local j index
    ! in the real space layout.

    use constants
    implicit none
    integer :: j

    yy=ly*(jglobal(j)-1.-nly/2)/nly

  end function yy


!*********************************************************  
  real function zz(i)
    ! returns z value corresponding to i index
    ! in the real space layout.
    use constants
    implicit none
    integer::i
    zz=lz*(kglobal(i)-1.-nlz/2)/nlz

  end function zz


!*********************************************************  
  integer function jglobal(j)
    use constants, only: nly_par,NPE
    use mp, only: iproc
    implicit none
    integer::j
    
!    jglobal=j+iproc*nly_par
    jglobal=j+mod(iproc,NPE)*nly_par

    end function jglobal


!*********************************************************  
  integer function kglobal(k)
    use constants, only: nlz_par, npe
    use mp, only: iproc
    implicit none
    integer::k

    kglobal = k+iproc/NPE*nlz_par

  end function kglobal


!*********************************************************  
  real function y_glob(j)
    use constants
    implicit none
    integer::j
    y_glob=ly/nly*(j-1-nly/2)  
    
  end function y_glob


!*********************************************************  
  real function kx(i)
    !
    ! returns k number corresponding to a local i index
    ! in the k-space layout.
    !
    ! Presumably the real kx = 2*pi*k/Lx
    !
    use constants, only: npe, nkx_par
    use mp, only: iproc
    implicit none
    integer::i, iglobal

    iglobal = mod(iproc, npe)*nkx_par + i

    !   if(i <= nlx/2+1)then   
    kx = (iglobal-1.)
    !   else
    !      kx=(i-nlx-1.)
    !   endif

  end function kx


!TTR (not used: meant for the case with the XOR transpose not implemented in Viriato)
!*********************************************************  
  real function kx_new(i)
    !
    ! returns k number corresponding to a local i index
    ! in the k-space layout.
    !
    ! Presumably the real kx = 2*pi*k/Lx
    !
    use constants, only: npe, nkx_par
    use mp, only: iproc
    implicit none
    integer::i, iglobal

    iglobal = mod(iproc, npe)*(nkx_par-1) + i

    kx_new = (iglobal-1.)

  end function kx_new


!*********************************************************  
  logical function kkeep (i)

    use constants, only: Nkx_par, Nkx, NPE
    use mp, only: iproc
    implicit none
    integer, intent (in) :: i
    integer :: iglobal
    
    iglobal = mod(iproc,NPE)*Nkx_par + i

    kkeep = iglobal <= Nkx

  end function kkeep


!*********************************************************  
  real function ky(j)
    !
    ! returns k number corresponding to a j index
    ! in the k-space layout.
    !
    ! Presumably the real ky = 2*pi*k/Ly
    !
    use constants
    implicit none
    integer::j
    if(j<=nky/2+1)then
       ky=(j-1.)*lx/ly
    else
       ky=(j-nky-1.)*lx/ly
    endif
  end function ky
  

!*********************************************************  
  !AVK: kz for postproc
  real function kz(j)
    !
    ! returns k number corresponding to a j index
    ! in the k-space layout.
    !
    ! Presumably the real kz = 2*pi*k/Lz
    !
    use constants
    implicit none
    integer::j
    if(j<=nkz/2+1)then
       kz=(j-1.)*lx/lz
    else
       kz=(j-nkz-1.)*lx/lz
    endif

  end function kz


!*********************************************************  
  real function gama0(x)  
    implicit none
!taken from numerical recipes; calculates the modified bessel funtion
   ! of order 0
    real:: x,ax
    real:: p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
    save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
    data p1,p2,p3,p4,p5,p6,p7/1.0d0, 3.5156229d0, 3.0899424d0,&
         & 1.2067492d0, 0.2659732d0, 0.360768d-1, 0.45813d-2/
    data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0, 0.1328592d-1,&
         & 0.225319d-2, -0.157565d-2, 0.916281d-2, -0.2057706d-1,&
         & 0.2635537d-1,-0.1647633d-1, 0.392377d-2/
    if (abs(x).lt.3.75) then
       y=(x/3.75)**2
       gama0=exp(-x)*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
    else
       ax=abs(x)
       y=3.75/ax
       gama0=(exp(ax-x)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y &
            *(q7+y*(q8+y*q9))))))))
    endif
    
  end function gama0


!*********************************************************  
  function proc_id_r (r_variable, j, zk)

! returns the number of the processor that has the 
! global indices j, zk (in the (y,z) plane, that is) for a real-space variable
    use constants, only: nlz_par, NPE
    implicit none
    integer :: proc_id_r
    type (r_layout_type), intent (in) :: r_variable
    integer, intent (in) :: j, zk

    proc_id_r = (j-1)/r_variable%NLy_par + (zk-1)/nlz_par*NPE

  end function proc_id_r
  

!*********************************************************  
  function proc_id_k (k_variable, i, zk)
! returns the number of the processor that has the 
! global index i for a k-space variable
    use constants, only: nlz_par, NPE
    implicit none
    integer :: proc_id_k
    type (k_layout_type), intent (in) :: k_variable
    integer, intent (in) :: i, zk

    proc_id_k = (i-1)/k_variable%Nkx_par + (zk-1)/nlz_par*NPE

  end function proc_id_k


!*********************************************************  
! function proc_id_z(k)
!   use constants, only: nlz_par
!   implicit none
!   integer :: proc_id_z
!   integer, intent (in):: k
!
!   proc_id_z = (k-1)/nlz_par
!
! end function proc_id_z


!*********************************************************  
  function idx_local_r (r, j, zk)

    implicit none
    logical :: idx_local_r
    type (r_layout_type), intent (in) :: r
    integer, intent (in) :: j, zk

    idx_local_r = r%iproc == proc_id(r, j, zk)

  end function idx_local_r


!*********************************************************  
  function idx_local_k (k, i, zk)

    implicit none
    logical :: idx_local_k
    type (k_layout_type), intent (in) :: k
    integer, intent (in) :: i, zk

    idx_local_k = k%iproc == proc_id(k, i, zk)

  end function idx_local_k

end module grid
