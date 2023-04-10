module Functions

  implicit none

contains

!**********
  real function exp_nu(j, i, niu2, dti)

    use constants
    use grid
    implicit none

    integer, intent(in) :: i, j
    real,    intent(in) :: niu2, dti
  
    exp_nu = exp( -( niu*kperp(j,i)**2 + niu2*kperp(j,i)**(2*hyper_order) ) * dti )

  end function exp_nu

!**********
  real function exp_eta(j,i,res2,dti)

    use constants
    use grid

    implicit none

    integer :: i, j
    real    :: res2, dti

    exp_eta=exp(-(res*kperp(j,i)**2+res2*kperp(j,i)**(2*hyper_order))*dti/(1.0+kperp(j,i)**2*de**2))

  end function exp_eta

!**********
  real function exp_ng(ng,hyper_nuei,dti)

    use constants

    implicit none

    integer::ng
    real::hyper_nuei,dti

    exp_ng = exp(-(ng*nu_ei+ng**(2*hyper_morder)*hyper_nuei)*dti)

  end function exp_ng

!**********
  real function nu_func1(j, i, niu2, dti)

    use constants
    use grid

    implicit none

    integer, intent(in) :: i,j
    real,    intent(in) :: niu2, dti

    real    :: dd
  
    dd = niu*kperp(j,i)**2 + niu2*kperp(j,i)**(2*hyper_order)

    if (dd*dti < 0.01) then
       nu_func1 = dti
    else
       nu_func1 = 1./dd*( 1. - exp(-dd*dti) )
    end if

  end function nu_func1

!**********
  real function nu_func2(j, i, niu2, ng, hyper_nuei, dti)

    use constants
    use grid

    implicit none

    integer, intent(in) :: i, j, ng
    real,    intent(in) :: niu2, hyper_nuei, dti

    real :: dd, dc
  
    dd = niu*kperp(j,i)**2 + niu2*kperp(j,i)**(2*hyper_order)
    dc = ng*nu_ei + ng**(2*hyper_morder)*hyper_nuei

    nu_func2 = 1./(dd+dc)*( 1.-exp(-(dd+dc)*dti) )

  end function nu_func2

!**********
  real function nu_func3(j, i, niu2, dti)

    use constants
    use grid

    implicit none

    integer, intent(in) :: i, j
    real,    intent(in) :: niu2, dti

    real :: dd

    dd = niu*kperp(j,i)**2 + niu2*kperp(j,i)**(2*hyper_order)

    if ( dd*dti<0.01 ) then
       nu_func3 = 0.5*dti
    else
       nu_func3 = ( 1.-exp(-dd*dti)*(1+dd*dti) )/(dd**2*dti)
    end if

  end function nu_func3

!**********
  real function nu_func4(j, i, niu2, dti)

    use constants
    use grid

    implicit none

    integer, intent(in) :: i, j
    real,    intent(in) :: niu2, dti

    real :: dd

    dd = niu*kperp(j,i)**2 + niu2*kperp(j,i)**(2*hyper_order)

    if ( dd*dti<0.01 ) then
       nu_func4 = 0.5*dti
    else
       nu_func4 = ( dd*dti-1.+exp(-dd*dti) )/(dd**2*dti)
    end if

  end function nu_func4

!**********
  real function nu_func5(j, i, niu2, ng, hyper_nuei, dti)

    use constants
    use grid

    implicit none

    integer, intent(in) :: i, j, ng
    real,    intent(in) :: niu2, hyper_nuei, dti

    real :: dd, dc

    dd = niu*kperp(j,i)**2 + niu2*kperp(j,i)**(2*hyper_order)
    dc = ng*nu_ei + ng**(2*hyper_morder)*hyper_nuei

    dd = dd+dc
    nu_func5 = ( 1.-exp(-dd*dti)*(1+dd*dti) )/(dd**2*dti)

  end function nu_func5

!**********
  real function nu_func6(j, i, niu2, ng, hyper_nuei, dti)

    use constants
    use grid

    implicit none

    integer, intent(in) :: i, j, ng
    real,    intent(in) :: niu2, hyper_nuei, dti
  
    real :: dd, dc

    dd = niu*kperp(j,i)**2 + niu2*kperp(j,i)**(2*hyper_order)
    dc = ng*nu_ei + ng**(2*hyper_morder)*hyper_nuei
  
    nu_func6 = ( (dd+dc)*dti-1.+exp(-(dd+dc)*dti) )/( (dd+dc)**2*dti )
  
  end function nu_func6

!**********
  real function anj_kron(m)
  
    implicit none

    integer, intent(in) :: m
  
    if (m==1) then 
       anj_kron = 1.0
    else
       anj_kron = 0.0
    end if

  end function anj_kron

end module Functions
