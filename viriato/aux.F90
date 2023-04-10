module Aux

  implicit none

contains

!*********************************************************  
  subroutine dtnext(RELATIVE_ERROR,x,noinc,dti)
    use constants, only: epsilon
    implicit none
    real, intent(in) :: RELATIVE_ERROR,x
    real, intent(in out)::dti
    real::inc_fac!,low,high
    logical,intent (in out)::noinc
 
    if (noinc) then
       inc_fac=1.0
       noinc=.false.
    else
       inc_fac=1.08
    end if


    if (RELATIVE_ERROR .LT. 0.8*epsilon) then
       if (x .LT. inc_fac*dti) then
          dti=x
       else
          dti=inc_fac*dti
       end if
    else
       dti=min(x,dti)
    end if
    
  end subroutine dtnext


!****************************************************
  subroutine increase_fac(relative_error,inc_fac)
    use constants, only: epsilon
    implicit none
    real, intent(in)  :: relative_error
    real, intent(out) :: inc_fac
    if (RELATIVE_ERROR .LT. 0.01*epsilon) inc_fac=4.0
    if ((RELATIVE_ERROR .GE. 0.01*epsilon) .AND. (RELATIVE_ERROR .LT. 0.1*epsilon)) inc_fac=2.0
    if ((RELATIVE_ERROR .GE. 0.1*epsilon) .AND. (RELATIVE_ERROR .LT. 0.5*epsilon)) inc_fac=1.5
    if ((RELATIVE_ERROR .GE. 0.5*epsilon) .AND. (RELATIVE_ERROR .LT. 0.8*epsilon)) inc_fac=1.2

  end subroutine increase_fac


!****************************************************
  subroutine hyper_diff (bperp_max,omega_kaw,niu2,res2)
    use constants,only: nky, nkx, hyper_coef, rhos,rhoi, small_rhoi, de
    use grid, only:gama0,ky
    implicit none

    real,intent(in)::bperp_max
    real,intent(out):: omega_kaw,niu2,res2
    real::omega_kaw_nl
    real,save::kperp_dum
    logical,save::first=.true.

    if(first) then
       kperp_dum=sqrt(ky(nky/2+1)**2+(nkx*1.)**2)
       first=.false.
    end if

    if (rhoi.LT.small_rhoi) then
       omega_kaw=sqrt(1.+kperp_dum**2*(3./4.*rhoi**2+rhos**2))*&
            & ky(nky/2+1)*bperp_max/(1.0+kperp_dum**2*de**2)
       omega_kaw_nl=sqrt(1.+kperp_dum**2*(3./4.*rhoi**2+rhos**2))*&
            & kperp_dum*bperp_max/(1.0+kperp_dum**2*de**2)
    else
       omega_kaw=kperp_dum*&
            & sqrt(rhos**2-rhoi**2/(gama0(0.5*kperp_dum**2*rhoi**2)-1.))*&
            & ky(nky/2+1)*bperp_max/(1.0+kperp_dum**2*de**2)
       omega_kaw_nl=kperp_dum**2*&
            & sqrt(rhos**2-rhoi**2/(gama0(0.5*kperp_dum**2*rhoi**2)-1.))*&
            & bperp_max/(1.0+kperp_dum**2*de**2)
    end if

!    niu2=hyper_coef*omega_kaw_nl/kperp_dum**4
    niu2=hyper_coef*omega_kaw_nl/kperp_dum**4
    res2=niu2

  end subroutine hyper_diff


!*********************************************************  
  subroutine omegakaw (bperp_max,omega_kaw)

    use constants, only: nky, nkx, rhos, rhoi, small_rhoi, de, dz, three_D, pi
    use grid, only: gama0, ky
    implicit none

    real, intent(in)  :: bperp_max
    real, intent(out) :: omega_kaw
    real :: omega_kaw_nl
    real, save :: kperp_dum
    logical, save :: first=.true.

    if(first) then
       kperp_dum=sqrt(ky(nky/2+1)**2+(nkx*1.)**2)
       first=.false.
    end if

    if (rhoi.LT.small_rhoi) then
       omega_kaw=sqrt(1.+kperp_dum**2*(3./4.*rhoi**2+rhos**2))*&
            & ky(nky/2+1)*bperp_max/(1.0+kperp_dum**2*de**2)
       omega_kaw_nl=sqrt(1.+kperp_dum**2*(3./4.*rhoi**2+rhos**2))*&
            & kperp_dum*bperp_max/(1.0+kperp_dum**2*de**2)
    else
       if (three_D) then
          ! omega_kaw=max(kperp_dum*&
          !      & sqrt(rhos**2-rhoi**2/(gama0(0.5*kperp_dum**2*rhoi**2)-1.))*&
          !      & ky(nky/2+1)*bperp_max/sqrt(1.0+kperp_dum**2*de**2),kperp_dum*&
          !      & sqrt(rhos**2-0.5*rhoi**2/(gama0(0.5*kperp_dum**2*rhoi**2)-1.))&
          !      & *2.0*Pi/(dz*sqrt(1.0+kperp_dum**2*de**2)))
          omega_kaw=kperp_dum*&
               & sqrt(rhos**2-0.5*rhoi**2/(gama0(0.5*kperp_dum**2*rhoi**2)-1.))&
               & *2.0*Pi/(dz*sqrt(1.0+kperp_dum**2*de**2))
       else
          omega_kaw=kperp_dum*&
               & sqrt(rhos**2-rhoi**2/(gama0(0.5*kperp_dum**2*rhoi**2)-1.))*&
               & ky(nky/2+1)*bperp_max/sqrt(1.0+kperp_dum**2*de**2)
       endif
       omega_kaw_nl=kperp_dum**2*&
            & sqrt(rhos**2-rhoi**2/(gama0(0.5*kperp_dum**2*rhoi**2)-1.))*&
            & bperp_max/(1.0+kperp_dum**2*de**2)
    end if

  end subroutine omegakaw


!*********************************************************  
  subroutine z_diffusion(coeff,dt,field,newfield)
    !NFL, 25/11/2014
    !laplacian diffusion in z: can either be just for numerical purposes or
    ! to represent the d^2/dz^2 term in the asymptotic closure.
    use constants
    use mp
    implicit none

    complex, dimension(nky, nkx_par,nlz_par):: field, newfield
    complex,dimension(nky,nkx_par,1):: arrayk
    real::coeff,dt
    integer::k,ip

    do k=2,nlz_par-1
       newfield(:,:,k) = field(:,:,k) + coeff*dt/dz**2*&
            (field(:,:,k+1)+field(:,:,k-1)-2.*field(:,:,k))
    end do
    
    !pass k=nlz_par array forward:
    
    arrayk(:,:,1) = field(:,:,nlz_par)
    
    do ip=0, npe*npez-1-npe
       if (iproc==ip .and. iproc /= ip+npe) then
          call zsend(arrayk,ip+npe)
       end if
    end do
    do ip=npe*npez-npe, npe*npez-1
       if (iproc==ip .and. iproc /= mod(ip,npe)) then
          call zsend(arrayk,mod(ip,npe))
       end if
    end do
    do ip=npe,npe*npez-1
       if (iproc==ip .and. iproc /= ip-npe) then
          call zreceive(arrayk,ip-npe)
       end if
    end do
    do ip=0, npe-1
       if (iproc==ip .and. iproc /= npe*npez-npe+ip) then
          call zreceive(arrayk,npe*npez-npe+ip)
       end if
    end do
    call barrier
   
    newfield(:,:,1) = newfield(:,:,1) + coeff*dt/dz**2*&
         (field(:,:,2)+arrayk(:,:,1)-2.*field(:,:,1))
   
    !pass k=1 array backward:

    arrayk(:,:,1) = field(:,:,1)
    do ip=npe, npe*npez-1
       if (iproc==ip .and. iproc /= ip-npe) then
          call zsend(arrayk,ip-npe)
       end if
    end do
    do ip=0, npe-1
       if (iproc==ip .and. iproc /= npe*npez-npe+ip ) then
          call zsend(arrayk,npe*npez-npe+ip)
       end if
    end do

    do ip=0,npe*npez-NPE-1
       if (iproc==ip .and. iproc /= ip+npe) then
          call zreceive(arrayk,ip+NPE)
       end if
    end do
    do ip=npe*npez-NPE, npe*npez-1
       if (iproc==ip .and. iproc /= mod(ip,npe)) then
          call zreceive(arrayk,mod(ip,NPE))
       end if
    end do
    call barrier

    newfield(:,:,nlz_par) = newfield(:,:,nlz_par) + coeff*dt/dz**2*&
         (arrayk(:,:,1)+field(:,:,nlz_par-1)-2.*field(:,:,nlz_par))

  end subroutine z_diffusion


!**********************************************************
  subroutine PHI_POT(nek,phiK)
  
    use mp, only: iproc
    use grid, only: kperp, gama0
    use constants, only: rhoi, nky, nkx_par, nlz_par, small_rhoi, npe
    implicit none
  
    complex,dimension(nky,nkx_par,nlz_par),intent(in)::nek
    complex,dimension(nky,nkx_par,nlz_par),intent(out)::phiK
    integer::i,j

    if (rhoi.LT. small_rhoi)then
       if (mod(iproc,NPE)==0) then
          phiK(1,1,:)=0.0
          do j=2,nky 
             phiK(j,1,:)=-neK(j,1,:)/kperp(j,1)**2
          enddo
      
          do i=2,nkx_par
             do j=1,nky
                phiK(j,i,:)=-neK(j,i,:)/kperp(j,i)**2
             enddo
          enddo
       else
          do i=1,nkx_par
             do j=1,nky
                phiK(j,i,:)=-neK(j,i,:)/kperp(j,i)**2
             enddo
          enddo
       endif
    else
       if (mod(iproc,NPE)==0) then
          phiK(1,1,:)=0.0
          do j=2,nky 
             phiK(j,1,:) = rhoi**2*0.5D0/(gama0(kperp(j,1)**2*rhoi**2*0.5D0)-1.)*neK(j,1,:)
!t.          phiK(j,1,:) = rhoi**2*0.5D0/(gama0(kperp(j,1)**2*rhoi**2*0.5)-1.)*neK(j,1,:)
          enddo

          do i=2,nkx_par
             do j=1,nky
                phiK(j,i,:) = rhoi**2*0.5D0/(gama0(kperp(j,i)**2*rhoi**2*0.5D0)-1.)*neK(j,i,:) 
!t.             phiK(j,i,:) = rhoi**2*0.5D0/(gama0(kperp(j,i)**2*rhoi**2*0.5)-1.)*neK(j,i,:) 
             enddo
          enddo
       else
          do i=1,nkx_par
             do j=1,nky
                phiK(j,i,:) = rhoi**2*0.5D0/(gama0(kperp(j,i)**2*rhoi**2*0.5D0)-1.)*neK(j,i,:) 
!t.             phiK(j,i,:) = rhoi**2*0.5D0/(gama0(kperp(j,i)**2*rhoi**2*0.5)-1.)*neK(j,i,:) 
             enddo
          enddo
       endif
    endif
  end subroutine PHI_POT


!*********************************************************      
  subroutine SEMI_IMP_OP(dti,bperp_max,aa0,SI_oper)
    use mp, only: iproc
    use grid, only: kperp, gama0
    use constants,only:nky,nkx_par,nlz_par,rhos,rhoi,small_rhoi,de,npe
    implicit none

    !real,dimension(nky,nkx_par)::aa0
    complex,dimension(nky,nkx_par,nlz_par),intent(out)::SI_oper !,SI_oper_ne
    real,intent(in)::bperp_max,dti,aa0
    integer::i,j

    !aa0=theta*(nky*1./2.+1.)/sqrt((nky*1./2.+1.)**2*1.+nkx**2*1.)

    !SI_oper(j,i)=aa0**2*bperp_max**2*(rhos**2+rhoi**2)*kperp(j,i)**4*dti**2

    !bperp_max=1.0

    !aa0=theta
    if (rhoi.LT. small_rhoi)then
       if (mod(iproc,NPE)==0) then
          SI_oper(1,1,:)=0.0
          do j=2,nky 
             SI_oper(j,1,:) = aa0**2*(1.+kperp(j,1)**2*(3./4.*rhoi**2+rhos**2))* &
                  & kperp(j,1)**2*bperp_max**2* &
                  & dti**2/(1.0+kperp(j,1)**2*de**2) 
          enddo

          do i=2,nkx_par
             do j=1,nky
                SI_oper(j,i,:) = aa0**2*(1.+kperp(j,i)**2*(3./4.*rhoi**2+rhos**2))* &
                     & kperp(j,i)**2*bperp_max**2* &
                     & dti**2/(1.0+kperp(j,i)**2*de**2)  
             enddo
          enddo
       else
          do i=1,nkx_par
             do j=1,nky
                SI_oper(j,i,:) = aa0**2*(1.+kperp(j,i)**2*(3./4.*rhoi**2+rhos**2))* &
                     & kperp(j,i)**2*bperp_max**2* &
                     & dti**2/(1.0+kperp(j,i)**2*de**2)  
             enddo
          enddo
       endif
    else
       if (mod(iproc,npe)==0) then
          SI_oper(1,1,:)=0.0
          do j=2,nky 
             SI_oper(j,1,:) = aa0**2*(3.0*rhos**2-rhoi**2/(gama0(0.5*kperp(j,1)**2*rhoi**2)-1.))* &
                  & kperp(j,1)**4*bperp_max**2* &
                  & dti**2/(1.0+kperp(j,1)**2*de**2)  
          enddo

          do i=2,nkx_par
             do j=1,nky
                SI_oper(j,i,:) = aa0**2*(3.0*rhos**2-rhoi**2/(gama0(0.5*kperp(j,i)**2*rhoi**2)-1.))* &
                     & kperp(j,i)**4*bperp_max**2* &
                     & dti**2/(1.0+kperp(j,i)**2*de**2)  
             enddo
          enddo
       else
          do i=1,nkx_par
             do j=1,nky
                SI_oper(j,i,:) = aa0**2*(3.0*rhos**2-rhoi**2/(gama0(0.5*kperp(j,i)**2*rhoi**2)-1.))* &
                     & kperp(j,i)**4*bperp_max**2* &
                     & dti**2/(1.0+kperp(j,i)**2*de**2)  
             enddo
          enddo
       endif
    endif
  end subroutine SEMI_IMP_OP


 !***************************************
  subroutine resolution_check(delta,savetime)
    use constants
    implicit none
   
    real::delta,savetime
   
    ! resol_x=real(Lx)/real(Nkx)

    if (delta .LT. 2*dx) then
       write(221,*) 'not enough resolution at t=', savetime
    end if
   
  end subroutine resolution_check
 

!*********************************************************  
  subroutine double_resolution(old_array,new_array)
    !to double the resolution in restart files
    !Nuno Loureiro, 27/05/04
    
    use constants
    use mp, only: proc0, iproc, send, receive!, barrier
    implicit none
    
    real,dimension(nlx,nlypar_old,nlz_par)::old_array
    real,dimension(nlx,nly_par,nlz_par)::new_array
    real,dimension(nlx,nlz_par)::first_element,last_element
    integer::i,j,k,n
 
    ! case x (newdim_y=oldim_y) !double x-resolution only
    
    !do j=1, newdim_y
!   do i=1, newdim_x-1, i+2
    !      new_array(i,j)=old_array((i+1)/2,j)
    !   end do
    !   do i=2, newdim_x-2, i+2
    !      new_array(i,j)=(old_array(i/2,j)+old_array(i/2+1),j)/2.
    !   end do
    !   new_array(newdim_x,j)=(old_array(i/2,j)+old_array(1,j))/2.
    !end do
    
    !case y (newdim_x=oldim_x)  !double y-resolution only
   
    first_element(:,:)=old_array(:,1,:)
    do k=1,nlz_par
       if (proc0) call send(first_element(:,k),NPE-1)
       if (iproc==NPE-1) call receive(last_element(:,k),0)
  
       do n=1, NPE-1
          if (iproc==n) call send(first_element(:,k),n-1)
       end do
       do n=0, NPE-2
          if (iproc==n) call receive(last_element(:,k),n+1)
       end do
       
       do i=1, nlx
          do j=1, nly_par-1, 2
             new_array(i,j,k)=old_array(i,(j+1)/2,k)
          end do
          do j=2, nly_par-2, 2
             new_array(i,j,k)=(old_array(i,j/2,k)+old_array(i,j/2+1,k))/2.
          end do
          new_array(i,nly_par,k)=(old_array(i,nlypar_old,k)+last_element(i,k))/2.
       end do
    end do
    
  end subroutine double_resolution


!*********************************************************  
  subroutine comp_eigen(A,WR,WI,VL,VR)
    use mp
    use constants, only: dim_vec
    implicit none
    
!    integer,parameter::N=dim_vec
!    integer,parameter::LDA=N
!    integer,parameter::LDVL=N
!    integer,parameter::LDVR=N
    integer::N,LDA,LDVL,LDVR
    integer,parameter::LWMAX=1000
    integer::INFO, LWORK
    !
    !     .. Local Arrays ..
!    real(8) A( LDA, N ), VL( LDVL, N ), VR( LDVR, N ),&
!         WR( N ), WI( N ), WORK( LWMAX )
    real(8) A( dim_vec, dim_vec ), VL( dim_vec, dim_vec ), VR( dim_vec,dim_vec ),&
         WR( dim_vec ), WI( dim_vec ), WORK( LWMAX )
    external DGEEV

    N=dim_vec
    LDA=N
    LDVL=N
    LDVR=N

    !
    !     Query the optimal workspace.
    !
    LWORK = -1
    call DGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL,&
       VR, LDVR, WORK, LWORK, INFO )
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    !
    !     Solve eigenproblem.
    !
    call DGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL,&
         VR, LDVR, WORK, LWORK, INFO )
    
    !
    !     Check for convergence.
    !
!    if( INFO.gt.0 ) then
!       write(*,*)'The algorithm failed to compute eigenvalues.'
!       stop
!    end if
    
  end subroutine comp_eigen


!*********************************************************  
  subroutine inv(A,Ainv)
    use constants, only:dim_vec
    implicit none
    real(8), dimension(dim_vec,dim_vec), intent(in) :: A
    real(8), dimension(size(A,1),size(A,2)) :: Ainv
    
    real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI
    
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)
    
!    if (info /= 0) then
!       stop 'Matrix is numerically singular!'
!    end if
    
    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)
    
!    if (info /= 0) then
!       stop 'Matrix inversion failed!'
!    end if
    !end function inv
  end subroutine inv


!*********************************************************  
  subroutine calc_leftover(dxapar, dyapar, dxgMp1, dygMp1, gkM, sumleftover)

    use constants
    use transforms, only: FFT2d_inv
    use mp

    implicit none

    integer::i,j,k
    real, dimension(nlx, nly_par, nlz_par) :: dxapar, dyapar, dxgMp1, dygMp1, gM, leftover
    complex, dimension(nky, nkx_par, nlz_par) :: gkM
    real::sumleftover

    !note that ngtot = M+1
    sumleftover=0.0
    !TTR
#   ifdef gasca2d
    do k = 1, nlz_par
       call FFT2d_inv(gkM(:, :, k), gM(:, :, k))
    end do
#   elif defined(gasca3d)
    call FFT2d_inv(gkM(:, :, :), gM(:, :, :))
#   endif
    leftover = rhos**3/de*sqrt(ngtot*1.0) * gM * (dxapar*dygMp1 - dyapar*dxgMp1)
    
    do k = 1, nlz_par
       do j = 1, nly_par
          do i = 1, nlx
             sumleftover = sumleftover + leftover(i, j, k)
          end do
       end do
    end do
    sumleftover = sumleftover / (nlx * nly * nlz * 1.0)
    call sum_allreduce(sumleftover)

  end subroutine calc_leftover


!*********************************************************  
  subroutine closure(gkMm1,gkM,M,hyper_nuei,gkMp1)
    use constants
    use grid
    implicit none
    complex, dimension(nky, nkx_par)::gkM,gkMp1,gkMm1
    integer::i,j,M
    real::hyper_nuei

    !M closure:
    
    !      gk_star(:,:,ngtot)=  4.0*gk_star(:,:,ngtot-1) - 6.0*gk_star(:,:,ngtot-2) &
    !           + 4.0*gk_star(:,:,ngtot-3) - gk_star(:,:,ngtot-4)
    !      gk_star(:,:,ngtot) = 2.0*gk_star(:,:,ngtot-1) - gk_star(:,:,ngtot-2)
!    gkMp1= gkM
    gkMp1=0.0
    do i=1,nkx_par
       do j=2,nky
          gkMp1(j,i) = -((0.0,1.0)*(M*nu_ei+M**(2*hyper_morder)*hyper_nuei)*gkM(j,i) + &
               rhos_de*sqrt(M*1.0)*ky(j)*gkMm1(j,i))/rhos_de*1.0/sqrt(M+1.)*1.0/ky(j)
       end do
    end do
!    gk_star(:,:,ngtot) = 0.0

  end subroutine closure

end module Aux
