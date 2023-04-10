module Brackets

  !TTR
  !. FPP macro to simplify RZG PERFLIB calls (performance-counter)
# include "perfmacro.h"

  implicit none

  private

  public :: Funcne_i, FuncAkpar_i, Funcg2, Funcgm, Func_lastg, bracket_3

  !TTR
  interface Funcgm
     module procedure Funcgm_3, Funcgm_4
  end interface

  !TTR
  interface Bracket
     module procedure Bracket_3, Bracket_4
  end interface


contains


!TTR
!*********************************************************  
  subroutine Bracket_3(dxF, dyF, dxG, dyG, braxyk)

    !this subroutine calculates the "bracket" of two quantities, i.e., 
    !dxk/dx*dyk/dy-dxk/dy*dyk/dx. inputs are xk and yk, the quantities in
    ! k space. returns
    !braxyk, the value of the bracket also in k-space.

    use mp, only: proc0, iproc
    use constants, only: nkx_par, nky, nlx, nly_par, nlz_par, npe, linear
    use transforms, only: FFT2d_direct

    implicit none

    !...  intent IN
    !... size is nlx * nly_par * nlz_par
    real,    dimension(:, :, :), intent(in) :: dxF, dyF, dxG, dyG

    !...  intent OUT
    complex, dimension(nky, nkx_par, nlz_par), intent(out) :: braxyk

    !... local vars
    real, dimension(nlx, nly_par, nlz_par) :: braxy
    integer :: i, j, k

    if(linear) then  !testing on this additional feature is to be completed
    braxyk = cmplx(0.0, 0.0) !sets the results of the Poisson brackets to zero if the 'linear test' option is selected
    else 
    braxy = dxF*dyG - dyF*dxG

    !TTR
#   ifdef gasca2d
    do k = 1, nlz_par
       call FFT2d_direct (braxy(:, :, k), braxyk(:, :, k))
    end do
#   elif defined(gasca3d)
    call FFT2d_direct (braxy(:, :, :), braxyk(:, :, :))
#   endif

    if (mod(iproc, npe) == 0) then
       braxyk(1, 1, :) = cmplx(0.0, 0.0) !insures that no zeroth mode is created
    end if
    
    end if

  end subroutine Bracket_3


!TTR
!*********************************************************  
  subroutine Bracket_4(dxF, dyF, dxG, dyG, braxyk)

    !this subroutine calculates the "bracket" of two quantities, i.e., 
    !dxk/dx*dyk/dy-dxk/dy*dyk/dx. inputs are xk and yk, the quantities in
    ! k space. returns
    !braxyk, the value of the bracket also in k-space.

    use mp, only: proc0, iproc
    use constants, only: nkx_par, nky, nlx, nly_par, nlz_par, npe, gmin, ngtot, linear
    use transforms, only: FFT2d_direct
    implicit none

    !...  intent IN
    !... size is nlx * nly_par * nlz_par
    real,    dimension(:, :, :), intent(in) :: dxF, dyF
    !... size is nlx * nly_par * nlz_par * (ngtot-gmin+1)
    real,    dimension(:, :, :, gmin:), intent(in) :: dxG, dyG

    !...  intent OUT
    !... size is nky * nkx_par * nlz_par * (ngtot-gmin+1)
    complex, dimension(:, :, :, gmin:), intent(out) :: braxyk

    !... local vars
    real, dimension(nlx, nly_par, nlz_par, gmin:ngtot) :: braxy
    integer :: i, j, k, ng

    if(linear) then  !tests on this feature have to be completed
    braxyk = cmplx(0.0,0.0) ! 
    else 

    do ng = gmin+1, ngtot-1
       do k = 1, nlz_par
          do j = 1, nly_par
             do i = 1, nlx
                braxy(i, j, k, ng) = dxF(i, j, k)*dyG(i, j, k, ng) - dyF(i, j, k)*dxG(i, j, k, ng)
             end do
          end do
       end do
    end do

    call FFT2d_direct (braxy(:, :, :, :), braxyk(:, :, :, :))

    if (mod(iproc, npe) == 0) then
       !braxyk(1, 1, :, gmin+1:ngtot-1) = 0.0 !insures that no zeroth mode is created
       braxyk(1, 1, :, gmin:ngtot) = cmplx(0.0, 0.0) !insures that no zeroth mode is created
    end if
    end if
  end subroutine Bracket_4


!***************************************************
  subroutine Funcne_i(Dxphi,Dyphi,Dxne,Dyne,DxApar,DyApar,Dxuepar,Dyuepar,fne,braakparuekpar)

    use constants, only: nky, nkx_par, nlz_par, nlx, nly_par
!for dnek/dt=fne
!aa=phik
!bb=nek
!cc=akpar
!dd=uekpar
!ee=nek_eq

    implicit none
    real, dimension(nlx, nly_par,nlz_par),intent(in)::dxphi,dyphi,dxne,dyne,dxapar,dyapar,dxuepar,dyuepar
    complex, dimension(nky, nkx_par,nlz_par)::braphiknek !,braakparuekpar
    complex, dimension(nky, nkx_par,nlz_par),intent(out)::fne,braakparuekpar
    integer::i,j,k

    call bracket(dxphi,dyphi,dxne,dyne,braphiknek)
    call bracket(dxapar,dyapar,dxuepar,dyuepar,braakparuekpar)
    do k=1,nlz_par
       do i=1,nkx_par
          do j=1,nky
             fne(j,i,k)=-braphiknek(j,i,k)+braakparuekpar(j,i,k) !&
             !&-(niu*kperp(j,i)**2+niu2*kperp(j,i)**4)*bb(j,i)
          end do
       end do
    end do
  end subroutine Funcne_i


!**************************************************
  subroutine FuncAkpar_i(DxApar,DyApar,Dxphi,Dyphi,Dxne,Dyne,Dxuepar,Dyuepar,Dxg2,Dyg2,akpar,FApar,t)

    use constants  !,only:rhos,nkx_par,nky,nlx,nly_par,de,nlz_par
    use grid, only: kperp, zz, ky
    use mp, only: iproc

    implicit none
    real, dimension(nlx,nly_par,nlz_par),intent(in)::dxphi,dyphi,dxne,dyne,dxapar,dyapar
    real, dimension(nlx,nly_par,nlz_par),intent(in)::dxuepar,dyuepar,dxg2,dyg2
    real, dimension(nlx,nly_par,nlz_par)::tempx,tempy
    complex, dimension(nky,nkx_par,nlz_par)::braakparphik,brauekparphik,akpar
    complex, dimension(nky,nkx_par,nlz_par),intent(out)::fapar
    integer::i,j,k
    real::t

    tempx =dxphi -rhos**2*(dxne +sqrt(2.0)*dxg2 )
    tempy =dyphi -rhos**2*(dyne +sqrt(2.0)*dyg2 )

    call bracket(dxapar,dyapar,tempx,tempy,braakparphik)
    call bracket(dxuepar,dyuepar,dxphi,dyphi,brauekparphik)

    do k=1,nlz_par
       do i=1,nkx_par
          do j=1,nky
             fapar(j,i,k)=1.0/(1.0+kperp(j,i)**2*de**2)*&
                  (braakparphik(j,i,k)-de**2*brauekparphik(j,i,k)&
                  -notanj*1.0/sqrt(2.0)*rhos_de*rhoe_LTe*(0.0,1.0)*ky(j)*Akpar(j,i,k))
          end do
       end do
    end do

!NFL, 05/08/2012:  this is the antenna, but I never checked that it is
!working properly!!!
!   if(mod(iproc,npe)==0) then
!      do k=1,Nlz_par
!         fapar(j1,2,k)= fapar(j1,2,k) + amplitude*exp(-(0.0,1.0)*omega0*t)*cos(kpar0*2*pi*zz(k)/Lz)
!!         fapar(j2,1,k)= fapar(j2,1,k) + amplitude*exp((0.0,1.0)*omega0*t)*cos(kpar0*2*pi*zz(k)/Lz)
!      end do
!   end if
  end subroutine FuncAkpar_i


!*********************************************************      
  subroutine Funcg2(Dxg2, Dyg2, Dxphi, Dyphi, Dxapar, Dyapar, Dxg3, Dyg3, braakparuekpar, phik, fg2)

    use constants
    use grid, only: ky

    implicit none

    !... size is nlx * nly_par * nlz_par
    real, dimension(:, :, :), intent(in) ::  Dxapar, Dyapar, Dxphi, Dyphi
    real, dimension(:, :, :), intent(in) ::  Dxg2, Dyg2, Dxg3, Dyg3

    !... size is nky * nkx_par * nlz_par
    complex, dimension(:, :, :), intent(out) :: fg2

    complex, dimension(nky, nkx_par, nlz_par) :: brag2phik, braakparuekpar, braakparg3, phik
    integer :: i, j, k

    call bracket(Dxg2,   Dyg2,   Dxphi, Dyphi, brag2phik)
    call bracket(Dxapar, Dyapar, Dxg3,  Dyg3,  braakparg3)

    do k = 1, nlz_par
       do i = 1, nkx_par
          do j = 1, nky
             fg2(j,i,k) = brag2phik(j,i,k) + sqrt(gmin+1.0) * rhos_de * braakparg3(j,i,k) &
                  + notanj * sqrt(2.0) * braakparuekpar(j,i,k) &
                  - notanj * 1.0/(2.0*rhos*de) * rhoe_LTe*(0.0,1.0) * ky(j) * Phik(j,i,k)
!                  -notanj*0.5*rhos_de*rhoe_LTe*(0.0,1.0)*ky(j)*Phik(j,i,k)
!NFL: 24/05/13 commented line above had wrong normalizations
          end do
       end do
    end do

  end subroutine Funcg2


!TTR
!*********************************************************      
  subroutine Funcgm_3(m, Dxgm, Dygm, Dxg, Dyg, Dxgp, Dygp, Dxphi, Dyphi, Dxapar, Dyapar, akpar, fgm, t)

    use constants, only: gmin, ngtot, npe, nlx, nly, nly_par, nlz_par, nky, nkx_par, &
         &               kpar0, j1, omega0, facpm, amplitude, Lz, pi, lambda, de, &
         &               notanj, rhos_de, rhoe_lte
    use grid, only: ky, zz
    use mp, only: iproc, proc0
    use Functions, only: anj_kron
  
    implicit none

    !...  intent IN
    integer, intent(in) :: m
    !...  size is nlx * nly_par * nlz_par
    real, dimension(:, :, :), intent(in) :: DxApar, DyApar, Dxphi, Dyphi
    real, dimension(:, :, :), intent(in) :: Dxgm, Dygm, Dxg, Dyg, Dxgp, Dygp

    !...  intent OUT
    !... size is nky * nkx_par * nlz_par
    complex, dimension(:, :, :), intent(out) :: fgm

    !... Local vars
    complex, dimension(nky, nkx_par, nlz_par) :: braakpargpm, braphikg, Akpar
    real,    dimension(gmin+1:ngtot) :: lte_kron

    real :: t
    integer :: i, j, k

    lte_kron = 0.0
    lte_kron(gmin+1) = 1.0

    call Bracket(Dxphi,  Dyphi, Dxg, Dyg, braphikg)

    call Bracket(DxApar, DyApar, &
         sqrt((m+1)*1.0)*Dxgp + sqrt(m*1.0)*(1.-1./lambda*anj_kron(m)*(1.-notanj))*Dxgm, &
         sqrt((m+1)*1.0)*Dygp + sqrt(m*1.0)*(1.-1./lambda*anj_kron(m)*(1.-notanj))*Dygm, &
  	   braakpargpm)

    do k = 1, nlz_par
       do i = 1, nkx_par
          do j = 1, nky
             fgm(j, i, k) = &
                  -braphikg(j, i, k) + rhos_de * braakpargpm(j, i, k) &
                  + notanj * lte_kron(m) * sqrt(3.0)/(2.0*de**2) &
                  * rhoe_LTe*(0.0, 1.0) * ky(j) * Akpar(j, i, k)
!                  +notanj*lte_kron(m)*sqrt(3.0/2.0)*rhos_de*rhoe_LTe*(0.0,1.0)*ky(j)*Akpar(j,i,k)
!NFL: 24/05/13 commented line above had wrong normalizations
          end do
       end do
    end do

    !AVK: 20/03/2013: dbpar antenna
    !AVK: Needs (b dot grad) operating on the antenna
    if(mod(iproc, npe) == 0) then
       do k = 1, nlz_par
          fgm(j1, 2, k) = fgm(j1, 2, k) - anj_kron(m) * sqrt(nlx*nly*1.) &
               &          * sqrt(omega0) * facpm*amplitude &
               &          * exp(-(0.0, 1.0)*omega0*t) * cos(kpar0*2.*pi*zz(k)/Lz)
       end do
    end if

  end subroutine Funcgm_3


!TTR
!*********************************************************
  subroutine Funcgm_4(Dxg1, Dyg1, Dxg2, Dyg2, Dxphi, Dyphi, Dxapar, Dyapar, Akpar, fgm, t)

    use constants, only: gmin, ngtot, npe, nlx, nly, nly_par, nlz_par, nky, nkx_par, &
         &               kpar0, j1, omega0, facpm, amplitude, Lz, pi, lambda, de, &
         &               notanj, rhos_de, rhoe_lte
    use grid, only: ky, zz
    use mp, only: iproc, proc0
    use Functions, only: anj_kron
  
    implicit none

    !...  intent IN
    !...  size is nlx * nly_par * nlz_par * (ngtot-gmin+1)
    real, dimension(:, :, :, gmin:), intent(in) :: Dxg1, Dyg1, Dxg2, Dyg2
    !...  size is nlx * nly_par * nlz_par
    real, dimension(:, :, :), intent(in) :: DxApar, DyApar, Dxphi, Dyphi

    !...  intent OUT
    !... size is nky * nkx_par * nlz_par * (ngtot-gmin+1)
    complex, dimension(:, :, :, gmin:), intent(out) :: fgm

    !... Local vars
    complex, dimension(nky, nkx_par, nlz_par, gmin:ngtot) :: braakpargpm, braphikg
    real,    dimension(nlx, nly_par, nlz_par, gmin:ngtot) :: fDxg, fDyg
    complex, dimension(nky, nkx_par, nlz_par) :: Akpar
    real,    dimension(gmin+1:ngtot) :: lte_kron

    real :: t
    integer :: i, j, m, k

    fDxg(:, :, :, gmin) = 0.0
    fDxg(:, :, :, ngtot) = 0.0
    fDyg(:, :, :, gmin) = 0.0
    fDyg(:, :, :, ngtot) = 0.0
    do m = gmin+1, ngtot-1
       fDxg(:, :, :, m) = sqrt((m+1)*1.0)*Dxg2(:, :, :, m+1) &
            &           + sqrt(m*1.0)*(1.-1./lambda*anj_kron(m)*(1.-notanj)) * Dxg1(:, :, :, m-1)
       fDyg(:, :, :, m) = sqrt((m+1)*1.0)*Dyg2(:, :, :, m+1) &
            &           + sqrt(m*1.0)*(1.-1./lambda*anj_kron(m)*(1.-notanj)) * Dyg1(:, :, :, m-1)
    end do

    lte_kron = 0.0
    lte_kron(gmin+1) = 1.0

    call Bracket(Dxphi(:, :, :),   Dyphi(:, :, :), &
         &       Dxg2(:, :, :, :), Dyg2(:, :, :, :), braphikg(:, :, :, :))

    call Bracket(DxApar(:, :, :),   DyApar(:, :, :), &
         &       fDxg(:, :, :, :),  fDyg(:, :, :, :), braakpargpm(:, :, :, :))

    do m = gmin+1, ngtot-1
       do k = 1, nlz_par
          do i = 1, nkx_par
             do j = 1, nky
                fgm(j, i, k, m) = &
                     -braphikg(j, i, k, m) + rhos_de * braakpargpm(j, i, k, m) &
                     + notanj * lte_kron(m) * sqrt(3.0)/(2.0*de**2) &
                     * rhoe_LTe*(0.0, 1.0) * ky(j) * Akpar(j, i, k)
!                     +notanj*lte_kron(m)*sqrt(3.0/2.0)*rhos_de*rhoe_LTe*(0.0,1.0)*ky(j)*Akpar(j,i,k)
!NFL: 24/05/13 commented line above had wrong normalizations
             end do
          end do
       end do
    end do

    !AVK: 20/03/2013: dbpar antenna
    !AVK: Needs (b dot grad) operating on the antenna
    if(mod(iproc, npe) == 0) then
       do m = gmin+1, ngtot-1
          do k = 1, nlz_par
             fgm(j1, 2, k, m) = fgm(j1, 2, k, m) - anj_kron(m) * sqrt(nlx*nly*1.) &
                  &             * sqrt(omega0) * facpm*amplitude &
                  &             * exp(-(0.0, 1.0)*omega0*t) * cos(kpar0*2.*pi*zz(k)/Lz)
          end do
       end do
    end if

  end subroutine Funcgm_4


!****************************************************
  subroutine Func_lastg(hyper_nuei,niu2,Dxgm,Dygm,Dxg,Dyg,Dxphi,Dyphi,Dxapar,Dyapar,fglast)
    !Calculates g_M using the analytical nonlinearclosure
    !As is, only does the perp part

    use constants
    use mp
    use diag, only: Convol, Convol2
    use grid

    implicit none

    complex, dimension(nky, nkx_par,nlz_par), intent(out) :: fglast

    real,    dimension(nlx, nly_par, nlz_par) :: dxapar, dyapar, dxphi, dyphi
    real,    dimension(nlx, nly_par, nlz_par) :: dxgm, dygm, dxg, dyg, dxbra, dybra
    complex, dimension(nky, nkx_par, nlz_par) :: braakparg, brafikg, totbra

    real :: hyper_nuei, niu2
    integer :: i, j, k
    
    call Bracket(Dxphi,  Dyphi,  Dxg, Dyg, brafikg)
    call Bracket(Dxapar, Dyapar, Dxg, Dyg, braakparg)

    do k = 1, nlz_par
       do i = 1, nkx_par
          do j = 1, nky
             braakparg(j,i,k) = rhos_de**2 * (ngtot+1)/((ngtot+1) * nu_ei &
                  + (ngtot+1)**(2*hyper_morder) * hyper_nuei &
                  + niu * kperp(j,i)**2 + niu2*kperp(j,i)**(2*hyper_order)) * braakparg(j,i,k)
          end do
       end do
    end do
    !TTR
#   ifdef gasca2d
    call Convol(braakparg, dxbra, dybra)
#   elif defined(gasca3d)
    call Convol2(braakparg, dxbra, dybra)
#   endif

    call Bracket(Dxapar, Dyapar, Dxbra + rhos_de * sqrt(ngtot*1.0) * Dxgm, &
         Dybra + rhos_de * sqrt(ngtot*1.0) * Dygm, totbra)

    fglast = -brafikg + totbra

  end subroutine Func_lastg

end module Brackets
