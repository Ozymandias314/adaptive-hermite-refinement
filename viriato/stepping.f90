module Stepping

  implicit none

contains


!**************** Z-STEP *****************************************
  subroutine z_step(phik,nek,uekpar,akpar,gk,dti)

    use constants
    !    use diag
    use mp
    use grid, only: kperp
    use Functions, only: anj_kron
    use Aux,    only: PHI_POT

    implicit none
    
    integer :: i, j, k, ip, ng
    real    :: dti, new_etaz, new_etaz_g
    complex, dimension(nky, nkx_par, nlz_par) :: phik, nek, uekpar, akpar
    complex, dimension(nky, nkx_par, nlz_par) :: phik_star, nek_star, uekpar_star, akpar_star
    complex, dimension(nky, nkx_par, nlz_par) :: phiknew, neknew, uekparnew, akparnew
    complex, dimension(nky, nkx_par, nlz_par, gmin:ngtot) :: gk!,gk_star,gknew
    complex, allocatable, dimension(:,:,:,:) :: gk_star,gknew !AVK gknew might be unnecessary
    complex, allocatable, dimension(:,:,:)   :: array_zparbc
    !    complex, dimension(:,:,:,gmin:) :: gk
    !complex, dimension(:,:,:,:) :: gk
    
    if(g_inc) then
       allocate (gk_star(nky,nkx_par,nlz_par,gmin:ngtot))
       allocate (gknew(nky,nkx_par,nlz_par,gmin:ngtot))
       allocate (array_zparbc(nky,nkx_par,1:ngtot-gmin+5))
    else
       allocate (gk_star(1,1,1,1))
       allocate (gknew(1,1,1,1))
       allocate (array_zparbc(nky,nkx_par,4))
    end if
    
    !if(proc0) print*, etaz,etaz_g,dti
    
    new_etaz=etaz/dti
    new_etaz_g=etaz_g/dti
    
    !    if(proc0) print*, new_etaz,new_etaz_g
    
    !Using MacCormack scheme.
    !predictor step:
    !AVK: 04/04/13 adding z diffusion term
    
    do k=2,Nlz_par-1
       do i=1,nkx_par
          do j=1,nky
             nek_star(j,i,k) = nek(j,i,k) - dti/dz*(uekpar(j,i,k+1)-uekpar(j,i,k))&
                  +new_etaz*dti*(nek(j,i,k+1)-2.*nek(j,i,k)+nek(j,i,k-1))/(dz**2)
             
             if (g_inc) then
                akpar_star(j,i,k) = akpar(j,i,k) - dti/(dz*(1.+kperp(j,i)**2*de**2))*&
                     (phik(j,i,k+1)-phik(j,i,k)-rhos**2*(nek(j,i,k+1)-nek(j,i,k)) &
                     -rhos**2*sqrt(2.)*(gk(j,i,k+1,gmin)-gk(j,i,k,gmin))) &
                     +new_etaz*dti*(akpar(j,i,k+1) - 2.*akpar(j,i,k) + akpar(j,i,k-1))&
                     /(dz**2*(1.+kperp(j,i)**2*de**2))
                
                gk_star(j,i,k,gmin) = gk(j,i,k,gmin) &
                     -dti/dz*sqrt(gmin+1.)*rhos_de*(gk(j,i,k+1,gmin+1)-gk(j,i,k,gmin+1)) &
                     -dti/dz*sqrt(2.)*notanj*(uekpar(j,i,k+1)-uekpar(j,i,k)) &
                     +new_etaz_g*dti*(gk(j,i,k+1,gmin) - 2.*gk(j,i,k,gmin) + gk(j,i,k-1,gmin))&
                     /(dz**2)
                
                do ng=gmin+1, ngtot-1
                   gk_star(j,i,k,ng) = gk(j,i,k,ng) &
                        -dti/dz*(rhos_de)*sqrt(ng+1.0)*(gk(j,i,k+1,ng+1)-gk(j,i,k,ng+1)) &
                        -dti/dz*(rhos_de)*sqrt(ng*1.)*(1.-1./lambda*anj_kron(ng)*(1.-notanj)) &
                        *(gk(j,i,k+1,ng-1)-gk(j,i,k,ng-1)) &
                        + new_etaz_g*dti*(gk(j,i,k+1,ng) - 2.*gk(j,i,k,ng) + gk(j,i,k-1,ng))&
                        /(dz**2)
                end do
             else
                akpar_star(j,i,k) = akpar(j,i,k) - dti/(dz*(1.+kperp(j,i)**2*de**2))*&
                     (phik(j,i,k+1)-phik(j,i,k) - rhos**2*(nek(j,i,k+1)-nek(j,i,k)))&
                     +new_etaz*dti*(akpar(j,i,k+1) - 2.*akpar(j,i,k) + akpar(j,i,k-1))&
                     /(dz**2*(1.+kperp(j,i)**2*de**2))
             end if
          end do
       end do
    end do
    
    !impose periodic BCs, i.e., stuff(nlz+1) = stuff(1)
    !   if(proc0) then
    !      phizero_star(nlz) = phizero(nlz) - dti/dz*(akpar(0,0,1)-akpar(0,0,nlz_par))
    !   end if
    
    array_zparbc(:,:,1)=nek(:,:,1)
    array_zparbc(:,:,2)=phik(:,:,1)
    array_zparbc(:,:,3)=akpar(:,:,1)
    array_zparbc(:,:,4)=uekpar(:,:,1)
    if(g_inc) then
       array_zparbc(:,:,5:ngtot-gmin+5)=gk(:,:,1,gmin:ngtot)
    end if
    
    do ip=npe, npe*npez-1
       if (iproc==ip .and. iproc /= ip-npe) then
          !         print*, iproc, array_zparbc(1,1,1), nek(1,1,1)
          !         print*, '              '
          !write(*,*) "Before first zsend", ip, ip-npe
          call zsend(array_zparbc,ip-npe)
       end if
    end do
    do ip=0, npe-1
       if (iproc==ip .and. iproc /= npe*npez-npe+ip ) then
          !         print*, iproc, array_zparbc(1,1,1), nek(1,1,1)
          !         print*, '              '
          call zsend(array_zparbc,npe*npez-npe+ip)
       end if
    end do
    !AVK debugging
    !call barrier
    
    do ip=0,npe*npez-NPE-1
       if (iproc==ip .and. iproc /= ip+npe) then
          call zreceive(array_zparbc,ip+NPE)
       end if
    end do
    do ip=npe*npez-NPE, npe*npez-1
       if (iproc==ip .and. iproc /= mod(ip,npe)) then
          call zreceive(array_zparbc,mod(ip,NPE))
       end if
    end do
    call barrier
    
    
    do i=1,nkx_par
       do j=1,nky
          nek_star(j,i,nlz_par) = nek(j,i,nlz_par) - dti/dz*(array_zparbc(j,i,4)-uekpar(j,i,nlz_par))&
               +new_etaz*dti*(array_zparbc(j,i,1) - 2.*nek(j,i,nlz_par) + nek(j,i,nlz_par-1))&
               /(dz**2)
          if (g_inc) then
             akpar_star(j,i,nlz_par) = akpar(j,i,nlz_par)-dti/(dz*(1.+kperp(j,i)**2*de**2))*&
                  (array_zparbc(j,i,2)-phik(j,i,nlz_par) - rhos**2*(array_zparbc(j,i,1)-nek(j,i,nlz_par)) &
                  -rhos**2*sqrt(2.)*(array_zparbc(j,i,5)-gk(j,i,nlz_par,gmin)))&
                  + new_etaz*dti*(array_zparbc(j,i,3) - 2.*akpar(j,i,nlz_par) +akpar(j,i,nlz_par-1))&
                  /(dz**2*(1.+kperp(j,i)**2*de**2))
             
             gk_star(j,i,nlz_par,gmin) = gk(j,i,nlz_par,gmin) &
                  -dti/dz*sqrt(gmin+1.)*rhos_de*(array_zparbc(j,i,6)-gk(j,i,nlz_par,gmin+1)) &
                  -dti/dz*sqrt(2.)*notanj*(array_zparbc(j,i,4)-uekpar(j,i,nlz_par)) &
                  + new_etaz_g*dti*(array_zparbc(j,i,5) - 2.*gk(j,i,nlz_par,gmin) + gk(j,i,nlz_par-1,gmin)) &
                  /(dz**2)
             
             do ng=gmin+1, ngtot-1
                gk_star(j,i,nlz_par,ng) = gk(j,i,nlz_par,ng) &
                     -dti/dz*(rhos_de)*sqrt(ng+1.0)*(array_zparbc(j,i,ng-gmin+6)-gk(j,i,nlz_par,ng+1)) &
                     -dti/dz*(rhos_de)*sqrt(ng*1.)*(1.-1./lambda*anj_kron(ng)*(1.-notanj)) &
                     *(array_zparbc(j,i,ng-gmin+4)-gk(j,i,nlz_par,ng-1)) &
                     +new_etaz_g*dti*(array_zparbc(j,i,ng-gmin+5) - 2.*gk(j,i,nlz_par,ng) + gk(j,i,nlz_par-1,ng)) &
                     /(dz**2)
             end do
          else
             akpar_star(j,i,nlz_par) = akpar(j,i,nlz_par) -dti/(dz*(1.+kperp(j,i)**2*de**2))*&
                  (array_zparbc(j,i,2)-phik(j,i,nlz_par) - rhos**2*(array_zparbc(j,i,1)-nek(j,i,nlz_par)))&
                  + new_etaz*dti*(array_zparbc(j,i,3) - 2.*akpar(j,i,nlz_par) +akpar(j,i,nlz_par-1))&
                  /(dz**2*(1.+kperp(j,i)**2*de**2))
          end if
       end do
    end do
    
    !AVK k = 1: diffusion term needs stuff(0) = stuff(nlz)
    !NFL, 16/12/2014: the following lines used to say nek_star, etc, which is wrong.
    array_zparbc(:,:,1)=nek(:,:,nlz_par)
    array_zparbc(:,:,2)=phik(:,:,nlz_par)
    array_zparbc(:,:,3)=akpar(:,:,nlz_par)
    array_zparbc(:,:,4)=uekpar(:,:,nlz_par)
    if(g_inc) then
       array_zparbc(:,:,5:ngtot-gmin+5)=gk(:,:,nlz_par,gmin:ngtot)
    end if
    
    do ip=0, npe*npez-1-npe
       if (iproc==ip .and. iproc /= ip+npe) then
          call zsend(array_zparbc,ip+npe)
       end if
    end do
    do ip=npe*npez-npe, npe*npez-1
       if (iproc==ip .and. iproc /= mod(ip,npe)) then
          call zsend(array_zparbc,mod(ip,npe))
       end if
    end do
    do ip=npe,npe*npez-1
       if (iproc==ip .and. iproc /= ip-npe) then
          call zreceive(array_zparbc,ip-npe)
       end if
    end do
    do ip=0, npe-1
       if (iproc==ip .and. iproc /= npe*npez-npe+ip) then
          call zreceive(array_zparbc,npe*npez-npe+ip)
       end if
    end do
    call barrier

    do i=1,nkx_par
       do j=1,nky
          nek_star(j,i,1) = nek(j,i,1) - dti/dz*(uekpar(j,i,2)-uekpar(j,i,1)) &
               +new_etaz*dti*(nek(j,i,2) - 2.*nek(j,i,1) + array_zparbc(j,i,1))/(dz**2)
          
          if (g_inc) then
             akpar_star(j,i,1) = akpar(j,i,1) - dti/(dz*(1+kperp(j,i)**2*de**2))*&
                  (phik(j,i,2)-phik(j,i,1)-rhos**2*(nek(j,i,2)-nek(j,i,1)) &
                  -rhos**2*sqrt(2.)*(gk(j,i,2,gmin)-gk(j,i,1,gmin)))&
                  +new_etaz*dti*(akpar(j,i,2) - 2.*akpar(j,i,1) + array_zparbc(j,i,3))&
                  /(dz**2*(1.+kperp(j,i)**2*de**2))
             
             gk_star(j,i,1,gmin) = gk(j,i,1,gmin) &
                  -dti/dz*sqrt(gmin+1.)*rhos_de*(gk(j,i,2,gmin+1)-gk(j,i,1,gmin+1)) &
                  -dti/dz*sqrt(2.)*notanj*(uekpar(j,i,2)-uekpar(j,i,1)) &
                  +new_etaz_g*dti*(gk(j,i,2,gmin) - 2.*gk(j,i,1,gmin) + array_zparbc(j,i,5))&
                  /(dz**2)
             
             do ng=gmin+1, ngtot-1
                gk_star(j,i,1,ng) = gk(j,i,1,ng) &
                     -dti/dz*(rhos_de)*sqrt(ng+1.0)*(gk(j,i,2,ng+1)-gk(j,i,1,ng+1)) &
                     -dti/dz*(rhos_de)*sqrt(ng*1.)*(1.-1./lambda*anj_kron(ng)*(1.-notanj)) &
                     *(gk(j,i,2,ng-1)-gk(j,i,1,ng-1))&
                     +new_etaz_g*dti*(gk(j,i,2,ng) - 2.*gk(j,i,1,ng) + array_zparbc(j,i,ng-gmin+5))&
                     /(dz**2)
             end do
          else
             akpar_star(j,i,1) = akpar(j,i,1) - dti/(dz*(1.+kperp(j,i)**2*de**2))*&
                  (phik(j,i,2)-phik(j,i,1) - rhos**2*(nek(j,i,2)-nek(j,i,1)))&
                  +new_etaz*dti*(akpar(j,i,2) - 2.*akpar(j,i,1) + array_zparbc(j,i,3))&
                  /(dz**2*(1.+kperp(j,i)**2*de**2))
          end if
       end do
    end do
    
    call PHI_POT(nek_star,phiK_star)
    
    do k=1,nlz_par
       do i=1,nkx_par
          do j=1,nky
             uekpar_star(j,i,k) = -kperp(j,i)**2*Akpar_star(j,i,k)
          end do
       end do
    end do
    !M closure:
    if(g_inc) then
       gk_star(:,:,:,ngtot) = gk_star(:,:,:,ngtot-1)
    end if
    
    !Now the corrector step of the MacCormack scheme:
    !AVK: 04/04/13, adding z diffusion term
    
    do k=2,nlz_par-1
       do i=1,nkx_par
          do j=1,nky
             neknew(j,i,k) = (nek(j,i,k) + nek_star(j,i,k))/2. &
                  - dti/(2.*dz)*(uekpar_star(j,i,k) - uekpar_star(j,i,k-1))&
                  +new_etaz*dti*(nek_star(j,i,k+1)-2.*nek_star(j,i,k)+nek_star(j,i,k-1))/(dz**2)
             
             if (g_inc) then
                akparnew(j,i,k) = (akpar(j,i,k) + akpar_star(j,i,k))/2. &
                     - dti/(2.*dz*(1.+kperp(j,i)**2*de**2))* &
                     (phik_star(j,i,k)-phik_star(j,i,k-1) &
                     - rhos**2*(nek_star(j,i,k)-nek_star(j,i,k-1)) &
                     - rhos**2*sqrt(2.)*(gk_star(j,i,k,gmin)-gk_star(j,i,k-1,gmin)))&
                     +new_etaz*dti*(akpar_star(j,i,k+1) - 2.*akpar_star(j,i,k) + akpar_star(j,i,k-1))&
                     /(dz**2*(1.+kperp(j,i)**2*de**2))
                
                gknew(j,i,k,gmin) = (gk(j,i,k,gmin) + gk_star(j,i,k,gmin))/2. &
                     - dti/(2.*dz)*sqrt(gmin+1.)*rhos_de*(gk_star(j,i,k,gmin+1)-gk_star(j,i,k-1,gmin+1)) &
                     - dti/(2.*dz)*sqrt(2.)*notanj*(uekpar_star(j,i,k)-uekpar_star(j,i,k-1))&
                     +new_etaz_g*dti*(gk_star(j,i,k+1,gmin) - 2.*gk_star(j,i,k,gmin) + gk_star(j,i,k-1,gmin))&
                     /(dz**2)
                
                do ng=gmin+1, ngtot-1
                   gknew(j,i,k,ng) = (gk(j,i,k,ng) + gk_star(j,i,k,ng))/2. &
                        - dti/(2.*dz)*(rhos_de)*sqrt(ng+1.0)*(gk_star(j,i,k,ng+1)-gk_star(j,i,k-1,ng+1)) &
                        - dti/(2.*dz)*(rhos_de)*sqrt(ng*1.)*(1.-1./lambda*anj_kron(ng)*(1.-notanj)) &
                        *(gk_star(j,i,k,ng-1)-gk_star(j,i,k-1,ng-1))&
                        + new_etaz_g*dti*(gk_star(j,i,k+1,ng) - 2.*gk_star(j,i,k,ng) + gk_star(j,i,k-1,ng))&
                        /(dz**2)
                end do
             else
                akparnew(j,i,k) = (akpar(j,i,k) + akpar_star(j,i,k))/2. &
                     - dti/(2.*dz*(1.+kperp(j,i)**2*de**2))* &
                     (phik_star(j,i,k)-phik_star(j,i,k-1) &
                     - rhos**2*(nek_star(j,i,k)-nek_star(j,i,k-1)))&
                     +new_etaz*dti*(akpar_star(j,i,k+1) - 2.*akpar_star(j,i,k) + akpar_star(j,i,k-1))&
                     /(dz**2*(1.+kperp(j,i)**2*de**2))
             end if
          end do
       end do
    end do
    
    !  !impose periodic BCs, i.e., stuff(0) = stuff(nlz)
    
    array_zparbc(:,:,1)=nek_star(:,:,nlz_par)
    array_zparbc(:,:,2)=phik_star(:,:,nlz_par)
    array_zparbc(:,:,3)=akpar_star(:,:,nlz_par)
    array_zparbc(:,:,4)=uekpar_star(:,:,nlz_par)
    if(g_inc) then
       array_zparbc(:,:,5:ngtot-gmin+5)=gk_star(:,:,nlz_par,gmin:ngtot)
    end if
    
    do ip=0, npe*npez-1-npe
       if (iproc==ip .and. iproc /= ip+npe) then
          call zsend(array_zparbc,ip+npe)
       end if
    end do
    do ip=npe*npez-npe, npe*npez-1
       if (iproc==ip .and. iproc /= mod(ip,npe)) then
          call zsend(array_zparbc,mod(ip,npe))
       end if
    end do
    do ip=npe,npe*npez-1
       if (iproc==ip .and. iproc /= ip-npe) then
          call zreceive(array_zparbc,ip-npe)
       end if
    end do
    do ip=0, npe-1
       if (iproc==ip .and. iproc /= npe*npez-npe+ip) then
          call zreceive(array_zparbc,npe*npez-npe+ip)
       end if
    end do
    call barrier
    
    do i=1,nkx_par
       do j=1,nky
          neknew(j,i,1) = (nek(j,i,1) + nek_star(j,i,1))/2. &
               - dti/(2.*dz)*(uekpar_star(j,i,1) - array_zparbc(j,i,4)) &
               +new_etaz*dti*(nek_star(j,i,2) - 2.*nek_star(j,i,1) + array_zparbc(j,i,1))&
               /(dz**2)
          
          if (g_inc) then
             akparnew(j,i,1) = (akpar(j,i,1) + akpar_star(j,i,1))/2. &
                  - dti/(2.*dz*(1.+kperp(j,i)**2*de**2))* &
                  (phik_star(j,i,1)-array_zparbc(j,i,2) &
                  - rhos**2*(nek_star(j,i,1)-array_zparbc(j,i,1)) &
                  - rhos**2*sqrt(2.)*(gk_star(j,i,1,gmin)-array_zparbc(j,i,5)))&
                  +new_etaz*dti*(akpar_star(j,i,2) - 2.*akpar_star(j,i,1) + array_zparbc(j,i,3))&
                  /(dz**2*(1.+kperp(j,i)**2*de**2)) 
             
             gknew(j,i,1,gmin) = (gk(j,i,1,gmin) + gk_star(j,i,1,gmin))/2. &
                  - dti/(2.*dz)*sqrt(gmin+1.)*rhos_de*(gk_star(j,i,1,gmin+1)-array_zparbc(j,i,6)) &
                  - dti/(2.*dz)*sqrt(2.)*notanj*(uekpar_star(j,i,1)-array_zparbc(j,i,4)) &
                  +new_etaz_g*dti*(gk_star(j,i,2,gmin) - 2.*gk_star(j,i,1,gmin) + array_zparbc(j,i,5))&
                  /(dz**2)
             
             do ng=gmin+1, ngtot-1
                gknew(j,i,1,ng) = (gk(j,i,1,ng) + gk_star(j,i,1,ng))/2. &
                     - dti/(2.*dz)*(rhos_de)*sqrt(ng+1.0)*(gk_star(j,i,1,ng+1)-array_zparbc(j,i,ng-gmin+6)) &
                     - dti/(2.*dz)*(rhos_de)*sqrt(ng*1.)*(1.-1./lambda*anj_kron(ng)*(1.-notanj)) & 
                     *(gk_star(j,i,1,ng-1)-array_zparbc(j,i,ng-gmin+4)) &
                     +new_etaz_g*dti*(gk_star(j,i,2,ng) - 2.*gk_star(j,i,1,ng) + array_zparbc(j,i,ng-gmin+5))&
                     /(dz**2)
             end do
          else
             akparnew(j,i,1) = (akpar(j,i,1) + akpar_star(j,i,1))/2. &
                  - dti/(2.*dz*(1.+kperp(j,i)**2*de**2))* &
                  (phik_star(j,i,1)-array_zparbc(j,i,2) &
                  - rhos**2*(nek_star(j,i,1)-array_zparbc(j,i,1))) &
                  +new_etaz*dti*(akpar_star(j,i,2) - 2.*akpar_star(j,i,1) + array_zparbc(j,i,3))&
                  /(dz**2*(1.+kperp(j,i)**2*de**2))
          end if
       end do
    end do

    !AVK k = nlz, z diffusion term requires stuff(nlz+1) = stuff(1)
    array_zparbc(:,:,1)=nek_star(:,:,1)
    array_zparbc(:,:,2)=phik_star(:,:,1)
    array_zparbc(:,:,3)=akpar_star(:,:,1)
    array_zparbc(:,:,4)=uekpar_star(:,:,1)
    if(g_inc) then
       array_zparbc(:,:,5:ngtot-gmin+5)=gk_star(:,:,1,gmin:ngtot)
    end if
    
    do ip=npe, npe*npez-1
       if (iproc==ip .and. iproc /= ip-npe) then
          !         print*, iproc, array_zparbc(1,1,1), nek(1,1,1)
          !         print*, '              '
          !write(*,*) "Before first zsend", ip, ip-npe
          call zsend(array_zparbc,ip-npe)
       end if
    end do
    do ip=0, npe-1
       if (iproc==ip .and. iproc /= npe*npez-npe+ip ) then
!          print*, iproc, array_zparbc(1,1,1), nek(1,1,1)
!          print*, '              '
          call zsend(array_zparbc,npe*npez-npe+ip)
       end if
    end do
    !AVK debugging
    !call barrier

    do ip=0,npe*npez-NPE-1
       if (iproc==ip .and. iproc /= ip+npe) then
          call zreceive(array_zparbc,ip+NPE)
       end if
    end do
    do ip=npe*npez-NPE, npe*npez-1
       if (iproc==ip .and. iproc /= mod(ip,npe)) then
          call zreceive(array_zparbc,mod(ip,NPE))
       end if
    end do
    call barrier

    do i=1,nkx_par
       do j=1,nky
          neknew(j,i,nlz_par) = (nek(j,i,nlz_par) + nek_star(j,i,nlz_par))/2. &
               - dti/(2.*dz)*(uekpar_star(j,i,nlz_par) - uekpar_star(j,i,nlz_par-1))&
               +new_etaz*dti*(array_zparbc(j,i,1) - 2.*nek_star(j,i,nlz_par) + nek_star(j,i,nlz_par-1))&
               /(dz**2)
         
          if (g_inc) then
             akparnew(j,i,nlz_par) = (akpar(j,i,nlz_par) + akpar_star(j,i,nlz_par))/2. &
                  - dti/(2.*dz*(1.+kperp(j,i)**2*de**2))* &
                  (phik_star(j,i,nlz_par)-phik_star(j,i,nlz_par-1) &
                  - rhos**2*(nek_star(j,i,nlz_par)-nek_star(j,i,nlz_par-1)) &
                  - rhos**2*sqrt(2.)*(gk_star(j,i,nlz_par,gmin)-gk_star(j,i,nlz_par-1,gmin)))&
                  + new_etaz*dti*(array_zparbc(j,i,3) - 2.*akpar_star(j,i,nlz_par) + akpar_star(j,i,nlz_par-1))&
                  /(dz**2*(1.+kperp(j,i)**2*de**2))

             gknew(j,i,nlz_par,gmin) = (gk(j,i,nlz_par,gmin) + gk_star(j,i,nlz_par,gmin))/2. &
                  - dti/(2.*dz)*sqrt(gmin+1.)*rhos_de*(gk_star(j,i,nlz_par,gmin+1)-gk_star(j,i,nlz_par-1,gmin+1)) &
                  - dti/(2.*dz)*sqrt(2.)*notanj*(uekpar_star(j,i,nlz_par)-uekpar_star(j,i,nlz_par-1))&
                  + new_etaz_g*dti*(array_zparbc(j,i,5) - 2.*gk_star(j,i,nlz_par,gmin) + gk_star(j,i,nlz_par-1,gmin))&
                  /(dz**2)
            
            
             do ng=gmin+1, ngtot-1
                gknew(j,i,nlz_par,ng) = (gk(j,i,nlz_par,ng) + gk_star(j,i,nlz_par,ng))/2. &
                     - dti/(2.*dz)*(rhos_de)*sqrt(ng+1.0)*(gk_star(j,i,nlz_par,ng+1)-gk_star(j,i,nlz_par-1,ng+1)) &
                     - dti/(2.*dz)*(rhos_de)*sqrt(ng*1.)*(1.-1./lambda*anj_kron(ng)*(1.-notanj)) &
                     *(gk_star(j,i,nlz_par,ng-1)-gk_star(j,i,nlz_par-1,ng-1))&
                     + new_etaz_g*dti*(array_zparbc(j,i,ng-gmin+5) - 2.*gk_star(j,i,nlz_par,ng) + gk_star(j,i,nlz_par-1,ng))&
                     /(dz**2)
             end do
          else
             akparnew(j,i,nlz_par) = (akpar(j,i,nlz_par) + akpar_star(j,i,nlz_par))/2. &
                  - dti/(2.*dz*(1.+kperp(j,i)**2*de**2))* &
                  (phik_star(j,i,nlz_par)-phik_star(j,i,nlz_par-1) &
                  - rhos**2*(nek_star(j,i,nlz_par)-nek_star(j,i,nlz_par-1))) &
                  + new_etaz*dti*(array_zparbc(j,i,3) - 2.*akpar_star(j,i,nlz_par) + akpar_star(j,i,nlz_par-1))&
                  /(dz**2*(1.+kperp(j,i)**2*de**2))
          end if
       end do
    end do

    call PHI_POT(neKnew,phiKnew)
    !   if(proc0) phiknew(0,0,:)=phizero_new(:)
    do k=1,nlz_par
       do i=1,nkx_par
          do j=1,nky
             uekparnew(j,i,k) = -kperp(j,i)**2*Akparnew(j,i,k)
          end do
       end do
    end do
    !M closure:
    if(g_inc) then 
       gknew(:,:,:,ngtot) = gknew(:,:,:,ngtot-1)
    end if


!update the variables:
    !   dadtk=(Akparnew-AKpar)/dti
    neK=neKnew
    AKpar=AKparnew 
    phiK=phiKnew
    ueKpar=uekparnew
    if(g_inc) gk=gknew


    deallocate(gk_star)
    deallocate(gknew)
    deallocate(array_zparbc)

  end subroutine z_step


!**************** Z-STEP *****************************************
  subroutine z_step2(phik,nek,uekpar,akpar,gk,dti)
    !NFL, 14/05/13
    !Do the z-step in characteristics form
    !Use TVDRK3 for the discretization of the d/dt op.
    !use UW7 (Pirozolli, JCP 02) for discretizing the fluxes

    !General idea: the z-step involves solving a set of eqs. in the form
    ! 
    !  du/dt = A du/dz, 
    !
    ! where u=(ne, Apar, g2, ..., gM)^T is the solution vector 
    ! and A is a non-diagonal matrix.
    ! To be able to use upwind schemes, we need to first diagonalize A, as follows:
    !
    ! P^{-1} du/dt = P^{-1} A P P^{-1} du/dz.
    !
    ! We define w = P^{-1}u and solve for P requiring that 
    !
    ! P^{-1}AP = D, a diagonal matrix. 
    !
    ! The equation for w is in characteristics form:
    !
    !   dw/dt  = D dw/dz ---- if D(j)>0 wave moves right, else moves left.
    !
    ! the diagonal of D are the eigenvalues of A; P is the matrix whose column vectors
    !are the eigenvectors of A.
    !To calculate the eigenvalues and eigenvectors, I use the lapack routine dgeev.
    !To invert P I use the lapack routines dgetri and dgetrf, see below.

    use constants
    !    use diag
    use mp
    use grid,   only:kperp,gama0,zz
    use Fluxes, only: leftflux, rightflux
    use Aux,    only: PHI_POT, comp_eigen, inv

    implicit none

    integer::i,j,k,ip,ng,l1,l2,i1,pp
    complex, dimension(nky,nkx_par,nlz_par) :: phik,nek,uekpar,akpar
!    complex,dimension(:,:,:,:)::gk
    complex, dimension(nky,nkx_par,nlz_par,gmin:ngtot)::gk
    complex,dimension(nky,nkx_par,nlz_par,dim_vec)::u
    complex,dimension(nky,nkx_par,nlz_par,dim_vec,0:3)::w
    complex,dimension(nky,nkx_par,nlz_par,dim_vec,0:2)::flux
    real,save,allocatable,dimension(:,:,:,:) :: A,P,invP,VL
    real,save,allocatable,dimension(:,:,:)::D,WR,WI
    real::dti
!    real::anj_kron
    logical,save::first=.true.
        

!    allocate(u(nky,nkx_par,nlz_par,dim_vec))
!    allocate(w(nky,nkx_par,nlz_par,dim_vec,0:3))
!    allocate(flux(nky,nkx_par,nlz_par,dim_vec,0:2))

    u=0.0
    w=0.0
    flux=0.0

    if(first) then

!!!!!!!!!!DEBUG !!!!!!!!!!!!!!
       
!       if(proc0) then
!          open(unit=655,file='debug.txt')
!          do k=1,nlz_par
!             write(655,'(4g16.8)') zz(k), real(Akpar(1,:,k))
!          end do
!          close(unit=655)
!       end if
!       call barrier
!       
!       do i=npe,npe*npez-1,npe
!          if(iproc==i) then
!             open(unit=655,file='debug.txt',position='append')
!             do k=1,nlz_par
!                write(655,'(4g16.8)') zz(k), real(Akpar(1,:,k))
!             end do
!             close(unit=655)
!          end if
!          call barrier
!       end do
       
!       if(proc0) then
!          open(unit=655,file='debug.txt',position='append')
!          write(655,*) '               '
!          write(655,*) '               '
!          close(unit=655)
!       end if
!***********************************
       !AVK: dim_vec is set appropriately for g_inc values
       !if(g_inc) then
       allocate(A(nky,nkx_par,dim_vec,dim_vec))
       allocate(P(nky,nkx_par,dim_vec,dim_vec))
       allocate(invP(nky,nkx_par,dim_vec,dim_vec))
       allocate(VL(nky,nkx_par,dim_vec,dim_vec))
       allocate(D(nky,nkx_par,dim_vec))
       allocate(WR(nky,nkx_par,dim_vec))
       allocate(WI(nky,nkx_par,dim_vec))

       A=0.0
       P=0.0
       invP=0.0
       VL=0.0
       D=0.0
       WR=0.0
       WI=0.0

       !else
       !   allocate(A(nky,nkx_par,2,2))
       !   allocate(P(nky,nkx_par,2,2))
       !   allocate(invP(nky,nkx_par,2,2))
       !   allocate(VL(nky,nkx_par,2,2))
       !   allocate(D(nky,nkx_par,2))
       !   allocate(WR(nky,nkx_par,2))
       !   allocate(WI(nky,nkx_par,2))
       !end if
       
       !form matrix A: there's a matrix A for each different value of kperp
!       A=0.0
       !AVK: A construction for anjor = false > 
       if(.not. anjor) then
          if (mod(iproc,npe)==0) then  !the kperp=0 cases need to be handled separately:
             if(g_inc) then
                do j=2,nky
                   A(j,1,1,2)=kperp(j,1)**2
                   A(j,1,2,1)=1./(1.+kperp(j,1)**2*de**2)*&
                        (-rhoi**2/(2*(gama0(kperp(j,1)**2*rhoi**2/2.)-1.)) + rhos**2)
                   A(j,1,2,3)=sqrt(2.0)*rhos**2/(1.+kperp(j,1)**2*de**2)
                   A(j,1,3,2)=sqrt(2.0)*kperp(j,1)**2
                   
                   ! AVK: Changed rhos/de to rhos_de
                   do i1=4,dim_vec
                      A(j,1,i1-1,i1)=-sqrt((i1-1)*1.0)*rhos_de
                      A(j,1,i1,i1-1)=-sqrt((i1-1)*1.0)*rhos_de  !A(j,1,i1-1,i1)
                   end do
                end do
                
                do i=2,nkx_par
                   do j=1,nky
                      A(j,i,1,2)=kperp(j,i)**2
                      A(j,i,2,1)=1./(1.+kperp(j,i)**2*de**2)*&
                           (-rhoi**2/(2*(gama0(kperp(j,i)**2*rhoi**2/2.)-1.)) + rhos**2)
                      A(j,i,2,3)=sqrt(2.0)*rhos**2/(1.+kperp(j,i)**2*de**2)
                      A(j,i,3,2)=sqrt(2.0)*kperp(j,i)**2
                      
                      !AVK: Changed rhos/de to rhos_de
                      do i1=4,dim_vec
                         A(j,i,i1-1,i1)=-sqrt((i1-1)*1.0)*rhos_de
                         A(j,i,i1,i1-1)=-sqrt((i1-1)*1.0)*rhos_de  !A(j,i,i1-1,i1)
                      end do
                      
                   end do
                end do
             else
                do j=2,nky
                   A(j,1,1,2)=kperp(j,1)**2
                   A(j,1,2,1)=1./(1.+kperp(j,1)**2*de**2)*&
                        (-rhoi**2/(2*(gama0(kperp(j,1)**2*rhoi**2/2.)-1.)) + rhos**2)
                end do
                
                do i=2,nkx_par
                   do j=1,nky
                      A(j,i,1,2)=kperp(j,i)**2
                      A(j,i,2,1)=1./(1.+kperp(j,i)**2*de**2)*&
                           (-rhoi**2/(2*(gama0(kperp(j,i)**2*rhoi**2/2.)-1.)) + rhos**2)
                   end do
                end do
             end if
             !compute the eigenvalues and eigenvectors
             do j=2,nky
                call comp_eigen(A(j,1,:,:),WR(j,1,:),WI(j,1,:),VL(j,1,:,:),P(j,1,:,:))
                call inv(P(j,1,:,:),invP(j,1,:,:))
             end do
             
             do i=2,nkx_par
                do j=1,nky
                   call comp_eigen(A(j,i,:,:),WR(j,i,:),WI(j,i,:),VL(j,i,:,:),P(j,i,:,:))
                   call inv(P(j,i,:,:),invP(j,i,:,:))
                end do
             end do
             
          else
             if(g_inc) then
                do i=1,nkx_par
                   do j=1,nky
                      A(j,i,1,2)=kperp(j,i)**2
                      A(j,i,2,1)=1./(1.+kperp(j,i)**2*de**2)*&
                           (-rhoi**2/(2*(gama0(kperp(j,i)**2*rhoi**2/2.)-1.)) + rhos**2)
                      A(j,i,2,3)=sqrt(2.0)*rhos**2/(1.+kperp(j,i)**2*de**2)
                      A(j,i,3,2)=sqrt(2.0)*kperp(j,i)**2
                      
                      !AVK: Changed rhos/de to rhos_de
                      do i1=4,dim_vec
                         A(j,i,i1-1,i1)=-sqrt((i1-1)*1.0)*rhos_de
                         A(j,i,i1,i1-1)=-sqrt((i1-1)*1.0)*rhos_de  !A(j,i,i1-1,i1)
                      end do
                      
                   end do
                end do
             else
                do i=1,nkx_par
                   do j=1,nky
                      A(j,i,1,2)=kperp(j,i)**2
                      A(j,i,2,1)=1./(1.+kperp(j,i)**2*de**2)*&
                           (-rhoi**2/(2*(gama0(kperp(j,i)**2*rhoi**2/2.)-1.)) + rhos**2)
                   end do
                end do
             end if
             !compute the eigenvalues and eigenvectors
             do i=1,nkx_par
                do j=1,nky
                   call comp_eigen(A(j,i,:,:),WR(j,i,:),WI(j,i,:),VL(j,i,:,:),P(j,i,:,:))
                   call inv(P(j,i,:,:),invP(j,i,:,:))
                end do
             end do
          end if
          
          !AVK: A construction for anjor = false ends 
          
          !AVK: A construction for anjor = true begins
          !Warning! the following chunk of code assumes rhoi<small_rhoi for anjor=true
       else 
          
          if (mod(iproc,npe)==0) then  !the kperp=0 cases need to be handled separately:
             if(g_inc) then
                do j=2,nky
                   A(j,1,1,2)=kperp(j,1)**2
                   A(j,1,2,1)= 1./(kperp(j,1)**2)
                   
                   ! AVK: Changed rhos/de to rhos_de
                   do i1=4,dim_vec
                      A(j,1,i1-1,i1)=-sqrt((i1-3)*1.0)*rhos_de
                      A(j,1,i1,i1-1)=-sqrt((i1-3)*1.0)*rhos_de  !A(j,1,i1-1,i1)
                   end do
                   !AVK: m = 1 has an additional term with lambda in it
                   A(j,1,4,3) = A(j,1,4,3) + rhos_de/lambda
                   
                end do
                
                do i=2,nkx_par
                   do j=1,nky
                      A(j,i,1,2)=kperp(j,i)**2
                      A(j,i,2,1)=1./(kperp(j,i)**2)
                      
                      !AVK: Changed rhos/de to rhos_de
                      do i1=4,dim_vec
                         A(j,i,i1-1,i1)=-sqrt((i1-3)*1.0)*rhos_de
                         A(j,i,i1,i1-1)=-sqrt((i1-3)*1.0)*rhos_de  !A(j,i,i1-1,i1)
                      end do
                      
                      !AVK: m = 1 has an additional term with lambda in it
                      A(j,i,4,3) = A(j,i,4,3) + rhos_de/lambda
                      
                   end do
                end do
             else
                do j=2,nky
                   A(j,1,1,2)=kperp(j,1)**2
                   A(j,1,2,1)=1./(kperp(j,1)**2)
                end do
                
                do i=2,nkx_par
                   do j=1,nky
                      A(j,i,1,2)=kperp(j,i)**2
                      A(j,i,2,1)=1./(kperp(j,i)**2)
                   end do
                end do
             end if
             !compute the eigenvalues and eigenvectors
             do j=2,nky
                call comp_eigen(A(j,1,:,:),WR(j,1,:),WI(j,1,:),VL(j,1,:,:),P(j,1,:,:))
                call inv(P(j,1,:,:),invP(j,1,:,:))
             end do
             
             do i=2,nkx_par
                do j=1,nky
                   call comp_eigen(A(j,i,:,:),WR(j,i,:),WI(j,i,:),VL(j,i,:,:),P(j,i,:,:))
                   call inv(P(j,i,:,:),invP(j,i,:,:))
                end do
             end do
             
          else
             if(g_inc) then
                do i=1,nkx_par
                   do j=1,nky
                      A(j,i,1,2)=kperp(j,i)**2
                      A(j,i,2,1)=1./(kperp(j,i)**2)
                      
                      !AVK: Changed rhos/de to rhos_de
                      do i1=4,dim_vec
                         A(j,i,i1-1,i1)=-sqrt((i1-3)*1.0)*rhos_de
                         A(j,i,i1,i1-1)=-sqrt((i1-3)*1.0)*rhos_de  !A(j,i,i1-1,i1)
                      end do
                      !AVK: m = 1 has an additional term with lambda in it
                      A(j,i,4,3) = A(j,i,4,3) + rhos_de/lambda
                      
                   end do
                end do
             else
                do i=1,nkx_par
                   do j=1,nky
                      A(j,i,1,2)=kperp(j,i)**2
                      A(j,i,2,1)=1./(kperp(j,i)**2)
                   end do
                end do
             end if
             !compute the eigenvalues and eigenvectors
             do i=1,nkx_par
                do j=1,nky
                   call comp_eigen(A(j,i,:,:),WR(j,i,:),WI(j,i,:),VL(j,i,:,:),P(j,i,:,:))
                   call inv(P(j,i,:,:),invP(j,i,:,:))
                end do
             end do
          end if
       end if
       
!       if(proc0) print*,cfl_frac*dz/maxval(abs(WR)), dti
       !       if(proc0) print*, '                '
       !       if(proc0) print*, 'check inverseP', matmul(P(2,1,:,:),invP(2,1,:,:))
       
       first=.false.
       
 !      if(proc0) then
       do ip=0,npe*npez-1
          if(iproc==ip) then
             open(unit=654, file='eigenvalues.dat',position='append')
             do l1=1,dim_vec
                do i=1,nkx_par
                   do j=1,nky
                      write(654,'(i4,g16.8,i4,g16.8)') iproc, kperp(j,i), l1, WR(j,i,l1)
                   end do
                end do
             end do
             close(unit=654)
          end if
          call barrier
       end do


    end if

    !now transform fields into characteristics:
    !the characteristics arrays is w(nky,nkx,nz,field_index,t_iter)
    !w=inv_P u
    u=0.0
    u(:,:,:,1)=nek(:,:,:)
    u(:,:,:,2)=Akpar(:,:,:)
    if(g_inc) then
       do ng=gmin,ngtot
          pp=ng-gmin+1
          u(:,:,:,2+pp)=gk(:,:,:,ng)
       end do
    end if

    w=0.0

!    do l1=1,dim_vec
!       do l2=1,dim_vec
!          do k=1,Nlz_par
!             do i=1,nkx_par
!                do j=1,nky
!                   w(j,i,k,l1,0)=w(j,i,k,l1,0)+invP(j,i,l1,l2)*u(j,i,k,l2)
!                end do
!             end do
!          end do
!       end do
!    end do

    do k=1,Nlz_par
       do i=1,nkx_par
          do j=1,nky
             w(j,i,k,:,0)=matmul(invP(j,i,:,:),u(j,i,k,:))
!NFL, 03/09/2013:
!couldn't get the following to work, so let's stick with the above line for now
!             call zgemm('N','N',dim_vec,1,dim_vec,1D0,invP(j,i,:,:),dim_vec,u(j,i,k,:),dim_vec,0D0,w(j,i,k,:,0),dim_vec)
          end do
       end do
    end do


!*********** test conversion u -> w -> u ****************
!    if(proc0) then
!       print*,'u before transf ', u(1,2,1,2)
!    end if
!    u=0.0
!    do l1=1,dim_vec
!       do l2=1,dim_vec
!          do k=1,Nlz_par
!             do i=1,nkx_par
!                do j=1,nky
!                   u(j,i,k,l1)=u(j,i,k,l1)+P(j,i,l1,l2)*w(j,i,k,l2,0)
!                end do
!             end do
!          end do
!       end do
!    end do
!    if(proc0) then
!       print*,'u after transf ', u(1,2,1,2)
!    end if
!    call barrier
!    stop

!*******************************************************
      
    do l1=1,dim_vec
       do i=1,nkx_par
          do j=1,nky
             if(WR(j,i,l1)>=0.0) then
                call rightflux(WR(j,i,l1),w(j,i,:,l1,0),flux(j,i,:,l1,0))

                w(j,i,:,l1,1) = w(j,i,:,l1,0)+dti*flux(j,i,:,l1,0)

                call rightflux(WR(j,i,l1),w(j,i,:,l1,1),flux(j,i,:,l1,1))

                w(j,i,:,l1,2)=0.75D0*w(j,i,:,l1,0)+0.25D0*w(j,i,:,l1,1)+&
                     0.25D0*dti*flux(j,i,:,l1,1)

                call rightflux(WR(j,i,l1),w(j,i,:,l1,2),flux(j,i,:,l1,2))

                w(j,i,:,l1,3)=1./3.*w(j,i,:,l1,0)+2./3.*w(j,i,:,l1,2)+&
                  2./3.*dti*flux(j,i,:,l1,2)             

             else

                call leftflux(WR(j,i,l1),w(j,i,:,l1,0),flux(j,i,:,l1,0))

                w(j,i,:,l1,1) = w(j,i,:,l1,0)+dti*flux(j,i,:,l1,0)

                call leftflux(WR(j,i,l1),w(j,i,:,l1,1),flux(j,i,:,l1,1))

                w(j,i,:,l1,2)=0.75D0*w(j,i,:,l1,0)+0.25D0*w(j,i,:,l1,1)+&
                     0.25D0*dti*flux(j,i,:,l1,1)

                call leftflux(WR(j,i,l1),w(j,i,:,l1,2),flux(j,i,:,l1,2))
                
                w(j,i,:,l1,3)=1./3.*w(j,i,:,l1,0)+2./3.*w(j,i,:,l1,2)+&
                     2./3.*dti*flux(j,i,:,l1,2)             
             end if

          end do
       end do
    end do

    !now compute the new u:
    u=0.0

!    do l1=1,dim_vec
!       do l2=1,dim_vec
!          do k=1,Nlz_par
!             do i=1,nkx_par
!                do j=1,nky
!                   u(j,i,k,l1)=u(j,i,k,l1)+P(j,i,l1,l2)*w(j,i,k,l2,3)
!                end do
!             end do
!          end do
!       end do
!    end do
    
    do k=1,Nlz_par
       do i=1,nkx_par
          do j=1,nky
             u(j,i,k,:)=matmul(P(j,i,:,:),w(j,i,k,:,3))
          end do
       end do
    end do
    

    !now need to decompose u into the fields:
    nek(:,:,:)=u(:,:,:,1)
    Akpar(:,:,:)=u(:,:,:,2)
    if(g_inc) then
       do ng=gmin,ngtot
          pp=ng-gmin+1
          gk(:,:,:,ng)=u(:,:,:,2+pp)
       end do
    end if

    do k=1,nlz_par
       do i=1,nkx_par
          do j=1,nky
             uekpar(j,i,k) = -kperp(j,i)**2*Akpar(j,i,k)
          end do
       end do
    end do

    call PHI_POT(nek,phik)

!    do i=0,npe*npez-1,npe
!       if(iproc==i) then
!          open(unit=655,file='debug.txt',position='append')
!          do k=1,nlz_par
!             write(655,'(4g16.8)') zz(k), real(Akpar(1,:,k))
!          end do
!          close(unit=655)
!       end if
!       call barrier
!    end do
!
!    if(proc0) then
!       open(unit=655,file='debug.txt',position='append')
!       write(655,*) '               '
!       write(655,*) '               '
!       close(unit=655)
!    end if

  end subroutine z_step2


end module Stepping
