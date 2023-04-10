module Initialize

  implicit none

contains

!****************************************************
  subroutine init_X_point(proc_X, y_loc_loc)

    use constants, only: nly_par, y_loc, z_loc, npe
    use grid,      only: r_variable, proc_id
    implicit none
    integer, intent(out) :: proc_X, y_loc_loc

    !!!proc_X=proc_id(r_variable,y_loc) + proc_id(z_loc)*NPE   !proc num that contains the X-point
    proc_X = proc_id(r_variable, y_loc, z_loc)
    y_loc_loc = y_loc - mod(proc_X, npe)*nly_par  !local j-index of the x-point
    
  end subroutine init_X_point

!****************************************************
  subroutine init_perturb(field)

    use constants
    use transforms, only: FFT2d_inv
    use forcing,    only: ran1
    use grid, only: xx, yy, zz
    use mp,   only: iproc

    implicit none
    
    real, intent(out), dimension(nlx, nly_par, nlz_par) :: field
    complex, dimension(nky, nkx_par, nlz_par) :: fieldk
    integer :: i, j, k, idum1, n, nmax
    nmax = nlz/4 !this is used to seed multiple kz's in 'allk' and 'alkp' routines
    
    fieldk=0.0
    field=0.0

 
    if(perturb_type=='none') then
    end if

    if(perturb_type=='gaus') then
       do k=1,nlz_par
          do j=1,NLy_par     
             do i=1,NLx
                field(i,j,k) = exp(-(zz(k)*2*pi/lz)**2)*&
                     exp(-(yy(j)*2*pi/ly)**2)*&
                     exp(-(xx(i)*2*pi/lx)**2)
             end do
          end do
       end do
    end if

    if (perturb_type == 'allk') then
    !TTR
#   ifdef gasca2d
       do k = 1, nlz_par
          do i = 1, nkx_par
             do j = 1, nky
                fieldk(j,i,k) = (1.0,1.0) * cos(2*pi*zz(k)/Lz)
             end do
          end do
          if(mod(iproc,npe)==0) fieldk(1,1,k) = 0.0
          call FFT2d_inv(fieldk(:, :, k), field(:, :, k))
       end do
#   elif defined(gasca3d)
       do k = 1, nlz_par
       idum1 = -1*( mod(iproc,npe)+1 )
          do i = 1, nkx_par
             do j = 1, nky
                do n = 1, nmax
                fieldk(j,i,k) = fieldk(j,i,k) + (ran1(idum1)-0.5d0 + (0.0,1.0)*(ran1(idum1)-0.5d0))*cos(2*n*pi*zz(k)/Lz+2*pi*ran1(idum1))
                end do
             end do
          end do
          if(mod(iproc, npe) == 0) fieldk(1, 1, k) = 0.0
       end do
       call FFT2d_inv(fieldk(:, :, :), field(:, :, :))
#   endif
    end if

    if (perturb_type == 'alkp') then
      !TTR
#   ifdef gasca2d
         do k = 1, nlz_par
            do i = 1, nkx_par
               do j = 1, nky
                  fieldk(j,i,k) = (1.0,1.0) * cos(2*pi*zz(k)/Lz)
               end do
            end do
            if(mod(iproc,npe)==0) fieldk(1,1,k) = 0.0
            call FFT2d_inv(fieldk(:, :, k), field(:, :, k))
         end do
#   elif defined(gasca3d)
         do k = 1, nlz_par
         idum1 = -1*( mod(iproc,npe)+1 )
            do i = 1, nkx_par
               do j = 1, nky
                  do n = 1, nmax
                  fieldk(j,i,k) = fieldk(j,i,k) + (ran1(idum1)+0.5d0 + (0.0,1.0)*(ran1(idum1)-0.5d0))*cos(2*n*pi*zz(k)/Lz+2*pi*ran1(idum1))
                  end do
               end do
            end do
            if(mod(iproc, npe) == 0) fieldk(1, 1, k) = 0.0
         end do
         call FFT2d_inv(fieldk(:, :, :), field(:, :, :))
         do k=1,nlz_par
            do j=1,NLy_par     
               do i=1,NLx
                  field(i,j,k) = field(i,j,k) + amp_ratio*cos(2*ky1*pi*yy(j)/Ly)*cos(2*kz1*pi*zz(k)/Lz)*cos(2*kx1*pi*xx(i)/Lx)
               end do
            end do
         end do
      if(twomodes) then
         do k=1,nlz_par
            do j=1,NLy_par     
               do i=1,NLx
                  field(i,j,k) = field(i,j,k) + amp_ratio*cos(2*ky2*pi*yy(j)/Ly)*cos(2*kz2*pi*zz(k)/Lz)*cos(2*kx2*pi*xx(i)/Lx)
               end do
            end do
         end do
      end if
# endif
  end if

    if (perturb_type == 'onmo') then
    !TTR
#   ifdef gasca2d
       do k = 1, nlz_par
          do i = 1, nkx_par
             do j = 1, nky
                field(i,j,k) = (1.0,1.0)*(2*pi*yy(j)/Ly)*(2*pi*xx(i)/Ly)*(20*pi*zz(k)/Ly)
             end do
          end do
          if(mod(iproc,npe)==0) fieldk(1,1,k) = 0.0
          call FFT2d_inv(fieldk(:, :, k), field(:, :, k))
       end do
#   elif defined(gasca3d)
       do k = 1, nlz_par
          do i = 1, nkx_par
             do j = 1, nky
                fieldk(1,1,10) =  (1.0,1.0)
             end do
          end do
          if(mod(iproc, npe) == 0) fieldk(1, 1, k) = 0.0
       end do
       call FFT2d_inv(fieldk(:, :, :), field(:, :, :))
#   endif
    end if


    if (perturb_type == 'onem') then
    !TTR
#   ifdef gasca2d
       do k = 1, nlz_par
          do i = 1, nkx_par
             do j = 1, nky
                fieldk(1,1,10) = (1.0,1.0)
             end do
          end do
          if(mod(iproc,npe)==0) fieldk(1,1,k) = 0.0
          call FFT2d_inv(fieldk(:, :, k), field(:, :, k))
       end do
#   elif defined(gasca3d)
       do k = 1, nlz_par
          do i = 1, nkx_par
             do j = 1, nky
                fieldk(5,5,k) =  (1.0,1.0)*cos(2*pi*zz(k)/Lz)
               !fieldk(-10,10,k) = -(1.0,1.0)*cos(2*pi*zz(k)/Lz)
             end do
          end do
          if(mod(iproc, npe) == 0) fieldk(1, 1, k) = 0.0
       end do
       call FFT2d_inv(fieldk(:, :, :), field(:, :, :))
#   endif
    end if


    if (perturb_type == 'cosi') then
       do k=1,nlz_par
          do j=1,NLy_par     
             do i=1,NLx
                !field(i,j,k) = cos(20*pi*xx(i)/Lx+20*pi*yy(j)/Ly)*cos(4*pi*zz(k)/Lz)+sin(20*pi*xx(i)/Lx+20*pi*yy(j)/Ly)*cos(4*pi*zz(k)/Lz)+cos(20*pi*xx(i)/Lx-20*pi*yy(j)/Ly)*cos(4*pi*zz(k)/Lz)+sin(20*pi*xx(i)/Lx-20*pi*yy(j)/Ly)*cos(4*pi*zz(k)/Lz)
                field(i,j,k) = cos(1*2*pi*yy(j)/Ly)*cos(24*2*pi*zz(k)/Lz)*cos(8*2*pi*xx(i)/Lx)
             end do
          end do
       end do
    end if

    if (perturb_type == 'rand') then
    !TTR
#   ifdef gasca2d
       do k = 1, nlz_par
          idum1 = -1*( mod(iproc,npe)+1 )
          do i = 1, nkx_par
             do j = 1, nky
                fieldk(j,i,k) = -(ran1(idum1)-0.5d0 + (0.0,1.0)*(ran1(idum1)-0.5d0))*cos(2*pi*zz(k)/Lz)
             end do
          end do
          call FFT2d_inv(fieldk(:, :, k), field(:, :, k))
       end do
#   elif defined(gasca3d)
       do k = 1, nlz_par
          idum1 = -1*( mod(iproc,npe)+1 )
          do i = 1, nkx_par
             do j = 1, nky
                fieldk(j,i,k) = -(ran1(idum1)-0.5d0 + (0.0,1.0)*(ran1(idum1)-0.5d0))*cos(2*pi*zz(k)/Lz)
             end do
          end do
       end do
       call FFT2d_inv(fieldk(:, :, :), field(:, :, :))
#   endif
    end if

    !this was used to test the TVDRK3UW7 for the JCP paper:
    !              Apar_perturb(i,j,k)=perturb_A*cos(xx(i))*&
    !                   0.5*tanh(100*(zz(k)-0.0)-tanh(100*(zz(k)-0.5)))

    

    ! A. Velb 6/9/22. Perturbs single m,n. m is poloidal mode number, n is z mode number. ONLY for 3D!
    if (perturb_type == 'polo') then
       do k=1,nlz_par
          do j=1,NLy_par
             do i=1,NLx                                                     
                ! In the following line, can set m and n (poloidal and 'toroidal' mode number). m is the first number in the first cos, n is the first
                ! number in the second cos
                field(i,j,k) = cos(1*atan2(yy(j),xx(i)))*cos(8*2*pi*zz(k)/Lz)
             end do
          end do
       end do
    end if

    if (perturb_type == 'pol2') then
       do k=1,nlz_par
          do j=1,NLy_par
             do i=1,NLx
                                                                                              
                field(i,j,k) = cos(1*atan2(yy(j),xx(i)))*cos(1*2*pi*zz(k)/Lz) + cos(1*atan2(yy(j),xx(i)))*cos(3*2*pi*zz(k)/Lz)
             end do
          end do
       end do
    end if
    
    ! A. Velberg 3/14/23. Perturbs phi to begin island coalescence using the fdev equilibrium
    if (perturb_type == 'coal') then
       do k=1,nlz_par
          do j=1,NLy_par
             do i=1,NLx
                                                                                              
                field(i,j,k) = -sin(xx(i))*sin(yy(j))
             end do
          end do
       end do
    end if
    

    field = perturb_amp*field

  end subroutine init_perturb

!*****************************************
  subroutine equilibrium(Apar_eq,  Akpar_eq, uepar_eq, uekpar_eq, &
       &                 Apar_eq_double_prime, Akpar_eq_double_prime, phi_eq)

    use constants,  only: nlx, nly_par, nlz_par, nky, nkx_par, equilib_type, pi, &
         &                Lx, Ly, Lz, phi0, a0, mm, nn
    use grid,       only: xx, yy, zz, kperp, kx
    use transforms, only: FFT2d_direct, FFT2d_inv

    implicit none
  
    real,    intent(out), dimension(nlx, nly_par, nlz_par) :: Apar_eq, uepar_eq, Apar_eq_double_prime
    real,    intent(out), dimension(nlx, nly_par, nlz_par) :: phi_eq
    complex, intent(out), dimension(nky, nkx_par, nlz_par) :: AKpar_eq, ueKpar_eq,Akpar_eq_double_prime

    integer :: i, j, k

   !AVK 14/7/13 >
    apar_eq = 0.0
    Akpar_eq = 0.0
    uepar_eq = 0.0
    uekpar_eq = 0.0
    apar_eq_double_prime = 0.0
    akpar_eq_double_prime = 0.0
    phi_eq = 0.0
  
    !<
  
    if(equilib_type=='tear') then
       !    j0=2./(Leq*Sqrt(Pi))
       do k = 1, nlz_par
          do j = 1, NLy_par
             do i = 1, NLx
                !this is the original one
                !           Apar_eq(i,j)=A0/cosh(xx(i)*2*Pi/lx)**2*(-1.+tanh(XX(i)+Lx/2.)**2+tanh(-XX(i)+Lx/2.)**2)
                !Numata's modification
                Apar_eq(i,j,k) = -A0/cosh(xx(i))**2*&
                     1./(2.*tanh(Pi)**2-tanh(2*Pi))*(tanh(xx(i)-Pi)**2+Tanh(xx(i)+Pi)**2-Tanh(2*Pi)**2)
                phi_eq(i,j,k)=0.0
                
             end do
          end do
       end do
    end if
  
    if(equilib_type=='tsin') then
       do k=1,nlz_par
          do j=1,NLy_par
             do i=1,NLx
  !              Apar_eq(i,j,k)=-A0/4.*cos(2*pi/Lx*xx(i)*4)
                Apar_eq(i,j,k)=-A0*cos(2*pi/Lx*xx(i))
                phi_eq(i,j,k)=0.0
             end do
          end do
       end do
    end if
  
    if(equilib_type=='kh00') then
       do k=1,nlz_par
          do j=1,NLy_par
             do i=1,NLx
                Apar_eq(i,j,k)=0.0
                phi_eq(i,j,k)=-PHI0/cosh(xx(i))**2*&
                     1./(2.*tanh(Pi)**2-tanh(2*Pi))*(tanh(xx(i)-Pi)**2+Tanh(xx(i)+Pi)**2-Tanh(2*Pi)**2)              
             end do
          end do
       end do
    end if

    if(equilib_type=='OT00') then
       !NFL, 20/03/13
       !Orszag-Tang equilib -- from Numata, JCP 10
       do k=1,nlz_par
          do j=1,NLy_par
             do i=1,NLx
!                Apar_eq(i,j,k)=cos(4*pi*xx(i)/Lx)+2.*cos(2*pi*yy(j)/ly)
!                phi_eq(i,j,k)=-2.0*(cos(2*pi*xx(i)/Lx)+cos(2*pi*yy(j)/ly))
                Apar_eq(i,j,k) = 0.5*cos(4*pi*xx(i)/Lx) + cos(2*pi*yy(j)/ly)
                phi_eq(i,j,k) = -cos(2*pi*xx(i)/Lx) - cos(2*pi*yy(j)/ly)
             end do
          end do
       end do
    end if

    if(equilib_type=='OTPD') then
       !NFL, 04/07/13
       !Orszag-Tang equilib -- from Paul Dellar, JCP 2013
       do k=1,nlz_par
          do j=1,NLy_par
             do i=1,NLx
                Apar_eq(i,j,k)=2.0*cos(2*pi*xx(i)/Lx) - cos(4*pi*yy(j)/ly)
                phi_eq(i,j,k)=2.0*cos(2*pi*xx(i)/Lx) - 2.0*sin(2*pi*yy(j)/ly)
             end do
          end do
       end do
    end if

    if(equilib_type=='OTBW') then
       !NFL, 04/07/13
       !Orszag-Tang equilib -- from Biskamp-Welter, PoF B 1989
       do k=1,nlz_par
          do j=1,NLy_par
             do i=1,NLx
                Apar_eq(i,j,k)=cos(4*pi*xx(i)/Lx + 2.3) + cos(2*pi*yy(j)/ly + 4.1)
                phi_eq(i,j,k)=cos(2*pi*xx(i)/Lx + 1.4) + cos(2*pi*yy(j)/ly + 0.5)
             end do
          end do
       end do
    end if

    if(equilib_type=='OT3D') then
       !NFL, 04/07/13                                                                                                          
       !Orszag-Tang equilib -- from Biskamp-Welter, PoF B 1989                                                    
       !modulated by a z perturb
       do k=1,nlz_par
          do j=1,NLy_par
             do i=1,NLx
                Apar_eq(i,j,k)=(cos(4*pi*xx(i)/Lx + 2.3) + cos(2*pi*yy(j)/ly + 4.1))*cos(2*pi*zz(k)/lz)
                phi_eq(i,j,k)=(cos(2*pi*xx(i)/Lx + 1.4) + cos(2*pi*yy(j)/ly + 0.5))*sin(2*pi*zz(k)/lz)
             end do
          end do
       end do
    end if
  
    if(equilib_type=='OT01') then
       !NFL, 20/03/13
       do k=1,nlz_par
          do j=1,NLy_par
             do i=1,NLx
                Apar_eq(i,j,k)=(cos(4*pi*xx(i)/Lx)+2.*cos(2*pi*yy(j)/ly))*cos(2*pi*zz(k)/lz)
                phi_eq(i,j,k)=-2.0*(cos(2*pi*xx(i)/Lx)+cos(2*pi*yy(j)/ly))*cos(2*pi*zz(k)/lz)
             end do
          end do
       end do
    end if

    !AVK 14/07/13>
    if(equilib_type/='none') then
    !TTR
#   ifdef gasca2d
       do k=1, nlz_par
          call FFT2d_direct (Apar_eq(:, :, k), AKpar_eq(:, :, k))
          !     print*,'equilil',maxval(abs(Apar_eq(:,:,k))),maxval(abs(Akpar_eq(:,:,k)))
          do i = 1, nkx_par
             do j = 1, nky
                uekpar_eq(j, i, k) = -kperp(j, i)**2 * AKpar_eq(j, i, k)
                AKpar_eq_double_prime(j, i, k) = -kx(i)**2 * AKpar_eq(j, i, k)
             end do
          end do
          call FFT2d_inv (AKpar_eq_double_prime(:, :, k), Apar_eq_double_prime(:, :, k))
          call FFT2d_inv (uekpar_eq(:, :, k), uepar_eq(:, :, k))
       end do
#   elif defined(gasca3d)
       call FFT2d_direct (Apar_eq(:, :, :), AKpar_eq(:, :, :))
       do i = 1, nkx_par
          do j = 1, nky
             uekpar_eq(j, i, :) = -kperp(j, i)**2 * AKpar_eq(j, i, :)
             AKpar_eq_double_prime(j, i, :) = -kx(i)**2 * AKpar_eq(j, i, :)
          end do
       end do
       call FFT2d_inv (AKpar_eq_double_prime(:, :, :), Apar_eq_double_prime(:, :, :))
       call FFT2d_inv (uekpar_eq(:, :, :), uepar_eq(:, :, :))
#   endif
    end if

    if(equilib_type=='isd1')then
    !Islands merging MZ                                                                                      
       do k=1,nlz_par
          do j=1, NLy_par
             do i=1, NLx
                Apar_eq(i,j,k)=a0*cos(mm*2*pi*xx(i)/Lx)*cos(nn*2*pi*yy(j)/ly)
                phi_eq(i,j,k)=0
             end do
          end do
       end do
    end if

    if(equilib_type=='gaus') then
       do k=1,nlz_par
          do j=1,NLy_par     
             do i=1,NLx
                Apar_eq(i,j,k) = a0*exp(-(yy(j)*2*pi*2./ly)**2)*&
                     exp(-(xx(i)*2*pi*2./lx)**2)
                phi_eq(i,j,k)=0
             end do
          end do
       end do

       call FFT2d_direct (Apar_eq(:, :, :), AKpar_eq(:, :, :)) ! This is needed to calc the perturbed spectrum. Not sure this will work if not done in 3D, so be mindful of this. A. Velb 5/26/22
    end if

    if(equilib_type=='coal')then
    ! Two flux ropes/islands coalescing (not an equilibrium)
        do k=1,nlz_par
            do j=1,NLy_par
                do i=1,NLx
                    Apar_eq(i,j,k) = a0*exp(-(2.*pi*1.5/lx*(xx(i)+.25*lx))**2 - (2.*pi*1.5/ly*yy(j))**2)+&
                    a0*exp(-(2.*pi*1.5/lx*(xx(i)-.25*lx))**2 - (2.*pi*1.5/ly*yy(j))**2)
                    phi_eq(i,j,k) = 0
                end do
            end do
        end do
    end if
    
    if(equilib_type=='fdev')then
    ! Fadeev equilibrium, modified to enforce perodicity in x direction, like the "tear" equil above
    ! the number inside the cosh and cos is k, number of islands. The factor outside cos is eps, measure of island width.
    !  AAV 3/2/23
        do k = 1, nlz_par
          do j = 1, NLy_par
             do i = 1, NLx
                !Apar_eq(i,j,k) = A0* (1/log(cosh(2*Pi)))*(log(cosh(2*Pi))-log(cosh(2.*xx(i)) + 0.5*cos(2.*yy(j))))*&
                    ! 1./(2.*tanh(Pi)**2-tanh(2*Pi))*(tanh(xx(i)-Pi)**2+tanh(xx(i)+Pi)**2-tanh(2*Pi)**2)
                !Apar_eq(i,j,k) = A0* (1/log(cosh(2*Pi)))*(log(cosh(2*Pi))-log(cosh(2.*xx(i)) + 0.5*cos(2.*yy(j))))*&                                                                                           
                     !(tanh(2*(xx(i)-Pi))**2+tanh(2*(xx(i)+Pi))**2-1) 
                Apar_eq(i,j,k) = A0* (1/log(cosh(2*Pi)))*(log(cosh(2*Pi))-log(cosh(2.*xx(i)) + 0.5*cos(2.*yy(j))))*&
                    (tanh((xx(i)-Pi))**2*tanh((xx(i)+Pi))**2)/tanh(Pi)**4
                phi_eq(i,j,k)=0.0
                
             end do
          end do
       end do
    end if
    
        
  end subroutine equilibrium


end module Initialize
