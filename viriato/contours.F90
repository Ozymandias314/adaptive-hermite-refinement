program contours

!Nuno Loureiro, 23/01/10
!post-processing program for REGK
!This takes the code output files containing the fields and gfields and generates 
!new files with fewer data points (set by nx_points and ny_points) so that 
!they can easily be handled by a visualization program. 
!The new files contain information about the whole grid.

!reads data times from manually created file "contours_files.txt"


!  use constants
  use mp, only: init_mp, finish_mp, proc0, broadcast  !, iproc, nproc, max_allreduce
  use diag
  implicit none

  integer::t,howlong
  integer,parameter::nt=3   !how many files to read
  real, dimension(nt)::file_times
 
  call init_mp

  if(proc0) then
     open(13,File='contours_files.txt')
     do t=1,nt
        read(13,*) file_times(t)
     end do
  end if
  call broadcast(file_times)
  
  do t=1,nt
     call predatasave(file_times(t),howlong)
     call filesname(file_times(t),howlong)
  end do
  call finish_mp
  stop
end program contours


!*********************************************************      
subroutine PHI_POT(nek,FIK)

use mp,only:proc0
use grid, only:kx, ky, kperp, gama0
use constants
implicit none

complex,dimension(nky,nkx_par)::nek,FIK
integer::i,j

if (roi.LT.10e-5)then
   if (proc0) then
      FIK(1,1)=0.0
      do j=2,nky 
         FIK(j,1)=-neK(j,1)/kperp(j,1)**2
      enddo
      
      do i=2,nkx_par
         do j=1,nky
            FIK(j,i)=-neK(j,i)/kperp(j,i)**2
         enddo
      enddo
   else
      do i=1,nkx_par
         do j=1,nky
            FIK(j,i)=-neK(j,i)/kperp(j,i)**2
         enddo
      enddo
   endif
else
   if (proc0) then
      FIK(1,1)=0.0
      do j=2,nky 
         FIK(j,1) = roi**2/2./(gama0(kperp(j,1)**2*roi**2/2.)-1)*neK(j,1) 
      enddo
      
      do i=2,nkx_par
         do j=1,nky
            FIK(j,i) = roi**2/2./(gama0(kperp(j,i)**2*roi**2/2.)-1)*neK(j,i) 
         enddo
      enddo
   else
      do i=1,nkx_par
         do j=1,nky
            FIK(j,i) = roi**2/2./(gama0(kperp(j,i)**2*roi**2/2.)-1)*neK(j,i) 
         enddo
      enddo
   endif
endif
end subroutine PHI_POT

!*********************************************************   
subroutine equilibrium(Apar_eq,Akpar_eq, uepar_eq,uekpar_eq, &
     &Apar_eq_double_prime, Akpar_eq_double_prime)

  use constants, only: nlx, nly_par, nky, nkx_par, A0, Pi, lx
  use grid, only: xx, kperp, kx
  use transforms, only: FFT2d_inv

  implicit none
  
  integer :: i, j
  real, dimension(nlx, nly_par), intent(out) :: Apar_eq, uepar_eq, Apar_eq_double_prime
  complex, dimension(nky, nkx_par), intent(out) :: AKpar_eq, ueKpar_eq,Akpar_eq_double_prime
  
  do j = 1, nly_par
     do i = 1, nlx
!        Apar_eq(i,j)=A0/cosh(xx(i)*2*Pi/lx)**2*&
!             &(-1+tanh(XX(i)+Lx/2.)**2+tanh(-XX(i)+Lx/2.)**2)
           Apar_eq(i,j)=-A0/cosh(xx(i))**2*&
                &1./(2.*tanh(Pi)**2-tanh(2*Pi))*(tanh(xx(i)-Pi)**2+Tanh(xx(i)+Pi)**2-Tanh(2*Pi)**2)
        !       Apar_eq(i,j)=cos(xx(i)*2*Pi/lx)
     end do
  end do
  
  !t.call FFT1d_direct (Apar_eq, AKpar_eq)
  call FFT2d_direct (Apar_eq, AKpar_eq)
  
  DO i = 1, NKx_par
     DO j = 1, NKy
        uekpar_eq(j, i) = -kperp(j, i)**2 * AKpar_eq(j, i)
        AKpar_eq_double_prime(j, i) = -kx(i)**2 * AKpar_eq(j, i)
     END DO
  END DO
  !t.call FFT1d_inv (AKpar_eq_double_prime, Apar_eq_double_prime)
  !t.call FFT1d_inv (uekpar_eq,uepar_eq)
  call FFT2d_inv (AKpar_eq_double_prime, Apar_eq_double_prime)
  call FFT2d_inv (uekpar_eq,uepar_eq)
  
end subroutine equilibrium


!***********************************************************
subroutine Boozer_Int_y(uepar_pert,file6,file7)
use constants
use mp, only: proc0,sum_reduce,iproc
use transforms
use grid, only:yy,kx,y_glob
implicit none

integer::i,j,k
real, DIMENSION(NLx,NLy_par) :: uepar_pert
COMPLEX, DIMENSION(NKy,NKx_par) :: uepar_pert_k 
COMPLEX, DIMENSION(Nly/2+1)::Int_ky
real,dimension(nly)::Int_y
character(len=*)::file6,file7

Int_ky=0.0
Int_y=0.0

do j=1, nly_par
   k=j+iproc*nly_par
   do i=2, nlx-1
      !integrate using trapezoidal rule:
      Int_y(k)=Int_y(k)+Lx/(Nlx*1.0)*(0.5*(uepar_pert(1,j)+uepar_pert(nlx,j))+&
           & uepar_pert(i,j))
   end do
end do

call sum_reduce(Int_y,0)

if (proc0) then
   call oneDfourndirect(Int_y, Int_ky, nly) 
   open (unit=6, file= file7)
   open (unit=7, file= file6)
   do j=1,nly/2+1
      write(6,40)  j-1, abs(Int_ky(j))
   end do
   do k=1,nly
      write(7,41) y_glob(k), Int_y(k)
   end do
40 format(i4, f16.8)
41 format(2f16.8)

   close(unit=6)
   close(unit=7)
end if

end subroutine Boozer_Int_y


!*********************************************************   
!subroutine filesname(time,length,F1,F2,F3,F4,F5,Fg,Fg1)
subroutine filesname(time,length) !,MMF1,F2,F3,F4,F5,Fg,Fg1)

  use constants
  use mp
  use transforms, only: FFT2d_direct, FFT2d_inv
  use diag
  use grid, only: kperp, xx, yy, kx, ky, gama0, proc_id
  use grid, only: r_variable, k_variable, init_grid

  implicit none 
  
  integer,intent(in)::length
  real, intent (in)::time
  character(len=len3+7+length+4+4)::F1
  character(len=len3+6+length+4+4)::F2
  character(len=len3+8+length+4+4)::F3
  character(len=len3+8+length+4+4)::F4
  character(len=len3+6+length+4+4)::F5
  character(len=len3+8+length+4+4)::Fg
  character(len=len3+9+length+4+4)::Fg1, Fg2, Fg3
  character(len=200)::MF1,MF2,MF3,MF4,MF5,MFg,MFg1, MFg2, MFg3
  character(len=200)::MMF1,MMF2,MMF3,MMF4,MMF5,MMFg,MMFg1,MMFg2, MMFg3
  integer:: i, j, k, ip, iaux, jaux, x_fact, y_fact, jglob, n, t, howlong
  complex, dimension(nky, nkx_par) :: nek, FIK, uekpar, Akpar, Akpar_eq
  complex, dimension(nky, nkx_par, 2:ngtot) :: Akpar_eq_double_prime, uekpar_eq, coll_diss_k
  real, dimension(nlx, nly_par) :: ne, vex, vey, vrosx, vrosy, tot_coll_diss
  real, dimension(nlx, nly_par) :: Apar, uepar, FI, Bx, By, epar
  real, dimension(nlx, nly_par) :: dxfi, dyfi, dxne, dyne, dxapar, dyapar
  real, dimension(nlx, nly_par, 2:ngtot) :: gm, coll_diss
  real, dimension(nlx, nly_par) :: Apar_eq, uepar_eq, Apar_eq_double_prime
  integer:: nx_points=384
  integer:: ny_points=384
  integer::nvpar,m
!  logical::boozer=.false.
  real::bxmax,bymax,bperp_max,x1,x2,x3
  real,dimension(2:ngtot)::xgm
!  real::vpar(ngtot)
!  real, dimension(nlx,nly_par,2:ngtot)::g
  complex,dimension(nky,nkx_par,2:ngtot)::gk
  real::gksq(2:ngtot)

  if (time .lt. 10) then
     write(MF1,'(a100,a7,f5.3,a4)') PATH,fieldsfile,time,".dat"
     MMF1=adjustl(MF1)
     F1=trim(MMF1)
     write(MF2,'(a100,a6,f5.3,a4)') PATH,"autot_",time,".dat"
     MMF2=adjustl(MF2)
     F2=trim(MMF2)
     write(MF3,'(a100,a8,f5.3,a4)') PATH,"finetot_",time,".dat"
     MMF3=adjustl(MF3)
     F3=trim(MMF3)
     write(MF4,'(a100,a8,f5.3,a4)') PATH,"flowtot_",time,".dat"
     MMF4=adjustl(MF4)
     F4=trim(MMF4)
     write(MF5,'(a100,a6,f5.3,a4)') PATH,"EBtot_",time,".dat"
     MMF5=adjustl(MF5)
     F5=trim(MMF5)
     write(MFg,'(a100,a8,f5.3,a4)') PATH,gfieldsfile,time,".dat"
     MMFg=adjustl(MFg)
     Fg=trim(MMFg)
     write(MFg1,'(a100,a5,f5.3,a4)') PATH,"gtot_",time,".dat"
     MMFg1=adjustl(MFg1)
     Fg1=trim(MMFg1)
     write(MFg2,'(a100,a5,f5.3,a4)') PATH,"gksp_",time,".dat"
     MMFg2=adjustl(MFg2)
     Fg2=trim(MMFg2)
     write(MFg3,'(a100,a9,f5.3,a4)') PATH,"distfunc_",time,".dat"
     MMFg3=adjustl(MFg3)
     Fg3=trim(MMFg3)
  end if
  if (time .ge. 10 .AND. time .lt. 100) then
     write(MF1,'(a100,a7,f6.3,a4)') PATH,"fields_",time,".dat"
     MMF1=adjustl(MF1)
     F1=trim(MMF1)
     write(MF2,'(a100,a6,f6.3,a4)') PATH,"autot_",time,".dat"
     MMF2=adjustl(MF2)
     F2=trim(MMF2)
     write(MF3,'(a100,a8,f6.3,a4)') PATH,"finetot_",time,".dat"
     MMF3=adjustl(MF3)
     F3=trim(MMF3)
     write(MF4,'(a100,a8,f6.3,a4)') PATH,"flowtot_",time,".dat"
     MMF4=adjustl(MF4)
     F4=trim(MMF4)
     write(MF5,'(a100,a6,f6.3,a4)') PATH,"EBtot_",time,".dat"
     MMF5=adjustl(MF5)
     F5=trim(MMF5)
     write(MFg,'(a100,a8,f6.3,a4)') PATH,gfieldsfile,time,".dat"
     MMFg=adjustl(MFg)
     Fg=trim(MMFg)
     write(MFg1,'(a100,a5,f6.3,a4)') PATH,"gtot_",time,".dat"
     MMFg1=adjustl(MFg1)
     Fg1=trim(MMFg1)
     write(MFg2,'(a100,a5,f6.3,a4)') PATH,"gksp_",time,".dat"
     MMFg2=adjustl(MFg2)
     Fg2=trim(MMFg2)
     write(MFg3,'(a100,a9,f6.3,a4)') PATH,"distfunc_",time,".dat"
     MMFg3=adjustl(MFg3)
     Fg3=trim(MMFg3)
  end if
  if (time .ge. 100 .AND. time .lt. 1000) then
     write(MF1,'(a100,a7,f7.3,a4)') PATH,"fields_",time,".dat"
     MMF1=adjustl(MF1)
     F1=trim(MMF1)
     write(MF2,'(a100,a6,f7.3,a4)') PATH,"autot_",time,".dat"
     MMF2=adjustl(MF2)
     F2=trim(MMF2)
     write(MF3,'(a100,a8,f7.3,a4)') PATH,"finetot_",time,".dat"
     MMF3=adjustl(MF3)
     F3=trim(MMF3)
     write(MF4,'(a100,a8,f7.3,a4)') PATH,"flowtot_",time,".dat"
     MMF4=adjustl(MF4)
     F4=trim(MMF4)
     write(MF5,'(a100,a6,f7.3,a4)') PATH,"EBtot_",time,".dat"
     MMF5=adjustl(MF5)
     F5=trim(MMF5)
     write(MFg,'(a100,a8,f7.3,a4)') PATH,gfieldsfile,time,".dat"
     MMFg=adjustl(MFg)
     Fg=trim(MMFg)
     write(MFg1,'(a100,a5,f7.3,a4)') PATH,"gtot_",time,".dat"
     MMFg1=adjustl(MFg1)
     Fg1=trim(MMFg1)
     write(MFg2,'(a100,a5,f7.3,a4)') PATH,"gksp_",time,".dat"
     MMFg2=adjustl(MFg2)
     Fg2=trim(MMFg2)
     write(MFg3,'(a100,a9,f7.3,a4)') PATH,"distfunc_",time,".dat"
     MMFg3=adjustl(MFg3)
     Fg3=trim(MMFg3)
  end if
  if (time .ge. 1000 .AND. time .lt. 10000) then
     write(MF1,'(a100,a7,f8.3,a4)') PATH,"fields_",time,".dat"
     MMF1=adjustl(MF1)
     F1=trim(MMF1)
     write(MF2,'(a100,a6,f8.3,a4)') PATH,"autot_",time,".dat"
     MMF2=adjustl(MF2)
     F2=trim(MMF2)
     write(MF3,'(a100,a8,f8.3,a4)') PATH,"finetot_",time,".dat"
     MMF3=adjustl(MF3)
     F3=trim(MMF3)
     write(MF4,'(a100,a8,f8.3,a4)') PATH,"flowtot_",time,".dat"
     MMF4=adjustl(MF4)
     F4=trim(MMF4)
     write(MF5,'(a100,a6,f8.3,a4)') PATH,"EBtot_",time,".dat"
     MMF5=adjustl(MF5)
     F5=trim(MMF5)
     write(MFg,'(a100,a8,f8.3,a4)') PATH,gfieldsfile,time,".dat"
     MMFg=adjustl(MFg)
     Fg=trim(MMFg)
     write(MFg1,'(a100,a5,f8.3,a4)') PATH,"gtot_",time,".dat"
     MMFg1=adjustl(MFg1)
     Fg1=trim(MMFg1)
     write(MFg2,'(a100,a5,f8.3,a4)') PATH,"gksp_",time,".dat"
     MMFg2=adjustl(MFg2)
     Fg2=trim(MMFg2)
     write(MFg3,'(a100,a9,f8.3,a4)') PATH,"distfunc_",time,".dat"
     MMFg3=adjustl(MFg3)
     Fg3=trim(MMFg3)
  end if
  if (time .ge. 10000 .AND. time .lt. 100000) then
     write(MF1,'(a100,a7,f9.3,a4)') PATH,"fields_",time,".dat"
     MMF1=adjustl(MF1)
     F1=trim(MMF1)
     write(MF2,'(a100,a6,f9.3,a4)') PATH,"autot_",time,".dat"
     MMF2=adjustl(MF2)
     F2=trim(MMF2)
     write(MF3,'(a100,a8,f9.3,a4)') PATH,"finetot_",time,".dat"
     MMF3=adjustl(MF3)
     F3=trim(MMF3)
     write(MF4,'(a100,a8,f9.3,a4)') PATH,"flowtot_",time,".dat"
     MMF4=adjustl(MF4)
     F4=trim(MMF4)
     write(MF5,'(a100,a6,f9.3,a4)') PATH,"EBtot_",time,".dat"
     MMF5=adjustl(MF5)
     F5=trim(MMF5)
     write(MFg,'(a100,a8,f9.3,a4)') PATH,gfieldsfile,time,".dat"
     MMFg=adjustl(MFg)
     Fg=trim(MMFg)
     write(MFg1,'(a100,a5,f9.3,a4)') PATH,"gtot_",time,".dat"
     MMFg1=adjustl(MFg1)
     Fg1=trim(MMFg1)
     write(MFg2,'(a100,a5,f9.3,a4)') PATH,"gksp_",time,".dat"
     MMFg2=adjustl(MFg2)
     Fg2=trim(MMFg2)
     write(MFg3,'(a100,a9,f9.3,a4)') PATH,"distfunc_",time,".dat"
     MMFg3=adjustl(MFg3)
     Fg3=trim(MMFg3)
  end if


!  call init_constants
  call init_grid
  call init_transforms

  if (nlx <= nx_points) nx_points=nlx
  if (nly <= ny_points) ny_points=nly

  x_fact=nlx/nx_points
  y_fact=nly/ny_points

  open(unit=1, file= F1)
  if(g_inc)  open(unit=11, file= Fg)
  !*****************************
  !READ FROM F
  if (g_inc) then
     do k=1,iproc*nly_par
        do i=1, nlx
           READ (1,*) x1,x2,x3
           READ (11,*) xgm(:)
        end do
     end do
     do k=iproc*nly_par+1, (iproc+1)*nly_par
        do i=1, nlx
           READ (1,*) Apar(i,k-iproc*nly_par), ne(i,k-iproc*nly_par), &
                epar(i,k-iproc*nly_par)
           READ (11,*) gm(i,k-iproc*nly_par,:)
        end do
     end do
     close (unit=1)
     close (unit=11)
  else
     do k=1,iproc*nly_par
        do i=1, nlx
           READ (1,*) x1,x2,x3
        end do
     end do
     do k=iproc*nly_par+1, (iproc+1)*nly_par
        do i=1, nlx
           READ (1,*) Apar(i,k-iproc*nly_par), ne(i,k-iproc*nly_par), &
                epar(i,k-iproc*nly_par)
        end do
     end do
     close (unit=1)
  end if
     !*****************

     !calculate secondary fields:
     !t.call FFT1d_direct (Apar, AKpar)
     call FFT2d_direct (Apar, AKpar)
     call convol(akpar,dxapar,dyapar)
     call b_field(dxapar,dyapar,bx,by,bxmax,bymax,bperp_max)

     !t.call FFT1d_direct (ne, nek)
     call FFT2d_direct (ne, nek)
     do i = 1, nkx_par
        do j = 1, nky
           uekpar(j, i) = -kperp(j, i)**2 * Akpar(j, i)
        end do
     end do
     !t.call FFT1d_inv (uekpar, uepar)
     call FFT2d_inv (uekpar, uepar)

     call PHI_POT(neK, FIK)
     call convol(fik,dxfi,dyfi)
     call convol(nek,dxne,dyne)
     call flows(dxfi,dyfi,dxne,dyne,vex,vey,vrosx,vrosy)
     !t.call FFT1d_inv(FIK, FI)
     call FFT2d_inv(FIK, FI)
     
     if (g_inc) then
!        call compute_g(0.0,gm,vpar,g)
        call compute_g(2.3,gm,Fg3)
        do nvpar = 2, ngtot
           !t.call FFT1d_direct(gm(:, :, nvpar), gk(:, :, nvpar))
           call FFT2d_direct(gm(:, :, nvpar), gk(:, :, nvpar))
        end do
!     
!        if (proc0) then
!           do j=1,nky
!              gksq(:)=gksq(:)+0.5*gk(j,1,:)*conjg(gk(j,1,:))
!           end do
!           do i=2,nkx_par
!              do j=1,nky
!                 gksq(:)=gksq(:)+gk(j,1,:)*conjg(gk(j,1,:))
!              end do
!        end do
!     else
!        do i=1,nkx_par
!           do j=1,nky
!              gksq(:)=gksq(:)+gk(j,1,:)*conjg(gk(j,1,:))
!           end do
!        end do
     end if
!     gksq=gksq/(nlx*nly*1.0)
!     do nvpar=2,ngtot
!        gksq(nvpar)=gksq(nvpar)/exp(-vpar(nvpar)**2)*pi**0.5
!     end do
!     call sum_reduce(gksq,0)
!     if (proc0) then
!        do nvpar=2,ngtot
!           write(765,*) vpar(nvpar),gksq(nvpar)
!        end do
!     end if
!  end if

  call mkspectra(time,Fg2,gk,coll_diss_k)
  do m = 2, ngtot
     !t.call FFT1d_inv(coll_diss_k(:, :, m), coll_diss(:, :, m))
     call FFT2d_inv(coll_diss_k(:, :, m), coll_diss(:, :, m))
  end do
  coll_diss(:, :, :) = coll_diss(:, :, :) * gm(:, :, :)
  tot_coll_diss = 0.0
  do j = 1, nly_par
     do i = 1, nlx
        do m = 2, ngtot
           tot_coll_diss(i, j) = tot_coll_diss(i, j) + coll_diss(i, j, m)
        end do
     end do
  end do


!     if (boozer) then
!        call equilibrium(Apar_eq,Akpar_eq, uepar_eq,uekpar_eq, &
!             & Apar_eq_double_prime, Akpar_eq_double_prime)
!        call Boozer_Int_y(uepar-uepar_eq,F6,F7)
!     end if
     
     do n=0,NPE-1
        if(iproc==n) then
           open (unit=2, file= F2, position = 'APPEND')
           open (unit=3, file= F3, position = 'APPEND')
           open (unit=4, file= F4, position = 'APPEND')
           open (unit=5, file= F5, position = 'APPEND')
           open (unit=6, file= Fg1, position = 'APPEND')
           do j=1,nly_par, y_fact
              do i=1,nlx, x_fact
                 write(2,35)  xx(i), yy(j), apar(i,j), uepar(i,j), tot_coll_diss(i,j)
                 write(3,33)  xx(i), yy(j), FI(i,j), ne(i,j)
                 write(4,34)  xx(i), yy(j), vex(i,j), vey(i,j), vrosx(i,j), vrosy(i,j)
                 write(5,35)  xx(i), yy(j), epar(i,j), bx(i,j),by(i,j)
                 write(6,35)  xx(i), yy(j), gm(i,j,2), gm(i,j,3), gm(i,j,4)
              end do
           end do
33         format(4f16.8)
34         format(6f16.8)
35         format(5f16.8)
           close(unit=2)
           close(unit=3)
           close(unit=4)
           close(unit=5)
           close(unit=6)
        end if
        call barrier
     end do

  
end subroutine filesname


!*********************************************************   
subroutine compute_g(loc,gm,filename)
  use constants
  use grid  !, only: r_variable, proc_id, xx
  use mp, only: iproc,proc0
  implicit none

  integer::ng,nvpar,i
  real::loc  !the y position of the x-cut
  real,dimension(ngtot+1)::vpar,w
  real, dimension(nlx,nly_par,2:ngtot)::gm
  real, dimension(nlx,ngtot+1)::g
  real::h(ngtot+1,-1:ngtot)
  integer::yi_glob,yi_loc,proc_loc
  character(len=*)::filename

  yi_glob=int(loc*nly/ly+1+nly/2)
  proc_loc=proc_id(r_variable,yi_glob)   
  yi_loc=yi_glob-mod(proc_loc,NPE)*nly_par  !local j-index 

  if (iproc==proc_loc) then
     open(unit=122,file=filename)
     call gauher2(vpar,w,ngtot+1)
     call hermite(vpar/sqrt(2.)*de/ros,h)

     g(:,:)=0.0
     do nvpar=1,ngtot+1
        do i=1,nlx
           do ng=2,ngtot
              g(i,nvpar)=g(i,nvpar) + gm(i,yi_loc,ng)*h(nvpar,ng)
           end do
        end do
        g(:,nvpar)=g(:,nvpar)*exp(-vpar(nvpar)**2)/pi**0.5
     end do

     do nvpar=1,ngtot+1
        do i=1, nlx
           write(122,'(3g16.8)') xx(i), vpar(nvpar), g(i,nvpar)
        end do
     end do
     close(unit=122)
  end if

end subroutine compute_g


!*********************************************************   
SUBROUTINE gauher(x,w,n)
!NFL:from Numerical Recipes;
!does not work for n> ~150 so I switched to gauher2
      INTEGER n,MAXIT
      REAL w(n),x(n)
      DOUBLE PRECISION EPS,PIM4
      PARAMETER (EPS=3.D-14,PIM4=.7511255444649425D0,MAXIT=10)
      INTEGER i,its,j,m
      DOUBLE PRECISION p1,p2,p3,pp,z,z1
      m=(n+1)/2
      do 13 i=1,m
        if(i.eq.1)then
          z=sqrt(float(2*n+1))-1.85575*(2*n+1)**(-.16667)
        else if(i.eq.2)then
          z=z-1.14*n**.426/z
        else if (i.eq.3)then
          z=1.86*z-.86*x(1)
        else if (i.eq.4)then
          z=1.91*z-.91*x(2)
        else
          z=2.*z-x(i-2)
        endif
        do 12 its=1,MAXIT
          p1=PIM4
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=z*sqrt(2.d0/j)*p2-sqrt(dble(j-1)/dble(j))*p3
11        continue
          pp=sqrt(2.d0*n)*p2
          z1=z
          z=z1-p1/pp
          if(abs(z-z1).le.EPS)goto 1
12      continue
        pause 'too many iterations in gauher'
1       x(i)=z
        x(n+1-i)=-z
        w(i)=2.d0/(pp*pp)
        w(n+1-i)=w(i)
13    continue
      return
 END


!*********************************************************   
subroutine hermite (x,h)
  use constants
  implicit none
  integer::i,j
  real::x(ngtot+1)
  real::h(ngtot+1,-1:ngtot)

  h(:,-1)=0.0
  h(:,0)=1.0/pi**(0.25)

  do j=0,ngtot-1
     do i=1,ngtot+1
        h(i,j+1)=x(i)*Sqrt(2.0/(j+1.0))*h(i,j)-sqrt(j/(j+1.0))*h(i,j-1)
     end do
  end do
  h = h*pi**(0.25)

end subroutine hermite


!*********************************************************   
subroutine mkspectra(time,Fg2,gk,coll_diss_k)
  use constants
  use mp
  use grid
  use transforms
  implicit none

  integer::m,mm, i,j,kp,t
  complex,dimension(nky,nkx_par,2:ngtot)::gk
  complex,dimension(nky,nkx_par,2:ngtot)::coll_diss_k
  real:: kpmax, dh, hyper_colls
  real:: dt,time,dummy(8)
  real::error=5e-3
  real, dimension(2:ngtot)::ek_gm
  real, dimension(:,:),allocatable::kshell
  real, dimension(:),allocatable::dum1,dum2
  character(len=*)::Fg2
  logical::found

  open(unit=12,file=Fg2)
  open(unit=14,file='timestep.dat')
  found=.false.
  do while (found==.false.)
     read(14,'(8g16.8)') dummy(:)
!     if (proc0) print*, dummy(8), time, abs(dummy(1)-time)/time
     if (abs(dummy(1)-time)/time<error) then
        dt=dummy(8)
        found=.true.
     end if
  end do

!  if(proc0) print*, time, dt
  dh = hyper_nu/dt/(NLx/2.)**(2*hyper_order)
  hyper_colls = hyperm_coef/dt/ngtot**(2*hyper_morder)
!  if (proc0) print*,dh,hyper_colls

  do m=2,ngtot
     if(proc0) then
        do j=1, nky
           ek_gm(m)=ek_gm(m)+0.5*abs(gk(j,1,m))**2
        end do
        do i=2,NKx_par
           do j=1,NKy
              ek_gm(m)=ek_gm(m)+abs(gk(j,i,m))**2
           end do
        end do
     else
        do i=1,NKx_par
           do j=1,NKy
              ek_gm(m)=ek_gm(m)+abs(gk(j,i,m))**2
           end do
        end do
     end if
  end do
  ek_gm(:)=ek_gm(:)/(nlx*nly*1.0)

  call sum_reduce(ek_gm,0)  
  if(mod(ngtot-1,2)==0)then 
     mm=ngtot
  else
     mm=ngtot-1
  end if

  kpmax=sqrt(ky(nky)**2.+((nkx-1.)*1.0)**2.)
  allocate(dum1(ceiling(kpmax)))
  allocate(dum2(ceiling(kpmax)))
!  allocate(kshell(ceiling(kpmax),2:ngtot))
  dum1=0.0
  dum2=0.0

  coll_diss_k=0.0

  do m=2,mm-1,2
     call k_energy(ceiling(kpmax),gk(:,:,m),gk(:,:,m+1),dum1,dum2)
!     kshell(:,m)=dum1(:)
!     kshell(:,m+1)=dum2(:)

     if(proc0) then 
        do kp=1,ceiling(kpmax)
           write (12,'(2i4,2g16.8)') m,kp, dum1(kp), &
                (dh*kp**(2*hyper_order) + nu_ei*m + &
                hyper_colls*m**(2*hyper_morder))*dum1(kp)
        end do
!        write(12,*) '                  '
        do kp=1,ceiling(kpmax)
           write (12,'(2i4,2g16.8)') m+1,kp, dum2(kp), &
                (dh*kp**(2*hyper_order) + nu_ei*m + &
                hyper_colls*m**(2*hyper_morder))*dum2(kp) 
        end do
!        write(12,*) '                  '
     end if

  end do

  do i=1,NKx_par
     do j=1,NKy
        coll_diss_k(j,i,2)=0.0  ! only has hyper-diffusion
        do m=3,ngtot
           coll_diss_k(j,i,m)=ros**2*(-nu_ei*m - &
                hyper_colls*m**(2*hyper_morder))*gk(j,i,m)  !*conjg(gk(j,i,m)
        end do
     end do
  end do

!  deallocate(dum1,dum2,kshell)
  deallocate(dum1,dum2)
  close(unit=12)
  close(unit=14)
end subroutine mkspectra


!*********************************************************   
subroutine k_energy(kpmax,akpar,fik,kshell_b, kshell_u)
    !energy per unit volume
    use constants
    use mp
    use grid
    implicit none
    
    COMPLEX, DIMENSION(NKy,NKX_PAR)::akpar,fik
    complex,dimension(nky,nkx)::sutot,sbtot
    complex, dimension (nky*nkx_par*NPE) :: work_u,work_b
    real, dimension(nky,nkx_par)::energy_u,energy_b
 !   real, dimension(:),allocatable::kshell_u,kshell_b
    real, dimension(kpmax)::kshell_u,kshell_b
    real:: ek_u,ek_b,time
    integer::i,j,p,kp,m,kpmax
    
    ek_u=0.0
    ek_b=0.0
    !the kx=0 modes only get added once
    !the kx diff 0 modes get added twice because of the reality condition
    if(proc0) then
       do j=1, nky
!          energy(j,1)=0.5*abs(kperp(j,1)*akpar(j,1))**2+abs(kperp(j,1)*fik(j,1))**2      
          energy_b(j,1)=0.5*abs(kperp(j,1)*akpar(j,1))**2     
          energy_u(j,1)=0.5*abs(kperp(j,1)*fik(j,1))**2     
          ek_b=ek_b+energy_b(j,1)
          ek_u=ek_u+energy_u(j,1)
       end do
       do i=2,NKx_par
          do j=1,NKy
!             energy(j,i)=abs(kperp(j,i)*akpar(j,i))**2+abs(kperp(j,i)*fik(j,i))**2
!             ek=ek+energy(j,i)
             energy_b(j,i)=abs(kperp(j,i)*akpar(j,i))**2     
             energy_u(j,i)=abs(kperp(j,i)*fik(j,i))**2     
             ek_b=ek_b+energy_b(j,i)
             ek_u=ek_u+energy_u(j,i)
          end do
       end do
    else
       do i=1,NKx_par
          do j=1,NKy
!             energy(j,i)=abs(kperp(j,i)*akpar(j,i))**2+abs(kperp(j,i)*fik(j,i))**2
!             ek=ek+energy(j,i)
             energy_b(j,i)=abs(kperp(j,i)*akpar(j,i))**2     
             energy_u(j,i)=abs(kperp(j,i)*fik(j,i))**2     
             ek_b=ek_b+energy_b(j,i)
             ek_u=ek_u+energy_u(j,i)
          end do
       end do
    end if
    ek_b=ek_b/(nlx*nly*1.0)
    ek_u=ek_u/(nlx*nly*1.0)

    energy_b=energy_b/(nlx*nly*1.0)
    energy_u=energy_u/(nlx*nly*1.0)
    call sum_reduce(ek_b,0)
    call sum_reduce(ek_u,0)

    work_u = 0.
    work_b = 0.
    m = 1+iproc*nky*nkx_par
    do i=1,nkx_par
       do j=1,nky
          work_b(m  )=energy_b(j,i)
          work_u(m  )=energy_u(j,i)
          m = m + 1
       end do
    end do
    call sum_reduce(work_b,0)
    call sum_reduce(work_u,0)
       
    if (proc0) then
       m=1
       do i=1,nkx
          do j=1,nky
             !          k=1+mod(k-1,nky*nkx_par)
             sutot(j,i)=work_u(m)
             sbtot(j,i)=work_b(m)
             m=m+1
          end do
       end do
       
!       kpmax=sqrt(ky(nky)**2.+((nkx-1.)*1.0)**2.)
!       allocate(kshell_u(ceiling(kpmax)))
!       allocate(kshell_b(ceiling(kpmax)))
       kshell_u=0.0
       kshell_b=0.0
       
!       do kp=1,ceiling(kpmax)
       do kp=1,kpmax
          do i=1,nkx
             do j=1,nky
                if (sqrt(((i-1.)*1.0)**2+ky(j)**2).GE. kp-1 .AND. &
                     & sqrt(((i-1.)*1.0)**2+ky(j)**2).LT.kp+1) then
                   kshell_b(kp)=kshell_b(kp)+sbtot(j,i)
                   kshell_u(kp)=kshell_u(kp)+sutot(j,i)
                end if
             end do
          end do
       end do
!       deallocate(kshell_u)
!       deallocate(kshell_b)
    end if
    call broadcast(kshell_u)
    call broadcast(kshell_b)
  end subroutine k_energy


!******************************
subroutine gauher2(x,w,order)
!NFL
!got from here:
!http://people.sc.fsu.edu/~jburkardt/f_src/gen_hermite_rule/gen_hermite_rule.html

!*****************************************************************************80
!
!
!  Discussion:
!
!    This program computes a generalized Gauss-Hermite quadrature rule 
!    and writes it to a file.
!
!    The user specifies:
!    * the ORDER (number of points) in the rule;
!    * ALPHA, the exponent of X;
!    * A, the center point;
!    * B, a scale factor;
!    * FILENAME, the root name of the output files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 )::  a=0.0D0
  real ( kind = 8 )::  alpha=0.0D0
  real ( kind = 8 )::  b=1.0D0
  real ( kind = 8 )  beta
  integer   ( kind = 4 )  iarg
  integer   ( kind = 4 )  iargc
  integer   ( kind = 4 )  ierror
  integer   ( kind = 4 )  kind
  integer   ( kind = 4 )  last
  integer   ( kind = 4 )  order
  real ( kind = 8 ),  dimension ( order ) :: w
  real ( kind = 8 ),  dimension ( order ) :: x

!
!  Initialize parameters.
!
  beta = 0.0D+00

  kind = 6

  call cgqf ( order, kind, alpha, beta, a, b, x, w )

end subroutine gauher2


!*********************************************************   
subroutine cdgqf ( nt, kind, alpha, beta, t, wts )

!*****************************************************************************80
!
!! CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
!
!  Discussion:
!
!    This routine computes all the knots and weights of a Gauss quadrature
!    formula with a classical weight function with default values for A and B,
!    and only simple knots.
!
!    There are no moments checks and no printing is done.
!
!    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Output, real ( kind = 8 ) T(NT), the knots.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) aj(nt)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) kind
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
  real ( kind = 8 ) zemu

  call parchk ( kind, 2 * nt, alpha, beta )
!
!  Get the Jacobi matrix and zero-th moment.
!
  call class_matrix ( kind, nt, alpha, beta, aj, bj, zemu )
!
!  Compute the knots and weights.
!
  call sgqf ( nt, aj, bj, zemu, t, wts )

  return
end


!*********************************************************   
subroutine cgqf ( nt, kind, alpha, beta, a, b, t, wts )

!*****************************************************************************80
!
!! CGQF computes knots and weights of a Gauss quadrature formula.
!
!  Discussion:
!
!    The user may specify the interval (A,B).
!
!    Only simple knots are produced.
!
!    Use routine EIQFS to evaluate this quadrature formula.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints, or
!    other parameters.
!
!    Output, real ( kind = 8 ) T(NT), the knots.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  integer ( kind = 4 ), allocatable :: mlt(:)
  integer ( kind = 4 ), allocatable :: ndx(:)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
!
!  Compute the Gauss quadrature formula for default values of A and B.
!
  call cdgqf ( nt, kind, alpha, beta, t, wts )
!
!  Prepare to scale the quadrature formula to other weight function with 
!  valid A and B.
!
  allocate ( mlt(1:nt) )

  mlt(1:nt) = 1

  allocate ( ndx(1:nt) )

  do i = 1, nt 
    ndx(i) = i
  end do

  call scqf ( nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b )

  deallocate ( mlt )
  deallocate ( ndx )

  return
end


!*********************************************************   
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character              c
  integer   ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end


!*********************************************************   
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

  character              c
  integer   ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end


!*********************************************************   
subroutine class_matrix ( kind, m, alpha, beta, aj, bj, zemu )

!*****************************************************************************80
!
!! CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
!
!  Discussion:
!
!    This routine computes the diagonal AJ and sub-diagonal BJ
!    elements of the order M tridiagonal symmetric Jacobi matrix
!    associated with the polynomials orthogonal with respect to
!    the weight function specified by KIND.
!
!    For weight functions 1-7, M elements are defined in BJ even
!    though only M-1 are needed.  For weight function 8, BJ(M) is
!    set to zero.
!
!    The zero-th moment of the weight function is returned in ZEMU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
!
!    Input, integer ( kind = 4 ) M, the order of the Jacobi matrix.
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Output, real ( kind = 8 ) AJ(M), BJ(M), the diagonal and subdiagonal
!    of the Jacobi matrix.
!
!    Output, real ( kind = 8 ) ZEMU, the zero-th moment.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a2b2
  real ( kind = 8 ) ab
  real ( kind = 8 ) aba
  real ( kind = 8 ) abi
  real ( kind = 8 ) abj
  real ( kind = 8 ) abti
  real ( kind = 8 ) aj(m)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) apone
  real ( kind = 8 ) beta
  real ( kind = 8 ) bj(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kind
  real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) temp
  real ( kind = 8 ) temp2
  real ( kind = 8 ) zemu

  temp = epsilon ( temp )

  call parchk ( kind, 2 * m - 1, alpha, beta )

  temp2 = 0.5D+00

  if ( 500.0D+00 * temp < abs ( ( r8_gamma ( temp2 ) )**2 - pi ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLASS_MATRIX - Fatal error!'
    write ( *, '(a)' ) '  Gamma function does not match machine parameters.'
    stop
  end if

  if ( kind == 1 ) then

    ab = 0.0D+00

    zemu = 2.0D+00 / ( ab + 1.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      abi = i + ab * mod ( i, 2 )
      abj = 2 * i + ab
      bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 2 ) then

    zemu = pi

    aj(1:m) = 0.0D+00

    bj(1) =  sqrt ( 0.5D+00 )
    bj(2:m) = 0.5D+00

  else if ( kind == 3 ) then

    ab = alpha * 2.0D+00
    zemu = 2.0D+00**( ab + 1.0D+00 ) * r8_gamma ( alpha + 1.0D+00 )**2 &
      / r8_gamma ( ab + 2.0D+00 )

    aj(1:m) = 0.0D+00
    bj(1) = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )
    do i = 2, m
      bj(i) = i * ( i + ab ) / ( 4.0D+00 * ( i + alpha )**2 - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 4 ) then

    ab = alpha + beta
    abi = 2.0D+00 + ab
    zemu = 2.0D+00**( ab + 1.0D+00 ) * r8_gamma ( alpha + 1.0D+00 ) &
      * r8_gamma ( beta + 1.0D+00 ) / r8_gamma ( abi )
    aj(1) = ( beta - alpha ) / abi
    bj(1) = 4.0D+00 * ( 1.0 + alpha ) * ( 1.0D+00 + beta ) &
      / ( ( abi + 1.0D+00 ) * abi * abi )
    a2b2 = beta * beta - alpha * alpha

    do i = 2, m
      abi = 2.0D+00 * i + ab
      aj(i) = a2b2 / ( ( abi - 2.0D+00 ) * abi )
      abi = abi**2
      bj(i) = 4.0D+00 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) &
        / ( ( abi - 1.0D+00 ) * abi )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 5 ) then

    zemu = r8_gamma ( alpha + 1.0D+00 )

    do i = 1, m
      aj(i) = 2.0D+00 * i - 1.0D+00 + alpha
      bj(i) = i * ( i + alpha )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 6 ) then

    zemu = r8_gamma ( ( alpha + 1.0D+00 ) / 2.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      bj(i) = ( i + alpha * mod ( i, 2 ) ) / 2.0D+00
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 7 ) then

    ab = alpha
    zemu = 2.0D+00 / ( ab + 1.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      abi = i + ab * mod ( i, 2 )
      abj = 2 * i + ab
      bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( kind == 8 ) then

    ab = alpha + beta
    zemu = r8_gamma ( alpha + 1.0D+00 ) * r8_gamma ( - ( ab + 1.0D+00 ) ) &
      / r8_gamma ( - beta )
    apone = alpha + 1.0D+00
    aba = ab * apone
    aj(1) = - apone / ( ab + 2.0D+00 )
    bj(1) = - aj(1) * ( beta + 1.0D+00 ) / ( ab + 2.0D+00 ) / ( ab + 3.0D+00 )
    do i = 2, m
      abti = ab + 2.0D+00 * i
      aj(i) = aba + 2.0D+00 * ( ab + i ) * ( i - 1 )
      aj(i) = - aj(i) / abti / ( abti - 2.0D+00 )
    end do

    do i = 2, m - 1
      abti = ab + 2.0D+00 * i
      bj(i) = i * ( alpha + i ) / ( abti - 1.0D+00 ) * ( beta + i ) &
        / ( abti**2 ) * ( ab + i ) / ( abti + 1.0D+00 )
    end do

    bj(m) = 0.0D+00
    bj(1:m) =  sqrt ( bj(1:m) )

  end if

  return
end


!*********************************************************   
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical              lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end


!*********************************************************   
subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to 
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine. 
!
!    It has been modified to produce the product Q' * Z, where Z is an input 
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
!    The changes consist (essentially) of applying the orthogonal 
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = 8 ) E(N), the subdiagonal entries of the 
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = 8 ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ), parameter :: itn = 30
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  real ( kind = 8 ) prec
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end


!*********************************************************   
subroutine parchk ( kind, m, alpha, beta )

!*****************************************************************************80
!
!! PARCHK checks parameters ALPHA and BETA for classical weight functions. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
!
!    Input, integer ( kind = 4 ) M, the order of the highest moment to
!    be calculated.  This value is only needed when KIND = 8.
!
!    Input, real ( kind = 8 ) ALPHA, BETA, the parameters, if required
!    by the value of KIND.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) m
  real ( kind = 8 ) tmp

  if ( kind <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARCHK - Fatal error!'
    write ( *, '(a)' ) '  KIND <= 0.'
    stop
  end if
!
!  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
!
  if ( 3 <= kind .and. alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARCHK - Fatal error!'
    write ( *, '(a)' ) '  3 <= KIND and ALPHA <= -1.'
    stop
  end if
!
!  Check BETA for Jacobi.
!
  if ( kind == 4 .and. beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARCHK - Fatal error!'
    write ( *, '(a)' ) '  KIND == 4 and BETA <= -1.0.'
    stop
  end if
!
!  Check ALPHA and BETA for rational.
!
  if ( kind == 8 ) then
    tmp = alpha + beta + m + 1.0D+00
    if ( 0.0D+00 <= tmp .or. tmp <= beta ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PARCHK - Fatal error!'
      write ( *, '(a)' ) '  KIND == 8 but condition on ALPHA and BETA fails.'
      stop
    end if
  end if

  return
end
function r8_epsilon ( )

!*****************************************************************************80
!
!! R8_EPSILON returns the R8 roundoff unit.
!
!  Discussion:
!
!    The roundoff unit is a number R which is a power of 2 with the
!    property that, to the precision of the computer's arithmetic,
!      1 < 1 + R
!    but
!      1 = ( 1 + R / 2 )
!
!    FORTRAN90 provides the superior library routine
!
!      EPSILON ( X )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_EPSILON, the R8 round-off unit.
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) d_test
  real ( kind = 8 ) r8_epsilon

  d = 1.0D+00
  d_test = 1.0D+00 + d / 2.0D+00

  do while ( 1.0D+00 < d_test )
    d = d / 2.0D+00
    d_test = 1.0D+00 + d / 2.0D+00
  end do

  r8_epsilon = d

  return
end
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
!
  implicit none

  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real ( kind = 8 ), parameter :: eps = 2.22D-16
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) sum
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 171.624D+00
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.0D+30
  real ( kind = 8 ), parameter :: xminin = 2.23D-308
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z

  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= 0.0D+00 ) then

      if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + 1.0D+00

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = 1.0D+00 / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < 1.0D+00 ) then

      z = y
      y = y + 1.0D+00
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - 1.0D+00

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + 1.0D+00
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - 0.5D+00 ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= 1.0D+00 ) then
    res = fact / res
  end if

  r8_gamma = res

  return
end


!*********************************************************   
subroutine scqf ( nt, t, mlt, wts, nwts, ndx, swts, st, kind, alpha, beta, a, &
  b )

!*****************************************************************************80
!
!! SCQF scales a quadrature formula to a nonstandard interval.
!
!  Discussion:
!
!    The arrays WTS and SWTS may coincide.
!
!    The arrays T and ST may coincide.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) T(NT), the original knots.
!
!    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
!
!    Input, real ( kind = 8 ) WTS(NWTS), the weights.
!
!    Input, integer ( kind = 4 ) NWTS, the number of weights.
!
!    Input, integer ( kind = 4 ) NDX(NT), used to index the array WTS.  
!    For more details see the comments in CAWIQ.
!
!    Output, real ( kind = 8 ) SWTS(NWTS), the scaled weights.
!
!    Output, real ( kind = 8 ) ST(NT), the scaled knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints.
!
  implicit none

  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts

  real ( kind = 8 ) a
  real ( kind = 8 ) al
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ) be
  real ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kind
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mlt(nt)
  integer ( kind = 4 ) ndx(nt)
  real ( kind = 8 ) p
  real ( kind = 8 ) shft
  real ( kind = 8 ) slp
  real ( kind = 8 ) st(nt)
  real ( kind = 8 ) swts(nwts)
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) temp
  real ( kind = 8 ) tmp
  real ( kind = 8 ) wts(nwts)

  temp = epsilon ( temp )

  call parchk ( kind, 1, alpha, beta )

  if ( kind == 1 ) then

    al = 0.0D+00
    be = 0.0D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 2 ) then

    al = -0.5D+00
    be = -0.5D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 3 ) then

    al = alpha
    be = alpha

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 4 ) then

    al = alpha
    be = beta

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 5 ) then

    if ( b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  B <= 0'
      stop
    end if

    shft = a
    slp = 1.0D+00 / b
    al = alpha
    be = 0.0D+00

  else if ( kind == 6 ) then

    if ( b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  B <= 0.'
      stop
    end if

    shft = a
    slp = 1.0D+00 / sqrt ( b )
    al = alpha
    be = 0.0D+00

  else if ( kind == 7 ) then

    al = alpha
    be = 0.0D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( kind == 8 ) then

    if ( a + b <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  A + B <= 0.'
      stop
    end if

    shft = a
    slp = a + b
    al = alpha
    be = beta

  else if ( kind == 9 ) then

    al = 0.5D+00
    be = 0.5D+00

    if ( abs ( b - a ) <= temp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCQF - Fatal error!'
      write ( *, '(a)' ) '  |B - A| too small.'
      stop
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  end if

  p = slp**( al + be + 1.0D+00 )

  do k = 1, nt

    st(k) = shft + slp * t(k)
    l = abs ( ndx(k) )

    if ( l /= 0 ) then
      tmp = p
      do i = l, l + mlt(k) - 1
        swts(i) = wts(i) * tmp
        tmp = tmp * slp
      end do
    end if

  end do

  return
end


!*********************************************************   
subroutine sgqf ( nt, aj, bj, zemu, t, wts )

!*****************************************************************************80
!
!! SGQF computes knots and weights of a Gauss Quadrature formula.
!
!  Discussion:
!
!    This routine computes all the knots and weights of a Gauss quadrature
!    formula with simple knots from the Jacobi matrix and the zero-th
!    moment of the weight function, using the Golub-Welsch technique.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( kind = 8 ) AJ(NT), the diagonal of the Jacobi matrix.
!
!    Input/output, real ( kind = 8 ) BJ(NT), the subdiagonal of the Jacobi 
!    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
!
!    Input, real ( kind = 8 ) ZEMU, the zero-th moment of the weight function.
!
!    Output, real ( kind = 8 ) T(NT), the knots.
!
!    Output, real ( kind = 8 ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) aj(nt)
  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
  real ( kind = 8 ) zemu
!
!  Exit if the zero-th moment is not positive.
!
  if ( zemu <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGQF - Fatal error!'
    write ( *, '(a)' ) '  ZEMU <= 0.'
    stop
  end if
!
!  Set up vectors for IMTQLX.
!
  t(1:nt) = aj(1:nt)

  wts(1) = sqrt ( zemu )
  wts(2:nt) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( nt, t, bj, wts )

  wts(1:nt) = wts(1:nt)**2

  return
end
