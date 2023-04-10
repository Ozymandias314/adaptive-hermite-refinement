module diag

  !TTR
  !. FPP macro to simplify RZG PERFLIB calls (performance-counter)
# include "perfmacro.h"
#define WITH_HAC_AGGREGATE 1
#define WITH_HAC_USE_TIMESTR 1
#if WITH_HAC
  !< MMH
  use mpi ! hac
  use hlst_adios_checkpoint ! hac
  use grid,       only: kglobal, jglobal
  use constants,  only: use_hac_checkpoint 
#endif

  implicit none

  private

  public :: Convol, Convol2
  public :: proc0_reads_to_all, proc0_reads_to_all_nogs, flows, b_field, injected
  public :: dataloadtest
#if WITH_HAC
  public :: dataloadtest_hac, datasavetest_hac
#endif
  public :: predatasave, datasavetest, kdatasave3d, test_of_test, cuts, k_energy, k_energy_pert, gamma_ky, helicity_test
  public :: diagnostics, zsteploss, gm_spectrum, gm_timetrace, zonalflows, zderivgk2
#if WITH_HAC
  !< MMH: HLST ADIOS CHECKPOINT AUXILIARY VARIABLES ****************
  !< The prefix 'hv_' stays for the initials of HAC-VIRIATO.
  integer,parameter :: hv_ad = 4 !< Array dimensions
  integer*4,public :: hv_GADA(hv_ad) !< Global Array Dimensions Array
  integer*4,public :: hv_GLBA(hv_ad) !< Globally indexed Local Bounds Array
  integer*4,public :: hv_err_ !< Error status variable
  character(len=3), parameter :: hv_fe='.bp'   !< hac file extension; "BP" is the name of the ADIOS format
  character(len=5), parameter :: hv_fl='_hac_' !< hac file label
  real,public,save :: hv_tstime = 0.0 !< total save time
  real :: hv_stime1, hv_stime2 !< save times: first and second intermediate
  character(len=*),parameter,public:: hv_lbl = "viriato_hac_info: "
  integer*4,public,save :: hv_vv = HAC_QUIET !< HAC Verbosity variable 
  logical :: debugging_d = .true.    
#endif

  !TTR
  interface Convol2
     module procedure Convol2_3
     module procedure Convol2_4
  end interface

contains
!*********************************************
   subroutine zderivgk2(gk,dzgk2)
  ! LMM 
  ! This subroutine computes the derivative (in real space) of the g2 moment of the distribution function, corresponding to 
  ! temperature perturbations. It is called in the k_energy routine; it is not called if g_inc is set to false
  ! The treatment of the information exchange between processors is based on the 'z_step' subroutine
   use constants
   use mp
   use grid, only:kperp

   integer :: i,j,k,ip
   complex, dimension(nky, nkx_par, nlz_par, gmin:ngtot):: gk
   complex, dimension(nky, nkx_par, nlz_par)::dzgk2
   complex, allocatable, dimension(:,:,:) :: array_zparbc
   !the first and only entry of the zparbc array is for gk2 

   allocate(array_zparbc(nky,nkx_par,1))
   !performing the for loop for which to information exchange with other processors is NOT necessary
   do k=2,Nlz_par-1
       do i=1,nkx_par
          do j=1,nky
             dzgk2(j,i,k)=(gk(j,i,k+1,gmin)-gk(j,i,k-1,gmin))/dz
          end do
       end do
   end do 
   !exchanging information with other processors, array to send information 'backward' 
    array_zparbc(:,:,1)=gk(:,:,1,gmin)
   
    do ip=npe, npe*npez-1
       if (iproc==ip .and. iproc /= ip-npe) then
          call zsend(array_zparbc,ip-npe)
       end if 
    end do
    do ip=0, npe-1
       if (iproc==ip .and. iproc /= npe*npez-npe+ip) then
          call zsend(array_zparbc,npe*npez-npe+ip)
       end if
    end do
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
    !basically right now array_zparbc corresponds to the datapoint nlz_par+1
    do i=1,nkx_par
       do j=1,nky
          dzgk2(j,i,nlz_par)=(array_zparbc(j,i,1) - gk(j,i,nlz_par-1,gmin))/dz
       end do 
    end do     
    !sending information forward
    array_zparbc(:,:,1)=gk(:,:,1,nlz_par)
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
   !Right now array_zparbc corresponds to the datapoint '0'
    do i=1,nkx_par
       do j=1,nky
          dzgk2(j,i,1)=(gk(j,i,2,gmin)-array_zparbc(j,i,1))/dz
       end do
    end do

end subroutine zderivgk2


!*********************************************
  subroutine kdatasave(iternumber,unitnumber,s1,s2,s3,s5,runname)
    
    use constants
    use mp
!s1=akpar
!s2=fik
!s3=nek
!s4=nik
!s5=uekpar
!s6=gyroavakpar
!s7=gyroavfik
    implicit none
    integer::iternumber,unitnumber,i,j,m
!  character(len=51)::filename 
    character(len=100)::filename 
    character(len=100), intent(in)::runname 
    complex,dimension(nky,nkx_par)::s1,s2,s3,s5
    complex,dimension(nky,nkx)::s1tot,s2tot,s3tot,s5tot
    complex,dimension (nky*nkx_par*4*NPE) :: work
!    complex, pointer, dimension(:) :: s1p, s2p, s3p, s4p, s5p


    if (proc0) then
       write(filename,'(a8,i6.6)') ktimefile,iternumber	
       filename = trim(runname)//"_"//trim(filename)
       open (unit=unitnumber, file=trim(filename))
       write(unitnumber,'(a5,a5,a15,a15,a15,a15,a15,a15,a15,a15)') 'i'&
            &,'j','reakpaar', 'refik','imfik','renek','imnek',&
            & 'reuekpar','imuekpar' 
    end if

    work = 0.
    m = 1+iproc*nky*nkx_par*4
    do j=1,nkx_par
       do i=1,nky
          work(m  )=s1(i,j)
          work(m+1)=s2(i,j)
          work(m+2)=s3(i,j)
          work(m+3)=s5(i,j)
          m = m + 4
       end do
    end do
    call sum_reduce(work,0)

    if (proc0) then
       m=1
       do j=1,nkx
          do i=1,nky
!          k=1+mod(k-1,nky*nkx_par)
          s1tot(i,j)=work(m)
          s2tot(i,j)=work(m+1)
          s3tot(i,j)=work(m+2)
          s5tot(i,j)=work(m+3)
          m=m+4
          

          write(unitnumber,321) i,j,real(s1tot(i,j)),aimag(s1tot(i,j)),&
               real(s2tot(i,j)),aimag(s2tot(i,j))!, real(s3tot(i,j)),&
              ! aimag(s3tot(i,j)),real(s5tot(i,j)),aimag(s5tot(i,j))

          end do
       end do
       close (unit=unitnumber) 	 
    end if

321 format (i5,    i5,    es15.6,    es15.6,   es15.6,   es15.6,  &
         & es15.6, es15.6,   es15.6,    es15.6)
  end subroutine kdatasave

!*****************************
subroutine test_of_test()
character(len=100)::test_name
test_name = 'this_is_a_test_of_a_test'
open(unit=1231, file=trim(test_name))
write(1231,*) 'test_completed'
close(unit=1231)
end subroutine test_of_test
!*********************************************!*********************************************
subroutine kdatasave3d(runname, time,length, akpar, nek, phik) !LMM
!(runname, time,length,apar,ne,epar,g)
use constants
use mp,only:iproc,barrier,proc0

implicit none
integer,intent(in)::length
integer::n,i,j,k
!character(len=len3+7+length+4+4)::filename1
!character(len=len3+8+length+4+4)::filename2
character(len=100)::filename1, filename2
character(len=100), intent(in) :: runname
real,intent(in)::time
complex,dimension(nkx_par,nky,nlz_par),intent(in)::akpar,nek,phik
character(len=20) :: time_str !MMH

if (time .lt. 10) then
        write(time_str,'(f5.3)') time
else if (time .ge. 10 .AND. time .lt. 100) then
        write(time_str,'(f6.3)') time
else if (time .ge. 100 .AND. time .lt. 1000) then
        write(time_str,'(f7.3)') time
else if (time .ge. 1000 .AND. time .lt. 10000) then
        write(time_str,'(f8.3)') time
else if (time .ge. 10000 .AND. time .lt. 100000) then
        write(time_str,'(f9.3)') time
end if

filename1 = trim(runname)//'_DEBUGGING_'//trim(time_str)//'.dat'

if(proc0) print*, 'inside kdatasave3d'

if(iproc==0) then
open (unit=1342, file= trim(filename1))
write(1342, '(8a15)') 'nky', 'nkx_par' ,'nlz_par', 'RE[akpar]','IM[akpar]','RE[nek]','IM[nek]', 'RE[fik]', 'IM[fik]'

do k=1, nlz_par
    do j=1, nkx_par
        do i=1, nky   
        write(1342,'(3i4,6g16.8)') i, j, k, REAL(akpar(i,j,k)), AIMAG(akpar(i,j,k)), REAL(nek(i,j,k)), AIMAG(nek(i,j,k)), REAL(phik(i,j,k)), AIMAG(phik(i,j,k))
        end do
    end do
end do
close(unit=1342)
end if
call barrier

do n=1,NPE*npez-1
if(iproc==n) then
open (unit=1342, file= trim(filename1), position = 'APPEND')!,form='unformatted')
write(1342, '(8a15)') 'nky', 'nkx_par', 'nlz_par', 'RE[akpar]', 'IM[akpar]', 'RE[nek]', 'IM[nek]', 'RE[fik]', 'IM[fik]'
do k=1, nlz_par
    do j=1,nkx_par
        do i=1, nky 
        write(1342,'(3i4,6g16.8)') i, j, k, REAL(akpar(i,j,k)), AIMAG(akpar(i,j,k)), REAL(nek(i,j,k)), AIMAG(nek(i,j,k)), REAL(phik(i,j,k)), AIMAG(phik(i,j,k))
        end do
    end do
end do
close(unit=1342)
end if
call barrier
end do

end subroutine kdatasave3d



!*********************************************!*********************************************
subroutine helicity_test(runname, time,length, akpar, nek) !LMM

use constants
use mp,only:iproc,barrier,proc0
use grid

implicit none
integer,intent(in)::length
integer::n,i,j,k
character(len=100)::filename1, filename2
character(len=100), intent(in) :: runname
real,intent(in)::time
complex,dimension(nkx_par,nky,nlz_par),intent(in)::akpar,nek
complex,dimension(nkx_par,nky,nlz_par)::ndya,dxapar,dyapar
character(len=20) :: time_str !MMH

if (time .lt. 10) then
        write(time_str,'(f5.3)') time
else if (time .ge. 10 .AND. time .lt. 100) then
        write(time_str,'(f6.3)') time
else if (time .ge. 100 .AND. time .lt. 1000) then
        write(time_str,'(f7.3)') time
else if (time .ge. 1000 .AND. time .lt. 10000) then
        write(time_str,'(f8.3)') time
else if (time .ge. 10000 .AND. time .lt. 100000) then
        write(time_str,'(f9.3)') time
end if

filename1 = trim(runname)//'_helicity_test_'//trim(time_str)//'.dat'

if(iproc==0) then
open (unit=1342, file= trim(filename1))
write(1342, '(a15,a15,a15,a15)') 'nky', 'nkx_par' ,'nlz_par', 'ndya'
do k=1, nlz_par
        do j=1, nky   
        ndya(j,1,k) = REAL(nek(j,1,k)*cmplx(0.0, 1.0)*ky(j)*akpar(j,1,k))
        write(1342,'(3i4,1g16.8)') j, 1, k, ndya(j,1,k)
        end do
    do i=2, nkx_par
        do j=1, nky   
        ndya(j,i,k) = 2*REAL(nek(j,i,k)**cmplx(0.0, 1.0)*ky(j)*akpar(j,i,k))
        write(1342,'(3i4,1g16.8)') j, i, k, ndya(j,i,k)
        end do
    end do
end do
close(unit=1342)
end if
call barrier

do n=1,NPE*npez-1
if(iproc==n) then
open (unit=1342, file= trim(filename1), position = 'APPEND')!,form='unformatted')
write(1342, '(a15,a15,a15,a15,a15,a15,a15)') 'nky', 'nkx_par', 'nlz_par', 'RE[akpar]', 'IM[akpar]', 'RE[nek]', 'IM[nek]'
do k=1, nlz_par
    do j=1,nkx_par
        do i=1, nky 
        write(1342,'(3i4,4g16.8)') i, j, k, REAL(akpar(i,j,k)), AIMAG(akpar(i,j,k)), REAL(nek(i,j,k)), AIMAG(nek(i,j,k))
        end do
    end do
end do
close(unit=1342)
end if
call barrier
end do

end subroutine helicity_test

!*********************************************





!*********************************************
  subroutine predatasave(time,howlong)
    implicit none
    real,intent(in)::time
    integer,intent(out)::howlong

!    if(time.lt.1) then
!       howlong=0
!    end if
!    if (time.ge.1 .AND. time .lt. 10) then
    if (time .lt. 10) then
       howlong=1
    end if
    if (time .ge. 10 .AND. time .lt. 100) then
       howlong=2
    end if
    if (time .ge. 100 .AND. time .lt. 1000) then
       howlong=3
    end if
    if (time .ge. 1000 .AND. time .lt. 10000) then
       howlong=4
    end if
    if (time .ge. 10000 .AND. time .lt. 100000) then
       howlong=5
    end if
  
  end subroutine predatasave


!*********************************************
#if WITH_HAC
  subroutine dataloadtest_hac(runname, apar,Apar_buff,ne,ne_buff,g,gdummy,epar,suffix_)

    use constants
    use mp, only: iproc, proc0, broadcast, barrier

    implicit none
    integer::n,i,j,k
    !character(len=len3+8+length+4+4)::filename2
    character(len=100)::bpfilename
    character(len=100), intent(in) :: runname
    real,dimension(:,:,:),intent(inout)::apar,ne,epar
    real::aa0,dti,x1,x2,x3,ek_u,ek_b, ek_e, sumleftover !  FIXME: why x1, x2, x3 ?
    real, DIMENSION(:,:,:):: Apar_buff, ne_buff
    real, DIMENSION(:,:,:,:):: g,gdummy
    character(len=20) :: suffix
    character(len=20), intent(in), optional :: suffix_

#if WITH_HAC
     if_hacload: if ( use_hac_checkpoint ) then
       !< MMH: Restart.
       hv_GLBA=(/1, jglobal(1), kglobal(1), gmin/)
       suffix=''
#if WITH_HAC_USE_TIMESTR
       if ( present ( suffix_ ) )  suffix=trim(suffix_)
#endif
#if WITH_HAC_AGGREGATE 
       bpfilename=trim(runname)//hv_fl//"epar+apar+ne+g"//trim(suffix)//hv_fe
#else
       bpfilename=trim(runname)//hv_fl//"apar"//trim(suffix)//hv_fe
#endif
       call hac_info (bpfilename, GADA_=hv_GADA, ierr_=hv_err_)
  hacr:if (hv_err_ == HAC_OK ) then
         if( hv_vv /= HAC_QUIET .and. proc0 ) then
           write(*,'(a)') hv_lbl//"Reading checkpoint file "//trim(bpfilename)//" with HAC ..."
           write(*,'(2(a,4(i0,x)),4(a,i0,x))') &
           &hv_lbl//" Dims: ", hv_GADA, " Bounds: ", hv_GLBA, " Procs: ", nly/nly_par, " x ", nlz/nlz_par,", ",gmin," : ",ngtot
         end if
#if WITH_HAC_AGGREGATE 
         call hv_hac_readm (epar, apar, ne, g, hv_GLBA(1:4), bpfilename)
#else
         if(g_inc) then
           call hac_read (g, hv_GLBA, trim(runname)//hv_fl//"g"//hv_fe)
         else
           call hv_hac_read3 (epar, hv_GLBA(1:3), trim(runname)//hv_fl//"epar"//hv_fe)
         end if
         call hv_hac_read3 (apar, hv_GLBA(1:3), trim(runname)//hv_fl//"apar"//hv_fe)
         call hv_hac_read3 (ne  , hv_GLBA(1:3), trim(runname)//hv_fl//"ne"//hv_fe)
#endif
       elseif (hv_err_ == HAC_ERR_NOFILE ) then
         if(proc0) write(*,'(a)') hv_lbl//"Skipping reading checkpoint file "//trim(bpfilename)//": it did not exist."
       else
         if(proc0) write(*,'(a)') hv_lbl//"hac_info() encountered a critical problem: Stopping here."
         stop
       end if hacr
       return
     endif if_hacload
#endif  
! WITH_HAC
  end subroutine dataloadtest_hac
#endif

!*********************************************
  subroutine dataloadtest(runname, apar,Apar_buff,ne,ne_buff,g,gdummy,epar,suffix)

    use constants
    use mp, only: iproc, proc0, broadcast, barrier

    implicit none
    integer::n,i,j,k
    !character(len=len3+8+length+4+4)::filename2
    character(len=100)::bpfilename
    character(len=100), intent(in) :: runname
    real,dimension(nlx,nly_par,nlz_par),intent(inout)::apar,ne,epar
    real::aa0,dti,x1,x2,x3,ek_u,ek_b, ek_e, sumleftover !  FIXME: why x1, x2, x3 ?
    real, DIMENSION(nlx,nly_par,nlz_par):: Apar_buff, ne_buff
    real, DIMENSION(nlx,nly_par,nlz_par,gmin:ngtot):: g,gdummy
    character(len=20), intent(in) :: suffix

#if WITH_HAC
    call dataloadtest_hac(runname,apar,Apar_buff,ne,ne_buff,g,gdummy,epar,suffix)
    return
#endif  
! WITH_HAC
     if (g_inc) then

       call proc0_reads_to_all(restart_file1,&
             restart_file2,Apar,Apar_buff,ne,ne_buff,g,gdummy)
     else
        call proc0_reads_to_all_nogs(restart_file1,&
             restart_file2,Apar,Apar_buff,ne,ne_buff)

     end if
  end subroutine dataloadtest


!*********************************************
#if WITH_HAC
!> MMH: Writes a few global 3D arrays using hac_write.
!       Assumes size(GAS1) == size(GAS2) == size(GAS3) == size(GAS4(:,:,:)).
!       FIXME: this is experimental, and coupled to hv_hac_readm.
  subroutine hv_hac_writem (GAS1, GAS2, GAS3, GAS4, LBA, fname, ierr_)
   use constants
   implicit none
!
   real,dimension(:,:,:),intent(in) :: GAS1, GAS2, GAS3 !< Global Array Slices
   integer*4,intent(in),dimension(:) :: LBA !< Local Bounds Array
   character(len=*),intent(in) :: fname !< File name
   integer*4,optional :: ierr_ !< Error variable
   real, allocatable, dimension(:,:,:,:) :: hv_b34d ! buffer from 3D to 4D
   real,dimension(:,:,:,:),intent(in) :: GAS4 !< Global Array Slice
!
   integer*4 :: tdims !< 
!
   tdims = 3
   if(g_inc) tdims = tdims + size(GAS4,4)
   allocate(hv_b34d(size(GAS1,1), size(GAS1,2), size(GAS1,3),tdims))
   hv_b34d(:,:,:,1) = GAS1
   hv_b34d(:,:,:,2) = GAS2
   hv_b34d(:,:,:,3) = GAS3
   if(g_inc) hv_b34d(:,:,:,4:) = GAS4
   call hac_write (hv_b34d, (/LBA(1), LBA(2), LBA(3), LBA(4)/), fname)
   deallocate(hv_b34d)
!
  end subroutine hv_hac_writem
!
!> MMH: Writes a global 3D array using hac_write.
  subroutine hv_hac_write3 (GAS, LBA, fname, ierr_)
   implicit none
!
   real,dimension(:,:,:),intent(in) :: GAS !< Global Array Slice
   integer*4,intent(in),dimension(:) :: LBA !< Local Bounds Array
   character(len=*),intent(in) :: fname !< File name
   integer*4,optional :: ierr_ !< Error variable
!
   call hac_write (reshape(GAS,(/size(GAS,1), size(GAS,2), size(GAS,3),1/)), (/LBA(1), LBA(2), LBA(3), 1/), fname)
!
  end subroutine hv_hac_write3
!
!> MMH: Reads a global 3D array using hac_read.
  subroutine hv_hac_read3(GAS, LBA, fname, ierr_)
   implicit none
!
   real,dimension(:,:,:),intent(inout) :: GAS !< Global Array Slice
   integer*4,intent(in),dimension(:) :: LBA !< Local Bounds Array
   character(len=*),intent(in) :: fname !< File name
   integer*4,optional :: ierr_ !< Error variable
   real, allocatable, dimension(:,:,:,:) :: hv_b34d ! buffer from 3D to 4D
!
   allocate(hv_b34d(size(GAS,1), size(GAS,2), size(GAS,3),1))
   call hac_read (hv_b34d, (/LBA(1), LBA(2), LBA(3), 1/), fname)
   GAS = hv_b34d(:,:,:,1)
   deallocate(hv_b34d)
!
  end subroutine hv_hac_read3
!
!> MMH: Reads all the global 3D arrays using hac_read.
!  FIXME: this is experimental, and coupled to hv_hac_writem.
  subroutine hv_hac_readm(GAS1, GAS2, GAS3, GAS4, LBA, fname, ierr_)
   use constants
   implicit none
!
   real,dimension(:,:,:),intent(inout) :: GAS1, GAS2, GAS3 !< Global Array Slices
   real,dimension(:,:,:,:),intent(inout) :: GAS4 !< Global Array Slice
   integer*4,intent(in),dimension(:) :: LBA !< Local Bounds Array
   character(len=*),intent(in) :: fname !< File name
   integer*4,optional :: ierr_ !< Error variable
   real, allocatable, dimension(:,:,:,:) :: hv_b34d ! buffer from 3D to 4D
!
   integer*4 :: tdims !< 
!
   tdims = 3
   if(g_inc) then
     tdims = tdims + size(GAS4,4)
   end if
   allocate(hv_b34d(size(GAS1,1), size(GAS1,2), size(GAS1,3), tdims))
   call hac_read (hv_b34d, (/LBA(1), LBA(2), LBA(3), LBA(4)/), fname)
   GAS2 = hv_b34d(:,:,:,2)
   GAS3 = hv_b34d(:,:,:,3)
   if(g_inc) then
     GAS4 = hv_b34d(:,:,:,4:)
   else
     GAS1 = hv_b34d(:,:,:,1)
   end if
   deallocate(hv_b34d)
!
  end subroutine hv_hac_readm
#endif

!!*********************************************
#if WITH_HAC
  subroutine datasavetest_hac(runname, time,length,apar,ne,epar,g,suffix_)

    use constants
    use mp,only:iproc,barrier,proc0

    implicit none
    integer,intent(in)::length
    integer::n,i,j,k
    !character(len=len3+7+length+4+4)::filename1
    !character(len=len3+8+length+4+4)::filename2
    character(len=100)::filename1, filename2
    character(len=100), intent(in) :: runname
    real,intent(in)::time
    real,dimension(:,:,:),intent(in)::apar,ne,epar
    real,dimension(:,:,:,:)::g
    character(len=20), intent(in), optional :: suffix_
    character(len=20) :: suffix

#if WITH_HAC
    if_hacsave: if ( use_hac_checkpoint ) then
      !< MMH: Checkpoint code
      call cpu_time(hv_stime1)
      suffix = ''
#if WITH_HAC_USE_TIMESTR
      if ( present ( suffix_ ) ) suffix = suffix_
#endif
#if WITH_HAC_AGGREGATE 
      filename1=trim(runname)//hv_fl//"epar+apar+ne+g"//trim(suffix)//hv_fe
#else
      filename1=trim(runname)//hv_fl//"apar"//trim(suffix)//hv_fe
#endif
      if( hv_vv /= HAC_QUIET .and. proc0 ) then
        write(*,'(a)') hv_lbl//"Writing checkpoint file "//trim(filename1)//" and others with HAC ..."
      end if
      hv_GLBA=(/1, jglobal(1), kglobal(1), gmin/)
#if WITH_HAC_AGGREGATE 
      call hv_hac_writem (epar, apar, ne, g, hv_GLBA(1:4), trim(filename1))
#else
      if(g_inc) then
        call hac_write (g, hv_GLBA, trim(runname)//hv_fl//"g"//hv_fe)
      endif
      call hv_hac_write3 (epar, hv_GLBA(1:3), trim(runname)//hv_fl//"epar"//hv_fe)
      call hv_hac_write3 (apar, hv_GLBA(1:3), trim(runname)//hv_fl//"apar"//hv_fe)
      call hv_hac_write3 (ne  , hv_GLBA(1:3), trim(runname)//hv_fl//"ne"//hv_fe)
#endif
      call cpu_time(hv_stime2)
      hv_tstime=hv_tstime+hv_stime2-hv_stime1
      return
    endif if_hacsave
#endif
! WITH_HAC
  end subroutine datasavetest_hac
#endif
! WITH_HAC


!!*********************************************
  subroutine datasavetest(runname, time,length,apar,ne,epar,g)

    use constants
    use mp,only:iproc,barrier,proc0

    implicit none
    integer,intent(in)::length
    integer::n,i,j,k
    !character(len=len3+7+length+4+4)::filename1
    !character(len=len3+8+length+4+4)::filename2
    character(len=100)::filename1, filename2
    character(len=100), intent(in) :: runname
    real,intent(in)::time
    real,dimension(nlx,nly_par,nlz_par),intent(in)::apar,ne,epar
    real,dimension(nlx,nly_par,nlz_par,gmin:ngtot)::g
    character(len=20) :: time_str !MMH
    if(proc0) write(1729,*) time
#if WITH_HAC
    time_str=''
    call file_name(runname, time,length,filename1,filename2,time_str)
    call datasavetest_hac(runname, time,length,apar,ne,epar,g,time_str)
    return
#endif
! WITH_HAC
    if(g_inc) then
       call file_name(runname, time,length,filename1,filename2)
       call proc0_recvs_and_writes(filename1,filename2,Apar,ne,g)

!        if(iproc==0) then
!           open (unit=12, file= trim(filename1))
!           open (unit=13, file= trim(filename2))
!           do k=1,nlz_par
!              do j=1,nly_par
!                 do i=1,nlx
!                    write(12,*)   apar(i,j,k), ne(i,j,k), epar(i,j,k)
!                    write(13,*)   g(i,j,k,:)
!                 end do
!              end do
!           end do
!           close(unit=12)
!           close(unit=13)
!        end if
!        call barrier
       
!        do n=1,NPE*npez-1
!           if (iproc==n) then
!              open (unit=12, file= trim(filename1), position = 'APPEND')!,form='unformatted')
!              open (unit=13, file= trim(filename2), position = 'APPEND')!,form='unformatted')
!              do k=1,nlz_par
!                 do j=1,nly_par
!                    do i=1,nlx
!                       write(12,*)   apar(i,j,k), ne(i,j,k), epar(i,j,k)
!                       write(13,*)   g(i,j,k,:)
!                    end do
!                 end do
!              end do
!              close(unit=12)
!              close(unit=13)
!           end if
!           call barrier
!        end do
    else
       call file_name(runname, time,length,filename1,filename2)
       !print*, 'about to output data with the new subroutine' - DEBUGGING LMM
       call proc0_recvs_and_writes_nogs(filename1,filename2,Apar,ne)
    end if 
  end subroutine datasavetest


!*********************************************
  subroutine cuts(unitnumber, time, proc_X,y_loc_loc,uepar,apar,ux,uy,runname)

    use constants
    use mp
    use grid, only: xx,yy, zz, idx_local, r_variable, k_variable, proc_id, y_glob
    use transforms
    
    implicit none
    integer,intent(in)::unitnumber,proc_X,y_loc_loc
    integer::i,ip,m,j,k,procz
    real, dimension(nlx,nly_par,nlz_par),intent (in)::uepar,apar,ux,uy
    real, dimension (nly) :: workuepar, workapar, workuy
    real, dimension (nlz) :: zworkuepar, zworkapar, zworkuy
    real,intent(in) ::time
    character(len=100)::filename1,filename2,filename3
    character(len=100),intent(in)::runname
    

    !cut in the x-direction:
    if (iproc==proc_X) then
       if (time.lt.10) then
          write(filename1,'(a5,f5.3,a4)') cutx,time,".dat"
       end if
       if (time .ge. 10 .AND. time .lt. 100) then
          write(filename1,'(a5,f6.3,a4)') cutx,time,".dat"
       end if
       if (time .ge. 100 .AND. time .lt. 1000) then
          write(filename1,'(a5,f7.3,a4)') cutx,time,".dat"
       end if
       if (time .ge. 1000 .AND. time .lt. 10000) then
          write(filename1,'(a5,f8.3,a4)') cutx,time,".dat"
       end if
	  
       filename1 = trim(runname)//"_"//trim(filename1)

       open (unit=unitnumber, file= trim(filename1))
       write(unitnumber, '(a10,4a15)') 'x' ,'uepar_x', 'apar', 'ux', 'uy'

       do i=1, nlx
          write(unitnumber,'(4g16.8)') xx(i),uepar(i,y_loc_loc,nlz_par), &
               & apar(i,y_loc_loc,nlz_par), ux(i,y_loc_loc,nlz_par)
       end do
       close (unit=unitnumber) 
    end if


!    !Now cut in the y-direction:
    
    workuepar=0.
    workapar=0.
    workuy=0.

    if(iproc<NPE) then
       m=1+iproc*nly_par
       do j=1,nly_par
          workuepar(m)=uepar(x_loc,j,1)
          workapar(m)=apar(x_loc,j,1)
          workuy(m)=uy(x_loc,j,1)
          m=m+1
       end do
    end if

    call sum_reduce(workuepar,0)
    call sum_reduce(workapar,0)
    call sum_reduce(workuy,0)
    if(proc0) then
       if (time.lt.10) then
          write(filename2,'(a5,f5.3,a4)') cuty,time,".dat"
       end if
       if (time .ge. 10 .AND. time .lt. 100) then
          write(filename2,'(a5,f6.3,a4)') cuty,time,".dat"
       end if
       if (time .ge. 100 .AND. time .lt. 1000) then
          write(filename2,'(a5,f7.3,a4)') cuty,time,".dat"
       end if
       if (time .ge. 1000 .AND. time .lt. 10000) then
          write(filename2,'(a5,f8.3,a4)') cuty,time,".dat"
       end if
       
       filename2 = trim(runname)//"_"//trim(filename2)

       open (unit=unitnumber+1, file= filename2)
       write(unitnumber+1, '(a10,a15)') 'y' ,'uepar_y', 'apar', 'uy'
       
       do j=1, nly
          write(unitnumber+1,'(4g16.8)') y_glob(j),workuepar(j), &
               & workapar(j), workuy(j)
       end do
       close (unit=unitnumber+1)
    end if

!now cut in the z-direction, at (0,0,z)

    zworkuepar=0.
    zworkapar=0.
    zworkuy=0.

    procz=proc_id(r_variable,y_loc,1)
    if(mod(iproc,npe)-procz==0) then
       m=(iproc-procz)/npe*nlz_par+1
       do k=1,nlz_par
          zworkuepar(m)=uepar(x_loc,y_loc_loc,k)
          zworkapar(m)=apar(x_loc,y_loc_loc,k)
          zworkuy(m)=uy(x_loc,y_loc_loc,k)
          m=m+1
       end do
    end if

    call sum_reduce(zworkuepar,0)
    call sum_reduce(zworkapar,0)
    call sum_reduce(zworkuy,0)

    if(proc0) then
       if (time.lt.10) then
          write(filename3,'(a5,f5.3,a4)') cutz,time,".dat"
       end if
       if (time .ge. 10 .AND. time .lt. 100) then
          write(filename3,'(a5,f6.3,a4)') cutz,time,".dat"
       end if
       if (time .ge. 100 .AND. time .lt. 1000) then
          write(filename3,'(a5,f7.3,a4)') cutz,time,".dat"
       end if
       if (time .ge. 1000 .AND. time .lt. 10000) then
          write(filename3,'(a5,f8.3,a4)') cutz,time,".dat"
       end if

       filename3 = trim(runname)//"_"//trim(filename3)

       open (unit=unitnumber+2, file= filename3)
       write(unitnumber+2, '(a10,4a15)') 'z' ,'uepar_z', 'apar', 'ux', 'uy'

       do k=1, nlz
          write(unitnumber+2,'(4g16.8)') lz/nlz*(k-1-nlz/2),zworkuepar(k), &
               zworkapar(k), zworkuy(k)
       end do
       close (unit=unitnumber+2) 
    end if

    call barrier
  end subroutine cuts


!*********************************************
  subroutine b_field(dxapar,dyapar,bx,by)
    !calculates the magnetic field from the Apar vector    
    use constants
    implicit none

    real, dimension(nlx,nly_par,nlz_par),intent(out):: bx, by
    real, dimension(nlx,nly_par,nlz_par)::dxapar,dyapar
    
    bx = dyapar
    by = -dxapar

  end subroutine b_field


!*********************************************
  subroutine epar_field(dadtk, epar)
!this calculates the parallel electric field.
    
    use constants
    use transforms, only: FFT2d_inv

    implicit none

    integer::e1,e2
    complex, dimension(nky, nkx_par) ::  ekpar, dadtk
    real, dimension(nlx, nly_par) :: epar, dadt
    
    do e1 = 1, nkx_par
       do e2 = 1, nky
            ekpar(e2, e1) = -dadtk(e2, e1)
        end do
    end do
    
    !TTR
    call FFT2d_inv (ekpar, epar)

  end subroutine epar_field


!*********************************************
  subroutine average_epar(epar,a,b,epar_avg)

    !returns spatial average of the electric field over some rectangular region
    !of width 2a * 2b
    use constants
    use mp
    implicit none

    real, dimension(nlx,nly_par):: epar
    real::a,b
    real::epar_avg
    real,allocatable,dimension(:)::work
    real,allocatable,dimension(:,:)::temp
    integer::i1,i2,j1_glob,j2_glob,i,j,m
 !   logical,save::first=.true.
    !first find the i and global j indices of lx,ly

    i1=nlx/lx*a+1+nlx/2
    i2=nlx/lx*(-a)+1+nlx/2

    j1_glob=nly/ly*b+1+nly/2
    j2_glob=nly/ly*(-b)+1+nly/2

    allocate (work((i1-i2+1)*nly))

    work = 0.
    m = 1+iproc*nly_par*(i1-i2+1)       
    do j=1,nly_par
       do i=i2,i1
          work(m)=epar(i,j)
          m = m + 1
       end do
    end do
    call sum_reduce(work,0)

    if (proc0) then
       allocate(temp(i1-i2+1,nly))
       temp=0.0
       epar_avg=0.0
       m=1
       do j=1,nly
          do i=1,i1-i2+1
             temp(i,j)=work(m)
             m=m+1
          end do
       end do
       do j=j2_glob,j1_glob
          do i=1,i1-i2+1
             epar_avg=epar_avg+temp(i,j)
          end do
       end do
       epar_avg=epar_avg/((i1-i2+1)*(j1_glob-j2_glob+1))
       deallocate(temp)
    end if
    deallocate(work)
  
  end subroutine average_epar


!*********************************************
  subroutine flows(dxfi,dyfi,dxne,dyne, vex, vey, vrhosx, vrhosy)
    !this calculates the e*b drift and the rhos drift
    use constants
    implicit none

    real,dimension(nlx,nly_par,nlz_par)::dxfi,dyfi,dxne,dyne
    real, dimension(nlx,nly_par,nlz_par),intent(out):: vex, vey, vrhosx, vrhosy
    
    if (rhoi .LT. small_rhoi) then
       vex = -dyfi
       vey = dxfi
       vrhosx = 0.0
       vrhosy = 0.0
    else
       vex = -dyfi
       vey = dxfi
       vrhosx = -rhos**2*dxne
       vrhosy = rhos**2*dyne
    end if
  end subroutine flows


!*********************************************
  subroutine tot_energy(vex,vey,bx,by,energy)
    !energy per unit volume
    use constants
    use mp
    use grid
    implicit none

    real, DIMENSION(NLx,NLy_PAR):: vex,vey,bx,by
    real::energy
    integer::i,j

    energy=0.0

    do j=1,nly_par
       do i=1,nlx
          energy=vex(i,j)**2+vey(i,j)**2+bx(i,j)**2+by(i,j)**2+energy
       end do
    end do

    call sum_reduce(energy,0)
    if (proc0) then 
       energy=0.5/(nlx*nly*1.0)*energy
    end if

  end subroutine tot_energy


!*********************************************
  subroutine k_energy(unitnumber,time,p,akpar,fik,nek,gk,dzgk2,braakpargk2,ek_u,ek_b,ek_e,hk_a,hk_t,ik_t,ik_a,ek_n,hk_aa,hk_ta,ik_ta,ik_aa,hk_1,hk_2,chk,runname)
    !energy per unit volume
    !helicity and energy injection - LMM 
    use constants
    use mp
    use grid
    implicit none
    
    COMPLEX, DIMENSION(NKy,NKX_PAR,nlz_par)::akpar,fik,nek,dzgk2,braakpargk2
    complex, dimension(nky,nkx_par,nlz_par,gmin:ngtot)::gk
    complex,dimension(nky,nkx)::sutot,sbtot,setot,shatot,shttot,sittot,siatot,sntot,shaatot,shtatot,sitatot,siaatot,shk1tot, shk2tot, schtot
    complex, dimension (nky*nkx_par*NPE) :: work_u,work_b, work_e, work_ha, work_ht, work_it, work_ia, work_n, work_haa, work_hta, work_ita, work_iaa, work_hk1, work_hk2, work_ch
    real, dimension(nky,nkx_par)::energy_u,energy_b, energy_e,energy_ha, energy_ht, energy_it, energy_ia, energy_n, energy_haa, energy_hta, energy_ita, energy_iaa, energy_hk1, energy_hk2, energy_ch ! LUIS, 28/4/15: kinetic, magnetic, electric energy
    real, dimension(:),allocatable::kshell_u,kshell_b,kshell_e,kshell_ha, kshell_ht, kshell_it, kshell_ia, kshell_n,kshell_haa, kshell_hta, kshell_ita, kshell_iaa, kshell_hk1, kshell_hk2, kshell_ch
    real:: ek_u,ek_b,ek_e,hk_a,hk_t,ik_t,ik_a,ek_n,hk_aa,hk_ta,ik_ta,ik_aa,hk_1,hk_2,time,kpmax,chk!,de,l
    integer::i,j,k,p,unitnumber,kp,m
    character(len=100)::filename1
    character(len=100),intent(in)::runname

    if (proc0 .AND. p .GT. 1  .AND. mod(p,save_energyfiles)==0) then
       if (time.lt.10) then
          write(filename1,'(a8,f5.3,a4)') energyk,time,".dat"
       end if
       if (time .ge. 10 .AND. time .lt. 100) then
          write(filename1,'(a8,f6.3,a4)') energyk,time,".dat"
       end if
       if (time .ge. 100 .AND. time .lt. 1000) then
          write(filename1,'(a8,f7.3,a4)') energyk,time,".dat"
       end if
       if (time .ge. 1000 .AND. time .lt. 10000) then
          write(filename1,'(a8,f8.3,a4)') energyk,time,".dat"
       end if
       if (time .ge. 10000 .AND. time .lt. 100000) then
          write(filename1,'(a8,f9.3,a4)') energyk,time,".dat"
       end if

       filename1 = trim(runname)//"_"//trim(filename1)
!NFL: 02/07/13 the following line is wrong!       
!       kpmax=sqrt(ky(nky)**2.+((nkx-1.)*1.0)**2.)
!       kpmax=sqrt((nky/2.*lx/ly)**2.+((nkx-1.)*1.0)**2.)
       kpmax=sqrt(ky(nky/2+1)**2.+((nkx-1.)*1.0)**2.)
       allocate(kshell_u(ceiling(kpmax)))
       allocate(kshell_b(ceiling(kpmax)))
       allocate(kshell_e(ceiling(kpmax)))
       allocate(kshell_ha(ceiling(kpmax)))
       allocate(kshell_ht(ceiling(kpmax)))
       allocate(kshell_it(ceiling(kpmax)))   
       allocate(kshell_ia(ceiling(kpmax)))    
       allocate(kshell_n(ceiling(kpmax)))
       allocate(kshell_haa(ceiling(kpmax)))
       allocate(kshell_hta(ceiling(kpmax)))
       allocate(kshell_ita(ceiling(kpmax)))   
       allocate(kshell_iaa(ceiling(kpmax)))     
       allocate(kshell_hk1(ceiling(kpmax)))   
       allocate(kshell_hk2(ceiling(kpmax)))   
       allocate(kshell_ch(ceiling(kpmax)))
       open (unit=unitnumber, file= trim(filename1))
       write(unitnumber, '(20a15)') 'kglobal(k)', 'kp' ,'kinetic', 'magnetic', 'electric', 'helicity_t', 'helicity_a', 'injection_t', 'injection_a', 'helicity_ta', 'helicity_aa', 'injection_ta', 'injection_aa', 'density', 'helicity_real', 'helicity_abs','cross_helicity'
     end if
    
    do k=1,nlz_par
       ek_u=0.0
       ek_e=0.0
       ek_b=0.0
       hk_a=0.0
       hk_t=0.0
       ik_t=0.0
       ik_a=0.0
       ek_n=0.0
       hk_aa=0.0
       hk_ta=0.0
       ik_ta=0.0
       ik_aa=0.0
       hk_1=0.0
       hk_2=0.0
       chk=0.0
       !the kx=0 modes only get added once
       !the kx diff 0 modes get added twice because of the reality condition
       if(mod(iproc,npe)==0) then
          do j=1, nky
             !          energy(j,1)=0.5*abs(kperp(j,1)*akpar(j,1))**2+abs(kperp(j,1)*fik(j,1))**2      
             energy_b(j,1)=0.5*abs(kperp(j,1)*akpar(j,1,k))**2     
             energy_e(j,1)=0.5*abs(kperp(j,1)*fik(j,1,k))**2     
             energy_u(j,1)=-1.0/rhoi**2*(gama0(kperp(j,1)**2*rhoi**2/2.)-1.)*fik(j,1,k)*conjg(fik(j,1,k))     
             energy_ha(j,1)= -real((rhoe_LTe/(de*rhos))* nek(j,1,k)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,1,k)))
             energy_ht(j,1)= real(nek(j,1,k)*conjg(sqrt(2.0)*dzgk2(j,1,k)))+real(nek(j,1,k)*sqrt(2.0)*conjg(braakpargk2(j,1,k)))
             energy_it(j,1) = -real(gk(j,1,k,gmin)/(sqrt(2.0))*conjg(cmplx(0.0,1.0)*ky(j)*fik(j,1,k)))
             energy_ia(j,1) = real((sqrt(3.0)*rhos)/(2.0*de)*gk(j,1,k,gmin+1)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,1,k)) )
             energy_n(j,1) = 0.5*abs(nek(j,1,k))**2
             energy_haa(j,1)= abs((rhoe_LTe/(de*rhos))* nek(j,1,k)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,1,k)))
             energy_hta(j,1)= abs(nek(j,1,k)*conjg(sqrt(2.0)*dzgk2(j,1,k)))+abs(nek(j,1,k)*sqrt(2.0)*conjg(braakpargk2(j,1,k)))
             energy_ita(j,1) = abs(gk(j,1,k,gmin)/(sqrt(2.0))*conjg(cmplx(0.0,1.0)*ky(j)*fik(j,1,k)))
             energy_iaa(j,1) = abs((sqrt(3.0)*rhos)/(2.0*de)*gk(j,1,k,gmin+1)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,1,k)))
             energy_hk1(j,1)= real((1+de**2*kperp(j,1)**2)*(fik(j,1,k)*conjg(akpar(j,1,k)))*(2*rhos**2/(rhoi**2))*(1.-gama0(kperp(j,1)**2*rhoi**2/2.)))
             energy_hk2(j,1)= abs((1+de**2*kperp(j,1)**2)*(fik(j,1,k)*conjg(akpar(j,1,k)))*(2*rhos**2/(rhoi**2))*(1.-gama0(kperp(j,1)**2*rhoi**2/2.)))    
             energy_ch(j,1) = real(-kperp(j,1)**2*fik(j,1,k)*(1+de**2*kperp(j,1)**2)*conjg(akpar(j,1,k)))
             ek_b=ek_b+energy_b(j,1)
             ek_u=ek_u+energy_u(j,1)
             ek_e=ek_e+energy_e(j,1)
             hk_a=hk_a+energy_ha(j,1)
             hk_t=hk_t+energy_ht(j,1)
             ik_t=ik_t+energy_it(j,1)
             ik_a=ik_a+energy_ia(j,1)
             ek_n=ek_n+energy_n(j,1)
             hk_aa=hk_aa+energy_haa(j,1)
             hk_ta=hk_ta+energy_hta(j,1)
             ik_ta=ik_ta+energy_ita(j,1)
             ik_aa=ik_aa+energy_iaa(j,1)
             hk_1=hk_1+energy_hk1(j,1)
             hk_2=hk_2+energy_hk2(j,1)
             chk=chk+energy_ch(j,1)
          end do
          do i=2,NKx_par
             do j=1,NKy
                !             energy(j,i)=abs(kperp(j,i)*akpar(j,i))**2+abs(kperp(j,i)*fik(j,i))**2
                !             ek=ek+energy(j,i)
                energy_b(j,i)=abs(kperp(j,i)*akpar(j,i,k))**2     
                energy_e(j,i)=abs(kperp(j,i)*fik(j,i,k))**2     
                energy_u(j,i)=-2./rhoi**2*(gama0(kperp(j,i)**2*rhoi**2/2.)-1.)*fik(j,i,k)*conjg(fik(j,i,k))     
                energy_ha(j,i)=-2.0*real((rhoe_LTe/(de*rhos))* nek(j,i,k)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k)))
                energy_ht(j,i)=2.0*(real(nek(j,i,k)*conjg(sqrt(2.0)*dzgk2(j,i,k)))+real(nek(j,i,k)*sqrt(2.0)*conjg(braakpargk2(j,i,k))))
                energy_it(j,i)=-2.0*(real(gk(j,i,k,gmin)/(sqrt(2.0))*conjg(cmplx(0.0,1.0)*ky(j)*fik(j,i,k))))
                energy_ia(j,i)=2.0*(real((sqrt(3.0)*rhos)/(2.0*de)*gk(j,i,k,gmin+1)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k))))
                energy_n(j,i) = abs(nek(j,i,k))**2
                energy_haa(j,i)=2.0*abs((rhoe_LTe/(de*rhos))* nek(j,i,k)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k)))
                energy_hta(j,i)=2.0*(abs(nek(j,i,k)*conjg(sqrt(2.0)*dzgk2(j,i,k)))+abs(nek(j,i,k)*sqrt(2.0)*conjg(braakpargk2(j,i,k))))
                energy_ita(j,i)=2.0*(abs(gk(j,i,k,gmin)/(sqrt(2.0))*conjg(cmplx(0.0,1.0)*ky(j)*fik(j,i,k))))
                energy_iaa(j,i)=2.0*(abs((sqrt(3.0)*rhos)/(2.0*de)*gk(j,i,k,gmin+1)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k))))
                energy_hk1(j,i)= 2.0*real((1+de**2*kperp(j,i)**2)*(fik(j,i,k)*conjg(akpar(j,i,k)))*(2*rhos**2/(rhoi**2))*(1.-gama0(kperp(j,i)**2*rhoi**2/2.)))
                energy_hk2(j,i)= 2.0*abs((1+de**2*kperp(j,i)**2)*(fik(j,i,k)*conjg(akpar(j,i,k)))*(2*rhos**2/(rhoi**2))*(1.-gama0(kperp(j,i)**2*rhoi**2/2.)))
                energy_ch(j,i) = 2.0*real(-kperp(j,i)**2*fik(j,i,k)*(1+de**2*kperp(j,i)**2)*conjg(akpar(j,i,k)))
                ek_b=ek_b+energy_b(j,i)
                ek_u=ek_u+energy_u(j,i)
                ek_e=ek_e+energy_e(j,i)
                hk_a=hk_a+energy_ha(j,i)
                hk_t=hk_t+energy_ht(j,i)
                ik_t=ik_t+energy_it(j,i)
                ik_a=ik_a+energy_ia(j,i)
                hk_aa=hk_aa+energy_haa(j,i)
                hk_ta=hk_ta+energy_hta(j,i)
                ik_ta=ik_ta+energy_ita(j,i)
                ik_aa=ik_aa+energy_iaa(j,i)
                hk_1=hk_1+energy_hk1(j,i)
                hk_2=hk_2+energy_hk2(j,i)
                chk=chk+energy_ch(j,i)

             end do
          end do
       else
          do i=1,NKx_par
             do j=1,NKy
                !             energy(j,i)=abs(kperp(j,i)*akpar(j,i))**2+abs(kperp(j,i)*fik(j,i))**2
                !             ek=ek+energy(j,i)
                energy_b(j,i)=abs(kperp(j,i)*akpar(j,i,k))**2     
                energy_e(j,i)=abs(kperp(j,i)*fik(j,i,k))**2
                energy_u(j,i)=-2./rhoi**2*(gama0(kperp(j,i)**2*rhoi**2/2.)-1.)*fik(j,i,k)*conjg(fik(j,i,k))
                energy_ha(j,i)=-2.0*real((rhoe_LTe/(de*rhos))* nek(j,i,k)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k)))
                energy_ht(j,i)=2.0*(real(nek(j,i,k)*conjg(sqrt(2.0)*dzgk2(j,i,k)))+real(nek(j,i,k)*sqrt(2.0)*conjg(braakpargk2(j,i,k))))
                energy_it(j,i)=-2.0*(real(gk(j,i,k,gmin)/(sqrt(2.0))*conjg(cmplx(0.0,1.0)*ky(j)*fik(j,i,k))))
                energy_ia(j,i)=2.0*(real((sqrt(3.0)*rhos)/(2.0*de)*gk(j,i,k,gmin+1)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k))))
                energy_n(j,i) = abs(nek(j,i,k))**2
                energy_haa(j,i)=2.0*abs((rhoe_LTe/(de*rhos))* nek(j,i,k)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k)))
                energy_hta(j,i)=2.0*(abs(nek(j,i,k)*conjg(sqrt(2.0)*dzgk2(j,i,k)))+abs(nek(j,i,k)*sqrt(2.0)*conjg(braakpargk2(j,i,k))))
                energy_ita(j,i)=2.0*(abs(gk(j,i,k,gmin)/(sqrt(2.0))*conjg(cmplx(0.0,1.0)*ky(j)*fik(j,i,k))))
                energy_iaa(j,i)=2.0*(abs((sqrt(3.0)*rhos)/(2.0*de)*gk(j,i,k,gmin+1)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k))))
                energy_hk1(j,i)= 2.0*real((1+de**2*kperp(j,i)**2)*(fik(j,i,k)*conjg(akpar(j,i,k)))*(2*rhos**2/(rhoi**2))*(1.-gama0(kperp(j,i)**2*rhoi**2/2.)))
                energy_hk2(j,i)= 2.0*abs((1+de**2*kperp(j,i)**2)*(fik(j,i,k)*conjg(akpar(j,i,k)))*(2*rhos**2/(rhoi**2))*(1.-gama0(kperp(j,i)**2*rhoi**2/2.)))
                energy_ch(j,i) = 2.0*real(-kperp(j,i)**2*fik(j,i,k)*(1+de**2*kperp(j,i)**2)*conjg(akpar(j,i,k)))
                ek_b=ek_b+energy_b(j,i)
                ek_u=ek_u+energy_u(j,i)
                ek_e=ek_e+energy_e(j,i)
                hk_a=hk_a+energy_ha(j,i)
                hk_t=hk_t+energy_ht(j,i)
                ik_t=ik_t+energy_it(j,i)
                ik_a=ik_a+energy_ia(j,i)
                ek_n=ek_n+energy_n(j,i)
                hk_aa=hk_aa+energy_haa(j,i)
                hk_ta=hk_ta+energy_hta(j,i)
                ik_ta=ik_ta+energy_ita(j,i)
                ik_aa=ik_aa+energy_iaa(j,i)
                hk_1=hk_1+energy_hk1(j,i)
                hk_2=hk_2+energy_hk2(j,i)
                chk=chk+energy_ch(j,i)
             end do
          end do
       end if
       ! ek_b=ek_b/(nlx*nly*1.0)
       ! ek_u=ek_u/(nlx*nly*1.0)
       ! ek_e=ek_e/(nlx*nly*1.0)
       ! hk_a=hk_a/(nlx*nly*1.0)
       ! hk_t=hk_t/(nlx*nly*1.0)
       ! ik_t=ik_t/(nlx*nly*1.0)
       ! ik_a=ik_a/(nlx*nly*1.0)
       ! ek_n=ek_n/(nlx*nly*1.0)
       ! hk_aa=hk_aa/(nlx*nly*1.0)
       ! hk_ta=hk_ta/(nlx*nly*1.0)
       ! ik_ta=ik_ta/(nlx*nly*1.0)
       ! ik_aa=ik_aa/(nlx*nly*1.0)
       ! hk_1=hk_1/(nlx*nly*1.0)
       ! hk_2=hk_2/(nlx*nly*1.0)
       ! chk=chk/(nlx*nly*1.0)
       ! energy_b=energy_b/(nlx*nly*1.0)
       ! energy_u=energy_u/(nlx*nly*1.0)
       ! energy_e=energy_e/(nlx*nly*1.0)
       ! energy_ha=energy_ha/(nlx*nly*1.0)
       ! energy_ht=energy_ht/(nlx*nly*1.0)
       ! energy_it=energy_it/(nlx*nly*1.0)
       ! energy_ia=energy_ia/(nlx*nly*1.0)
       ! energy_n=energy_n/(nlx*nly*1.0)
       ! energy_haa=energy_haa/(nlx*nly*1.0)
       ! energy_hta=energy_hta/(nlx*nly*1.0)
       ! energy_ita=energy_ita/(nlx*nly*1.0)
       ! energy_iaa=energy_iaa/(nlx*nly*1.0)
       ! energy_hk1=energy_hk1/(nlx*nly*1.0)
       ! energy_hk2=energy_hk2/(nlx*nly*1.0)
       ! energy_ch=energy_ch/(nlx*nly*1.0)
       call sum_reduce(ek_b,0)
       call sum_reduce(ek_u,0)
       call sum_reduce(ek_e,0)
       call sum_reduce(hk_a,0)
       call sum_reduce(hk_t,0)
       call sum_reduce(ik_t,0)
       call sum_reduce(ik_a,0)
       call sum_reduce(ek_n,0)
       call sum_reduce(hk_aa,0)
       call sum_reduce(hk_ta,0)
       call sum_reduce(ik_ta,0)
       call sum_reduce(ik_aa,0)
       call sum_reduce(hk_1,0)
       call sum_reduce(hk_2,0)
       call sum_reduce(chk,0)

       If(p .GT. 1  .AND. mod(p,save_energyfiles)==0)THEN
          
          work_u = 0.
          work_b = 0.
          work_e = 0.
          work_ha = 0.
          work_ht = 0.
          work_it = 0.
          work_ia = 0.
          work_n = 0.
          work_haa = 0.
          work_hta = 0.
          work_ita = 0.
          work_iaa = 0.
          work_hk1 = 0.
          work_hk2 = 0.
          work_ch = 0.
          m = 1+mod(iproc,NPE)*nky*nkx_par
          do i=1,nkx_par
             do j=1,nky
                work_b(m  )=energy_b(j,i)
                work_u(m  )=energy_u(j,i)
                work_e(m  )=energy_e(j,i)
                work_ha(m )=energy_ha(j,i)
                work_ht(m )=energy_ht(j,i)
                work_it(m )=energy_it(j,i)
                work_ia(m )=energy_ia(j,i)
                work_n(m ) =energy_n(j,i)
                work_haa(m )=energy_haa(j,i)
                work_hta(m )=energy_hta(j,i)
                work_ita(m )=energy_ita(j,i)
                work_iaa(m )=energy_iaa(j,i)
                work_hk1(m )=energy_hk1(j,i)
                work_hk2(m )=energy_hk2(j,i)
                work_ch(m )=energy_ch(j,i)
                m = m + 1
             end do
          end do
          call sum_reduce(work_b,0)
          call sum_reduce(work_u,0)
          call sum_reduce(work_e,0)
          call sum_reduce(work_ha,0)
          call sum_reduce(work_ht,0)
          call sum_reduce(work_it,0)
          call sum_reduce(work_ia,0)
          call sum_reduce(work_n,0)
          call sum_reduce(work_haa,0)
          call sum_reduce(work_hta,0)
          call sum_reduce(work_ita,0)
          call sum_reduce(work_iaa,0)
          call sum_reduce(hk_1,0)
          call sum_reduce(hk_2,0)
          call sum_reduce(chk,0)
          if(proc0) then
             m=1
             do i=1,nkx
                do j=1,nky
                   !          k=1+mod(k-1,nky*nkx_par)
                   sutot(j,i)=work_u(m)
                   sbtot(j,i)=work_b(m)
                   setot(j,i)=work_e(m)
                   shatot(j,i)=work_ha(m)
                   shttot(j,i)=work_ht(m)
                   sittot(j,i)=work_it(m)
                   siatot(j,i)=work_ia(m)
                   sntot(j,i)=work_n(m)
                   shaatot(j,i)=work_haa(m)
                   shtatot(j,i)=work_hta(m) 
                   sitatot(j,i)=work_ita(m)
                   siaatot(j,i)=work_iaa(m)
                   shk1tot(j,i)=work_hk1(m)
                   shk2tot(j,i)=work_hk2(m)
                   schtot(j,i)=work_ch(m)
                   m=m+1
                end do
             end do
             
          
             kshell_u=0.0
             kshell_b=0.0
             kshell_e=0.0
             kshell_ha=0.0
             kshell_ht=0.0
             kshell_it=0.0
             kshell_ia=0.0
             kshell_n=0.0
             kshell_haa=0.0
             kshell_hta=0.0
             kshell_ita=0.0
             kshell_iaa=0.0
             kshell_hk1=0.0
             kshell_hk2=0.0
             kshell_ch=0.0
             do kp=1,ceiling(kpmax)
                do i=1,nkx
                   do j=1,nky
                      if (sqrt(((i-1.)*1.0)**2+ky(j)**2).GE. kp-1 .AND. &
                           & sqrt(((i-1.)*1.0)**2+ky(j)**2).LT.kp+1) then
                         kshell_b(kp)=kshell_b(kp)+sbtot(j,i)
                         kshell_u(kp)=kshell_u(kp)+sutot(j,i)
                         kshell_e(kp)=kshell_e(kp)+setot(j,i)
                         kshell_ha(kp)=kshell_ha(kp)+shatot(j,i)
                         kshell_ht(kp)=kshell_ht(kp)+shttot(j,i)
                         kshell_it(kp)=kshell_it(kp)+sittot(j,i)
                         kshell_ia(kp)=kshell_ia(kp)+siatot(j,i)
                         kshell_n(kp)=kshell_n(kp)+sntot(j,i)
                         kshell_haa(kp)=kshell_haa(kp)+shaatot(j,i)
                         kshell_hta(kp)=kshell_hta(kp)+shtatot(j,i)
                         kshell_ita(kp)=kshell_ita(kp)+sitatot(j,i)
                         kshell_iaa(kp)=kshell_iaa(kp)+siaatot(j,i)
                         kshell_hk1(kp)=kshell_hk1(kp)+shk1tot(j,i)
                         kshell_hk2(kp)=kshell_hk2(kp)+shk2tot(j,i)
                         kshell_ch(kp)=kshell_ch(kp)+schtot(j,i)
                      end if
                   end do
                end do
                !write(unitnumber, '(2i5,2g16.8)') kglobal(k), kp ,kshell_u(kp), kshell_b(kp)
				!AVK: we were missing a factor of nlz
                ! write(unitnumber, '(2i5,15g16.8)') kglobal(k), kp ,kshell_u(kp)/(nlz*1.), &
		!		&kshell_b(kp)/(nlz*1.), kshell_e(kp)/(nlz*1.), kshell_ht(kp)/(nlz*1.), kshell_ha(kp)/(nlz*1.), kshell_it(kp)/(nlz*1.0), kshell_ia(kp)/(nlz*1.0), kshell_hta(kp)/(nlz*1.), kshell_haa(kp)/(nlz*1.), kshell_ita(kp)/(nlz*1.0), kshell_iaa(kp)/(nlz*1.0), kshell_n(kp)/(nlz*1.0), kshell_hk1(kp)/(nlz*1.0), kshell_hk2(kp)/(nlz*1.0), kshell_ch(kp)/(nlz*1.0)
                write(unitnumber, '(2i5,18g16.8)') kglobal(k), kp, kshell_u(kp), &
				&kshell_b(kp), kshell_e(kp), kshell_ht(kp), kshell_ha(kp), kshell_it(kp), kshell_ia(kp), kshell_hta(kp), kshell_haa(kp), kshell_ita(kp), kshell_iaa(kp), kshell_n(kp)
             end do
             write(unitnumber, *) '              '   !useful?
          end if
       end if
    end do
    if (proc0 .AND. p .GT. 1  .AND. mod(p,save_energyfiles)==0) then
       close(unit=unitnumber)
       deallocate(kshell_u)
       deallocate(kshell_b)
       deallocate(kshell_e)
       deallocate(kshell_ha)
       deallocate(kshell_ht)
       deallocate(kshell_it)
       deallocate(kshell_ia)
       deallocate(kshell_n)
       deallocate(kshell_haa)
       deallocate(kshell_hta)
       deallocate(kshell_ita)
       deallocate(kshell_iaa)
       deallocate(kshell_ch)
    end if

  end subroutine k_energy
  
!**********************************************
! DEPRECATED 2/28/23 A. Velberg perturbed energy output 5/20/22. All this does is create a new file called "_energykp_" using the same code, but called with the perturbed quantities as the input... 
! To inlcude again, uncomment call in REGK.F90, add k_energy_pert from diag in REGK, and add k_energy_pert as public on line 26 here
!**********************************************
subroutine k_energy_pert(unitnumber,time,p,akpar,fik,nek,gk,dzgk2,braakpargk2,ek_u,ek_b,ek_e,hk_a,hk_t,ik_t,ik_a,ek_n,hk_aa,hk_ta,ik_ta,ik_aa,hk_1,hk_2,chk,runname)
    !energy per unit volume
    !helicity and energy injection - LMM 
    use constants
    use mp
    use grid
    implicit none
    
    COMPLEX, DIMENSION(NKy,NKX_PAR,nlz_par)::akpar,fik,nek,dzgk2,braakpargk2
    complex, dimension(nky,nkx_par,nlz_par,gmin:ngtot)::gk
    complex,dimension(nky,nkx)::sutot,sbtot,setot,shatot,shttot,sittot,siatot,sntot,shaatot,shtatot,sitatot,siaatot,shk1tot, shk2tot, schtot
    complex, dimension (nky*nkx_par*NPE) :: work_u,work_b, work_e, work_ha, work_ht, work_it, work_ia, work_n, work_haa, work_hta, work_ita, work_iaa, work_hk1, work_hk2, work_ch
    real, dimension(nky,nkx_par)::energy_u,energy_b, energy_e,energy_ha, energy_ht, energy_it, energy_ia, energy_n, energy_haa, energy_hta, energy_ita, energy_iaa, energy_hk1, energy_hk2, energy_ch ! LUIS, 28/4/15: kinetic, magnetic, electric energy
    real, dimension(:),allocatable::kshell_u,kshell_b,kshell_e,kshell_ha, kshell_ht, kshell_it, kshell_ia, kshell_n,kshell_haa, kshell_hta, kshell_ita, kshell_iaa, kshell_hk1, kshell_hk2, kshell_ch
    real:: ek_u,ek_b,ek_e,hk_a,hk_t,ik_t,ik_a,ek_n,hk_aa,hk_ta,ik_ta,ik_aa,hk_1,hk_2,time,kpmax,chk!,de,l
    integer::i,j,k,p,unitnumber,kp,m
    character(len=100)::filename1
    character(len=100),intent(in)::runname

    if (proc0 .AND. p .GT. 1  .AND. mod(p,save_energyfiles)==0) then
       if (time.lt.10) then
          write(filename1,'(a8,f5.3,a4)') energykp,time,".dat"
       end if
       if (time .ge. 10 .AND. time .lt. 100) then
          write(filename1,'(a8,f6.3,a4)') energykp,time,".dat"
       end if
       if (time .ge. 100 .AND. time .lt. 1000) then
          write(filename1,'(a8,f7.3,a4)') energykp,time,".dat"
       end if
       if (time .ge. 1000 .AND. time .lt. 10000) then
          write(filename1,'(a8,f8.3,a4)') energykp,time,".dat"
       end if
       if (time .ge. 10000 .AND. time .lt. 100000) then
          write(filename1,'(a8,f9.3,a4)') energykp,time,".dat"
       end if

       filename1 = trim(runname)//"_"//trim(filename1)
!NFL: 02/07/13 the following line is wrong!       
!       kpmax=sqrt(ky(nky)**2.+((nkx-1.)*1.0)**2.)
!       kpmax=sqrt((nky/2.*lx/ly)**2.+((nkx-1.)*1.0)**2.)
       kpmax=sqrt(ky(nky/2+1)**2.+((nkx-1.)*1.0)**2.)
       allocate(kshell_u(ceiling(kpmax)))
       allocate(kshell_b(ceiling(kpmax)))
       allocate(kshell_e(ceiling(kpmax)))
       allocate(kshell_ha(ceiling(kpmax)))
       allocate(kshell_ht(ceiling(kpmax)))
       allocate(kshell_it(ceiling(kpmax)))   
       allocate(kshell_ia(ceiling(kpmax)))    
       allocate(kshell_n(ceiling(kpmax)))
       allocate(kshell_haa(ceiling(kpmax)))
       allocate(kshell_hta(ceiling(kpmax)))
       allocate(kshell_ita(ceiling(kpmax)))   
       allocate(kshell_iaa(ceiling(kpmax)))     
       allocate(kshell_hk1(ceiling(kpmax)))   
       allocate(kshell_hk2(ceiling(kpmax)))   
       allocate(kshell_ch(ceiling(kpmax)))
       open (unit=unitnumber, file= trim(filename1))
       write(unitnumber, '(20a15)') 'kglobal(k)', 'kp' ,'kinetic', 'magnetic', 'electric', 'helicity_t', 'helicity_a', 'injection_t', 'injection_a', 'helicity_ta', 'helicity_aa', 'injection_ta', 'injection_aa', 'density', 'helicity_real', 'helicity_abs','cross_helicity'
     end if
    
    do k=1,nlz_par
       ek_u=0.0
       ek_e=0.0
       ek_b=0.0
       hk_a=0.0
       hk_t=0.0
       ik_t=0.0
       ik_a=0.0
       ek_n=0.0
       hk_aa=0.0
       hk_ta=0.0
       ik_ta=0.0
       ik_aa=0.0
       hk_1=0.0
       hk_2=0.0
       chk=0.0
       !the kx=0 modes only get added once
       !the kx diff 0 modes get added twice because of the reality condition
       if(mod(iproc,npe)==0) then
          do j=1, nky
             !          energy(j,1)=0.5*abs(kperp(j,1)*akpar(j,1))**2+abs(kperp(j,1)*fik(j,1))**2      
             energy_b(j,1)=0.5*abs(kperp(j,1)*akpar(j,1,k))**2     
             energy_e(j,1)=0.5*abs(kperp(j,1)*fik(j,1,k))**2     
             energy_u(j,1)=-1.0/rhoi**2*(gama0(kperp(j,1)**2*rhoi**2/2.)-1.)*fik(j,1,k)*conjg(fik(j,1,k))     
             energy_ha(j,1)= -real((rhoe_LTe/(de*rhos))* nek(j,1,k)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,1,k)))
             energy_ht(j,1)= real(nek(j,1,k)*conjg(sqrt(2.0)*dzgk2(j,1,k)))+real(nek(j,1,k)*sqrt(2.0)*conjg(braakpargk2(j,1,k)))
             energy_it(j,1) = -real(gk(j,1,k,gmin)/(sqrt(2.0))*conjg(cmplx(0.0,1.0)*ky(j)*fik(j,1,k)))
             energy_ia(j,1) = real((sqrt(3.0)*rhos)/(2.0*de)*gk(j,1,k,gmin+1)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,1,k)) )
             energy_n(j,1) = 0.5*abs(nek(j,1,k))**2
             energy_haa(j,1)= abs((rhoe_LTe/(de*rhos))* nek(j,1,k)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,1,k)))
             energy_hta(j,1)= abs(nek(j,1,k)*conjg(sqrt(2.0)*dzgk2(j,1,k)))+abs(nek(j,1,k)*sqrt(2.0)*conjg(braakpargk2(j,1,k)))
             energy_ita(j,1) = abs(gk(j,1,k,gmin)/(sqrt(2.0))*conjg(cmplx(0.0,1.0)*ky(j)*fik(j,1,k)))
             energy_iaa(j,1) = abs((sqrt(3.0)*rhos)/(2.0*de)*gk(j,1,k,gmin+1)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,1,k)))
             energy_hk1(j,1)= real((1+de**2*kperp(j,1)**2)*(fik(j,1,k)*conjg(akpar(j,1,k)))*(2*rhos**2/(rhoi**2))*(1.-gama0(kperp(j,1)**2*rhoi**2/2.)))
             energy_hk2(j,1)= abs((1+de**2*kperp(j,1)**2)*(fik(j,1,k)*conjg(akpar(j,1,k)))*(2*rhos**2/(rhoi**2))*(1.-gama0(kperp(j,1)**2*rhoi**2/2.)))    
             energy_ch(j,1) = real(-kperp(j,1)**2*fik(j,1,k)*(1+de**2*kperp(j,1)**2)*conjg(akpar(j,1,k)))
             ek_b=ek_b+energy_b(j,1)
             ek_u=ek_u+energy_u(j,1)
             ek_e=ek_e+energy_e(j,1)
             hk_a=hk_a+energy_ha(j,1)
             hk_t=hk_t+energy_ht(j,1)
             ik_t=ik_t+energy_it(j,1)
             ik_a=ik_a+energy_ia(j,1)
             ek_n=ek_n+energy_n(j,1)
             hk_aa=hk_aa+energy_haa(j,1)
             hk_ta=hk_ta+energy_hta(j,1)
             ik_ta=ik_ta+energy_ita(j,1)
             ik_aa=ik_aa+energy_iaa(j,1)
             hk_1=hk_1+energy_hk1(j,1)
             hk_2=hk_2+energy_hk2(j,1)
             chk=chk+energy_ch(j,1)
          end do
          do i=2,NKx_par
             do j=1,NKy
                !             energy(j,i)=abs(kperp(j,i)*akpar(j,i))**2+abs(kperp(j,i)*fik(j,i))**2
                !             ek=ek+energy(j,i)
                energy_b(j,i)=abs(kperp(j,i)*akpar(j,i,k))**2     
                energy_e(j,i)=abs(kperp(j,i)*fik(j,i,k))**2     
                energy_u(j,i)=-2./rhoi**2*(gama0(kperp(j,i)**2*rhoi**2/2.)-1.)*fik(j,i,k)*conjg(fik(j,i,k))     
                energy_ha(j,i)=-2.0*real((rhoe_LTe/(de*rhos))* nek(j,i,k)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k)))
                energy_ht(j,i)=2.0*(real(nek(j,i,k)*conjg(sqrt(2.0)*dzgk2(j,i,k)))+real(nek(j,i,k)*sqrt(2.0)*conjg(braakpargk2(j,i,k))))
                energy_it(j,i)=-2.0*(real(gk(j,i,k,gmin)/(sqrt(2.0))*conjg(cmplx(0.0,1.0)*ky(j)*fik(j,i,k))))
                energy_ia(j,i)=2.0*(real((sqrt(3.0)*rhos)/(2.0*de)*gk(j,i,k,gmin+1)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k))))
                energy_n(j,i) = abs(nek(j,i,k))**2
                energy_haa(j,i)=2.0*abs((rhoe_LTe/(de*rhos))* nek(j,i,k)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k)))
                energy_hta(j,i)=2.0*(abs(nek(j,i,k)*conjg(sqrt(2.0)*dzgk2(j,i,k)))+abs(nek(j,i,k)*sqrt(2.0)*conjg(braakpargk2(j,i,k))))
                energy_ita(j,i)=2.0*(abs(gk(j,i,k,gmin)/(sqrt(2.0))*conjg(cmplx(0.0,1.0)*ky(j)*fik(j,i,k))))
                energy_iaa(j,i)=2.0*(abs((sqrt(3.0)*rhos)/(2.0*de)*gk(j,i,k,gmin+1)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k))))
                energy_hk1(j,i)= 2.0*real((1+de**2*kperp(j,i)**2)*(fik(j,i,k)*conjg(akpar(j,i,k)))*(2*rhos**2/(rhoi**2))*(1.-gama0(kperp(j,i)**2*rhoi**2/2.)))
                energy_hk2(j,i)= 2.0*abs((1+de**2*kperp(j,i)**2)*(fik(j,i,k)*conjg(akpar(j,i,k)))*(2*rhos**2/(rhoi**2))*(1.-gama0(kperp(j,i)**2*rhoi**2/2.)))
                energy_ch(j,i) = 2.0*real(-kperp(j,i)**2*fik(j,i,k)*(1+de**2*kperp(j,i)**2)*conjg(akpar(j,i,k)))
                ek_b=ek_b+energy_b(j,i)
                ek_u=ek_u+energy_u(j,i)
                ek_e=ek_e+energy_e(j,i)
                hk_a=hk_a+energy_ha(j,i)
                hk_t=hk_t+energy_ht(j,i)
                ik_t=ik_t+energy_it(j,i)
                ik_a=ik_a+energy_ia(j,i)
                hk_aa=hk_aa+energy_haa(j,i)
                hk_ta=hk_ta+energy_hta(j,i)
                ik_ta=ik_ta+energy_ita(j,i)
                ik_aa=ik_aa+energy_iaa(j,i)
                hk_1=hk_1+energy_hk1(j,i)
                hk_2=hk_2+energy_hk2(j,i)
                chk=chk+energy_ch(j,i)

             end do
          end do
       else
          do i=1,NKx_par
             do j=1,NKy
                !             energy(j,i)=abs(kperp(j,i)*akpar(j,i))**2+abs(kperp(j,i)*fik(j,i))**2
                !             ek=ek+energy(j,i)
                energy_b(j,i)=abs(kperp(j,i)*akpar(j,i,k))**2     
                energy_e(j,i)=abs(kperp(j,i)*fik(j,i,k))**2
                energy_u(j,i)=-2./rhoi**2*(gama0(kperp(j,i)**2*rhoi**2/2.)-1.)*fik(j,i,k)*conjg(fik(j,i,k))
                energy_ha(j,i)=-2.0*real((rhoe_LTe/(de*rhos))* nek(j,i,k)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k)))
                energy_ht(j,i)=2.0*(real(nek(j,i,k)*conjg(sqrt(2.0)*dzgk2(j,i,k)))+real(nek(j,i,k)*sqrt(2.0)*conjg(braakpargk2(j,i,k))))
                energy_it(j,i)=-2.0*(real(gk(j,i,k,gmin)/(sqrt(2.0))*conjg(cmplx(0.0,1.0)*ky(j)*fik(j,i,k))))
                energy_ia(j,i)=2.0*(real((sqrt(3.0)*rhos)/(2.0*de)*gk(j,i,k,gmin+1)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k))))
                energy_n(j,i) = abs(nek(j,i,k))**2
                energy_haa(j,i)=2.0*abs((rhoe_LTe/(de*rhos))* nek(j,i,k)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k)))
                energy_hta(j,i)=2.0*(abs(nek(j,i,k)*conjg(sqrt(2.0)*dzgk2(j,i,k)))+abs(nek(j,i,k)*sqrt(2.0)*conjg(braakpargk2(j,i,k))))
                energy_ita(j,i)=2.0*(abs(gk(j,i,k,gmin)/(sqrt(2.0))*conjg(cmplx(0.0,1.0)*ky(j)*fik(j,i,k))))
                energy_iaa(j,i)=2.0*(abs((sqrt(3.0)*rhos)/(2.0*de)*gk(j,i,k,gmin+1)*conjg(cmplx(0.0,1.0)*ky(j)*akpar(j,i,k))))
                energy_hk1(j,i)= 2.0*real((1+de**2*kperp(j,i)**2)*(fik(j,i,k)*conjg(akpar(j,i,k)))*(2*rhos**2/(rhoi**2))*(1.-gama0(kperp(j,i)**2*rhoi**2/2.)))
                energy_hk2(j,i)= 2.0*abs((1+de**2*kperp(j,i)**2)*(fik(j,i,k)*conjg(akpar(j,i,k)))*(2*rhos**2/(rhoi**2))*(1.-gama0(kperp(j,i)**2*rhoi**2/2.)))
                energy_ch(j,i) = 2.0*real(-kperp(j,i)**2*fik(j,i,k)*(1+de**2*kperp(j,i)**2)*conjg(akpar(j,i,k)))
                ek_b=ek_b+energy_b(j,i)
                ek_u=ek_u+energy_u(j,i)
                ek_e=ek_e+energy_e(j,i)
                hk_a=hk_a+energy_ha(j,i)
                hk_t=hk_t+energy_ht(j,i)
                ik_t=ik_t+energy_it(j,i)
                ik_a=ik_a+energy_ia(j,i)
                ek_n=ek_n+energy_n(j,i)
                hk_aa=hk_aa+energy_haa(j,i)
                hk_ta=hk_ta+energy_hta(j,i)
                ik_ta=ik_ta+energy_ita(j,i)
                ik_aa=ik_aa+energy_iaa(j,i)
                hk_1=hk_1+energy_hk1(j,i)
                hk_2=hk_2+energy_hk2(j,i)
                chk=chk+energy_ch(j,i)
             end do
          end do
       end if
       
       call sum_reduce(ek_b,0)
       call sum_reduce(ek_u,0)
       call sum_reduce(ek_e,0)
       call sum_reduce(hk_a,0)
       call sum_reduce(hk_t,0)
       call sum_reduce(ik_t,0)
       call sum_reduce(ik_a,0)
       call sum_reduce(ek_n,0)
       call sum_reduce(hk_aa,0)
       call sum_reduce(hk_ta,0)
       call sum_reduce(ik_ta,0)
       call sum_reduce(ik_aa,0)
       call sum_reduce(hk_1,0)
       call sum_reduce(hk_2,0)
       call sum_reduce(chk,0)

       If(p .GT. 1  .AND. mod(p,save_energyfiles)==0)THEN
          
          work_u = 0.
          work_b = 0.
          work_e = 0.
          work_ha = 0.
          work_ht = 0.
          work_it = 0.
          work_ia = 0.
          work_n = 0.
          work_haa = 0.
          work_hta = 0.
          work_ita = 0.
          work_iaa = 0.
          work_hk1 = 0.
          work_hk2 = 0.
          work_ch = 0.
          m = 1+mod(iproc,NPE)*nky*nkx_par
          do i=1,nkx_par
             do j=1,nky
                work_b(m  )=energy_b(j,i)
                work_u(m  )=energy_u(j,i)
                work_e(m  )=energy_e(j,i)
                work_ha(m )=energy_ha(j,i)
                work_ht(m )=energy_ht(j,i)
                work_it(m )=energy_it(j,i)
                work_ia(m )=energy_ia(j,i)
                work_n(m ) =energy_n(j,i)
                work_haa(m )=energy_haa(j,i)
                work_hta(m )=energy_hta(j,i)
                work_ita(m )=energy_ita(j,i)
                work_iaa(m )=energy_iaa(j,i)
                work_hk1(m )=energy_hk1(j,i)
                work_hk2(m )=energy_hk2(j,i)
                work_ch(m )=energy_ch(j,i)
                m = m + 1
             end do
          end do
          call sum_reduce(work_b,0)
          call sum_reduce(work_u,0)
          call sum_reduce(work_e,0)
          call sum_reduce(work_ha,0)
          call sum_reduce(work_ht,0)
          call sum_reduce(work_it,0)
          call sum_reduce(work_ia,0)
          call sum_reduce(work_n,0)
          call sum_reduce(work_haa,0)
          call sum_reduce(work_hta,0)
          call sum_reduce(work_ita,0)
          call sum_reduce(work_iaa,0)
          call sum_reduce(hk_1,0)
          call sum_reduce(hk_2,0)
          call sum_reduce(chk,0)
          if(proc0) then
             m=1
             do i=1,nkx
                do j=1,nky
                   !          k=1+mod(k-1,nky*nkx_par)
                   sutot(j,i)=work_u(m)
                   sbtot(j,i)=work_b(m)
                   setot(j,i)=work_e(m)
                   shatot(j,i)=work_ha(m)
                   shttot(j,i)=work_ht(m)
                   sittot(j,i)=work_it(m)
                   siatot(j,i)=work_ia(m)
                   sntot(j,i)=work_n(m)
                   shaatot(j,i)=work_haa(m)
                   shtatot(j,i)=work_hta(m) 
                   sitatot(j,i)=work_ita(m)
                   siaatot(j,i)=work_iaa(m)
                   shk1tot(j,i)=work_hk1(m)
                   shk2tot(j,i)=work_hk2(m)
                   schtot(j,i)=work_ch(m)
                   m=m+1
                end do
             end do
             
          
             kshell_u=0.0
             kshell_b=0.0
             kshell_e=0.0
             kshell_ha=0.0
             kshell_ht=0.0
             kshell_it=0.0
             kshell_ia=0.0
             kshell_n=0.0
             kshell_haa=0.0
             kshell_hta=0.0
             kshell_ita=0.0
             kshell_iaa=0.0
             kshell_hk1=0.0
             kshell_hk2=0.0
             kshell_ch=0.0
             do kp=1,ceiling(kpmax)
                do i=1,nkx
                   do j=1,nky
                      if (sqrt(((i-1.)*1.0)**2+ky(j)**2).GE. kp-1 .AND. &
                           & sqrt(((i-1.)*1.0)**2+ky(j)**2).LT.kp+1) then
                         kshell_b(kp)=kshell_b(kp)+sbtot(j,i)
                         kshell_u(kp)=kshell_u(kp)+sutot(j,i)
                         kshell_e(kp)=kshell_e(kp)+setot(j,i)
                         kshell_ha(kp)=kshell_ha(kp)+shatot(j,i)
                         kshell_ht(kp)=kshell_ht(kp)+shttot(j,i)
                         kshell_it(kp)=kshell_it(kp)+sittot(j,i)
                         kshell_ia(kp)=kshell_ia(kp)+siatot(j,i)
                         kshell_n(kp)=kshell_n(kp)+sntot(j,i)
                         kshell_haa(kp)=kshell_haa(kp)+shaatot(j,i)
                         kshell_hta(kp)=kshell_hta(kp)+shtatot(j,i)
                         kshell_ita(kp)=kshell_ita(kp)+sitatot(j,i)
                         kshell_iaa(kp)=kshell_iaa(kp)+siaatot(j,i)
                         kshell_hk1(kp)=kshell_hk1(kp)+shk1tot(j,i)
                         kshell_hk2(kp)=kshell_hk2(kp)+shk2tot(j,i)
                         kshell_ch(kp)=kshell_ch(kp)+schtot(j,i)
                      end if
                   end do
                end do
                !write(unitnumber, '(2i5,2g16.8)') kglobal(k), kp ,kshell_u(kp), kshell_b(kp)
				!AVK: we were missing a factor of nlz
                ! write(unitnumber, '(2i5,15g16.8)') kglobal(k), kp ,kshell_u(kp)/(nlz*1.), &
		!		&kshell_b(kp)/(nlz*1.), kshell_e(kp)/(nlz*1.), kshell_ht(kp)/(nlz*1.), kshell_ha(kp)/(nlz*1.), kshell_it(kp)/(nlz*1.0), kshell_ia(kp)/(nlz*1.0), kshell_hta(kp)/(nlz*1.), kshell_haa(kp)/(nlz*1.), kshell_ita(kp)/(nlz*1.0), kshell_iaa(kp)/(nlz*1.0), kshell_n(kp)/(nlz*1.0), kshell_hk1(kp)/(nlz*1.0), kshell_hk2(kp)/(nlz*1.0), kshell_ch(kp)/(nlz*1.0)
                write(unitnumber, '(2i5,18g16.8)') kglobal(k), kp, kshell_u(kp), &
				&kshell_b(kp), kshell_e(kp), kshell_ht(kp), kshell_ha(kp), kshell_it(kp), kshell_ia(kp), kshell_hta(kp), kshell_haa(kp), kshell_ita(kp), kshell_iaa(kp), kshell_n(kp)
             end do
             write(unitnumber, *) '              '   !useful?
          end if
       end if
    end do
    if (proc0 .AND. p .GT. 1  .AND. mod(p,save_energyfiles)==0) then
       close(unit=unitnumber)
       deallocate(kshell_u)
       deallocate(kshell_b)
       deallocate(kshell_e)
       deallocate(kshell_ha)
       deallocate(kshell_ht)
       deallocate(kshell_it)
       deallocate(kshell_ia)
       deallocate(kshell_n)
       deallocate(kshell_haa)
       deallocate(kshell_hta)
       deallocate(kshell_ita)
       deallocate(kshell_iaa)
       deallocate(kshell_ch)
    end if

  end subroutine k_energy_pert


!*********************************************
!  subroutine k_energy(akpar,fik,ek)
!    !energy per unit volume
!    use constants
!    use mp
!    use grid
!    implicit none!
!
!    COMPLEX, DIMENSION(NKy,NKX_PAR)::akpar,fik    
!    real:: ek
!    integer::i,j
!
!    ek=0.0
!    !the kx=0 modes only get added once
!    !the kx diff 0 modes get added twice because of the reality condition
!    if(proc0) then
!       do j=1, nky
!          ek=ek+0.5*abs(kperp(j,1)*akpar(j,1))**2+abs(kperp(j,1)*fik(j,1))**2
!       end do
!       do i=2,NKx_par
!          do j=1,NKy
!             ek=ek+abs(kperp(j,i)*akpar(j,i))**2+abs(kperp(j,i)*fik(j,i))**2
!          end do
!       end do
!    else
!       do i=1,NKx_par
!          do j=1,NKy
!             ek=ek+abs(kperp(j,i)*akpar(j,i))**2+abs(kperp(j,i)*fik(j,i))**2
!          end do
!       end do
!    end if
!    ek=ek/(nlx*nly*1.0)
!    call sum_reduce(ek,0)
!
!    end subroutine k_energy
!*********************************************


!*********************************************
  subroutine tot_diss(uepar,uepar_eq,fi,ne,lost)
    !dissipated energy per unit volume
    use constants
    use mp
    use grid
    implicit none
    
    real, DIMENSION(NLx,NLy_PAR):: uepar,uepar_eq,fi,ne
    real::lost
    integer::i,j

    lost=0.0

    do j=1,nly_par
       do i=1,nlx
          lost=lost+res*uepar(i,j)*(uepar(i,j)-uepar_eq(i,j))+niu*fi(i,j)*ne(i,j)
       end do
    end do
    
    call sum_reduce(lost,0)
    if (proc0) then 
       lost=lost/(nlx*nly*1.0)
    end if

  end subroutine tot_diss


!*********************************************
   !subroutine injected(dummy_real,fi,inj)! Z.Liu 6/29/2020
    subroutine injected(dummy_real,fi,ne,inj)
    use constants
    use mp
    use grid
    implicit none
    
    !real, DIMENSION(NLx,NLy_PAR,NLz_PAR):: dummy_real,fi
    real, DIMENSION(NLx,NLy_PAR,NLz_PAR):: dummy_real,fi,ne
    real::inj
    integer::i,j,k

    inj=0.0
    do k=1,nlz_par
       do j=1,nly_par
          do i=1,nlx
             !inj=inj + dummy_real(i,j,k)*fi(i,j,k)
             !inj=inj + rhos_diag**2*ne(i,j,k)*dummy_real(i,j,k)
            inj = inj + (rhos_diag**2*ne(i,j,k) - fi(i,j,k))*dummy_real(i,j,k) !Z.Liu 7/3/2020
            !inj = inj + (rhos_diag**2*(gama0(kperp(j,1)**2*rhoi**2/2.)-1.)/rhoi**2 -1)*fi(i,j,k)*dummy_real(i,j,k)
            !inj = inj + (rhos_diag**2 - 1./(2./rhoi**2*(gama0(kperp(j,i)**2*rhoi**2/2.)-1.)))*ne(i,j,k)*dummy_real(i,j,k)
          end do
       end do
    end do
    
    call sum_reduce(inj,0)
    if (proc0) then 
       inj=Lx/(1.0*nlx)*ly/(1.0*nly)*lz/(1.0*nlz)*inj
    end if

  end subroutine injected


!*********************************************
  subroutine zsteploss(FIk,nek,uekpar,Akpar,gk,dt,proc_X,energyzloss)
    !original subroutine by AM to quantify energy loss in z-step
    !NFL, 15/08/2012 -- changed to allow for no g's and to be compatible with anjor=.true.
	!AVK 26/03/13: Changed rhos to rhos_diag. It is set in set_anjor. 
	! rhos_diag = rhos if anjor=false, otherwise rhos_diag = 1
    ! Luis, Nov 2014, fixed a small bug: b_energy, de_energy, ne_energy, phine_energy,sumgsquaretot 
    ! were being called with size "nlz", but only "nlz_par" elements were then used; 
    !also improved z_integration step (new subroutine) 
    use mp
    USE constants
    use grid
!    use transforms
    implicit none
    
    COMPLEX, DIMENSION(NKy,NKX_PAR,nlz_par):: FIK, nek, uekpar, Akpar
    !complex, dimension(nlz/2+1):: bkz_energy
    COMPLEX, DIMENSION(NKy,NKX_PAR,nlz_par,gmin:ngtot)::gk 
!    COMPLEX, DIMENSION(:,:,:,:)::gk 
    real,dimension(nlz_par)::b_energy, de_energy, ne_energy, phine_energy
    real::b_energytot, de_energytot, ne_energytot, phine_energytot
    real,dimension(nlz_par)::sumgsquare
    real::sumgsquaretot
    real::energyz,dt,energyzloss
    real,save::energyz_old
    logical,save::first=.true.
    integer::i,j,k, proc_X,ng

    if(first) then
       energyz_old=0.0
       first=.false.
    end if

    b_energy=0.0
    b_energytot=0.0
    de_energy=0.0
    de_energytot=0.0
    ne_energy=0.0
    ne_energytot=0.0
    phine_energy=0.0
    phine_energytot=0.0
    sumgsquare=0.0
    sumgsquaretot=0.0

    if (mod(iproc,npe)==0) then
       do k=1,nlz_par
          do j=1, nky
             if(g_inc) then
                do ng=gmin,ngtot
                   sumgsquare(k)=sumgsquare(k)+0.5*rhos_diag**2*gk(j,1,k,ng)*conjg(gk(j,1,k,ng))
                end do
             else
                sumgsquare=0.0
             end if
             b_energy(k)=b_energy(k)+0.5*kperp(j,1)**2*Akpar(j,1,k)*conjg(akpar(j,1,k))
             de_energy(k)=de_energy(k)+0.5*de**2*uekpar(j,1,k)*conjg(uekpar(j,1,k))
             ne_energy(k)=ne_energy(k)+0.5*rhos_diag**2*nek(j,1,k)*conjg(nek(j,1,k))
             phine_energy(k)=phine_energy(k)-1.0/rhoi**2*(gama0(kperp(j,1)**2*rhoi**2/2.)-1.)*&
                  FIK(j,1,k)*conjg(FIK(j,1,k))
             
          end do
          do i=2, nkx_par
             do j=1, nky
                if(g_inc) then
                   do ng=gmin,ngtot
                      sumgsquare(k)=sumgsquare(k)+rhos_diag**2*gk(j,i,k,ng)*conjg(gk(j,i,k,ng))
                   end do
                else
                   sumgsquare=0.0
                end if
                b_energy(k)=b_energy(k)+kperp(j,i)**2*Akpar(j,i,k)*conjg(akpar(j,i,k))
                de_energy(k)=de_energy(k)+de**2*uekpar(j,i,k)*conjg(uekpar(j,i,k))
                ne_energy(k)=ne_energy(k)+rhos_diag**2*nek(j,i,k)*conjg(nek(j,i,k))
                phine_energy(k)=phine_energy(k)-2./rhoi**2*(gama0(kperp(j,i)**2*rhoi**2/2.)-1.)*&
                     FIK(j,i,k)*conjg(FIK(j,i,k))
                
             end do
          end do
       end do
    else
       do k=1,nlz_par
          do i=1,NKx_par
             do j=1,NKy
                if(g_inc) then
                   do ng=gmin,ngtot
                      sumgsquare(k)=sumgsquare(k)+rhos_diag**2*gk(j,i,k,ng)*conjg(gk(j,i,k,ng))
                   end do
                else
                   sumgsquare=0.0
                end if
                b_energy(k)=b_energy(k)+kperp(j,i)**2*Akpar(j,i,k)*conjg(akpar(j,i,k))
                de_energy(k)=de_energy(k)+de**2*uekpar(j,i,k)*conjg(uekpar(j,i,k))
                ne_energy(k)=ne_energy(k)+rhos_diag**2*nek(j,i,k)*conjg(nek(j,i,k))
                phine_energy(k)=phine_energy(k)-2./rhoi**2*(gama0(kperp(j,i)**2*rhoi**2/2.)-1.)*&
                     FIK(j,i,k)*conjg(FIK(j,i,k))
                
             end do
          end do
       end do
    end if
    b_energy(:)=b_energy(:)/(nlx*nly*1.0)
    de_energy(:)=de_energy(:)/(nlx*nly*1.0)
    ne_energy(:)=ne_energy(:)/(nlx*nly*1.0)
    phine_energy(:)=phine_energy(:)/(nlx*nly*1.0)
    sumgsquare=sumgsquare/(nlx*nly*1.0)

    if(three_D) then
       
          call z_integration(b_energy,b_energytot,dz,nlz_par) 
          call z_integration(de_energy,de_energytot,dz,nlz_par)
          call z_integration(ne_energy,ne_energytot,dz,nlz_par)
          call z_integration(phine_energy,phine_energytot,dz,nlz_par)
          call z_integration(sumgsquare,sumgsquaretot,dz,nlz_par)
           
    else
       b_energytot = b_energy(1)
       de_energytot = de_energy(1)
       ne_energytot = ne_energy(1)
       phine_energytot = phine_energy(1)
       sumgsquaretot = sumgsquare(1)
    end if

!NFL: commented the following lines, 18/03/13    
!    do k=1, nlz_par
!       b_energytot = b_energytot + b_energy(k)
!       de_energytot = de_energytot + de_energy(k)
!       ne_energytot = ne_energytot + ne_energy(k)
!       phine_energytot = phine_energytot + phine_energy(k)
!       sumgsquaretot = sumgsquaretot + sumgsquare(k)
!    end do
!    b_energytot=b_energytot/(nlz*1.0)
!    de_energytot=de_energytot/(nlz*1.0)
!    ne_energytot=ne_energytot/(nlz*1.0)
!    phine_energytot=phine_energytot/(nlz*1.0)
!    sumgsquaretot=sumgsquaretot/(nlz*1.0)

    call sum_reduce (b_energytot, proc_X)
    call sum_reduce (de_energytot, proc_X)
    call sum_reduce (ne_energytot, proc_X)
    call sum_reduce (phine_energytot, proc_X)
    call sum_reduce (sumgsquaretot, proc_X)
    
    if (iproc==proc_X) then
       energyz=b_energytot+de_energytot+ne_energytot+phine_energytot+sumgsquaretot
  !     energyzloss=0.5*(energyz-energyz_old)/dt
       !NFL, 24/03/2013: why the 0.5 factor in the above line?
       energyzloss=(energyz-energyz_old)/dt
       !write(*,*) energyzloss
       energyz_old=energyz  ! This makes sense because of "save" LUIS
       !write(*,*) 'incalc'
       !write(*,*) energyzloss
    end if
!    deallocate(gk)

  end subroutine zsteploss


!*********************************************
  subroutine z_integration(array, integral, step, a_size)
     ! Subroutine to compute the integral of an array with elements
     !    array(1), array(2)...., array (a_size), in the z-direction 
     ! is assuming that "integral=0.0" before being called - Luis, Nov 2014
    use mp
    USE constants
    use grid
!    use transforms
    implicit none
    integer:: k, a_size
    real,dimension(a_size)::array
    real:: integral, step
    
    integral=0.0

    ! LUIS: This may be useful in the future, for non-periodic BC's
!      do k=1, a_size-1 ! finite difference scheme
!           integral = integral + step/2.*(array(k)+array(k+1))
!      end do
! 
!      if (npez.GT.1) then ! now add the bits missing from the integration boundaries -Luis, 17/11/14
!         integral = integral + step/2.*(array(1)+array(a_size))
! 
!         if (iproc.LT.npe) then ! fix the "overshoot" at 1st z-boundary - Luis
!            integral = integral - step/2.*array(1)
!         else if (iproc.GT.(npe*npez-npe-1)) then! fix the "overshoot" at 2nd z-boundary - Luis
!            integral = integral - step/2.*array(a_size)
!         endif
!      endif
       !LUIS: assuming periodic BC's
     do k=1,a_size
         integral = integral + step*array(k)
     end do

  end subroutine z_integration


!*********************************************
  subroutine diagnostics(FIK,nek,uekpar,uekpar_eq, &
       Akpar,Akpar_eq,gk,sumleftover,savetime,res2,niu2,nu_g,hyper_nuei,dt,proc_X,inj,energyzloss)
	!AVK 26/03/13: Changed rhos to rhos_diag. It is set in set_anjor. 
	! rhos_diag = rhos if anjor=false, otherwise rhos_diag = 1
    use mp
    USE constants
    use grid
!    use transforms
    implicit none

    
    COMPLEX, DIMENSION(NKy,NKX_PAR,nlz_par):: FIK, nek, uekpar, Akpar, uekpar_eq, Akpar_eq
    COMPLEX, DIMENSION(NKy,NKX_PAR,nlz_par):: temp2k
    !complex, dimension(nlz/2+1):: bkz_energy
    COMPLEX, DIMENSION(NKy,NKX_PAR,nlz_par,gmin:ngtot)::gk 
!    complex, allocatable, dimension(:,:,:,:)::gk
!    real, dimension(nlx,nly_par,nlz_par)::temp2
!    real, DIMENSION(NLx,NLy_PAR,nlz_par):: Apar, Apar_eq, uepar_eq, uepar
    real,dimension(nlz_par)::b_energy, de_energy, ne_energy, phine_energy, ohmic_diss
    real::b_energytot, de_energytot, ne_energytot, phine_energytot, ohmic_disstot
    real,dimension(nlz_par)::visc_diss, gm_diss, hyper_eta_diss, hyper_nu_diss, hyper_gm_diss, &
         sumgsquare, hyper_gm_kdiss
    real::visc_disstot, gm_disstot, hyper_eta_disstot, hyper_nu_disstot, hyper_gm_disstot, &
         sumgsquaretot, hyper_gm_kdisstot
    real::energy,savetime,res2,niu2,nu_g,dt,dwdt,hyper_nuei,inj,energyzloss,sumleftover
    real,save::energy_old,savetime_old
    logical,save::first=.true.
    integer::i,j,k, proc_X,ng

    if(first) then
       energy_old=0.0
       savetime_old=0.0
       first=.false.
    end if

    b_energy=0.0
    b_energytot=0.0
    de_energy=0.0
    de_energytot=0.0
    ne_energy=0.0
    ne_energytot=0.0
    phine_energy=0.0
    phine_energytot=0.0
    ohmic_diss=0.0
    ohmic_disstot=0.0
    visc_diss=0.0
    visc_disstot=0.0
    gm_diss=0.0
    gm_disstot=0.0
    hyper_eta_diss=0.0
    hyper_eta_disstot=0.0
    hyper_nu_diss=0.0
    hyper_nu_disstot=0.0
    hyper_gm_diss=0.0
    hyper_gm_disstot=0.0
    hyper_gm_kdiss=0.0
    hyper_gm_kdisstot=0.0
    sumgsquare=0.0
    sumgsquaretot=0.0

    do k=1,nlz_par
       do i=1, nkx_par
          do j=1, nky 
             !          temp2k(j,i)=kperp(j,i)**2*(uekpar(j,i)-uekpar_eq(j,i))
             temp2k(j,i,k)=kperp(j,i)**(2*hyper_order)*(Akpar(j,i,k)-Akpar_eq(j,i,k))
          end do
       end do
    end do
    
    if (mod(iproc,npe)==0) then
       do k=1,nlz_par
          do j=1, nky
             b_energy(k)=b_energy(k)+0.5*kperp(j,1)**2*Akpar(j,1,k)*conjg(akpar(j,1,k))
             de_energy(k)=de_energy(k)+0.5*de**2*uekpar(j,1,k)*conjg(uekpar(j,1,k))
             ne_energy(k)=ne_energy(k)+0.5*rhos_diag**2*nek(j,1,k)*conjg(nek(j,1,k))
             phine_energy(k)=phine_energy(k)-1.0/rhoi**2*(gama0(kperp(j,1)**2*rhoi**2/2.)-1.)*&
                  FIK(j,1,k)*conjg(FIK(j,1,k))
             
             !dissipation:
             ohmic_diss(k)=ohmic_diss(k)+res*(uekpar(j,1,k)-uekpar_eq(j,1,k))*&
                  CONJG(uekpar(j,1,k)-uekpar_eq(j,1,k))+&
                  res*uekpar_eq(j,1,k)*conjg((uekpar(j,1,k)-uekpar_eq(j,1,k)))
             visc_diss(k)=visc_diss(k)+niu*kperp(j,1)**2*(rhos_diag**2*nek(j,1,k)-FIK(j,1,k))*&
                  conjg(nek(j,1,k))

             if(g_inc) then
                do ng=gmin+1,ngtot
                   gm_diss(k)=gm_diss(k)+nu_ei*rhos_diag**2*ng*gk(j,1,k,ng)*conjg(gk(j,1,k,ng))
                   hyper_gm_diss(k)=hyper_gm_diss(k) + &
                        hyper_nuei*rhos_diag**2*ng**(2*hyper_morder)*gk(j,1,k,ng)*conjg(gk(j,1,k,ng))
                end do
                do ng=gmin,ngtot
                   sumgsquare(k)=sumgsquare(k)+0.5*rhos_diag**2*gk(j,1,k,ng)*conjg(gk(j,1,k,ng))
                   hyper_gm_kdiss(k)=hyper_gm_kdiss(k)+nu_g*rhos_diag**2*kperp(j,1)**(2*hyper_order_g)*&
                        gk(j,1,k,ng)*conjg(gk(j,1,k,ng))
                end do
             end if

             !          hyper_eta_diss(k)=hyper_eta+0.5*res2*kperp(j,1,k)**2*&
             !               (uekpar(j,1,k)-uekpar_eq(j,1,k))*conjg(uekpar(j,1,k)-uekpar_eq(j,1,k))&
             !               +0.5*res2*temp2k(j,1,k)
             hyper_eta_diss(k)=hyper_eta_diss(k)-res2*uekpar(j,1,k)*conjg(temp2k(j,1,k))
             hyper_nu_diss(k)=hyper_nu_diss(k) + niu2*kperp(j,1)**(2*hyper_order)*&
                  (rhos_diag**2*nek(j,1,k)-FIK(j,1,k))*conjg(nek(j,1,k))
             
          end do
          do i=2, nkx_par
             do j=1, nky
                b_energy(k)=b_energy(k)+kperp(j,i)**2*Akpar(j,i,k)*conjg(akpar(j,i,k))
                de_energy(k)=de_energy(k)+de**2*uekpar(j,i,k)*conjg(uekpar(j,i,k))
                ne_energy(k)=ne_energy(k)+rhos_diag**2*nek(j,i,k)*conjg(nek(j,i,k))
                phine_energy(k)=phine_energy(k)-2./rhoi**2*(gama0(kperp(j,i)**2*rhoi**2/2.)-1.)*&
                     FIK(j,i,k)*conjg(FIK(j,i,k))
                !    energy=energy-kperp(j,i)**2*FIK(j,i)&
                !         &*CONJG(FIK(j,i))+ uekpar(j,i)*CONJG(AKpar(j,i))
                
                !dissipation:
                ohmic_diss(k)=ohmic_diss(k)+2.0*res*(uekpar(j,i,k)-uekpar_eq(j,i,k))*&
                     CONJG(uekpar(j,i,k)-uekpar_eq(j,i,k))+&
                     2.0*res*uekpar_eq(j,i,k)*conjg((uekpar(j,i,k)-uekpar_eq(j,i,k)))
                !             visc_diss(k)=visc_diss(k)+niu*kperp(j,i)**2*nek(j,i)*conjg(nek(j,i))
                visc_diss(k)=visc_diss(k)+2.0*niu*kperp(j,i)**2*(rhos_diag**2*nek(j,i,k)-FIK(j,i,k))*&
                     conjg(nek(j,i,k))

                if(g_inc) then
                   do ng=gmin+1,ngtot
                      gm_diss(k)=gm_diss(k)+2.0*nu_ei*rhos_diag**2*ng*gk(j,i,k,ng)*conjg(gk(j,i,k,ng))
                      hyper_gm_diss(k)=hyper_gm_diss(k) + &
                           2.0*hyper_nuei*rhos_diag**2*ng**(2*hyper_morder)*gk(j,i,k,ng)*conjg(gk(j,i,k,ng))
                   end do
                   do ng=gmin,ngtot
                      sumgsquare(k)=sumgsquare(k)+rhos_diag**2*gk(j,i,k,ng)*conjg(gk(j,i,k,ng))
                      hyper_gm_kdiss(k)=hyper_gm_kdiss(k)+2.0*nu_g*rhos_diag**2*kperp(j,i)**(2*hyper_order_g)*&
                           gk(j,i,k,ng)*conjg(gk(j,i,k,ng))
                   end do
                end if
                
                !             hyper_eta_diss(k)=hyper_eta+res2*kperp(j,i)**2*&
                !                  (uekpar(j,i)-uekpar_eq(j,i))*conjg(uekpar(j,i)-uekpar_eq(j,i))&
                !                  +res2*temp2k(j,i)
                hyper_eta_diss(k)=hyper_eta_diss(k)-2.0*res2*uekpar(j,i,k)*conjg(temp2k(j,i,k))
                hyper_nu_diss(k)=hyper_nu_diss(k) + 2.0*niu2*kperp(j,i)**(2*hyper_order)*&
                     (rhos_diag**2*nek(j,i,k)-FIK(j,i,k))*conjg(nek(j,i,k))
             end do
          end do
       end do
    else
       do k=1,nlz_par
          do i=1,NKx_par
             do j=1,NKy
                b_energy(k)=b_energy(k)+kperp(j,i)**2*Akpar(j,i,k)*conjg(akpar(j,i,k))
                de_energy(k)=de_energy(k)+de**2*uekpar(j,i,k)*conjg(uekpar(j,i,k))
                ne_energy(k)=ne_energy(k)+rhos_diag**2*nek(j,i,k)*conjg(nek(j,i,k))
                phine_energy(k)=phine_energy(k)-2./rhoi**2*(gama0(kperp(j,i)**2*rhoi**2/2.)-1.)*&
                     FIK(j,i,k)*conjg(FIK(j,i,k))
                !    energy=energy-kperp(j,i)**2*FIK(j,i)&
                !         &*CONJG(FIK(j,i))+ uekpar(j,i)*CONJG(AKpar(j,i))
                
                !dissipation:
                ohmic_diss(k)=ohmic_diss(k)+2.0*res*(uekpar(j,i,k)-uekpar_eq(j,i,k))*&
                     CONJG(uekpar(j,i,k)-uekpar_eq(j,i,k))+&
                     2.0*res*uekpar_eq(j,i,k)*conjg((uekpar(j,i,k)-uekpar_eq(j,i,k)))
                !             visc_diss(k)=visc_diss(k)+niu*kperp(j,i)**2*nek(j,i)*conjg(nek(j,i))
                visc_diss(k)=visc_diss(k)+2.0*niu*kperp(j,i)**2*(rhos_diag**2*nek(j,i,k)-FIK(j,i,k))*&
                     conjg(nek(j,i,k))

                if(g_inc) then ! gmin=2
                   do ng=gmin+1,ngtot
                      gm_diss(k)=gm_diss(k)+2.0*nu_ei*rhos_diag**2*ng*gk(j,i,k,ng)*conjg(gk(j,i,k,ng))
                      hyper_gm_diss(k)=hyper_gm_diss(k) + &
                           2.0*hyper_nuei*rhos_diag**2*ng**(2*hyper_morder)*gk(j,i,k,ng)*conjg(gk(j,i,k,ng))
                   end do
                   do ng=gmin,ngtot
                      sumgsquare(k)=sumgsquare(k)+rhos_diag**2*gk(j,i,k,ng)*conjg(gk(j,i,k,ng))
                      hyper_gm_kdiss(k)=hyper_gm_kdiss(k)+2.0*nu_g*rhos_diag**2*kperp(j,i)**(2*hyper_order_g)*&
                           gk(j,i,k,ng)*conjg(gk(j,i,k,ng))
                   end do
                end if  

                !             hyper_eta_diss(k)=hyper_eta+res2*kperp(j,i)**2*&
                !                  (uekpar(j,i)-uekpar_eq(j,i))*conjg(uekpar(j,i)-uekpar_eq(j,i))&
                !                  +res2*temp2k(j,i)
                hyper_eta_diss(k)=hyper_eta_diss(k)-2.0*res2*uekpar(j,i,k)*conjg(temp2k(j,i,k))
                hyper_nu_diss(k)=hyper_nu_diss(k) + 2.0*niu2*kperp(j,i)**(2*hyper_order)*&
                     (rhos_diag**2*nek(j,i,k)-FIK(j,i,k))*conjg(nek(j,i,k))
             end do
          end do
       end do
    end if

    ! Z.Liu 7/2/2020
    !b_energy(:)=b_energy(:)/(nlx*nly*1.0)
    !de_energy(:)=de_energy(:)/(nlx*nly*1.0)
    !ne_energy(:)=ne_energy(:)/(nlx*nly*1.0)
    !phine_energy(:)=phine_energy(:)/(nlx*nly*1.0)
    !sumgsquare=sumgsquare/(nlx*nly*1.0)
    !ohmic_diss(:)=ohmic_diss(:)/(nlx*nly*1.0)
    !visc_diss(:)=visc_diss(:)/(nlx*nly*1.0)
    !gm_diss(:)=gm_diss(:)/(nlx*nly*1.0)
    !hyper_eta_diss(:)=hyper_eta_diss(:)/(nlx*nly*1.0)
    !hyper_nu_diss(:)=hyper_nu_diss(:)/(nlx*nly*1.0)
    !hyper_gm_diss(:)=hyper_gm_diss(:)/(nlx*nly*1.0)
    !hyper_gm_kdiss(:)=hyper_gm_kdiss(:)/(nlx*nly*1.0)

!    print*,iproc,b_energy
     
    if(three_D) then
          call z_integration(b_energy,b_energytot,dz,nlz_par) 
          call z_integration(de_energy,de_energytot,dz,nlz_par)
          call z_integration(ne_energy,ne_energytot,dz,nlz_par)
          call z_integration(phine_energy,phine_energytot,dz,nlz_par)
          call z_integration(sumgsquare,sumgsquaretot,dz,nlz_par)
          call z_integration(ohmic_diss,ohmic_disstot,dz,nlz_par) 
          call z_integration(visc_diss,visc_disstot,dz,nlz_par)
          call z_integration(gm_diss,gm_disstot,dz,nlz_par)
          call z_integration(hyper_eta_diss,hyper_eta_disstot,dz,nlz_par)
          call z_integration(hyper_nu_diss,hyper_nu_disstot,dz,nlz_par)
          call z_integration(hyper_gm_diss,hyper_gm_disstot,dz,nlz_par)
          call z_integration(hyper_gm_kdiss,hyper_gm_kdisstot,dz,nlz_par)
    else
       b_energytot = b_energy(1)
       de_energytot = de_energy(1)
       ne_energytot = ne_energy(1)
       phine_energytot = phine_energy(1)
       sumgsquaretot = sumgsquare(1)
       ohmic_disstot = ohmic_diss(1)
       visc_disstot = visc_diss(1)
       gm_disstot = gm_diss(1)
       hyper_eta_disstot = hyper_eta_diss(1)
       hyper_nu_disstot = hyper_nu_diss(1)
       hyper_gm_disstot = hyper_gm_diss(1)
       hyper_gm_kdisstot = hyper_gm_kdiss(1)
    end if

!    b_energytot=b_energytot/(nlz*1.0)
!    de_energytot=de_energytot/(nlz*1.0)
!    ne_energytot=ne_energytot/(nlz*1.0)
!    phine_energytot=phine_energytot/(nlz*1.0)
!    sumgsquaretot=sumgsquaretot/(nlz*1.0)
!    ohmic_disstot=ohmic_disstot/(nlz*1.0)
!    visc_disstot=visc_disstot/(nlz*1.0)
!    gm_disstot=gm_disstot/(nlz*1.0)
!    hyper_eta_disstot=hyper_eta_disstot/(nlz*1.0)
!    hyper_nu_disstot=hyper_nu_disstot/(nlz*1.0)
!    hyper_gm_disstot=hyper_gm_disstot/(nlz*1.0)
!    hyper_gm_kdisstot=hyper_gm_kdisstot/(nlz*1.0)

    call sum_reduce (b_energytot, proc_X)
    call sum_reduce (de_energytot, proc_X)
    call sum_reduce (ne_energytot, proc_X)
    call sum_reduce (phine_energytot, proc_X)
    call sum_reduce (sumgsquaretot, proc_X)
    call sum_reduce (ohmic_disstot,proc_X)
    call sum_reduce (visc_disstot, proc_X)
    call sum_reduce (gm_disstot,proc_X)
    call sum_reduce (hyper_eta_disstot, proc_X)
    call sum_reduce (hyper_nu_disstot, proc_X)
    call sum_reduce (hyper_gm_disstot, proc_X)
    call sum_reduce (hyper_gm_kdisstot, proc_X)
    
    if (iproc==proc_X) then
       energy=b_energytot+de_energytot+ne_energytot+phine_energytot+sumgsquaretot
       dwdt=(energy-energy_old)/(savetime-savetime_old)
       write(88,39) savetime, b_energytot,de_energytot,ne_energytot,phine_energytot,&
            sumgsquaretot,dwdt,inj,ohmic_disstot,visc_disstot,gm_disstot,hyper_eta_disstot,&
            hyper_nu_disstot,hyper_gm_disstot, hyper_gm_kdisstot,energyzloss,sumleftover
       energy_old=energy    !This makes sense because of "save", in variable declaration LUIS
       savetime_old=savetime
    end if
39  format(17g16.8)

  end subroutine diagnostics


!*********************************************
  subroutine isl_width(Apar,proc_X,y_loc_loc,a_X,W2,pos)
    use constants
    use grid
    use mp
    implicit none

    real, dimension(NLx,NLy_par)::Apar
    real::W2
    real::a_X,Y
    real::epslon
    integer::i,pos,proc_X,y_loc_loc
    
    epslon=1e-5

    Y=0.0
    i=x_loc    
    do while (abs((Y-a_X)/a_X).GT.epslon)
       Y=Apar(i,1)
       i=i+1
       if (i==Nlx) then
          epslon=2*epslon
          i=x_loc
       end if
    end do
    
    W2=2*(XX(i-1)-XX(x_loc))
    pos=i-1-x_loc
  end subroutine isl_width


!*********************************************
  subroutine SP(vey,by,p,proc_X,y_loc_loc,outflow_max,bup_max,bup_jemella)
    use constants
    use mp
    use grid
    implicit none
    
    real, DIMENSION(NLx,NLy_par) ::vey,by
    real,dimension(NLx)::dum1
    real,dimension(NLy_par)::dum2
    real,dimension(nly)::work
    real::bup_max, bup_jemella, outflow_max
    integer::p,proc_X,y_loc_loc,m,j
    
    if (iproc==proc_X) then
       dum1(:)=by(:,y_loc_loc)
       bup_max=maxval(dum1)
       bup_jemella=by(p,y_loc_loc)
    end if
   
    work=0.
    m=1+iproc*nly_par
    do j=1,nly_par
       work(m)=vey(x_loc,j)
       m=m+1
    end do

    call sum_reduce(work,0)
    if (proc0) outflow_max=maxval(work)
    
  end subroutine SP


!*********************************************
!  subroutine layer_width(uepar_perturb,proc_X,y_loc_loc,delta,p)
!    use constants
!    use mp
!    use grid
!    implicit none
    
!    real,dimension(nlx,nly_par),intent(in)::uepar_perturb
!    real,dimension(nlx/2)::temp1,temp2
!    real::half_max,min1,min2
!    real,intent(out)::delta
!    type(r_layout_type)::where
!    integer::pos(1),i,n,proc_X,y_loc_loc
!    integer,intent(out)::p
!    
!    
! !   if (iproc==proc_X) then
!       half_max=uepar_perturb(x_loc,y_loc_loc)/2.
!       i=1
!       do while((uepar_perturb(x_loc-i,y_loc_loc) .LT. 0) .OR. &
!            & (uepar_perturb(x_loc+i,y_loc_loc) .LT. 0))
!          temp1(i)=abs(uepar_perturb(x_loc-i,y_loc_loc)-half_max)
!          temp2(i)=abs(uepar_perturb(x_loc+i,y_loc_loc)-half_max)
!          i=i+1
!       end do
!       
!       do n=i, nlx/2
!          temp1(n)=1000
!          temp2(n)=1000
!       end do
!       
!       min1=minval(temp1)
!       min2=minval(temp2)
!       
!       if (min1 .LE. min2) then
!          pos=minloc(temp1)
!       else
!          pos=minloc(temp2)
!       end if
!       
!       delta=abs(xx(x_loc+pos(1))-xx(x_loc))*2.
!       p=x_loc-4*pos(1)
!  !  end if
!    
!  end subroutine layer_width


!*********************************************
  subroutine layer_width(uepar_perturb,y_loc_loc,delta,p)
    use constants
    use mp
    use grid
    implicit none
    
    real,dimension(nlx,nly_par),intent(in)::uepar_perturb
    real,dimension(nlx)::work
    real::half_max,min1,min2,Y,epslon
    real,intent(out)::delta
    integer::pos(1),i,n,y_loc_loc
    integer,intent(out)::p
    
    half_max=uepar_perturb(x_loc,y_loc_loc)/2.

    epslon=1e-5
    Y=0.0
    i=x_loc    
    do while (abs((Y-half_max)/half_max).GT.epslon)
       Y=uepar_perturb(i,y_loc_loc)
       i=i-1
!       if (i==1) then
       if ((Y/half_max .LT. 0).OR. i==0) then
          epslon=2*epslon
          i=x_loc
       end if
    end do

    delta=2*(xx(x_loc)-xx(i+1))
    p=i+1
 
  end subroutine layer_width

!*********************************************  
!  subroutine sheet_length(uepar_perturb,L_sheet)
!    use constants
!    use mp
!    use grid, only: yy
!    
!    implicit none
!    real,dimension(nlx,nly_par),intent(in)::uepar_perturb
!    real,dimension(nly/2)::temp1,temp2
!    real,dimension(nly)::work
!    real::half_max,min1,min2
!    real,intent(out)::L_sheet
!    integer::pos(1),j,n,m
!    
!    work=0.
!    m=1+iproc*nly_par
!    do j=1,nly_par
!       work(m)=uepar_perturb(x_loc,j)
!       m=m+1
!    end do
!    
!    call sum_reduce(work,0)
!    
!    if (proc0) then
!       half_max=work(y_loc)/2.
!       j=1
!       do while((work(y_loc-j) .LT. 0) .OR. &
!            &(work(y_loc+j) .LT. 0))
!          temp1(j)=abs(work(y_loc-j)-half_max)
!          temp2(j)=abs(work(y_loc+j)-half_max)
!          j=j+1
!       end do
!       
!       do n=j, nly/2
!          temp1(n)=1000
!          temp2(n)=1000
!       end do
!       
!       min1=minval(abs(temp1))
!       min2=minval(abs(temp2))
!       
!       If(min1 .LE. min2) Then
!          pos=minloc(temp1)
!       else
!          pos=minloc(temp2)
!       end if
!       
!       L_sheet=abs(yy(y_loc+pos(1))-yy(y_loc))*2.
!    end if
!    
!  end subroutine sheet_length
!************************************
  subroutine sheet_length(uepar_perturb,L_sheet)
    use constants
    use mp
    use grid, only: yy
    
    implicit none
    real,dimension(nlx,nly_par),intent(in)::uepar_perturb
    real,dimension(nly)::temp1,temp2
    real,dimension(nly)::work
    real::half_max,min1,min2,Y,epslon
    real,intent(out)::L_sheet
    integer::pos(1),j,n,m,i
    
    work=0.
    m=1+iproc*nly_par
    do j=1,nly_par
       work(m)=uepar_perturb(x_loc,j)
       m=m+1
    end do
    
    call sum_reduce(work,0)
    
    if (proc0) then
       half_max=work(y_loc)/2.
       work(:)=abs(work(:)-half_max)
       pos=minloc(work)
       L_sheet=2.*abs(yy(y_loc)-yy(pos(1)))
    end if
    
  end subroutine sheet_length

!************************************
  subroutine sheet_length2(vey,LCS2)
    use constants
    use mp
    use grid
    implicit none
    real,dimension(nlx,nly_par)::vey
    real,dimension(NPE)::vout_max
    real,dimension(nly_par)::dum_vey
    integer::k
    integer,dimension(1)::voutmax_loc,jloc
    real::LCS2

    dum_vey(:)=abs(vey(x_loc,:))
    vout_max=0.0
    do k=0,NPE-1
       if (iproc==k) then
          vout_max(k+1)=maxval(dum_vey)
       end if
    end do
    call sum_allreduce(vout_max)
    voutmax_loc=maxloc(vout_max)
    if (iproc==voutmax_loc(1)-1) then
       jloc=maxloc(dum_vey)
       LCS2=2.0*yy(jloc(1))
       call send(LCS2,0)   
    end if
    if (proc0) call receive(LCS2,voutmax_loc(1)-1)
   
  end subroutine sheet_length2


!*********************************************
  subroutine X_point(apar,uepar_perturb,y_loc_loc,a_X,u_X)
    !X-point quantities and initialization
    use constants
    use mp
    implicit none
    
    real,dimension(nlx,nly_par),intent(in)::apar,uepar_perturb
    real,intent(out)::a_X,u_X
    integer,intent(in)::y_loc_loc
    
    a_X=apar(x_loc,y_loc_loc)
    u_X=uepar_perturb(x_loc,y_loc_loc)
        
  end subroutine X_point


!*********************************************
  subroutine file_name(runname, time,length,F3,Fg,time_str_)

    use constants, only: len3, PATH, fieldsfile, gfieldsfile, g_inc
!    use mp

    implicit none 

    integer,intent(in)::length
    real, intent (in)::time
    !    character(len=len3+3+length+4+4)::F1
    !    character(len=len3+5+length+4+4)::F2
    !character(len=len3+7+length+4+4),intent(out)::F3
    !character(len=len3+8+length+4+4),intent(out)::Fg
    character(len=100),intent(out)::F3
    character(len=100),intent(out)::Fg
!    character(len=len3+5+length+4+4)::F4
    character(len=200)::MF3,MMF3,MFg,MMFg
	character(len=20)::time_str
	character(len=20),intent(out),optional::time_str_
	character(len=100), intent(in) :: runname
    
    if (time .lt. 10) then
       !write(MF3,'(a100,a7,f5.3,a4)') PATH,fieldsfile,time,".dat"
	   write(time_str,'(f5.3)') time
	   MF3 = trim(PATH)//trim(runname)//trim(fieldsfile)//trim(time_str)//".dat"
       MMF3=adjustl(MF3)
       F3=trim(MMF3)
       if(g_inc) then
          !write(MFg,'(a100,a8,f5.3,a4)') PATH,gfieldsfile,time,".dat"
		  MFg = trim(PATH)//trim(runname)//trim(gfieldsfile)//trim(time_str)//".dat"
          MMFg=adjustl(MFg)
          Fg=trim(MMFg)
       end if
    end if
    if (time .ge. 10 .AND. time .lt. 100) then
       !write(MF3,'(a100,a7,f6.3,a4)') PATH,fieldsfile,time,".dat"
	   write(time_str,'(f6.3)') time
	   MF3 = trim(PATH)//trim(runname)//trim(fieldsfile)//trim(time_str)//".dat"
       MMF3=adjustl(MF3)
       F3=trim(MMF3)
       if(g_inc) then
          !write(MFg,'(a100,a8,f6.3,a4)') PATH,gfieldsfile,time,".dat"
		  MFg = trim(PATH)//trim(runname)//trim(gfieldsfile)//trim(time_str)//".dat"
          MMFg=adjustl(MFg)
          Fg=trim(MMFg)
       end if
    end if
    if (time .ge. 100 .AND. time .lt. 1000) then
       !write(MF3,'(a100,a7,f7.3,a4)') PATH,fieldsfile,time,".dat"
	   write(time_str,'(f7.3)') time
	   MF3 = trim(PATH)//trim(runname)//trim(fieldsfile)//trim(time_str)//".dat"
       MMF3=adjustl(MF3)
       F3=trim(MMF3)
       if(g_inc) then
          !write(MFg,'(a100,a8,f7.3,a4)') PATH,gfieldsfile,time,".dat"
		  MFg = trim(PATH)//trim(runname)//trim(gfieldsfile)//trim(time_str)//".dat"
          MMFg=adjustl(MFg)
          Fg=trim(MMFg)
       end if
    end if
    if (time .ge. 1000 .AND. time .lt. 10000) then
       !write(MF3,'(a100,a7,f8.3,a4)') PATH,fieldsfile,time,".dat"
	   write(time_str,'(f8.3)') time
	   MF3 = trim(PATH)//trim(runname)//trim(fieldsfile)//trim(time_str)//".dat"
       MMF3=adjustl(MF3)
       F3=trim(MMF3)
       if(g_inc) then
          !write(MFg,'(a100,a8,f8.3,a4)') PATH,gfieldsfile,time,".dat"
		  MFg = trim(PATH)//trim(runname)//trim(gfieldsfile)//trim(time_str)//".dat"
          MMFg=adjustl(MFg)
          Fg=trim(MMFg)
       end if
    end if
    if (time .ge. 10000 .AND. time .lt. 100000) then
       !write(MF3,'(a100,a7,f9.3,a4)') PATH,fieldsfile,time,".dat"
	   write(time_str,'(f9.3)') time
	   MF3 = trim(PATH)//trim(runname)//trim(fieldsfile)//trim(time_str)//".dat"
       MMF3=adjustl(MF3)
       F3=trim(MMF3)
       if(g_inc) then
          !write(MFg,'(a100,a8,f9.3,a4)') PATH,gfieldsfile,time,".dat"
		  MFg = trim(PATH)//trim(runname)//trim(gfieldsfile)//trim(time_str)//".dat"
          MMFg=adjustl(MFg)
          Fg=trim(MMFg)
       end if
    end if
    if( present ( time_str_ ) ) then
      time_str_ = time_str
    end if

  end subroutine file_name

!*********************************************
subroutine kfile_name(runname, time,length,F3,Fg,time_str_)

use constants, only: len3, PATH, kfieldsfile, gfieldsfile, g_inc
!    use mp

implicit none

integer,intent(in)::length
real, intent (in)::time
!    character(len=len3+3+length+4+4)::F1
!    character(len=len3+5+length+4+4)::F2
!character(len=len3+7+length+4+4),intent(out)::F3
!character(len=len3+8+length+4+4),intent(out)::Fg
character(len=100),intent(out)::F3
character(len=100),intent(out)::Fg
!    character(len=len3+5+length+4+4)::F4
character(len=200)::MF3,MMF3,MFg,MMFg
character(len=20)::time_str
character(len=20),intent(out),optional::time_str_
character(len=100), intent(in) :: runname

if (time .lt. 10) then
!write(MF3,'(a100,a7,f5.3,a4)') PATH,fieldsfile,time,".dat"
write(time_str,'(f5.3)') time
MF3 = trim(PATH)//trim(runname)//trim(kfieldsfile)//trim(time_str)//".dat"
MMF3=adjustl(MF3)
F3=trim(MMF3)
if(g_inc) then
!write(MFg,'(a100,a8,f5.3,a4)') PATH,gfieldsfile,time,".dat"
MFg = trim(PATH)//trim(runname)//trim(gfieldsfile)//trim(time_str)//".dat"
MMFg=adjustl(MFg)
Fg=trim(MMFg)
end if
end if
if (time .ge. 10 .AND. time .lt. 100) then
!write(MF3,'(a100,a7,f6.3,a4)') PATH,fieldsfile,time,".dat"
write(time_str,'(f6.3)') time
MF3 = trim(PATH)//trim(runname)//trim(kfieldsfile)//trim(time_str)//".dat"
MMF3=adjustl(MF3)
F3=trim(MMF3)
if(g_inc) then
!write(MFg,'(a100,a8,f6.3,a4)') PATH,gfieldsfile,time,".dat"
MFg = trim(PATH)//trim(runname)//trim(gfieldsfile)//trim(time_str)//".dat"
MMFg=adjustl(MFg)
Fg=trim(MMFg)
end if
end if
if (time .ge. 100 .AND. time .lt. 1000) then
!write(MF3,'(a100,a7,f7.3,a4)') PATH,fieldsfile,time,".dat"
write(time_str,'(f7.3)') time
MF3 = trim(PATH)//trim(runname)//trim(kfieldsfile)//trim(time_str)//".dat"
MMF3=adjustl(MF3)
F3=trim(MMF3)
if(g_inc) then
!write(MFg,'(a100,a8,f7.3,a4)') PATH,gfieldsfile,time,".dat"
MFg = trim(PATH)//trim(runname)//trim(gfieldsfile)//trim(time_str)//".dat"
MMFg=adjustl(MFg)
Fg=trim(MMFg)
end if
end if
if (time .ge. 1000 .AND. time .lt. 10000) then
!write(MF3,'(a100,a7,f8.3,a4)') PATH,fieldsfile,time,".dat"
write(time_str,'(f8.3)') time
MF3 = trim(PATH)//trim(runname)//trim(kfieldsfile)//trim(time_str)//".dat"
MMF3=adjustl(MF3)
F3=trim(MMF3)
if(g_inc) then
!write(MFg,'(a100,a8,f8.3,a4)') PATH,gfieldsfile,time,".dat"
MFg = trim(PATH)//trim(runname)//trim(gfieldsfile)//trim(time_str)//".dat"
MMFg=adjustl(MFg)
Fg=trim(MMFg)
end if
end if
if (time .ge. 10000 .AND. time .lt. 100000) then
!write(MF3,'(a100,a7,f9.3,a4)') PATH,fieldsfile,time,".dat"
write(time_str,'(f9.3)') time
MF3 = trim(PATH)//trim(runname)//trim(kfieldsfile)//trim(time_str)//".dat"
MMF3=adjustl(MF3)
F3=trim(MMF3)
if(g_inc) then
!write(MFg,'(a100,a8,f9.3,a4)') PATH,gfieldsfile,time,".dat"
MFg = trim(PATH)//trim(runname)//trim(gfieldsfile)//trim(time_str)//".dat"
MMFg=adjustl(MFg)
Fg=trim(MMFg)
end if
end if
if( present ( time_str_ ) ) then
time_str_ = time_str
end if

end subroutine kfile_name

!*********************************************

  subroutine Convol(Fk, DxF, DyF)

    use constants,  only: nkx_par, nky, nlx, nly_par, nlz_par
    use grid,       only: kx, ky
    use transforms, only: FFT2d_inv

    implicit none

    complex, intent(in),  dimension(nky, nkx_par, nlz_par) :: Fk
    real,    intent(out), dimension(nlx, nly_par, nlz_par) :: DxF, DyF

    complex, dimension(nky, nkx_par, nlz_par) :: Fk_ikx, Fk_iky
    integer :: i, j, k

    Myperfon('Convol22')

    do k = 1, nlz_par
       do i = 1, nkx_par
          do j = 1, nky
             Fk_ikx(j, i, k) = cmplx(0.0, 1.0) * kx(i) * Fk(j, i, k)
             Fk_iky(j, i, k) = cmplx(0.0, 1.0) * ky(j) * Fk(j, i, k)
          end do
       end do

       !TTR
       call FFT2d_inv (Fk_ikx(:, :, k), DxF(:, :, k))
       call FFT2d_inv (Fk_iky(:, :, k), DyF(:, :, k))

    end do

    Myperfoff ! Convol22
    
  end subroutine Convol


!TTR
!*********************************************
  subroutine Convol2_3(Fk, DxF, DyF)

    use constants,  only: nkx_par, nky, nlz_par
    use grid,       only: kx, ky
    use transforms, only: FFT2d_inv

    implicit none

    !. size is (nky, nkx_par, nlz_par)
    complex, intent(in),  dimension(:, :, :) :: Fk
    !. size is (nlx, nly_par, nlz_par)
    real,    intent(out), dimension(:, :, :) :: DxF, DyF

    complex, dimension(nky, nkx_par, nlz_par) :: Fk_ikx, Fk_iky
    integer :: i, j, k

    Myperfon('Convol23')

    do k = 1, nlz_par
       do i = 1, nkx_par
          do j = 1, nky
             Fk_ikx(j, i, k) = cmplx(0.0, 1.0) * kx(i) * Fk(j, i, k)
             Fk_iky(j, i, k) = cmplx(0.0, 1.0) * ky(j) * Fk(j, i, k)
          end do
       end do
    end do

    call FFT2d_inv (Fk_ikx, DxF)
    call FFT2d_inv (Fk_iky, DyF)

    Myperfoff ! Convol23

  end subroutine Convol2_3


!TTR
!*********************************************
  subroutine Convol2_4(Fk, DxF, DyF)

    use constants,  only: nkx_par, nky, nlz_par, gmin, ngtot
    use grid,       only: kx, ky
    use transforms, only: FFT2d_inv

    implicit none
    
    !. size is (nky, nkx_par, nlz_par, gmin:ngtot)
    complex, intent(in),  dimension(:, :, :, gmin:) :: Fk
    !. size is (nlx, nly_par, nlz_par, gmin:)
    real,    intent(out), dimension(:, :, :, gmin:) :: DxF, DyF

    complex, dimension(nky, nkx_par, nlz_par, gmin:ngtot) :: Fk_ikx, Fk_iky
    integer :: i, j, k

    Myperfon('Convol24')

    do k = 1, nlz_par
       do i = 1, nkx_par
          do j = 1, nky
             Fk_ikx(j, i, k, :) = cmplx(0.0, 1.0) * kx(i) * Fk(j, i, k, :)
             Fk_iky(j, i, k, :) = cmplx(0.0, 1.0) * ky(j) * Fk(j, i, k, :)
          end do
       end do
    end do

    call FFT2d_inv (Fk_ikx, DxF)
    call FFT2d_inv (Fk_iky, DyF)

    Myperfoff ! Convol24
    
  end subroutine Convol2_4


!*********************************************
  subroutine gm_spectrum(unitnumber, time, gk, runname)
!NFL, 31/01/10; outputs gm spectrum
    use constants
    use mp
    implicit none

    character(len=100), intent(in) :: runname

    complex, dimension(nky, nkx_par, nlz_par, gmin:ngtot) :: gk
    real, dimension(nlz, gmin:ngtot) :: ek_gm
    real :: time, ek_gm_tot(gmin:ngtot)
    character(len=100) :: filename1
    integer :: i, j, k, m, unitnumber

    ek_gm=0.0
    ek_gm_tot=0.0
    !the kx=0 modes only get added once
    !the kx diff 0 modes get added twice because of the reality condition
    do m = gmin, ngtot
       if(proc0) then
          do k = 1, nlz_par
             do j = 1, nky
                ek_gm(k, m) = ek_gm(k, m) + 0.5*abs(gk(j, 1, k, m))**2*rhos_diag**2 !Z.Liu 7/8/3030
             end do
             do i = 2, nkx_par
                do j = 1, nky
                   ek_gm(k, m) = ek_gm(k, m) + abs(gk(j, i, k, m))**2*rhos_diag**2 !Z.Liu 7/8/3030
                end do
             end do
          end do
       else
          do k = 1, nlz_par
             do i = 1, nkx_par
                do j = 1, nky
                   ek_gm(k, m) = ek_gm(k, m) + abs(gk(j, i, k, m))**2*rhos_diag**2 !Z.Liu 7/8/3030
                end do
             end do
          end do
       end if
       ek_gm = ek_gm/(nlx*nly*1.0)
       !now integrate in z:
       if(three_D) then 
          do k = 1, nlz_par-1
             ek_gm_tot(m) = ek_gm_tot(m) + dz/2.*(ek_gm(k,m)+ek_gm(k+1,m))
          end do
       else
          ek_gm_tot(m) = ek_gm(1, m)
       end if
    end do
    call sum_reduce(ek_gm_tot, 0)

    if (proc0) then

       if (time.lt.10) then
          write(filename1,'(a9,f5.3,a4)') energy_gm,time,".dat"
       end if
       if (time .ge. 10 .AND. time .lt. 100) then
          write(filename1,'(a9,f6.3,a4)') energy_gm,time,".dat"
       end if
       if (time .ge. 100 .AND. time .lt. 1000) then
          write(filename1,'(a9,f7.3,a4)') energy_gm,time,".dat"
       end if
       if (time .ge. 1000 .AND. time .lt. 10000) then
          write(filename1,'(a9,f8.3,a4)') energy_gm,time,".dat"
       end if

       filename1 = trim(runname)//"_"//trim(filename1)
       open (unit=unitnumber, file= trim(filename1))
       
       do m = gmin, ngtot
          write(unitnumber,*) m, ek_gm_tot(m)
       end do

    end if

  end subroutine gm_spectrum


!*********************************************
  subroutine gm_timetrace(time,gk)
    !AVK, 13/09/12; outputs gm timetrace for slow mode runs
    use constants
    use mp
    implicit none

    complex, dimension(nky,nkx_par,nlz_par,gmin:ngtot)::gk
    real, dimension(nlz)::ek_gm
    real::time,ek_gm_tot, gm_inv
    !character(len=21)::filename1
    integer::i,j,k,m
    
    if (proc0) then
       write(99, 123, advance='no')  time
    end if

    gm_inv = 0.0
    
    !gmin done separately for the total energy
    ek_gm=0.0
    ek_gm_tot=0.0
    !the kx=0 modes only get added once
    !the kx diff 0 modes get added twice because of the reality condition
    if(proc0) then
       do k=1,nlz_par
          do j=1, nky
             ek_gm(k)=ek_gm(k)+0.5*abs(gk(j,1,k,gmin))**2*rhos_diag**2
          end do
          do i=2,NKx_par
             do j=1,NKy
                ek_gm(k)=ek_gm(k)+abs(gk(j,i,k,gmin))**2*rhos_diag**2
             end do
          end do
       end do
    else
       do k=1,nlz_par
          do i=1,NKx_par
             do j=1,NKy
                ek_gm(k)=ek_gm(k)+abs(gk(j,i,k,gmin))**2*rhos_diag**2
             end do
          end do
       end do
    end if
    ek_gm=ek_gm/(nlx*nly*1.0) ! AVK: why not divide by nlz as well?
    !now integrate in z:
    if(three_d) then
       do k=1,nlz_par-1
          ek_gm_tot = ek_gm_tot + dz/2.*(ek_gm(k)+ek_gm(k+1))
       end do
    else
       ek_gm_tot = ek_gm(1)
    end if
    call sum_reduce(ek_gm_tot,0)
    
    if (proc0) then
       write(99, 123, advance='no') ek_gm_tot
    end if
    
    gm_inv = gm_inv + ek_gm_tot*(1.0-1.0/lambda)
    
    do m = gmin+1, ngtot
       ek_gm=0.0
       ek_gm_tot=0.0
       !the kx=0 modes only get added once
       !the kx diff 0 modes get added twice because of the reality condition
       if(proc0) then
          do k=1,nlz_par
             do j=1, nky
                ek_gm(k)=ek_gm(k)+0.5*abs(gk(j,1,k,m))**2*rhos_diag**2 !Z.LIU 
             end do
             do i=2,NKx_par
                do j=1,NKy
                   ek_gm(k)=ek_gm(k)+abs(gk(j,i,k,m))**2*rhos_diag**2
                end do
             end do
          end do
       else
          do k=1,nlz_par
             do i=1,NKx_par
                do j=1,NKy
                   ek_gm(k)=ek_gm(k)+abs(gk(j,i,k,m))**2*rhos_diag**2
                end do
             end do
          end do
       end if
       ek_gm=ek_gm/(nlx*nly*1.0) ! AVK: why not divide by nlz as well?
       !now integrate in z:
       if(three_d) then
          do k=1,nlz_par-1
             ek_gm_tot = ek_gm_tot + dz/2.*(ek_gm(k)+ek_gm(k+1))
          end do
       else
          ek_gm_tot = ek_gm(1)
       end if
       call sum_reduce(ek_gm_tot,0)
       
       if (proc0) then
          write(99, 123, advance='no') ek_gm_tot
       end if
       
       gm_inv = gm_inv + ek_gm_tot
    end do

    if(proc0 .and. anjor) then 
       write(99, 123) gm_inv
    end if
    !write(99, *) '  '
123 format(es15.6)
  end subroutine gm_timetrace


!*********************************************
  subroutine layer_xy_widths (uekpar, uepar_eq, y_loc_loc, proc_X, delta_x, delta_y)

    use constants
    use grid
    use mp
    use transforms, only: FFT2d_inv
    implicit none

    complex, dimension(nky, nkx_par), intent(in) :: uekpar
    real,    dimension(nlx, nly_par), intent(in) :: uepar_eq
    integer, intent(in) :: y_loc_loc, proc_X

    real,    dimension(nlx, nly_par) :: uepar_double_prime_x, uepar_double_prime_y
    complex, dimension(nky, nkx_par) :: uekpar_double_prime_x, uekpar_double_prime_y

    real, intent(out) :: delta_x, delta_y
    integer:: i, j

    do i = 1, nkx_par
       do j = 1, nky
          uekpar_double_prime_x(j, i) = -kx(i)**2 * uekpar(j, i)
          uekpar_double_prime_y(j, i) = -ky(j)**2 * uekpar(j, i)
       end do
    end do
    !TTR
    call FFT2d_inv (uekpar_double_prime_x, uepar_double_prime_x)
    call FFT2d_inv (uekpar_double_prime_y, uepar_double_prime_y)

    if(iproc==proc_X) then 
       delta_x = sqrt(abs(uepar_eq(x_loc, y_loc_loc)/uepar_double_prime_x(x_loc,y_loc_loc)))
       delta_y = sqrt(abs(uepar_eq(x_loc, y_loc_loc)/uepar_double_prime_y(x_loc,y_loc_loc)))
    end if

    call broadcast (delta_x, proc_X)
    call broadcast (delta_y, proc_X)
  end subroutine layer_xy_widths


!*********************************************
  subroutine gamma_ky(time,nek,apar,runname)
    use constants
    use grid
    use mp
    use transforms
    implicit none

    real::time,gamma1(nky),gamma2(nky),omega(nky)
!   real, dimension (nly) :: workapar
    real, dimension(nlx,nly_par,nlz_par),intent (in)::apar
    real,save::time_old
    real::dt
    integer::j,i,k,m,ip
    complex, dimension(nky,nkx_par,nlz_par):: nek
    complex,dimension(nky,nkx)::s1tot
    complex, dimension (nky*nkx_par*NPE) :: work
    complex, dimension(:,:,:),allocatable,save:: nek_old
!   complex,dimension(nly/2+1)::ak
!   complex,dimension(:),allocatable,save::ne_old1,ne_old2,ak_old
    logical, save::first=.true.
    character(len=100),intent(in)::runname

    if (first) then
       grates=trim(runname)//trim(grates)
       open (unit=65, file=trim(grates))
!      open (unit=66, file='grates_raw.txt')
!      open (unit=67, file='grates_ak.txt')
!      open (unit=68, file='grates_ak_raw.txt')
       first=.false.
       time_old=0.0
!      allocate(ne_old1(nky))
!      allocate(ne_old2(nky))
!      allocate(ak_old(nly/2+1))
!      ne_old1=1e-20
!      ne_old2=1e-20
!      ak_old=1.e-20
       allocate(nek_old(nky,nkx_par,nlz_par))
       nek_old=0.0
    end if

    dt=time-time_old

!   if(proc0) then
!      do j=1,nky
!         gamma1(j) = (log(abs(nek(j,2,1)))-log(abs(ne_old1(j))))/(time-time_old)
!         gamma2(j) = (log(abs(nek(j,3,1)))-log(abs(ne_old2(j))))/(time-time_old)
!         write(65,'(a2,g16.8,a7,g16.8,a11,2g16.8,a7,2g16.8)') &
!              't=', time, ' ky=', ky(j), ' ne_energy=', abs(nek(j,2,1)), abs(nek(j,3,1)),&
!              ' gamma=', gamma1(j), gamma2(j)
!         write(66,'(6g16.8)') time,ky(j),abs(nek(j,2,1)),abs(nek(j,3,1)),gamma1(j),gamma2(j)
!      end do
!      ne_old1(:)=nek(:,2,1)
!      ne_old2(:)=nek(:,3,1)
!
!   end if
!
!   workapar=0.0
!   m=1+iproc*nly_par
!    do j=1,nly_par
!       workapar(m)=apar(x_loc,j,1)
!       m=m+1
!    end do
!
!    call sum_reduce(workapar,0)
!    if(proc0) then
!       call oneDfourndirect(workapar,ak,nly)
!       do j=1,nly/2+1
!          gamma1(j) =  (log(abs(ak(j)))-log(abs(ak_old(j))))/(time-time_old)
!          write(67,'(a2,g16.8,a4,g16.8,a8,g16.8,a7,g16.8)') &
!              't=', time, ' ky=', ky(j), ' abs_ak=', abs(ak(j)),&
!              ' gamma=', gamma1(j)
!          write(68,'(4g16.8)') time, ky(j), abs(ak(j)),gamma1(j)
!       end do
!       ak_old(:)=ak(:)
!    end if
    
    !> NFL 20/05/2014
    work = 0.
    if(iproc<npe) then
       m = 1+iproc*nky*nkx_par
       do i=1,nkx_par
          do j=1,nky
             work(m)=nek(j,i,1)
             m = m + 1
          end do
       end do
    end if
    call sum_reduce(work,0)
   
    if (proc0) then
       m=1
       do i=1,nkx
          do j=1,nky
             write(65,'(6g16.8)') time, ky(j), i-1., real(work(m)), aimag(work(m)), iproc
             m=m+1
          end do
       end do
    end if
   
    !<
   
    nek_old=nek
    time_old=time

  end subroutine gamma_ky


!*********************************************
  subroutine zonalflows(unitnumber,time,p,phik,ak,runname)
    use mp
    use constants
    use transforms
    use grid
    implicit none

    integer::i,m,p,unitnumber
    complex,dimension(nkx_par)::phik,ak
    complex,dimension(nkx)::workphi,workapar
    real,dimension(nlx)::phi,apar
    real::time
    !logical, save::first=.true.
    character(len=100)::filename1
    character(len=100),intent(in)::runname

   if (proc0 .AND. p .GT. 1  .AND. mod(p,save_energyfiles)==0) then
       if (time.lt.10) then
          write(filename1,'(a6,f5.3,a4)') ,time,".dat"
       end if
       if (time .ge. 10 .AND. time .lt. 100) then
          write(filename1,'(a6,f6.3,a4)') zonal,time,".dat"
       end if
       if (time .ge. 100 .AND. time .lt. 1000) then
          write(filename1,'(a6,f7.3,a4)') zonal,time,".dat"
       end if
       if (time .ge. 1000 .AND. time .lt. 10000) then
          write(filename1,'(a6,f8.3,a4)') zonal,time,".dat"
       end if
       if (time .ge. 10000 .AND. time .lt. 100000) then
          write(filename1,'(a6,f9.3,a4)') zonal,time,".dat"
       end if

   filename1 = trim(runname)//"_"//trim(filename1)

       open (unit=unitnumber, file= trim(filename1))
       write(unitnumber, '(a15,a15,a15)') '', 'xx' ,'phi', 'apar'
    end if

    workphi=0.0
    workapar=0.0

    if(iproc<npe) then
       m=1+iproc*nkx_par
       do i=1,nkx_par
          workphi(m)=phik(i)
          workapar(m)=ak(i)
          m=m+1
       end do
    end if
    call sum_reduce(workphi,0)
    call sum_reduce(workapar,0)
    if(proc0) then
       call oneDfourninv(workphi,phi,nlx)
       call oneDfourninv(workapar,apar,nlx)
       do i=1,nlx
          write(unitnumber,'(3g16.8)') xx(i),phi(i),apar(i)
       end do
       write(unitnumber,*) '                   '
    end if
      
  end subroutine zonalflows


!*********************************************
  subroutine proc0_recvs_and_writes(file1,file2,Apar,ne,g)
    !LF, 8/3/2015
    !write checkpoints: proc0 receives data from all others and writes to a file
    !Motivation: cut down on IO time

    use constants
    use mp,only:proc0,iproc
    use mpi_mod
    use HDF5
    implicit none
    !include 'mpif.h'

    integer:: ierror
    integer, dimension (MPI_STATUS_SIZE) :: status_mp  
    integer::i,j,k,n
    real, DIMENSION(nlx,nly_par,nlz_par):: Apar, ne
    !real, DIMENSION(nlx,nly_par,nlz_par):: epar
    real, DIMENSION(nlx,nly_par,nlz_par,gmin:ngtot):: g
    real:: trash
    character(len=100):: file1, file2
   ! --- hdf5 ---
    character(len=4) :: Apar_name = "Apar"
    character(len=2) :: ne_name = "ne"
    character(len=1) :: g_name = "g"
    integer(HID_T) :: file_id
    integer(HID_T) :: file_id_g
    integer(HSIZE_T), dimension(1) :: data_dims_total
    integer(HID_T) :: dataspace_Apar, dataspace_ne, dset_id_Apar, dset_id_ne, memspace
    integer(HSIZE_T), dimension(1) :: data_dims_single
    integer(HSIZE_T), dimension(1) :: data_dims_total_g
    integer(HID_T) :: dataspace_g, dset_id_g
    integer(HSIZE_T), dimension(1) :: data_dims_single_g

    integer(4) :: error
    integer(4) :: rank =1
    !integer(HSIZE_T), dimension(1) :: slab_size
    integer(HSIZE_T), dimension(1) :: offset
    integer(HSIZE_T), dimension(1) :: offset_g
    integer(HSIZE_T), dimension(1) :: stride = (/1/)
    integer(HSIZE_T), dimension(1) :: block_size = (/1/)

    real, allocatable, dimension(:,:,:) :: Apar_buff, ne_buff !epar_buff
    real, allocatable, dimension(:,:,:,:) :: gdummy
    
    allocate (Apar_buff(nlx,nly_par,nlz_par))
    allocate (ne_buff(nlx,nly_par,nlz_par))
    !allocate (epar_buff(nlx,nly_par,nlz_par))
    allocate (gdummy(nlx,nly_par,nlz_par,gmin:ngtot))

    trash = 0.0
    data_dims_single(1) = nlx*nly_par*nlz_par
    data_dims_total(1) = nlx*nly*nlz
    data_dims_single_g(1) = nlx*nly_par*nlz_par*(ngtot-gmin+1)
    data_dims_total_g(1) = nlx*nly*nlz*(ngtot-gmin+1)
    
    if (iproc==0) then
       !write(*,*) 'Hello world1'
      !  open (unit=16,file=trim(file1))
      !  open (unit=17,file=trim(file2))
       !write(*,*) 'Hello world2'
      !---------------hdf5-------------
       call h5open_f(error)
       call h5fcreate_f(trim(file1),H5F_ACC_TRUNC_F, file_id,error) ! create an HDF5 file
       call h5fcreate_f(trim(file2),H5F_ACC_TRUNC_F, file_id_g,error)

       call h5screate_simple_f(rank,data_dims_total,dataspace_Apar,error)
       call h5screate_simple_f(rank,data_dims_total,dataspace_ne,error)
       call h5screate_simple_f(rank,data_dims_total_g,dataspace_g,error)

       call h5dcreate_f(file_id, Apar_name, H5T_NATIVE_DOUBLE, dataspace_Apar, dset_id_Apar, error)
       call h5dcreate_f(file_id, ne_name, H5T_NATIVE_DOUBLE, dataspace_ne, dset_id_ne, error)
       call h5dcreate_f(file_id_g, g_name, H5T_NATIVE_DOUBLE, dataspace_g, dset_id_g, error)

       offset(1) = 0
       offset_g(1) = 0
      ! ---------------
       do k=1,nlz_par  !writes its own data
          do j=1,nly_par
             do i=1,nlx
                ! write(16,*) Apar(i,j,k), ne(i,j,k)
                ! write(16,*) Apar(i,j,k), ne(i,j,k), trash
                ! write(13,*) Apar(i,j,k), ne(i,j,k), epar(i,j,k)
                ! write(17,*) g(i,j,k,:)
               Apar_buff(i,j,k) = Apar(i,j,k)
               ne_buff(i,j,k) = ne(i,j,k)
               gdummy(i,j,k,:) = g(i,j,k,:)
             end do
          end do
       end do
       !write(*,*) 'Hello world3'

       !----hdf5-----
       call h5sselect_hyperslab_f(dataspace_Apar, H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
       call h5screate_simple_f(rank,data_dims_single,memspace,error)
       call h5dwrite_f(dset_id_Apar, H5T_NATIVE_DOUBLE, Apar_buff, data_dims_single, error, memspace, dataspace_Apar)
       call h5sclose_f(memspace, error)

       call h5sselect_hyperslab_f(dataspace_ne, H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
       call h5screate_simple_f(rank,data_dims_single,memspace,error)
       call h5dwrite_f(dset_id_ne, H5T_NATIVE_DOUBLE, ne_buff, data_dims_single, error, memspace, dataspace_ne)
       call h5sclose_f(memspace, error)

       call h5sselect_hyperslab_f(dataspace_g, H5S_SELECT_SET_F, offset_g, data_dims_single_g, error, stride, block_size)
       call h5screate_simple_f(rank,data_dims_single_g,memspace,error)
       call h5dwrite_f(dset_id_g, H5T_NATIVE_DOUBLE, gdummy,data_dims_single_g, error, memspace, dataspace_g)
       call h5sclose_f(memspace, error)

       do n=1,NPE*npez-1 
          call mpi_recv (Apar_buff, size(Apar_buff), MPI_DOUBLE_PRECISION, n, 1,& ! recv n-process's data
               MPI_COMM_WORLD, status_mp, ierror)
          call mpi_recv (ne_buff, size(ne_buff), MPI_DOUBLE_PRECISION, n, 2,&
               MPI_COMM_WORLD, status_mp, ierror)
          !call mpi_recv (epar_buff, size(epar_buff), MPI_DOUBLE_PRECISION, n, 2,&
          !             MPI_COMM_WORLD, status_mp, ierror)
          call mpi_recv (gdummy, size(gdummy), MPI_DOUBLE_PRECISION, n, 3,&
               MPI_COMM_WORLD, status_mp, ierror)
          !write(*,*) 'Hello world4-', n 
         !  do k=1,nlz_par  !write n-process's data
         !     do j=1,nly_par
         !        do i=1,nlx
         !           write(16,*) Apar_buff(i,j,k), ne_buff(i,j,k)
         !           !write(13,*) Apar_buff(i,j,k), ne_buff(i,j,k), epar(i,j,k)
         !           write(17,*) gdummy(i,j,k,:)
         !        end do
         !     end do
         !  end do

          offset(1) = n * data_dims_single(1)
          offset_g(1) = n * data_dims_single_g(1)

          call h5sselect_hyperslab_f(dataspace_Apar, H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
          call h5screate_simple_f(rank,data_dims_single,memspace,error)
          call h5dwrite_f(dset_id_Apar, H5T_NATIVE_DOUBLE, Apar_buff, data_dims_single,error, memspace, dataspace_Apar)
          call h5sclose_f(memspace,error)

          call h5sselect_hyperslab_f(dataspace_ne, H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
          call h5screate_simple_f(rank,data_dims_single,memspace,error)
          call h5dwrite_f(dset_id_ne, H5T_NATIVE_DOUBLE, ne_buff, data_dims_single, error, memspace, dataspace_ne)
          call h5sclose_f(memspace, error)
   
          call h5sselect_hyperslab_f(dataspace_g, H5S_SELECT_SET_F, offset_g, data_dims_single_g, error, stride, block_size)
          call h5screate_simple_f(rank,data_dims_single_g,memspace,error)
          call h5dwrite_f(dset_id_g, H5T_NATIVE_DOUBLE, gdummy,data_dims_single_g, error, memspace, dataspace_g)
          call h5sclose_f(memspace, error)
          !write(*,*) 'Hello world5-', n 
       end do
      !  close(16)
      !  close(17)
       call h5sclose_f(dataspace_Apar,error)
       call h5sclose_f(dataspace_ne,error)
       call h5sclose_f(dataspace_g,error)

       call h5dclose_f(dset_id_Apar, error)
       call h5dclose_f(dset_id_ne, error)
       call h5dclose_f(dset_id_g, error)

       call h5fclose_f(file_id, error)
       call h5fclose_f(file_id_g, error)
       call h5close_f(error)
    end if

    do n=1,NPE*npez-1
       if (iproc==n) then
          do k=1,nlz_par
             do j=1,nly_par
                do i=1,nlx
                   Apar_buff(i,j,k) = Apar(i,j,k)
                   ne_buff(i,j,k) = ne(i,j,k)
                   !epar_buff(i,j,k) = epar(i,j,k)
                   gdummy(i,j,k,:) = g(i,j,k,:)
                end do
             end do
          end do
          call mpi_send (Apar_buff, size(Apar), MPI_DOUBLE_PRECISION, 0, 1,&
               MPI_COMM_WORLD, ierror)
          call mpi_send (ne_buff, size(ne_buff), MPI_DOUBLE_PRECISION, 0, 2,&
               MPI_COMM_WORLD, ierror)
          !call mpi_send (epar_buff, size(epar_buff), MPI_DOUBLE_PRECISION, 0, 2,&
          !     MPI_COMM_WORLD, ierror)
          call mpi_send (gdummy, size(gdummy), MPI_DOUBLE_PRECISION, 0, 3,&
               MPI_COMM_WORLD, ierror)
       end if
    end do

    deallocate(Apar_buff)
    deallocate(ne_buff)
    !deallocate(epar_buff)
    deallocate(gdummy)

  end subroutine proc0_recvs_and_writes

!*********************************************
  subroutine proc0_recvs_and_writes_nogs(file1,file2,Apar,ne)
    !LMM: 28/10/2018 
    !write checkpoints: proc0 receives data from all others and writes to a file
    !Motivation: cut down on IO time

    use constants
    use mp,only:proc0,iproc
    use mpi_mod
    use HDF5
    implicit none
    integer:: ierror
    integer, dimension (MPI_STATUS_SIZE) :: status_mp  
    integer::i,j,k,n
    real, DIMENSION(nlx,nly_par,nlz_par):: Apar, ne
    real:: trash
    character(len=100):: file1, file2

    !hdf5 file and data
    character(len=4):: Apar_name = "Apar"
    character(len=2):: ne_name = "ne"
    integer(HID_T) :: file_id ! file identifier
    integer(HSIZE_T), dimension(1) :: data_dims_total ! size of entire domain dataset
    integer(HID_T):: dataspace_Apar, dataspace_ne, memspace, dset_id_Apar, dset_id_ne
    integer(4)     ::   error ! Error flag
    integer(4):: rank = 1  ! save data as 1d array -- should be possible to make this a 3d array, but do this case first.
    integer(HSIZE_T), dimension(1):: data_dims_single !size of single processor dataset
    INTEGER(HSIZE_T), DIMENSION(1) :: slab_size  ! Size of hyperslab = single processor dataset size
    INTEGER(HSIZE_T), DIMENSION(1) :: offset ! Hyperslab offset
    INTEGER(HSIZE_T), DIMENSION(1) :: stride = (/1/) ! Hyperslab stride 
    INTEGER(HSIZE_T), DIMENSION(1) :: block_size = (/1/)  ! Hyperslab block size 


    real, allocatable, dimension(:,:,:) :: Apar_buff, ne_buff
    allocate (Apar_buff(nlx,nly_par,nlz_par))
    allocate (ne_buff(nlx,nly_par,nlz_par))
    !print*, 'inside new subroutine'
    trash = 0.0

    data_dims_single(1) = nlx*nly_par*nlz_par
    data_dims_total(1) = nlx*nly*nlz

    if (iproc==0) then
       !open (unit=16,file=trim(file1))
       ! ---hdf5---
       call h5open_f(error)
       call h5fcreate_f(trim(file1), H5F_ACC_TRUNC_F, file_id, error)
       ! create dataspace for A and n
       call h5screate_simple_f(rank, data_dims_total, dataspace_Apar, error)
       call h5screate_simple_f(rank, data_dims_total, dataspace_ne, error)
       ! create datasets for A and n 
       call h5dcreate_f(file_id, Apar_name, H5T_NATIVE_DOUBLE, dataspace_Apar, dset_id_Apar, error)
       call h5dcreate_f(file_id, ne_name, H5T_NATIVE_DOUBLE, dataspace_ne, dset_id_ne, error)
       offset(1) = 0
       ! ----------
       do k=1,nlz_par  !writes its own data
          do j=1,nly_par
             do i=1,nlx
                !write(16,*) Apar(i,j,k), ne(i,j,k)
                Apar_buff(i,j,k) = Apar(i,j,k)
                ne_buff(i,j,k) = ne(i,j,k)
             end do
          end do
       end do
       ! select hyperslab for proc 0 for apar
       call h5sselect_hyperslab_f(dataspace_Apar,H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
       ! cretate memory dataspace for apar subset
       call h5screate_simple_f(rank, data_dims_single, memspace, error)
       ! write apar
       call h5dwrite_f(dset_id_Apar, H5T_NATIVE_DOUBLE, Apar_buff, data_dims_single, error, memspace, dataspace_Apar)
       ! close memspace (not sure if necessary, put here to be space) 
       call h5sclose_f(memspace, error)
       ! repeat for ne
       call h5sselect_hyperslab_f(dataspace_ne,H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
       ! cretate memory dataspace for ne
       call h5screate_simple_f(rank, data_dims_single, memspace, error)
       ! write ne
       call h5dwrite_f(dset_id_ne, H5T_NATIVE_DOUBLE, ne_buff, data_dims_single, error, memspace, dataspace_ne)
       ! close ne
       call h5sclose_f(memspace, error)

       do n=1,NPE*npez-1 
          call mpi_recv (Apar_buff, size(Apar_buff), MPI_DOUBLE_PRECISION, n, 1,& ! recv n-process's data
               MPI_COMM_WORLD, status_mp, ierror)
          call mpi_recv (ne_buff, size(ne_buff), MPI_DOUBLE_PRECISION, n, 2,&
               MPI_COMM_WORLD, status_mp, ierror)
         !  do k=1,nlz_par  !write n-process's data
         !     do j=1,nly_par
         !        do i=1,nlx
         !           write(16,*) Apar_buff(i,j,k), ne_buff(i,j,k)
         !        end do
         !     end do
         !  end do 
         !---hdf5---
          offset(1) = n*data_dims_single(1)
          call h5sselect_hyperslab_f(dataspace_Apar,H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
          call h5screate_simple_f(rank, data_dims_single, memspace, error)
          ! write apar
          call h5dwrite_f(dset_id_Apar, H5T_NATIVE_DOUBLE, Apar_buff, data_dims_single, error, memspace, dataspace_Apar)
          ! close memspace (not sure if necessary, put here to be space) 
          call h5sclose_f(memspace, error)
          ! repeat for ne
          call h5sselect_hyperslab_f(dataspace_ne,H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
          ! cretate memory dataspace for ne
          call h5screate_simple_f(rank, data_dims_single, memspace, error)
       ! write ne
          call h5dwrite_f(dset_id_ne, H5T_NATIVE_DOUBLE, ne_buff, data_dims_single, error, memspace, dataspace_ne)
          ! close ne
          call h5sclose_f(memspace, error)
       end do
       !close(16)
       call h5sclose_f(dataspace_Apar,error)
       call h5sclose_f(dataspace_ne,error)
       call h5dclose_f(dset_id_Apar, error)
       call h5dclose_f(dset_id_ne, error)
       call h5fclose_f(file_id, error)
       call h5close_f(error)
    end if

    do n=1,NPE*npez-1
       if (iproc==n) then
          do k=1,nlz_par
             do j=1,nly_par
                do i=1,nlx
                   Apar_buff(i,j,k) = Apar(i,j,k)
                   ne_buff(i,j,k) = ne(i,j,k)
                end do
             end do
          end do
          call mpi_send (Apar_buff, size(Apar), MPI_DOUBLE_PRECISION, 0, 1,&
               MPI_COMM_WORLD, ierror)
          call mpi_send (ne_buff, size(ne_buff), MPI_DOUBLE_PRECISION, 0, 2,&
               MPI_COMM_WORLD, ierror)
       end if
    end do

    deallocate(Apar_buff)
    deallocate(ne_buff)

  end subroutine proc0_recvs_and_writes_nogs


!*********************************************



  subroutine proc0_reads_to_all(restart1,&
       restart2,Apar,Apar_buff,ne,ne_buff,g,gdummy)
    !LF, 7/3/2015
    !restart process: proc0 reads the files and sends to others, sequentially
    !Motivation: try to circunvent having files of many gigabytes read by many
    !processes

    use constants
    use mp,only:proc0,iproc,send,receive
    use mpi_mod
    use HDF5
    implicit none
    !include 'mpif.h'

   ! Z. Liu, 4/8/2020
    integer:: ierror
    integer, dimension (MPI_STATUS_SIZE) :: status_mp

    integer::i,j,k,n
    real, DIMENSION(nlx,nly_par,nlz_par):: Apar, Apar_buff, ne, ne_buff
    real, DIMENSION(nlx,nly_par,nlz_par,gmin:ngtot):: g, gdummy
    real:: trash
    character(len=100):: restart1, restart2

    ! hdf5 file and data properties, zliu 04/06/2020
    character(len=4) :: Apar_name = "Apar"
    character(len=2) :: ne_name = "ne"
    character(len=1) :: g_name = "g"

    integer(HID_T) :: file_id ! file identifier1
    integer(HID_T) :: file_id_g ! file identifier for g file (restart2)
    integer(HSIZE_T), dimension(1) :: data_dims_total ! size of entire domain dataset
    integer(HID_T) :: dataspace_Apar, dataspace_ne, memspace, dset_id_Apar, dset_id_ne !???
    integer(HID_T) :: dataspace_g, dset_id_g ! dset_id means dataset identifier
    integer(4) :: error
    integer(4) :: rank = 1 ! save data as 1d array    
    !! hyperslab is portions of datasets, here it might be refer to a portion of dataset in ONE PROCESSOR
    integer(HSIZE_T), dimension(1) :: data_dims_single !size of single processor dataset
    integer(HSIZE_T), dimension(1) :: slab_size !? = single processor dataset size
    !! stride and block_size are parameters when selecting hyperslab
    integer(HSIZE_T), dimension(1) :: stride = (/1/)
    integer(HSIZE_T), dimension(1) :: block_size = (/1/)
    integer(HSIZE_T), dimension(1) :: offset
    integer(HSIZE_T), dimension(1) :: offset_g
    ! two dimensions but for g
    integer(HSIZE_T), dimension(1) :: data_dims_single_g
    integer(HSIZE_T), dimension(1) :: data_dims_total_g

    real, allocatable, dimension(:,:,:) :: Apar_buff2, ne_buff2
    real, allocatable, dimension(:,:,:,:):: gdummy2

    allocate (Apar_buff2(nlx,nly_par,nlz_par))
    allocate (ne_buff2(nlx,nly_par,nlz_par))
    allocate (gdummy2(nlx,nly_par,nlz_par,gmin:ngtot))    
    
    !hdf5 data dimensions
    data_dims_total(1) = nlx*nly*nlz
    data_dims_single(1) = nlx*nly_par*nlz_par
    ! two dimensions but for g
    data_dims_single_g(1) = nlx*nly_par*nlz_par*(ngtot-gmin+1)
    data_dims_total_g(1) = nlx*nly*nlz*(ngtot-gmin+1)
    
    do n=1,NPE*npez-1
         if (iproc==n) then
           call mpi_recv (Apar_buff, size(Apar_buff), MPI_DOUBLE_PRECISION, 0,1,&
                        MPI_COMM_WORLD, status_mp, ierror)
           call mpi_recv (ne_buff, size(ne_buff), MPI_DOUBLE_PRECISION, 0, 2,&
                        MPI_COMM_WORLD, status_mp, ierror)
           !call mpi_recv (epar_buff, size(epar_buff), MPI_DOUBLE_PRECISION, 0, 2,&
           !             MPI_COMM_WORLD, status_mp, ierror)
           call mpi_recv (gdummy, size(gdummy), MPI_DOUBLE_PRECISION, 0, 3,&
                        MPI_COMM_WORLD, status_mp, ierror)


           do k=1,nlz_par
             do j=1,nly_par
               do i=1,nlx
                 Apar(i,j,k) = Apar_buff(i,j,k)
                 ne(i,j,k) = ne_buff(i,j,k)
                 !epar(i,j,k) = epar_buff(i,j,k)
                 g(i,j,k,:) = gdummy(i,j,k,:)
               end do
             end do
           end do
         end if
    end do

    if (iproc==0) then
      !   open (unit=12,file=trim(restart2),status='OLD')
      !   open (unit=13,file=trim(restart1),status='OLD')

      !!!!!!!!!!! read hdf5 files !!!!!!!!!!!!!!!!!!!!
      offset(1) = 0
      offset_g(1) = 0
      !open file
      call h5open_f(error) ! initialize the library
      call h5fopen_f(trim(restart1),H5F_ACC_RDONLY_F,file_id,error) ! ACC_RDONLY means READ_ONLY
      call h5fopen_f(trim(restart2),H5F_ACC_RDONLY_F,file_id_g,error)
      !open an existing dataset, returns an identifier for the dataset
      call h5dopen_f(file_id, Apar_name, dset_id_Apar, error)
      call h5dopen_f(file_id, ne_name, dset_id_ne,error)
      call h5dopen_f(file_id_g, g_name,dset_id_g,error)
      !make a copy of dataspace of the dataset, returns an identifier for a copy of the dataspace of the dataset
      !! dataspace id a part of metadata of a dataset, it describes the layout of a dataset's data elements
      !! two roll of a dataspce:
      !!! 1. contains spatial information (rand and dimensions)
      !!! 2. it can be used to select a portion or subset of a dataset by describing data buffers participating in I/O
      call h5dget_space_f(dset_id_Apar, dataspace_Apar, error)
      call h5dget_space_f(dset_id_ne, dataspace_ne, error)
      call h5dget_space_f(dset_id_g, dataspace_g, error)

      !Apar initial read for proc0
      call h5sselect_hyperslab_f(dataspace_Apar, H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
      call h5screate_simple_f(rank, data_dims_single, memspace, error)
      call h5dread_f(dset_id_Apar, H5T_NATIVE_DOUBLE, Apar_buff2, data_dims_single, error, memspace, dataspace_Apar)
      call h5sclose_f(memspace, error)
      !ne initial read for proc0
      call h5sselect_hyperslab_f(dataspace_ne, H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
      call h5screate_simple_f(rank, data_dims_single, memspace, error)
      call h5dread_f(dset_id_ne, H5T_NATIVE_DOUBLE, ne_buff2, data_dims_single, error, memspace, dataspace_ne)
      call h5sclose_f(memspace, error)

      !g initial read for proc0
      !! select a hyperslab region to add to the current selected region
      call h5sselect_hyperslab_f(dataspace_g, H5S_SELECT_SET_F, offset_g, data_dims_single_g, error, stride, block_size)
      !! creates a new dataspace and opens it for access
      !!! rank is the number of dimension, data_dims is the current dimension
      !!! returns a dataspace identifier (memory dataspace) ("memspace" here)
      call h5screate_simple_f(rank, data_dims_single_g, memspace, error)
      !! reads raw dataset specified by data_id from the file into an application memory buffer
      !!! NATIVE_DOUBLE is a memtype id, gdummy2 is buffer
      !!! The part of the dataset to read is defined by memspace_id and filespace_id
      !!!! memspace_id specifies the momory dataspace and the selection within it, filespace_id specifies the selection within the file dataset's dataspace
      !!!! File is refered to original dataset, it has data_id and dataspace_id, memspace is the target
      !!!! I think dataspace_g here is actually H5S_ALL ? NO, bcause HYPER_SELECTION!
      call h5dread_f(dset_id_g, H5T_NATIVE_DOUBLE, gdummy2, data_dims_single_g, error, memspace, dataspace_g)
      call h5sclose_f(memspace, error)


         do k=1,nlz_par  !reads its own data
            do j=1,nly_par
               do i=1,nlx
                  Apar(i,j,k) = Apar_buff2(i,j,k)
                  ne(i,j,k) = ne_buff2(i,j,k)
                  g(i,j,k,:) = gdummy2(i,j,k,:)
                 !read(13,*) Apar(i,j,k), ne(i,j,k)
                 !read(13,*) Apar(i,j,k), ne(i,j,k), epar(i,j,k)
                 !read(12,*) g(i,j,k,:)
               end do
            end do 
         end do

         do n=1,NPE*npez-1
            offset(1) = n*data_dims_single(1)
            offset_g(1) = n*data_dims_single_g(1)
            ! Apar read for proc n
            call h5sselect_hyperslab_f(dataspace_Apar, H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
            call h5screate_simple_f(rank, data_dims_single, memspace, error)
            call h5dread_f(dset_id_Apar, H5T_NATIVE_DOUBLE, Apar_buff2, data_dims_single, error, memspace, dataspace_Apar)
            call h5sclose_f(memspace, error)
            ! ne read
            call h5sselect_hyperslab_f(dataspace_ne, H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
            call h5screate_simple_f(rank, data_dims_single, memspace, error)
            call h5dread_f(dset_id_ne, H5T_NATIVE_DOUBLE, ne_buff2, data_dims_single, error, memspace,dataspace_ne)
            call h5sclose_f(memspace, error)
            ! g read
            call h5sselect_hyperslab_f(dataspace_g, H5S_SELECT_SET_F, offset_g, data_dims_single_g, error, stride, block_size)
            call h5screate_simple_f(rank, data_dims_single_g, memspace, error)
            call h5dread_f(dset_id_g, H5T_NATIVE_DOUBLE, gdummy2, data_dims_single_g, error, memspace, dataspace_g)
            call h5sclose_f(memspace, error)


         !  do k=1,nlz_par  !reads n-process's data
         !     do j=1,nly_par
         !       do i=1,nlx
         !         read(13,*) Apar_buff(i,j,k), ne_buff(i,j,k)
         !         !read(13,*) Apar_buff(i,j,k), ne_buff(i,j,k), epar(i,j,k)
         !         read(12,*) gdummy(i,j,k,:)
         !       end do
         !     end do
         !  end do

         !  call mpi_send (Apar_buff, size(Apar), MPI_DOUBLE_PRECISION, n, 1,&
         !          MPI_COMM_WORLD, ierror)
         !  call mpi_send (ne_buff, size(ne_buff), MPI_DOUBLE_PRECISION, n, 2,&
         !         MPI_COMM_WORLD, ierror)
         !  !call mpi_send (epar_buff, size(epar_buff), MPI_DOUBLE_PRECISION, n, 2,&
         !  !       MPI_COMM_WORLD, ierror)
         !  call mpi_send (gdummy, size(gdummy), MPI_DOUBLE_PRECISION, n, 3,&
         !         MPI_COMM_WORLD, ierror)
            call mpi_send(Apar_buff2, size(Apar_buff2), MPI_DOUBLE_PRECISION,n,1,&
                   MPI_COMM_WORLD, ierror)
            call mpi_send(ne_buff2, size(ne_buff2), MPI_DOUBLE_PRECISION,n,2,&
                   MPI_COMM_WORLD, ierror)
            call mpi_send(gdummy2, size(gdummy2), MPI_DOUBLE_PRECISION,n,3,&
                   MPI_COMM_WORLD, ierror)
         end do

         ! hdf5 close
         call h5sclose_f(dataspace_Apar,error)
         call h5sclose_f(dataspace_ne, error)
         call h5sclose_f(dataspace_g, error)
         call h5dclose_f(dset_id_Apar,error)
         call h5dclose_f(dset_id_ne,error)
         call h5dclose_f(dset_id_g,error)
         call h5fclose_f(file_id,error)
         call h5fclose_f(file_id_g,error)
         call h5close_f(error)
   end if

 end subroutine proc0_reads_to_all

!*********************************************
  subroutine proc0_reads_to_all_nogs(restart1,&
       restart2,Apar,Apar_buff,ne,ne_buff)
    !LF, 7/3/2015
    !restart process: proc0 reads the files and sends to others, sequentially
    !Motivation: try to circunvent having files of many gigabytes read by many
    !processes

    use constants
    use mp,only:proc0,iproc,send,receive
    use mpi_mod
    use HDF5
    implicit none
    !include 'mpif.h'

    integer:: ierror
    integer, dimension (MPI_STATUS_SIZE) :: status_mp

    integer::i,j,k,n
    real, DIMENSION(nlx,nly_par,nlz_par):: Apar, Apar_buff, ne, ne_buff, epar, epar_buff
    real:: trash
    character(len=100):: restart1, restart2

    ! hdf5 file and data propertoes

    character(len=4):: Apar_name = "Apar"
    character(len=2):: ne_name = "ne"
    integer(HID_T) :: file_id ! file identifier
    integer(HSIZE_T), dimension(1) :: data_dims_total ! size of entire domain dataset
    integer(HID_T):: dataspace_Apar, dataspace_ne, memspace, dset_id_Apar, dset_id_ne
    integer(4)     ::   error ! Error flag
    integer(4):: rank = 1  ! save data as 1d array -- should be possible to make this a 3d array, but do this case first.
    integer(HSIZE_T), dimension(1):: data_dims_single !size of single processor dataset
    !! hyperslab is portions of datasets, here it might be refer to a portion of dataset in ONE PROCESSOR
    INTEGER(HSIZE_T), DIMENSION(1) :: slab_size  ! Size of hyperslab = single processor dataset size
    INTEGER(HSIZE_T), DIMENSION(1) :: offset ! Hyperslab offset
    INTEGER(HSIZE_T), DIMENSION(1) :: stride = (/1/) ! Hyperslab stride 
    INTEGER(HSIZE_T), DIMENSION(1) :: block_size = (/1/)  ! Hyperslab block size 

    real, allocatable, dimension(:,:,:) :: Apar_buff2, ne_buff2 !epar_buff
    allocate (Apar_buff2(nlx,nly_par,nlz_par))
    allocate (ne_buff2(nlx,nly_par,nlz_par))
    trash = 0.0
   ! hdf5 data dimensions
    data_dims_single(1) = nlx*nly_par*nlz_par
    data_dims_total(1) = nlx*nly*nlz

    !print*, 'loading within nogs'
    do n=1,NPE*npez-1
         if (iproc==n) then
           call mpi_recv (Apar_buff, size(Apar_buff), MPI_DOUBLE_PRECISION, 0,1,&
                        MPI_COMM_WORLD, status_mp, ierror)
           call mpi_recv (ne_buff, size(ne_buff), MPI_DOUBLE_PRECISION, 0, 2,&
                        MPI_COMM_WORLD, status_mp, ierror)
           !call mpi_recv (epar_buff, size(epar_buff), MPI_DOUBLE_PRECISION, 0, 2,&
           !             MPI_COMM_WORLD, status_mp, ierror)
           !call mpi_recv (gdummy, size(gdummy), MPI_DOUBLE_PRECISION, 0, 3,&
           !             MPI_COMM_WORLD, status_mp, ierror)


           do k=1,nlz_par
             do j=1,nly_par
               do i=1,nlx
                 Apar(i,j,k) = Apar_buff(i,j,k)
                 ne(i,j,k) = ne_buff(i,j,k)
                 !epar(i,j,k) = epar_buff(i,j,k)
                 !g(i,j,k,:) = gdummy(i,j,k,:)
               end do
             end do
           end do
         end if
    end do
    !print*, 'loading2'
    if (iproc==0) then
      !  ! open (unit=12,file=trim(restart2),status='OLD')
      !   open (unit=13,file=trim(restart1),status='OLD')
      !    do k=1,nlz_par  !reads its own data
      !       do j=1,nly_par
      !          do i=1,nlx
      !             read(13,*) Apar(i,j,k), ne(i,j,k)
      !            !read(13,*) Apar(i,j,k), ne(i,j,k), epar(i,j,k)
      !            !read(12,*) g(i,j,k,:)
      !          end do
      !       end do
      !    end do
      !    !print*, 'getting there'
      !    do n=1,NPE*npez-1
      !     do k=1,nlz_par  !reads n-process's data
      !        do j=1,nly_par
      !          do i=1,nlx
      !            read(13,*) Apar_buff(i,j,k), ne_buff(i,j,k)
      !            !read(13,*) Apar_buff(i,j,k), ne_buff(i,j,k), epar_buff(i,j,k)
      !            !read(12,*) gdummy(i,j,k,:)
      !          end do
      !        end do
      !     end do
      !     !print*, 'almost done'
      !     call mpi_send (Apar_buff, size(Apar), MPI_DOUBLE_PRECISION, n, 1,&
      !             MPI_COMM_WORLD, ierror)
      !     call mpi_send (ne_buff, size(ne_buff), MPI_DOUBLE_PRECISION, n, 2,&
      !            MPI_COMM_WORLD, ierror)
      !     !call mpi_send (epar_buff, size(epar_buff), MPI_DOUBLE_PRECISION, n, 2,&
      !     !       MPI_COMM_WORLD, ierror)
      !     !call mpi_send (gdummy, size(gdummy), MPI_DOUBLE_PRECISION, n, 3,&
      !     !       MPI_COMM_WORLD, ierror)
      !    end do

       ! open (unit=12,file=trim(restart2),status='OLD')
       !        open (unit=13,file=trim(restart1),status='OLD')
       !       read hdf5 file
       offset(1) = 0
       !open file
       call h5open_f(error)
       call h5fopen_f(trim(restart1), H5F_ACC_RDONLY_F, file_id, error)
       !open dataset
       call h5dopen_f(file_id, Apar_name, dset_id_Apar, error)
       call h5dopen_f(file_id, ne_name, dset_id_ne, error)
       !get dataspace
       call h5dget_space_f(dset_id_Apar, dataspace_Apar, error)
       call h5dget_space_f(dset_id_ne, dataspace_ne, error)
       ! Apar initial read for  proc0
       call h5sselect_hyperslab_f(dataspace_Apar, H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
       call h5screate_simple_f(rank, data_dims_single, memspace, error)
       call h5dread_f(dset_id_Apar, H5T_NATIVE_DOUBLE, Apar_buff2, data_dims_single, error, memspace, dataspace_Apar)
       call h5sclose_f(memspace, error)
       ! ne initial read
       call h5sselect_hyperslab_f(dataspace_ne, H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
       call h5screate_simple_f(rank, data_dims_single, memspace, error)
       call h5dread_f(dset_id_ne, H5T_NATIVE_DOUBLE, ne_buff2, data_dims_single, error, memspace, dataspace_ne)
       call h5sclose_f(memspace, error)
      
         do k=1,nlz_par  !reads its own data
            do j=1,nly_par
               do i=1,nlx
                  Apar(i,j,k) = Apar_buff2(i,j,k)
                  ne(i,j,k) = ne_buff2(i,j,k)
!                  read(13,*) Apar(i,j,k), ne(i,j,k)
                !read(13,*) Apar(i,j,k), ne(i,j,k), epar(i,j,k)
                !read(12,*) g(i,j,k,:)
               end do
            end do
         end do
        !print*, 'getting there'
         do n=1,NPE*npez-1
            offset(1) = n*data_dims_single(1)
            ! Apar read for proc n
            call h5sselect_hyperslab_f(dataspace_Apar, H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
            call h5screate_simple_f(rank, data_dims_single, memspace, error)
            call h5dread_f(dset_id_Apar, H5T_NATIVE_DOUBLE, Apar_buff2, data_dims_single, error, memspace, dataspace_Apar)
            call h5sclose_f(memspace, error)
            ! ne read
            call h5sselect_hyperslab_f(dataspace_ne, H5S_SELECT_SET_F, offset, data_dims_single, error, stride, block_size)
            call h5screate_simple_f(rank, data_dims_single, memspace, error)
            call h5dread_f(dset_id_ne, H5T_NATIVE_DOUBLE, ne_buff2, data_dims_single, error, memspace,dataspace_ne)
            call h5sclose_f(memspace, error)
           !          do k=1,nlz_par  !reads n-process's data
!            do j=1,nly_par
 !             do i=1,nlx
  !              read(13,*) Apar_buff(i,j,k), ne_buff(i,j,k)
                !read(13,*) Apar_buff(i,j,k), ne_buff(i,j,k), epar_buff(i,j,k)
                !read(12,*) gdummy(i,j,k,:)
   !           end do
    !        end do
     !    end do
         !print*, 'almost done'
          call mpi_send (Apar_buff2, size(Apar_buff2), MPI_DOUBLE_PRECISION, n, 1,&
                  MPI_COMM_WORLD, ierror)
          call mpi_send (ne_buff2, size(ne_buff2), MPI_DOUBLE_PRECISION, n, 2,&
                 MPI_COMM_WORLD, ierror)
         !call mpi_send (epar_buff, size(epar_buff), MPI_DOUBLE_PRECISION, n, 2,&
         !       MPI_COMM_WORLD, ierror)
         !call mpi_send (gdummy, size(gdummy), MPI_DOUBLE_PRECISION, n, 3,&
         !       MPI_COMM_WORLD, ierror)
      end do

      ! close hdf5 stuff
      call h5sclose_f(dataspace_Apar,error)
      call h5sclose_f(dataspace_ne,error)
      call h5dclose_f(dset_id_Apar, error)
      call h5dclose_f(dset_id_ne, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)      
   end if

 end subroutine proc0_reads_to_all_nogs

end module diag


