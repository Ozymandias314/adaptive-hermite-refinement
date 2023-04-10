module forcing

  implicit none

contains

  subroutine  force(kfz,dtau,kfp1,kfp2,feps,fieldk)
    !Tarek's forcing routine, edited by NFL
    ! subroutine  force(kfz,feps,dtau,field)
    ! White noise in time forcing in the range between kp1 and kp2 
    ! set in run.in.  Work (per volume) done is feps.  dtau is the size if time
    ! step. zp is forced if field='p', zm if field='m', or both zp and zm  if
    ! field='b'.
    ! use global       ! kfp2, kfp1, nx,ny,  nky
    
    use mp,only:iproc,proc0
    use constants, only: pi, nlx, nly, nkx_par, nky, nlz_par, lx, ly, lz, npe, rhoi, rhos
    use grid	
    !  use mp,only:proc0
    implicit none
    
    real, intent(in) :: dtau, kfz, kfp1,kfp2, feps
    complex, dimension (nky,nkx_par,nlz_par)::fieldk
    real, dimension(nky,nkx_par,nlz_par)::re_fieldk,im_fieldk
    integer, dimension(:,:) , allocatable :: kfs_big
    integer, dimension(:,:) , allocatable, save :: kfs
    real :: kfp, kav=0
    real :: amp, phi , phiz
    integer ::  i,j,k,iloc,ik,nk=0,sgn,proc_force !,ip,i_rl,i_im
    logical, save :: lfirst=.true.

    fieldk=0.0

    ! First time calculate possible k_perp for forcing
    ! kfp1 <|k_perp|< kfp2  and |kperp| .ne. 0
    ! kfp1 and kfp2 are set in run.in
    
    if (lfirst) then  
       allocate(kfs_big(nkx_par*nky,2))
       do i=-ceiling(kfp2),ceiling(kfp2)    !because kx=0....nkx; ky=-nky/2-1...nky/2
          do j=-ceiling(kfp2),ceiling(kfp2)
             kfp=sqrt(i**2.+j**2.)
             ! can not force zero mode
             !.ttr. this should be checked:
             !....  if ( (kfp .ge. kfp1) .and. (kfp .le. kfp2) .and. (abs(kfp) .ge. 1.e-40) ) then 
             if ((kfp.ge.kfp1).and.(kfp.le.kfp2).and.(kfp.ne.0.)) then 
             
                nk=nk+1
                kfs_big(nk,1)=i
                kfs_big(nk,2)=j
                kav=kav+kfp
             endif
          enddo
       enddo
       kav=kav/nk
       allocate(kfs(nk,2))
       kfs=kfs_big(1:nk,1:2)
       deallocate(kfs_big)
       !   if ((ip.ge.10) .and.proc0) then
       if(proc0) then
          write(*,'(A,I3,G10.3,A)',advance='no') 'Forcing:,nk, kav =', nk, kav,' modes:'
          do i=1,nk 
             write(*,'(A,I3,A,I3,A)',advance='no') '(', kfs(i,1),',',kfs(i,2),')'
          enddo
          write(*,*)
       endif
       
       lfirst=.false. 
       
    endif
    
    ! Choose a random element from kfs
    ik=nk*.9999*uniran()+1
    
    kfp=sqrt(kfs(ik,1)**2d0+kfs(ik,2)**2d0)
    !  kfp=sqrt(kx(kfs(ik,1))**2+ky(kfs(ik,2))**2)
    
    ! amplitude and phase random phase 
    amp=1.0/abs(kfp)*sqrt(-feps/dtau*log(uniran()))*sqrt(nlx*nly*1.0)
    
    !NFL: this amp depends on the field that's actually being forced;  
    
    phi=pi*(2d0*uniran()-1d0)
    
    ! if kx < 0 then use u(-\vec k)=u^*(\vec k)
    !NFL: now we convert the above i,j to array indices:
    sgn=sign(1,kfs(ik,1))
    !  i_rl=2*sgn*kfs(ik,1)+1
    !  i_im=2*sgn*kfs(ik,1)+2
    i=sgn*kfs(ik,1)+1
    !  j=sgn*kfs(ik,2)+1
	! AVK, 17/03/2013: where did the factor of ly/lx come from?
    j=sgn*kfs(ik,2)*ly/lx+1
!NFL, 11/03/2013: the following line may be wrong
!AVK, 17/03/2013: it is correct if ly=lx
    if (j.le.0) j=j+nky   
    !identify the proc num that contains the k to be forced:
    !alfred: this will need editing with new parallelization.
    proc_force=(i-1)/nkx_par
    iloc=i-proc_force*nkx_par
    
    phi=sgn*phi
    
    !.ttr. this should be checked:
    !...   if (abs(kfz) .ge. 1.e-40) then 
    if (kfz .ne. 0.) then 
       ! normalize
       amp=amp*sqrt(2d0)
       ! random phase
       phiz=pi*(2d0*uniran()-1d0) 
    else
       phiz=0d0
    endif
    
    !  if (ip.ge.20) print*,'# Forcing ',field,' kfx,kfy,kfz,|,amp,phi,phiz', &
    !              kfs(ik,1),kfs(ik,2),kfz,'|',amp,phi,phiz 
    
    !  ! force zp
    !  if ((field.eq.'b').or.(field.eq.'p')) then
    !     Sp(i_rl,j,:)=Sp(i_rl,j,:)+amp*cos(phi)*cos(kfz*z+phiz)
    !     Sp(i_im,j,:)=Sp(i_im,j,:)+amp*sin(phi)*cos(kfz*z+phiz)
    !  endif 
    
    !  ! force zm
    !  if ((field.eq.'b').or.(field.eq.'m')) then
    !     Sm(i_rl,j,:)=Sm(i_rl,j,:)+amp*cos(phi)*cos(kfz*z+phiz)
    !     Sm(i_im,j,:)=Sm(i_im,j,:)+amp*sin(phi)*cos(kfz*z+phiz)
    !  endif
    re_fieldk(:,:,:)=REAL(fieldk(:,:,:))  
    im_fieldk(:,:,:)=AIMAG(fieldk(:,:,:))
    
    !  re_fieldk(j,i)=re_fieldk(j,i)-kperp(j,i)**2*amp*cos(phi)                              
    !  im_fieldk(j,i)=im_fieldk(j,i)-kperp(j,i)**2*amp*sin(phi)

!    do ip=0, NPE*npez-1
    if (mod(iproc,NPE)==proc_force) then
       do k=1, nlz_par
          !re_fieldk(j,iloc,k)=re_fieldk(j,iloc,k)-kfp**2*amp*cos(phi)*cos(2*pi*kfz*zz(k)/lz+phiz)   
          !im_fieldk(j,iloc,k)=im_fieldk(j,iloc,k)-kfp**2*amp*sin(phi)*cos(2*pi*kfz*zz(k)/lz+phiz)
          !re_fieldk(j,iloc,k)=re_fieldk(j,iloc,k)+2./rhoi**2*(gama0(kfp**2*rhoi**2/2.)-1.)*amp*cos(phi)*cos(2*pi*kfz*zz(k)/lz+phiz)
          !im_fieldk(j,iloc,k)=im_fieldk(j,iloc,k)+2./rhoi**2*(gama0(kfp**2*rhoi**2/2.)-1.)*amp*sin(phi)*cos(2*pi*kfz*zz(k)/lz+phiz)
          !Z.Liu 7/2/2020 
          !The modified normalization factor does NOT mean energy injected is independent of rho_i, it is just logically consistent with inplementations above 
          !re_fieldk(j,iloc,k)=re_fieldk(j,iloc,k)+1./(rhos**2-1./(2./rhoi**2*(gama0(kfp**2*rhoi**2/2.)-1.)))*amp*cos(phi)*cos(2*pi*kfz*zz(k)/lz+phiz)
          !im_fieldk(j,iloc,k)=im_fieldk(j,iloc,k)+1./(rhos**2-1./(2./rhoi**2*(gama0(kfp**2*rhoi**2/2.)-1.)))*amp*sin(phi)*cos(2*pi*kfz*zz(k)/lz+phiz)   
          !Z.Liu 7/3/2020 
          !"sqrt" added here is to seek a way which may make injected energy similar for different rhoi. Adding "sqrt" assumes that not only S, but ne will also be divided by this factor.
          re_fieldk(j,iloc,k)=re_fieldk(j,iloc,k)+1./sqrt(rhos**2-1./(2./rhoi**2*(gama0(kfp**2*rhoi**2/2.)-1.)))*amp*cos(phi)*cos(2*pi*kfz*zz(k)/lz+phiz)
          im_fieldk(j,iloc,k)=im_fieldk(j,iloc,k)+1./sqrt(rhos**2-1./(2./rhoi**2*(gama0(kfp**2*rhoi**2/2.)-1.)))*amp*sin(phi)*cos(2*pi*kfz*zz(k)/lz+phiz) 
       enddo
       !NFL: note: if FLR effects are included, then the factor is not kperp^2,
       !but something involving Gamma_0(b)
       
       if (i .eq. 1) then 
          ! ensure that f*(k)=f(-k) on Im axis 
          !    fieldk(nky-j+2,i)=fieldk(j,i)

!NFL, 11/03/2013: the followings line may be wrong if j=1 (see above)
          re_fieldk(nky-j+2,i,:)=re_fieldk(j,i,:)
          im_fieldk(nky-j+2,i,:)=-im_fieldk(j,i,:)
          
          !    Sp(i_rl,nky-j+2,:)= Sp(i_rl,j,:)
          !    Sp(i_im,nky-j+2,:)=-Sp(i_im,j,:)
          !    Sm(i_rl,nky-j+2,:)= Sm(i_rl,j,:)
          !    Sm(i_im,nky-j+2,:)=-Sm(i_im,j,:)
       endif
    end if
    fieldk(:,:,:)=re_fieldk(:,:,:)+(0.0,1.0)*im_fieldk(:,:,:)
 !end do
  end subroutine force


  function uniran(init)
  !
  ! 26-sep-02/wolf: Adapted from `Numerical Recipes for F90' ran() routine
  ! 23-mai-06/tay: Ripped from Pencil-code
  !
  ! "Minimal" random number generator of Park and Miller combined
  ! with a Marsaglia shift sequence. Returns a uniform random deviate
  ! between 0.0 and 1.0 (exclusive of the endpoint values).
  ! Call with (INIT=ival) to initialize.
  ! The period of this generator is supposed to be about 3.1� 10^18.
  !
  implicit none
  !
  real :: uniran
  integer, parameter :: mseed=256
 ! integer, dimension(mseed), save :: seed=0
  integer, dimension(mseed), save :: rstate=0
  real, parameter :: impossible=3.9085e37
  real, save :: am=impossible  ! will be constant on a given platform
  integer, optional, intent(in) :: init
  integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
  integer :: k,init_ts=1812   ! default value
  logical, save :: first_call=.true.

  !ajw This doesn't appear to always get set!
  if (first_call) then
    am=nearest(1.0,-1.0)/im
    first_call=.false.
  endif
  if (present(init) .or. rstate(1)==0 .or. rstate(2)<=0) then
    !
    ! initialize
    !
    if (present(init)) init_ts = init
    am=nearest(1.0,-1.0)/im
    rstate(1)=ieor(777755555,abs(init_ts))
    rstate(2)=ior(ieor(888889999,abs(init_ts)),1)
  endif
  !
  ! Marsaglia shift sequence with period 2^32-1
  !
  rstate(1)=ieor(rstate(1),ishft(rstate(1),13))
  rstate(1)=ieor(rstate(1),ishft(rstate(1),-17))
  rstate(1)=ieor(rstate(1),ishft(rstate(1),5)) 
  !
  ! Park-Miller sequence by Schrage's method, period 2^31-2
  !
  k=rstate(2)/iq
  rstate(2)=ia*(rstate(2)-k*iq)-ir*k
  if (rstate(2) < 0) rstate(2)=rstate(2)+im
  !
  ! combine the two generators with masking to ensure nonzero value
  !
  uniran=am*ior(iand(im,ieor(rstate(1),rstate(2))),1)
  endfunction uniran

  FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
           NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
   END function ran1

end module forcing
