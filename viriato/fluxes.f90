module Fluxes

  implicit none

contains


!*********************************************************  
  subroutine leftflux(v,w,lw)
    !7th order upwind from Pirozolli / Samtaney
    use constants, only:dz,nlz_par
    implicit none
    real(8)::v
    complex(8),dimension(nlz_par)::w,lw
    complex(8),dimension(-3:nlz_par+3)::wext
    complex(8), dimension(-3:3) :: array_zparbc
    integer::i
    
    call pass_z_array_left_flux(w,array_zparbc)

    !extended index:
    wext(-3)=array_zparbc(-3)
    wext(-2)=array_zparbc(-2)
    wext(-1)=array_zparbc(-1)
    wext(0)=array_zparbc(0)
    wext(nlz_par+1)=array_zparbc(1)
    wext(nlz_par+2)=array_zparbc(2)
    wext(nlz_par+3)=array_zparbc(3)
    
    do i=1,nlz_par
       wext(i)=w(i)
    end do
    
    do i=1,nlz_par
       lw(i)=v/dz*((1./105.*wext(i+3)-19./210.*wext(i+2)+107./210.*wext(i+1)+&
            319./420.*wext(i)-101./420.*wext(i-1)+5./84.*wext(i-2)-1./140.*wext(i-3))-&
            (1./105.*wext(i+2)-19./210.*wext(i+1)+107./210.*wext(i)+&
            319./420.*wext(i-1)-101./420.*wext(i-2)+5./84.*wext(i-3)-1./140.*wext(i-4)))
    end do
    
  end subroutine leftflux


!*********************************************************  
  subroutine rightflux(v,w,rw)
    !7th order upwind from Pirozolli / Samtaney
    use constants, only:dz,nlz_par
    implicit none
    real(8)::v
    complex(8),dimension(nlz_par)::w,rw
    complex(8),dimension(-2:nlz_par+4)::wext
    complex(8), dimension(-2:4) :: array_zparbc
    integer::i
    
    call pass_z_array_right_flux(w,array_zparbc)

    wext(-2)=array_zparbc(-2)
    wext(-1)=array_zparbc(-1)
    wext(0)=array_zparbc(0)
    wext(nlz_par+1)=array_zparbc(1)
    wext(nlz_par+2)=array_zparbc(2)
    wext(nlz_par+3)=array_zparbc(3)
    wext(nlz_par+4)=array_zparbc(4)

    do i=1,nlz_par
       wext(i)=w(i)
    end do
    
    do i=1,nlz_par
       rw(i)=v/dz*((1./105.*wext(i-2)-19./210.*wext(i-1)+107./210.*wext(i)+&
            319./420.*wext(i+1)-101./420.*wext(i+2)+5./84.*wext(i+3)-1./140.*wext(i+4))-&
            (1./105.*wext(i-3)-19./210.*wext(i-2)+107./210.*wext(i-1)+&
            319./420.*wext(i)-101./420.*wext(i+1)+5./84.*wext(i+2)-1./140.*wext(i+3)))
    end do
    
  end subroutine rightflux


!*********************************************************  
 subroutine pass_z_array_left_flux(field,array_zparbc)
    use constants
    use mp!!!, ONLY : send, receive, ipro
    implicit none
    
    integer:: ip!, iz
    complex(8), dimension(nlz_par)::field
    complex(8), dimension(-3:3) :: array_zparbc

!    do iz=1,3
       do ip=npe, npe*npez-1
          if (iproc==ip .and. iproc /= ip-npe) then
             !             call zsend(array_zparbc,ip-npe)
!             call send(field(iz),ip-npe)
             call send(field(1:3),ip-npe)
          end if
       end do
    
       do ip=0,npe*npez-NPE-1
          if (iproc==ip .and. iproc /= ip+npe) then
!             call receive(array_zparbc(iz),ip+NPE)
             call receive(array_zparbc(1:3),ip+NPE)
          end if
       end do

       do ip=0, npe-1
          if (iproc==ip .and. iproc /= npe*npez-npe+ip ) then
             !          call zsend(array_zparbc,npe*npez-npe+ip)
!             call send(field(iz),npe*npez-npe+ip)
             call send(field(1:3),npe*npez-npe+ip)
          end if
       end do

       do ip=npe*npez-NPE, npe*npez-1
          if (iproc==ip .and. iproc /= mod(ip,npe)) then
!             call receive(array_zparbc(iz),mod(ip,NPE))
             call receive(array_zparbc(1:3),mod(ip,NPE))
          end if
       end do
!       call barrier
!    end do

!    do iz=0,3
       do ip=0, npe*npez-1-npe
          if (iproc==ip .and. iproc /= ip+npe) then
!             call send(field(nlz_par-iz),ip+npe)
             call send(field(nlz_par-3:nlz_par),ip+npe)
          end if
       end do

       do ip=npe,npe*npez-1
          if (iproc==ip .and. iproc /= ip-npe) then
!             call receive(array_zparbc(-iz),ip-npe)
             call receive(array_zparbc(-3:0),ip-npe)
          end if
       end do

       do ip=npe*npez-npe, npe*npez-1
          if (iproc==ip .and. iproc /= mod(ip,npe)) then
!             call send(field(nlz_par-iz),mod(ip,npe))
             call send(field(nlz_par-3:nlz_par),mod(ip,npe))
          end if
       end do

       do ip=0, npe-1
          if (iproc==ip .and. iproc /= npe*npez-npe+ip) then
!             call receive(array_zparbc(-iz),npe*npez-npe+ip)
             call receive(array_zparbc(-3:0),npe*npez-npe+ip)
          end if
       end do
!       call barrier

    
  end subroutine pass_z_array_left_flux


!*********************************************************  
  subroutine pass_z_array_right_flux(field,array_zparbc)
    use constants
    use mp
    implicit none
    
    integer :: ip!, iz
    complex(8), dimension(nlz_par)::field
    complex(8), dimension(-2:4) :: array_zparbc

!    do iz=1,4
       do ip=npe, npe*npez-1
          if (iproc==ip .and. iproc /= ip-npe) then
             !             call zsend(array_zparbc,ip-npe)
!             call send(field(iz),ip-npe)
             call send(field(1:4),ip-npe)
          end if
       end do
       do ip=0,npe*npez-NPE-1
          if (iproc==ip .and. iproc /= ip+npe) then
!             call receive(array_zparbc(iz),ip+NPE)
             call receive(array_zparbc(1:4),ip+NPE)
          end if
       end do

       do ip=0, npe-1
          if (iproc==ip .and. iproc /= npe*npez-npe+ip ) then
             !          call zsend(array_zparbc,npe*npez-npe+ip)
!             call send(field(iz),npe*npez-npe+ip)
             call send(field(1:4),npe*npez-npe+ip)
          end if
       end do
    
       do ip=npe*npez-NPE, npe*npez-1
          if (iproc==ip .and. iproc /= mod(ip,npe)) then
!             call receive(array_zparbc(iz),mod(ip,NPE))
             call receive(array_zparbc(1:4),mod(ip,NPE))
          end if
       end do
!       call barrier
!    end do

!    do iz=0,2
       do ip=0, npe*npez-1-npe
          if (iproc==ip .and. iproc /= ip+npe) then
!             call send(field(nlz_par-iz),ip+npe)
             call send(field(nlz_par-2:nlz_par),ip+npe)
          end if
       end do

       do ip=npe,npe*npez-1
          if (iproc==ip .and. iproc /= ip-npe) then
!             call receive(array_zparbc(-iz),ip-npe)
             call receive(array_zparbc(-2:0),ip-npe)
          end if
       end do

       do ip=npe*npez-npe, npe*npez-1
          if (iproc==ip .and. iproc /= mod(ip,npe)) then
!             call send(field(nlz_par-iz),mod(ip,npe))
             call send(field(nlz_par-2:nlz_par),mod(ip,npe))
          end if
       end do

       do ip=0, npe-1
          if (iproc==ip .and. iproc /= npe*npez-npe+ip) then
!             call receive(array_zparbc(-iz),npe*npez-npe+ip)
             call receive(array_zparbc(-2:0),npe*npez-npe+ip)
          end if
       end do
!       call barrier
!    end do
    
  end subroutine pass_z_array_right_flux

end module Fluxes
