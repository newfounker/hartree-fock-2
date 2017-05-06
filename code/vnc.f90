!!$ vnc(r) = VLambda(r)
!!$  r -> 0,    vnc(r) -> -2*Za/r  (for R=0)
!!$  r -> inf,  vcentre(r) ->  -Zas/r
!!$  in this case Zas = 2*Za = 2


subroutine construct_vnc(nr, gridr)
  use input_data
  use vnc_module
  implicit none

  integer, intent(in):: nr
  real*8, dimension(nr), intent(in):: gridr
  !
  integer:: lamtop_vc_set

  lamtop_vc = data_in%ltmax  ! change in future if required...
  allocate(vnc(nr,0:lamtop_vc),minvnc(0:lamtop_vc),maxvnc(0:lamtop_vc))
  vnc(:,:) = 0d0

  call  VLambdaR(nr,gridr, lamtop_vc, data_in%Z1, data_in%Z2, data_in%Rd, data_in%origin, vnc,minvnc,maxvnc,lamtop_vc_set)

  if( lamtop_vc_set .ne. lamtop_vc ) then

    lamtop_vc = lamtop_vc_set   ! interesting :)

  endif

  !    print*, 'module vnc: test vc: ', vnc(1,0),vnc(1,0)*gridr(1), Za


end subroutine construct_vnc
!---------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine VLambdaR(maxr,gridr, lambda_max, Z1, Z2, Rd, origin, vlr,minvnc,maxvnc, lamtop_vc_set)
  implicit none

  integer, intent(in):: maxr
  real*8, dimension(maxr),intent(in):: gridr
  integer, intent(in):: lambda_max
  real*8, intent(in)  :: Z1, Z2
  real*8,  intent(in) :: Rd, origin
  real*8, dimension(maxr,0:lambda_max), intent(out):: vlr
  integer, dimension(0:lambda_max), intent(out):: minvnc,maxvnc
  integer, intent(out):: lamtop_vc_set

  integer::  i, lambda, iRd, i1, i2, lambda_max_local
  real*8:: const1, r, Rdh, Za
  real*8 :: R1, R2
  real*8, dimension(maxr):: temp

  lambda_max_local = lambda_max

  vlr(:,:) = 0d0
  maxvnc(:) = 1

  write (*, "(a)") "> nuclear potential"

  if ((abs(Z1 - Z2) < 1.0e-5) .and. (2* origin - 1.0 < 1.0e-5)) then

    write (*, *) 'symmetric diatomic molecule'
    write (*, *)

    Rdh = Rd / 2d0  !  distance from COM to a nuclear

    if(Rd .eq. 0d0) then ! joint nuclear, spherically symmetric case, only lam=0 term

      iRd = 1
      lambda_max_local = 0

    else
      do i = 1, maxr
        r = gridr(i)
        if(r .gt. Rdh) exit
      enddo

      iRd = i

    endif

    Za = Z1 + Z2

    !??         print*, 'lambda_max=',lambda_max, lambda_max_local
    do lambda = 0, lambda_max_local, 2
      temp(:) = 0d0

      do i = 1, iRd
        r = gridr(i)
        temp(i) = -Za * (r**lambda/Rdh**(lambda + 1))
      end do
      do i = iRd, maxr
        r = gridr(i)
        temp(i) = -Za * (Rdh**lambda/r**(lambda + 1))
      end do

      call minmaxi(temp,maxr,i1,i2)
      vlr(i1:i2,lambda) = temp(i1:i2)
      minvnc(lambda) = i1
      maxvnc(lambda) = i2
      !??            print*, '>> VLambdaR:  lam,i1,i2=', lambda,i1,i2
    end do

    lamtop_vc_set = lambda_max_local

  else

    if ((abs(Z1) < 1.0e-5) .or. (abs(Z2) < 1.0e-5)) then

      write (*, *) 'atom'
      write (*, *)

    else

      write (*, *) 'asymmetric diatomic molecule'
      write (*, *)

    end if

    R1 = Rd * origin
    R2 = Rd * (1 - origin)

    do lambda = 0, lambda_max_local
      temp(:) = 0d0

      do i = 1, maxr

        r = gridr(i)

        temp(i) = - ( (Z1 * (min(r, R1) ** lambda) / (max(r, R1) ** (lambda + 1))) + &
            (((-1) ** lambda) * Z2 * (min(r, R2) ** lambda) / (max(r, R2) ** (lambda + 1))))


      end do

      call minmaxi(temp,maxr,i1,i2)
      vlr(i1:i2,lambda) = temp(i1:i2)
      minvnc(lambda) = i1
      maxvnc(lambda) = i2

    end do

    lamtop_vc_set = lambda_max_local

  end if


end subroutine VLambdaR

!-----------------------------------------------------------------------------------------
!!$  the potential has rank lam (up to max value lamtop_vc) and its projection=0
!!$  so that the magnetic sublevel projections of pi and pj must be the same
function VLambdaR_ME(pi,pj,ma)
  use input_data
  use grid_radial
  use sturmian_class
  use vnc_module

  implicit none

  real*8:: VLambdaR_ME
  type(sturmian_nr), intent(in):: pi, pj
  integer, intent(in):: ma

  real*8, pointer, dimension(:):: fi, fj, weight, gridr

  integer:: minfi, maxfi, minfj, maxfj, li,lj, minf, maxf, lam_min, lam_max
  real*8:: tmpres, tmp, tmp1,res
  real*8:: Yint
  integer:: lam, lamtmp, i1, i2

  VLambdaR_ME = 0d0

  weight => grid%weight


  tmpres = 0d0

  fi => fpointer(pi)
  minfi = get_minf(pi)
  maxfi = get_maxf(pi)
  li = get_ang_mom(pi)

  fj => fpointer(pj)
  minfj = get_minf(pj)
  maxfj = get_maxf(pj)
  lj = get_ang_mom(pj)

  minf = max(minfi,minfj)
  maxf = min(maxfi,maxfj)


  !      if(data_in%iSlater .eq. 0) then
  !         call fcexch_nr(pi,pj,res)
  !         tmpres = tmpres + res
  !      endif

  lam_min = abs(li-lj)
  lam_max = min(abs(li+lj),lamtop_vc)

  if((-1)**lam_min .eq. -1) lam_min = lam_min + 1

  tmp = 0d0
  do lam=lam_min,lam_max,2
    tmp1 = Yint(dble(li),dble(ma),dble(lam),dble(0),dble(lj),dble(ma))
    !                   print*, i,j, lam, tmp1
    if(tmp1 .eq. 0) cycle
    i1 = max(minvnc(lam),minf)
    i2 = min(maxvnc(lam),maxf)
    tmp = tmp + tmp1 * SUM(weight(i1:i2)*fi(i1:i2)*fj(i1:i2)*vnc(i1:i2,lam))
  enddo
  tmpres = tmpres + tmp

  VLambdaR_ME = tmpres

end function VLambdaR_ME
