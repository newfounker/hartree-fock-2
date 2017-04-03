! This is matrix element of V(1,2) electron-electron potential for 4 functions:
! <n1 n2 | V(1,2) | n1p n2p >, where V(1,2) = sum_{lam} V_{lam}(1,2)
! Note: potential V(12) can not change the total ang.mom. projection M of a configuration
subroutine V12me(pn1,pn2,pn1p,pn2p,m1,m2,m1p,m2p,result)
  use input_data
  use grid_radial
  use sturmian_class

  type(sturmian_nr), intent(in):: pn1,pn2,pn1p,pn2p
  integer, intent(in):: m1, m1p, m2, m2p
  real*8, intent(out):: result

  real*8:: Yint

  integer:: l1, l2, l1p, l2p
  integer:: maxfm, minfm, i1, i2, minfun, maxfun, minfun1, maxfun1, lam, lammin, lammax
  real*8:: rlam, rq, reslam
  real*8:: tmp, tmp1, tmp2, sum1, sum2
  integer:: minf1,maxf1, minf1p,maxf1p,minf2,maxf2, minf2p,maxf2p
  real*8, pointer, dimension(:):: f1, f1p, f2, f2p
  real*8, dimension(grid%nr):: temp, fun, fun1, fun11
  real*8:: rl1,rl2,rl1p,rl2p,rm1, rm2,rm1p,rm2p

  result = 0d0

  l1 = get_ang_mom(pn1)
  l2 = get_ang_mom(pn2)
  l1p = get_ang_mom(pn1p)
  l2p = get_ang_mom(pn2p)
  rl1 = l1
  rl2 = l2
  rl1p = l1p
  rl2p = l2p

  rm1 = m1
  rm2 = m2
  rm1p = m1p
  rm2p = m2p

  f1 => fpointer(pn1)
  f1p => fpointer(pn1p)

  f2 => fpointer(pn2)
  f2p => fpointer(pn2p)

  minf1 = get_minf(pn1)
  maxf1 = get_maxf(pn1)
  minf1p = get_minf(pn1p)
  maxf1p = get_maxf(pn1p)

  minf2 = get_minf(pn2)
  maxf2 = get_maxf(pn2)
  minf2p = get_minf(pn2p)
  maxf2p = get_maxf(pn2p)


  maxfm = min(maxf1,maxf1p)
  minfm = max(minf1,minf1p)

  fun1(minfm:maxfm) = f1(minfm:maxfm)*f1p(minfm:maxfm)*grid%weight(minfm:maxfm)

  minfun = max(minf2,minf2p)
  maxfun = min(maxf2,maxf2p)

  fun(minfun:maxfun) = f2(minfun:maxfun)*f2p(minfun:maxfun) * grid%weight(minfun:maxfun)

  lammin=max(abs(l1-l1p),abs(l2-l2p))
  lammax=min(l1+l1p,l2+l2p)

  do lam=lammin,lammax

    rlam = lam
    rq = m1p - m1
    tmp1 = (-1)**(nint(rq))*Yint(rl1,rm1,rlam,-rq,rl1p,rm1p)
    tmp2 = Yint(rl2,rm2,rlam,rq,rl2p,rm2p)
    tmp = tmp1 * tmp2

    !     print*, 'lam, q, tmp', lam, rq,tmp1,tmp2

    if( tmp .eq. 0d0) then
      cycle
    endif

    call form(lam,fun,minfun,maxfun,maxfm,temp,i1,i2)

    minfun1 = max(i1,minfm)
    maxfun1 = min(i2,maxfm)

    fun11(minfun1:maxfun1) = fun1(minfun1:maxfun1) * temp(minfun1:maxfun1)

    reslam = SUM(fun11(minfun1:maxfun1))

    result = result + reslam * tmp

  enddo

end subroutine V12me
