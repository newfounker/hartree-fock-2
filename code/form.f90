!  This subroutine is used to calculate e-e Coulomb potential.
!  This routine returns  temp(r) = int dr' fun(r') * f(r<) * g(R>), where
!  the range of integration is zero to infinity. 
! fun(r) : is passed via argument list together with its start and last points, it already contains the Simson's integration weights.
! f(r) = r**l, and g(r) = 1/r**(l+1): are passed via module rpow.

! l : angular momentum
! nr : size of r-grid
! maxfm : index to the maximum r value for which temp(i) need to be calculated.

subroutine form(l,fun,minfun,maxfun,maxfm,temp,i1,i2)
  use grid_radial 

  implicit none
  
  integer, intent(in):: l
  real*8, dimension(grid%nr), intent(in):: fun
  integer, intent(in):: minfun, maxfun, maxfm
  real*8, dimension(grid%nr), intent(out):: temp
  integer, intent(out):: i1, i2
  
  real*8, dimension(0:grid%nr):: t1, t2
  integer:: i
  integer:: if1, if2, istop, it1max, it2min

  if(l .gt. grid%ltmax) then
     i1=2
     i2=1     
     return
  endif

! Set correct limits for integration: \int dr' fun(r') f(r')
  if1 = max(minfun,grid%i1_f(l))
  if2 = min(maxfun,grid%i2_g(l))


  if (if2 .le. if1 + 2) then
     i1=2
     i2=1
     return
  end if

! do not initialise array temp to zero for faster calculation:
!  temp = 0.0

  istop=min(if2,maxfm)

!  Find the integral for fun(r') * f(r') from 0 to istop. The function fun 
!  already contains the Simson's integration weights. 
! limits for which integral is required are:  if1:istop
         t1(if1 - 1) = 0.0
         do i=if1,istop
            t1(i) = t1(i-1) + fun(i)*grid%rpow_f(i,l)
         end do 
! account for the case when:   maxfun < it1max=min(maxfm,i2_g),
! note that in this case istop = maxfun
         it1max = min(maxfm,grid%i2_g(l))
         if(maxfun .lt. it1max) then
            t1(istop+1:it1max) = t1(istop)      
         else
            it1max = istop
         endif
! limits now became:  if1:it1max
!         write(*,'("t1=",3F15.8)') t1(maxfm),t1(maxfm-1), t1(maxfm-10)


!  Find the integral of fun(r') * g(r') from infinity to if1
! limits for which integral should be taken:  if1:if2
         t2(if2) = 0.0
         do i=if2,if1,-1
            t2(i-1) = t2(i) + fun(i)*grid%rpow_g(i,l)
         end do 
! account for the case when:   minfun > i1_f(l)
         if(minfun .gt. grid%i1_f(l)) then
            t2(grid%i1_f(l):if1-1) = t2(if1)
            it2min = grid%i1_f(l)
         else
            it2min = if1
         endif
! limits now became:  it2min:if2, but the upper limit is istop, so the final limits are:
!  it2min:istop

!  Make the form factor by summing two parts

         temp(if1:it1max) = t1(if1:it1max)*grid%rpow_g(if1:it1max,l)

! if array tmp is initialised to zero then:
!         temp(it2min:istop) = temp(it2min:istop) +  t2(it2min:istop)*rpow_f(it2min:istop,l)
! if not:
         if(if1 .le. it2min) then
            if(it1max .ge. istop) then
               temp(it2min:istop) = temp(it2min:istop) +  t2(it2min:istop)*grid%rpow_f(it2min:istop,l)
            else
               temp(it2min:it1max) = temp(it2min:it1max) +  t2(it2min:it1max)*grid%rpow_f(it2min:it1max,l)
               temp(it1max:istop) =  t2(it1max:istop)*grid%rpow_f(it1max:istop,l)
            end if
         else
            if(it1max .ge. istop) then
               temp(if1:istop) = temp(if1:istop) +  t2(if1:istop)*grid%rpow_f(if1:istop,l)
               temp(it2min:if1) =  t2(it2min:if1)*grid%rpow_f(it2min:if1,l)
            else
               temp(if1:it1max) = temp(if1:it1max) +  t2(if1:it1max)*grid%rpow_f(if1:it1max,l)
               temp(it2min:if1) =  t2(it2min:if1)*grid%rpow_f(it2min:if1,l)
               temp(it1max:istop) =  t2(it1max:istop)*grid%rpow_f(it1max:istop,l)
            end if
         endif
! limits: 
         i1 = min(if1,it2min)
         i2 = max(it1max,istop)

end subroutine form

