module sturmian_class

  public::  new_basis, copy_basis, init_function, copy, multiply_byscalar, construct, destruct, basis_size, get_k, get_ang_mom, get_ang_mom_proj, &
      value, get_minf, get_maxf, fpointer, me_1el, ovlp3_nr, get_energy, print_energy, sort_by_energy, print_wf, ovlp, set_k, set_ang_mom_proj, get_alpha

  private:: construct_nr, destruct_nr, value_nr
!
  type, public:: sturmian_nr
     private
     real*8:: en  ! enrgy of one-electron state, for Sturmians en is set to zero: en=0.0
     integer:: l  ! angular momentum of sturmian function
     integer:: m  ! projection m
     integer:: k  ! order of sturmian function, or principal q. number, or redundant for one-e.states
     integer:: minf, maxf
     real*8, pointer, dimension(:) :: f => NULL()
     real*8:: alpha

  end type sturmian_nr
!
!
  type, public:: basis_sturmian_nr
!     private
     type(sturmian_nr), pointer, dimension(:) :: b  => NULL()
     integer:: n, nExt   ! Number of functions in regular & extended basis.
     integer, dimension(:), allocatable :: iArray, lArray
     real*8, pointer, dimension(:,:):: ortint  => NULL()
  end type basis_sturmian_nr
!
!
  interface init_function
     module procedure init_function_nr,  init_function_nr_m, init_function_nr_st, init_function_nr_st_m
  end interface
!
  interface copy
     module procedure copy_nr
  end interface
!
  interface multiply_byscalar
     module procedure multiply_byscalar_nr
  end interface
!
interface new_basis
     module procedure new_basis_nr
  end interface
!
interface copy_basis
     module procedure copy_basis_nr
  end interface
!
  interface construct
     module procedure construct_nr,  construct_wflocalpot_nr, construct_all_nr
  end interface
!
  interface destruct
     module procedure destruct_nr,  destruct_nr_bf
  end interface
!
  interface basis_size
     module procedure basis_size_nr
  end interface
!
  interface ext_basis_size
     module procedure ext_basis_size_nr
  end interface
!
  interface value
     module procedure  value_nr
  end interface
!
  interface get_energy
     module procedure  get_energy_nr
  end interface
!
  interface get_ang_mom
     module procedure get_ang_mom_nr
  end interface
!
  interface set_ang_mom
     module procedure set_ang_mom_nr
  end interface
!
  interface get_ang_mom_proj
     module procedure get_ang_mom_proj_nr
  end interface
!
  interface set_ang_mom_proj
     module procedure set_ang_mom_proj_nr
  end interface
!
  interface get_k
     module procedure get_k_nr
  end interface
!
  interface set_k
     module procedure set_k_nr
  end interface
!
  interface fpointer
     module procedure  fpointer_nr
  end interface
!
  interface get_minf
     module procedure  get_minf_nr
  end interface
!
  interface get_maxf
     module procedure  get_maxf_nr
  end interface
!
  interface print_energy
     module procedure  print_energy_nr
  end interface
!
  interface print_wf
     module procedure  print_wf_nr
  end interface
!
  interface sort_by_energy
     module procedure  sort_by_energy_nr
  end interface
!
!
  interface gsort
     module procedure gsort_nr
  end interface
!
  interface me_1el
     module procedure  me_1el_nr
  end interface
!
!
  interface ovlp
     module procedure  ovlp_nr
  end interface
!
!
contains

  function get_calc_type()
    get_calc_type = calc_type
  end function get_calc_type
!
!
  subroutine init_function_nr(self,l,k,i1,i2,temp,nr,alpha)
    type(sturmian_nr), intent(inout):: self
    integer, intent(in):: l, k, i1, i2, nr
    real*8, dimension(nr), intent(in):: temp
    real*8:: alpha

    self%en = 0d0
    self%l = l
    self%m = 0
    self%k = k
    self%minf = i1
    self%maxf = i2
    ! create array for one-electron function
    allocate(self%f(i2))
    self%f(1:i1) = 0.0
    self%f(i1:i2) = temp(i1:i2)
    self%alpha = 0d0
    self%alpha = alpha
  end subroutine init_function_nr
!
  subroutine init_function_nr_m(self,l,m,k,i1,i2,temp,nr)
    type(sturmian_nr), intent(inout):: self
    integer, intent(in):: l, k, i1, i2, nr
    real*8, dimension(nr), intent(in):: temp

    self%en = 0d0
    self%l = l
    self%m = m
    self%k = k
    self%minf = i1
    self%maxf = i2
    ! create array for one-electron function
    allocate(self%f(i2))
    self%f(1:i1) = 0d0
    self%f(i1:i2) = temp(i1:i2)
    self%alpha = 0d0

  end subroutine init_function_nr_m
!
  subroutine init_function_nr_st(self,en,l,k,i1,i2,temp,nr)
    type(sturmian_nr), intent(inout):: self
    real*8, intent(in):: en
    integer, intent(in):: l, k, i1, i2, nr
    real*8, dimension(nr), intent(in):: temp

    self%en = en
    self%l = l
    self%m = 0
    self%k = k
    self%minf = i1
    self%maxf = i2
    ! create array for one-electron function
    allocate(self%f(i2))
    self%f(1:i1) = 0.0
    self%f(i1:i2) = temp(i1:i2)
    self%alpha = 0d0
  end subroutine init_function_nr_st
!
  subroutine init_function_nr_st_m(self,en,l,m,k,i1,i2,temp,nr)
    type(sturmian_nr), intent(inout):: self
    real*8, intent(in):: en
    integer, intent(in):: l, m, k, i1, i2, nr
    real*8, dimension(nr), intent(in):: temp

    self%en = en
    self%l = l
    self%m = m
    self%k = k
    self%minf = i1
    self%maxf = i2
    ! create array for one-electron function
    allocate(self%f(i2))
    self%f(1:i1) = 0.0
    self%f(i1:i2) = temp(i1:i2)
    self%alpha = 0d0
  end subroutine init_function_nr_st_m
!
!
  subroutine copy_nr(state_l,state_r)
    type(sturmian_nr), intent(out):: state_l
    type(sturmian_nr), intent(in):: state_r
    integer:: i1,i2

    state_l%en = state_r%en
    state_l%l = state_r%l
    state_l%m = state_r%m
    state_l%k = state_r%k
    i1 = state_r%minf
    i2 = state_r%maxf
    state_l%minf = i1
    state_l%maxf = i2
    if(associated(state_l%f)) then
       deallocate(state_l%f)
    endif
    allocate(state_l%f(1:i2))
    state_l%f(i1:i2) = state_r%f(i1:i2)
    state_l%alpha = state_r%alpha
  end subroutine copy_nr
!
!
  subroutine multiply_byscalar_nr(scalar,state)
    real*8, intent(in):: scalar
    type(sturmian_nr), intent(inout):: state
    integer:: i1,i2

    i1 = state%minf
    i2 = state%maxf
    if(associated(state%f)) then
       state%f(i1:i2) = scalar*state%f(i1:i2)
    endif
  end subroutine multiply_byscalar_nr
!
!
  subroutine new_basis_nr(self,n)
    implicit none

    type(basis_sturmian_nr), intent(inout):: self
    integer, intent(in):: n
    integer :: m, noid

    self%n = n

    ! create array of n one-electron functions
    if(n .ne. 0)  then
       allocate( self%b(n), self%ortint(self%n,self%n) )
       self%ortint(:,:) = 0.0d0
    endif

  end subroutine new_basis_nr
!
!
  subroutine copy_basis_nr(basis_l,basis_r)
    type(basis_sturmian_nr), intent(inout):: basis_l,basis_r
!!$    integer, intent(in):: n

    if(basis_l%n .ne. basis_r%n) then
       print*, 'sturmian.f90: Can not copy basis of unequal size'
       stop
    endif

    if(basis_l%n .eq. 0) then
       print*, 'sturmian.f90: basis size is zero'
       stop
    endif

    n = basis_l%n
    do i = 1, n
       call copy(basis_l%b(i), basis_r%b(i))
    end do
    basis_l%ortint(:,:) =  basis_r%ortint(:,:)


  end subroutine copy_basis_nr
!
! create basis of n Sturmian functions and allocate space for each radial part of the function in the basis
  subroutine construct_nr(self,n,l,al)
    use grid_radial
    implicit none
    type(basis_sturmian_nr), intent(inout):: self
    integer, intent(in):: n
    integer, intent(in):: l
    real*8, intent(in):: al

    integer:: i, k, i1, i2, la
    real*8:: f8(grid%nr,n), x, tmp
    real*8, dimension(grid%nr):: temp
    real*8:: lambda
    real*8, pointer, dimension(:):: weight, gridr
    integer:: ortog, m

    m = 0   !   NOTE to be corrected...

!
    weight => grid%weight
!     weight => grid%bool
    gridr => grid%gridr

! create array of n one-electron functions for given l
    call new_basis(self,n)
    lambda = 2d0 * al
    ortog = 0
    if(ortog .eq. 1) then
       call  lagpol8(dble(2*l+2),lambda,f8,grid%nr,n,dble(gridr),grid%nr)
    else
       call  lagpol8(dble(2*l+1),dble(2.0*al),f8,grid%nr,n,dble(gridr),grid%nr)
    endif

    do k=1,n
       if(ortog .eq. 1) then
          tmp = k
          do i=1,2*l+1
             tmp = tmp*(k+i)
          enddo
          tmp = sqrt(lambda/tmp)
       else
          tmp =  2.0*(k+l)
          do i=0,2*l
             tmp = tmp*(k+i)
          enddo
          tmp = sqrt(tmp)
          tmp = sqrt(2d0*al)/tmp
       endif

       do i=1,grid%nr
          x = gridr(i) * 2d0*al
          temp(i) = tmp*(x**(l+1))*exp(-x/2.0D0)*f8(i,k)
       end do


       call minmaxi(temp,grid%nr,i1,i2)
       call init_function(self%b(k),l,k,i1,i2,temp,grid%nr,al)
!       print*, 'k=',k,', size=',SIZE(self%b(k)%f)
    enddo



!!$ make overlap matrix
    la = l
    do k=1,n
       self%ortint(k,k) = 1d0
    enddo
    do k=1,n-1
          self%ortint(k,k+1) =  -0.5*sqrt(1.0-dble(la*(la+1))/dble((la+k)*(la+k+1)))
          self%ortint(k+1,k) = self%ortint(k,k+1)
    enddo


!??    write(*,'("Created nonrel. Sturmian basis with n=",i5,"  functions for angular momentum l=",i5)') n, l
  end subroutine construct_nr
!
!-------------------------------------------------------------------------------------------------!
! Make a radial basis of all Sturmians.
! Functions are in the form N * exp(-x/2) * x^a/2 * L_n^k(x). If a = k the basis is orthogonal.
! Not originally written, but generalised to (non)orthogonal spherical/spheroidal on 28/9/11 by JS.
!
! INPUT: self (type basis_sturmian_nr)
!
! OUTPUT: self
!
!
  subroutine construct_all_nr(self, dataARG)

    use grid_radial
    use input_data

    implicit none

    type(basis_sturmian_nr), intent(inout):: self
    type(input), intent(in):: dataARG

    integer:: i, j, k, i1, i2, kn, knext, nall, next, l, labot, latop, n, la, nsp, ni, nj, li, lj, mi, mj
    real*8:: f8(grid%nr,MAXVAL(dataARG%nps(dataARG%labot:dataARG%latop))), x, tmp
    real*8, dimension(grid%nr):: temp
    real*8:: lambda, lagarg
    real*8, pointer, dimension(:):: weight, gridr
    integer:: ortog, m

    ortog = 0   ! Orthogonality of the basis functions. 0: nonorthogonal (2l+1 type); 1: orthogonal (2l+2 type); 2: spheroidal orthogonal (m type).

    weight => grid%weight
    gridr => grid%gridr

    ! create array of nall one-electron functions for given all l (or m in spheroidal case)
    labot = dataARG%labot
    latop = dataARG%latop

    nall = SUM(dataARG%nps(labot:latop))
    call new_basis(self, nall)

    kn = 0   ! Counter of functions calculated.
    knext = 0   ! Counter of functions in the extended basis.
    do l=labot,latop

       n = dataARG%nps(l)

       lambda = 2d0 * dataARG%alpha(l)

       ! Calculate the Laguerre polynomials up to order n on the given r-grid.
       f8(:,:) = 0d0
       if(ortog .eq. 0) then   ! Nonorthogonal spherical.
          lagarg = dble(2*l+1)
       elseif(ortog .eq. 1) then   ! Orthogonal spherical.
          lagarg = dble(2*l+2)
       endif

       call lagpol8( lagarg, lambda, f8, grid%nr, n, dble(gridr), grid%nr )

       ! Determine the normalisation constant.
!       print*, '*** ortog=', ortog
       do k=1,n
          if(ortog .eq. 0) then   ! Nonorthogonal--extra factor of 1/(2k+2l).
             tmp = 2.0d0*(k+l)
          elseif(ortog .eq. 1) then   ! Orthogonal.
             tmp = 1.0d0
          ! elseif(ortog .eq. 2) then   ! Spheroidal.
          !    tmp = 4.0*acos(0.0)   ! 2pi
          endif
          do i=0, nint(lagarg)-1   ! Factorial quotient (k-1+lagarg)!/(k-1)!
             tmp = tmp*(k+i)
          enddo
          tmp = sqrt( lambda / tmp )

          ! Multiply the Laguerre polynomials by the normalisation, exponential and polynomial factors.
          if(ortog .eq. 0) then   ! Nonorthogonal--correction to the power of x.
             lagarg = lagarg + 1.0d0   ! 2l+1 -> 2(l+1)
          endif
          do i=1,grid%nr
             x = gridr(i) * lambda
             temp(i) = tmp * x**(lagarg/2.0d0) * exp(-x/2.0d0) * f8(i,k)
          end do
          if( ortog .eq. 0 ) then   ! Uncorrection correction.
             lagarg = lagarg - 1.0d0   ! 2(l+1) -> 2l+1
          endif

          call minmaxi(temp,grid%nr,i1,i2)

          kn = kn + 1

          ! Initialise each Sturmian function (type sturmian_nr) within the basis.
          call init_function(self%b(kn),l,k,i1,i2,temp,grid%nr,lambda/2d0)
!          print*, 'k=',k,', size=',SIZE(self%b(k)%f)

!!$open(unit=1336+kn)
!!$write(1336+kn, '(5(A,I4))') '# k=', k, '   l=', l, '   m=', m, '   i1=', i1, '   i2=', i2
!!$do i=i1,i2
!!$   write(1336+kn, '(100F20.10)') grid%gridr(i), temp(i)
!!$enddo
!!$close(unit=1336+kn)

       enddo ! k -- each loop creates a single basis function.
       write(*,'(2(A,I2))'), 'Created radial basis with n=',n, ' Sturmian functions for angular momentum=',l
    enddo ! l -- one loop for each labot<=l<=latop.
!print*, 'iArray:', self%iArray
!print*, 'lArray:', self%lArray
    ! Check if all requested functions have been created.
    if(nall .ne. kn) then
       print*,'sturmian.f90: nall != kn:', nall,  kn
       stop
    endif

    self%n = nall

    ! Overlap matrix.
    if(ortog.eq.0 .or. ortog.eq.1) then
       do i=1,nall   ! nallxnall identity matrix.
          self%ortint(i,i) = 1d0
       enddo
       if(ortog .eq. 0) then   ! Nonorthogonal--off-diagonal overlap matrix elements.
          do nsp=1,nall-1
             la = get_ang_mom(self%b(nsp))
             i = get_k(self%b(nsp))
             if( la .eq.  get_ang_mom(self%b(nsp+1))) then
                self%ortint(nsp,nsp+1) =  -0.5*sqrt(1.0-dble(la*(la+1))/dble((la+i)*(la+i+1)))
                self%ortint(nsp+1,nsp) = self%ortint(nsp,nsp+1)
             endif
          enddo
       endif
    endif

    ! Code for verifying the overlap matrix if one wishes to do so.
    if( 2 + 2 == 5 ) then
       do n=1,nall
          do i=1,n
             call ovlp_nr( self%b(n), self%b(i), tmp )
             write(*,'(2I5, 2F20.10)') n, i, tmp, self%ortint(n,i)
          enddo
          print*
       enddo
    endif


    write(*,'(3(A,I3))'), 'Created radial basis with n=',nall, ' Sturmian functions for angular momentum from ',labot, ' to ',latop


  end subroutine construct_all_nr
!
!-------------------------------------------------------------------------------------------------!
!
  subroutine destruct_nr_bf(self)
    implicit none
    type(sturmian_nr), intent(inout):: self
    integer:: stat

    deallocate(self%f)
  end subroutine destruct_nr_bf
!
!
  subroutine destruct_nr(self)
    implicit none
    type(basis_sturmian_nr), intent(inout):: self
    integer:: n
    integer:: i
    integer:: stat

    n= self%n
    do i=1,n
!       print*,'dealocate: i=', i, ', size=', SIZE(self%b(i)%f)
       deallocate(self%b(i)%f)
    enddo
    deallocate(self%b, STAT=stat)
  end subroutine destruct_nr
!
!
  function basis_size_nr(self)
    implicit none
    integer:: basis_size_nr
    type(basis_sturmian_nr), intent(in):: self

    basis_size_nr = self%n
  end function basis_size_nr
!
!
  function ext_basis_size_nr(self)
    implicit none
    integer :: ext_basis_size_nr
    type(basis_sturmian_nr), intent(in) :: self

    ext_basis_size_nr = self%nExt
  end function ext_basis_size_nr
!
!
  function value_nr(self,i)
    implicit none
    real*8:: value_nr
    type(sturmian_nr), intent(in):: self
    integer, intent(in):: i  !

    if( i .gt. self%maxf .or.  i .lt. self%minf) then
       value_nr = 0.0
    else
       value_nr = self%f(i)
    endif

  end function value_nr
!
!
  function get_ang_mom_nr(self)
    implicit none
    integer:: get_ang_mom_nr
    type(sturmian_nr), intent(in):: self

    get_ang_mom_nr = self%l
  end function get_ang_mom_nr
!
  subroutine set_ang_mom_nr(self,L)
    implicit none
    type(sturmian_nr):: self
    integer, intent(in):: L

    self%l = L
  end subroutine set_ang_mom_nr
!
!
  function get_ang_mom_proj_nr(self)
    implicit none
    integer:: get_ang_mom_proj_nr
    real*8:: get_ang_mom_nr
    type(sturmian_nr), intent(in):: self

    get_ang_mom_proj_nr = self%m
  end function get_ang_mom_proj_nr
!
  subroutine set_ang_mom_proj_nr(self,M)
    implicit none
    type(sturmian_nr):: self
    integer, intent(in):: M

    self%m = M
  end subroutine set_ang_mom_proj_nr
!
!
function   get_min_L(self)
  implicit none
  integer:: get_min_L
  type(basis_sturmian_nr), intent(in):: self
  integer:: n, l, ltmp, Nmax

  Nmax = basis_size(self)
  if(Nmax .le. 0) then
     print*,"sturmian.f90:  Nmax <=0 in get_min_L(self)"
     stop
  endif

  l = get_ang_mom(self%b(1))
  do n=2,Nmax
     ltmp = get_ang_mom(self%b(n))
     l = min(l,ltmp)
  enddo

  get_min_L = l

  return
end function get_min_L
!
function   get_max_L(self)
  implicit none
  integer:: get_max_L
  type(basis_sturmian_nr), intent(in):: self
  integer:: n, l, ltmp, Nmax

  Nmax = basis_size(self)
  if(Nmax .le. 0) then
     print*,"sturmian.f90:  Nmax <=0 in get_max_L(self)"
     stop
  endif

  l = get_ang_mom(self%b(1))
  do n=2,Nmax
     ltmp = get_ang_mom(self%b(n))
     l = max(l,ltmp)
  enddo

  get_max_l = l

  return
end function get_max_L
!
!  Finds how many orbitals with given L are in the basis
function   get_max_nL(self,l)
  implicit none
  integer:: get_max_nL
  type(basis_sturmian_nr), intent(in):: self
  integer, intent(in):: l
  integer:: n, ltmp, Nmax, numk

  Nmax = basis_size(self)
  if(Nmax .le. 0) then
     print*,"sturmian.f90:  Nmax <=0 in get_maxnL(self,l)"
     stop
  endif

  numk = 0
  do n=1,Nmax
     ltmp = get_ang_mom(self%b(n))
     if(l .eq. ltmp) then
        numk = numk + 1
     endif
  enddo

  get_max_nL = numk

  return
end function get_max_nL
!
!  Finds how many orbitals with given L are in the basis up to a basis function M
function   get_nL_uptoM(self,l,M)
  implicit none
  integer:: get_nL_uptoM
  type(basis_sturmian_nr), intent(in):: self
  integer, intent(in):: l, M
  integer:: n, ltmp, Nmax, numk

  Nmax = basis_size(self)
  if(Nmax .le. 0) then
     print*,"sturmian.f90:  Nmax <=0 in get_nL_uptoM(self)"
     stop
  endif

  if(M .le. 0) then
     print*,"sturmian.f90:  M <=0 in get_min_kappa_uptoM(self,kappa,M)"
     stop
  endif

  if(M .gt. Nmax) then
     print*,"sturmian.f90:  M <= Nmax in get_nL_upto(self,kappa,M)"
     stop
  endif


  numk = 0
  do n=1,M
     ltmp = get_ang_mom(self%b(n))
     if(l .eq. ltmp) then
        numk = numk + 1
     endif
  enddo

  get_nL_uptoM = numk

  return
end function get_nL_uptoM
!
!  Finds the maximum value of orbitals for given L over all possible L are in the basis
function   get_maxall_nL(self)
  implicit none
  integer:: get_maxall_nL
  type(basis_sturmian_nr), intent(in):: self
  integer:: l, l_min, l_max, numk, numktmp

  l_max = get_max_L(self)
  l_min = get_min_L(self)

  numk = 0
  do l=l_min,l_max
     numktmp =  get_max_nl(self,l)
     numk = max(numk,numktmp)
  enddo

  get_maxall_nL = numk

  return
end function get_maxall_nL
!
!
  function get_k_nr(self)
    implicit none
    integer:: get_k_nr
    type(sturmian_nr), intent(in):: self

    get_k_nr = self%k
  end function get_k_nr
!
  subroutine set_k_nr(self,k)
    implicit none
    type(sturmian_nr):: self
    integer, intent(in):: k

    self%k = k
  end subroutine set_k_nr
!
!
  function fpointer_nr(self)
    implicit none
    real*8, pointer, dimension(:):: fpointer_nr
    type(sturmian_nr), intent(in):: self

    fpointer_nr => self%f
  end function fpointer_nr
!
!
  function get_minf_nr(self)
    implicit none
    integer:: get_minf_nr
    type(sturmian_nr), intent(in):: self

    get_minf_nr = self%minf
  end function get_minf_nr
!
!
  function get_maxf_nr(self)
    implicit none
    integer:: get_maxf_nr
    type(sturmian_nr), intent(in):: self

    get_maxf_nr = self%maxf
  end function get_maxf_nr
!
!
  function get_alpha(self)
    implicit none
    real*8:: get_alpha
    type(sturmian_nr), intent(in):: self

    get_alpha = self%alpha
  end function get_alpha
!
!
 function get_energy_nr(self)
    implicit none
    real*8:: get_energy_nr
    type(sturmian_nr), intent(in):: self

    get_energy_nr = self%en
  end function get_energy_nr
!
  subroutine  print_energy_nr(self)
    implicit none
    type(basis_sturmian_nr), intent(in):: self
    integer:: n
    real*8:: ioniz_en, exit_en

     write(*,'("    N     J     exitation en.   ionization en.")')
    do n=1,self%N
       exit_en = (self%b(n)%en-self%b(1)%en)*27.2116
       ioniz_en = (self%b(n)%en)*27.2116  ! in eV
       write(*,'(i5,I6,2F17.5)') n, get_ang_mom(self%b(n)), exit_en, ioniz_en
    enddo
    print*

  end subroutine print_energy_nr
!
!
  subroutine  print_wf_nr(self,N,Rmax,filename)
    use grid_radial

    implicit none
    type(basis_sturmian_nr), intent(in):: self  ! basis
    integer, intent(in):: N                     ! number of functions to be printed
    real*8, intent(in):: Rmax                   ! max value of R
    character(LEN=40), intent(in):: filename    ! name of file where to print functions
    integer:: i,m

    open(150,file=filename,status='REPLACE')

    write(150,'("# n:",12X,100(I10,5X))') (get_k(self%b(m)), m=1,N)
    write(150,'("# kappa:",8X,100(I10,5X))') (get_ang_mom(self%b(m)), m=1,N)
    write(150,'("# energy:",8X,100(10X,ES15.6,5X))') (get_energy(self%b(m)), m=1,N)
    do i=1,grid%nr
       if(grid%gridr(i) .le. grid%rmax .and. grid%gridr(i) .le. Rmax) then
          write(150,'(F15.5,100(1X,E14.5))') grid%gridr(i), (value(self%b(m),i), m=1,N)
       endif
    enddo
    close(150)

  end subroutine print_wf_nr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shellsort algorithm
  subroutine sort_by_energy_nr(self)
    type(basis_sturmian_nr), intent(inout):: self
    integer:: gap, i, j, N
    type(sturmian_nr):: Tmp

    N = self%n
    gap = N/2
    do
       if(gap .le. 0) exit
       do i=gap+1,N
          call copy(Tmp,self%b(i))
          do j=i,gap+1,-gap
             if(Tmp%en .lt. self%b(j-gap)%en) then
                call copy(self%b(j),self%b(j-gap))
             else
                exit
             endif
             call copy(self%b(j-gap),Tmp)
          enddo
       enddo
       if ( gap .eq. 2) then
          gap = 1
       else
          gap = gap/2.2
       endif
    enddo

  end subroutine sort_by_energy_nr
!
!******************************************************************
!
   subroutine ovlp3_nr(pn,pnp,v,i1,i2,result)
      use grid_radial
      implicit none

      type(sturmian_nr), intent(in):: pn,pnp !  one-electron states
      real*8, dimension(grid%nr), intent(in):: v
      integer, intent(in)::  i1, i2
      real*8, intent(out):: result
!
      real*8, dimension(grid%nr):: temp, fun
      real*8, pointer, dimension(:):: f, fp
      integer:: minf,maxf, minfp,maxfp
      integer::  minfun, maxfun
!
      result = 0d0

      if(get_ang_mom(pn) .ne. get_ang_mom(pnp)) return
      if(get_ang_mom_proj(pn) .ne. get_ang_mom_proj(pnp)) return

      f  => fpointer(pn)
      fp => fpointer(pnp)
      minf = get_minf(pn)
      maxf = get_maxf(pn)
      minfp = get_minf(pnp)
      maxfp = get_maxf(pnp)


      result = 0.0
      minfun = max(max(minf,minfp),i1)
      maxfun = min(min(maxf,maxfp),i2)

      fun(minfun:maxfun) = f(minfun:maxfun)*fp(minfun:maxfun)*grid%weight(minfun:maxfun)*v(minfun:maxfun)
      result = SUM(fun(minfun:maxfun))

    end subroutine ovlp3_nr
!
   subroutine ovlp_nr(pn,pnp,result)
      use grid_radial
      implicit none

      type(sturmian_nr), intent(in):: pn,pnp !  one-electron states
      real*8, intent(out):: result
!
      real*8, dimension(grid%nr):: temp, fun
      real*8, pointer, dimension(:):: f, fp
      integer:: minf,maxf, minfp,maxfp
      integer::  minfun, maxfun
!

      result = 0.0
      if(get_ang_mom(pn) .ne. get_ang_mom(pnp)) return

      f  => fpointer(pn)
      fp => fpointer(pnp)
      minf = get_minf(pn)
      maxf = get_maxf(pn)
      minfp = get_minf(pnp)
      maxfp = get_maxf(pnp)


      minfun = max(minf,minfp)
      maxfun = min(maxf,maxfp)

      fun(minfun:maxfun) = f(minfun:maxfun)*fp(minfun:maxfun)*grid%weight(minfun:maxfun)
      result = SUM(fun(minfun:maxfun))

    end subroutine ovlp_nr
!
!
   subroutine me_1el_nr(lam,pn,pnp,v,i1,i2,result)
      use grid_radial
      implicit none

      integer, intent(in):: lam
      type(sturmian_nr), intent(in):: pn,pnp !  one-electron states
      real*8, dimension(grid%nr), intent(in):: v   ! v is a scalar
      integer, intent(in)::  i1, i2
      real*8, intent(out):: result
!
      real*8, dimension(grid%nr):: temp, fun
      real*8, pointer, dimension(:):: f, fp
      integer:: minf,maxf, minfp,maxfp
      integer:: l, lp
      real*8:: j, jp, rlam
      integer:: k
      real*8:: CLAM
      real*8:: tmp
      integer::  minfun, maxfun, i1old, i2old
!


      if(get_ang_mom(pn) .ne. get_ang_mom(pnp)) return
      if(get_ang_mom_proj(pn) .ne. get_ang_mom_proj(pnp)) return

      rlam = lam

      l  =  get_ang_mom(pn)
      lp =  get_ang_mom(pnp)
      f  => fpointer(pn)
      fp => fpointer(pnp)
      minf = get_minf(pn)
      maxf = get_maxf(pn)
      minfp = get_minf(pnp)
      maxfp = get_maxf(pnp)


      result = 0.0
      minfun = max(minf,minfp)
      maxfun = min(maxf,maxfp)
      tmp = sqrt(2*lp +1d0)*CLAM(dble(l),rlam,dble(lp))

      if(tmp .ne. 0d0) then

         fun(minfun:maxfun) = f(minfun:maxfun)*fp(minfun:maxfun)*grid%weight(minfun:maxfun)*v(minfun:maxfun)
         result = tmp * SUM(fun(minfun:maxfun))
      endif

    end subroutine me_1el_nr

!
!!$----------------------------------------------------------------------------------------

subroutine gsort_nr(self)
  use  grid_radial

  implicit none

  type(basis_sturmian_nr), intent(inout):: self

  real*8, dimension(grid%nr):: v
  real*8, dimension(self%n,self%n):: ortint
  real*8, pointer, dimension(:):: f1, f2
  integer:: nspm, n, m, i1, i2, i, nr
  integer:: l1, l2, m1, m2
  real*8:: tmp, sum1
  real*8::  sum2
  integer:: j,k

  write(*,'(" G - S orthogonalization")')

  nr = grid%nr
  nspm= self%n
  v(:) = 0.0

  do n=1,nspm
     f1 => fpointer(self%b(n))
     l1 = get_ang_mom(self%b(n))
     m1 = get_ang_mom_proj(self%b(n))
     i1 = get_minf(self%b(n))
     i2 = get_maxf(self%b(n))
     do m=1,n
        f2 => fpointer(self%b(m))
        l2 = get_ang_mom(self%b(m))
        m2 = get_ang_mom_proj(self%b(m))
        i1 = max(i1,get_minf(self%b(m)))
        i2 = min(i2,get_maxf(self%b(m)))
        tmp = 0.0
        if(l1 .eq. l2 .and. m1 .eq. m2) then
           tmp = SUM( (f1(i1:i2)* f2(i1:i2) )*grid%weight(i1:i2) )
        endif
        ortint(n,m) = tmp
        ortint(m,n) = tmp
     enddo
  enddo

!!$     form overlap array <u_j|v_k>, u_j - old nonorthogonal set,
!!$     v_k - new set but not yet normalized
!!$     Only elements with  j >= k  required.
  do j=1,nspm
     do k=1,j
        sum1 = 0d0
        do n=1,k-1
           sum1 = sum1 + ortint(k,n)*ortint(j,n)/ortint(n,n)
        end do
        ortint(j,k) = ortint(j,k) - sum1
!!$            write(20,'("j,k =",2I3,", <j|k> =",F10.5)')  j, k, real(ortint(j,k))
     end do
  end do


!!$    form new orthonomal set vb_k, f_k = vb_k
  do k=1,nspm
     f1 => fpointer(self%b(k))
     v(:) = 0.0
     do i=1,nr
        sum2 = 0d0
        do n=1,k-1
           f2 => fpointer(self%b(n))
           i1 = get_minf(self%b(n))
           i2 = get_maxf(self%b(n))
           if(i .ge. i1 .and. i .le. i2) then
              sum2 = sum2 + dble(f2(i))*ortint(k,n)/dsqrt(ortint(n,n))
           endif
        end do
        v(i) = (value(self%b(k),i) - sum2)/dsqrt(ortint(k,k))
     end do
     call minmaxi2(v,nr,i1,i2)
     self%b(k)%maxf = i2
     self%b(k)%minf = i1
     deallocate(self%b(k)%f)
     allocate(self%b(k)%f(1:i2))

     tmp = sqrt(SUM( (v(i1:i2)*v(i1:i2))*grid%weight(i1:i2)))
!     print*,k,', tmp=', tmp, ortint(k,k)
     v(i1:i2) = v(i1:i2)/tmp

     self%b(k)%f(1:i2) = v(1:i2)
  end do

  return
end subroutine gsort_nr
!!$------------------------------------------------------------------------------------
! This subroutine make use nonrelativistic Sturmians
!!   V is a scalar
!!   the basis  self  is built from functions with the same angular momentum la  and its projection
  subroutine construct_wflocalpot_nr(self,Nwf,V,nd,la,al)
    use grid_radial
!    use sturmian_class

    implicit none

    type(basis_sturmian_nr), intent(out):: self
    integer, intent(in):: Nwf
    real*8, dimension(grid%nr), intent(in):: V  !!!  V is a scalar
    integer, intent(in):: nd, la
    real*8, intent(in):: al


    type(basis_sturmian_nr):: bst
    real*8, pointer, dimension(:):: weight, gridr
    real*8, dimension(:,:), allocatable:: H, b, CI
    real*8, dimension(:), allocatable:: w
    integer:: i, j, m, n, nr, nrj,  imax, imin, im1,im2, N_oneel, kk
    integer:: matz, ierr
    real*8, pointer, dimension(:):: fi, fj, fm
    integer:: minfi, maxfi, minfj, maxfj, minf, maxf
    real*8:: tmpC
    real*8, dimension(grid%nr):: temp
    real*8:: energy, overlap, tmp, tmp1

    if(Nwf .gt. nd) then
       print*,'sturmian.f90: construct_wflocalpot_nr(): Nwf>nd:',Nwf , nd
    endif

    weight => grid%weight
    gridr => grid%gridr

    call new_basis(self,Nwf)

!!$ For given (la,nd,al) construct Sturmian basis.
    call construct(bst,nd, la, al)

!!$ Temporary arrays
    allocate(H(nd,nd))
    allocate(b(nd,nd))
    allocate(CI(nd,nd))
    allocate(w(nd))

!!$ Calculate H matrix
    b = 0d0
    H = 0d0
!!$ Overlap matrix
    b(1:nd,1:nd) = bst%ortint(1:nd,1:nd)

!!$ get kinetic energy: (using special properties of Lag func)
    H(1:nd,1:nd) = -al*al*b(1:nd,1:nd)/2d0
    do i = 1,nd
       H(i,i) = H(i,i) + al*al
    end do


!!$ Hamiltonian  matrix
    do i = 1,nd
       fi => fpointer(bst%b(i))
       minfi = get_minf(bst%b(i))
       maxfi = get_maxf(bst%b(i))
       do j = 1, i
          fj => fpointer(bst%b(j))
          minfj = get_minf(bst%b(j))
          maxfj = get_maxf(bst%b(j))
          minf = max(minfi,minfj)
          maxf = min(maxfi,maxfj)
          tmp = SUM(weight(minf:maxf)*fi(minf:maxf)*fj(minf:maxf)*V(minf:maxf))
          H(i,j) = H(i,j)  + tmp
          H(j,i) = H(i,j)
       end do
    end do

!!$          write(*,'(1P,5E15.5)') b
!!$          print*
!!$          print*
!!$          write(*,'(1P,5E15.5)') H


    matz=2
    call  rsg(nd,nd,H,b,w,matz,CI,ierr)
    if(ierr .ne. 0) then
       write(*,'("ierr =",I3)') ierr
       stop 'construct_wflocalpot: ERROR IN DIAGONALIZATION'
    endif

!??    write(*,'("Diagonalizing with L, nd,alpha:",2I5,F15.5)')la,nd,al
!??    do j=1, nd
!??       write(*,'(I5,1P,E15.5)') j,w(j)
!??    enddo
!??    print*

    N_oneel = 0
    do j=1,Nwf
       temp = 0d0
       do m=1,nd
          fm => fpointer(bst%b(m))
          tmpC =  CI(m,j)
          imax = get_maxf(bst%b(m))
          imin = get_minf(bst%b(m))
          do i=imin,imax
             temp(i) = temp(i) + tmpC*fm(i)
          enddo
       enddo
       call minmaxi(temp(1:grid%nr),grid%nr,im1,im2)
!!$ Fix sign of the one-electron functions
       if( sign(1d0,temp(im1+1)) .lt. 0)  then
          temp = -temp
       end if
       im1 = min(im1,imin)
       im2 = max(im2,imax)
       N_oneel = N_oneel + 1
       if(N_oneel .gt. Nwf) then
          print*,'N_oneel=',N_oneel, ', Nwf=', Nwf, ' nd=', nd
          stop 'N_oneel .gt. Nwf'
       endif
       ! record one-electron function
       energy = w(j)
       kk = j
       !              print*, '*** n,kk:', N_oneel, kk
       call init_function(self%b(N_oneel),energy,la,kk,im1,im2,temp,grid%nr)

    end do

    deallocate(H)
    deallocate(b)
    deallocate(CI)
    deallocate(w)
    call destruct(bst)

    if(N_oneel .ne. Nwf) then
       print*,'N_oneel=',N_oneel, ', Nwf=', Nwf
       stop 'They should be equal!!!'
    endif

    call sort_by_energy(self)
!!$    call print_energy(self)

  end subroutine construct_wflocalpot_nr

!
end module sturmian_class



!===================================================================
!                                     m
!   Laguerre's polinomials  f(i,n) = L   (dlambda * grid(i))
!                                     n-1
!===================================================================
subroutine lagpol8 (m, dlambda, f, nload, nmax, grid, nr)
  implicit none
  integer:: nload, nmax, nr
  real*8:: m, dlambda, f(nload, nmax), grid(nr)
  real*8:: pl(nmax)
  real*8:: L2, r, x, pl1, pl2, pl3
  integer:: i, n
!
! INPUT:
! -----
!  m - parameter of Laguerre polinomial.
!  nmax  - the max number of n.
!  nload - dimension of f.
!  grid  - "R" - grid of "NR" points.
!  pl(nmax) - work space for this program.
!
! OUTPUT:
! ------
!  f(i,n) - Laguerre polinomials.
!
  L2 = m - 2
!
!     Loop by  r-grid
!
  do i = 1, nr
     r = grid(i)
     x = dlambda * r
!
!        Define Laguerre's polinomials
!        and store them in array 'f'
!
     pl1    = 1.0d0
     pl(1)  = pl1
     f(i,1) = pl1
!
     pl2    = dble(L2 + 3) - x
     if (nmax.gt.1) then
        pl(2)  = pl2
        f(i,2) = pl2
     endif
!
     do n = 3, nmax
        pl3 = ((dble(2*n-1+L2)-x)*pl2 - dble(n+L2)*pl1) / dble(n-1)
        pl(n) = pl3
        f(i,n) = pl3
        pl1 = pl2
        pl2 = pl3
     end do
  end do
  return
end subroutine lagpol8
!
!==================================================================
!
subroutine minmaxi(f,nr,i1,i2)
  use grid_radial
  implicit none
  integer, intent(in):: nr
  real*8, intent(in), dimension(nr)::  f
  integer, intent(out):: i1, i2

  integer:: i

  !< (tom ross) case where f(x) = 0
  !< modified i.gt.1 -> i.gt.i1
  i=1
  do while (i.lt.nr.and.abs(f(i)).lt.grid%regcut)
     i=i+1
  end do
  i1=i
  i=nr
  do while (i.gt.i1.and.abs(f(i)).lt.grid%expcut)
     i=i-1
  end do
  i2=i

  return
end subroutine minmaxi
!
subroutine minmaxi2(f,nr,i1,i2)
  use grid_radial
  implicit none
  integer, intent(in):: nr
  real*8, intent(in), dimension(nr,2)::  f
  integer, intent(out):: i1, i2

  integer:: i, i1s,i1l,i2s,i2l

  i=1
  do while (i.lt.nr.and.abs(f(i,1)).lt.grid%regcut)
     i=i+1
  end do
  i1l=i
  i=nr
  do while (i.gt.1.and.abs(f(i,1)).lt.grid%expcut)
     i=i-1
  end do
  i2l=i

  i=1
  do while (i.lt.nr.and.abs(f(i,2)).lt.grid%regcut)
     i=i+1
  end do
  i1s=i
  i=nr
  do while (i.gt.1.and.abs(f(i,2)).lt.grid%expcut)
     i=i-1
  end do
  i2s=i

  i1 = min(i1s,i1l)
  i2 = max(i2s,i2l)

  return
end subroutine minmaxi2
!
!==================================================================
!
!==================================================================
