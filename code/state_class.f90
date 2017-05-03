module state_class

   private
  type, public:: state
     character*4 :: label ! Label representing the quantum numbers (e.g. ground state is 1sSg).
     real*8:: M ! total angular momentum projection of the state
     integer:: parity
     real*8:: spin  ! spin = 0,1 for two-electron targets, 1/2 for one-electron targets
     real*8:: energy
     real*8, dimension(:), allocatable:: CI
     integer, dimension(:), allocatable:: na
     integer, dimension(:), allocatable:: ma ! magnetic sublevels
     integer:: nam   !  size of arrays na(:), ma(:), CI(:)

     integer:: l   ! major config l value
     integer:: n   ! major config n value

     integer:: inum ! state index for given target symmetry

  end type state
!
! NOTE that JS has added some variables to both classes (label and n for state, l for basis_state) without updating the methods. If things aren't working (such as states not copying properly), this is why. Come see me !
!
  type, public:: basis_state
     real*8:: Mmax
     type(state), dimension(:), pointer:: b  ! array of states
     integer:: Nmax    ! number of states (size of the basis)
     integer:: Nstates    ! last initialized state
     integer:: nicm  ! number of core orbitals, set in H12.f90
     integer, dimension(:), allocatable:: ncore, mncore    ! index to core orbitals, allocated in H12.f90
     real*8:: en_ion   ! energy of one-el. ion ground state
     integer:: Nmax_bound

  end type basis_state


  public::  new_basis_st, destruct_basis_st, sort_by_energy_basis_st, &
      basis_size_st, print_energy_basis_st, calc_spectro_factors
  public:: setupCI, construct_st, destruct_st, copy_st, get_angmom_proj, &
      get_par_st, get_energy_st, get_inum_st, get_nam, get_na, get_ma, get_CI, &
      ovlp_st, H1el_st,  osc_st, modify_CI_state, get_l_majconf, set_l_majconf


  interface get_nam
     module procedure get_nam_st
  end interface

  interface get_na
     module procedure get_na_st
  end interface

  interface get_ma
     module procedure get_ma_st
  end interface

  interface print_energy
     module procedure print_energy_basis_st
  end interface

 interface new_basis_st
     module procedure new_basis_st
  end interface


contains
!
!
  subroutine  construct_st(self,m,parity,spin,energy,inum,ncm,CI,no1,mo1)

   implicit none

    type(state), intent(inout):: self
    real*8, intent(in):: m
    integer, intent(in):: parity
    real*8, intent(in):: spin
    real*8, intent(in):: energy
    integer, intent(in):: inum
    integer, intent(in):: ncm
    real*8, dimension(ncm), intent(in) :: CI
    integer, dimension(ncm), intent(in) :: no1, mo1

    self%label = "  "
    self%m = m
    self%parity = parity
    self%spin = spin
    self%energy = energy

    self%l = -1
    self%n = -1

    self%inum = inum

    allocate(self%CI(ncm),self%na(ncm),self%ma(ncm))
    self%nam=ncm
    self%na(1:ncm) = no1(1:ncm)
    self%ma(1:ncm) = mo1(1:ncm)
    self%CI(1:ncm) = CI(1:ncm)

!    print*,'created a state with energy:', energy

  end   subroutine  construct_st
  !
  !  only for one-electron states
  subroutine modify_CI_state(self,nam_in,na_in,C_in, ma_in)

    implicit none

    type(state), intent(inout):: self
    integer, intent(in):: nam_in
    integer, dimension(nam_in), intent(in):: na_in
    real*8, dimension(nam_in), intent(in):: C_in
    integer, dimension(nam_in), intent(in):: ma_in
    integer:: nam

    if(nam_in .ne. self%nam) then
       if(allocated(self%CI)) then
          deallocate(self%CI,self%na, self%ma)
       endif
       nam = nam_in
       allocate(self%CI(nam),self%na(nam),self%ma(nam))
    endif
    self%nam = nam_in
    self%na(:) = na_in(:)
    self%CI(:) = C_in(:)
    self%ma(:) = ma_in(:)

  end subroutine modify_CI_state
!
! Shellsort algorithm for integer array na(:)
  subroutine sort_by_value(N,na)
    implicit none

    integer, intent(in):: N
    integer, dimension(N), intent(inout):: na
    integer:: gap, i, j
    integer:: Tmp

    gap = N/2
    do
       if(gap .le. 0) exit
       do i=gap+1,N
          Tmp = na(i)
          do j=i,gap+1,-gap
             if(Tmp.lt. na(j-gap)) then
                na(j) = na(j-gap)
             else
                exit
             endif
             na(j-gap) = Tmp
          enddo
       enddo
       if ( gap .eq. 2) then
          gap = 1
       else
          gap = gap/2.2
       endif
    enddo

  end subroutine sort_by_value
  !
  !
  subroutine  destruct_st(self)
    implicit none
    type(state), intent(inout):: self

    if(allocated(self%CI)) deallocate(self%CI)
    if(allocated(self%na)) deallocate(self%na)
    if(allocated(self%ma)) deallocate(self%ma)

  end subroutine destruct_st
  !
  subroutine copy_st(state_l,state_r)
    implicit none

    type(state), intent(out):: state_l
    type(state), intent(in):: state_r
    integer:: i1,i2, nam

    nam = state_r%nam

    state_l%energy = state_r%energy
    state_l%M = state_r%M
    state_l%parity = state_r%parity
    state_l%spin = state_r%spin

    state_l%l = state_r%l
    state_l%n = state_r%n

    state_l%inum = state_r%inum

    if(allocated(state_l%CI)) then
       deallocate(state_l%CI,state_l%na,state_l%ma)
    endif
    state_l%nam = state_r%nam

    allocate(state_l%na(nam))
    allocate(state_l%CI(nam))
    allocate(state_l%ma(nam))

    state_l%CI = state_r%CI
    state_l%na = state_r%na
    state_l%ma = state_r%ma

  end subroutine copy_st
!
  function get_energy_st(self)
    implicit none
    real*8:: get_energy_st
    type(state), intent(in):: self

    get_energy_st = self%energy
  end function get_energy_st
!
!
  function get_angmom_proj(self)
    implicit none
    real*8:: get_angmom_proj
    type(state), intent(in):: self

    get_angmom_proj = self%M

  end function get_angmom_proj
!
!
  function get_par_st(self)
    implicit none
    integer:: get_par_st
    type(state), intent(in):: self

    get_par_st = self%parity

  end function get_par_st
!
!
  function get_inum_st(self)
    implicit none
    integer:: get_inum_st
    type(state), intent(in):: self

    get_inum_st = self%inum

  end function get_inum_st
!
!
  function get_l_majconf(self)
    implicit none
    integer:: get_l_majconf
    type(state), intent(in):: self

    get_l_majconf = self%l

  end function get_l_majconf
!
!
  subroutine set_l_majconf(self,l)
    implicit none
    type(state):: self
    integer, intent(in):: l

    self%l = l

  end subroutine set_l_majconf
!
!
  function get_nam_st(self)
    implicit none
    integer:: get_nam_st
    type(state), intent(in):: self

    get_nam_st = self%nam

  end function get_nam_st
!
  function get_na_st(self,n)
    implicit none
    integer:: get_na_st
    type(state), intent(in):: self
    integer, intent(in):: n


    if(n .le. self%nam .and. n.ge. 1) then
        get_na_st = self%na(n)
    else
       print*,'Error: state_class.f90, get_na_st(self,n): value of n is out of bounds: n=', n
       print*, 'self%nam=', self%nam
       stop
    endif

  end function get_na_st
!
!
  function get_ma_st(self,n)
    implicit none
    integer:: get_ma_st
    type(state), intent(in):: self
    integer, intent(in):: n

    if(n .le. self%nam .and. n.ge. 1) then
        get_ma_st = self%ma(n)
    else
       print*,'Error: state_class.f90, get_ma_st(self,n): value of n is out of bounds: n=', n
       stop
    endif

  end function get_ma_st
!
  function get_CI(self,i)
    implicit none
    real*8:: get_CI
    type(state), intent(in):: self
    integer, intent(in):: i

    get_CI = self%CI(i)

  end function get_CI
!
!
  subroutine new_basis_st(self,n)
    implicit none

    type(basis_state), intent(inout):: self
    integer, intent(in):: n  ! number of states
    !

    self%Nstates = 0
    self%Nmax = n
    ! create array of n states
    if(n.ne. 0) allocate( self%b(n))

  end subroutine new_basis_st

 subroutine destruct_basis_st(self)
    implicit none
    type(basis_state), intent(inout):: self
    integer:: n
    integer:: i
    integer:: stat

    n= self%Nmax
    do i=1,n
       !       print*,'dealocate: i=', i
       call  destruct_st(self%b(i))
    enddo
    deallocate(self%b, STAT=stat)
  end subroutine destruct_basis_st
!
!
  function basis_size_st(self)
    implicit none
    integer:: basis_size_st
    type(basis_state), intent(in):: self

    basis_size_st = self%Nmax

  end function basis_size_st
!
!
  subroutine  print_energy_basis_st(self)
    use data_targets

    implicit none
    type(basis_state), intent(in):: self
    real*8:: en_ion
    integer:: nc
    real*8:: ioniz_en, ioniz_en_au, exit_en, two_electron_en

    write(*,'("#********* Energy levels: **********")')

!    print*, 'hlike=',self%hlike

    en_ion = 0d0

    write(*,'(4X,"N    m   par state label    excitation energy(eV) ionization energy(au) ionization energy(eV)" )')
    do nc=1,self%Nmax
        exit_en =  real(27.2114*( get_energy_st(self%b(nc)) -  get_energy_st(self%b(1))))
        ioniz_en_au = real(( get_energy_st(self%b(nc)) - en_ion))
        ioniz_en = real(27.2116 * ioniz_en_au )

        write(*,'(3I5,A10,4F22.5)') nc, int(self%b(nc)%M), self%b(nc)%label, exit_en, ioniz_en_au, ioniz_en

    end do

    print*

  end subroutine print_energy_basis_st
  !
  ! Shellsort algorithm
  subroutine sort_by_energy_basis_st(self)
    implicit none

    type(basis_state), intent(inout):: self
    integer:: gap, i, j, N
    type(state):: Tmp

    N = self%Nmax
    gap = N/2
    do
       if(gap .le. 0) exit
       do i=gap+1,N
          call copy_st(Tmp,self%b(i))
          do j=i,gap+1,-gap
             if(Tmp%energy .lt. self%b(j-gap)%energy) then
                call copy_st(self%b(j),self%b(j-gap))
             else
                exit
             endif
             call copy_st(self%b(j-gap),Tmp)
          enddo
       enddo
       if ( gap .eq. 2) then
          gap = 1
       else
          gap = gap/2!.2
       endif
    enddo

  end subroutine sort_by_energy_basis_st
  !
  !
  !
  !
  function make_label(n,l,m,par)
    !
    ! Converts a set of quantum numbers (n,l,m,parity) into a 4-character label.
    !
    ! First is n; second is l (in lowercase s,p,d,f notation);
    ! third is m (in uppercase S,P,D,F notation to symbolise sigma,pi,delta,phi); and
    ! fourth is parity (g = gerade = even = +1, u = ungerade = odd = -1).
    !
    ! For example, 1sSg is the ground state. The next energetic is 2pSu.
    !
    ! Note that the absolute value of m is taken, so m-degenerate states are labelled the same.
    !
    !
    ! JS 6/12/11
    !
    implicit none

    character*4 :: make_label
    integer, intent(in) :: n, l, m, par

!    write(make_label,'(I1)') n
    make_label(1:2) = char(n+48)   ! 48 -> 0, 49 -> 1, etc. Google "ASCII table" if you want.

    select case(l)
       case(0)
          make_label(2:3) = 's'
       case(1)
          make_label(2:3) = 'p'
       case(2)
          make_label(2:3) = 'd'
       case default
          make_label(2:3) = char(l+99)   ! 102 -> f, 103 -> g, etc. Meaningless beyond l=12 or so.
    end select

    select case( abs(m) )   ! Degenerate +/- m states are treated as the same.
       case(0)
          make_label(3:4) = 'S'
       case(1)
          make_label(3:4) = 'P'
       case(2)
          make_label(3:4) = 'D'
       case default
          make_label(3:4) = char(abs(m)+67)   ! 70 -> F, 71 -> G, etc. Meaningless beyond m=12.
    end select

    if ( par == 1 ) then
       make_label(4:) = 'g'
    elseif ( par == -1 ) then
       make_label(4:) = 'u'
    endif


  end function make_label


end module state_class
