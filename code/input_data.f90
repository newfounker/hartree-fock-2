module input_data

  ! forked from ~dmitry/Structure_code_Mol/input_data.f90 by Dmitry Fursa
  ! modified 09/01/2017 by Tom Ross
  ! assuming atomic units
  !     length: bohr
  !     energy: hartree

  public:: readin

  type input
      character*10     :: target           ! label of target

      real*8           :: Z1, Z2, Zasym    ! diatomic charges, total charges
      real*8           :: Rd               ! inter-nuclear distance
      real*8           :: origin           ! 0: on Z1, 0.5: middle, 1: on Z2

      integer          :: labot, latop     ! min, max atomic ang. mom. (l)
      integer, pointer :: nps(:)           ! number of basis functions per l
      real*8,  pointer :: alpha(:)         ! exponential fall-off for each l

      real*8           :: rmax             ! max radial value in grid
      real*8           :: qmax             ! max momentum integrable over grid
      integer          :: ndouble          ! number of doubling
      integer          :: npdbl            ! points per tier in grid
      integer          :: npwave           ! points per oscillation
      integer          :: ltmax            ! max l in v(1, 2) expansion
      real*8           :: formcut
      real*8           :: regcut
      real*8           :: expcut

      integer          :: Mt_min , Mt_max  ! min, max values of Mt
      integer, pointer :: nst(:, :)        ! 1e states per m / +- parity

      integer          :: n_e              ! number of electrons in molecule
      integer          :: iter_max         ! maximum number of hf iterations
      real*8           :: tolerance        ! convergence tolerance of density matrix

  end type input

  type(input):: data_in     ! contains all input data

contains

  subroutine readin(self, input_dir, lwrite)
    implicit none
    type(input)                     :: self
    character(len = *) , intent(in) :: input_dir
    logical            , intent(in) :: lwrite      ! T: write to screen, F: dont
    character(len = 1000)           :: filepath
    logical                         :: file_exists
    integer                         :: unitno
    integer                         :: io_stat
    integer                         :: l, m

    filepath = input_dir//"/data.in"

    inquire(file = filepath, exist = file_exists)

    if (.not. file_exists) then

      write (*, *) filepath, " does not exist"

      stop

    end if

    !< open file
    unitno = 1000
    write (*, *) "opening ", filepath

    open(unitno, file = filepath, action = 'read', iostat = io_stat)

    if (io_stat /= 0) then

      write (*, *) filepath, " cannot be opened"

      stop

    end if

    write (*, '(a)') 'reading data.in'

    ! reading target
    read (unitno, '(a10)') self%target

    if (lwrite) then
      write (*, '(a, a)') ' target: ', self%target
      write (*, *)
    end if

    ! reading charges
    read (unitno, *) self%Z1, self%Z2, self%Zasym

    if (lwrite) then
      write (*, '(a, f10.5)') ' Z1:    ', self%Z1
      write (*, '(a, f10.5)') ' Z2:    ', self%Z2
      write (*, '(a, f10.5)') ' Zasym: ', self%Zasym
      write (*, *)
    end if

    ! reading inter-nuclear distance, origin
    read (unitno, '(2f10.5)') self%Rd, self%origin

    if (lwrite) then
      write (*, '(a, f10.5)') ' inter-nuclear distance: ', self%Rd
      write (*, '(a, f10.5)') ' origin:                 ', self%origin
      write (*, *)
    end if

    ! reading l bottom and top
    read (unitno, '(2i6)') self%labot, self%latop

    if (lwrite) then
      write (*, '(a, i4)') ' l bottom: ', self%labot
      write (*, '(a, i4)') ' l top:    ', self%latop
      write (*, *)
    end if

    allocate(self%alpha(0:self%latop))
    allocate(self%nps(0:self%latop))

    ! reading nps(:), alpha(:)
    read(unitno,*) (self%nps(l), self%alpha(l), l = self%labot, self%latop)

    if (lwrite) then
      write (*,'(a, 100(i5,f10.4,1X))') ' nps(:), alpha(:):', (self%nps(l), self%alpha(l), l = self%labot, self%latop)
      write (*, *)
    end if

    ! reading rmax, qmax, ndouble, npdbl, npwave, ltmax
    read (unitno, *) self%rmax, self%qmax, self%ndouble, self%npdbl, self%npwave, self%ltmax

    if (lwrite) then
      write (*, '(a, f10.5)') ' rmax:    ', self%rmax
      write (*, '(a, f10.5)') ' qmax:    ', self%qmax
      write (*, '(a, i4)')    ' ndouble: ', self%ndouble
      write (*, '(a, i4)')    ' npdbl:   ', self%npdbl
      write (*, '(a, i4)')    ' npwave:  ', self%npwave
      write (*, '(a, i4)')    ' ltmax:   ', self%ltmax
      write (*, *)
    end if

    ! reading formcut, regcut, expcut
    read (unitno, '(3e10.2)') self%formcut, self%regcut, self%expcut

    if (lwrite) then
      write (*, '(a, es10.3)') ' formcut: ', self%formcut
      write (*, '(a, es10.3)') ' regcut:  ', self%regcut
      write (*, '(a, es10.3)') ' expcut:  ', self%expcut
      write (*, *)
    end if

    ! reading min, max values of Mt
    read (unitno, '(2i6)') self%Mt_min, self%Mt_max

    if (lwrite) then
      write (*, '(a, i4)') ' Mt_min: ', self%Mt_min
      write (*, '(a, i4)') ' Mt_max: ', self%Mt_max
      write (*, *)
    end if

    allocate(self%nst(0:self%Mt_max, -1:1))

    self%nst(:, :) = 0

    ! reading number of 1e states per m
    read (unitno, *) self%nst(self%Mt_min:self%Mt_max,  1)
    read (unitno, *) self%nst(self%Mt_min:self%Mt_max, -1)

    if (lwrite) then
      write (*, '(a, 100i4)') &
          ' nst(m, +1):', self%nst(self%Mt_min:self%Mt_max,  1)
      write (*, '(a, 100i4)') &
          ' nst(m, -1):', self%nst(self%Mt_min:self%Mt_max, -1)
      write (*, *)
    end if

    ! reading n_e, iter_max, tolerance
    read (unitno, *) self%n_e, self%iter_max, self%tolerance

    if (lwrite) then
      write (*, '(a, i4)')     ' n_e:         ', self%n_e
      write (*, '(a, i4)')     ' iter_max:    ', self%iter_max
      write (*, '(a, es10.3)') ' tolerance:   ', self%tolerance
      write (*, *)
    end if

    write (*, *)

  end subroutine readin

end module input_data
