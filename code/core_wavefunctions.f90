module core_wavefunctions

!> - coulomb_pw -> form subroutine

!> - symmetry check in spectrum
!> - in potential_me adjust for nuclear term
!> - modify for multiple electronic orbitals
!> - draw potential energy curve

  use input_data
  use grid_radial
  use sturmian_class
  use vnc_module
  use target_states
  implicit none

  private
  public core_state

  !< core state
  !<  pw(l, s) is the l-th partial wave (radial function) of the s-th spatial
  !<    orbital (there are (n_e + 1)/2 spatial orbitals.)
  !<  coulomb_pw(:, l) is the l-th partial wave of the e-e (electron-electron)
  !<    coulomb potential
  !<  potential_pw(:, l) is the l-th partial wave of the e-e and nuclear
  !<    potentials summed together
  type core_state
    private

    integer                         :: m
    integer                         :: n_e
    integer                         :: labot, latop
    real*8                          :: radial_distance
    real*8                          :: electronic_energy
    real*8                          :: nuclear_energy
    type(sturmian_nr) , allocatable :: pw(:, :)

    real*8            , allocatable :: coulomb_pw(:, :)
    real*8            , allocatable :: potential_pw(:, :)

  contains

    procedure , pass :: construct_core_state

    procedure , pass :: get_radial_distance
    procedure , pass :: get_electronic_energy
    procedure , pass :: get_nuclear_energy

    procedure , pass :: get_pw_functions
    procedure , pass :: get_coulomb_pw
    procedure , pass :: get_potential_pw

    procedure , pass :: coulomb_me
    procedure , pass :: exchange_me
    procedure , pass :: potential_me

    procedure , pass :: spectrum

    procedure , pass :: read_from
    procedure , pass :: write_to

    procedure , pass :: write_pw_to
    procedure , pass :: write_orbital_pw_to
    procedure , pass :: write_coulomb_pw_to
    procedure , pass :: write_potential_pw_to

  end type core_state

contains

!> construct_core_state
  subroutine construct_core_state (core, m, n_e, labot, latop, &
      radial_distance, electronic_energy, nuclear_energy, pw)
    class(core_state) , intent(inout) :: core
    integer           , intent(in)    :: m
    integer           , intent(in)    :: n_e
    integer           , intent(in)    :: labot, latop
    real*8            , intent(in)    :: radial_distance
    real*8            , intent(in)    :: electronic_energy
    real*8            , intent(in)    :: nuclear_energy
    real*8            , intent(in)    :: pw(:, labot:, :)
    integer                           :: s, l, kk, lambda
    integer                           :: mini, maxi

    write (*, "(a)") &
        "> core-wavefunctions :: &
        constructing"

    !< record state information
    core%m = m
    core%n_e = n_e
    core%labot = labot
    core%latop = latop
    core%radial_distance = radial_distance
    core%electronic_energy = electronic_energy
    core%nuclear_energy = nuclear_energy

    !< orbital partial waves
    allocate(core%pw(core%labot:core%latop, 1:(core%n_e + 1) / 2))

    grid%nr = size(pw(:, :, :), 1)

    kk = 1
    do s = 1, (core%n_e + 1) / 2

      do l = core%labot, core%latop

        call minmaxi(pw(:, l, s), grid%nr, mini, maxi)

        call init_function(core%pw(l, s), l, core%m, kk, mini, maxi, &
            pw(:, l, s), grid%nr)

        kk = kk + 1

      end do

    end do

    !< examine spectroscopic factors
    call calc_spectroscopic_factors(core)

    !< coulomb partial waves
    allocate(core%coulomb_pw(1:grid%nr, 0:2*core%latop))
    call calc_coulomb_potential(core)

    !< potential (e-e coulomb + nuclear coulomb) partial waves
    allocate(core%potential_pw(1:grid%nr, 0:max(ubound(vnc, 2), 2*core%latop)))

    core%potential_pw(:, :) = 0.0

    do lambda = lbound(vnc, 2), ubound(vnc, 2)

      core%potential_pw(:, lambda) = core%potential_pw(:, lambda) + vnc(:, lambda)

    end do

    do lambda = 0, 2*core%latop

      core%potential_pw(:, lambda) = core%potential_pw(:, lambda) + core%coulomb_pw(:, lambda)

    end do

  end subroutine construct_core_state

!> get_radial_distance
!>   energy associated with frozen core
  function get_radial_distance (core) result (radial_distance)
    class(core_state) , intent(in) :: core
    real*8                         :: radial_distance

    radial_distance = core%radial_distance

  end function get_radial_distance

!> get_electronic_energy
!>   energy associated with frozen core
  function get_electronic_energy (core) result (electronic_energy)
    class(core_state) , intent(in) :: core
    real*8                         :: electronic_energy

    electronic_energy = core%electronic_energy

  end function get_electronic_energy

!> get_nuclear_energy
!>   energy associated with frozen core
  function get_nuclear_energy (core) result (nuclear_energy)
    class(core_state) , intent(in) :: core
    real*8                         :: nuclear_energy

    nuclear_energy = core%nuclear_energy

  end function get_nuclear_energy

!> get_pw_functions
  subroutine get_pw_functions (core, pw)
    class(core_state) , intent(in)  :: core
    real*8            , intent(out) :: pw(:, core%labot:, 1:)
    integer                         :: mini, maxi
    real*8            , allocatable :: temp(:)
    integer                         :: s, l

    allocate(temp(1:grid%nr))

    do s = 1, (core%n_e + 1) / 2

      do l = core%labot, core%latop

        mini = get_minf(core%pw(l, s))
        maxi = get_maxf(core%pw(l, s))

        temp(:) = fpointer(core%pw(l, s))

        pw(:, l, s) = 0d0
        pw(mini:maxi, l, s) = temp(mini:maxi)

      end do

    end do

  end subroutine get_pw_functions

!> get_coulomb_pw
  function get_coulomb_pw (core) result (coulomb_pw)
    class(core_state) , intent(in)  :: core
    real*8            , allocatable :: coulomb_pw(:, :)

    allocate(coulomb_pw(1:grid%nr, 0:core%latop))

    coulomb_pw(:, :) = core%coulomb_pw(:, :)

  end function get_coulomb_pw

!> get_potential_pw
  function get_potential_pw (core) result (potential_pw)
    class(core_state) , intent(in)  :: core
    real*8            , allocatable :: potential_pw(:, :)

    allocate(potential_pw(1:grid%nr, 0:core%latop))

    potential_pw(:, :) = core%potential_pw(:, :)

  end function get_potential_pw

!> coulomb_me
!>  e-e coulomb matrix element
  function coulomb_me (core, pi, pj) result (V_ij)
    class(core_state)           , intent(in) :: core
    type(sturmian_nr) , pointer , intent(in) :: pi, pj
    real*8                                   :: V_ij
    real*8            , pointer              :: fi(:), fj(:)
    integer                                  :: li, lj
    integer                                  :: mi, mj
    integer                                  :: mini, maxi
    integer                                  :: lambda
    integer                                  :: lam_mini, lam_maxi
    real*8                                   :: Yint

    fi => fpointer(pi)
    li = get_ang_mom(pi)
    mi = get_ang_mom_proj(pi)

    fj => fpointer(pj)
    lj = get_ang_mom(pj)
    mj = get_ang_mom_proj(pj)

    do lambda = lbound(core%coulomb_pw, 2), ubound(core%coulomb_pw, 2)

      call minmaxi(core%coulomb_pw(:, lambda), size(core%coulomb_pw, 1), &
          lam_mini, lam_maxi)

      mini = max(get_minf(pi), get_minf(pj), lam_mini)
      mini = min(get_maxf(pi), get_maxf(pj), lam_maxi)

      V_ij = sum(fi(mini:maxi) * fj(mini:maxi) * &
          core%coulomb_pw(mini:maxi, lambda) * grid%weight(mini:maxi)) * &
          Yint(&
          dble(li), dble(mi), &
          dble(lambda), dble(0), &
          dble(lj), dble(mj))

    end do

  end function coulomb_me

!> exchange_me
!>  e-e exchange matrix element
  function exchange_me (core, pi, pj) result (K_ij)
    class(core_state)           , intent(in) :: core
    type(sturmian_nr) , pointer , intent(in) :: pi, pj
    real*8                                   :: K_ij
    type(sturmian_nr)                        :: pk, pl
    integer                                  :: mi, mj, mk, ml
    integer                                  :: s, l_1, l_2
    real*8                                   :: temp

    mi = get_ang_mom_proj(pi)
    mj = get_ang_mom_proj(pj)

    K_ij = 0d0

    do s = 1, (core%n_e + 1) / 2

    do l_1 = core%labot, core%latop

      pk = core%pw(l_1, s)
      mk = get_ang_mom_proj(pk)

      do l_2 = core%labot, core%latop

        pl = core%pw(l_2, s)
        ml = get_ang_mom_proj(pl)

        call V12me(pi, pk, pl, pj, mi, mk, ml, mj, temp)

        K_ij = K_ij + temp

      end do

    end do

    end do

  end function exchange_me

!> potential_me
!>  nuclear + e-e coulomb + e-e exchange
  function potential_me (core, pi, pj) result (V_ij)
    class(core_state)           , intent(in) :: core
    type(sturmian_nr) , pointer , intent(in) :: pi, pj
    real*8                                   :: V_ij
    real*8                                   :: VLambdaR_ME

    V_ij = core%coulomb_me(pi, pj) + core%exchange_me(pi, pj)
        ! + VLambdaR_ME(pi, pj, core%m)

  end function potential_me

!> diagonalise a basis with frozen core electron potentials
!>  sa_i(:) are the indexes of the symmetry-adapted basis functions (w.r.t. the
!>    larger basis)
!>  sa_n is the number of symmetry-adapted basis functions
  subroutine spectrum (core, basis, sa_n, sa_i, T, S)
    class(core_state)       , intent(in)  :: core
    type(basis_sturmian_nr) , intent(in)  :: basis
    integer                 , intent(in)  :: sa_n
    integer                 , intent(in)  :: sa_i(:)
    real*8                  , intent(in)  :: T(:, :)
    real*8                  , intent(in)  :: S(:, :)
    real*8                  , allocatable :: V(:, :)
    real*8                  , allocatable :: C(:, :), w(:)
    type(sturmian_nr)       , pointer     :: pi, pj
    integer                               :: ierr
    integer                               :: ii, jj

    write (*, *) "diagonalising with nuclear, coulomb and exchange potentials"

    !< (coulomb + exchange) potential matrix
    allocate(V(1:sa_n, 1:sa_n))
    V(:, :) = 0d0

    do ii = 1, sa_n

      do jj = 1, sa_n

        pi => basis%b(sa_i(ii))
        pj => basis%b(sa_i(jj))

        V(ii, jj) = core%potential_me(pi, pj)

      end do

    end do

    !< diagonalise
    allocate(C(1:sa_n, 1:sa_n))
    allocate(w(1:sa_n))

    call rsg(sa_n, sa_n, T + V, S, w, 2, C, ierr)

    if (ierr /= 0) then
      write (*, "(a)") "rsg failed to diagonalise the system"
    end if

    !< write eigenvalues
    do ii = 1, sa_n
      write (*, '(i4, f10.3)') ii, w(ii)
    end do

    write (*, *)

  end subroutine spectrum

!> read_from
!>  Reads core radial wavefunctions (partial wave expansions) from a given file.
!>  The wave functions are then interpolated to be plotted on a new
!>  grid.
!>  The partial wave pw(:, l, s) is the l-th partial wave of the s-th spatial
!>  orbiatl.
  subroutine read_from (core, filepath)
    class(core_state)  , intent(inout) :: core
    character(len = *) , intent(in)    :: filepath
    integer                            :: unitno, io_stat
    logical                            :: file_exists
    integer                            :: temp_nr
    real*8             , allocatable   :: temp_grid(:)
    real*8             , allocatable   :: temp_pw(:, :, :)
    real*8             , allocatable   :: intpl_pw(:, :, :)
    integer                            :: mini, maxi
    logical                            :: grid_cutoff_located
    integer                            :: grid_cutoff
    integer                            :: s, l, lambda, ii, kk

    write (*, "(a)") &
        "> core-wavefunctions :: &
        reading"

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

    !< read file
    write (*, *) "reading ", filepath

    !< read state information
    read (unitno, *) core%m
    read (unitno, *) core%n_e
    read (unitno, *) core%labot, core%latop
    read (unitno, *) core%radial_distance
    read (unitno, *) core%electronic_energy
    read (unitno, *) core%nuclear_energy

    !< read original grid
    read (unitno, *) temp_nr

    allocate(temp_grid(1:temp_nr))

    read (unitno, *) temp_grid(1:temp_nr)

    !< locate point where original grid ends
    grid_cutoff = grid%nr
    grid_cutoff_located = .false.
    ii = 1

    do while ((.not. grid_cutoff_located) .and. (ii < grid%nr))

      if (grid%gridr(ii) >= temp_grid(temp_nr)) then

        grid_cutoff = ii

        grid_cutoff_located = .true.

      else

        ii = ii + 1

      end if

    end do

    !< allocate partial wave arrays
    allocate(temp_pw(1:temp_nr, core%labot:core%latop, 1:(core%n_e + 1) / 2))
    allocate(intpl_pw(1:grid%nr, core%labot:core%latop, 1:(core%n_e + 1) / 2))

    !< partial waves (plotted on original grid)
    do s = 1, (core%n_e + 1) / 2

      do l = core%labot, core%latop

        read (unitno, *) temp_pw(:, l, s)

      end do

    end do

    !< plot partial waves on current grid by interpolating
    do s = 1, (core%n_e + 1) / 2

      do l = core%labot, core%latop

        call intrpl(temp_nr, temp_grid, temp_pw(:, l, s), &
            grid_cutoff, grid%gridr(1:grid_cutoff), &
            intpl_pw(1:grid_cutoff, l, s))

        intpl_pw(grid_cutoff:grid%nr, l, s) = 0.0

      end do

    end do

    !< store partial waves as sturmian_nr type
    allocate(core%pw(core%labot:core%latop, 1:(core%n_e + 1) / 2))

    kk = 1
    do s = 1, (core%n_e + 1) / 2

      do l = core%labot, core%latop

        call minmaxi(intpl_pw(:, l, s), grid%nr, mini, maxi)
        ! write (*, '(2i2, 2i6)') s, l, mini, maxi

        call init_function(core%pw(l, s), l, core%m, kk, mini, maxi, &
            intpl_pw(:, l, s), grid%nr)

        kk = kk + 1

      end do

    end do

    !< close file
    write (*, *) "closing ", filepath
    close(unitno)

    !< empty line
    write (*, *)

    !< examine spectroscopic factors
    call calc_spectroscopic_factors(core)

    !< coulomb partial waves
    allocate(core%coulomb_pw(1:grid%nr, 0:2*core%latop))
    call calc_coulomb_potential(core)

    !< potential (e-e coulomb + nuclear coulomb) partial waves
    allocate(core%potential_pw(1:grid%nr, 0:max(ubound(vnc, 2), 2*core%latop)))

    core%potential_pw(:, :) = 0.0

    do lambda = lbound(vnc, 2), ubound(vnc, 2)

      core%potential_pw(:, lambda) = core%potential_pw(:, lambda) + vnc(:, lambda)

    end do

    do lambda = 0, 2*core%latop

      core%potential_pw(:, lambda) = core%potential_pw(:, lambda) + core%coulomb_pw(:, lambda)

    end do

  end subroutine read_from


!> write_to
!>  writes core state information to a file (same form as it would expect to
!>  read_from).
  subroutine write_to (core, filepath)
    class(core_state)  , intent(inout) :: core
    character(len = *) , intent(in)    :: filepath
    real*8             , allocatable   :: pw(:, :, :)
    integer                            :: unitno
    integer                            :: s, l

    write (*, "(a)") &
        "> core-wavefunctions :: &
        writing"

    !< open file
    unitno = 1000
    write (*, *) "opening ", filepath

    open(unitno, file = filepath)

    !< write to file
    write (*, *) "writing ", filepath

    !< write state information
    write (unitno, *) core%m
    write (unitno, *) core%n_e
    write (unitno, *) core%labot, core%latop
    write (unitno, *) core%radial_distance
    write (unitno, *) core%electronic_energy
    write (unitno, *) core%nuclear_energy

    !< write grid
    write (unitno, *) grid%nr
    write (unitno, *) grid%gridr(:)

    !< write orbital partial waves
    allocate(pw(1:grid%nr, core%labot:core%latop, 1:(core%n_e + 1) / 2))
    call core%get_pw_functions(pw)

    do s = 1, (core%n_e + 1) / 2

      do l = core%labot, core%latop

        write (unitno, *) pw(:, l, s)

      end do

    end do

    !< close file
    write (*, *) "closing ", filepath
    close(unitno)

    !< empty line
    write (*, *)

  end subroutine write_to

!> write_pw_to
!>  calls write_pw_to, write_coulomb_pw_to, write_potential_pw_to, writing to
!>  data files in a given directory
  subroutine write_pw_to (core, directory)
    class(core_state)  , intent(inout) :: core
    character(len = *) , intent(in)    :: directory

    call core%write_orbital_pw_to(directory//"core_plots.dat")
    call core%write_coulomb_pw_to(directory//"core_coulomb.dat")
    call core%write_potential_pw_to(directory//"core_potential.dat")

  end subroutine write_pw_to

!> write_orbital_pw_to
!>  write partial waves of core orbitals to a file, in a format suitable for
!>  gnuplot
  subroutine write_orbital_pw_to (core, filepath)
    class(core_state)  , intent(inout) :: core
    character(len = *) , intent(in)    :: filepath
    real*8             , allocatable   :: pw(:, :, :)
    integer                            :: unitno
    integer                            :: s, l, ii

    write (*, "(a)") &
        "> core-wavefunctions :: &
        writing spatial orbital partial waves"

    !< open file
    unitno = 1000
    write (*, *) "opening ", filepath

    open(unitno, file = filepath)

    !< write partial waves to file
    allocate(pw(1:grid%nr, core%labot:core%latop, 1:(core%n_e + 1) / 2))
    call core%get_pw_functions(pw)

    write (*, *) "writing partial waves"
    do ii = 1, grid%nr

      write(unitno, *) grid%gridr(ii), pw(ii, core%labot:core%latop, &
          1:(core%n_e + 1) / 2)

    end do

    !< close file
    write (*, *) "closing ", filepath
    close(unitno)

    !< empty line
    write (*, *)

  end subroutine write_orbital_pw_to

!> write_coulomb_pw_to
!>  write partial waves of core, coulomb potential to a file, in a format
!>  suitable for gnuplot
  subroutine write_coulomb_pw_to (core, filepath)
    class(core_state)  , intent(inout) :: core
    character(len = *) , intent(in)    :: filepath
    integer                            :: unitno
    integer                            :: ii


    write (*, "(a)") &
        "> core-wavefunctions :: &
        writing e-e coulomb potential partial waves"

    !< open file
    write (*, *) "opening ", filepath
    unitno = 1000
    open(unitno, file = filepath)

    !< write partial waves to file
    write (*, *) "writing partial waves"
    do ii = 1, grid%nr

      write(unitno, *) grid%gridr(ii), core%coulomb_pw(ii, 0:2*core%latop)

    end do

    !< close file
    write (*, *) "closing ", filepath
    close(unitno)

    !< empty line
    write (*, *)

  end subroutine write_coulomb_pw_to

!> write_potential_pw_to
!>  write partial waves of core, electron + nuclear potential to a file, in a
!>  format suitable for gnuplot
  subroutine write_potential_pw_to (core, filepath)
    class(core_state)  , intent(inout) :: core
    character(len = *) , intent(in)    :: filepath
    integer                            :: unitno
    integer                            :: ii

    write (*, "(a)") &
        "> core-wavefunctions :: &
        writing (nuclear + e-e coulomb) potential partial waves"

    !< open file
    write (*, *) "opening ", filepath
    unitno = 1000
    open(unitno, file = filepath)

    !< write partial waves to file
    write (*, *) "writing partial waves"
    do ii = 1, grid%nr

      write(unitno, *) grid%gridr(ii), core%potential_pw(ii, 0:2*core%latop)

    end do

    !< close file
    write (*, *) "closing ", filepath
    close(unitno)

    !< empty line
    write (*, *)

  end subroutine write_potential_pw_to

!> calculate spectroscopic factors & print to screen
  subroutine calc_spectroscopic_factors (core)
    type(core_state) , intent(inout) :: core
    real*8           , allocatable   :: spectroscopic(:, :)
    real*8           , allocatable   :: pw(:, :, :)
    integer                          :: s, l

    write (*, "(a)") &
        "> core-wavefunctions :: &
        calculating spectroscopic factors"

    !< extract partial waves
    allocate(pw(1:grid%nr, core%labot:core%latop, 1:(core%n_e + 1) / 2))
    call core%get_pw_functions(pw)

    !< spectroscopic factors
    allocate(spectroscopic(core%labot:core%latop, 1:(core%n_e + 1) / 2))
    spectroscopic(:, :) = 0.0

    do s = 1, (core%n_e + 1) / 2

      write (*, "(a, i2)") &
          " orbtial: ", s

      do l = core%labot, core%latop

        spectroscopic(l, s) = sum(pw(:, l, s) * pw(:, l, s) * grid%weight(:))

          if (abs(spectroscopic(l, s)) > 1.0e-5) then

            write (*, '(i4, f10.5)') l, spectroscopic(l, s)

          end if

      end do

      write (*, '(a, f10.5)') " sum", sum(spectroscopic(:, s))
      write (*, *)

    end do

  end subroutine calc_spectroscopic_factors

!> calc_coulomb_potential
!>  l_t       angular quantum number of target
!>  lambda    angular quantum number summed over in laplace expansion
!>  laplace   radial component of laplace expansion
!>  V_pw      V_pw(:, l) is the l-th partial wave of the core potential
  subroutine calc_coulomb_potential (core)
    type(core_state) , intent(inout) :: core
    real*8           , allocatable   :: pw(:, :, :)
    integer                          :: mini1, maxi1, mini2, maxi2
    integer                          :: minf1, maxf1, minf2, maxf2
    real*8           , allocatable   :: integral(:)
    real*8           , allocatable   :: temp(:)
    integer                          :: lambda, lambda_min, lambda_max
    integer                          :: s, l_1, l_2
    real*8                           :: Yint, harmonic
    integer                          :: ii, jj

    write (*, "(a)") &
        "> core-wavefunctions :: &
        calculating e-e coulomb potential partial waves"

    !< summed over radial potential
    lambda_min = 0
    lambda_max = 2 * core%latop

    !< extract partial waves
    allocate(pw(1:grid%nr, core%labot:core%latop, 1:(core%n_e + 1) / 2))
    call core%get_pw_functions(pw)

    !< temporary arrays
    allocate(integral(1:grid%nr))
    allocate(temp(1:grid%nr))

    !< lambda potential
    core%coulomb_pw(:, :) = 0.0

    do s = 1, (core%n_e + 1) / 2

      do lambda = lambda_min, lambda_max

        do l_1 = core%labot, core%latop

          call minmaxi(pw(:, l_1, s), grid%nr, mini1, maxi1)

          do l_2 = core%labot, core%latop

            call minmaxi(pw(:, l_2, s), grid%nr, mini2, maxi2)

            minf1 = max(mini1, mini2)
            maxf1 = min(maxi1, maxi2)

            temp(:) = 0.0
            integral(:) = 0.0

            temp(minf1:maxf1) = pw(minf1:maxf1, l_1, s) * &
                pw(minf1:maxf1, l_2, s) * grid%weight(minf1:maxf1)

            harmonic = Yint(&
                dble(l_1), dble(core%m), &
                dble(lambda), dble(0), &
                dble(l_2), dble(core%m))

            call form(lambda, temp(:), minf1, maxf1, grid%nr, integral(:), minf2, maxf2)

            core%coulomb_pw(minf2:maxf2, lambda) = &
                core%coulomb_pw(minf2:maxf2, lambda) + &
                (2.0 * integral(minf2:maxf2) * harmonic)

          end do

        end do

      end do

    end do

    !< empty line
    write (*, *)

  end subroutine calc_coulomb_potential

end module core_wavefunctions
