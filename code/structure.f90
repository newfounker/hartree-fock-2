module structure

  use input_data
  use grid_radial
  use sturmian_class
  use vnc_module
  use one_electron_func_mod
  use target_states
  use ovlpste1me
  use hartree_fock
  use core_wavefunctions
  implicit none

  private
  public structure_simple, structure_core

contains

!> Performs the Hartree-Fock procedure for the system specified in the global
!> variable data_in.
  subroutine structure_simple ()
    type(basis_sturmian_nr)               :: basis
    real*8                  , allocatable :: H(:, :)
    real*8                  , allocatable :: integrals(:, :, :, :)
    type(core_state)                      :: core
    integer                               :: basis_n, ii

    !< construct basis
    call construct_diagonalised(basis)

    basis_n = basis_size(basis)

    !< construct two-electron integrals
    allocate(integrals(1:basis_n, 1:basis_n, 1:basis_n, 1:basis_n))
    call two_electron_integrals(basis, integrals)

    !< allocate nuclear hamiltonian
    allocate(H(1:basis_n, 1:basis_n))

    !< construct nuclear hamiltonian matrix
    call nuclear_hamiltonian(basis, H)

    !< perform hartree_fock procedure (assuming m = 0)
    call hf_procedure(basis, H, integrals, core)

  end subroutine structure_simple

!> Performs the Hartree-Fock procedure for the system specified in the global
!>  variable data_in, but across a range of radial distance values.
!> Stores the results in an array of core_states, one for each radial distance.
!> Note that this subroutine modifies the value of data_in%Rd.
  subroutine structure_core (core_grid, core_states)
    real*8                  , intent(in)  :: core_grid(:)
    type(core_state)        , intent(out) :: core_states(:)
    type(basis_sturmian_nr)               :: basis
    real*8                  , allocatable :: H(:, :)
    real*8                  , allocatable :: integrals(:, :, :, :)
    integer                               :: basis_n, core_n, ii

    !< construct basis
    call construct_diagonalised(basis)

    basis_n = basis_size(basis)

    !< construct two-electron integrals
    allocate(integrals(1:basis_n, 1:basis_n, 1:basis_n, 1:basis_n))
    call two_electron_integrals(basis, integrals)

    !< allocate nuclear hamiltonian
    allocate(H(1:basis_n, 1:basis_n))

    !< calculate core states for varying radial-separation
    core_n = size(core_grid)

    do ii = 1, core_n

      !< construct nuclear potential (note. alters data_in%Rd)
      call nuclear_potential(core_grid(ii))

      !< construct nuclear hamiltonian matrix
      call nuclear_hamiltonian(basis, H)

      !< perform hartree_fock procedure (assuming m = 0)
      call hf_procedure(basis, H, integrals, core_states(ii))

    end do

    do ii = 1, core_n

      write (*, "(3f10.5)") &
          core_grid(ii), &
          core_states(ii)%get_electronic_energy(), &
          core_states(ii)%get_nuclear_energy()

    end do

  end subroutine structure_core

!> Constructs a diagonalised basis, specified by the global data_in variable,
!> which is diagonalised with regard to the kinetic energy operator.
!> After this basis is constructed, the kinetic energy matrix can be built by
!> T(ii, ii) = get_energy(bst%b(ii))
  subroutine construct_diagonalised (basis)
    type(basis_sturmian_nr) , intent(out) :: basis
    type(basis_sturmian_nr) , pointer     :: basis_sets(:)
    type(sturmian_nr)       , pointer     :: p
    real*8                  , allocatable :: potential(:)
    integer                               :: size, n, l, k, m

    write (*, "(a)") "> establishing basis"

    !< allocate entire basis (indexing over all k, l, m values)
    !< modified for only m = 0 basis functions
    size = 0
    do l = data_in%labot, data_in%latop

      ! size = size + (data_in%nps(l) * ((2 * l) + 1))
      size = size + data_in%nps(l)

    end do

    call new_basis_nr(basis, size)

    !< zero-potential, used in basis diagonalisation
    allocate(potential(1:grid%nr))
    potential = 0.0

    !< entire basis index
    n = 1

    !< loop through basis for given l value
    allocate (basis_sets(data_in%labot:data_in%latop))
    do l = data_in%labot, data_in%latop

      !< construct diagonal basis (w.r.t kinetic energy) with m = 0, for given l
      call construct_wflocalpot_nr(basis_sets(l), data_in%nps(l), potential, &
          data_in%nps(l), l, data_in%alpha(l))

    end do

    !< loop through basis for given l value
    write (*, "(4a4)") "n", "k", "l", "m"

    !< use only m = 0 values
    ! do m = -data_in%latop, data_in%latop, 1
    do m = 0, 0, 1

      do l = abs(m), data_in%latop

        do k = 1, data_in%nps(l)

          call copy(basis%b(n), basis_sets(l)%b(k))

          call set_ang_mom_proj(basis%b(n), m)

          basis%ortint(n, n) = 1.0

          p => basis%b(n)
          write (*, "(4i4)") n, get_k(p), get_ang_mom(p), get_ang_mom_proj(p)

          n = n + 1

        end do

      end do

    end do


    !< destruct temporary basis sets
    do l = data_in%labot, data_in%latop

      call destruct(basis_sets(l))

    end do

    write (*, *)

  end subroutine construct_diagonalised

!> Construct nuclear potential.
!>  Note that data_in%Rd is modified.
  subroutine nuclear_potential (rd)
    real*8  , intent(in) :: rd
    integer              :: lamtop_vc_set

    !< reset data_in%Rd to given value
    data_in%Rd = rd

    !< recalculate nuclear potential partial waves
    call VLambdaR(grid%nr, grid%gridr, data_in%ltmax, data_in%Z1, data_in%Z2, &
        data_in%Rd, data_in%origin, vnc, minvnc, maxvnc, lamtop_vc_set)

  end subroutine nuclear_potential

!> Constructs the hamiltonian (kinetic + nuclear potential) for a given basis,
!> assuming that is has already been diagonalised with regard to the kinetic
!> energy operator.
  subroutine nuclear_hamiltonian (basis, H)
    type(basis_sturmian_nr) , intent(in)  :: basis
    real*8                  , intent(out) :: H(:, :)
    type(sturmian_nr)       , pointer     :: pi, pj
    integer                               :: m
    integer                               :: n, ii, jj
    real*8                                :: VLambdaR_ME

    n = basis_size(basis)

    H(:, :) = 0.0

    !< kinetic
    do ii = 1, n

      pi => basis%b(ii)

      H(ii, ii) = get_energy(pi)

    end do

    do ii = 1, n

      pi => basis%b(ii)

      !< magnetic quantum number
      m = get_ang_mom_proj(pi)

      !< nuclear potential
      do jj = 1, ii

        pj => basis%b(jj)

        if (m == get_ang_mom_proj(pj)) then

          H(ii, jj) = H(ii, jj) + VLambdaR_ME(pi, pj, m)

          !< symmetric matrix
          H(jj, ii) = H(ii, jj)

        end if

      end do

    end do

  end subroutine nuclear_hamiltonian

!> Calculates two-electron integrals for a given basis of functions.
  subroutine two_electron_integrals (basis, integrals)
    type(basis_sturmian_nr) , intent(in)  :: basis
    real*8                  , intent(out) :: integrals(:, :, :, :)
    logical                 , allocatable :: calculated(:, :, :, :)
    type(sturmian_nr)       , pointer     :: pi, pj, pk, pl
    integer                               :: mi, mj, mk, ml
    integer                               :: ii, jj, kk, ll
    integer                               :: n
    !$ real*8                             :: start, finish

    !$ interface
    !$  double precision function omp_get_wtime()
    !$  end function omp_get_wtime
    !$ end interface

    n = basis_size(basis)

    write (*, "(a)") "> two-electron integrals"
    write (*, '(a, es10.3)') &
        " basis size:        ", 1.0 * n
    write (*, '(a, es10.3)') &
        " integrals:         ", (n ** 4) / 1.0
    write (*, '(a, es10.3)') &
        " unique integrals:  ", (n ** 4) / 4.0
    write (*, '(a)') &
        " estimated time (s) "
    write (*, '(a, i6)') &
        "  for 1.0e-7 s/int: ", nint(1.0e-7 * (n ** 4))
    write (*, '(a, i6)') &
        "  for 1.0e-6 s/int: ", nint(1.0e-6 * (n ** 4))
    write (*, '(a, i6)') &
        "  for 1.0e-5 s/int: ", nint(1.0e-5 * (n ** 4))
    write (*, '(a, i6)') &
        "  for 1.0e-4 s/int: ", nint(1.0e-4 * (n ** 4))

    !$ start = omp_get_wtime()

    allocate(calculated(1:n, 1:n, 1:n, 1:n))

    calculated = .false.

    !$omp parallel do &
    !$omp& private(ii, jj, kk, ll, pi, pj, pk, pl) &
    !$omp& shared(basis, integrals, calculated)
    do ii = 1, n

      pi => basis%b(ii)
      mi = get_ang_mom_proj(pi)

      do jj = 1, n

        pj => basis%b(jj)
        mj = get_ang_mom_proj(pj)

        do kk = 1, n

          pk => basis%b(kk)
          mk = get_ang_mom_proj(pk)

          do ll = 1, n

            pl => basis%b(ll)
            ml = get_ang_mom_proj(pl)

            if (.not. calculated(ii, jj, kk, ll)) then

              calculated(ii, jj, kk, ll) = .true.
              calculated(kk, jj, ii, ll) = .true.
              calculated(ii, ll, kk, jj) = .true.
              calculated(kk, ll, ii, jj) = .true.

              call V12me(pi, pj, pk, pl, mi, mj, mk, ml, &
                  integrals(ii, jj, kk, ll))

              integrals(kk, jj, ii, ll) = integrals(ii, jj, kk, ll)
              integrals(ii, ll, kk, jj) = integrals(ii, jj, kk, ll)
              integrals(kk, ll, ii, jj) = integrals(ii, jj, kk, ll)

            end if

          end do

        end do

      end do

    end do
    !$omp end parallel do

    !$ finish = omp_get_wtime()

    !$ write (*, '(a, f10.3)') &
    !$     " time taken:        ", finish - start
    !$ write (*, '(a, es10.3)') &
    !$     " time per integral: ", (finish - start) / (n ** 4)
    write (*, *)

  end subroutine two_electron_integrals

!> Record elements of a basis which satisfy given symmetries of parity and
!> angular momentum projection.
  subroutine symmetry_adaption (basis, parity, m, sa_n, sa_i)
    type(basis_sturmian_nr) , intent(in)  :: basis
    integer                 , intent(in)  :: parity
    integer                 , intent(in)  :: m
    integer                 , intent(out) :: sa_n
    integer       , pointer , intent(out) :: sa_i(:)
    type(sturmian_nr)       , pointer     :: p
    integer                               :: sa_n_temp
    integer                 , allocatable :: sa_i_temp(:)
    integer                               :: n, ii

    !< get basis size
    n = basis_size(basis)

    !< set up temporary symmetry-adapted indexes
    allocate(sa_i_temp(1:n))

    sa_i_temp(:) = 0
    sa_n_temp = 0

    !< record symmetry-adapted basis functions
    do ii = 1, n

      p => basis%b(ii)

      if ((((-1) ** get_ang_mom(p)) == parity) &
          .and. (get_ang_mom_proj(p) == m)) then

        sa_n_temp = sa_n_temp + 1

        sa_i_temp(sa_n_temp) = ii

      end if

    end do

    !< transfer temporary values across
    sa_n = sa_n_temp

    allocate(sa_i(1:sa_n))

    sa_i(:) = sa_i_temp(:)

  end subroutine symmetry_adaption

end module structure
