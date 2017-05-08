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

contains

!> Performs the Hartree-Fock procedure for the system specified in the global
!> variable data_in.
!> Writes the results of the HF procedure to the specified output file
  subroutine hf_structure (m, hf_file)
    integer                 , intent(in)  :: m
    character(len = *)      , intent(in)  :: hf_file
    type(basis_sturmian_nr)               :: basis
    real*8                  , allocatable :: H(:, :)
    integer                               :: n, ii
    type(core_state)                      :: core

    !< construct basis
    call construct_diagonalised(basis)

    n = basis_size(basis)

    !< construct nuclear hamiltonian matrix
    allocate(H(1:n, 1:n))
    call nuclear_hamiltonian(basis, H)

    !< perform hartree_fock procedure (assuming m = 0)
    call hf_procedure(basis, H, m, hf_file)

    !< core wavefunctions
    ! call core%read_from("../output/hf_results.dat")
    ! call core%write_pw_to("../output/core_plots.dat")
    ! call core%write_coulomb_pw_to("../output/core_coulomb.dat")
    ! call core%write_potential_pw_to("../output/core_potential.dat")
    ! call core_spectrum(core, basis, n, no, H, S)

  end subroutine hf_structure

!> Constructs a diagonalised basis, specified by the global data_in variable,
!> which is diagonalised with regard to the kinetic energy operator.
!> After this basis is constructed, the kinetic energy matrix can be built by
!> T(ii, ii) = get_energy(bst%b(ii))
  subroutine construct_diagonalised (basis)
    type(basis_sturmian_nr) , intent(out) :: basis
    type(basis_sturmian_nr)               :: basis_l
    type(sturmian_nr)       , pointer     :: p
    real*8                  , allocatable :: potential(:)
    integer                               :: size, n, l, k, m

    write (*, "(a)") "> establishing basis"

    !< allocate entire basis (indexing over all k, l, m values)
    size = 0
    do l = data_in%labot, data_in%latop
      size = size + (data_in%nps(l) * ((2 * l) + 1))
    end do

    call new_basis_nr(basis, size)

    !< zero-potential, used in basis diagonalisation
    allocate(potential(1:grid%nr))
    potential = 0.0

    !< entire basis index
    n = 1

    !< loop through basis for given l value
    write (*, "(4a4)") "n", "k", "l", "m"
    do l = data_in%labot, data_in%latop

      !< construct diagonal basis (w.r.t kinetic energy) with m = 0, for given
      !< l.
      call construct_wflocalpot_nr(basis_l, data_in%nps(l), potential, &
          data_in%nps(l), l, data_in%alpha(l))

      do k = 1, data_in%nps(l)

        do m = -l, l, 1

          call copy(basis%b(n), basis_l%b(k))

          call set_ang_mom_proj(basis%b(n), m)

          basis%ortint(n, n) = 1.0

          p => basis%b(n)
          write (*, "(4i4)") n, get_k(p), get_ang_mom(p), get_ang_mom_proj(p)

          n = n + 1

        end do

      end do

      call destruct(basis_l)

    end do

    write (*, *)

  end subroutine construct_diagonalised

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

end module structure
