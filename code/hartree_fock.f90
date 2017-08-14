module hartree_fock
!> possible problem with odd number of electrons - check expression for G

  use input_data
  use grid_radial
  use sturmian_class
  use core_wavefunctions
  use omp_lib
  implicit none

  private
  public hf_procedure

contains

!> Performs Hartree-Fock procedure assuming that the basis has already been
!>  diagonalised, and that the hamiltonian contains the nuclear potential.
!> Records results of procedure in core_state.
  subroutine hf_procedure (basis, H, integrals, core)
    type(basis_sturmian_nr) , intent(in)  :: basis
    real*8                  , intent(in)  :: H(:, :)
    real*8                  , intent(in)  :: integrals(:, :, :, :)
    type(core_state)        , intent(out) :: core
    real*8                  , allocatable :: C(:, :)
    real*8                  , allocatable :: P(:, :)
    real*8                  , allocatable :: G(:, :)
    real*8                  , allocatable :: F(:, :)
    real*8                  , allocatable :: w(:)
    real*8                                :: hf_energy
    logical                               :: converged
    integer                               :: n, ii

    write (*, "(a)") "> hartree-fock procedure"
    write (*, "(a, a)") &
        " target: ", trim(data_in%target)
    write (*, *)

    n = basis_size(basis)

    allocate(C(1:n, 1:n))
    allocate(P(1:n, 1:n))
    allocate(G(1:n, 1:n))
    allocate(F(1:n, 1:n))
    allocate(w(1:n))

    C(:, :) = 0.0
    P(:, :) = 0.0
    G(:, :) = 0.0
    F(:, :) = 0.0
    w(:)    = 0.0

    converged = .true.
    ii        = 0

    do while ((ii <= data_in%n_e) .and. (converged))

      call hf_iteration(integrals, ii, H, C, P, G, F, w, converged, basis)

      call write_energy(H, integrals, ii, C, w)

      ii = ii + 2

    end do

    hf_energy = (data_in%Z1 * data_in%Z2 / data_in%Rd) + &
        calc_electronic_energy (data_in%n_e, C, P, H, F)

    call hf_core(basis, C, hf_energy, core)

  end subroutine hf_procedure


!> calculation schemes

!> hf_iteration
!>  H           core hamiltonian (kinetic + nuclear potential) matrix
!>  integrals   4-d matrix of two electron integrals
!>  C           basis expansion coefficients
!>  S           overlap matrix
!>  P           charge density matrix
!>  P_iter      iter-th iteration of charge density matrix
!>  G           two electron interation matrix
!>  F           fock matrix
!>  w           fock eigen-values
!>  converged   flags if hf_iteration procedure has converged
!>  K           dimension of basis
!>  iter        iteration count
  subroutine hf_iteration (integrals, n_e, H, C, P, G, F, w, converged, basis)
    type(basis_sturmian_nr) , intent(in)    :: basis
    real*8  , intent(in)    :: integrals(:, :, :, :)
    integer , intent(in)    :: n_e
    real*8  , intent(in)    :: H(:, :)
    real*8  , intent(inout) :: C(:, :)
    real*8  , intent(inout) :: P(:, :)
    real*8  , intent(inout) :: G(:, :)
    real*8  , intent(inout) :: F(:, :)
    real*8  , intent(inout) :: w(:)
    logical , intent(inout) :: converged
    real*8  , allocatable   :: P_iter(:, :, :)
    integer                 :: K
    integer                 :: iter
    integer                 :: ii, jj ! for printing matrix elements

    write (*, "(a, i3)") ">> hartree-fock iterative scheme"
    write (*, "(a, i4)") &
        " electrons: ", n_e

    K = size(C, 1)

    allocate(P_iter(1:data_in%iter_max, 1:K, 1:K))

    P_iter(:, :, :) = 0.0

    call set_density(n_e, C, P_iter(1, :, :))

    converged = .false.
    iter      = 1

    do while ((.not. converged) .and. (iter < data_in%iter_max))

      iter = iter + 1

      call set_two_electron(integrals, P_iter(iter - 1, :, :), C, n_e, G)

      call set_fock(H, G, F)

      call diagonalise(F, C, w)

      call set_density(n_e, C, P_iter(iter, :, :))

      call check_convergence(P_iter, iter, converged)

    end do

    if (converged) then

      P(:, :) = P_iter(iter, :, :)

      write (*, "(a, i5)") &
          " iterations: ", iter

      ! call print_matrix(C)
      ! call print_matrix(P)
      ! call print_matrix(G)
      ! call print_matrix(F)

    else

      write (*, *) ' failed to converge'

      write (*, '(a)') '> stopping'

      stop

    end if

    write (*, *)

  end subroutine hf_iteration


!> calculation procedures

  subroutine set_density (n_e, C, P)
    integer , intent(in)  :: n_e
    real*8  , intent(in)  :: C(:, :)
    real*8  , intent(out) :: P(:, :)
    integer               :: K, ii, jj, kk

    K = size(C, 1)

    !$omp parallel do private(ii, jj) shared (P, C, n_e, K)
    do ii = 1, K

      do jj = 1, K

        P(ii, jj) = 0.0

        do kk = 1, (n_e / 2)

          P(ii, jj) = P(ii, jj) + &
              (2.0 * C(ii, kk) * C(jj, kk))

        end do

      end do

    end do
    !$omp end parallel do

  end subroutine set_density


  subroutine set_two_electron (integrals, P, C, n_e, G)
    real*8  , intent(in)  :: integrals(:, :, :, :)
    real*8  , intent(in)  :: P(:, :)
    real*8  , intent(in)  :: C(:, :)
    integer , intent(in)  :: n_e
    real*8  , intent(out) :: G(:, :)
    integer               :: K
    integer               :: s
    integer               :: ii, jj, kk, ll

    K = size(P, 1)

    !$omp parallel do private(ii, jj, kk, ll) shared(G, P, integrals, K)
    do ii = 1, K

      do jj = 1, K

        G(ii, jj) = 0.0

        do kk = 1, K

          do ll = 1, K

            G(ii, jj) = G(ii, jj) + &
                (P(ll, kk) * (integrals(ii, kk, jj, ll) - &
                (0.5 * integrals(ii, kk, ll, jj))))

          end do

        end do

      end do

    end do
    !$omp end parallel do

    !< additional sum for odd-number of electrons
    if (mod(n_e, 2) == 1) then

      s = n_e / 2

      !$omp parallel do private(ii, jj, kk, ll) shared(G, C, integrals, K, S)
      do ii = 1, K

        do jj = 1, K

          do kk = 1, K

            do ll = 1, K

              G(ii, jj) = G(ii, jj) + &
                  (C(ll, s+1) * C(kk, s+1) * (integrals(ii, kk, jj, ll) - &
                  integrals(ii, kk, ll, jj)))

            end do

          end do

        end do

      end do
      !$omp end parallel do

    end if

  end subroutine set_two_electron


  subroutine set_fock (H, G, F)
    real*8 , intent(in)  :: H(:, :)
    real*8 , intent(in)  :: G(:, :)
    real*8 , intent(out) :: F(:, :)

    F(:, :) = H(:, :) + G(:, :)

  end subroutine set_fock


  subroutine diagonalise (F, C, w)
    real*8  , intent(in)  :: F(:, :)
    real*8  , intent(out) :: C(:, :)
    real*8  , intent(out) :: w(:)
    real*8  , allocatable :: S(:, :)
    integer               :: K
    integer               :: error
    integer               :: ii

    K = size(F, 1)

    allocate(S(1:K, 1:K))

    S(:, :) = 0.0
    !$omp parallel do private(ii) shared(S, K)
    do ii = 1, K

      S(ii, ii) = 1.0

    end do
    !$omp end parallel do

    call rsg(K, K, F(:, :), S(:, :), w(:), 1, C(:, :), error)

    if (error /= 0) then

      write (*, "(a)") "RSG failed to diagonalise the system"

    end if

  end subroutine diagonalise


  subroutine check_convergence (P, iter, converged)
    real*8  , intent(in)  :: P(:, :, :)
    integer , intent(in)  :: iter
    logical , intent(out) :: converged

    converged = (check_std_dev(P, iter) < data_in%tolerance)

  end subroutine check_convergence


  function check_std_dev (P, iter) result (std_dev)
    real*8  , intent(in)  :: P(:, :, :)
    integer , intent(in)  :: iter
    real*8                :: std_dev

    std_dev = sqrt(sum((P(iter, :, :) - P(iter - 1, :, :)) ** 2)) / size(P, 2)
    write (*, '(a, es10.2)') &
        ' std. dev.: ', std_dev

  end function check_std_dev


!> post-calculation procedures

  function calc_electronic_energy (n_e, C, P, H, F) result (energy)
    integer , intent(in) :: n_e
    real*8  , intent(in) :: C(:, :)
    real*8  , intent(in) :: P(:, :)
    real*8  , intent(in) :: H(:, :)
    real*8  , intent(in) :: F(:, :)
    real*8               :: energy
    integer              :: K
    integer              :: s
    integer              :: ii, jj, kk

    K = size(P, 1)

    energy = 0.0

    do ii = 1, K

      do jj = 1, K

        energy = energy + &
            (0.5 * P(jj, ii) * (H(ii, jj) + F(ii, jj)))

      end do

    end do

    !< odd number of electrons
    if (mod(n_e, 2) == 1) then

      s = n_e / 2

      do ii = 1, K

        do jj = 1, K

          energy = energy + &
              (C(ii, s+1) * C(jj, s+1) * ((1.5 * F(ii, jj)) - (0.5 * H(ii, jj))))

        end do

      end do

    end if

  end function calc_electronic_energy


  subroutine write_energy (H, integrals, n_e, C, w)
    real*8  , intent(in)  :: H(:, :)
    real*8  , intent(in)  :: integrals(:, :, :, :)
    integer , intent(in)  :: n_e
    real*8  , intent(in)  :: C(:, :)
    real*8  , intent(in)  :: w(:)
    real*8  , allocatable :: P(:, :)
    real*8  , allocatable :: G(:, :)
    real*8  , allocatable :: F(:, :)
    real*8                :: energy_nuclear
    real*8                :: energy_electronic
    real*8                :: energy_hf
    integer               :: K, ii

    K = size(C, 1)

    allocate(P(1:K, 1:K))
    allocate(G(1:K, 1:K))
    allocate(F(1:K, 1:K))

    call set_density(n_e, C, P)
    call set_two_electron(integrals, P, C, n_e, G)
    call set_fock(H, G, F)

    !< print energies
    write (*, '(a)') '>> energies'

    !< nuclear
    if ((data_in%Z1 < 1.0e-5) .or. (data_in%Z2 < 1.0e-5) &
        .or. (data_in%Rd < 1.0e-5)) then

      energy_nuclear = 0.0

    else

      energy_nuclear = data_in%Z1 * data_in%Z2 / data_in%Rd

    end if

    call write_title('nuclear')
    call write_value(energy_nuclear)
    write (*, *)

    !< orbitals
    call write_title('orbitals')

    do ii = 1, (n_e / 2)
      call write_value(w(ii))
      call write_value(w(ii))
    end do

    if (mod(n_e, 2) == 1) then
      call write_value(w((n_e / 2) + 1))
    end if

    write (*, *)

    !< electronic
    energy_electronic = calc_electronic_energy(n_e, C, P, H, F)
    call write_title('electronic')
    call write_value(energy_electronic)
    write (*, *)

    !< hf
    energy_hf = energy_nuclear + energy_electronic
    call write_title('hartree-fock')
    call write_value(energy_hf)
    write (*, *)

  end subroutine write_energy


  subroutine write_title (title)
    character(len = *) , intent(in) :: title

    write (*, '(a14, a)') &
        title, '           a.u.             eV             Ry'

  end subroutine write_title


  subroutine write_value (energy)
    real*8 , intent(in) :: energy

    write (*, '(a14, 3f15.5)') &
        ' ', energy, (27.2116 * energy), (2 * energy)

  end subroutine write_value


  subroutine print_matrix (matrix)
    real*8  , intent(in) :: matrix(:, :)
    integer              :: ii, jj, max_size

    max_size = 25

    if (size(matrix, 1) > max_size) then

      write (*, "(a)") "(truncated)"

    end if

    do ii = 1, min(max_size, size(matrix, 1))

      do jj = 1, min(max_size, size(matrix, 2))

        if (abs(matrix(ii, jj)) < 1.0e-5) then

          write (*, "(a6)", advance = 'no') "."

        else

          write (*, "(f6.2)", advance = 'no') matrix(ii, jj)

        end if

      end do

      write (*, *)

    end do

    write (*, *)

  end subroutine print_matrix


!> writes hartree_fock results to file
!>  temp(n, l, :) is the plot of the n-th spatial orbital's l-th partial wave.
  subroutine hf_write_results (basis, C, hf_energy, filepath)
    type(basis_sturmian_nr) , intent(in)  :: basis
    real*8                  , intent(in)  :: C(:, :)
    real*8                  , intent(in)  :: hf_energy
    character(len = *)      , intent(in)  :: filepath
    real*8                  , allocatable :: temp(:, :, :)
    integer                               :: s, l, ii
    integer                               :: unitno
    real*8                                :: spectroscopic

    ! allocate(temp(1:grid%nr, -data_in%latop:data_in%latop, 0:data_in%latop, &
    !     1:((data_in%n_e + 1)/2)))
    allocate(temp(1:grid%nr, 0:data_in%latop, 1:((data_in%n_e + 1)/2)))

    !< plot partial wave expansions of spatial orbitals
    temp(:, :, :) = 0.0

    do s = 1, (data_in%n_e + 1) / 2

      do l = 0, data_in%latop

        do ii = 1, basis_size(basis)

          if (get_ang_mom(basis%b(ii)) == l) then

            temp(:, l, s) = temp(:, l, s) + &
                (C(ii, s) * fpointer(basis%b(ii)))

          end if

        end do

      end do

    end do

    !< remove underflow errors
    where (abs(temp) > 1.0e5)

      temp = 0.0

    end where

    !< write plots to file
    unitno = 1000
    open (unitno, file = filepath)

    write(unitno, *) 0 ! core%m
    write(unitno, *) data_in%n_e
    write(unitno, *) data_in%labot, data_in%latop
    write(unitno, *) data_in%Rd
    write(unitno, *) hf_energy
    write(unitno, *) data_in%Z1 * data_in%Z2 / data_in%Rd

    write(unitno, *) grid%nr
    write(unitno, *) grid%gridr(:)

    do s = 1, (data_in%n_e + 1) / 2

      write (*, "(a, i3)") ">> spatial orbtial:", s
      write (*, "(a3, a8)") "l", "<l|l>"

      do l = data_in%labot, data_in%latop

          write (unitno, *) temp(:, l, s)

          spectroscopic = sum((temp(:, l, s) ** 2) * grid%weight(:))

          if (abs(spectroscopic) > 1.0e-5) then

            write (*, "(i3, f8.5)") &
                l, spectroscopic

          end if

      end do

      write (*, *)

    end do

    close (unitno)

  end subroutine hf_write_results


!> records hartree_fock results in a core_state variable.
!>  pw(:, l, s) is the plot of the s-th spatial orbital's l-th partial wave.
  subroutine hf_core (basis, C, hf_energy, core)
    type(basis_sturmian_nr) , intent(in)  :: basis
    real*8                  , intent(in)  :: C(:, :)
    real*8                  , intent(in)  :: hf_energy
    type(core_state)        , intent(out) :: core
    real*8                  , allocatable :: pw(:, :, :)
    integer                               :: s, l, ii
    real*8                                :: spectroscopic

    allocate(pw(1:grid%nr, 0:data_in%latop, 1:((data_in%n_e + 1)/2)))

    !< plot partial wave expansions of spatial orbitals
    pw(:, :, :) = 0.0

    do s = 1, (data_in%n_e + 1) / 2

      do l = 0, data_in%latop

        do ii = 1, basis_size(basis)

          if (get_ang_mom(basis%b(ii)) == l) then

            pw(:, l, s) = pw(:, l, s) + &
                (C(ii, s) * fpointer(basis%b(ii)))

          end if

        end do

      end do

    end do

    !< remove underflow errors
    where (abs(pw) > 1.0e5)

      pw = 0.0

    end where

    do s = 1, (data_in%n_e + 1) / 2

      write (*, "(a, i3)") ">> spatial orbtial:", s
      write (*, "(a3, a8)") "l", "<l|l>"

      do l = data_in%labot, data_in%latop

          spectroscopic = sum((pw(:, l, s) ** 2) * grid%weight(:))

          if (abs(spectroscopic) > 1.0e-5) then

            write (*, "(i3, f8.5)") &
                l, spectroscopic

          end if

      end do

      write (*, *)

    end do

    !< record state information and partial waves in core_state
    call core%construct_core_state(0, data_in%n_e, data_in%labot, data_in%latop,&
        data_in%Rd, hf_energy, data_in%Z1 * data_in%Z2 / data_in%Rd, pw)

  end subroutine hf_core

end module hartree_fock
