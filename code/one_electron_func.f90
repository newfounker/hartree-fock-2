subroutine construct_1el_basis_nr(Number_one_electron_func, output_dir)
    use input_data
    use sturmian_class
    use vnc_module
    use one_electron_func_mod
    use target_states
    use ovlpste1me
    use hartree_fock
    use core_wavefunctions

    implicit none

    integer, intent(in):: Number_one_electron_func
    !
    real*8, dimension(:,:), allocatable:: H, b, CI
    real*8, dimension(:), allocatable:: w
    integer:: i, j, n, nd, N_oneel, jstart, jstop, ipar, l_n, ipar_n, li, lj, Nmax, Ntmp, lorb, itmp
    real*8:: al
    integer:: matz, ierr
    integer:: ma, ni, nj
    real*8:: tmp, res, tmp1, rval,  tmp_sign_CI
    real*8:: energy
    integer, dimension(:), allocatable:: no, mo
    logical:: hlike
    type(sturmian_nr), pointer:: pi, pj
    real*8::  VLambdaR_ME

    !< added by tom ross
    character(len = *) , intent(in) :: output_dir
    logical :: prediag
    real*8 :: hf_energy
    type(core_state) :: core

    !
    !
    hlike = .true.

    ! Create space for basis of  1-electron pseudostates
    call new_basis_st(TargetStates,Number_one_electron_func)
!??    write(*,'("allocated space for one-electron target states Nmax=",I5)') Number_one_electron_func


    N_oneel = 0

    ! For given input file data.in construct Sturmian basis for all values of l: is kept in module one_electron_func
!??    print*, 'Start making nonrrel  Sturmian basis'

    prediag = .true.

    if (prediag) then

        call construct_prediag(bst, data_in)

    else

        call construct(bst, data_in)

    end if


    nspm = basis_size(bst) ! number of sturmian finctions
!??    print*, 'Size of Sturmian basis: nspm=', nspm
    allocate(no(nspm),mo(nspm))
    no (:) = 0
    mo (:) = 0
    !

    do ma = data_in%Mt_min, data_in%Mt_max
      do ipar = 1, -1, -2

        if(data_in%nst(ma, ipar) .le. 0) cycle

!!$   form set of (one-electron) configurations to be used to digonalise H
        i = 0
        do n=1,nspm
          l_n = get_ang_mom(bst%b(n))
          ipar_n = (-1)**(l_n)
            if( l_n .ge. ma .and. ipar .eq. ipar_n) then
                i = i + 1
                no(i) = n
!??                write (*, *) n, l_n, get_energy(bst%b(n))
            endif
        end do

        nd = i    ! number of sturmian orbitals to be used in diagonlization
        if(nd .eq. 0) then
          write (*, '(a, i3)') &
              'zero symmetry-adapted basis functions for m, par: ', ma, ipar
            stop
        endif

        if(data_in%nst(ma, ipar) .gt. nd) then
          print*, 'one_electron_func.f:  data_in%nst(ma, ipar) .gt. nd :', &
              data_in%nst(ma, ipar), nd
          stop
        endif


        mo(1:nd) = ma

        ! Temporary arrays
        allocate(H(nd,nd))
        allocate(b(nd,nd))
        allocate(CI(nd,nd))
        allocate(w(nd))

        H(:,:) = 0d0
        b(:,:) = 0d0

!!$ Calculate H matrix: see DF notes for details
!!$ Here we rely on the special form of Laguerre functions of order (l+1)
!!$ in order to  calculate overlaps and 1/r integrals analytically
        b(1:nd,1:nd) = bst%ortint(no(1:nd),no(1:nd))

!> prediag: true -> kinetic energy already diagonalised
        if (prediag) then

            do ni = 1, nd

                H(ni, ni) = get_energy(bst%b(no(ni)))

            end do

        else

!!$   get kinetic energy: (using special properties of Lag func)
            do ni =1,nd
                i = no(ni)
                li = get_ang_mom(bst%b(i))
                al = data_in%alpha(li)
                do nj=1,nd
                    j = no(nj)
                    lj = get_ang_mom(bst%b(j))
                    if(li.eq. lj) then
                        H(ni,nj) = -al*al*b(ni,nj)/2d0
                    endif
                enddo
            enddo

!!$ Laguerre functions of order (l+1) are diagonal with 1/r potential weight.
            do ni = 1,nd
                i = no(ni)
                li = get_ang_mom(bst%b(i))
                al = data_in%alpha(li)
                H(ni,ni) = H(ni,ni) + al*al
            end do

        end if

        do ni = 1,nd
            i = no(ni)
            pi => bst%b(i)

            do nj = 1,ni
                j = no(nj)
                pj => bst%b(j)

                tmp =  VLambdaR_ME(pi,pj,ma)

                ! tmp1 = data_in%Z1*data_in%Z2/data_in%Rd * bst%ortint(i,j)
                tmp1 = 0.0

                H(ni,nj) = H(ni,nj) + tmp + tmp1
                H(nj,ni) = H(ni,nj)
                !                call Hlagorb(bst,i,j,ma,result)
                !             print'(4i5,2E15.5)', i, j, ni, nj, H(ni,nj), result

            end do
          end do

        write (*, '(a)') 'symmetry'
        write (*, '(a, i5)') ' m:      ', ma
        write (*, '(a, i5)') ' parity: ', ipar
        write (*, '(a, i5)') ' size:   ', nd


!??          matz=2
!??          call rsg(nd, nd, H, b, w, 2, CI, ierr)
!??          write(*,'("ierr =",I3)') ierr
!??          print*, ' Energies in a.u.'
!??          write(*,'(5F15.5)') (real(w(i)), i=1,nd)
!??          print*
!??          print*, ' Energies in eV'
!??          write(*,'(5F15.5)') (real(27.2116*w(i)), i=1,nd)

!>> hartree fock procedure
        write (*, *)

        call rsg(nd, nd, H, b, w, 2, CI, ierr)

        ! call hf_procedure(bst, nd, no, ma, ipar, H, CI, w, output_dir//"/hf_results.dat")
        ! call core%read_from(output_dir//"/hf_results.dat")
        ! call core%write_pw_to(output_dir//"/core_plots.dat")
        ! call core%write_coulomb_pw_to(output_dir//"/core_coulomb.dat")
        ! call core%write_potential_pw_to(output_dir//"/core_potential.dat")
        ! call core_spectrum(core, bst, nd, no, H, b)

!!$ Create basis of 1-electron pseudostates from Laguere basis and CI coef.

        jstart = 1
        jstop = min(nd, jstart+data_in%nst(ma, ipar) - 1)
        do j=jstart,jstop
            energy = w(j)
            N_oneel = N_oneel + 1
            tmp_sign_CI = SUM(CI(1:nd,j))
            if(  tmp_sign_CI .lt. 0d0 ) then
                CI(1:nd,j) = -CI(1:nd,j)
            endif
            TargetStates%Nstates = N_oneel
            call  construct_st(TargetStates%b(N_oneel),dble(ma),ipar,0.5d0,energy,j,nd,CI(1:nd,j),no(1:nd),mo(1:nd))
!!$ NOTE here we might need to acount for degeneracy of the molecule energy levels for ma .ne. 0
!!$ this could lead to introdicing addtonal state with -ma
            if( ma .ne. 0 ) then
                N_oneel = N_oneel + 1
                TargetStates%Nstates = N_oneel
                call  construct_st(TargetStates%b(N_oneel),-dble(ma),ipar,05d0,energy,j,nd,CI(1:nd,j),no(1:nd),-1*mo(1:nd))
            endif


        end do

        deallocate(H)
        deallocate(b)
        deallocate(CI)
        deallocate(w)
      enddo
    enddo

    if( TargetStates%Nstates .ne. TargetStates%Nmax) then
       print*, 'one_electron_func.f90: TargetStates%Nstates .ne. TargetStates%Nmax:', TargetStates%Nstates, TargetStates%Nmax
       stop
    endif


    print*
    call sort_by_energy_basis_st(TargetStates)
    ! call calc_spectro_factors(TargetStates)
!??    call print_energy_basis_st(TargetStates)

!!$  populate arary e1me to be used later in H12 diagonalization
    Nmax = basis_size_st(TargetStates)
    if(allocated(e1me)) deallocate(e1me)
    allocate(e1me(Nmax, Nmax))
    allocate(ovlpst(Nmax, Nmax))
    e1me(:,:) = 0d0
    ovlpst(:,:) = 0d0
    do j=1,Nmax
       e1me(j,j) = get_energy_st(TargetStates%b(j))
       ovlpst(j,j) = 1d0
    enddo


  end subroutine construct_1el_basis_nr


  !$$  It is assumed that basis bst is made from Laguerre function of
!$$  order  (2l+1)  that are orthogonal with weight 1/r
!$$  these functions have definite value of orbital angular momentum l
 subroutine Hlagorb(bst,ist,jst,ma,result)

    use sturmian_class

    implicit none

    type(basis_sturmian_nr), intent(in) :: bst   ! this is Sturmian basis
    integer, intent(in):: ist,jst
    integer, intent(in):: ma
    real*8, intent(out):: result

    type(sturmian_nr), pointer:: pi, pj
    integer:: li, lj
    real*8:: al, tmpovlp
    real*8::  VLambdaR_ME

    result = 0d0

    pi => bst%b(ist)
    pj => bst%b(jst)

    li = get_ang_mom(pi)
    lj = get_ang_mom(pj)

    if(abs(ma) .gt. li  .or.  abs(ma) .gt. lj)  return

    if(lj .eq. li) then

    al = get_alpha(bst%b(ist))    ! it is assumed that it is Laguerre functions, and made for the same alpha forgiven l
                                  ! term with  alpha  appears only in matrix elements of the kinetic energy operator where li=lj

       tmpovlp = bst%ortint(ist,jst)

       result = -al*al*tmpovlp/2d0
       if(ist .eq. jst) then
          result = result + al*al
       endif
    endif

!    print*, result, VLambdaR_ME(pi,pj,ma)
    result = result + VLambdaR_ME(pi,pj,ma)


  end subroutine Hlagorb
!
!
!!$  It is assumed that H2+ gs orbital has m=0 and it is the first in TargetStates array.


subroutine construct_prediag(bst_orth, dataARG)
  use input_data
  use grid_radial
  use sturmian_class
  use vnc_module
  use one_electron_func_mod
  use target_states
  use ovlpste1me

  implicit none

  type(input)             , intent(in)  :: dataARG
  type(basis_sturmian_nr) , intent(out) :: bst_orth
  type(basis_sturmian_nr)               :: bst_temp
  real*8                  , allocatable :: potential(:)
  integer                               :: ll, kk, n

  allocate(potential(1:grid%nr))

  potential = 0.0

  call new_basis_nr(bst_orth, sum(dataARG%nps(dataARG%labot:dataARG%latop)))

  n = 1

  do ll = dataARG%labot, dataARG%latop

      call construct_wflocalpot_nr(bst_temp, dataARG%nps(ll), potential, dataARG%nps(ll), ll, dataARG%alpha(ll))

      do kk = 1, dataARG%nps(ll)

          call copy(bst_orth%b(n), bst_temp%b(kk))

          bst_orth%ortint(n, n) = 1d0

          n = n + 1

      end do

      call destruct(bst_temp)

  end do

end subroutine construct_prediag
