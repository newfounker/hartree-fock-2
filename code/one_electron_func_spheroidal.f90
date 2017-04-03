!!$  This subroutine performs one-electron diagonalization for diatomic molecule like H2+ 
!!$  in spheroidal coordinates and a single centre scheme
!!$  for given number of conserved quantum numbers: magnetic sublevel and parity


subroutine construct_1el_basis_spheroidal(numStates)
  !
  ! JS
  !
  use input_data
  use grid_radial
  use one_electron_func
  use sturmian_class
  use state_class
  use  target_states

  implicit none

  integer, intent(in) :: numStates

  integer :: basisType, stateNum, stateCount, func, basisSize, numFunctions, countFunc, i,j, ni,nj, li,lj, m,par, ierr
  integer, dimension(:), allocatable :: nArray, mArray
  real*8 :: R, lambda, factor, numerical
  real*8, dimension(:), allocatable :: energies
  type(sturmian_nr), pointer :: fi, fj
  real*8, dimension(:,:), allocatable :: overlapMatrix, HMatrix, CIMatrix
  logical :: hlike, debug=.FALSE.


  basisType = data_in%calculation_type
  hlike = data_in%hlike
  R = data_in%Rd

  ! Create space for the target states.
  call new_basis_st( TargetStates, numStates, hlike, basisType )
  stateCount = 0
  write(*,'(A,I4,A)') 'Created space for', numStates, ' one-electron target states.'
  print*

  ! Construct the radial basis functions.
  call construct(bst,data_in)
  basisSize = basis_size(bst)
  numFunctions = ext_basis_size(bst)
  write(*,'(A,I4,A,I4,A)') 'Created a spheroidal basis of size', basisSize, ' extended to', numFunctions, ' functions.'
  print*
  print*

  allocate( nArray(basisSize), mArray(basisSize) )


  ! Loop through each (m,par) symmetry and do stuff.
  do m = data_in%Mt_min, data_in%Mt_max   ! Loop through each m as specified by the input file.
     do par = 1,-1,-2   ! Loop through positive and negative parity.
        write(*,'(3(A,I2))') 'Symmetry (m,par) = (', m, ',', par, '). nst =', data_in%nst(m,par)

        if( data_in%nst(m,par) .gt. 0 ) then   ! States will be constructed for this symmetry.

           !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
           ! This block of code determines which of the Sturmian functions will be used in
           ! constructing the particular (m,par) state.

           nArray = -1   ! The indices of the extended basis functions that will construct a state.
           mArray = -1   ! The corresponding values of m.
           
           countFunc = 0   ! He sounds cool. I'd like to meet him.
           do func = 1, numFunctions   ! Loop through the functions.

              if( get_ang_mom_proj( bst%b(bst%iArray(func)) ).eq.m .and. (-1)**bst%lArray(func).eq.par ) then   ! Function has the correct values of m and parity.
                 countFunc = countFunc + 1
                 
                 nArray(countFunc) = func
                 mArray(countFunc) = m
              endif

           enddo ! func

if(debug) then
   print*, 'nArray =', nArray(1:countFunc)
   print*, 'mArray =', mArray(1:countFunc)
   print*
endif

           !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
           ! This block of code builds the Hamiltonian matrix with the chosen functions.

           allocate( overlapMatrix(countFunc,countFunc), HMatrix(countFunc,countFunc), CIMatrix(countFunc,countFunc), energies(countFunc) )
           HMatrix = 0d0

           overlapMatrix = bst%ortint( nArray(1:countFunc), nArray(1:countFunc) )

           lambda = 2.0 * data_in%alpha(m)
           do i = 1, countFunc   ! Loop through the relevant functions.
              fi => bst%b( bst%iArray(nArray(i)) )
              ni = get_k(fi)
              li = bst%lArray( nArray(i) )

              factor = R / (2.0*li+1.0)
              do j = li-m+1, li+m   ! (l+m)! / (l-m)!
                 factor = factor * j
              enddo

              do j = 1, i   ! Fill out the lower triangle, diagonal inclusive.
                 fj => bst%b( bst%iArray(nArray(j)) )
                 nj = get_k(fj)
                 lj = bst%lArray( nArray(j) )

                 if( ni.eq.nj+2 .and. li.eq.lj ) then
                    HMatrix(i,j) = factor * -0.125 * sqrt( (ni-2.0)*(ni-1.0) * (ni+m-2.0)*(ni+m-1.0) )
                 elseif( ni.eq.nj+1 .and. li.eq.lj ) then
                    HMatrix(i,j) = factor * 0.25 * sqrt( (ni-1.0)*(ni+m-1.0) ) * ( lambda + 4.0*R/lambda )
                 elseif( ni.eq.nj .and. li.eq.lj ) then
                     HMatrix(i,i) = factor * 0.25 * ( ni**2 - ni + ni*m - m/2.0 + 1.0 + lambda*(2.0*ni+m-1.0) - 4.0*R*(2.0*ni+m-1.0+lambda)/lambda + 2.0*li*(li+1.0) )
                 endif

                 ! Numerically integrate the extra m-term in the Hamiltonian.
                 call H_mterm_num( fi, fj, numerical )
                 HMatrix(i,j) = HMatrix(i,j) - factor * m**2 * acos(0d0) * numerical

                 HMatrix(j,i) = Hmatrix(i,j)   ! Take advantage of symmetry.

              enddo ! j

           enddo ! i

if(debug) then
   print*, 'Overlap matrix elements'
   do i = 1, countFunc
      write(*,'(100F10.5)') ( overlapMatrix(i,j), j = 1,countFunc )
   enddo
   print*
   print*, 'Hamiltonian matrix elements'
   do i = 1, countFunc
      write(*,'(100F10.5)') ( HMatrix(i,j), j = 1,countFunc )
   enddo
   print*
endif

           !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
           ! This block of code solves the generalised eigenvalue problem Hv=abv
           ! and builds the target states.

           call rsg( countFunc, countFunc, HMatrix, overlapMatrix, energies, 1, CIMatrix, ierr)

if(debug) then
   print*, 'CI coefficients'
   do i = 1, countFunc
      write(*,'(100F10.5)') ( CIMatrix(i,j), j = 1,countFunc )
   enddo
   print*
endif

           print*
           print*, 'Energies in a.u. :'
           write(*,'(5F20.5)') energies
           print*, 'Energies in eV   :'
           write(*,'(5F20.5)') energies * 27.211
           print*
           print*

           ! Build the target states requested in data.in.
           do stateNum = 1, data_in%nst(m,par)
              stateCount = stateCount + 1

              if( sum( CIMatrix(:,stateNum) ) .lt. 0d0 ) CIMatrix(:,stateNum) = -CIMatrix(:,stateNum)   ! Make the state predominately positive, I think ?

              call construct_st(TargetStates%b(stateCount), hlike, dble(m), par, 0.5d0, energies(stateNum), stateNum, countFunc, CIMatrix(:,stateNum), nArray(1:stateCount), mArray(1:stateCount) )

              if( m.gt.0 ) then   ! We have degeneracy.
                 stateCount = stateCount + 1
                 call construct_st(TargetStates%b(stateCount), hlike, dble(-m), par, 0.5d0, energies(stateNum), stateNum, countFunc, CIMatrix(:,stateNum), nArray(1:stateCount), mArray(1:stateCount) )
              endif
           enddo


           deallocate( overlapMatrix, HMatrix, CIMatrix, energies )


        endif
     enddo ! par
  enddo ! m


  print*
  print*
  write(*,'(2(A,I4))') 'Number of states created:', stateCount, ' and requested:', numStates
  print*
  print*

  if( stateCount .ne. numStates ) stop

  call sort_by_energy_basis_st(TargetStates)
  call calc_spectro_factors(TargetStates)
  call print_energy_basis_st(TargetStates)


end subroutine construct_1el_basis_spheroidal



subroutine construct_1el_basis_spheroidal2(Number_one_electron_func)
    use input_data
    use grid_radial
    use sturmian_class
    use vnc_module
    use one_electron_func
    use state_class
    use  target_states
  
    implicit none

    integer, intent(in):: Number_one_electron_func
!
    real*8, dimension(:,:), allocatable:: H, b, CI
    real*8, dimension(:), allocatable:: w, HNum
    integer:: la, nd, i, j, k, m, i1, i2, N_oneel, jstart,Nst, lamtmp, n, nsp, li,lj, mi,mj, lam_min, lam_max, lam, ipar, ma_n, lmin, lmax, ndtot, offset
    real*8:: al
    integer:: matz, ierr
    integer:: maxfi, maxfj, minfi, minfj, minf, maxf, imax,imin
    integer:: ma, ni, nj, indi, indj
    real*8, pointer, dimension(:):: fi, fj, weight, gridr
    real*8, dimension(grid%nr):: temp
    integer:: minf_temp, maxf_temp, minf_bj, maxf_bj
    real*8:: tmp, res, tmp1, rval,  tmp_sign_CI
    real*8:: energy, Rd, lambda
    integer, dimension(:), allocatable:: no, mo
    logical:: hlike
    real*8:: Yint
    type(sturmian_nr), pointer:: pi, pj
    real*8::  VLambdaR_ME
!
!
    basis_type = data_in%calculation_type
    if(basis_type .ne. 2) then
       print*, 'Only calculations with a spheroidal basis should be here!'
       stop
    endif
    hlike = .true.

! Create space for basis of  1-electron pseudostates
    call new_basis_st(TargetStates,Number_one_electron_func,hlike,basis_type)
    write(*,'("allocated space for one-electron target states Nmax=",I5)') Number_one_electron_func
    

    N_oneel = 0
    weight => grid%weight
    
    
! For all given (la,nd,al) construct Sturmian basis: is kept in module one_electron_func
    print*, 'Start making nonrrel  Sturmian basis'
    call construct(bst,data_in)
    Nst = basis_size(bst) ! number of sturmian finctions
!    nspm =  basis_size(bst)
    print*, 'Size of Sturmian basis: Nst=', Nst 
    allocate(no(Nst),mo(Nst)) 
    no (:) = 0   ! Array of n of each state.
    mo (:) = 0   ! Array of m of each state.

    do ma = data_in%Mt_min, data_in%Mt_max   ! Loop through m as specified in the input file.
       do ipar = 1,-1,-2   ! Positive and negative parities.
          if(data_in%nst(ma,ipar) .le. 0) cycle   ! Skip non-positive values of nst.

          ! Form set of (one-electron) configurations to be used to diagonalise H2+.
          i = 0
          do n=1,Nst
             ma_n = get_ang_mom_proj(bst%b(n))
!             ipar_n = (-1)**(ma_n)
!             if( ma_n.eq.ma .and. ipar.eq.ipar_n ) then   ! l>=m and correct parity.
             if ( ma_n.eq.ma ) then
                i = i + 1
                no(i) = n
             endif
          end do

          nd = i    ! Number of sturmian orbitals to be used in diagonalisation.
          lmin = max( data_in%labot, ma )   ! lmin >= labot & ma
          lmin = lmin + mod( lmin + (1-ipar)/2, 2 )   ! Correct parity.
          lmax = data_in%latop - mod( data_in%latop + (1-ipar)/2, 2 )   ! lmax = latop or latop-1
          ndtot = nd * (lmax-lmin+2)/2   ! Total length of the matrices is n*l.

          if(ndtot .le. 0) then
             write(*,'(3(A,I2),A)'), 'ndtot =', ndtot, ' (ma =', ma, ' and ipar =', ipar, ')'
             stop
          endif
          mo(1:nd) = ma

          ! Temporary arrays
          allocate(H(ndtot,ndtot))
allocate( Hnum(ndtot) )
          allocate(b(ndtot,ndtot))
          allocate(CI(ndtot,ndtot))
          allocate(w(ndtot))
          
          H(:,:) = 0d0
          b(:,:) = 0d0

          Rd = data_in%Rd

          ! Copy overlap matrix.
!          b(1:nd,1:nd) = bst%ortint( no(1:nd), no(1:nd) )
          ! Calculate overlap matrix elements.
          do la = lmin, lmax, 2
             do indi=1,nd
                i = no(indi)   ! Re-use the same nd functions for each l.
                ni = get_k( bst%b(i) )

                lambda = 2.0 * data_in%alpha(la)

                do indj=1,indi
                   j = no(indj)
                   nj = get_k( bst%b(j) )

                   tmp = Rd**3 / ( 4.0 * lambda**2 * (2.0*la+1.0) )
                   do m = la-ma+1, la+ma
                      tmp = tmp * dble(m)
                   enddo

                   offset = nd * (la-lmin)/2   ! Progresses down the matrix as l increases.

                   if (ni .eq. nj+2) then
                      b(indi+offset,indj+offset) = sqrt( nj*(nj+1.0) * (nj+ma)*(nj+ma+1.0) ) * tmp
                   elseif (ni .eq. nj+1) then
                      b(indi+offset,indj+offset) = sqrt( nj * (nj+ma+0.0) ) * tmp * -4.0*(nj+ma/2.0+lambda/2.0)
                   elseif (ni .eq. nj) then
                      b(indi+offset,indj+offset) = tmp * ( 6.0*(nj**2-nj+nj*ma-ma/2.0+ma**2/6.0+1.0/3.0) + 4.0*lambda*(nj+ma/2.0-0.5) + 2.0*lambda**2*(la**2+la+ma**2-1.0)/(2.0*la-1.0)/(2.0*la+3.0) )
                      if ( la+2 .le. lmax ) then
                         b(indi+offset+nd,indj+offset) = tmp * -lambda**2 * (la+ma+1.0) * (la+ma+2.0) / (2.0*la+3.0) / (2.0*la+5.0)
                         b(indi+offset,indj+offset+nd) = b(indi+offset+nd,indj+offset)
                      endif
                   endif
                   b(indj+offset,indi+offset) = b(indi+offset,indj+offset)
                enddo ! indj
             enddo ! indi
          enddo ! la

          ! Calculate H matrix.
          do la = lmin, lmax, 2
             do indi=1,nd
                i = no(indi)
                ni = get_k( bst%b(i) )
             
                lambda = 2.0 * data_in%alpha(la)

                do indj=1,indi
                   j = no(indj)
                   nj = get_k( bst%b(j) )

                   tmp = 0.125 * Rd / (2.0*la+1.0)
                   do m = la-ma+1, la+ma
                      tmp = tmp * dble(m)
                   enddo

                   offset = nd * (la-lmin)/2   ! Progresses down the matrix as l increases.

                   if( ni .eq. nj+2 ) then
                      H(indi+offset,indj+offset) = -tmp * sqrt( nj*(nj+1.0) * (nj+ma)*(nj+ma+1.0) )
                   elseif( ni .eq. nj+1 ) then
                      H(indi+offset,indj+offset) = 2.0 * tmp * sqrt( nj * (nj+ma+0.0) ) * ( lambda + 4.0*Rd/lambda )
                   elseif( ni .eq. nj ) then
                      H(indi+offset,indj+offset) = 2.0 * tmp * ( nj**2-nj+nj*ma-ma/2.0+1.0 + lambda*(2.0*nj+ma-1.0) + 2.0*la*(la+1.0) - 4.0*Rd*(2.0*nj+ma-1.0+lambda)/lambda )
                   endif
                   
                   call H_mterm_num( bst%b(i), bst%b(j), tmp )   ! Numerically integrates the nasty term.
                   tmp = -acos(0.0) * Rd / (2.0*la+1.0) * ma**2 * tmp
                   do m = la-ma+1, la+ma
                      tmp = tmp * dble(m)
                   enddo

                   H(indi+offset,indj+offset) = H(indi+offset,indj+offset) + tmp
                   
                   H(indj+offset,indi+offset) = H(indi+offset,indj+offset)

                enddo   ! indj
             enddo   ! indi
          enddo   ! la
!          H = H + b / Rd   ! 1/R neutron-neutron repulsion term.

          print*
          write(*, '("Symmetry and size: ma, ipar, ndtot = ",2I3,I5)') ma, ipar, ndtot
          print*, 'Overlap matrix elements:'
          do i=1,ndtot
             write(*,'(100F10.5)') (b(i,j), j=1,ndtot)
          enddo
          print*, 'Hamiltonian matrix elements:'
          do i=1,ndtot
             write(*,'(100F10.5)') (H(i,j), j=1,ndtot)
          enddo
          print*
!!$   
!!$
!!$          do la=lmin,lmax,2
!!$             do i=1,nd
!!$                do matz=lmin,lmax,2
!!$                   do j=1,nd
!!$                      call Hnum_element_spheroidal( bst%b(no(i)), bst%b(no(j)), la, matz, tmp )
!!$                      Hnum( j + nd*(matz-lmin)/2 ) = tmp
!!$                   enddo
!!$                enddo
!!$                write(*,'(100F10.5)') HNum
!!$             enddo
!!$          enddo


          matz=2
          call  rsg(ndtot,ndtot,H,b,w,matz,CI,ierr)
          write(*,'("ierr =",I3)') ierr
          print*, ' Energies in a.u.'
          write(*,'(5F15.5)') (real(w(i)), i=1,nd)
          print*
          print*, ' Energies in eV'
          write(*,'(5F15.5)') (real(27.2116*w(i)), i=1,nd)

print*
do i=1,ndtot
   write(*,'(100F10.5)') (CI(i,j), j=1,ndtot)
enddo
print*


!!$temp = 0.0
!!$do i=1,nd
!!$   ni = no(i)
!!$   fi => fpointer( bst%b(ni) )
!!$   minf_temp = get_minf( bst%b(ni) )
!!$   maxf_temp = get_maxf( bst%b(ni) )
!!$   do li = data_in%labot, data_in%latop, 2
!!$      print*, ni, li, CI(i+li*nd/2,1)
!!$      temp(minf_temp:maxf_temp) = temp(minf_temp:maxf_temp) + CI(i+li*nd/2,1) * fi(minf_temp:maxf_temp)
!!$   enddo
!!$enddo
!!$open(unit=4414)
!!$do j = minf_temp, maxf_temp
!!$   write(4414,'(100F20.10)') grid%gridr(j), temp(j)
!!$enddo
!!$close(unit=4414)

!!$ Create basis of 1-electron pseudostates from Laguere basis and CI coef.
!!$          if(la .le. data_in%la_core) then
!!$             jstart = data_in%npr_core(la) - la + 1
!!$          else 
!!$             jstart = 1
!!$          endif

          jstart = 1
          do j=jstart,jstart+data_in%nst(ma,ipar)-1 ! nd             
             energy = w(j)
             N_oneel = N_oneel + 1
             tmp_sign_CI = SUM(CI(1:nd,j))
             if( tmp_sign_CI .lt. 0d0 ) then
                CI(1:nd,j) = -CI(1:nd,j)
             endif
             call construct_st(TargetStates%b(N_oneel),hlike,dble(ma),ipar,0.5d0,energy,j,nd,CI(1:nd,j),no(1:nd),mo(1:nd))
!!$ NOTE here we might need to acount for degeneracy of the molecule energy levels for ma .ne. 0
!!$ this could lead to introdicing addtonal state with -ma 
             if( ma .ne. 0 ) then
                N_oneel = N_oneel + 1
                call  construct_st(TargetStates%b(N_oneel),hlike,-dble(ma),ipar,0.5d0,energy,j,nd,CI(1:nd,j),no(1:nd),-1*mo(1:nd))
             endif

!!$             temp(:) = 0d0
!!$             maxf_temp = 0
!!$             minf_temp = grid%nr
!!$
!!$             do i=1,grid%nr
!!$                
!!$                do n=1,nd
!!$                   minf_bj = get_minf(bst%b(n))
!!$                   maxf_bj = get_maxf(bst%b(n))
!!$                   if(minf_bj .le. i .and. maxf_bj .ge. i) then
!!$                      maxf_temp = max(maxf_temp, maxf_bj)
!!$                      minf_temp = min(minf_temp, minf_bj)
!!$                      temp(i) = temp(i) + CI(n,j) * value(bst%b(n),i)
!!$                   endif
!!$                enddo
!!$
!!$             enddo
!!$             tmp_sign_CI = SUM(CI(1:nd,j))
!!$             if(  tmp_sign_CI .lt. 0d0 ) then
!!$                temp(:) = -temp(:)
!!$             endif
!!$             
!!$             close(1289)
!!$             open(1289, file='tmp_file_wf')
!!$             
!!$             do i=minf_temp, maxf_temp
!!$                rval =  grid%gridr(i)
!!$                
!!$                tmp = 2d0*rval*(1d0-rval)*exp(-rval)
!!$                
!!$                write(1289,*) rval, temp(i), tmp
!!$                
!!$             enddo
!!$             
!!$             
!!$             print*, "Enter return to continue"
!!$             read(*,*)
!!$             close(1289)
             

          end do
          
          deallocate(H)
deallocate(Hnum)
          deallocate(b)
          deallocate(CI)
          deallocate(w)
          
       enddo
    enddo

stop

    print*
    call sort_by_energy_basis_st(TargetStates)
    call calc_spectro_factors(TargetStates)
    call print_energy_basis_st(TargetStates)
  end subroutine construct_1el_basis_spheroidal2



subroutine H_mterm_num( bi, bj, result )
  !
  ! Numerically integrates the products of two radial basis functions divided by rho+2.
  !
  use sturmian_class
  use grid_radial
  implicit none

  type(sturmian_nr), intent(in) :: bi, bj
  real*8, intent(out) :: result

  integer :: mindex, indmax
  real*8, dimension(grid%nr) :: f
  real*8, dimension(:), pointer :: fi, fj


  fi => fpointer(bi)
  fj => fpointer(bj)
  
  mindex = max( get_minf(bi), get_minf(bj) )
  indmax = min( get_maxf(bi), get_maxf(bj) )

  f(mindex:indmax) = fi(mindex:indmax) * fj(mindex:indmax) * grid%weight(mindex:indmax) / (grid%gridr(mindex:indmax)+2.0)
  result = sum( f(mindex:indmax) )


end subroutine H_mterm_num



subroutine Hnum_element_spheroidal( bi, bj, li, lj, result )
  use input_data
  use sturmian_class
  use grid_radial
  implicit none

  type(sturmian_nr), intent(in) :: bi, bj
  integer, intent(in) :: li, lj
  real*8, intent(out) :: result

  integer :: mindex, indmax, ind, mj
  real*8 :: lambda, R, pi
  real*8, dimension(grid%nr) :: rho, f, dfj, d2fj
  real*8, dimension(:), pointer :: fi, fj


  pi = 2.0*acos(0.0)
  result = 0.0

  mj = get_ang_mom_proj(bj)
  if( li.eq.lj .and. get_ang_mom_proj(bi).eq.mj ) then
     lambda = 2.0 * data_in%alpha( get_ang_mom_proj(bi) )
     R = data_in%Rd

     fi => fpointer(bi)
     fj => fpointer(bj)
     mindex = max( get_minf(bi), get_minf(bj) )
     indmax = min( get_maxf(bi), get_maxf(bj) )
     rho(mindex:indmax) = grid%gridr(mindex:indmax)

     call diff( fj, mindex, indmax, dfj )
     dfj = ( rho**2 + 2.0*rho ) * dfj

     call diff( dfj, mindex, indmax, d2fj )
     d2fj = -pi*R/(2.0*lj+1.0) * ( d2fj + ( 2.0*R*(rho+1.0) - lj*(lj+1.0) - mj**2/(rho**2+2.0*rho) )*fj )
! + 0*0.5*pi*R**2 * (rho**2 + 2.0*rho + 2.0/3.0) * fj )

     f(mindex:indmax) = grid%weight(mindex:indmax) * fi(mindex:indmax) * d2fj(mindex:indmax)
!     f(mindex:indmax) = grid%weight(mindex:indmax) * fi(mindex:indmax) * fj(mindex:indmax) * ( rho(mindex:indmax)**2 + 2.0*rho(mindex:indmax) + 1.0 - (2.0*lj**2+2.0*lj-2.0*mj**2-1.0)/(2.0*lj-1.0)/(2.0*lj+3.0) ) * pi*R**3/2.0 / (2.0*lj+1.0)
     do ind = lj-mj+1,lj+mj
        f = f * ind
     enddo

     result = sum( f(mindex:indmax) )

!!$print*, result
!!$open(unit=1336)
!!$write(1336,'(5(A,I2))') '#   n=', get_k(bj), ' ->', get_k(bi), '   l=', lj, ' -> ', li, '   m=', mj
!!$do ind = mindex,indmax
!!$   write(1336,'(100F20.10)') rho(ind), fj(ind), dfj(ind), d2fj(ind)
!!$enddo
!!$close(unit=1336)
!!$if ( get_k(bj) .eq. 2 ) then
!!$   stop
!!$endif

  endif


end subroutine Hnum_element_spheroidal



subroutine diff( f, lower, upper, df )
! Numerically calculates the derivative of the array stored in sturm
! over the range [mindex, indmax] and stores the resulting array in df.
  use sturmian_class
  use grid_radial
  implicit none

  real*8, dimension(grid%nr), intent(in) :: f
  integer, intent(in) :: lower, upper
  real*8, dimension(grid%nr), intent(out) :: df

  integer :: j, k, mindex, indmax


  df = 0.0

  ! End-point finite difference derivatives. Less accurate.
  mindex = max( lower-1, 1 )   ! Pushes below lower bound if possible.
  indmax = min( upper+1, grid%nr )   ! Pushes above upper bound if possible.
  df(mindex) = ( f(mindex+1)-f(mindex) ) / ( grid%gridr(mindex+1)-grid%gridr(mindex) )
  df(indmax) = ( f(indmax)-f(indmax-1) ) / ( grid%gridr(indmax)-grid%gridr(indmax-1) )
df(indmax) = 0.0

  ! Central finite difference derivative. More accurate.
  df(mindex+1:indmax-1) = (/( (f(j+1)-f(j-1))/(grid%gridr(j+1)-grid%gridr(j-1)) , j=mindex+1,indmax-1 )/)

  ! Correct the points at which the interval is doubled (general formula).
  do j=1,grid%ndouble
     k = j*grid%npdbl
     df(k) = df(k) - ( df(k+1)-df(k-1) )*( grid%gridr(k+1)/2.0-grid%gridr(k)+grid%gridr(k-1)/2.0 )/( grid%gridr(k+1)-grid%gridr(k-1) )   ! Basically, interpolation to an off-centre point.
  enddo


end subroutine diff

