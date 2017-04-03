module state_class_extra

  subroutine calc_spectro_factors( self )
    !
    ! Calculates the spectroscopic factors of each state for each subset of basis functions
    ! and guesses the l (in TargetState%l array) and thus n (in each state%n) of the state.
    !
    ! JS 6/12/11
    !
    use input_data
    use sturmian_class
    use one_electron_func

    implicit none

    type(basis_state), intent(inout) :: self

    type(state), pointer :: stateObj
    integer :: stateNum, n, l, m, par, numFunctions, funcDown, funcRight, ind
    integer, dimension(:), allocatable :: funcIndices!, lArray
    real*8 :: energy
    real*8, dimension(0:data_in%latop) :: spectroVec
    real*8, dimension(:,:), allocatable :: overlapMatrix


    open( 1336, file='states.core_parts' )
    write(1336,*) 'state n   l   m  par  label  Energy    S(l)   sum(S)   S(0)    S(1)  . . .'


    do stateNum = 1, basis_size_st(self)   ! Loop through all the states.

       stateObj => self%b(stateNum)   ! Pointer to the current state object.
       energy = get_energy_st(stateObj)   ! Energy of the state for printing purposes.
       m = get_angmom_proj(stateObj)   ! Angular momentum projection of the state.
       par = get_par_st(stateObj)   ! Parity (+1 or -1) of the state.
       numFunctions = get_nam_st(stateObj)   ! Number of functions used to construct this state.
!print*
!print*, 'stateNum:', stateNum, 'numFunctions:', numFunctions
!print*, 'na:', stateObj%na
!print*, 'ma:', stateObj%ma
       allocate( funcIndices(numFunctions), overlapMatrix(numFunctions,numFunctions) )
       funcIndices = stateObj%na
       overlapMatrix = bst%ortint( funcIndices, funcIndices )   ! Overlap matrix of the relevant functions.

       spectroVec(:) = 0.0d0
       do funcDown = 1, numFunctions   ! Loop through the basis functions comprising the state
          if( data_in%calculation_type.eq.0 .or. data_in%calculation_type.eq.1 ) then
             l = get_ang_mom(bst%b(funcIndices(funcDown)))   ! Contributing angular momentum.
          elseif( data_in%calculation_type .eq. 2 ) then
             l = bst%lArray( funcIndices(funcDown) )
          endif

          do funcRight = 1, numFunctions   ! Loop through the functions again.

             spectroVec(l) = spectroVec(l) + get_CI(stateObj,funcDown) * get_CI(stateObj,funcRight) * overlapMatrix(funcDown,funcRight)   ! Each function's contribution to the spectroscopic factor.

          enddo ! funcRight
       enddo ! funcDown

       deallocate( funcIndices, overlapMatrix )


       l = maxloc( spectroVec, 1 ) - 1   ! Guess at l from the largest spectroscopic factor.
       n = l + 1   ! Ensures n is greater than l.
       self%b(stateNum)%l = l

       do ind = 1, stateNum-1   ! Search through previous states for similar l & m to generate n.
          if ( self%b(ind)%l == l .and. get_angmom_proj(self%b(ind)) == m ) n=n+1
       enddo

       stateObj%n = n

       stateObj%label = make_label(n,l,m,par)

       write( 1336, '(5I4,A7,F10.5,100F8.4)') stateNum, n, l, m, par, stateObj%label, energy, spectroVec(l), sum(spectroVec), spectroVec

    enddo ! stateNum

    close(1336)


  end subroutine calc_spectro_factors



!
!!$
!!$
!!$  state_l,state_r  are one-electron target states !>> to do make it work with two-electron states too
  function ovlp_st(state_l,state_r)

    use sturmian_class
    use one_electron_func
    use ovlpste1me

    implicit none

    real*8:: ovlp_st

    type(state), intent(in):: state_l, state_r

    real*8, dimension(:,:), pointer:: ortint
    real*8:: tmp_l, tmp_r, tmpsum, ovlpsp1, ovlpsp2
    integer:: i, j, m, n, ni, nm,  n1l, n2l, n1r, n2r, li, lm


    ovlp_st = 0d0



    if(nint(2*state_l%m) .ne. nint(2*state_r%m)) then
!       print*, 'm: l:', state_l%m
!       print*, '   r:', state_r%m
       return
    endif
    if(state_l%parity .ne. state_r%parity) then
!       print*, 'par: l:',  state_l%parity
!       print*, '     r:',  state_r%parity
       return
    endif
    if(nint(2*state_l%spin) .ne. nint(2*state_r%spin)) then
!       print*, 'spin: l:', state_l%spin
!       print*, '      r:', state_r%spin
       return
    endif

    if(state_r%hlike .ne. state_l%hlike) then
       print*,'state_class.f90:  ovlp_st(): state_r%hlike .ne. statel%hlike'
       print*,'stop'
       stop
    endif

!    if(state_r%hlike) then
       ortint => bst%ortint
!    endif


    tmpsum = 0d0

    if(state_r%hlike) then ! one-electron states
        do i=1,state_l%nam
          ni = state_l%na(i)   ! index to a Sturmian type function with fixe dvalue of angular moemntum
!          print*, 'ni=', ni
          li = get_ang_mom(bst%b(ni))
          tmp_l = get_CI(state_l,i)
          do m=1,state_r%nam
             nm = state_r%na(m)  ! index to a Sturmian type function with fixe dvalue of angular moemntum
!             print*, 'nm=', nm
             lm = get_ang_mom(bst%b(nm))
             if( li .ne. lm )  cycle
             ovlpsp1 = ortint(ni,nm)
             if(ovlpsp1 .eq. 0d0) cycle
             tmp_r = get_CI(state_r,m)
             tmpsum = tmpsum + tmp_l * tmp_r * ovlpsp1
          enddo
       enddo
    else  ! two-electron states
       do i=1,state_l%nam
          n1l = state_l%na(i)
          n2l = state_l%nb(i)
          tmp_l = get_CI(state_l,i)
          do j=1,state_r%nam
             n1r = state_r%na(j)
             n2r = state_r%nb(j)
             tmp_r = get_CI(state_r,j)

             tmpsum = tmpsum + ortint(n1l,n1r)*ortint(n2l,n2r)*tmp_l*tmp_r
             !    tmpsum = tmpsum + ovlpst(n1l,n1r)*ovlpst(n2l,n2r)*tmp_l*tmp_r
             !  print'(">>>",4E15.6)', ovlpst(n1l,n1r),ovlpst(n2l,n2r),tmp_l,tmp_r
!             print'(2i3,2(2i5,F15.5))', i,j, n1l, n1r,ortint(n1l,n1r), n2l, n2r,ortint(n2l,n2r)

          enddo
       enddo

    endif
    ovlp_st = tmpsum

  end function ovlp_st
!
!!$  this is a function that calculate matrix element for one-electron Hamiltonian (H2+)
!!$  state_l,state_r  are one-electron target states
!!$  it relies on special form of underlying Laguerre basis (2l+1)
  function H1el_st(state_l,state_r)

    use sturmian_class
    use one_electron_func

    implicit none

    real*8:: H1el_st
    type(state), intent(in):: state_l, state_r  ! it is one-electron states

    real*8, dimension(:,:), pointer:: ortint
    real*8:: tmp_l, tmp_r, tmpsum, result
    integer:: ma, i, m, n, ni, nm


    ortint => bst%ortint

    H1el_st = 0d0

    if(state_l%m .ne. state_r%m) return
    if(state_l%parity .ne. state_r%parity) return

    ma = state_l%m

    if(state_l%hlike) then
    else
       print*,'>>> state_class.f90: H1el_st(...)'
       print*,'Two-electron states have not been coded yet.'
       print*, 'stop'
       stop
    endif

    tmpsum = 0d0
    do i=1,state_l%nam
       ni = state_l%na(i)
       tmp_l = get_CI(state_l,i)
       do m=1,state_r%nam
          nm = state_r%na(m)
          tmp_r = get_CI(state_r,m)

          call Hlagorb(bst,ni,nm,ma,result)

          tmpsum = tmpsum +  result * tmp_l * tmp_r
!          print'(4i5,3E15.5)', i,ni,m,nm, tmp_l, tmp_r,result
       enddo
    enddo

     H1el_st = tmpsum

   end function H1el_st

end module state_class_extra
