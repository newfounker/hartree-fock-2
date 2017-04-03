!--------------------------------------------------------------------------------------------------------
subroutine rearrange(bst,LMAX,TargetStates)
!!$ OneState is a pointer to one target state, and bst is its single particle basis for  target states
!!$ will return back bst which contains the new single particle basis (an array of sturmian_nr structure)
  
  use grid_radial
  use sturmian_class
  use state_class

  implicit none
  
  type(basis_sturmian_nr), intent(inout) :: bst 
  integer, intent(in):: Lmax   ! max L valus in one-electron basis bst
  type(basis_state), intent(inout):: TargetStates

  integer, dimension(TargetStates%Nmax,LMAX+1):: av_l_newf ! an array that will store  l values 
  integer, dimension(TargetStates%Nmax,LMAX+1):: av_n_newf ! an array that will store  n values 
  integer, dimension(TargetStates%Nmax):: num_n_newf
  type(basis_sturmian_nr):: bst_rearr ! these are sturmian basis (depend on k and l)
  type(state), pointer:: OneState !( a pointer to) one target state only
  integer:: n
  integer:: l, lfn
  real*8, pointer, dimension(:):: fn
  real*8::rCI
  type(sturmian_nr), pointer:: pn, pnnew  !  one-electron orbitals
  integer:: nc, ncm
  integer:: maxl
  integer::minf,maxf, minf_all, maxf_all
  integer:: nst, nstmax
  integer:: tot_numf ! total number of new basis (function) to allocate (for all target states)
  integer:: icntr ! counter variable for indexing bst_rearr%b(bst_rearr)
  integer:: ictr ! counter variable for counting the number of basis to allocate
  integer:: i,j, i1, j1
  integer, dimension(:), allocatable:: ma_in, na_in
  real*8, dimension(:), allocatable:: CI_in
  integer:: nam_in, num ! num = num_n_newf(nst), updated in the loop over each target state
                        ! it is used to allocate no, na_in, CI_in with the size num
  real*8, dimension(grid%nr):: arrfnew
  character(20):: filename
  integer:: icl, icheck, licl, mst, nspm_loc, nicl
!!$------------------------------------------------
  integer:: ii, jj, nsti, nstf, ncm_i, ncm_f, msti, mstf
  integer:: icl_i, icl_f, licl_i, licl_f, nicl_i, nicl_f, nci, ncf
  integer:: n_i, n_f, lfn_i, lfn_f
  real*8::  rCI_i, rCI_f
  real*8:: sum_new
  type(state), pointer:: OneState_i, OneState_f
  type(sturmian_nr), pointer:: pn_i, pn_f
!!$-----------------------------------------------
  
  print*, 'start rearange() for one-electron states'

  tot_numf = 0 ! intialization

  nstmax = TargetStates%Nmax ! number of target states
!!$  print*,'nstmax = ', nstmax
  
  num_n_newf(:) = 0 ! initialiize to 0
  
!!$ Initialize av_l_newf array elements with -1
  av_l_newf(:,:) = -1
  
  
  av_n_newf(:,:) = 0
  
  do nst = 1, nstmax 

     ictr = 0 ! a counter when a new l is encountered in the loop over all the Laguerre basis of a particular target state, initialize ictr to zero in every target state

     OneState => TargetStates%b(nst)

     ncm = get_nam(OneState) ! ncm = nam which is the size of the array na(:), which is also the number of CI coeff
     do nc = 1,ncm ! go through every Lagueere basis (psi_kl(r)) of that target state
!!$---------------------------------------
        rCI = get_CI(OneState, nc)
        if(rCI .eq. 0.0) cycle
!!$----------------------------------------
        n = get_na(OneState,nc,1) ! which is na(nc), where nc = 1,2,...ncm
        pn => bst%b(n) 
        l = get_ang_mom(pn) 
!        print*, '>>>>> nst, n, l:', nst, n, l, rCI

        icheck = 0  ! check if this value of  l  is already counted
        do icl=1,ictr
           if(l .eq. av_l_newf(nst,icl)) then
              icheck = 1   ! found it counted
              exit
           endif
        enddo
        
        if(icheck .eq. 0) then !  found new l value
           ictr = ictr + 1
!           print*, ictr, l
           tot_numf = tot_numf + 1
           av_l_newf(nst,ictr) = l                         
           av_n_newf(nst,ictr) = tot_numf                         
        endif
        
     end do ! end nc loop
     num_n_newf(nst) = ictr ! record the total number of different l basis for that target state
     
  end do ! end nst loop
  
  write(*,'("rearrange.f90:  number of sp orbitals: tot_numf = ",I5)') tot_numf

!!$---------------------------------------------
  call new_basis_nr(bst_rearr,tot_numf) ! allocate array of structure sturmian_nr (for the new basis)
  ! for all the target states, with tot_numf as the array (basis) size
  nspm_loc = 0
  do nst = 1, nstmax ! go through every target state
     
     OneState => TargetStates%b(nst)
     mst = get_angmom_proj(Onestate)
     ncm = get_nam(OneState)
     
     
     do icl=1, num_n_newf(nst) ! go through every l of the target state
        licl = av_l_newf(nst, icl)
        nicl = av_n_newf(nst, icl)

        nspm_loc = nspm_loc + 1

        if(nicl .ne. nspm_loc) then
           print*, '*** rearrange.f90'
           print*, '*** nspm_loc .ne. nicl:', nspm_loc, nicl
           stop
        endif
        
        arrfnew(:) = 0d0 ! all elements are initialized to zero
        
        maxf_all = 1
        minf_all = grid%nr
        
        do nc = 1, ncm 
           n = get_na(OneState,nc,1)
           pn => bst%b(n)
           lfn = get_ang_mom(pn)
           
           if(lfn .eq. licl) then
              
              fn => fpointer(pn) ! means f points to bst%b(n)%f
              minf = get_minf(pn)
              maxf = get_maxf(pn)
              rCI = get_CI(OneState,nc)
              
              arrfnew(minf:maxf) =  arrfnew(minf:maxf) + rCI*fn(minf:maxf)
              
              maxf_all = max(maxf_all, maxf)
              minf_all = min(minf_all, minf)
              
           end if
        end do ! nc
        
!!$ copy arrfnew array for that particular lfn into bst_rearr%b(icntr)
        call init_function(bst_rearr%b(nspm_loc),licl,mst,nspm_loc,minf_all,maxf_all,arrfnew,grid%nr)
                
     end do
     
  end do ! nst
  
  if(nspm_loc .ne. tot_numf) then
     print*, '*** rearrange.f90'
     print*, '*** nspm_loc .ne. tot_numf:', nspm_loc, tot_numf 
     stop
  endif
  
!
!!$------------------------------------------------------------------------------------
!!$ code for ortint(:,:)

  do nsti = 1, nstmax
     OneState_i => TargetStates%b(nsti)
     ncm_i = get_nam(OneState_i)
     msti = get_angmom_proj(Onestate_i)
 
     do icl_i  = 1, num_n_newf(nsti) ! walk through all l for that target state
        licl_i = av_l_newf(nsti, icl_i)
        nicl_i = av_n_newf(nsti, icl_i)

           do nstf = 1, nstmax
              OneState_f => TargetStates%b(nstf)
              ncm_f = get_nam(OneState_f)
              mstf = get_angmom_proj(Onestate_f)

!!$ check the same M-state values
              if(msti .ne. mstf) then
                 do  icl_f  = 1, num_n_newf(nstf) ! walk through all l for that target state
                    licl_f = av_l_newf(nstf, icl_f)
                    nicl_f = av_n_newf(nstf, icl_f)

                    bst_rearr%ortint(nicl_f,nicl_i) = 0d0

                 end do
                 cycle  ! over this value of nstf
              end if

              do icl_f  = 1, num_n_newf(nstf) ! walk through all l for that target state
                 licl_f = av_l_newf(nstf, icl_f) !!
                 nicl_f = av_n_newf(nstf, icl_f) !!

                 if(nicl_f .gt. nicl_i) cycle


!! put code here 
!!$-------------------------------------------------------------------------------------------
                 sum_new = 0d0

                 do nci = 1, ncm_i
                    n_i = get_na(OneState_i,nci, 1)
                    pn_i => bst%b(n_i)
                    lfn_i = get_ang_mom(pn_i)

                    if(lfn_i .ne. licl_i) cycle
                    rCI_i = get_CI(OneState_i, nci)

                    do ncf = 1, ncm_f
                       n_f = get_na(OneState_f,ncf, 1)
                       pn_f => bst%b(n_f)
                       lfn_f = get_ang_mom(pn_f)

                       if(lfn_f .ne. licl_f) cycle
                       rCI_f = get_CI(OneState_f, ncf)

                       sum_new = sum_new + rCI_i * rCI_f * bst%ortint(n_i,n_f)

                    end do ! ncf

                 end do ! nci

!!$------------------------------------------------------------------------------------------ 
                 bst_rearr%ortint(nicl_f,nicl_i) = sum_new
                 bst_rearr%ortint(nicl_i,nicl_f) = sum_new
              end do ! icl_f

           end do ! nstf

  end do ! icl_i
end do ! nsti

!!$--------------------------------------------------------------------------------------
!!$ This is the code to overwride the bst (copying bst_rearr to bst)

  call destruct(bst) ! deallocate old bst (old basis) first
  call new_basis(bst, nspm_loc) ! allocate nspm_loc new space for a new created bst
  
  call copy_basis(bst,bst_rearr)
! nspm =  basis_size(bst)

  call destruct(bst_rearr)
!!$-------------------------------------------------------------------------------------
!!$ Check if it's correct
!   do ii = 1, nspm
!     do jj = 1, nspm
!        print*,'bst%ortint(',ii,',',jj,') = ', bst%ortint(ii,jj)
!     end do
!  end do
!!$ -------------------------------------------------------------------------------
!!$ This is code to update the TargetStates links to the new basis

  do nst = 1, nstmax ! go through every target state
     
     OneState => TargetStates%b(nst)

     mst = get_angmom_proj(Onestate)
     num = num_n_newf(nst)
     
     allocate(na_in(num), CI_in(num), ma_in(num)) 
     
     nam_in = num
     na_in(1:num) = av_n_newf(nst,1:num)
     CI_in(1:num) = 1d0
     ma_in(1:num) = mst
     
     call modify_CI_state(OneState,nam_in,na_in,CI_in, ma_in)
     
     deallocate(ma_in)
     deallocate(na_in)
     deallocate(CI_in) 
     
  end do ! nst

  print*, 'finish  rearrange' 

end subroutine rearrange


!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
