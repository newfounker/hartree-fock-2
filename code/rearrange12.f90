!--------------------------------------------------------------------------------------------------------
subroutine rearrange12_st(TargetStates1el, TargetStates2el,nicm,Nmax,Nmax1el)

  use sturmian_class
  use state_class
  use one_electron_func
  use grid_radial
  use ovlpste1me

  implicit none
  
  type(basis_state):: TargetStates1el, TargetStates2el 
  integer, intent(in):: nicm,Nmax,Nmax1el

  integer:: nst, ma, inc, n1, n2, i, istcore, nstcore, nam, nst1el, ipar, j, nd, ncm, is, nstp, nspbst
  real*8:: energy, rma, tmp
!!$ representation of f.c. 1el states via orginal one-electron states
  integer, dimension(:,:), allocatable:: igiven
  real*8, dimension(:,:), allocatable:: CIigiven
  integer, dimension(:), allocatable:: igiven_max

!!$ representation of 2el states via f.c. one-electron states
  integer, dimension(:), allocatable:: manst1el, iparnst1el, itmpmo
  integer, dimension(Nmax):: nstton1
  integer, dimension(Nmax,nicm):: nstton2
  integer, dimension(Nmax,nicm):: nstA
  real*8, dimension(Nmax,nicm):: CInstA
  logical iscore, islarger, itmplog
  type(basis_state):: tmp1elst

  integer, dimension(:), allocatable:: no1,no2,mo1,mo2, phase
  real*8, dimension(2*Nmax1el):: CI
  integer, dimension(:), allocatable:: nobst, mobst
  real*8, dimension(:), allocatable:: CInobst
  integer:: itmpmax
  

!  will need to change the size of arrays 2*Nmax1el to the number of f.c. states 
!  allocate(igiven(2*Nmax1el,Nmax1el), CIigiven(2*Nmax1el,Nmax1el), igiven_max(2*Nmax1el), manst1el(2*Nmax1el), iparnst1el(2*Nmax1el), itmpmo(Nmax1el))
!  allocate(nstton1(Nmax), nstton2(Nmax,nicm))



  nstton1(:) = 0 
  nstton2(:,:) = 0
  CI(:) = 0d0

  CInstA(:,:) = 0d0
  nstA(:,:) = 0

!!$  first determine the number of frozen-core orbitals

  nst1el = nicm
  
  do nst=1, Nmax
     
     do istcore=1,nicm  ! go through the core one-electron states
        
        nstcore = TargetStates2el%ncore(istcore)
        
        nam = get_nam(TargetStates2el%b(nst))
        
        i = 0 ! start the counter for f.c. 1el. state representation here
        j = 0  ! to deal with symetric configuration
        do inc=1,nam
           n1 = get_na(TargetStates2el%b(nst),inc,1)  ! index to one-electron TargetStates
           if(n1 .ne. nstcore) cycle

           n2 = get_na(TargetStates2el%b(nst),inc,2)  


           itmplog = islarger(n2,istcore)

           if( n1 .eq. n2) then
              j = n1
              cycle
           elseif(itmplog) then 
              ! contniue to build new f.c. state
           else
              cycle ! to avoid double counting of configurations 
           endif
           
           ! f.c. orbital representation is done below
           ! only states n2 with the same value of magnetic quantum number  and parity will be here           
           ! need only to count the first occurance to add a new f.c. 1el state
           if(i .eq. 0 ) then
              nst1el = nst1el + 1  
           endif

           i = i + 1

        enddo

       
     enddo

  enddo

  allocate(igiven_max(nst1el),CIigiven(nst1el,Nmax1el),igiven(nst1el,Nmax1el))
  igiven_max(:) = 0
  igiven(:,:) = 0
  CIigiven(:,:) = 0d0
  allocate(manst1el(nst1el), iparnst1el(nst1el), itmpmo(nst1el))
  manst1el(:) = 0
  iparnst1el(:) = 0
  allocate(no1(nst1el),no2(nst1el),mo1(nst1el),mo2(nst1el), phase(nst1el))
  no1(:) = 0
  no2(:) = 0
  mo1(:) = 0
  mo2(:) = 0
  phase(:) = 0



print '("rearrange12(): Nmax,Nmax1el,nicm,nst1el: ",4i7)',   Nmax, Nmax1el, nicm, nst1el


  do istcore=1,nicm  ! go through the core one-electron states
     
     nstcore = TargetStates2el%ncore(istcore)
     manst1el(istcore) = get_angmom_proj(TargetStates1el%b(nstcore))
     iparnst1el(istcore) = get_par_st(TargetStates1el%b(nstcore))
     
  enddo

  nst1el = nicm
 
  do nst=1, Nmax

     do istcore=1,nicm  ! go through the core one-electron states
!        print*
        nstcore = TargetStates2el%ncore(istcore)
        
        nam = get_nam(TargetStates2el%b(nst))
        
        i = 0 ! start the counter for f.c. 1el. state representation here
        j = 0  ! to deal with symetric configuration
        do inc=1,nam
           n1 = get_na(TargetStates2el%b(nst),inc,1)  ! index to one-electron TargetStates
           if(n1 .ne. nstcore) cycle

           n2 = get_na(TargetStates2el%b(nst),inc,2)  

!           print*, '--', n1, n2

           itmplog = islarger(n2,istcore)

           if( n1 .eq. n2) then
              j = n1
              nstA(nst,istcore) = istcore
              CInstA(nst,istcore) = get_CI(TargetStates2el%b(nst),inc)
              cycle
           elseif(itmplog) then 
              ! contniue to build new f.c. state
           else
              cycle ! to avoid double counting of configurations 
           endif
           
           ! f.c. orbital representation is done below
           ! only states n2 with the same value of magnetic quantum number  and parity will be here           
           ! need only to count the first occurance to add a new f.c. 1el state
           if(i .eq. 0 ) then
              nst1el = nst1el + 1  
              nstton2(nst,istcore) = nst1el
              manst1el(nst1el) = get_angmom_proj(TargetStates1el%b(n2))
              iparnst1el(nst1el) =   get_par_st(TargetStates1el%b(n2))
           endif

           i = i + 1
           igiven(nst1el,i) = n2
           CIigiven(nst1el,i) = get_CI(TargetStates2el%b(nst),inc)
!           print'(10i5)', nst, n1,nst1el,n2, i
           
           igiven_max(nst1el) = i

        enddo

        if( i .gt. 0 .and. j .ne. 0) then
! if j = 0 then there was no symmetric configuration
! if i = 0 then there was no f.c. orbital for thsi core orbital 
           i = i + 1
           igiven_max(nst1el) = i
           igiven(nst1el,i) = j
           CIigiven(nst1el,i) = CInstA(nst,istcore)/2d0
           nstA(nst,istcore) = 0  ! set to zero to avoid including it in two places
!           print*,'adding symmetric configuration:', j, j
!           print*, 'igiven_max(nst1el)=', igiven_max(nst1el)
        endif
     enddo

  enddo
!  print*

  call new_basis_st(tmp1elst,nst1el,.true.,0)

!  print*, 'nicm = ', nicm
  do i=1,nicm  ! go through the core one-electron states
     
     nstcore = TargetStates2el%ncore(i)
     call copy_st(tmp1elst%b(i),TargetStates1el%b(nstcore))
     tmp1elst%Nstates = i

  enddo


  nspbst = basis_size(bst)
!  print*, 'nspbst=',nspbst
  allocate(nobst(nspbst),mobst(nspbst),CInobst(nspbst))
  nobst(:) = 0 
  mobst(:) = 0
  CInobst(:) = 0d0


!
  do i=nicm+1,nst1el

!!$  need to transform from representation of frozen-core orbitals in terms of one-electron 
!!$  target states TargetStates1el to reprresentation via underlying one-electron functions with 
!!$  fixed orbital angular momentum given by  bst
!     print*,'i=',i
     nd = igiven_max(i)
     call make_nomo(Targetstates1el,bst,nd,CIigiven(i,1:nd),igiven(i,1:nd),itmpmo(i),nspbst,nobst,mobst,CInobst,itmpmax)
!!$--------
     ipar = iparnst1el(i)
     rma = manst1el(i)
     energy = 0d0
     j = 0
     nd = itmpmax
     tmp1elst%Nstates = i
     call  construct_st(tmp1elst%b(i),.true.,rma,ipar,0.5d0,energy,i,nd,CInobst(1:nd),nobst(1:nd),mobst(1:nd))
  enddo

!  print*, 'size of tmp1elst:',  basis_size_st(tmp1elst)


  call rearrange(bst,get_max_L(bst),tmp1elst)

  do nst=1, Nmax
     
!!$ determine how many antisymmetric configurations are used for this state
!!$ it is equal to the number of core orbital used in description of this state
     is = nint(TargetStates2el%b(nst)%spin)
     i = 0
     phase(:) = (-1)**(is)
     do istcore=1,nicm
        
        if(nstA(nst,istcore) .ne. 0) then   ! symmetric configuration
           i = i+ 1 
           no1(i) = istcore 
           no2(i) =  istcore 
           nstcore = TargetStates2el%ncore(istcore)
           mo1(i) = get_angmom_proj(TargetStates1el%b(nstcore))
           mo2(i) = mo1(i)
           CI(i) = CInstA(nst,istcore)           
        endif
        if(nstton2(nst,istcore) .ne. 0) then  ! f.c. configuration
           i = i+ 1 
           no1(i) = istcore 
           nst1el = nstton2(nst,istcore)
           no2(i) = nst1el
           nstcore = TargetStates2el%ncore(istcore)
           mo1(i) = get_angmom_proj(TargetStates1el%b(nstcore))
           mo2(i) = manst1el(nst1el)
           if(no1(i) .eq. no2(i)) then
              CI(i) = CIigiven(istcore,istcore)
           else
              CI(i) = sqrt(2d0)  ! as per  setupCI() 1d0
           endif
        endif
     enddo
     ncm = i

     call setupCI(TargetStates2el%b(nst),ncm,CI,no1,mo1,no2,mo2,phase)        
     
     
  enddo

  call destruct_basis_st(TargetStates1el)
  call new_basis_st(TargetStates1el,tmp1elst%Nstates,.true.,0)
  print*, 'tmp1elst%Nstates:', tmp1elst%Nstates
  TargetStates1el%Nmax = tmp1elst%Nmax
  TargetStates1el%Nstates = tmp1elst%Nstates
  do i=1,tmp1elst%Nstates
     call copy_st(TargetStates1el%b(i), tmp1elst%b(i))
  enddo

!!$ Note  arrays ovlpst(:,:) and e1me(:,:)  might need redefinition (change of diminsions) as number of one-elctron target states is changed and can be larger than previously defined.
!!$! 1el
!!$!  print*, '1el'
!!$  do nst=1,TargetStates1el%Nmax
!!$     do nstp=1,nst
!!$        tmp = ovlp_st(TargetStates1el%b(nst),TargetStates1el%b(nstp))
!!$        ovlpst(nst,nstp) = tmp
!!$        ovlpst(nstp,nst) = tmp
!!$!        print*, '3. ', nst, nstp, tmp
!!$     enddo
!!$  enddo


  return
end subroutine rearrange12_st
!!$
!!$ Check if this 1el state n is a core 1el state
function iscore(n)

  use state_class
  use target_states

  implicit none

  integer, intent(in):: n
  logical:: iscore
  integer:: ic, ncst

  iscore = .false.
  do ic=1,TargetStates2el%nicm
     ncst = TargetStates2el%ncore(ic)
     if(n .eq. ncst) then
        iscore = .true.
        exit
     endif
  enddo

  return
end function iscore
!!$
!!$ Check if this 1el state n is a core 1el state and if its  number is larger that core state number k
function islarger(n,k)

  use state_class
  use target_states

  implicit none
  
  logical:: islarger
  integer, intent(in):: n,k
  logical  iscore
  integer:: ic, ncst

  iscore = .false.
  islarger = .false.
  do ic=1,TargetStates2el%nicm
     ncst = TargetStates2el%ncore(ic)
     if(n .eq. ncst) then
        iscore = .true.
        exit
     endif
  enddo

  if(iscore) then 
     if(ic .gt. k) then
        islarger = .true.
!        print*, '!!', n, ic,k
     endif
  else
     islarger = .true.
!     print*, '>>', n, ic,k
  endif

  return
end function islarger

!!$-----------------------------------------------------------------------------------------------------------------
!!$  transform from repreesentation of frozen-core orbitals in terms of one-electron 
!!$  target states TargetStates1el to reprresentation via underlying one-electron functions with 
!!$  fixed orbital ongular momentum given by  bst
subroutine  make_nomo(TargetStates1el,bst,nd,CIigiven,igiven,itmpmo,nspbst,nobst,mobst,CInobst,itmpmax)

  use sturmian_class
  use state_class
  use ovlpste1me

  implicit none
  
  type(basis_state), intent(in):: TargetStates1el
  type(basis_sturmian_nr), intent(in):: bst   ! has fixed  angular momentum 
  integer, intent(in):: nd
  real*8, dimension(nd), intent(in):: CIigiven
  integer, dimension(nd), intent(in):: igiven
  integer, intent(in)::  itmpmo, nspbst
  integer, dimension(nspbst), intent(out):: nobst, mobst
  real*8, dimension(nspbst), intent(out):: CInobst
  integer, intent(out):: itmpmax

  integer:: i, n,nst, mst, nam, j, np, nstp
  real*8:: CI1el, CI2el, tmp, CI2elp
  real*8, dimension(nspbst):: CItmp

  nobst(:) = 0 
  mobst(:) = 0
  CInobst(:) = 0d0
  itmpmax = 0

  CItmp(:) = 0d0

  mobst(:) = itmpmo  !! they have all the same value m

  mst = itmpmo
  
  do n=1,nd

     nst = igiven(n)  ! index to one-electron target states  TargetStates1el
!     print*, 'n,nst=', n,nst
     CI2el = CIigiven(n)  ! this is coef. with which those one-electron states contribute this the current FC orbital
     
     nam = get_nam(TargetStates1el%b(nst))  ! representation of target state nst  via one-electron  bst  basis
     do i=1,nam
        
        j = get_na(TargetStates1el%b(nst),i)   ! index to one-electron basis function
        CI1el = get_CI(TargetStates1el%b(nst),i)

        CInobst(j) = CInobst(j) + CI1el*CI2el
!        write(*,'(4i5,2E15.5)') n,nst,i,j, CI1el,CI2el
     enddo

  enddo


  i = 0
  do j=1,nspbst

     if(CInobst(j) .eq. 0 ) cycle
     
     i = i + 1
     
     CItmp(i) = CInobst(j)
     nobst(i) = j
     
!     print*,'make_nomo: ***', j, get_angmom(bst%b(j)), CItmp(i)

  enddo
  
!  print*, '>>>', i, nspbst

  CInobst(:) = 0d0
  CInobst(1:i) = CItmp(1:i)

  itmpmax = i

  
  
!!$!test
!!$
!!$  tmp = 0d0
!!$
!!$ do n=1,nd
!!$
!!$     nst = igiven(n)  ! index to one-electron target states  TargetStates1el
!!$     CI2el = CIigiven(n)  ! this is coef. with which those one-electron states contribute this the current FC orbital
!!$
!!$     
!!$     do np=1,nd
!!$        
!!$        nstp = igiven(np)  ! index to one-electron target states  TargetStates1el
!!$        CI2elp = CIigiven(np)  ! this is coef. with which those one-electron states contribute this the current FC orbital
!!$        
!!$        tmp = tmp + CI2el*CI2elp*ovlpst(nst,nstp)
!!$
!!$     enddo
!!$
!!$  enddo
!!$
!!$  print*, 'tmp=',tmp
!!$
!!$
!!$
!!$
!!$  tmp = 0d0
!!$
!!$ do n=1, itmpmax
!!$
!!$     nst = nobst(n)
!!$     CI2el = CInobst(n)  ! this is coef. with which those one-electron states contribute this the current FC orbital
!!$
!!$     
!!$     do np=1, itmpmax
!!$        
!!$        nstp = nobst(np)
!!$        CI2elp = CInobst(np)  ! this is coef. with which those one-electron states contribute this the current FC orbital
!!$        
!!$        tmp = tmp + CI2el*CI2elp * bst%ortint(nst,nstp)
!!$
!!$     enddo
!!$
!!$  enddo
!!$
!!$  print*, 'tmp=',tmp



end subroutine make_nomo

