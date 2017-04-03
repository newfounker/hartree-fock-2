subroutine structure12

  use input_data
  use sturmian_class
  use state_class
  use target_states
  use ovlpste1me
  use one_electron_func_mod

  implicit none

  integer:: lorbmax
  parameter( lorbmax = 20 )

  character(LEN=20):: target
  integer:: Mmax
  integer, dimension(:,:,:), allocatable:: nstate

  integer:: ma, ip, is, itmp, inc, labot, latop, nsp, ncm, nc, ncp, l,m, Nmax, Nmax1el, nst, nstp
  real*8:: rma, ris,  en_ion
  integer, dimension(:), allocatable:: phase
  logical:: hlike
  real*8, allocatable,dimension(:):: alpha
  integer, allocatable,dimension(:):: nps
  type(input):: dataF5
  integer, allocatable,dimension(:):: no1,no2, mo1, mo2
  real*8:: tmp, resultb, resultH

  integer:: matz, ierr, i
  real*8, dimension(:,:), allocatable:: H, b, CI
  real*8, dimension(:), allocatable:: w
  integer:: nicm, ni, icheck
  integer, allocatable,dimension(:):: ncore, mncore
  integer:: nsp1,nsp2, nsp1p, nsp2p, mnsp1, mnsp2, mnsp1p, mnsp2p, Nmax_bound, n
  integer:: l12max, ltmp
  integer, dimension(0:lorbmax):: nk1, nk2
  integer:: l_ion_core
  integer, dimension(0:lorbmax):: n_ion_core
  integer:: i_st  ! =0 to use s.p. basis for 2 el. diagonaliation, =1 to use one-electron target states for diagonalization
  real*8:: tmp_sign_CI

  hlike = .false.

  nk1(:) = 0
  nk2(:) = 0
  n_ion_core(:) = 0

  open(3,file='F5')

  print*
  print*,'Start CI routines for quasi two-electron di-atomic molecule'
  print*

  read(3,*) target
  write(*,*) 'target'
  write(*,*) target

  read(3,*) Mmax
  write(*,*) 'Mmax abs value'
  write(*,*) Mmax

  allocate(nstate(0:Mmax,-1:1,0:1))
  nstate(:,:,:) = 0

  read(3,*) (nstate(l,1,0), l=0,Mmax) ! positive parity, singlet
  write(*,*) 'nstate(l,1,0) : singlet, positive parity'
  write(*,*) '     ', (nstate(l,1,0), l=0,Mmax)
  read(3,*) (nstate(l,-1,0), l=0,MMax) ! negative parity, singlet
  write(*,*) 'nstate(l,1,0) : singlet, negative parity'
  write(*,*) '     ', (nstate(l,-1,0), l=0,Mmax)
  read(3,*) (nstate(l,1,1), l=0,Mmax) ! positive parity, triplet
  write(*,*) 'nstate(l,-1,1) : triplet, positive parity'
  write(*,*) '     ', (nstate(l,1,1), l=0,Mmax)
  read(3,*) (nstate(l,-1,1), l=0,MMax) ! negative parity, triplet
  write(*,*) 'nstate(l,-1,1) : triplet, negative parity'
  write(*,*) '     ', (nstate(l,-1,1), l=0,Mmax)


  Nmax = 0
  do m=0,Mmax
     do ip=-1,1,2
        do is=0,1
           if(nstate(m,ip,is) .gt. 0) then
              itmp = 1
              if(m .ne. 0) itmp = 2  ! include both positive and negative m states
              Nmax = Nmax +  itmp * nstate(m,ip,is)
           endif
        enddo
     enddo
  enddo
  print*, 'Nmax = ', Nmax


  read(3,*) i_st   !  i_st =0 for sp.p.o. from F5, =1 for one-electrob target states
  write(*,*) i_st, '            */ i_st'


  basis_type = 0  !  data_in%calculation_type, only spherical basis is coded
  hlike = .false.

!!$  Target states were created for one-electron target: one-electon ion, get the ground state energy of the ion and then destroy the old target states
  en_ion = get_energy_st(TargetStates%b(1))
  print*, 'en_ion=', en_ion
  if( i_st .eq. 0 )  call destruct_basis_st(TargetStates)
  call new_basis_st(TargetStates2el,Nmax,hlike,basis_type)

  TargetStates2el%hlike = .false.
  TargetStates2el%en_ion = en_ion

!$$   make s.p. Laguerre basis
  dataF5%calculation_type = 0
  read(3,*) labot, latop
  write(*,'("labot,latop: ",2I5)') labot, latop
  allocate(alpha(0:latop))
  allocate(nps(0:latop))

  read(3,*) (nps(l),alpha(l), l=labot,latop)
  write(*,'("nps(l),alpha(l): ",20(I5,F10.4))') (nps(l),alpha(l), l=labot,latop)

  dataF5%latop = latop
  dataF5%labot = labot
  allocate(dataF5%alpha(0:latop))
  allocate(dataF5%nps(0:latop))
  dataF5%nps(:) = nps(:)
  dataF5%alpha(:) = alpha(:)

  if( i_st .eq. 0 ) then
     call destruct(bst)

     call construct(bst, dataF5)
     nspm =  basis_size(bst) ! number of sturmian finctions
     print*, 'Size of Sturmian basis: nspm=', nspm
  else
     print*
     print*, 'Use target states from one-electron diagonalization'
     print*, 's.p. basis in F5 is ignored'
     print*
  endif
!!$----------------------

  print*
  read(3,*) inc   ! inc=1 for call rearange
  write(*,*) inc, '            */ inc'

  read(3,*) l_ion_core,(n_ion_core(l), l=0,l_ion_core)
  print*,  l_ion_core, (n_ion_core(l), l=0,l_ion_core)
  if(l_ion_core .gt. 20) then
     print*, 'increase size of array n_ion_core(20) to at least l_ion_core=', l_ion_core
     stop
  endif

  print*

  ncm = nspm*(2*latop+1)
  allocate(ncore(ncm),mncore(ncm))


!!$   read CI model only once, might want to change later to read per symmetry as in atomic code.
  read(3,*) l12max   !
  write(*,*) l12max, '           */ loutmax'
  read(3,*) (nk2(i), i=0,l12max)   !
  write(*,*) (nk2(i), i=0,l12max), '           */ nkout'
  read(3,*) (nk1(i), i=0,l12max)   !
  write(*,*) (nk1(i), i=0,l12max), '           */ nkin'
!!$
  if(l_ion_core .gt. l12max) then
     print*,'l_ion_core .gt. l12max', l_ion_core, l12max
     l_ion_core = l12max
     print*, 'redefine l_ion_core to l12max value'
  endif

  nicm = 0
  ncm = -1

  do ma = -Mmax, Mmax     !  M of atom
     rma = ma
     do ip=1,-1,-2    ! parity
        do is = 0,1   ! spin
           ris = is

           if(nstate(abs(ma),ip,is) .le. 0) cycle

           print*

           write(*,'("symmetry: ma,ip,is=",3i5)') ma,ip,is

!!$ Set up list of configurations, the CI model in F5 file is read in config12()
!!$    first find number of configurations
           ncm = 0  ! number of configuration
           if(i_st .eq. 0 ) then
              call config12_tmp(bst,ma,ip,is,l12max,nk1(0:l12max),nk2(0:l12max),l_ion_core,n_ion_core(0:l12max),ncm)  ! call first time to find ncm
           else
              call config12_tmp_st(TargetStates,ma,ip,is,l12max,nk1(0:l12max),nk2(0:l12max),l_ion_core,n_ion_core(0:l12max),ncm)  ! call first time to find ncm
           endif
           print*, ' ncm =', ncm
           if(ncm .le. 0) then
              !              if(myid.eq.0)then
              print*,'Problem: no configurations have been built for this symmetry: ncm=',ncm
              !              endif
              cycle
           endif
           if(allocated(no1)) then
              deallocate(no1,no2,mo1,mo2)
           endif
           allocate(no1(ncm),mo1(ncm),no2(ncm),mo2(ncm))
           no1 (:) = 0
           mo1 (:) = 0
           no2 (:) = 0
           mo2 (:) = 0
           if(i_st .eq. 0 ) then
              call config12(bst,ma,ip,is,l12max,nk1(0:l12max),nk2(0:l12max),l_ion_core,n_ion_core(0:l12max),ncm,no1,mo1,no2,mo2) ! call second time to populate arrays
           else
              call config12_st(TargetStates,ma,ip,is,l12max,nk1(0:l12max),nk2(0:l12max),l_ion_core,n_ion_core(0:l12max),ncm,no1,mo1,no2,mo2) ! call second time to populate arrays
           endif
!!$ END Set up list of configurations


!!$  ------------------------------------------------
!!$     find core orbitals: they are in general given by array nk1(l2);
!!$     but not all of them will be used (due to selection rules). Therefore
!!$     it is better to find core orbitals after list of configurations is set.
!!$     This means that core orbitals is now correspond to array no1(nco)
           do nc=1,ncm
!              print*, no1(nc), no2(nc), mo1(nc), mo2(nc)
              icheck = 0  ! checking if this orbital wa sa;ready included
              do ni=1,nicm
                 if(ncore(ni).eq.no1(nc) .and. mncore(ni).eq.mo1(nc)) then
                    icheck = 1
                 end if
              end do
              if(icheck.eq.0) then
                 nicm = nicm + 1
                 ncore(nicm) = no1(nc)
                 mncore(nicm) = mo1(nc)
!            print*, '      nicm,ncore(nicm),mncore(nicm)', nicm,ncore(nicm), mncore(nicm)
              end if
           enddo
!!$  ------------------------------------------------

!!$   Temporary arrays
           if(allocated(H)) deallocate(H,b,CI,w)
           allocate(H(ncm,ncm),b(ncm,ncm),CI(ncm,ncm),w(ncm))

!!$ Form H matrix

           do nc=1,ncm
              nsp1 = no1(nc)
              mnsp1 = mo1(nc)
              nsp2 = no2(nc)
              mnsp2 = mo2(nc)

              do ncp=nc,ncm
                 nsp1p = no1(ncp)
                 mnsp1p = mo1(ncp)
                 nsp2p = no2(ncp)
                 mnsp2p = mo2(ncp)

                 b(nc,ncp) = 0d0
                 H(nc,ncp) = 0d0

                 if(i_st .eq. 0 ) then
                    call H12me(data_in%Z,data_in%Rd,is,bst,nsp1,nsp2,nsp1p,nsp2p,mnsp1,mnsp2,mnsp1p,mnsp2p,resultH,resultb)
                 else
                    if(data_in%inc .eq. 10 ) then
                       Nmax1el = basis_size_st(TargetStates)
                       call H12me_st_notortog(data_in%Z,data_in%Rd,is,TargetStates,Nmax1el,e1me,ovlpst,bst,nsp1,nsp2,nsp1p,nsp2p,mnsp1,mnsp2,mnsp1p,mnsp2p,resultH,resultb)
                    else
                       Nmax1el = basis_size_st(TargetStates)
                       call H12me_st_notortog(data_in%Z,data_in%Rd,is,TargetStates,Nmax1el,e1me,ovlpst,bst,nsp1,nsp2,nsp1p,nsp2p,mnsp1,mnsp2,mnsp1p,mnsp2p,resultH,resultb)
!                       call H12me_st(data_in%Z,data_in%Rd,is,TargetStates,bst,nsp1,nsp2,nsp1p,nsp2p,mnsp1,mnsp2,mnsp1p,mnsp2p,resultH,resultb)
                    endif
                 endif

                 b(nc,ncp) = resultb
                 b(ncp,nc) = resultb
                 H(nc,ncp) = resultH
                 H(ncp,nc) = resultH

              enddo
           enddo
!!$
!!$          print*, 'H matrix'
!!$          do nc=1,ncm
!!$             write(*,'(10E14.5)') (real(H(nc,i)), i=1,ncm)
!!$          enddo
!!$          print*, 'b matrix'
!!$          do nc=1,ncm
!!$             write(*,'(10E14.5)') (real(b(nc,i)), i=1,ncm)
!!$          enddo


           matz=2
           call  rsg(ncm,ncm,H,b,w,matz,CI,ierr)
           write(*,'("ierr =",I3)') ierr
           print*, ' Energies in a.u.'
           write(*,'(5F15.5)') (real(w(i)), i=1,ncm)
           print*, ' Energies in eV - en_ion'
           write(*,'(5F15.5)') (real(27.2116*(w(i)-en_ion)), i=1,ncm)

           if(allocated(phase)) deallocate(phase)
           allocate(phase(ncm))

           if(nstate(abs(ma),ip,is) .gt. ncm) then
              print*,'H12.f: nstate(ma,ip,is) .gt. ncm :', nstate(abs(ma),ip,is), ncm
              print*, 'increase number of s.p. orbitals or decrease number of states in F5'
              stop
           endif

           phase(:) = (-1)**(is)
           do nc=1,nstate(abs(ma),ip,is)

              ! Mark: Fix sign of coefficients
              tmp_sign_CI = SUM(CI(1:ncm,nc))
              if(  tmp_sign_CI .lt. 0d0 ) then
                 CI(1:ncm,nc) = -CI(1:ncm,nc)
              endif


              TargetStates2el%Nstates = TargetStates2el%Nstates + 1
              call construct_st(TargetStates2el%b(TargetStates2el%Nstates),hlike,rma,ip,ris,w(nc),nc,ncm,CI(1:ncm,nc),no1,mo1,no2,mo2,phase)
           enddo

        enddo ! end is loop
     enddo ! end ip loop
  enddo    ! end ma loop

  print*
  if(TargetStates2el%Nstates .ne. Nmax) then
     print*,'Number of calculated states is less than the number of ordered states'
     print*,'Check the CI configuration list, most likely ncm=0 for some symmetry'
     stop
  endif

  TargetStates2el%nicm = nicm   ! nicm was determined in above CI loops
  allocate(TargetStates2el%ncore(1:nicm))
  allocate(TargetStates2el%mncore(1:nicm))
  TargetStates2el%ncore(1:nicm) = ncore(1:nicm)
  TargetStates2el%mncore(1:nicm) = mncore(1:nicm)
  print*,'nicm=',nicm
  do l=1, TargetStates2el%nicm
     print*,l, TargetStates2el%ncore(l), TargetStates2el%mncore(l)
  enddo

  close(3)

  print*, " Sort states"
  call sort_by_energy_basis_st(TargetStates2el)
  print*, " Print energies"
  call print_energy_basis_st(TargetStates2el)

  Nmax_bound = 0
  do n=1,Nmax
     tmp = get_energy_st(TargetStates2el%b(n))
     if(tmp .lt. 0) then
        Nmax_bound = n
     endif
  enddo
  TargetStates2el%Nmax_bound = Nmax_bound

  do n=1,Nmax
     do m=1,n
!!$       write(*,'("overlap:",2i5,2X,es15.5)') n,m, ovlp_st(TargetStates2el%b(n),TargetStates2el%b(m))
     enddo
  enddo


!!$  print*, ' CI representation of two-electron target states'
!!$  do n=1,Nmax
!!$     do i=1,get_nam(TargetStates2el%b(n))
!!$        write(*,'("n,i, CI(i):",2i5,2X,e15.5)') n,i, get_CI(TargetStates2el%b(n),i)
!!$     enddo
!!$  enddo


  if(i_st .eq. 1 ) then
     if(inc .eq. 1) then
        call rearrange12_st(TargetStates,TargetStates2el,TargetStates2el%nicm,TargetStates2el%Nmax,TargetStates%Nmax)
     endif
  else
     !  need to write rearrange12 here
  endif

  call convert_from_st_to_sp(bst,TargetStates,TargetStates2el,basis_size(bst),TargetStates%Nmax,TargetStates2el%Nmax)

!!$! 2el
  print*
  print*, 'Check orthogonality of target states'
  do nst=1,Nmax
     do nstp=1,nst

        tmp = ovlp_st(TargetStates2el%b(nst),TargetStates2el%b(nstp))

        print'(2i5,E15.6)', nst, nstp, tmp

     enddo
  enddo
  print*


! test orthogonality of the target states and calculate osc.strength.
  do n=1,Nmax
     do m=1,n
!        call osc_st(TargetStates2el%b(n),TargetStates2el%b(m))
!        write(*,'(2i5,2X,es15.5)') n,m, ovlp_st(TargetStates2el%b(n),TargetStates2el%b(m))
     enddo
  enddo

  deallocate(H,b,CI,w,phase)

end subroutine structure12
!-----------------------------------------------------------------
subroutine config12(bstF5,ma,ip,is,l12max,nk1,nk2,l_ion_core,n_ion_core,ncm,no1,mo1,no2,mo2)

  use sturmian_class

  implicit none

  type(basis_sturmian_nr), intent(in):: bstF5
  integer, intent(in):: ma, ip, is
  integer, intent(in):: l12max
  integer, dimension(0:l12max), intent(in):: nk1, nk2
  integer, intent(in):: l_ion_core
  integer, dimension(0:l12max), intent(in):: n_ion_core
  integer, intent(inout):: ncm
  integer, dimension(ncm),intent(inout):: no1, mo1, no2, mo2

  integer ico, nsp1, nsp2, nspm, lnsp1, mnsp1, lnsp2, mnsp2, icheck
  integer:: k1,k2

  nspm =  basis_size(bstF5) ! number of sturmian finctions

  ico = 0
  do nsp1=1,nspm
     lnsp1 =  get_ang_mom(bstF5%b(nsp1))
     do mnsp1 = -lnsp1,lnsp1

        do nsp2=nsp1,nspm
           lnsp2 =  get_ang_mom(bstF5%b(nsp2))
           if((-1)**(lnsp1+lnsp2) .ne. ip) cycle

           do mnsp2 = -lnsp2,lnsp2

!!$--------------------- logic block: include or not config ------------
              if(mnsp1+mnsp2 .ne. ma) cycle
              if(nsp1 .eq. nsp2) then
                 if(mnsp1 .gt. mnsp2) cycle  ! to avoid the same config. counted twice
                 if(mnsp1 .eq. mnsp2) then
                    if(is .eq. 1) cycle  ! symmetry condition
                 endif
              endif

              k1 = get_k(bstF5%b(nsp1))   ! inner orbital
              k2 = get_k(bstF5%b(nsp2))
              if(k1 .gt. nk1(lnsp1) ) cycle
              if(k2 .gt. nk2(lnsp2) ) cycle

              if(lnsp1 .gt. l12max) cycle
              if(lnsp2 .gt. l12max) cycle

              if(l_ion_core .ge. 0) then
                 if( k2 .le. nk1(lnsp2)  .or. k1 .le. n_ion_core(lnsp1)) then
                    continue ! include this configuration
                 else
                    cycle   ! ! exclude this configuration
                 endif
              endif


              if(ncm .gt. 0) then
                 icheck = 1
                 call testsameconfig(ncm,no1,no2,mo1,mo2,ico,nsp1,mnsp1,nsp2,mnsp2,icheck)
                 if(icheck .eq. 0) cycle
              endif

!              print*, nk1(lnsp1), k1,  lnsp1, k2, lnsp2

!!$--------------------- end logic block: include or not config ------------

              ico = ico + 1

              if(ncm .gt. 0) then
                 no1(ico) = nsp1
                 mo1(ico) = mnsp1
                 no2(ico) = nsp2
                 mo2(ico) = mnsp2
!                 print*, nsp1, nsp2, mnsp1, mnsp2
              endif

           enddo
        enddo
     enddo
  enddo

  ncm = ico

return
end subroutine config12
!
!
subroutine config12_tmp(bstF5,ma,ip,is,l12max,nk1,nk2,l_ion_core,n_ion_core,ncm)

  use sturmian_class

  implicit none

  type(basis_sturmian_nr), intent(in):: bstF5
  integer, intent(in):: ma, ip, is
  integer, intent(in):: l12max! l12max
  integer, dimension(0:l12max), intent(in):: nk1, nk2
  integer, intent(in):: l_ion_core
  integer, dimension(0:l12max), intent(in):: n_ion_core
  integer, intent(inout):: ncm

  integer ico, nsp1, nsp2, nspm, lnsp1, mnsp1, lnsp2, mnsp2, icheck
  integer:: k1,k2

  nspm =  basis_size(bstF5) ! number of sturmian finctions

  ico = 0
  do nsp1=1,nspm
     lnsp1 =  get_ang_mom(bstF5%b(nsp1))
     do mnsp1 = -lnsp1,lnsp1

        do nsp2=nsp1,nspm
           lnsp2 =  get_ang_mom(bstF5%b(nsp2))
           if((-1)**(lnsp1+lnsp2) .ne. ip) cycle

           do mnsp2 = -lnsp2,lnsp2

!!$--------------------- logic block: include or not config ------------
              if(mnsp1+mnsp2 .ne. ma) cycle
              if(nsp1 .eq. nsp2) then
                 if(mnsp1 .gt. mnsp2) cycle  ! to avoid the same config. counted twice
                 if(mnsp1 .eq. mnsp2) then
                    if(is .eq. 1) cycle  ! symmetry condition
                 endif
              endif

              k1 = get_k(bstF5%b(nsp1))  ! inner orbital
              k2 = get_k(bstF5%b(nsp2))
              if(k1 .gt. nk1(lnsp1) ) cycle
              if(k2 .gt. nk2(lnsp2) ) cycle

              if(lnsp1 .gt. l12max) cycle
              if(lnsp2 .gt. l12max) cycle

              if(l_ion_core .ge. 0) then
                 if( k2 .le. nk1(lnsp2)  .or. k1 .le. n_ion_core(lnsp1)) then
                    continue ! include this configuration
                 else
                    cycle   ! ! exclude this configuration
                 endif
              endif

!!$--------------------- end logic block: include or not config ------------

              ico = ico + 1

           enddo
        enddo
     enddo
  enddo

  ncm = ico

return
end subroutine config12_tmp
!
!
subroutine config12_tmp_st(TargetStates1el,ma,ip,is,l12max,nk1,nk2,l_ion_core,n_ion_core,ncm)

  use state_class

  implicit none

  type(basis_state), intent(in):: TargetStates1el
  integer, intent(in):: ma, ip, is
  integer, intent(in):: l12max
  integer, dimension(0:l12max), intent(in):: nk1, nk2
  integer, intent(in):: l_ion_core
  integer, dimension(0:l12max), intent(in):: n_ion_core
  integer, intent(inout):: ncm

  integer Nmax, ico, nst1, nst2, ipar1, ipar2, mnst1, mnst2, icheck
  integer:: k1,k2, l1_majconf, l2_majconf

  Nmax =  basis_size_st(TargetStates1el) ! number of  one-el. target states

  ico = 0
  do nst1=1,Nmax
     ipar1 =  get_par_st(TargetStates1el%b(nst1))
     mnst1 = get_angmom_proj(TargetStates1el%b(nst1))
     if(abs(mnst1) .gt. l12max) cycle
     do nst2=nst1,Nmax
        ipar2 =  get_par_st(TargetStates1el%b(nst2))
        mnst2 = get_angmom_proj(TargetStates1el%b(nst2))
        if(abs(mnst2) .gt. l12max) cycle

             if(ipar1*ipar2 .ne. ip) cycle

!!$--------------------- logic block: include or not config ------------
        if(mnst1+mnst2 .ne. ma) cycle
        if(nst1 .eq. nst2) then
           if(is .eq. 1) cycle  ! symmetry condition
        endif

        k1 = get_inum_st(TargetStates1el%b(nst1))
        k2 = get_inum_st(TargetStates1el%b(nst2))
        if(k1 .gt. nk1(abs(mnst1)) ) cycle  ! inner orbital
        if(k2 .gt. nk2(abs(mnst2)) ) cycle

        l1_majconf = get_l_majconf(TargetStates1el%b(nst1))
        l2_majconf = get_l_majconf(TargetStates1el%b(nst2))
        if(l1_majconf .gt. l12max) cycle
        if(l2_majconf .gt. l12max) cycle

        if(l_ion_core .ge. 0) then
!           print*, l1_majconf, k1, n1(l1_majconf), n_ion_core(l1_majconf)
           if( k2 .le. nk1(l2_majconf)  .or. k1 .le. n_ion_core(l1_majconf)) then
              continue ! include this configuration
           else
              cycle   ! ! exclude this configuration
           endif
        endif

!!$--------------------- end logic block: include or not config ------------

        ico = ico + 1

     enddo
  enddo

  ncm = ico

return
end subroutine config12_tmp_st
!
!
subroutine config12_st(TargetStates1el,ma,ip,is,l12max,nk1,nk2,l_ion_core,n_ion_core,ncm,no1,mo1,no2,mo2)

  use state_class

  implicit none

  type(basis_state), intent(in):: TargetStates1el
  integer, intent(in):: ma, ip, is
  integer, intent(in):: l12max
  integer, dimension(0:l12max), intent(in):: nk1, nk2
  integer, intent(in):: l_ion_core
  integer, dimension(0:l12max), intent(in):: n_ion_core
  integer, intent(inout):: ncm
  integer, dimension(ncm),intent(inout):: no1, mo1, no2, mo2

  integer ico, Nmax, nst1, nst2, ipar1, ipar2, mnst1, mnst2, icheck
  integer:: k1,k2,  l1_majconf, l2_majconf

  Nmax =  basis_size_st(TargetStates1el) ! number of one-el. target states

  ico = 0
  do nst1=1,Nmax
     ipar1 =  get_par_st(TargetStates1el%b(nst1))
     mnst1 = get_angmom_proj(TargetStates1el%b(nst1))
     if(abs(mnst1) .gt. l12max) cycle
     do nst2=nst1,Nmax
        ipar2 =  get_par_st(TargetStates1el%b(nst2))
        mnst2 = get_angmom_proj(TargetStates1el%b(nst2))
        if(abs(mnst2) .gt. l12max) cycle

             if(ipar1*ipar2 .ne. ip) cycle

!!$--------------------- logic block: include or not config ------------
        if(mnst1+mnst2 .ne. ma) cycle
        if(nst1 .eq. nst2) then
           if(is .eq. 1) cycle  ! symmetry condition
        endif

        k1 = get_inum_st(TargetStates1el%b(nst1))
        k2 = get_inum_st(TargetStates1el%b(nst2))
        if(k1 .gt. nk1(abs(mnst1)) ) cycle   !  inner orbital
        if(k2 .gt. nk2(abs(mnst2)) ) cycle

        l1_majconf = get_l_majconf(TargetStates1el%b(nst1))
        l2_majconf = get_l_majconf(TargetStates1el%b(nst2))
        if(l1_majconf .gt. l12max) cycle
        if(l2_majconf .gt. l12max) cycle

        if(l_ion_core .ge. 0) then
           if( k2 .le. nk1(l2_majconf)  .or. k1 .le. n_ion_core(l1_majconf)) then
              continue ! include this configuration
           else
              cycle   ! ! exclude this configuration
           endif
        endif

        if(ncm .gt. 0) then
           icheck = 1
           call testsameconfig(ncm,no1,no2,mo1,mo2,ico,nst1,mnst1,nst2,mnst2,icheck)
           if(icheck .eq. 0) cycle
        endif
!!$--------------------- end logic block: include or not config ------------

        ico = ico + 1

        if(ncm .gt. 0) then
           no1(ico) = nst1
           mo1(ico) = mnst1
           no2(ico) = nst2
           mo2(ico) = mnst2
!          print'("st: ",8i5)', nst1, nst2, k1, k2, mnst1, mnst2, get_l_majconf(TargetStates1el%b(nst1)), get_l_majconf(TargetStates1el%b(nst2))
        endif

     enddo
  enddo

  ncm = ico

return
end subroutine config12_st
!
!
subroutine testsameconfig(ncm,no1,no2,mo1,mo2,ico,nsp1,mnsp1,nsp2,mnsp2,icheck)

  integer, dimension(ncm),intent(in):: no1, mo1, no2, mo2
  integer, intent(in):: nsp1,mnsp1,nsp2,mnsp2, ico
  integer, intent(inout):: icheck

  do n=1,ico

     if( (nsp1 .eq. no1(n) .and. mnsp1 .eq. mo1(n)) .and. (nsp2 .eq. no2(n) .and. mnsp2 .eq. mo2(n)) ) then
        icheck = 0
        exit
     endif

     if( (nsp2 .eq. no1(n) .and. mnsp2 .eq. mo1(n)) .and. (nsp1 .eq. no2(n) .and. mnsp1 .eq. mo2(n)) ) then
        icheck = 0
        exit
     endif

  enddo


  return
end subroutine testsameconfig
!$**************************************************************************************************
!$$   This subroutine is to calculate matrix elements of H for H_2 molecule for 4 one-electron functions
!$$   with given orbital angular momentum and its projection.
!$$   antisymmetric configurations

subroutine H12me(Znuc,Rd,is,bstF5,nsp1,nsp2,nsp1p,nsp2p,mnsp1,mnsp2,mnsp1p,mnsp2p,resultH,resultb)

  use sturmian_class

  implicit none

  real*8, intent(in):: Znuc, Rd
  integer, intent(in):: is
  type(basis_sturmian_nr), intent(in):: bstF5
  integer, intent(in):: nsp1,nsp2,nsp1p,nsp2p
  integer, intent(in):: mnsp1,mnsp2,mnsp1p,mnsp2p
  real*8, intent(out):: resultH,resultb

  type(sturmian_nr), pointer:: p1, p2, p1p, p2p
  real*8:: oneelme, result,  twoelme, tmp


  resultH = 0d0
  resultb = 0d0

  oneelme = 0d0
  twoelme = 0d0
  tmp = 0d0

  p1 => bstF5%b(nsp1)
  p2 => bstF5%b(nsp2)
  p1p => bstF5%b(nsp1p)
  p2p => bstF5%b(nsp2p)


!!$   overlap and one-electron operator matrix
  if(mnsp1 .eq. mnsp1p .and. mnsp2 .eq. mnsp2p) then
     resultb =  bstF5%ortint(nsp1,nsp1p) * bstF5%ortint(nsp2,nsp2p)

     call Hlagorb(bstF5,nsp1,nsp1p,mnsp1,result)
     oneelme =  result * bstF5%ortint(nsp2,nsp2p)
     call Hlagorb(bstF5,nsp2,nsp2p,mnsp2,result)
     oneelme =  oneelme + result * bstF5%ortint(nsp1,nsp1p)

  endif

  if((nsp1 .eq. nsp2 .and. mnsp1 .eq. mnsp2) .and. (nsp1p .ne. nsp2p .or. mnsp1p .ne. mnsp2p)) then
     ! <aa | ..| b1 b2>
     resultb = sqrt(2d0) * resultb
     oneelme = sqrt(2d0) * oneelme
  elseif((nsp1 .ne. nsp2 .or. mnsp1 .ne. mnsp2) .and. (nsp1p .eq. nsp2p .and. mnsp1p .eq. mnsp2p)) then
     ! <a1a2 | ..| bb>
     resultb = sqrt(2d0) * resultb
     oneelme = sqrt(2d0) * oneelme
  elseif((nsp1 .eq. nsp2 .and. mnsp1 .eq. mnsp2) .and. (nsp1p .eq. nsp2p .and. mnsp1p .eq. mnsp2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...
  else
     !<a1 a2 | ..| b1 b2>
     if(mnsp1 .eq. mnsp2p .and. mnsp2 .eq. mnsp1p) then
        resultb = resultb + (-1)**(is) * bstF5%ortint(nsp1,nsp2p) * bstF5%ortint(nsp2,nsp1p)

        call Hlagorb(bstF5,nsp1,nsp2p,mnsp1,result)
        oneelme =  oneelme + (-1)**(is) * result * bstF5%ortint(nsp2,nsp1p)

        call Hlagorb(bstF5,nsp2,nsp1p,mnsp2,result)
        oneelme =  oneelme + (-1)**(is) * result * bstF5%ortint(nsp1,nsp2p)
     endif
  endif


!!$   two-electron operator  matrix
  call V12me(p1,p2,p1p,p2p,mnsp1,mnsp2,mnsp1p,mnsp2p,result)
  twoelme  = result
  !                 print*, 'result', result
  if((nsp1 .eq. nsp2 .and. mnsp1 .eq. mnsp2) .and. (nsp1p .ne. nsp2p .or. mnsp1p .ne. mnsp2p)) then
     ! <aa | ..| b1 b2>
     twoelme  =  sqrt(2d0) * twoelme
  elseif((nsp1 .ne. nsp2 .or. mnsp1 .ne. mnsp2) .and. (nsp1p .eq. nsp2p .and. mnsp1p .eq. mnsp2p)) then
     ! <a1a2 | ..| bb>
     twoelme  =  sqrt(2d0) * twoelme
  elseif((nsp1 .eq. nsp2 .and. mnsp1 .eq. mnsp2) .and. (nsp1p .eq. nsp2p .and. mnsp1p .eq. mnsp2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...
  else
     !<a1 a2 | ..| b1 b2>
     call V12me(p1,p2,p2p,p1p,mnsp1,mnsp2,mnsp2p,mnsp1p,result)
     twoelme  = twoelme + (-1)**(is) * result
  endif

!!$  account for ZZ/R term
  tmp = 0d0
  if( Rd .ne. 0) then
     tmp = Znuc*Znuc/Rd  * resultb
  endif

  resultH = oneelme + tmp + twoelme


  return
end subroutine H12me
!$**************************************************************************************************
!
!$$  two-electron configurations are made from one-electron target states,
!$$  these states are not orthogonal
!$$  the idea here is to use H2+ 1s orbital together with Laguerre expansion (large exp. fall-off)
!$$  to capture electron-electron corellations for low lying states.

!$$   This subroutine is to calculate matrix elements of H for H_2 molecule for 4 one-electron target states
!$$   with given projection of orbital angular momentum and parity
!$$   antisymmetric configurations

subroutine H12me_st_notortog(Znuc,Rd,is,TargetStates1el,Nmax1el,e1me,ovlpst,bst,nst1,nst2,nst1p,nst2p,mnst1,mnst2,mnst1p,mnst2p,resultH,resultb)

  use sturmian_class
  use state_class

  implicit none

  real*8, intent(in):: Znuc, Rd
  integer, intent(in):: is
  type(basis_state), intent(in):: TargetStates1el
  integer, intent(in):: Nmax1el
  real*8, dimension(Nmax1el,Nmax1el), intent(in):: e1me, ovlpst
  type(basis_sturmian_nr), intent(in):: bst
  integer, intent(in):: nst1,nst2,nst1p,nst2p
  integer, intent(in):: mnst1,mnst2,mnst1p,mnst2p
  real*8, intent(out):: resultH,resultb

  type(sturmian_nr), pointer:: p1, p2, p1p, p2p
  integer:: nsp1,nsp2,nsp1p,nsp2p, i1,i2,i1p,i2p
  real*8:: tmp, tmpCI1, tmpCI1p, tmpCI2, tmpCI2p, ttt, result, oneelME, twoelME
  integer:: i1max, i2max, i1pmax, i2pmax

  real*8:: a1,a2, a3,a4, a5, a6

  resultH = 0d0
  resultb = 0d0

  oneelme = 0d0
  twoelme = 0d0
  tmp = 0d0


  i1max = get_nam(TargetStates1el%b(nst1))
  i2max = get_nam(TargetStates1el%b(nst2))
  i1pmax = get_nam(TargetStates1el%b(nst1p))
  i2pmax = get_nam(TargetStates1el%b(nst2p))




!!$   overlap and one-electron operator matrix
  if(mnst1 .eq. mnst1p .and. mnst2 .eq. mnst2p) then

     resultb = ovlpst(nst1,nst1p) * ovlpst(nst2,nst2p)
!!$     resultb = ovlp_st(TargetStates1el%b(nst1),TargetStates1el%b(nst1p)) * ovlpst(TargetStates1el%b(nst2),TargetStates1el%b(nst2p))
!     a1 =  ovlp_st(TargetStates1el%b(nst1),TargetStates1el%b(nst1p))
!     a2 = ovlp_st(TargetStates1el%b(nst2),TargetStates1el%b(nst2p))
!     print*, a1,ovlpst(nst1,nst1p)
!     print*, a2,ovlpst(nst2,nst2p)

     oneelme =   e1me(nst1,nst1p) * ovlpst(nst2,nst2p)
!!$     oneelme =   H1el_st(TargetStates1el%b(nst1),TargetStates1el%b(nst1p)) * ovlp_st(TargetStates1el%b(nst2),TargetStates1el%b(nst2p))
1     a3 = H1el_st(TargetStates1el%b(nst1),TargetStates1el%b(nst1p))
!     print*, a3,e1me(nst1,nst1p)

     oneelme =  oneelme + e1me(nst2,nst2p) * ovlpst(nst1,nst1p)
!!$     oneelme =  oneelme +   H1el_st(TargetStates1el%b(nst2),TargetStates1el%b(nst2p)) * ovlp_st(TargetStates1el%b(nst1),TargetStates1el%b(nst1p))
!     a5 = H1el_st(TargetStates1el%b(nst2),TargetStates1el%b(nst2p))
!     print*, a5,e1me(nst2,nst2p)

!     read*

  endif

  if((nst1 .eq. nst2) .and. (nst1p .ne. nst2p)) then
     ! <aa | ..| b1 b2>
     resultb = sqrt(2d0) * resultb
     oneelme = sqrt(2d0) * oneelme
  elseif((nst1 .ne. nst2) .and. (nst1p .eq. nst2p)) then
     ! <a1a2 | ..| bb>
     resultb = sqrt(2d0) * resultb
     oneelme = sqrt(2d0) * oneelme
  elseif((nst1 .eq. nst2) .and. (nst1p .eq. nst2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...
  else
     !<a1 a2 | ..| b1 b2>
     if(mnst1 .eq. mnst2p .and. mnst2 .eq. mnst1p) then

        resultb = resultb + (-1)**(is) * ovlpst(nst1,nst2p) * ovlpst(nst2,nst1p)
!!$        resultb = resultb + (-1)**(is) * ovlp_st(TargetStates1el%b(nst1),TargetStates1el%b(nst2p)) * ovlp_st(TargetStates1el%b(nst2),TargetStates1el%b(nst1p))

        oneelme =  oneelme + (-1)**(is) * e1me(nst1,nst2p) * ovlpst(nst2,nst1p) ! ovlp_st(TargetStates1el%b(nst2),TargetStates1el%b(nst1p)) *  H1el_st(TargetStates1el%b(nst1),TargetStates1el%b(nst2p))

        oneelme =  oneelme + (-1)**(is) * e1me(nst2,nst1p) *  ovlpst(nst1,nst2p) ! ovlp_st(TargetStates1el%b(nst1),TargetStates1el%b(nst2p)) *  H1el_st(TargetStates1el%b(nst2),TargetStates1el%b(nst1p))

     endif
  endif


!!$   two-electron operator  matrix
  ttt = 0d0
  do i1=1,i1max
     nsp1 = get_na(TargetStates1el%b(nst1),i1,1)
     tmpCI1 = get_CI(TargetStates1el%b(nst1),i1)
     p1 => bst%b(nsp1)
     do i2=1,i2max
        nsp2 = get_na(TargetStates1el%b(nst2),i2,1)
        tmpCI2 = get_CI(TargetStates1el%b(nst2),i2)
        p2 => bst%b(nsp2)

        do i1p=1,i1pmax
           nsp1p = get_na(TargetStates1el%b(nst1p),i1p,1)
           tmpCI1p = get_CI(TargetStates1el%b(nst1p),i1p)
           p1p => bst%b(nsp1p)
           do i2p=1,i2pmax
              nsp2p = get_na(TargetStates1el%b(nst2p),i2p,1)
              tmpCI2p = get_CI(TargetStates1el%b(nst2p),i2p)
              p2p => bst%b(nsp2p)

              tmp = tmpCI1 * tmpCI1p * tmpCI2 * tmpCI2p

              call V12me(p1,p2,p1p,p2p,mnst1,mnst2,mnst1p,mnst2p,result)

              ttt = ttt +  tmp * result

           enddo
        enddo
     enddo
  enddo

  twoelME = ttt



  !                 print*, 'result', result
  if((nst1 .eq. nst2) .and. (nst1p .ne. nst2p)) then
     ! <aa | ..| b1 b2>
     twoelme  =  sqrt(2d0) * twoelme
  elseif((nst1 .ne. nst2) .and. (nst1p .eq. nst2p)) then
     ! <a1a2 | ..| bb>
     twoelme  =  sqrt(2d0) * twoelme
  elseif((nst1 .eq. nst2) .and. (nst1p .eq. nst2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...
  else
     !<a1 a2 | ..| b1 b2>
     ttt = 0d0
     do i1=1,i1max
        nsp1 = get_na(TargetStates1el%b(nst1),i1,1)
        tmpCI1 = get_CI(TargetStates1el%b(nst1),i1)
        p1 => bst%b(nsp1)
        do i2=1,i2max
           nsp2 = get_na(TargetStates1el%b(nst2),i2,1)
           tmpCI2 = get_CI(TargetStates1el%b(nst2),i2)
           p2 => bst%b(nsp2)

           do i1p=1,i1pmax
              nsp1p = get_na(TargetStates1el%b(nst1p),i1p,1)
              tmpCI1p = get_CI(TargetStates1el%b(nst1p),i1p)
              p1p => bst%b(nsp1p)
              do i2p=1,i2pmax
                 nsp2p = get_na(TargetStates1el%b(nst2p),i2p,1)
                 tmpCI2p = get_CI(TargetStates1el%b(nst2p),i2p)
                 p2p => bst%b(nsp2p)

                 tmp = tmpCI1 * tmpCI1p * tmpCI2 * tmpCI2p

                 call V12me(p1,p2,p2p,p1p,mnst1,mnst2,mnst2p,mnst1p,result)

                 ttt = ttt + tmp * result

              enddo
           enddo
        enddo
     enddo
     twoelme  = twoelme + (-1)**(is) * ttt



  endif

!!$  account for ZZ/R term
  tmp = 0d0
  if( Rd .ne. 0) then
     tmp = Znuc*Znuc/Rd  * resultb
  endif

  resultH = oneelme + tmp + twoelme


  return
end subroutine H12me_st_notortog

!$**************************************************************************************************
!$$   two-electron configurations are made from one-electron target states,
!$$   This subroutine is to calculate matrix elements of H for H_2 molecule for 4 one-electron target states
!$$   with given projection of orbital angular momentum and parity
!$$   antisymmetric configurations

subroutine H12me_st(Znuc,Rd,is,TargetStates1el,bst,nst1,nst2,nst1p,nst2p,mnst1,mnst2,mnst1p,mnst2p,resultH,resultb)

  use sturmian_class
  use state_class

  implicit none

  real*8, intent(in):: Znuc, Rd
  integer, intent(in):: is
  type(basis_state), intent(in):: TargetStates1el
  type(basis_sturmian_nr), intent(in):: bst
  integer, intent(in):: nst1,nst2,nst1p,nst2p
  integer, intent(in):: mnst1,mnst2,mnst1p,mnst2p
  real*8, intent(out):: resultH,resultb

  type(sturmian_nr), pointer:: p1, p2, p1p, p2p
  integer:: nsp1,nsp2,nsp1p,nsp2p, i1,i2,i1p,i2p
  real*8:: tmp, tmpCI1, tmpCI1p, tmpCI2, tmpCI2p, ttt, result, oneelME, twoelME
  integer:: i1max, i2max, i1pmax, i2pmax


  resultH = 0d0
  resultb = 0d0


  i1max = get_nam(TargetStates1el%b(nst1))
  i2max = get_nam(TargetStates1el%b(nst2))
  i1pmax = get_nam(TargetStates1el%b(nst1p))
  i2pmax = get_nam(TargetStates1el%b(nst2p))
!  print*, 'i1max=', i1max
!  print*, 'CI=', get_CI(TargetStates1el%b(nst1),1)
!  print*, 'nst1 =',  get_na(TargetStates1el%b(nst1),1,1)


!!$  deal with one-electron ME: diagonal for one-electron target states

  oneelME = 0d0

  if(nst1 .eq. nst1p .and. nst2 .eq. nst2p) then

     resultb = 1d0

     oneelME = get_energy_st(TargetStates1el%b(nst1)) + get_energy_st(TargetStates1el%b(nst2))

  endif

!!$   two-electron operator  matrix
  twoelME = 0d0

  ttt = 0d0
  do i1=1,i1max
     nsp1 = get_na(TargetStates1el%b(nst1),i1,1)
     tmpCI1 = get_CI(TargetStates1el%b(nst1),i1)
     p1 => bst%b(nsp1)
     do i2=1,i2max
        nsp2 = get_na(TargetStates1el%b(nst2),i2,1)
        tmpCI2 = get_CI(TargetStates1el%b(nst2),i2)
        p2 => bst%b(nsp2)

        do i1p=1,i1pmax
           nsp1p = get_na(TargetStates1el%b(nst1p),i1p,1)
           tmpCI1p = get_CI(TargetStates1el%b(nst1p),i1p)
           p1p => bst%b(nsp1p)
           do i2p=1,i2pmax
              nsp2p = get_na(TargetStates1el%b(nst2p),i2p,1)
              tmpCI2p = get_CI(TargetStates1el%b(nst2p),i2p)
              p2p => bst%b(nsp2p)

              tmp = tmpCI1 * tmpCI1p * tmpCI2 * tmpCI2p

              call V12me(p1,p2,p1p,p2p,mnst1,mnst2,mnst1p,mnst2p,result)

              ttt = ttt +  tmp * result

           enddo
        enddo
     enddo
  enddo

  twoelME = ttt

  !                 print*, 'result', result
  if((nst1 .eq. nst2) .and. (nst1p .ne. nst2p)) then
     ! <aa | ..| b1 b2>
     twoelme  =  sqrt(2d0) * twoelme
  elseif((nst1 .ne. nst2) .and. (nst1p .eq. nst2p)) then
     ! <a1a2 | ..| bb>
     twoelme  =  sqrt(2d0) * twoelme
  elseif((nst1 .eq. nst2) .and. (nst1p .eq. nst2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...
  else
     !<a1 a2 | ..| b1 b2>
     ttt = 0d0
     do i1=1,i1max
        nsp1 = get_na(TargetStates1el%b(nst1),i1,1)
        tmpCI1 = get_CI(TargetStates1el%b(nst1),i1)
        p1 => bst%b(nsp1)
        do i2=1,i2max
           nsp2 = get_na(TargetStates1el%b(nst2),i2,1)
           tmpCI2 = get_CI(TargetStates1el%b(nst2),i2)
           p2 => bst%b(nsp2)

           do i1p=1,i1pmax
              nsp1p = get_na(TargetStates1el%b(nst1p),i1p,1)
              tmpCI1p = get_CI(TargetStates1el%b(nst1p),i1p)
              p1p => bst%b(nsp1p)
              do i2p=1,i2pmax
                 nsp2p = get_na(TargetStates1el%b(nst2p),i2p,1)
                 tmpCI2p = get_CI(TargetStates1el%b(nst2p),i2p)
                 p2p => bst%b(nsp2p)

                 tmp = tmpCI1 * tmpCI1p * tmpCI2 * tmpCI2p

                 call V12me(p1,p2,p2p,p1p,mnst1,mnst2,mnst2p,mnst1p,result)

                 ttt = ttt + tmp * result

              enddo
           enddo
        enddo
     enddo
     twoelme  = twoelme + (-1)**(is) * ttt
  endif

!!$  account for ZZ/R term
  tmp = 0d0
  if( Rd .ne. 0) then
     tmp = Znuc*Znuc/Rd  * resultb
  endif

  resultH = oneelme + tmp +  twoelme

  return
end subroutine H12me_st

!$**************************************************************************************************
!
! This is matrix element of V(1,2) electron-electron potential for 4 functions:
! <n1 n2 | V(1,2) | n1p n2p >, where V(1,2) = sum_{lam} V_{lam}(1,2)
! Note: potential V(12) can not change the total ang.mom. projection M of a configuration
subroutine V12me(pn1,pn2,pn1p,pn2p,m1,m2,m1p,m2p,result)

  use sturmian_class
  use grid_radial

  implicit none

  type(sturmian_nr), intent(in):: pn1,pn2,pn1p,pn2p
  integer, intent(in):: m1, m1p, m2, m2p
  real*8, intent(out):: result

  real*8:: Yint


  integer:: l1, l2, l1p, l2p
  integer:: maxfm, minfm, i1, i2, minfun, maxfun, minfun1, maxfun1, lam, lammin, lammax
  real*8:: rlam, rq, reslam
  real*8:: tmp, tmp1, tmp2, sum1, sum2
  integer:: minf1,maxf1, minf1p,maxf1p,minf2,maxf2, minf2p,maxf2p
  real*8, pointer, dimension(:):: f1, f1p, f2, f2p
  real*8, dimension(grid%nr):: temp, fun, fun1, fun11
  real*8:: rl1,rl2,rl1p,rl2p,rm1, rm2,rm1p,rm2p


  result = 0d0

  l1 = get_ang_mom(pn1)
  l2 = get_ang_mom(pn2)
  l1p = get_ang_mom(pn1p)
  l2p = get_ang_mom(pn2p)
  rl1 = l1
  rl2 = l2
  rl1p = l1p
  rl2p = l2p

  rm1 = m1
  rm2 = m2
  rm1p = m1p
  rm2p = m2p

  f1 => fpointer(pn1)
  f1p => fpointer(pn1p)

  f2 => fpointer(pn2)
  f2p => fpointer(pn2p)

  minf1 = get_minf(pn1)
  maxf1 = get_maxf(pn1)
  minf1p = get_minf(pn1p)
  maxf1p = get_maxf(pn1p)

  minf2 = get_minf(pn2)
  maxf2 = get_maxf(pn2)
  minf2p = get_minf(pn2p)
  maxf2p = get_maxf(pn2p)


  maxfm = min(maxf1,maxf1p)
  minfm = max(minf1,minf1p)

  fun1(minfm:maxfm) = f1(minfm:maxfm)*f1p(minfm:maxfm)*grid%weight(minfm:maxfm)

  minfun = max(minf2,minf2p)
  maxfun = min(maxf2,maxf2p)

  fun(minfun:maxfun) = f2(minfun:maxfun)*f2p(minfun:maxfun) * grid%weight(minfun:maxfun)


  lammin=max(abs(l1-l1p),abs(l2-l2p))
  lammax=min(l1+l1p,l2+l2p)

  do lam=lammin,lammax

     rlam = lam
     rq = m1p - m1
     tmp1 = (-1)**(nint(rq))*Yint(rl1,rm1,rlam,-rq,rl1p,rm1p)
     tmp2 = Yint(rl2,rm2,rlam,rq,rl2p,rm2p)
     tmp = tmp1 * tmp2

!     print*, 'lam, q, tmp', lam, rq,tmp1,tmp2

     if( tmp .eq. 0d0) then
        cycle
     endif

     call form(lam,fun,minfun,maxfun,maxfm,temp,i1,i2)

     minfun1 = max(i1,minfm)
     maxfun1 = min(i2,maxfm)

     fun11(minfun1:maxfun1) = fun1(minfun1:maxfun1) * temp(minfun1:maxfun1)

     reslam = SUM(fun11(minfun1:maxfun1))

     result = result + reslam * tmp

  enddo


end subroutine V12me




!$**************************************************************************************************

subroutine convert_from_st_to_sp(bst,TargetStates1el,TargetStates2el,nspm,Nmax1el,Nmax)

  use sturmian_class
  use state_class
  use grid_radial

  implicit none

  type(basis_sturmian_nr), intent(in):: bst
  type(basis_state), intent(in):: TargetStates1el
  type(basis_state):: TargetStates2el
  integer, intent(in):: nspm,Nmax,Nmax1el

  integer:: nst, namst, i, ist, n1st, n2st, nam1, nam2, i1, i2, n1, n2, ncm
  real*8:: CIst, CI1, CI2
!  real*8, dimension(nspm,nspm):: CIno
!  real*8, dimension(nspm*nspm):: CItmp
  real*8, dimension(:,:), allocatable:: CIno
  real*8, dimension(:), allocatable:: CItmp
  integer, dimension(:), allocatable:: no1, no2, mo1, mo2, phase
  integer, dimension(nspm):: nspar
  integer:: newnspm, j, jn1, jn2
  integer, dimension(:), allocatable:: arnsp

print*, 'start: convert_from_st_to_sp()'
  print*, 'nspm,Nmax1el,Nmax:', nspm,Nmax1el,Nmax

!!$ this loop can be paralalized
  do nst=1,Nmax

!     print*, '>> nst =', nst

     namst = get_nam(TargetStates2el%b(nst))
!     print*, 'namst=', namst


!!$ as nspm can be very large we need to avoid dealing with all one-electron orbitals
!!$ and use only those orbitals that are used for given state i
!!$  First create two arrays nspar(1:nspm) and arnsp(1:newnspm)
!!$
!!$ create array nspar(1:nspm) where we record 1 if the orbital is used in the description of the given two-electron state ist

     nspar(:) = 0

     do ist=1,namst

        n1st = get_na(TargetStates2el%b(nst),ist,1)
        nam1 = get_nam(TargetStates1el%b(n1st))
        do i1=1,nam1
           n1 = get_na(TargetStates1el%b(n1st),i1)
           nspar(n1) = 1
        enddo

        n2st = get_na(TargetStates2el%b(nst),ist,2)
        nam2 = get_nam(TargetStates1el%b(n2st))
        do i2=1,nam2
           n2 = get_na(TargetStates1el%b(n2st),i2)
           nspar(n2) = 1
        enddo
     enddo
!!$ check how many orbitals are used for decription of the state ist
     newnspm = 0
     do i=1,nspm
        if(nspar(i) .ne. 0) then
           newnspm = newnspm + 1
        endif
     enddo
     allocate(arnsp(newnspm))
     arnsp(:) = 0
!     print*, 'newnspm=', newnspm

     j = 0
     do i=1,nspm
        if(nspar(i) .ne. 0) then
           j = j + 1
           arnsp(j) = i
           nspar(i) = j
        endif
     enddo
!!$ at this stage arrays nspar() and arnsp() allow to move forward and backwards between old enumeration (1 to nspm)
!!$ and new enumeration (1 to newnspm)
!!$ Hopefully newnspm << nspm

     allocate(CIno(newnspm,newnspm))
     allocate(CItmp(newnspm*newnspm),no1(newnspm*newnspm),no2(newnspm*newnspm),mo1(newnspm*newnspm),mo2(newnspm*newnspm),phase(newnspm*newnspm))
     CIno(:,:) = 0d0

     do ist=1,namst

        n1st = get_na(TargetStates2el%b(nst),ist,1)
        n2st = get_na(TargetStates2el%b(nst),ist,2)

        CIst = get_CI(TargetStates2el%b(nst),ist)

        nam1 = get_nam(TargetStates1el%b(n1st))
        nam2 = get_nam(TargetStates1el%b(n2st))
!        print*, 'n1st, n2st=', n1st, n2st
        do i1=1,nam1

           n1 = get_na(TargetStates1el%b(n1st),i1)
           CI1 = get_CI(TargetStates1el%b(n1st),i1)
!           mo1(n1) = get_angmom_proj(TargetStates1el%b(n1st))

           jn1 = nspar(n1)

           do i2=1,nam2

              n2 = get_na(TargetStates1el%b(n2st),i2)
              CI2 = get_CI(TargetStates1el%b(n2st),i2)
!              mo2(n2) = get_angmom_proj(TargetStates1el%b(n2st))
              jn2 = nspar(n2)

              CIno(jn1,jn2) = CIno(jn1,jn2) + CIst * CI1 * CI2
!              print*, n1, n2, jn1, jn2,  CIno(jn1,jn2)
           enddo

        enddo

     enddo


     i = 0
     do jn1=1,newnspm
        do jn2=jn1,newnspm

           if(CIno(jn1,jn2) .eq. 0 ) cycle

           i = i + 1

           n1 = arnsp(jn1)
           no1(i) = n1
           n2 = arnsp(jn2)
           no2(i) = n2
           mo1(i) = get_ang_mom_proj(bst%b(n1))
           mo2(i) = get_ang_mom_proj(bst%b(n2))
           if(n1 .eq. n2)  then
              CItmp(i) = CIno(jn1,jn2)
              phase(i) = 1d0
           else
              CItmp(i) = CIno(jn1,jn2) * sqrt(2d0)
              phase(i) = nint(CIno(jn1,jn2)/CIno(jn2,jn1))
           endif
!           print*, i, n1, n2, CItmp(i), phase(i)
        enddo
     enddo

     ncm = i

     call setupCI(TargetStates2el%b(nst),ncm,CItmp(1:ncm),no1(1:ncm),mo1(1:ncm),no2(1:ncm),mo2(1:ncm),phase(1:ncm))
     deallocate(arnsp,CIno, CItmp,no1,no2,mo1,mo2,phase)

  enddo

  print*, 'finish convert_from_st_to_sp()'

end subroutine convert_from_st_to_sp
