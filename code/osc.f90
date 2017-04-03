subroutine oscstr1
  use input_data 
  use one_electron_func
  use sturmian_class
  use state_class
  use  target_states
  use data_targets

  implicit none  

  integer:: n, np, Nst, ido_osc ! Loops and input option controls
  real*8::  result_l, result_v,  result_p, result_pc ! Results
  real*8:: TranEnergy, ioniz_en, temp_p
  integer:: map, ma, dM

  open(301,file='osc.strength')
  open(302,file='dip.mom')
  
  ioniz_en =  real(27.2114 * get_energy_st(TargetStates%b(1)) )

  write(301,'("OSCILLATOR STRENGTHS  " )')
  write(301,'("Exact solutions for the Oscillator Strengths are for the R = 2.0 case" )')
  write(301,'("Ref[1]: , Ref[2]: , Ref[3]: " )')

  write(301,'("Initial (np) to Final States (n):" )') 
  write(301,'(3X,"np   Mp   parp   n    M    par    Excitation Energy(eV)  Length       Velocity       Exact       Ref" )')
  
  write(302,'("STATIC DIPOLE POLARIZABILITY") ')
  write(302,'("Exact solution for the Polarizability are for the R = 2.0 case." )')
  write(302,'("Polarizability calculated using oscillator strength length guages." )')
  write(302,'("Ref[1]: , Ref[2]: , Ref[3]: " )')

  write(302,'(3X,"np   Mp   parp   n    M    par    Excitation Energy(eV)    Polarizability(a_0^3)" )')

  Nst = basis_size_st(TargetStates)

  do np = 1, 1 ! Nst
     map = TargetStates%b(np)%M

     result_p = 0d0	
     result_pc = 0d0

     do n = np+1, Nst
        if(TargetStates%b(np)%parity * TargetStates%b(n)%parity .ne. -1) cycle
        ma = TargetStates%b(n)%M ! Get magnetic sublevels? M = m for one-electron
        dM = map - ma 
        if (abs(dM) .ge. 2) cycle

        TranEnergy =  get_energy_st(TargetStates%b(n)) -  get_energy_st(TargetStates%b(np))  ! In units of a.u., hence multiply in osc by 2 for Ryd.
        
        if( data_in%calculation_type.eq.0 .or. data_in%calculation_type.eq.1 ) then 
           call oscstr1_st(data_in%iosc,n,np,TargetStates%b(n),TargetStates%b(np),TranEnergy,result_l,result_v)           
        elseif( data_in%calculation_type .eq. 2 ) then
           ! Call Spheroidal Oscillator Strength Routine
        endif
        
        temp_p = 0d0
        if ( ABS(dM) == 1) then
           temp_p =  result_l / ( 2.0 * TranEnergy * TranEnergy)  ! Summing over f_n dm = +/-1, therefore divide by 2
        else if ( dM == 0 ) then
           temp_p =  result_l / ( TranEnergy * TranEnergy)
        end if
        if (  get_energy_st(TargetStates%b(n)) > 0 ) then ! For continuum contribution
           result_pc = result_pc + temp_p
        end if
        result_p = result_p + temp_p

        write(301,'(I5,F6.1,I5,I5,F6.1,I5,5X,F10.5,8X,3(2X,E11.4,1X),1X,I5)') np, real(map), TargetStates%b(np)%parity,n, real(ma), TargetStates%b(n)%parity, 27.2114 * TranEnergy,result_l,result_v, H2Iosc(np,n), Ref_H2Iosc(np,n)          

        ! file: dip.mom
        write(302,'(I5,F6.1,I5,I5,F6.1,I5,5X,F10.5,15X,E10.3)') np, real(map), TargetStates%b(np)%parity,n, real(ma), TargetStates%b(n)%parity, 27.2114 * TranEnergy, temp_p

     end do

     write(302,'(3X,"State:",I5,"   Discrete Spectrum Polarizability:",E11.4,"   Continuum Polarizability:",E11.4 )') np, result_p, result_pc
  end do
  
  close(301)
  close(302)

end subroutine oscstr1!
!
subroutine  oscstr1_st(iosc,n,np,pn,pnp,TranEnergy,result_l,result_v)
  use grid_radial 
  use sturmian_class
  use state_class
  use  target_states
  use one_electron_func


  implicit none  
  integer, intent(in):: iosc
  integer, intent(in):: n, np     ! Molecular states
  type(state), intent(in):: pn, pnp   ! One electron states
  real*8, intent(in):: TranEnergy
  real*8, intent(out):: result_l, result_v

  integer::  Nst, Nstp, l, lp, ma, map, ne, nep,  dL, dM  !quantum numbers and loops
  type(sturmian_nr), pointer:: nn, nnp  !  one-electron orbitals
  integer:: i, ipp
  integer::  minf, maxf, minfp, maxfp,  i1, i2 ! Radial integrations
  real*8:: Coeff, nnpRnn, temp, CI, CIp, temp_l
  real*8:: temp_v, ideriv,  vint
  real*8, pointer, dimension(:)::  f, fp ! f One-electron functions
  real*8, pointer, dimension(:):: weight, gridr

 
  weight => grid%weight
  gridr => grid%gridr
  result_l = 0d0
  temp_l = 0d0
  result_v = 0d0
  temp_v = 0d0

  Nstp = get_nam(pnp) 
  map = TargetStates%b(np)%M

  Nst = get_nam(pn)  ! ncm = nam which is the size of the array na(:), which is also the number of CI coeff
  ma = TargetStates%b(n)%M ! Get magnetic sublevels? M = m for one-electron
  
  dM = map - ma 

  if ( ABS(dM) <= 1 ) then ! Selections rules dm = +/- 1, 0

     do nep = 1, Nstp  ! Looping over the number of one-electron orbitals which make up each molecule state
        
        CIp = get_CI(pnp, nep)
        ipp = get_na(pnp,nep,1)
        nnp => bst%b(ipp)
        fp => fpointer(nnp)
        lp = get_ang_mom(nnp)
        
        do ne =  1, Nst
           
           CI = get_CI(pn, ne) ! Get CI coefficient for that one-electron function
           i = get_na(pn,ne,1) ! which is na(ne), where nc = 1,2,...Nst
           nn => bst%b(i)   ! Makes nn the one-electron states?
           f => fpointer(nn)   ! One electron functions?
           l = get_ang_mom(nn) ! Gets angular momentum from the one-electron orbitals
           
           dL = lp - l

           if ( dL == 1 .OR. dL == -1  ) then ! Selection rules dl = +/- 1 
              temp = (2 * l + 2 + dL) * (2 * l + dL)
              
              ! Extent of radial grid for integrations
              minf = get_minf(nn)
              maxf = get_maxf(nn)
              minfp = get_minf(nnp)
              maxfp = get_maxf(nnp)
              
              i1 = max(minf,minfp)
              i2 = min(maxf,maxfp)
              
              ! Radial Integration <fp|r|f>
              nnpRnn = SUM( fp(i1:i2) *  f(i1:i2)  * gridr(i1:i2) * weight(i1:i2) ) ! One-electron function normalized so that r^2 term is taken out of integrations.


              !Velocity Integrals
              ! Derivativie Integral Note: d/dr(phi(r)/r) = phi'(r)/r -  phi(r)/(r^2)
              ideriv =  SUM((( f(i1+1:i2) - f(i1:i2-1)) / (gridr(i1+1:i2) - gridr(i1:i2-1))) * fp(i1:i2-1) * weight(i1:i2-1))

              ! Varsholovich pg 147 and look at derivative operator second term.
              vint = SUM( fp(i1:i2) * f(i1:i2) * weight(i1:i2) / gridr(i1:i2) )
              if ( l < lp ) then
                 vint = ideriv - (l + 1.0) * vint
              else if ( l > lp) then
                 vint = ideriv + dble(l) * vint
              end if

              Coeff = 0d0 ! Length Guage Block Varsholovich pg 145 
              if ( dM == 0) then ! Parallel Transitions <z> 
                 Coeff = sqrt( ( dble(l) + 0.5 + dble(dL) * 0.5 ) ** 2.0 - dble(ma * ma)) 
              end if
              if ( ABS(dM) == 1) then
                 if (ma < map  .AND. l < lp   ) then  ! Perpendicular Transitions <x +/- iy>              
                    Coeff = sqrt(dble((l + ma + 2 ) * (l + ma + 1)))
                 else if (  ma > map  .AND. l < lp ) then
                    Coeff = - sqrt(dble((l - ma + 2) * (l - ma + 1)))  
                 else if (  ma < map  .AND. l > lp  ) then         
                    Coeff = - sqrt(dble((l - ma ) * (l - ma - 1)))  
                 else if (  ma > map  .AND. l > lp ) then
                    Coeff = sqrt(dble((l + ma ) * (l + ma - 1)))               
                 end if
                 Coeff = Coeff /  sqrt(2.0)                 
              end if

              Coeff = Coeff / sqrt(temp) 
              
              result_l = temp_l +  CIp * CI * Coeff * nnpRnn 
              result_v = temp_v +  CIp * CI * Coeff * vint
              temp_l = result_l
              temp_v = result_v

           end if ! End dL
        
        end do ! ne RHS one-electron functions
     
     end do    ! nep LHS one-electron functions

     if ( ABS(dM) == 1) then
        result_l = 4.0 * ABS(TranEnergy) * result_l * result_l / 3.0 ! Added orbital degeneracy g = 2, should test for ma=+/-2 etc 
        result_v = 4.0 * result_v * result_v / (3.0 * ABS(TranEnergy))
     else 
        result_l = 2.0 * ABS(TranEnergy) * result_l * result_l / 3.0 ! Need to add orbital degeneracy, prob best to do it in the coefficients
        result_v = 2.0 * result_v * result_v / ( 3.0  * ABS(TranEnergy) )
     end if

  end if ! End the selection rule ABS(dM) <= 1 
   
end subroutine oscstr1_st

