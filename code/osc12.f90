
subroutine oscstr_2e()

  use input_data 
  use data_targets
  use one_electron_func
  use sturmian_class
  use target_states

  implicit none  
  integer:: nstate, npstate, Nmax, ninitial_max         ! Loops and input option controls
  ! Results
  real*8::  result_l, result_v,  result_p, result_pc   
  real*8:: result_p_par, result_p_perp
  real*8:: TranEnergy,  temp_p, Energy, Osc 
  integer:: map, ma, dM, Spin, Spinp, parity, parityp   ! Quantum Numbers Selection Rules

  Nmax = basis_size_st(TargetStates2el)                 ! Number of Molecular States Generated

  ninitial_max = min(data_in%iosc,Nmax)

  print*,"******************************************"
  print*,"OSCILLATOR STRENGTHS AND POLARIZABILITY"
  print*,"******************************************"

  ! INITIAL Molecular State Loop
  do nstate = 1, ninitial_max                           
     ! Initial State number n
     ma = NINT(TargetStates2el%b(nstate)%M)       ! 2e Molecular State Angular Projection
     Spin = NINT(TargetStates2el%b(nstate)%spin)  ! 2e Molecular State Spin
     parity = TargetStates2el%b(nstate)%parity    ! 2e Molecular State Parity

     result_p = 0d0
     result_pc = 0d0
     result_p_par = 0d0
     result_p_perp = 0d0

     ! FINAL Molecular State Loop     
     do npstate = nstate + 1, Nmax   
        ! Final State number np
        map = NINT(TargetStates2el%b(npstate)%M)          
        Spinp = NINT(TargetStates2el%b(npstate)%spin)
        parityp = TargetStates2el%b(npstate)%parity
        
        dM = map - ma        
        
        TranEnergy = get_energy_st(TargetStates2el%b(npstate)) -  get_energy_st(TargetStates2el%b(nstate))

        result_l= 0d0
        result_v = 0d0
        ! MAKE SURE BOTH MOLECULAR STATES HAVE THE SAME SPIN and difference Angular Projections <= 1
        ! Parity must change
        if ( Spin /= Spinp ) cycle
        if (ABS(dM) > 1 ) cycle
        if ( parity == parityp) cycle

!        call  oscstrength_2e(nstate,npstate,TranEnergy,result_l,result_v)
! subroutine _config should be more efficient 
        call   oscstrength_2e_config(nstate,npstate,TranEnergy,result_l,result_v)

        print*,"ni",nstate,"  nf",npstate,"  dM",dM
        print*,"L",result_l,"V",result_v

        ! Polarizabilty calulated from Oscillator strengths
        temp_p = 0d0
        if ( ABS(dM) == 1) then
           ! Summing over f_n dm = +/-1, therefore divide by 2
           temp_p =  result_l / ( 2.0 * TranEnergy * TranEnergy)
           result_p_perp =  result_p_perp + 1.5 * temp_p
        else if ( dM == 0 ) then
           temp_p =  result_l / ( TranEnergy * TranEnergy )
           result_p_par =  result_p_par + 3.0 * temp_p
        end if
        ! For continuum contribution
        if (  get_energy_st(TargetStates2el%b(npstate)) > 0 ) then
           result_pc = result_pc + temp_p
        end if
        result_p = result_p + temp_p            
        


     end do
     if ( nstate == 1) then
        print*
        ! Oscillator Strengths: R-matrix S. Branchett 92
        ! Polarizability: W. Kolos 67  TABLE II R = 1.4
        print*,"Polarizability: Ref[1]:  W. Kolos 67 R = 1.4"

        print*,"Parallel:",result_p_par,"Ref[1]: 6.38049"
        print*,"Perpendicular:",result_p_perp,"Ref[1]: 4.57769"
        ! a_T = 1/3 * a_par + 2/3 * a_perp
        print*,"Total:",result_p,"Ref[1]: 5.17862"
     end if
  end do
  print*

end subroutine oscstr_2e




subroutine  oscstrength_2e(nstate,npstate,TranEnergy,result_l,result_v)
! Calculates oscillator strengths in the form of Length <np1,np2|L|n1,n2>
! f_L =  sum C_(np1,np2) C(n1,n2) <np1,np2|L|n1,n2> 

  use grid_radial 
  use sturmian_class
  use target_states
  use one_electron_func
  use state_class


  implicit none  
  integer, intent(in):: nstate, npstate               ! Molecular state index
  real*8, intent(in):: TranEnergy
  real*8, intent(out):: result_l, result_v
  ! Quantum numbers 
  integer::  la1, la2, lap1, lap2, ma, map, mp1, mp2, m1, m2,  dL, dM, Spin, Spinp
  integer:: parity, parityp
  integer::  ne_con,  nep_con, Ncon, Npcon                               ! Configuration Loops
  type(sturmian_nr), pointer:: tn1, tn2, tnp1, tnp2                      ! One-electron orbitals
  integer:: ind1, ind2, indp1, indp2
  integer::  minf, maxf, minfp, maxfp,  ir1, ir2                         ! Radial integrations
  real*8:: Coeff, nnpRnn, temp, CI, CIp, temp_l
  real*8:: temp_v, ideriv1,  vint1, r_overlap1, overlap2                 ! Length velocity Calcs
  real*8, pointer, dimension(:)::  f1, f2, fp1, fp2                      ! f One-electron functions
  real*8, pointer, dimension(:):: weight, gridr

 
  weight => grid%weight
  gridr => grid%gridr
  result_l = 0d0
  temp_l = 0d0
  result_v = 0d0
  temp_v = 0d0


! Final State number npstate
  Npcon =  TargetStates2el%b(npstate)%nam        ! Number of A.O. and Molecular ion configurations.        
  map = NINT(TargetStates2el%b(npstate)%M )      ! Molecular state Angular Projection
  Spinp = NINT(TargetStates2el%b(npstate)%spin)  ! 2e Molecular State Spin
  parityp = TargetStates2el%b(npstate)%parity    ! 2e Molecular State Parity

! Initial State number nstate
  Ncon = TargetStates2el%b(nstate)%nam         
  ma = NINT(TargetStates2el%b(nstate)%M)       
  Spin = NINT(TargetStates2el%b(nstate)%spin)
  parity = TargetStates2el%b(nstate)%parity     
  
  dM = map - ma 

! MAKE SURE BOTH MOLECULAR STATES HAVE THE SAME SPIN and difference Angular Projections <= 1
  ! Parity must change
  if ( Spin /= Spinp ) return
  if (ABS(dM) > 1 ) return
  if ( parity == parityp) return


! Looping over FINAL Molecular State orbitals.
  do nep_con =  1, Npcon

! Quantum numbers and functions for ALPHA
     indp1 = TargetStates2el%b(npstate)%na(nep_con)       ! Final state number np. nep A.O.
     tnp1 => bst%b(indp1)                              !           
     fp1 => fpointer(tnp1)                                ! One electron functions
     lap1 = get_ang_mom(tnp1)                             ! Gets angular momentum from the one-electron orbitals
     mp1 =  TargetStates2el%b(npstate)%ma(nep_con)        ! Get angular projection of A.O.   

! Quantum numbers and functions for BETA
     indp2 = TargetStates2el%b(npstate)%nb(nep_con)       
     tnp2 => bst%b(indp2)                                         
     fp2 => fpointer(tnp2)                                
     lap2 = get_ang_mom(tnp2)                             
     mp2 =  TargetStates2el%b(npstate)%mb(nep_con)           
     CIp = get_CI(TargetStates2el%b(npstate),nep_con)  

! looping over INITIAL Molecular State orbitals.
     do ne_con = 1, Ncon

! Quantum numbers and functions for GAMMA
        ind1 = TargetStates2el%b(nstate)%na(ne_con)       
        tn1 => bst%b(ind1)                                          
        f1 => fpointer(tn1)                               
        la1 = get_ang_mom(tn1)                            
        m1 =  TargetStates2el%b(nstate)%ma(ne_con)           


! Quantum numbers and functions for DELTA
        ind2 = TargetStates2el%b(nstate)%nb(ne_con)       
        tn2 => bst%b(ind2)                                         
        f2 => fpointer(tn2)                               
        la2 = get_ang_mom(tn2)                            
        m2 =  TargetStates2el%b(nstate)%mb(ne_con)           
        CI = get_CI(TargetStates2el%b(nstate),ne_con)  
        
        dL = lap1 - la1
        dM = mp1 - m1
        
        ! Selections Rules 
        if ( ABS(dL) > 1 .OR. ABS(dM) > 1 ) cycle
       
        if ( lap2 /= la2 .OR. mp2 /= m2  ) cycle
        ! DO OVERLAPOF COORDINATE SPACE 2  
        minf = get_minf(tn2)
        maxf = get_maxf(tn2)
        minfp = get_minf(tnp2)
        maxfp = get_maxf(tnp2)      
        ir1 = max(minf,minfp)
        ir2 = min(maxf,maxfp)  
        overlap2 = SUM( fp2(ir1:ir2) *  f2(ir1:ir2) * weight(ir1:ir2) )

        ! COORDINATE 1 RADIAL INTEGRALS
        minf = get_minf(tn1)
        maxf = get_maxf(tn1)
        minfp = get_minf(tnp1)
        maxfp = get_maxf(tnp1)
        ir1 = max(minf,minfp)
        ir2 = min(maxf,maxfp) 
        ! Radial Integration <fp1|r|f1>
         r_overlap1 = SUM( fp1(ir1:ir2) *  f1(ir1:ir2)  * gridr(ir1:ir2) * weight(ir1:ir2) ) 
        
        !Velocity Integrals
        ! Derivatiie Integral Note: d/dr(phi(r)/r) = phi'(r)/r -  phi(r)/(r^2)
        ideriv1 =  SUM((( f1(ir1 + 1 : ir2 ) - f1(ir1: ir2 - 1)) /(gridr(ir1 + 1 :ir2) - gridr(ir1: ir2 - 1))) * fp1(ir1: ir2 - 1) * weight(ir1:ir2 - 1))
        vint1 = SUM( fp1(ir1:ir2) * f1(ir1:ir2) * weight(ir1:ir2) / gridr(ir1:ir2) )


        if ( la1 < lap1 ) then
           vint1 = ideriv1 - (dble(la1) + 1.0) * vint1
        else if ( la1 > lap1) then
           vint1 = ideriv1 + dble(la1) * vint1
        end if
                    

        temp = (2.0 * la1 + 2.0 + dL) * (2.0 * la1 + dL)
        Coeff = 0d0 ! Length Guage Block Varsholovich pg 145 
        if ( dM == 0) then ! Parallel Transitions <z> 
           Coeff = sqrt( ( dble(la1) + 0.5 + dble(dL) * 0.5 ) ** 2.0 - dble(m1 * m1)) 
        end if
        if ( ABS(dM) == 1) then
           if ( la1 < lap1   ) then  ! Perpendicular Transitions <x +/- iy> 
              Coeff = dble(dM)*sqrt(dble(la1+dM*m1+2)*dble(la1+dM*m1+1))
           else if ( la1 > lap1) then
              Coeff = -dble(dM)*sqrt(dble(la1-dM*m1)*dble(la1-dM*m1-1))
           end if
           Coeff = Coeff /  sqrt(2.0)                 
        end if
              
        Coeff = Coeff / sqrt(temp) 
              
        ! Multiply by 2 for 2e target. Operator L: z_1 + z_2 = 2*z_1
        result_l = temp_l + 2.0 * CIp * CI * Coeff * r_overlap1 * overlap2
        result_v = temp_v + 2.0 * CIp * CI * Coeff * vint1 * overlap2
        temp_l = result_l
        temp_v = result_v


     end do ! Initial States Orbitals
     
  end do    ! Final States Orbitals

  dM = map - ma 
! TranEnergy is in a.u. converted to Rydbergs here.
  if ( ABS(dM) == 1) then
     result_l = 4.0 * ABS(TranEnergy) * result_l * result_l / 3.0 ! Added orbital degeneracy g = 2, should test for ma=+/-2 etc 
     result_v = 4.0 * result_v * result_v / (3.0 * ABS(TranEnergy))
  else 
     result_l = 2.0 * ABS(TranEnergy) * result_l * result_l / 3.0 
     result_v = 2.0 * result_v * result_v / ( 3.0  * ABS(TranEnergy) )
  end if


 
end subroutine oscstrength_2e







subroutine  oscstrength_2e_config(nstate,npstate,TranEnergy,result_l,result_v)
! Calculates oscillator strengths in the form of Length <np1,np2|L|n1,n2>
! f_L =  sum C_(np1,np2) C(n1,n2) <np1,np2|L|n1,n2> 

  use grid_radial 
  use sturmian_class
  use target_states
  use one_electron_func
  use state_class


  implicit none  
  integer, intent(in):: nstate, npstate               ! Molecular state index
  real*8, intent(in):: TranEnergy
  real*8, intent(out):: result_l, result_v
  ! Quantum numbers 
  integer::  la1, la2, lap1, lap2, ma, map, mp1, mp2, m1, m2,  dL, dM, Spin, Spinp
  integer:: m1_temp, mp1_temp, parity, parityp
  ! Configuration Loops
  integer::  ne_con,  nep_con, Ncon, Npcon,  n_orb, np_orb, n_orb_max, np_orb_max 
  integer:: nuse_norb, nuse_nporb
  type(sturmian_nr), pointer:: tn1, tn2, tnp1, tnp2                      ! One-electron orbitals
  integer:: ind1, ind2, indp1, indp2, indp1_2, ind1_2
  integer::  minf, maxf, minfp, maxfp,  ir1, ir2                         ! Radial integrations
  real*8:: Coeff, nnpRnn, temp, CI, CIp, temp_l
  real*8:: temp_v, ideriv1,  vint1, r_overlap1, temp_overlap2, overlap2  ! Length velocity Calcs
  real*8, pointer, dimension(:)::  f1, f2, fp1, fp2                      ! f One-electron functions
  real*8, pointer, dimension(:):: weight, gridr

 
  weight => grid%weight
  gridr => grid%gridr
  result_l = 0d0
  temp_l = 0d0
  result_v = 0d0
  temp_v = 0d0


  ! Final State number npstate
  Npcon =  TargetStates2el%b(npstate)%nam          ! Number of A.O. and Molecular ion configurations.        
  map = NINT(TargetStates2el%b(npstate)%M )        ! Molecular state Angular Projection
  Spinp = NINT(TargetStates2el%b(npstate)%spin)    ! 2e Molecular State Spin
  parityp = TargetStates2el%b(npstate)%parity      ! 2e Molecular State Parity
  np_orb_max = TargetStates2el%b(npstate)%nusemax  ! Number of orbitals used to describe this state

  ! Initial State number nstate
  Ncon = TargetStates2el%b(nstate)%nam         
  ma = NINT(TargetStates2el%b(nstate)%M)       
  Spin = NINT(TargetStates2el%b(nstate)%spin)
  parity = TargetStates2el%b(nstate)%parity
  n_orb_max = TargetStates2el%b(nstate)%nusemax     
  
  dM = map - ma 

  ! MAKE SURE BOTH MOLECULAR STATES HAVE THE SAME SPIN and difference Angular Projections <= 1
  ! Parity must change
  if ( Spin /= Spinp ) return
  if (ABS(dM) > 1 ) return
  if ( parity == parityp) return


! Below sums all the overlaps for COORDINATE 2  same configuratins(1s,2s,..) in COORDINATE 1
  ! FINAL State COORDINATE 1
  do np_orb = 1, np_orb_max
   
     nuse_nporb= TargetStates2el%b(npstate)%nuse(np_orb)
     
     ! INITIAL State COORDINATE 1
     do n_orb = 1, n_orb_max
      
        nuse_norb =  TargetStates2el%b(nstate)%nuse(n_orb)   
       
        overlap2 = 0d0

        ! Looping over FINAL Molecular State orbitals. COORDINATE 2
        do nep_con =  1, Npcon          
           
           indp1_2 = TargetStates2el%b(npstate)%na(nep_con)    ! Final state number np. nep A.O.
           
           if ( nuse_nporb /= indp1_2 ) cycle
         
           ! Quantum numbers for ALPHA
           tnp1 => bst%b(indp1_2)                           !        
           lap1 = get_ang_mom(tnp1)         ! Gets Angular momentum A.O.
           mp1 = get_ang_mom_proj(tnp1)     ! Get angular projection of A.O.   
           
           ! Quantum numbers and functions for BETA
           indp2 = TargetStates2el%b(npstate)%nb(nep_con)                
           tnp2 => bst%b(indp2)                                         
           fp2 => fpointer(tnp2)                                
           lap2 = get_ang_mom(tnp2)                             
           mp2 = get_ang_mom_proj(tnp2)          
           CIp = get_CI(TargetStates2el%b(npstate),nep_con)  
           
           ! Looping over INITIAL Molecular State orbitals.  COORDINATE 2
           do ne_con = 1, Ncon       
              
              ind1_2 = TargetStates2el%b(nstate)%na(ne_con)  
              
              if ( nuse_norb /= ind1_2 ) cycle
              
              ! Quantum numbers for GAMMA
              tn1 => bst%b(ind1_2)                                                                                     
              la1 = get_ang_mom(tn1)                            
              m1 = get_ang_mom_proj(tn1)   

              
              ! Quantum numbers and functions for DELTA
              ind2 = TargetStates2el%b(nstate)%nb(ne_con) 
              tn2 => bst%b(ind2)                                         
              f2 => fpointer(tn2)                               
              la2 = get_ang_mom(tn2)                            
              m2 = get_ang_mom_proj(tn2)            
              CI = get_CI(TargetStates2el%b(nstate),ne_con)  
              

              dL = lap1 - la1
              dM = mp1 - m1
              
              ! Selections Rules 
              if ( ABS(dL) > 1 .OR. ABS(dM) > 1 ) cycle    
              if ( lap2 /= la2 .OR. mp2 /= m2  ) then
                 cycle
              else
                 ! Saving Last Angular Projection used so can get mp1 and m1
                 m1_temp = m1
                 mp1_temp = mp1
              end if
                               
              temp_overlap2 = 0d0
              ! DO OVERLAPOF COORDINATE SPACE 2  
              minf = get_minf(tn2)
              maxf = get_maxf(tn2)
              minfp = get_minf(tnp2)
              maxfp = get_maxf(tnp2)      
              ir1 = max(minf,minfp)
              ir2 = min(maxf,maxfp)  
              temp_overlap2 = SUM( fp2(ir1:ir2) *  f2(ir1:ir2) * weight(ir1:ir2) )
              
              overlap2 = overlap2 +  CIp * CI * temp_overlap2
              
           end do  ! INITIAL STATE COORDINATE 2
           
        end do    ! FINAL STATE COORDINATE 2

        if ( overlap2 == 0d0 ) cycle


        ! COORDINATE 1 RADIAL INTEGRALS
        ! Quantum numbers and functions for ALPHA
        indp1 = nuse_nporb                                  ! Final state number np. nep A.O.
        tnp1 => bst%b(indp1)                             !           
        fp1 => fpointer(tnp1)                               ! One electron functions
        lap1 = get_ang_mom(tnp1)                            ! Gets Angular momentum A.O.
        mp1 = get_ang_mom_proj(tnp1)                        ! Get angular projection of A.O. 
        

        ! Quantum numbers and functions for GAMMA
        ind1 = nuse_norb               
        tn1 => bst%b(ind1)                                          
        f1 => fpointer(tn1)                               
        la1 = get_ang_mom(tn1)                            
        m1 = get_ang_mom_proj(tn1)     
        
        dL = lap1 - la1
        dM = mp1 - m1
        ! Selections Rules 
        if ( ABS(dL) > 1 .OR. ABS(dM) > 1 ) cycle    

        minf = get_minf(tn1)
        maxf = get_maxf(tn1)
        minfp = get_minf(tnp1)
        maxfp = get_maxf(tnp1)
        ir1 = max(minf,minfp)
        ir2 = min(maxf,maxfp) 
        ! Radial Integration <fp1|r|f1>
        r_overlap1 = SUM( fp1(ir1:ir2) *  f1(ir1:ir2)  * gridr(ir1:ir2) * weight(ir1:ir2) ) 
        
        !Velocity Integrals
        ! Derivatiie Integral Note: d/dr(phi(r)/r) = phi'(r)/r -  phi(r)/(r^2)
        ideriv1 =  SUM((( f1(ir1 + 1 : ir2 ) - f1(ir1: ir2 - 1)) /(gridr(ir1 + 1 :ir2) - gridr(ir1: ir2 - 1))) * fp1(ir1: ir2 - 1) * weight(ir1:ir2 - 1))
        vint1 = SUM( fp1(ir1:ir2) * f1(ir1:ir2) * weight(ir1:ir2) / gridr(ir1:ir2) )
        
        if ( la1 < lap1 ) then
           vint1 = ideriv1 - (dble(la1) + 1.0) * vint1
        else if ( la1 > lap1) then
           vint1 = ideriv1 + dble(la1) * vint1
        end if
        
        temp = (2.0 * la1 + 2.0 + dL) * (2.0 * la1 + dL)
        Coeff = 0d0 ! Length Guage Block Varsholovich pg 145 
        if ( dM == 0) then ! Parallel Transitions <z> 
           Coeff = sqrt( ( dble(la1) + 0.5 + dble(dL) * 0.5 ) ** 2.0 - dble(m1 * m1)) 
        end if
        if ( ABS(dM) == 1) then
           if ( la1 < lap1   ) then  ! Perpendicular Transitions <x +/- iy> 
              Coeff = dble(dM)*sqrt(dble(la1+dM*m1+2)*dble(la1+dM*m1+1))
           else if ( la1 > lap1) then
              Coeff = -dble(dM)*sqrt(dble(la1-dM*m1)*dble(la1-dM*m1-1))
           end if
           Coeff = Coeff /  sqrt(2.0)                 
        end if
        
        Coeff = Coeff / sqrt(temp) 
        
        ! Multiply by 2 for 2e target. Operator L: z_1 + z_2 = 2*z_1
        result_l = temp_l + 2.0 * Coeff * r_overlap1 * overlap2
        result_v = temp_v + 2.0 * Coeff * vint1 * overlap2
        temp_l = result_l
        temp_v = result_v
        
     end do ! INITIAL STATE COORDINATE 1
     
  end do   ! FINAL STATE COORDINATE 1

  dM = map - ma 
  ! TranEnergy is in a.u. converted to Rydbergs here.
  if ( ABS(dM) == 1) then
     result_l = 4.0 * ABS(TranEnergy) * result_l * result_l / 3.0 ! Added orbital degeneracy g = 2, should test for ma=+/-2 etc 
     result_v = 4.0 * result_v * result_v / (3.0 * ABS(TranEnergy))
  else 
     result_l = 2.0 * ABS(TranEnergy) * result_l * result_l / 3.0 
     result_v = 2.0 * result_v * result_v / ( 3.0  * ABS(TranEnergy) )
  end if


 
end subroutine oscstrength_2e_config
