module data_targets
  ! This data is for the molecule H^{+}_{2}. 
  ! Compare calculated values with exact calculations for Energies and Oscillator Strengths, which are stored here.
  
  
  !Need to learn modules mark if you want to allocate arrays in a subroutine
  ! This will allow me to call the arrays without reallocating them every time in the do loops
  
  implicit none 
  real:: Ex_osc, Ex_Energy
  integer:: Ref

contains
  real function H2Iosc(np,n)
    use input_data
    implicit none
    integer, intent(in):: np, n
    real*8, allocatable:: H2Iosc_exact(:,:)
    integer :: npmax, nmax
    
    npmax = 30
    nmax =  30
    
    allocate( H2Iosc_exact(1:npmax,2:nmax) )
    
    ! Write Once. Written up completely for up to n = 4, i.e. for Mt_max = n - 1
    if (  ABS(data_in%Mt_max) >= 1 ) then ! Buggy if i try to declar just once with condition np == 1 .AND. n == 2
       
       H2Iosc_exact(:,:) = 0d0
       H2Iosc_exact(1,2) = 0.319    ! 1s(m=0)-2p(m=0)
       H2Iosc_exact(1,3) = 0.46     ! 1s(m=0)-2p(m=-1)  
       H2Iosc_exact(1,4) = 0.46     ! 1s(m=0)-2p(m=1)
       H2Iosc_exact(1,6) = 8.24E-4  ! 1s(m=0)-3p(m=0) 
       H2Iosc_exact(2,5) = 1.36E-1  ! 2p(m=0)-2s(m=0)
       H2Iosc_exact(2,8) = 0.268    ! 2p(m=0)-3d(m=-1)
       H2Iosc_exact(2,9) = 0.268    ! 2p(m=0)-3d(m=1)
       H2Iosc_exact(2,7) = 2.18E-1  ! 2p(m=0)-3d(m=0)  
       H2Iosc_exact(3,8) = 0.275    ! 2p(m=-1)-3d(m=-1)
       H2Iosc_exact(4,9) = 0.275    ! 2p(m=1)-3d(m=1)
       
       if (  ABS(data_in%Mt_max) == 1 ) then
          
          H2Iosc_exact(1,13) = 5.54E-5  ! 1s(m=0)-4p(m=0)             
          H2Iosc_exact(1,15) = 4.07E-5  ! 1s(m=0)-4f(m=0)  
          H2Iosc_exact(2,12) = 1.46E-2  ! 2p(m=0)-3s(m=0)
          
       else if (  ABS(data_in%Mt_max) >= 2 ) then
          
          H2Iosc_exact(1,15) = 5.54E-5  ! 1s(m=0)-4p(m=0)             
          H2Iosc_exact(1,17) = 4.07E-5  ! 1s(m=0)-4f(m=0)  
          H2Iosc_exact(2,14) = 1.46E-2  ! 2p(m=0)-3s(m=0)
          
       end if
       
    end if
    
    H2Iosc = 0d0
    
    if ( np <= npmax .AND. n <= nmax  ) then
       H2Iosc = H2Iosc_exact(np,n)
    end if
    
  end function H2Iosc
  
  integer function Ref_H2Iosc(np,n)
    use input_data
    integer, intent(in):: np, n
    integer, allocatable:: Ref_H2Iosc_exact(:,:)
    integer :: npmax, nmax
    
    npmax = 30
    nmax =  30
    
    allocate ( Ref_H2Iosc_exact(1:npmax,2:nmax))
    
    ! Ref[0]: None, Ref[1]: D.R. Bates 1954, Ref[2]:, Ref[3]: D.R. Bates 1951

    if (  ABS(data_in%Mt_max) >= 1 ) then ! Buggy if i try to declar just once with condition np == 1 .AND. n == 2
       
       Ref_H2Iosc_exact(:,:) = 0
       Ref_H2Iosc_exact(1,2) = 3      ! 1s(m=0)-2p(m=0)
       Ref_H2Iosc_exact(1,3) = 2      ! 1s(m=0)-2p(m=-1)  
       Ref_H2Iosc_exact(1,4) = 2      ! 1s(m=0)-2p(m=1)
       Ref_H2Iosc_exact(1,6) = 1      ! 1s(m=0)-3p(m=0) 
       Ref_H2Iosc_exact(2,5) = 1      ! 2p(m=0)-2s(m=0)
       Ref_H2Iosc_exact(2,8) = 2      ! 2p(m=0)-3d(m=-1)
       Ref_H2Iosc_exact(2,9) = 2      ! 2p(m=0)-3d(m=1)
       Ref_H2Iosc_exact(2,7) = 1      ! 2p(m=0)-3d(m=0)  
       Ref_H2Iosc_exact(3,8) = 2      ! 2p(m=-1)-3d(m=-1)
       Ref_H2Iosc_exact(4,9) = 2      ! 2p(m=1)-3d(m=1)
       
       if (  ABS(data_in%Mt_max) == 1 ) then
          
          Ref_H2Iosc_exact(1,13) = 1  ! 1s(m=0)-4p(m=0)             
          Ref_H2Iosc_exact(1,15) = 1  ! 1s(m=0)-4f(m=0)  
          Ref_H2Iosc_exact(2,12) = 1  ! 2p(m=0)-3s(m=0)
          
       else if (  ABS(data_in%Mt_max) >= 2 ) then
          
          Ref_H2Iosc_exact(1,15) = 1  ! 1s(m=0)-4p(m=0)             
          Ref_H2Iosc_exact(1,17) = 1  ! 1s(m=0)-4f(m=0)  
          Ref_H2Iosc_exact(2,14) = 1  ! 2p(m=0)-3s(m=0)
          
       end if
       
    end if
    
    Ref_H2Iosc = 0
    
    if ( np <= npmax .AND. n <= nmax  ) then
       Ref_H2Iosc = Ref_H2Iosc_exact(np,n)
    end if
    
  end function Ref_H2Iosc


  real function H2Ieny(np)   
    use input_data
    implicit none
    integer, intent(in):: np
    real, allocatable:: H2Ieny_exact(:)
    integer:: npmax, Mt

    Mt = ABS(data_in%Mt_max)
    
    npmax = 30

    allocate( H2Ieny_exact(0:npmax) )
    
    ! Exact Values in atomic units. Written up completely for up to n = 4, i.e. for Mt_max = n - 1
    
    if (  Mt  >= 1 ) then ! Buggy when i try to declar array just once with np == 1
       H2Ieny_exact(:) = 0d0
       H2Ieny_exact(1) = -0.602635   ! 1s(m=0)
       H2Ieny_exact(2) = -0.167535   ! 2p(m=0)
       H2Ieny_exact(3) = 0.071229    ! 2p(m=1) 
       H2Ieny_exact(4) = 0.071229    ! 2p(m=1)
       H2Ieny_exact(5) = 0.139135    ! 2s(m=0)
       H2Ieny_exact(6) = 0.244586    ! 3p(m=0)
       H2Ieny_exact(7) = 0.264223    ! 3d(m=0) 
       H2Ieny_exact(8) = 0.2733      ! 3d(m=-1)
       H2Ieny_exact(9) = 0.2733      ! 3d(m=-1) 
       
       if (ABS(data_in%Mt_max) == 1) then  
          
          H2Ieny_exact(10) = 0.d0      ! 3p(m=-1) 
          H2Ieny_exact(11) = 0.d0      ! 3p(m=1)
          H2Ieny_exact(12) = 0.14464   ! 3s(m=0) Bates et al 1953 
          H2Ieny_exact(13) = 0.362685  ! 4p(m=0) Bates et al 1953
          H2Ieny_exact(14) = 0d0       ! 4d(m=0)
          H2Ieny_exact(15) = 0.373356  ! 4f(m=0)  
          H2Ieny_exact(16) = 0d0       ! 4d(m=-1)   
          H2Ieny_exact(17) = 0d0       ! 4d(m=1)
          H2Ieny_exact(18) = 0d0       ! 4f(m=-1)   
          H2Ieny_exact(19) = 0d0       ! 4f(m=1)  
          H2Ieny_exact(20) = 0.384085  ! 4p(m=-1)   
          H2Ieny_exact(21) = 0.384085  ! 4p(m=1)  
          
       else if ( ABS(data_in%Mt_max) == 2) then
          
          H2Ieny_exact(10) = 0.d0     ! 3d(m=-2) 
          H2Ieny_exact(11) = 0.d0     ! 3d(m=2)    
          H2Ieny_exact(12) = 0.d0     ! 3p(m=-1) 
          H2Ieny_exact(13) = 0.d0     ! 3p(m=1)
          H2Ieny_exact(14) = 0.32232  ! 3s(m=0) Bates et al 1953 
          H2Ieny_exact(15) = 0.362685 ! 4p(m=0) Bates et al 1953
          H2Ieny_exact(16) = 0d0      ! 4d(m=0)
          H2Ieny_exact(17) = 0.373356 ! 4f(m=0)   
          H2Ieny_exact(18) = 0d0      ! 4d(m=-1)   
          H2Ieny_exact(19) = 0d0      ! 4d(m=1)
          H2Ieny_exact(20) = 0d0      ! 4f(m=-1)   
          H2Ieny_exact(21) = 0d0      ! 4f(m=1)    
          H2Ieny_exact(22) = 0d0      ! 4f(m=-2)   
          H2Ieny_exact(23) = 0d0      ! 4f(m=2)
          H2Ieny_exact(24) = 0d0      ! 4d(m=-2)   
          H2Ieny_exact(25) = 0d0      ! 4d(m=2)   
          H2Ieny_exact(26) = 0.384085 ! 4p(m=-1)   
          H2Ieny_exact(27) = 0.384085 ! 4p(m=1)  
          H2Ieny_exact(28) = 0.394588 ! 4s(m=0)
            
       else if ( ABS(data_in%Mt_max) >= 3) then
          
          H2Ieny_exact(10) = 0.d0     ! 3d(m=-2) 
          H2Ieny_exact(11) = 0.d0     ! 3d(m=2)    
          H2Ieny_exact(12) = 0.d0     ! 3p(m=-1) 
          H2Ieny_exact(13) = 0.d0     ! 3p(m=1)
          H2Ieny_exact(14) = 0.14464  ! 3s(m=0) Bates et al 1953 
          H2Ieny_exact(15) = 0.362685 ! 4p(m=0) Bates et al 1953
          H2Ieny_exact(16) = 0d0      ! 4d(m=0)
          H2Ieny_exact(17) = 0.373356 ! 4f(m=0)   
          H2Ieny_exact(18) = 0d0      ! 4d(m=-1)   
          
          H2Ieny_exact(19) = 0d0      ! 4d(m=1)
          H2Ieny_exact(20) = 0d0      ! 4f(m=-1)   
          H2Ieny_exact(21) = 0d0      ! 4f(m=1)    
          H2Ieny_exact(22) = 0d0      ! 4f(m=-2)   
          H2Ieny_exact(23) = 0d0      ! 4f(m=2)
            
          H2Ieny_exact(24) = 0d0      ! 4f(m=-3) 
          H2Ieny_exact(25) = 0d0      ! 4f(m=3)  
          
          H2Ieny_exact(26) = 0d0      ! 4d(m=-2)   
          H2Ieny_exact(27) = 0d0      ! 4d(m=2)   
          H2Ieny_exact(28) = 0.384085 ! 4p(m=-1)   
          H2Ieny_exact(29) = 0.384085 ! 4p(m=1)  
          H2Ieny_exact(30) = 0.394588 ! 4s(m=0)
          
       end if
       
       H2Ieny_exact(0:npmax) = H2Ieny_exact(0:npmax) - 0.5
       
    end if

    H2Ieny = 0d0
    if ( np <= npmax ) then
       H2Ieny = H2Ieny_exact(np) 
    end if
  end function H2Ieny

  integer function Ref_H2Ieny(np)
    use input_data
    implicit none
    integer, intent(in):: np
    integer, allocatable:: H2Ieny_Ref(:)
    integer:: npmax
    
    npmax = 30

    allocate( H2Ieny_Ref(0:npmax) )
    
    !  Ref[1]: T.E. Sharp 1970. Ref[2]: D.R. Bates et al 1953. Need to get Ref for the energy havent found yet
    
    if (  data_in%Mt_max  >= 1 ) then ! Buggy when i try to declar array just once with np == 1
      
       H2Ieny_Ref(1:9) = 1 ! T.E Sharp 1970
       
       if (ABS(data_in%Mt_max) == 1) then            
          
          H2Ieny_Ref(10) = 0      ! 3p(m=-1) 
          H2Ieny_Ref(11) = 0      ! 3p(m=1)
          H2Ieny_Ref(12) = 2      ! 3s(m=0) Bates et al 1953 
          H2Ieny_Ref(13) = 2      ! 4p(m=0) Bates et al 1953
          H2Ieny_Ref(14) = 0      ! 4d(m=0)
          H2Ieny_Ref(15) = 1      ! 4f(m=0)  
          H2Ieny_Ref(16) = 0      ! 4d(m=-1)   
          H2Ieny_Ref(17) = 0      ! 4d(m=1)
          H2Ieny_Ref(18) = 0      ! 4f(m=-1)   
          H2Ieny_Ref(19) = 0      ! 4f(m=1)  
          H2Ieny_Ref(20) = 1      ! 4p(m=-1)   
          H2Ieny_Ref(21) = 1      ! 4p(m=1) 
          
       else if ( ABS(data_in%Mt_max) == 2) then          
          
          H2Ieny_Ref(10) = 0      ! 3d(m=-2) 
          H2Ieny_Ref(11) = 0      ! 3d(m=2)    
          H2Ieny_Ref(12) = 0      ! 3p(m=-1) 
          H2Ieny_Ref(13) = 0      ! 3p(m=1)
          H2Ieny_Ref(14) = 2      ! 3s(m=0) Bates et al 1953 
          H2Ieny_Ref(15) = 2      ! 4p(m=0) Bates et al 1953
          H2Ieny_Ref(16) = 0      ! 4d(m=0)
          H2Ieny_Ref(17) = 1      ! 4f(m=0)   
          H2Ieny_Ref(18) = 0      ! 4d(m=-1)   
          H2Ieny_Ref(19) = 0      ! 4d(m=1)
          H2Ieny_Ref(20) = 0      ! 4f(m=-1)   
          H2Ieny_Ref(21) = 0      ! 4f(m=1)    
          H2Ieny_Ref(22) = 0      ! 4f(m=-2)   
          H2Ieny_Ref(23) = 0      ! 4f(m=2)
          H2Ieny_Ref(24) = 0      ! 4d(m=-2)   
          H2Ieny_Ref(25) = 0      ! 4d(m=2)   
          H2Ieny_Ref(26) = 1      ! 4p(m=-1)   
          H2Ieny_Ref(27) = 1      ! 4p(m=1)  
          H2Ieny_Ref(28) = 1      ! 4s(m=0)

       else if ( ABS(data_in%Mt_max) >= 3) then
          

          H2Ieny_Ref(10) = 0      ! 3d(m=-2) 
          H2Ieny_Ref(11) = 0      ! 3d(m=2)    
          H2Ieny_Ref(12) = 0      ! 3p(m=-1) 
          H2Ieny_Ref(13) = 0      ! 3p(m=1)
          H2Ieny_Ref(14) = 2      ! 3s(m=0) Bates et al 1953 
          H2Ieny_Ref(15) = 2      ! 4p(m=0) Bates et al 1953
          H2Ieny_Ref(16) = 0      ! 4d(m=0)
          H2Ieny_Ref(17) = 1      ! 4f(m=0)   
          H2Ieny_Ref(18) = 0      ! 4d(m=-1)   
          
          H2Ieny_Ref(19) = 0      ! 4d(m=1)
          H2Ieny_Ref(20) = 0      ! 4f(m=-1)   
          H2Ieny_Ref(21) = 0      ! 4f(m=1)    
          H2Ieny_Ref(22) = 0      ! 4f(m=-2)   
          H2Ieny_Ref(23) = 0      ! 4f(m=2)
            
          H2Ieny_Ref(24) = 0      ! 4f(m=-3) 
          H2Ieny_Ref(25) = 0      ! 4f(m=3)  
          
          H2Ieny_Ref(26) = 0      ! 4d(m=-2)   
          H2Ieny_Ref(27) = 0      ! 4d(m=2)   
          H2Ieny_Ref(28) = 1      ! 4p(m=-1)   
          H2Ieny_Ref(29) = 1      ! 4p(m=1)  
          H2Ieny_Ref(30) = 1      ! 4s(m=0)
          
       end if
       
    end if
        
    Ref_H2Ieny = 0

    if ( np <= npmax) then
       Ref_H2Ieny = H2Ieny_Ref(np) 
    end if
        
  end function Ref_H2Ieny
  
end module data_targets



