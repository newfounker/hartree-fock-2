module grid_radial
 implicit none
 public

! This is declaration of the type rgrid.
 type rgrid
    integer:: npwave                      ! the number of points per oscillation
    integer:: npdbl                       ! the number of points with the same dx per interval
    real*8:: rmax                         ! maximum value of r
    real*8:: qmax                         ! maximum value of q for which integration is accurate
    integer:: ndouble                     ! number of doubling
    real*8:: formcut
    real*8:: regcut
    real*8:: expcut
! the above parameters are used to set up r-grid decribed below
    real*8, dimension(:), pointer:: gridr   ! r-grid values
    real*8, dimension(:), pointer:: weight  ! r-grid weights for Simpson rule
    real*8, dimension(:), pointer:: bool    ! r-grid weights for Bool rule
    integer:: nr                          ! number of points in r-grid
    integer, dimension(:), pointer:: jdouble ! j points where dx doubles

! rpow_f(r,l) = r**l, rpow_g(r,l) = 1/r**(l+1)
! i2_f = nr,  i1_g = 1
    integer:: ltmax                       ! used to set up arrays power_f(i,l) and power_g(i,l) for e-e Coulomb potential integration
    real*8, dimension(:,:), allocatable:: rpow_f, rpow_g
    integer, dimension(:), allocatable:: i1_f, i2_g
! note that i1_g = 1, and i2_f = nr

end type rgrid

!****  This is declaration of an object of the type rgrid.

 type(rgrid):: grid

!****

! This is public method for type rgrid
 public setgrids

contains

!======================================================================
!  This subroutine sets up radial grid:  grid. It is used to perform integrations
!  from zero to infinity using Simpson's rule. The structure is such
!  that a wave with momentum QMAX will have NPWAVE points per oscillation. This
!  determines HMAX, the largest step. At some point, XDBLE, the intervals,
!  dX, are progressively halved. This is done to accomodate irregular
!  solutions.
!  INPUT:
!    npwave  - number of points per oscillation
!    qmax  - the biggest momentum (a.u.) that can be calculated by this mesh.
!    ndoble - the number of doubling of intervals is NDOUBLE
!    rmax  - the largest "r=x" in the meshes.
!  OUTPUT:
!    gridr -  R grid
!    nr    - Number of X points
!    jdouble - j points where dx doubles
!    weight -  weights for Simpson rule
!    bool   -  weights for Bool rule
!
! It is assumed that wave function is exact zero at r=0 and at r=rmax (or at
! corrsponding maxf, minf points). Therefore the first (r=0) and last (r=rmax)
! terms in Simpson composed rule are exact zero and r grid starts from
! the second term and finishes at the one before last.
!======================================================================
subroutine setgrids(self)
  use input_data
  implicit none
  type(rgrid) self

  integer:: npwave, npdbl, nleft, j, nj, nd, max_nj, jb
  real*8:: hmax, hmin, rdble, rleft, r, dr
  integer:: lna

! This is r-grid parameters used in setgrids(grid)
  self%npwave = data_in%npwave
  self%npdbl = data_in%npdbl
  self%qmax = data_in%qmax
  self%ndouble = data_in%ndouble
  self%rmax = data_in%rmax
  self%formcut = data_in%formcut
  self%regcut = data_in%regcut
  self%expcut = data_in%expcut
  self%ltmax = data_in%ltmax



!  jdouble stores the J index for GRIDR where DR doubles, with the first
!  point set to 1 and the last to NR. Used in the numerov integrations.

!  Set up GRIDR

!  Make sure NPDBL is even
  self%npdbl=(self%npdbl/2) * 2
  if (self%npdbl.lt.4) then
     print*,'Need to have at least 4 equally spaced points in GRIDR'
     stop 'Need to have at least 4 equally spaced points in GRIDR'
  end if

  hmax = 3.14/float(self%npwave)/self%qmax
  hmin = hmax/float(2**self%ndouble)
!  The value of the R point from which dR is constant, = HMAX, is RDBLE
!  RDBLE = NPDBL * hmin * (2**NDOUBLE-1)
  rdble = float(self%npdbl) * hmin * (2**self%ndouble-1)
!  The remaining part from rdble to rmax is:
  rleft = self%rmax - rdble
!  nleft = int(rleft / hmax) / 2 * 2
  nleft = int(rleft / hmax) / 4 * 4
!  The total number of points:
  self%nr = nleft + self%npdbl * self%ndouble

!??  print*,'Grid R parameters:'
!??  print*,'NDOUBLE:',self%ndouble
!??  print*,'HMIN:',hmin
!??  print*,'HMAX:',hmax
!??  print*,'NPDBL:',self%npdbl
!??  print*,'RDBLE:',rdble, self%nr - nleft
!??  print*,'NR:', self%nr


  allocate(self%gridr(self%nr))
  allocate(self%weight(self%nr))
  allocate(self%bool(self%nr))
  allocate(self%jdouble(self%ndouble+2))


  self%jdouble(1) = 1
  do nd = 2, self%ndouble + 1
     self%jdouble(nd) = self%npdbl * (nd - 1)
  end do
  self%jdouble(self%ndouble+2) = self%nr

  dr = hmin
  r = 0.0
  j = 0
  jb = 0
  do nd = 1, self%ndouble+1
! For all intervals except for the last one max_nj should be equal to npdbl, for first interval it does not give corect result and is corrected in the line below, for the last interval it should give nleft.
     max_nj = self%jdouble(nd+1) - self%jdouble(nd)
     if(nd .eq. 1) max_nj = max_nj + 1
     do nj = 1, max_nj
        j = j + 1
        self%gridr(j) = r + float(nj) * dr
!  Simpson's rule weights
        self%weight(j) = float(mod(j,2) * 2 + 2) * dr / 3.0
     end do
     self%weight(j) = dr  !  this accounts for change of integration step at the boundary: (dr + 2*dr)/3 = dr
     r = self%gridr(j)
!
!  Bool's rule weights
     do nj=1,max_nj,4
        self%bool(jb+nj) = dr * 32.
        self%bool(jb+nj+1) = dr * 12.
        self%bool(jb+nj+2) = dr * 32.
        self%bool(jb+nj+3) = dr * 14.
     end do
     jb = self%jdouble(nd+1)
     self%bool(jb) = dr * 21.
!
     dr = dr * 2.0
  end do
  self%weight(j) = hmax/3.0  ! for j=nr
  self%bool(self%nr) = 7.
  self%bool(1:self%nr) = self%bool(1:self%nr) * 2.0 / 45.0

!??  print*,'Last R and NR:', self%gridr(self%nr), self%nr


!  call  rpow_construct(self%ltmax,self%nr,self%gridr,self%regcut,self%expcut)


!??  print*, 'Set rpow, ltmax =', self%ltmax

  allocate(self%rpow_f(1:self%nr,0:self%ltmax))
  allocate(self%rpow_g(1:self%nr,0:self%ltmax))
  allocate(self%i1_f(0:self%ltmax))
  allocate(self%i2_g(0:self%ltmax))

  self%i1_f(0) = 1
  self%i2_g(0) = self%nr
  self%rpow_f(1:self%nr,0)=1.0
  self%rpow_g(1:self%nr,0)=1.0/self%gridr(1:self%nr)

  do lna=1,self%ltmax

     self%i1_f(lna) = self%i1_f(lna-1)
     self%i2_g(lna) = self%i2_g(lna-1)
     self%rpow_f(1:self%nr,lna)=self%rpow_f(1:self%nr,lna-1)*self%gridr(1:self%nr)
     self%rpow_g(1:self%nr,lna)=self%rpow_g(1:self%nr,lna-1)/self%gridr(1:self%nr)

     do while (self%rpow_f(self%i1_f(lna),lna) .lt. self%regcut)
        self%i1_f(lna) = self%i1_f(lna)+1
     end do
     do while (self%rpow_g(self%i2_g(lna),lna) .lt. self%expcut*self%expcut)
        self%i2_g(lna) = self%i2_g(lna)-1
     end do

  end do


end subroutine setgrids

subroutine destruct_gridr(self)

  implicit none
  type(rgrid) self

  deallocate(self%rpow_f)
  deallocate(self%rpow_g)
  deallocate(self%i1_f)
  deallocate(self%i2_g)
  deallocate(self%gridr,self%weight,self%bool)

end subroutine destruct_gridr

end module grid_radial
