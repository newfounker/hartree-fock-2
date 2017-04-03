program main


  use input_data            ! definitions for input data type, construct routine, object data_in
  use grid_radial           ! set up r-grid, and keep it in object grid of type(gridr)
  use sturmian_class        ! defines one-electron functions (and basis of them) with operations on them
  use vnc_module        !  keeps central potential
  use one_electron_func_mod     ! keeps basis of one-electron functions

  implicit none
!
  integer:: Number_one_electron_func, ipar


!------------------------------------------------------------------------
! Input data: construct routine for input type.
  call readin(data_in, "data.in", .true.)
!
!-------------------------------------------------------------------------
! construct routine for radial grid
  call setgrids(grid)
!
!-------------------------------------------------------------------------
!
! Make non-central local potential
  call construct_vnc(grid%nr, grid%gridr)
!
!??  print*

!??  print*, 'Start structure calculation'
  ! Determine number of one-electron target states
  Number_one_electron_func = 0

  do ipar= -1, 1, 2
    Number_one_electron_func = Number_one_electron_func + &
        SUM(data_in%nst(data_in%Mt_min:data_in%Mt_max, ipar))

!!$ account for degeneracy of the states with nonzwero M
    Number_one_electron_func = Number_one_electron_func + &
        SUM(data_in%nst(max(1,data_in%Mt_min):data_in%Mt_max, ipar))
  end do


  call construct_1el_basis_nr(Number_one_electron_func)

!------------------------------------------



  stop

end program main
