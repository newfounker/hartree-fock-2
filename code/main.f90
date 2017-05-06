program main


  use input_data            ! definitions for input data type, construct routine, object data_in
  use grid_radial           ! set up r-grid, and keep it in object grid of type(gridr)
  use sturmian_class        ! defines one-electron functions (and basis of them) with operations on them
  use vnc_module        !  keeps central potential
  use one_electron_func_mod     ! keeps basis of one-electron functions
  use structure

  implicit none
!
  integer:: Number_one_electron_func, ipar

  character(len = 200) :: input_filepath, output_filepath, arg
  integer              :: ii

  write (*, "(a)") "> main"

  if (command_argument_count() < 2) then

    write (*, *) "missing (2) arguments: input directory and output directory"
    write (*, *) "assuming default values"
    write (*, *)

    input_filepath = "../input/input.dat"
    output_filepath = "../output/hf_results.dat"

  else

    call get_command_argument(1, input_filepath)
    call get_command_argument(2, output_filepath)

  end if

  call readin(data_in, trim(input_filepath), .true.)

  call setgrids(grid)

  call construct_vnc(grid%nr, grid%gridr)

  call hf_structure (0, trim(output_filepath))

!------------------------------------------------------------------------
! Input data: construct routine for input type.
!
!-------------------------------------------------------------------------
! construct routine for radial grid
!
!-------------------------------------------------------------------------
!
! Make non-central local potential
!
!??  print*

  !< assuming only m = 0 basis functions are used.

!< tom ross: removed for hf_calculation
! !??  print*, 'Start structure calculation'
!   ! Determine number of one-electron target states
!   Number_one_electron_func = 0

!   do ipar= -1, 1, 2
!     Number_one_electron_func = Number_one_electron_func + &
!         SUM(data_in%nst(data_in%Mt_min:data_in%Mt_max, ipar))

! !!$ account for degeneracy of the states with nonzwero M
!     Number_one_electron_func = Number_one_electron_func + &
!         SUM(data_in%nst(max(1,data_in%Mt_min):data_in%Mt_max, ipar))
!   end do


!   call construct_1el_basis_nr(Number_one_electron_func, trim(output_filepath))

! !------------------------------------------



  stop

end program main
