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

  character(len = 200) :: input_file, hf_file, curve_file
  integer              :: ii

  write (*, "(a)") "> main"

  if (command_argument_count() < 3) then

    write (*, *) "missing (< 3) arguments (input, hf and potential curve files)"
    write (*, *) "assuming default values"
    write (*, *)

    input_file = "../input/input.dat"
    hf_file = "../output/hf_results.dat"
    curve_file = "../output/curve.dat"

  else

    call get_command_argument(1, input_file)
    call get_command_argument(2, hf_file)
    call get_command_argument(3, curve_file)

  end if

  write (*, "(a, a)") &
      " input_file: ", input_file
  write (*, "(a, a)") &
      " hf_file:    ", hf_file
  write (*, "(a, a)") &
      " curve_file: ", curve_file


  call readin(data_in, trim(input_file), .true.)

  call setgrids(grid)

  call construct_vnc(grid%nr, grid%gridr)

  call hf_structure (trim(hf_file))

  ! Determine number of one-electron target states
  Number_one_electron_func = 0

  do ipar= -1, 1, 2
    Number_one_electron_func = Number_one_electron_func + &
        SUM(data_in%nst(data_in%Mt_min:data_in%Mt_max, ipar))

!!$ account for degeneracy of the states with nonzwero M
    Number_one_electron_func = Number_one_electron_func + &
        SUM(data_in%nst(max(1,data_in%Mt_min):data_in%Mt_max, ipar))
  end do

  call construct_1el_basis_nr(Number_one_electron_func, &
      trim(hf_file), trim(curve_file))

  stop

end program main
