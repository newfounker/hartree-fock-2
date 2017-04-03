module vnc_module
  use sturmian_class
  integer:: lamtop_vc ! limit for vnc
  real*8, dimension(:,:), allocatable:: vnc
  integer, dimension(:), allocatable:: minvnc, maxvnc

end module vnc_module
