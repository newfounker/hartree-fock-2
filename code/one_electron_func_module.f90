module one_electron_func_mod

  use sturmian_class

  integer:: basis_type ! =0 or 1 if  spherical basis is in use, or =2 if spheroidal.
  type(basis_sturmian_nr) :: bst   ! this is Sturmian basis
  integer:: nspm

end module one_electron_func_mod
