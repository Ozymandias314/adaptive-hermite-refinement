MODULE mpi_mod

#ifdef mpi1
  IMPLICIT NONE
  INCLUDE 'mpif.h'
#else
  USE mpi
#endif

END MODULE mpi_mod
