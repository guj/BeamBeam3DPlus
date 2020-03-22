program main
  use sensei_tunefoot
  use mpi
  implicit none
  !include 'mpif.h'                                                                                                                                                                                       
  double precision :: time
  integer :: ierr;
  integer :: xPos, yPos

  call MPI_INIT(ierr)

!---------------------------------------------------------------------
! NOTE: all of the following parameters are::: 0 based for C libraries
! 
!   format: 
!      function_name (start, count)
!---------------------------------------------------------------------

  !call tunefoot_set_particleRange(0, 100)
  !call tunefoot_set_bunch(0)
  !call tunefoot_setTurnRange(0, 1024)

  !xPos = 0; 
  !yPos = 2

  !call tunefoot_run(xPos, yPos) 
  call tunefoot_run_xml(MPI_COMM_WORLD, ierr)
  call MPI_Finalize(ierr);

end program main

