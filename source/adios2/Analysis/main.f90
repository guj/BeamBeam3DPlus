program main
  use adios2_tunefoot
  use mpi

  implicit none

  !!character(len=*), parameter :: version = '1.0'
  integer :: ierr;
  character(len=32) :: arg
  character(len=8) :: date
  character(len=10) :: time
  character(len=5) :: zone
  logical :: do_time = .false.
  integer :: i
  integer :: ptlStart=0, ptlCount=10
  integer :: tStart=0, tCount=2; !! turnStart/Count
  integer :: bunchID=0;
  integer :: attrPos1=0, attrPos2=1, hasMore;

  i = 1;
  !do i = 1, 
  
  do while (i<=command_argument_count())
     call get_command_argument(i, arg)

     select case (arg)
     case ('-h', '--help')
        call print_help()
        stop
     case ('-t', '--turnStartCount')
        call hasEnoughInput(i,2, arg);
        call get_command_argument(i+1, arg)

        READ(arg, *) tStart
        call get_command_argument(i+2, arg)
        READ(arg, *) tCount
        i = i+2
        write (*,*) "turn start/count=", tStart, tCount
     case ('-p', '--ptlStartCount')
        call hasEnoughInput(i, 2, arg);
        call get_command_argument(i+1, arg)

        READ(arg, *) ptlStart
        call get_command_argument(i+2, arg)
        READ(arg, *) ptlCount
        i = i+2
        write (*,*) "ptl  start/count=", ptlStart, ptlCount

     case ('-a', '--attrs')
        call hasEnoughInput(i,2, arg);
        call get_command_argument(i+1, arg)

        READ(arg, *) attrPos1
        call get_command_argument(i+2, arg)
        READ(arg, *) attrPos2

        call checkAttrPos();
        write (*,*) "The attrs ids (0 based) ", attrPos1, attrPos2
        i = i+2
        
     case ('-b', '--bunch')
        call hasEnoughInput(i,1, arg);
        call get_command_argument(i+1, arg)

        READ(arg, *) bunchID
        write (*,*) "bunchID=", bunchID
        i = i+1;

     case default
        print '(a,a,/)', 'Unrecognized command-line option: ', arg
        call print_help()
        stop
     end select
     i = i+1;
  end do

  call MPI_INIT(ierr)

  call tunefoot_init(MPI_COMM_WORLD, ierr)
  if (ierr .eq. 0) then  
!---------------------------------------------------------------------
! NOTE: all of the following parameters are::: 0 based for C libraries
! 
!   format: 
!      function_name (start, count)
!---------------------------------------------------------------------

     call tunefoot_set_particleRange(ptlStart, ptlCount)
     call tunefoot_set_bunch(bunchID)
     call tunefoot_setTurnRange(tStart, tCount)

     !call tunefoot_writer_init(MPI_COMM_WORLD, ierr);
     call tunefoot_writer_init(MPI_COMM_NULL, attrPos1, attrPos2, ierr);
     
     hasMore = 1
     do while (hasMore .gt. 0)
        call tunefoot_run(attrPos1, attrPos2, hasMore) 
     enddo
  endif
  
  call tunefoot_close()
  call MPI_Finalize(ierr)
contains

  subroutine hasEnoughInput(i, c, arg)
    integer, intent(in) :: i;
    integer, intent(in) ::  c;
    character (len=32), intent(in) :: arg

    if (i+c > command_argument_count()) then
       write (*,*) "Error: expects ", c, " numbers following ", arg
       stop
    end if
  end subroutine hasEnoughInput
  
  subroutine checkAttrPos()
    if ( (attrPos1 >= 0) .and. (attrPos1 < 6) .and. (attrPos2 >= 0) .and. (attrPos2 < 6)) then 
       if (attrPos1 == attrPos2) then
          write (*,*) "Error! expecting two different attrs."
          stop
       endif
    else 
       write (*,*) "Error! attrs ids needs to be 0-5, representing: x, Px, y, Py, z, Pz"
       stop
    endif
  end subroutine checkAttrPos

  subroutine print_help()
    write (*,*) command_argument_count()
    print '(a)', 'usage: cmdline [OPTIONS]'
    print '(a)', ''
    print '(a)', ''
    print '(a)', 'cmdline options:'
    print '(a)', ''
    print '(a)', '  -h, --help        print usage information and exit'
    print '(a)', '  -a,     --attrs             followed by two attr poss  (0-6)'
    print '(a)', '  -b,     --bunch             followed by bunch ID (0 or 1)'
    print '(a)', '  -p,     --ptlStartCount     followed by startId & count for the particles (0 based)'
    print '(a)', '  -t,     --turnStartCount    followed by turnId & count for the FFT (0 based)'
    print '(a)', '                              note: turnCount needs to be powers of 2, e.g. 8, 16, etc '
    print '(a)', 'for example: '
    print '(a)', '   ./tunefoot -a 1 2 -b 0 -p 0 10 -t 3 16'
  end subroutine print_help

end program main
