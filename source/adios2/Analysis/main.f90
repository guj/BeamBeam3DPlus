program main
  use adios2_tunefoot
  use mpi

  implicit none

  !!character(len=*), parameter :: version = '1.0'
  integer :: ierr, tmp;

  character(len=8) :: date
  character(len=10) :: time
  character(len=5) :: zone
  !logical :: do_time = .false.
  integer :: ptlStart=0, ptlCount=10000000
  integer :: tStart=0, tCount=32768; !! turnStart/Count
  integer :: bunchID=0;
  integer :: attrPos1=0, attrPos2=2, attrPos3=4
  integer, allocatable, dimension(:) :: tuneAttrPos; 
  integer :: tuneAttrPosSize=0;
  integer :: hasMore
  integer :: rank_id;

  call MPI_INIT(ierr)  
  call MPI_Comm_rank(MPI_COMM_WORLD, rank_id, ierr);

  call read_command_line(ierr)

  if (ierr .ne. 0) then
     call MPI_Finalize(ierr)
     return;
  endif

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

     !call tunefoot_writer_init(MPI_COMM_NULL, tuneAttrPos, tuneAttrPosSize, ierr);
     call tunefoot_writer_init(MPI_COMM_WORLD, tuneAttrPos, tuneAttrPosSize, ierr);
     
     hasMore = 1
     do while (hasMore .gt. 0)
        call tunefoot_run(tuneAttrPos, tuneAttrPosSize, hasMore) 
     enddo
  endif
  
  call tunefoot_close()
  call MPI_Finalize(ierr)
contains

!--------------------------------------------
  subroutine hasEnoughInput(i, c, arg, yesNo)
!--------------------------------------------
    integer, intent(in) :: i;
    integer, intent(in) ::  c;
    character (len=32), intent(in) :: arg
    logical, intent(out) :: yesNo

    yesNo = .true.;
    if (i+c > command_argument_count()) then
       if (rank_id .eq. 0) write (*,*) "Error: expects ", c, " numbers following ", arg
       yesNo = .false.
    end if
  end subroutine hasEnoughInput

!--------------------------------------------  
  subroutine checkAttrPos()
!--------------------------------------------
    if ( (attrPos1 >= 0) .and. (attrPos1 < 6) .and. (attrPos2 >= 0) .and. (attrPos2 < 6)) then 
       if (attrPos1 == attrPos2) then
          if (rank_id .eq. 0) write (*,*) "Error! expecting two different attrs."
          stop
       endif
    else 
       if (rank_id .eq. 0) write (*,*) "Error! attrs ids needs to be 0-5, representing: x, Px, y, Py, z, Pz"
       stop
    endif
  end subroutine checkAttrPos

!--------------------------------------------
  subroutine print_help()
!--------------------------------------------
    if (rank_id .eq. 0) then  
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
    endif
  end subroutine print_help

!--------------------------------------------
  subroutine  read_input_file(ierr)
    implicit none
    integer, intent(out) :: ierr;
    logical  :: file_exists
    character (len = 40) :: file_name  = "tunefoot.in"
    integer, dimension(6) :: attrPos
    integer :: counter, i;
    ierr = 0;

    INQUIRE(FILE=file_name, EXIST=file_exists)

    if  (.not. file_exists) then 
       if (rank_id .eq. 0) write(*,*)  "Expecting input file: ", file_name
       ierr = -1
       return;
    endif

    if (rank_id .eq. 0) write(*,*) "reading from input file:", file_name, file_exists       
    
    open(unit=13,file=file_name,status='old')

    read(13,*) bunchID
    read(13,*) ptlStart, ptlCount
    read(13,*) tStart, tCount
    read(13,*) attrPos
    close(13)

    counter = 0;
    do i=1,6
       if (attrPos(i) > 0) counter = counter +1
    enddo

    if (counter < 2) then 
       if (rank_id .eq. 0) write(*,*) "Please use >= 2 attributes for tunefoot study."
       ierr = -1
    endif
    
    tuneAttrPosSize = counter
    allocate(tuneAttrPos(counter));
    counter = 1;
    do i=1,6
       if (attrPos(i) > 0) then 
          tuneAttrPos(counter)=i-1;
          counter = counter+1;
       endif
    enddo
    !write(*,*) tuneAttrPos
  end subroutine read_input_file
!--------------------------------------------

!--------------------------------------------
  subroutine read_command_line(ierr)
!--------------------------------------------
    implicit none
    integer i;
    integer, intent(out) :: ierr
    character(len=32) :: arg
    logical :: yesNo
    
    ierr = 0; 
    i = 1

    if (command_argument_count() == 0) then
       call read_input_file(ierr)
       return
    endif

    do while (i<=command_argument_count())
     call get_command_argument(i, arg)

     select case (arg)
     case ('-h', '--help')
        call print_help()        
        ierr = -1;
        exit
     case ('-t', '--turnStartCount')        
        call hasEnoughInput(i,2, arg, yesNo);
        if (.not. yesNo) then 
           ierr = -1
           exit
        endif
        call get_command_argument(i+1, arg)

        READ(arg, *) tStart
        call get_command_argument(i+2, arg)
        READ(arg, *) tCount
        i = i+2
        if (rank_id .eq. 0) write (*,*) "turn start/count=", tStart, tCount
     case ('-p', '--ptlStartCount')
        call hasEnoughInput(i, 2, arg, yesNo);
        if (.not. yesNo) then 
           ierr = -1
           exit
        endif
        call get_command_argument(i+1, arg)

        READ(arg, *) ptlStart
        call get_command_argument(i+2, arg)
        READ(arg, *) ptlCount
        i = i+2
        if (rank_id .eq. 0) write (*,*) "ptl  start/count=", ptlStart, ptlCount

     case ('-a', '--attrs')
        call hasEnoughInput(i,2, arg, yesNo);
        if (.not. yesNo) then
           ierr = -1
           exit
        endif
        call get_command_argument(i+1, arg)

        READ(arg, *) attrPos1
        call get_command_argument(i+2, arg)
        READ(arg, *) attrPos2

        call checkAttrPos();
        if (rank_id .eq. 0) write (*,*) "The attrs ids (0 based) ", attrPos1, attrPos2
        i = i+2               
     case ('-b', '--bunch')
        call hasEnoughInput(i,1, arg, yesNo);
        if (.not. yesNo) then 
           ierr = -1
           exit
        endif
        call get_command_argument(i+1, arg)

        READ(arg, *) bunchID
        if (rank_id .eq. 0) write (*,*) "bunchID=", bunchID
        i = i+1;

     case default
        if (rank_id .eq. 0) print '(a,a,/)', 'Unrecognized command-line option: ', arg
        call print_help()
        ierr = -1
        exit
     end select
     i = i+1;
  end do
    
  end subroutine read_command_line

end program main
