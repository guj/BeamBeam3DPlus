program testSenseiBeam
  use mpi
  use sensei_beam_w_mod
  use sensei_beam_r_mod

  implicit none

  type(SenseiBeamWriter) :: senseiW
  integer rank, size
  integer :: ierr;
  !!character(len=:), allocatable :: var1_name

  integer*8:: total, count, start
  integer:: nbunchs, outRate, nTurns

  double precision, allocatable, dimension(:,:,:,:) :: data
  double precision, allocatable, dimension(:,:,:) :: dataPerTurn

  double precision, allocatable, dimension(:,:) :: temp2d, temp2dA, temp2dB
  double precision, allocatable, dimension(:) :: temp

  character*(20)                  :: attrName1, attrName2, attrName3, attrName4, attrName5, attrName6
  integer currTurn, attrID, bunchID
  integer*8 :: PID !! particle counter


  ! reader setup:
  type(SenseiBeamReader) :: senseiR
  integer*8 :: readInNumTurns, readInNptls;

  real*8, allocatable, dimension(:)   :: d1,    d2
  real*8, allocatable, dimension(:,:) :: data1, data2

  !
  ! Launch MPI
  !

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
  
  ! init
  total  = 20; count=total/size; start=rank*count;
  nbunchs = 3; outRate=5; nTurns = 20;
  
  attrName1 = "x";     attrName2 = "px" 
  attrName3 = "y";     attrName4 = "py"
  attrName5 = "z";     attrName6 = "pz"  

  allocate(data(6, count, nbunchs, nTurns));
   
    do attrID=1, 6
       do pID=1, count
          do bunchID=1, nbunchs
             do currTurn=1, nTurns
                if (attrID == 1) then
                   data(attrID, pID, bunchID, currTurn) = bunchID*1000 + pID*100 + rank*10 + currTurn*0.01;
                else
                   data(attrID, pID, bunchID, currTurn) = pID + 0.1*rank + bunchID*0.01 + currTurn*0.001 + attrID*0.0001;
                endif
             enddo
          enddo
       enddo
    enddo

  !!call sensei_w_init(senseiW, MPI_COMM_WORLD, "testFile", ierr)
  call sensei_w_init(senseiW, MPI_COMM_WORLD, char(0), ierr)



  call sensei_w_setupPtl(senseiW, nbunchs, start, count, outRate, ierr);
  do currTurn=1, nTurns
       allocate(dataPerTurn(6, count, nbunchs));
       dataPerTurn(:,:,:) = data(:,:,:, currTurn)

       call sensei_w_turn_start(senseiW, currTurn, ierr)

       allocate(temp(count));
       !! note: using "x" "px" directly instead of attrName* would result in 
       !! "xpx" in the attr name. So have to assign to attrName* first
       do bunchID=1,nbunchs
          temp(:) = dataPerTurn(1,:,bunchID)
          call sensei_writeAttr(senseiW, bunchID-1, attrName1, temp, ierr);
          temp(:) = dataPerTurn(2,:,bunchID)
          call sensei_writeAttr(senseiW, bunchID-1, attrName2, temp, ierr);
          temp(:) = dataPerTurn(3,:,bunchID)
          call sensei_writeAttr(senseiW, bunchID-1, attrName3, temp, ierr);
          temp(:) = dataPerTurn(4,:,bunchID)
          call sensei_writeAttr(senseiW, bunchID-1, attrName4, temp, ierr);
          temp(:) = dataPerTurn(5,:,bunchID)
          call sensei_writeAttr(senseiW, bunchID-1, attrName5, temp, ierr);
          temp(:) = dataPerTurn(6,:,bunchID)
          call sensei_writeAttr(senseiW, bunchID-1, attrName6, temp, ierr);
          !call sensei_writePtlData(bunchID-1,  temp2d)
       enddo
       call sensei_w_turn_end(senseiW, ierr)
       deallocate(temp);
       deallocate(dataPerTurn)
    enddo


    deallocate(data)
    call sensei_w_close(senseiW, ierr)


!
!
    
  call sensei_r_init(senseiR, MPI_COMM_WORLD, char(0), ierr)
  call sensei_r_get_params(senseiR, readInNumTurns, readInNptls);
  write (*,*) "num turns: ", readInNumTurns
  write (*,*) "num ptls: ", readInNptls

  allocate(d1(readInNumTurns * readInNptls))
  allocate(d2(readInNumTurns * readInNptls))

  call sensei_r_get_data(senseiR, d1, d2, ierr)

  write (*,'(4F10.4)') d1
  write (*,*)
  write (*,'(4F10.4)') d2
  
  !!verify d1, d2 and then put them to data1 data2 before passing to FFT
  deallocate(d1); 
  deallocate(d2); 

  call sensei_r_close(senseiR, ierr)
  call MPI_Finalize(ierr);
end program testSenseiBeam

  
