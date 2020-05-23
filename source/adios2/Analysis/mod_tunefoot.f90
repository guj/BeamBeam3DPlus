MODULE ADIOS2_TUNEFOOT
  use mpi
  use adios2 
  implicit  none
  public
  integer :: ParticleStart=0; !nbin
  integer :: ParticleSize=10; 

  integer :: TurnStart=0; !nline
  integer :: TurnSize=5; 
  
  integer :: BunchStart=0; !nfile
  integer :: BunchSize=1; 
  
  type(adios2_adios)      :: a2_handle
  type(adios2_io)         :: a2_io
  type(adios2_engine)     :: a2_Reader
  type(adios2_variable)   :: var_handle   

  type(adios2_io)         :: a2_io_tunefoot
  type(adios2_variable), allocatable, dimension(:)   :: tune_handles
  type(adios2_variable)   :: tune1_handle 
  type(adios2_variable)   :: tune2_handle 
  type(adios2_engine)     :: a2_Writer
  integer :: a2_store_tunefoot = 0

  integer :: rank, size;

  CONTAINS
    !----------------------------------------------------------------------------!
    SUBROUTINE tunefoot_init (comm, ierr)
    !----------------------------------------------------------------------------!          
      INCLUDE 'mpif.h'
      character*(20)                  :: fname
      integer, intent(out)            :: ierr
      integer, intent(in)             :: comm

      double precision :: start, end

      call MPI_Comm_rank(comm, rank, ierr);
      call MPI_Comm_size(comm, size, ierr);

      start  =  MPI_WTIME()
      ! Init adios2
      ! older adios version: call adios2_init_config (a2_handle, "adios2_config.xml", comm, &
      call adios2_init (a2_handle, "adios2_config.xml", comm, &
           adios2_debug_mode_on, ierr)

      call adios2_declare_io (a2_io, a2_handle, "SimulationOutput", ierr)

      fname = "beam3d.bp";
      call adios2_open (a2_Reader, a2_io, fname, adios2_mode_read, &
           comm, ierr)     

      end  = MPI_WTime()

    end SUBROUTINE tunefoot_init


    !----------------------------------------------------------------------------!
    SUBROUTINE  tunefoot_close()
    !----------------------------------------------------------------------------!
      INCLUDE 'mpif.h'

      integer :: ierr
      double precision :: start, end

      start =  MPI_WTime()
      call adios2_close(a2_Reader, ierr)
      if (rank == 0) call adios2_close(a2_Writer, ierr);
      call adios2_finalize(a2_handle, ierr)

      end =  MPI_WTime()
      if (rank  == 0) write (*, *) "tunefoot_close consumed:", end-start
    end SUBROUTINE tunefoot_close

    !----------------------------------------------------------------------------!
    SUBROUTINE tunefoot_set_particleRange(pStart, pSize)
    !----------------------------------------------------------------------------!
      integer, intent(in)::pStart, pSize

      ParticleStart = pStart
      ParticleSize = pSize

    end SUBROUTINE tunefoot_set_particleRange
    
    !----------------------------------------------------------------------------!
    SUBROUTINE tunefoot_set_bunch(b)
    !----------------------------------------------------------------------------!
      integer, intent(in)  ::  b

      BunchStart = b;
      BunchSize = 1
    end SUBROUTINE tunefoot_set_bunch

    !----------------------------------------------------------------------------!
    SUBROUTINE tunefoot_setTurnRange(tStart, tSize)
    !----------------------------------------------------------------------------!
      integer, intent(in)  ::  tStart, tSize

      TurnStart = tStart;
      TurnSize = tSize;
    end SUBROUTINE tunefoot_setTurnRange


    !----------------------------------------------------------------------------!
    SUBROUTINE tunefoot_display(pos1, pos2)
    !----------------------------------------------------------------------------!
      integer, intent(in):: pos1, pos2

      if (rank == 0) then
         write (*,*) "-----------------------------------"
         write (*,*) "Turns: ",TurnStart, TurnSize
         write (*,*) "Ptls:  ",ParticleStart, ParticleSize
         write (*,*) "Bunch: ",BunchStart
         write (*,*) "Attr:  ",pos1, pos2
      endif
    end SUBROUTINE tunefoot_display


    !----------------------------------------------------------------------------!
    SUBROUTINE tunefoot_display_all(attrPos, size)
    !----------------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: size
      integer, dimension(size), intent(in):: attrPos

      if (rank == 0) then
         write (*,*) "-----------------------------------"
         write (*,*) "Turns: ",TurnStart, TurnSize
         write (*,*) "Ptls:  ",ParticleStart, ParticleSize
         write (*,*) "Bunch: ",BunchStart
         write (*,*) "Attr:  ", attrPos
      endif
    end SUBROUTINE tunefoot_display_all

    !----------------------------------------------------------------------------!
    SUBROUTINE GetAttrName(pos,  name)
      integer,  intent(in)  :: pos
      character(2), intent(out) :: name
      write (*,*) "GetAttrName at: ",pos
      if (pos == 0)  then
         name = trim("x")
      else if (pos == 1) then
         name = "Px"
      else if  (pos == 2)  then
         name = "y"
      else if (pos == 3)  then 
         name  = "Py"
      else if (pos == 4) then  
         name  = "z"
      else
         name  = "Pz"
      endif
      write(*,*) name
    end SUBROUTINE GetAttrName
    !----------------------------------------------------------------------------!

    !----------------------------------------------------------------------------!
    SUBROUTINE  tunefoot_writer_init(comm, tuneAttrPos,  size, ierr)
    !----------------------------------------------------------------------------!
      implicit none
      character*(20)                  :: fname
      character*(2)                   :: attrName
      integer, intent(out)            :: ierr
      integer, intent(in)             :: comm, size
      integer, dimension(size),intent(in) :: tuneAttrPos
      integer*8, dimension(1)         :: vStart, vCount; !! one core, no need vTotal
      integer                         :: i
      Type(adios2_attribute)          :: attribute;


      if (rank > 0) return
      ! now write out tunefoot
      a2_store_tunefoot = 0;
      call adios2_declare_io (a2_io_tunefoot, a2_handle, "TuneFoot", ierr)

      if (ierr .ne. 0) return;

      fname = "tunefoot.bp"
      call adios2_open (a2_Writer, a2_io_tunefoot, fname, adios2_mode_write, &
           comm, ierr)      

      if (ierr .ne. 0) return;
      ! Define variables                                                                                                                         
      vStart(1) = 0;
      vCount(1) = ParticleSize; 
      !pNbunch = nbunchs

      allocate(tune_handles(size))

      do i=1,size
         call GetAttrName(tuneAttrPos(i), attrName)
         call adios2_define_variable(tune_handles(i), a2_io_tunefoot, trim(attrName), adios2_type_dp, 1, &
              vCount, vStart, vCount, .false., ierr)
         
         if (ierr .ne. 0) return;
      enddo

      call  adios2_define_attribute(attribute, a2_io_tunefoot,  "ParticleStart",  ParticleStart,    ierr);
      call  adios2_define_attribute(attribute, a2_io_tunefoot,  "TurnStart",  TurnStart,    ierr);
      call  adios2_define_attribute(attribute, a2_io_tunefoot,  "TurnSize",  TurnSize,    ierr);
      call  adios2_define_attribute(attribute, a2_io_tunefoot,  "Bunch",  BunchStart,    ierr);


      a2_store_tunefoot = 1;
    end SUBROUTINE TUNEFOOT_WRITER_INIT


    !----------------------------------------------------------------------------!
    SUBROUTINE  tunefoot_writer_init1(comm, pos1,  pos2, ierr)
    !----------------------------------------------------------------------------!
      implicit none
      character*(20)                  :: fname
      character*(2)                   :: attrName
      integer, intent(out)            :: ierr
      integer, intent(in)             :: comm, pos1, pos2
      integer*8, dimension(1)         :: vStart, vCount; !! one core, no need vTotal

      Type(adios2_attribute)          :: attribute;

      if (rank > 0) return
      ! now write out tunefoot
      a2_store_tunefoot = 0;
      call adios2_declare_io (a2_io_tunefoot, a2_handle, "TuneFoot", ierr)

      if (ierr .ne. 0) return;

      fname = "tunefoot.bp"
      call adios2_open (a2_Writer, a2_io_tunefoot, fname, adios2_mode_write, &
           comm, ierr)      

      if (ierr .ne. 0) return;
      ! Define variables                                                                                                                         
      vStart(1) = 0;
      vCount(1) = ParticleSize; 
      !pNbunch = nbunchs

      call GetAttrName(pos1, attrName)
      call adios2_define_variable(tune1_handle, a2_io_tunefoot, trim(attrName), adios2_type_dp, 1, &
           vCount, vStart, vCount, .false., ierr)

      if (ierr .ne. 0) return;
      call GetAttrName(pos2, attrName)
      call adios2_define_variable(tune2_handle, a2_io_tunefoot, trim(attrName), adios2_type_dp, 1, &
           vCount, vStart, vCount, .false., ierr)
                                       
      if (ierr .ne. 0) return;

      call  adios2_define_attribute(attribute, a2_io_tunefoot,  "ParticleStart",  ParticleStart,    ierr);
      call  adios2_define_attribute(attribute, a2_io_tunefoot,  "TurnStart",  TurnStart,    ierr);
      call  adios2_define_attribute(attribute, a2_io_tunefoot,  "TurnSize",  TurnSize,    ierr);
      call  adios2_define_attribute(attribute, a2_io_tunefoot,  "Bunch",  BunchStart,    ierr);


      a2_store_tunefoot = 1;
    end SUBROUTINE TUNEFOOT_WRITER_INIT1


    !----------------------------------------------------------------------------!
    SUBROUTINE  tunefoot_writer_put(tune, numTunes)
    !----------------------------------------------------------------------------!
      INCLUDE 'mpif.h'
      integer :: ierr,i;
      integer, intent(in) :: numTunes
      real*8, dimension(numTunes,ParticleSize), intent(in) :: tune 
      double precision :: start, end

      start = MPI_WTIME()
      if (a2_store_tunefoot .eq. 0) return
      
      call adios2_begin_step(a2_Writer, ierr);
      do i=1,numTunes
         call adios2_put (a2_Writer, tune_handles(i), tune(i,:), adios2_mode_sync, ierr)
      enddo
      call adios2_end_step(a2_Writer, ierr);

      end = MPI_WTIME()
      if (rank == 0)  write (*, *) " ... tunefoot_writer_put consumed:", end-start

    end SUBROUTINE TUNEFOOT_WRITER_PUT

    !----------------------------------------------------------------------------!
    SUBROUTINE  tunefoot_writer_put1(tune1, tune2)
    !----------------------------------------------------------------------------!
      INCLUDE 'mpif.h'
      integer :: ierr;
      real*8, dimension(ParticleSize), intent(in) :: tune1,tune2
      double precision :: start, end

      start = MPI_WTIME()
      if (a2_store_tunefoot .eq. 0) return
      
      call adios2_begin_step(a2_Writer, ierr);
      call adios2_put (a2_Writer, tune1_handle, tune1, adios2_mode_sync, ierr)
      call adios2_put (a2_Writer, tune2_handle, tune2, adios2_mode_sync, ierr)
      call adios2_end_step(a2_Writer, ierr);

      end = MPI_WTIME()
      if (rank == 0)  write (*, *) " ... tunefoot_writer_put consumed:", end-start

    end SUBROUTINE TUNEFOOT_WRITER_PUT1


    !----------------------------------------------------------------------------!
    SUBROUTINE tunefoot_run1(pos1, pos2, hasMore)      
    !----------------------------------------------------------------------------!
      implicit none
      INCLUDE 'mpif.h'
      integer :: step_status
      integer :: ierr, tmp,  i, perRankPtl;
      integer*8, dimension(3) :: selStart, selCount
      integer*8, allocatable, dimension(:) :: shape_in
      integer*8 :: current_step
      integer, intent(in):: pos1, pos2
      integer, intent(out):: hasMore
      integer :: turnCounter = 1;
      
      double precision, allocatable, dimension(:,:) :: array1, localArray1
      double precision, allocatable, dimension(:,:) :: array2, localArray2;

      real*8, allocatable, dimension(:) :: tune1,tune2

      double  precision :: start, end

      start = MPI_WTIME()

      turnCounter = 1;
      hasMore = 0;

      call tunefoot_display(pos1, pos2)

      selCount(1)= 1;      
      !selStart(2) = ParticleStart;   selCount(2) = ParticleSize;
      perRankPtl = ParticleSize/size;
      selStart(2) = ParticleStart+perRankPtl*rank; selCount(2) = perRankPtl
      selStart(3) = BunchStart;      selCount(3) = BunchSize;

      allocate(localArray1(TurnSize, perRankPtl)); 
      allocate(localArray2(TurnSize, perRankPtl)); 

      allocate(tune1(ParticleSize)); 
      allocate(tune2(ParticleSize)); 
     
      do
         if (turnCounter > TurnSize) exit;

         ierr = 1;
         !! note that need to use -1 to make sure wait for response.
         call adios2_begin_step(a2_Reader, adios2_step_mode_read, -1.0, step_status, ierr);
         if (ierr .ne. 0) exit;

         if(step_status == adios2_step_status_end_of_stream) exit
         if(step_status == adios2_step_status_not_ready) exit         

         call adios2_current_step(current_step, a2_Reader, ierr)
         !write(*,*) "  curr=", current_step, TurnStart
         if (current_step >= TurnStart) then
            !write(*,*) "curr step = ", current_step, "turncounter=", turnCounter, TurnStart, TurnSize, pos1, pos2
            call adios2_inquire_variable(var_handle,  a2_io, "particles", ierr)

            if( var_handle%valid .eqv. .true. ) then                     
               !call adios2_variable_shape(shape_in, var_handle%ndims, var_handle, ierr)               
               selStart(1) = pos1;
               call adios2_set_selection(var_handle, 3, selStart, selCount, ierr)                     
               
               call adios2_get(a2_Reader, var_handle, localArray1(turnCounter, 1:perRankPtl), adios2_mode_sync, ierr)
               
               selStart(1) = pos2;
               call adios2_set_selection(var_handle, 3, selStart, selCount, ierr)
               call adios2_get(a2_Reader, var_handle, localArray2(turnCounter, 1:perRankPtl), adios2_mode_sync, ierr)
            endif

            turnCounter = turnCounter + 1;
         endif
         call adios2_end_step(a2_Reader, ierr)
      enddo

      end  = MPI_WTIME()
      if (rank == 0)  write (*, *) " ... tunefoot_get_data consumed:", end-start


      tmp  = TurnSize * perRankPtl
      if (rank == 0) then
         allocate(array1(TurnSize, ParticleSize));       
         allocate(array2(TurnSize, ParticleSize));       
      endif         

      call MPI_Gather(localArray1, tmp, MPI_DOUBLE_PRECISION, &  
                      array1, tmp, MPI_DOUBLE_PRECISION, &
                      0, MPI_COMM_WORLD, ierr)

      call MPI_Gather(localArray2, tmp, MPI_DOUBLE_PRECISION, &  
                      array2, tmp, MPI_DOUBLE_PRECISION, &
                      0, MPI_COMM_WORLD, ierr)

      call MPI_Barrier(MPI_COMM_WORLD, ierr);

      !write (*,*) turnCounter, TurnSize, ParticleSize

         if (turnCounter == TurnSize+1) then
            start = end;
            if (rank .eq. 0) then
               !!call footprint(array1, array2,tune1,tune2, TurnSize, ParticleSize);      
               call footprintOne(array1, tune1,TurnSize, ParticleSize);      
               call footprintOne(array2, tune2,TurnSize, ParticleSize);      
               end =  MPI_WTIME()
               write (*, *) " ... tunefoot compute consumed:", end-start
         
               do i = 1, ParticleSize
                  write(4,*)tune1(i),tune2(i)
               enddo
               !!  temp disable
               call tunefoot_writer_put1(tune1, tune2)
            endif ! rank == 0

            hasMore = 1; !! can try to run more
         else if (turnCounter == 1) then
            if (rank == 0) write (*,*) "turn starts is too large. No actions taken"
            if (rank == 0) write (4,*) ""
         else
            if (rank == 0) write (*,*) "turn count is over the limit. no actions taken"
            write (4,*) ""
         endif

      deallocate(tune1);
      deallocate(tune2);

      deallocate(localArray1)
      deallocate(localArray2)

      if (rank == 0) then
         deallocate(array1)
         deallocate(array2)
      endif
      
    end SUBROUTINE tunefoot_run1
    

    !----------------------------------------------------------------------------!
    SUBROUTINE tunefoot_run(attrPos, numTunes, hasMore)      
    !----------------------------------------------------------------------------!
      implicit none
      INCLUDE 'mpif.h'
      integer :: step_status
      integer :: ierr, tmp,  i, perRankPtl;
      integer*8, dimension(3) :: selStart, selCount
      integer*8, allocatable, dimension(:) :: shape_in
      integer*8 :: current_step
      integer, dimension(numTunes), intent(in):: attrPos
      integer, intent(in) :: numTunes
      integer, intent(out):: hasMore
      integer :: turnCounter = 1;
      
      double precision, allocatable, dimension(:,:,:) :: array, localArray
      double precision, allocatable, dimension(:,:) :: tmpArray


      real*8, allocatable, dimension(:,:) :: tune 

      double  precision :: start, end

      start = MPI_WTIME()

      turnCounter = 1;
      hasMore = 0;

      call tunefoot_display_all(attrPos, numTunes)

      selCount(1)= 1;      
      !selStart(2) = ParticleStart;   selCount(2) = ParticleSize;
      perRankPtl = ParticleSize/size;
      selStart(2) = ParticleStart+perRankPtl*rank; selCount(2) = perRankPtl
      selStart(3) = BunchStart;      selCount(3) = BunchSize;

      allocate(localArray(numTunes, TurnSize, perRankPtl)); 
      allocate(tune(numTunes, ParticleSize)); 
     
      do
         if (turnCounter > TurnSize) exit;

         ierr = 1;
         !! note that need to use -1 to make sure wait for response.
         call adios2_begin_step(a2_Reader, adios2_step_mode_read, -1.0, step_status, ierr);
         if (ierr .ne. 0) exit;

         if(step_status == adios2_step_status_end_of_stream) exit
         if(step_status == adios2_step_status_not_ready) exit         

         call adios2_current_step(current_step, a2_Reader, ierr)
         !write(*,*) "  curr=", current_step, TurnStart
         if (current_step >= TurnStart) then
            !write(*,*) "curr step = ", current_step, "turncounter=", turnCounter, TurnStart, TurnSize, pos1, pos2
            call adios2_inquire_variable(var_handle,  a2_io, "particles", ierr)

            if( var_handle%valid .eqv. .true. ) then                     
               !call adios2_variable_shape(shape_in, var_handle%ndims, var_handle, ierr)               
               do i=1,numTunes
                  selStart(1) = attrPos(i) !selStart(1) = pos1;
                  call adios2_set_selection(var_handle, 3, selStart, selCount, ierr)                     
                  call adios2_get(a2_Reader, var_handle, localArray(i,turnCounter, 1:perRankPtl), adios2_mode_sync, ierr)
               enddo
            endif

            turnCounter = turnCounter + 1;
         endif
         call adios2_end_step(a2_Reader, ierr)
      enddo

      end  = MPI_WTIME()
      if (rank == 0)  write (*, *) " ... tunefoot_get_data consumed:", end-start


      tmp  = TurnSize * perRankPtl
      if (rank == 0) then
         !!allocate(array(numTunes, TurnSize, ParticleSize));       
         allocate(tmpArray(TurnSize, ParticleSize));       
      endif         

      do i=1,numTunes
         call MPI_Gather(localArray(i,:,:), tmp, MPI_DOUBLE_PRECISION, &  
              tmpArray, tmp, MPI_DOUBLE_PRECISION, &
              0, MPI_COMM_WORLD, ierr)

         call MPI_Barrier(MPI_COMM_WORLD, ierr);
         if (turnCounter == TurnSize+1) then
            if (rank .eq. 0) then
               call footprintOne(tmpArray, tune(i,:), TurnSize, ParticleSize)
            endif
         endif         
      enddo

      if (rank == 0)  write (*, *) " ... tunefoot_cal consumed:", MPI_WTime() - end;
      call MPI_Barrier(MPI_COMM_WORLD, ierr);


      if (turnCounter == TurnSize+1) then
            if (rank .eq. 0) then
               deallocate(tmpArray)         
               call tunefoot_writer_put(tune, numTunes)
            endif ! rank == 0

            hasMore = 1; !! can try to run more
         else if (turnCounter == 1) then
            if (rank == 0) write (*,*) "turn starts is too large. No actions taken"
            if (rank == 0) write (4,*) ""
         else
            if (rank == 0) write (*,*) "turn count is over the limit. no actions taken"
            write (4,*) ""
         endif

      deallocate(tune);

      deallocate(localArray)

      if (rank == 0) then
         !!deallocate(array)
      endif
      
    end SUBROUTINE tunefoot_run
    

    !----------------------------------------------------------------------------!
    ! in : one array of particle attributes, one of x/px/y/py/z/pz
    ! out: tune for that input
    subroutine footprintOne(input, tune, M,npt) 
    !----------------------------------------------------------------------------!
      implicit none
      integer :: M,npt
      real*8, dimension(M,npt) :: input 
      real*8, dimension(npt) :: tune
      real*8 data1(M)
      real power(M/2)
      double complex output(M)
      integer :: i, zero, sign,nn,isign
      real :: unit,pi,deltat,a,b,x1,x2,x3,x4,x5,x6,x7
      real :: y1,y2,y3,y4,y5,y6,y7,scale,zz 
      real*8, allocatable, dimension(:,:) :: btf
      integer :: Nloop,ilp
      real*8 :: pamp1,xfrqmin,xfrqmax,xfrq,frqkick
      real*8 :: xfrqmin1,xfrqmin2,xfrqmax1,xfrqmax2
      
      unit = 1.0
      zero = 0
      sign = 1
      scale = 1.0
      nn = 1
      isign = 1
      pi = 2.0*asin(1.0d0)

      deltat = 1.0

      Nloop = npt
      !tune range of the working points
      xfrqmin1 = 0.02
      xfrqmax1 = 0.5
      xfrqmin2 = 0.02
      xfrqmax2 = 0.5

      allocate(btf(2,Nloop))

      xfrqmin = xfrqmin1
      xfrqmax = xfrqmax1
      do ilp = 1, Nloop

        do i = 1, M
          data1(i) = input(i,ilp) 
        enddo


        call realft(data1,M,isign)
        !write(*,*) "out=>", data1(1:10)
        output = 0.0
        output(1) = dcmplx(data1(1),0)
        output(M/2+1) = dcmplx(data1(2),0)
        do i = 2, M/2
           output(i) = cmplx(data1(2*i-1),data1(2*i)) 
           !output(i) = dcmplx(data1(2*i-1),data1(2*i)) 
        enddo
        power(1) = abs(output(1))*abs(output(1))/(M*M);
        do i = 2, M/2
           power(i) = abs(output(i))*abs(output(i))/M/M
        enddo

        if(ilp.eq.2) then
          do i = 1, M/2
            write(2,*) (i-1)/(M*deltat),"  ",power(i)
          enddo
        endif

        pamp1 = 0.0
        do i = 2, M/2-1
          xfrq = (i-1)*1.0d0/(M*deltat)
          if(xfrq.ge.xfrqmin .and. xfrq.le.xfrqmax) then
             if(pamp1.le.power(i)) then
               pamp1 = power(i)
               btf(1,ilp) = xfrq
               btf(2,ilp) = pamp1
             endif
          endif
        enddo
        tune(ilp) = btf(1,ilp)
      enddo


      deallocate(btf)

      end subroutine footprintOne

    !----------------------------------------------------------------------------!
    !find the working points (tunes) for npt points with M turns.
    subroutine footprint(xtr,ytr,tunex,tuney,M,npt) 
    !----------------------------------------------------------------------------!
      implicit none
      integer :: M,npt
      real*8, dimension(M,npt) :: xtr,ytr
      real*8, dimension(npt) :: tunex,tuney
      real*8 data1(M)
      real power(M/2)
      double complex output(M)
      integer :: i, zero, sign,nn,isign
      real :: unit,pi,deltat,a,b,x1,x2,x3,x4,x5,x6,x7
      real :: y1,y2,y3,y4,y5,y6,y7,scale,zz 
      real*8, allocatable, dimension(:,:) :: btf
      integer :: Nloop,ilp
      real*8 :: pamp1,xfrqmin,xfrqmax,xfrq,frqkick
      real*8 :: xfrqmin1,xfrqmin2,xfrqmax1,xfrqmax2
      
      unit = 1.0
      zero = 0
      sign = 1
      scale = 1.0
      nn = 1
      isign = 1
      pi = 2.0*asin(1.0d0)

      deltat = 1.0

      Nloop = npt
      !tune range of the working points
      xfrqmin1 = 0.02
      xfrqmax1 = 0.5
      xfrqmin2 = 0.02
      xfrqmax2 = 0.5

      allocate(btf(2,Nloop))

      xfrqmin = xfrqmin1
      xfrqmax = xfrqmax1
      do ilp = 1, Nloop

        do i = 1, M
          data1(i) = xtr(i,ilp) 
        enddo


        call realft(data1,M,isign)
        !write(*,*) "out=>", data1(1:10)
        output = 0.0
        output(1) = dcmplx(data1(1),0)
        output(M/2+1) = dcmplx(data1(2),0)
        do i = 2, M/2
           output(i) = cmplx(data1(2*i-1),data1(2*i)) 
           !output(i) = dcmplx(data1(2*i-1),data1(2*i)) 
        enddo
        power(1) = abs(output(1))*abs(output(1))/(M*M);
        do i = 2, M/2
           power(i) = abs(output(i))*abs(output(i))/M/M
        enddo

        if(ilp.eq.2) then
          do i = 1, M/2
            write(2,*) (i-1)/(M*deltat),"  ",power(i)
          enddo
        endif

        pamp1 = 0.0
        do i = 2, M/2-1
          xfrq = (i-1)*1.0d0/(M*deltat)
          if(xfrq.ge.xfrqmin .and. xfrq.le.xfrqmax) then
             if(pamp1.le.power(i)) then
               pamp1 = power(i)
               btf(1,ilp) = xfrq
               btf(2,ilp) = pamp1
             endif
          endif
        enddo
        tunex(ilp) = btf(1,ilp)
      enddo

      xfrqmin = xfrqmin2
      xfrqmax = xfrqmax2
      do ilp = 1, Nloop
        do i = 1, M
          data1(i) = ytr(i,ilp) 
        enddo

        call realft(data1,M,isign)

        output = 0.0
        output(1) = dcmplx(data1(1),0)
        output(M/2+1) = dcmplx(data1(2),0)
        do i = 2, M/2
           output(i) = cmplx(data1(2*i-1),data1(2*i)) 
        enddo
        power(1) = abs(output(1))*abs(output(1))/(M*M);
        do i = 2, M/2
           power(i) = abs(output(i))*abs(output(i))/M/M
        enddo

        if(ilp.eq.2) then
          do i = 1, M/2
            write(1,*) (i-1)/(M*deltat),"  ",power(i)
          enddo
        endif

        pamp1 = 0.0
        do i = 2, M/2-1
          xfrq = (i-1)*1.0d0/(M*deltat)
          if(xfrq.ge.xfrqmin .and. xfrq.le.xfrqmax) then
             if(pamp1.le.power(i)) then
               pamp1 = power(i)
               btf(1,ilp) = xfrq
               btf(2,ilp) = pamp1
             endif
          endif
        enddo
        tuney(ilp) = btf(1,ilp)
      enddo

      deallocate(btf)

      end subroutine footprint

    !----------------------------------------------------------------------------!
    SUBROUTINE realft(data,n,isign)
    !----------------------------------------------------------------------------!
      INTEGER isign,n
      REAL*8 data(n)
      INTEGER i,i1,i2,i3,i4,n2p3
      REAL*8 c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n/2)
      c1=0.5
      if (isign.eq.1) then
        c2=-0.5
        call four1(data,n/2,+1)
      else
        c2=0.5
        theta=-theta
      endif
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=n+3
      do 11 i=2,n/4
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=wr
        wis=wi
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
11    continue
      if (isign.eq.1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call four1(data,n/2,-1)
      endif
      return
      END subroutine realft

      !----------------------------------------------------------------------------!
      SUBROUTINE four1(data,nn,isign)
      !----------------------------------------------------------------------------!
      INTEGER isign,nn
      REAL*8 data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL*8 tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=wr*data(j)-wi*data(j+1)
            tempi=wr*data(j+1)+wi*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
    END subroutine four1

END MODULE ADIOS2_TUNEFOOT
