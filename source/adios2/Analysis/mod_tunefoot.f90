MODULE ADIOS2_TUNEFOOT
  use adios2 
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
  type(adios2_variable)   :: var_handle   !sampled (ptFrac)

  integer :: rank, size;

  CONTAINS
    !----------------------------------------------------------------------------!
    SUBROUTINE tunefoot_init (comm, ierr)
    !----------------------------------------------------------------------------!          
      character*(20)                  :: fname
      integer, intent(out)            :: ierr
      integer, intent(in)             :: comm

      ! Init adios2
      call adios2_init_config (a2_handle, "adios2_config.xml", comm, &
           adios2_debug_mode_on, ierr)

      call adios2_declare_io (a2_io, a2_handle, "Tunefoot", ierr)

      fname = "beam3d.bp";
      call adios2_open (a2_Reader, a2_io, fname, adios2_mode_read, &
           comm, ierr)      
      write (*,*) ierr, a2_Reader


      call MPI_Comm_rank(comm, rank, ierr);
      call MPI_Comm_size(comm, size, ierr);
    end SUBROUTINE tunefoot_init

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
    SUBROUTINE tunefoot_run(pos1, pos2)      
    !----------------------------------------------------------------------------!
      integer :: step_status
      integer :: ierr;
      integer*8, dimension(3) :: selStart, selCount
      integer*8, allocatable, dimension(:) :: shape_in
      integer*8 :: current_step
      integer, intent(in):: pos1, pos2

      integer :: turnCounter = 1;

      double precision, allocatable, dimension(:,:) :: array1
      double precision, allocatable, dimension(:,:) :: array2;

      real*8, allocatable, dimension(:) :: tune1,tune2

      call tunefoot_display(pos1, pos2)

      selCount(1)= 1;
      selStart(2) = ParticleStart;   selCount(2) = ParticleSize;
      selStart(3) = BunchStart;      selCount(3) = BunchSize;

      allocate(array1(TurnSize, ParticleSize)); 
      allocate(array2(TurnSize, ParticleSize)); 

      allocate(tune1(ParticleSize)); 
      allocate(tune2(ParticleSize)); 
     
      do
         if (turnCounter > TurnSize) exit;

         !call adios2_begin_step_full(a2_Reader, adios2_step_mode_next_available, 0., &
         !     step_status, ierr)
         !!call adios2_begin_step(a2_Reader, ierr);
         ierr = 0;
         call adios2_begin_step(a2_Reader, adios2_step_mode_read, 0., step_status, ierr);
         if(step_status == adios2_step_status_end_of_stream) exit
         if (ierr .ne. 0) exit;

         call adios2_current_step(current_step, a2_Reader, ierr)

         if (current_step >= TurnStart) then
            !write(*,*) "curr step = ", current_step, "turncounter=", turnCounter, TurnStart, TurnSize, pos1, pos2
            call adios2_inquire_variable(var_handle,  a2_io, "particles", ierr)

            if( var_handle%valid .eqv. .true. ) then                     
               !call adios2_variable_shape(shape_in, var_handle%ndims, var_handle, ierr)               
               selStart(1) = pos1;
               call adios2_set_selection(var_handle, 3, selStart, selCount, ierr)                     
               
               call adios2_get(a2_Reader, var_handle, array1(turnCounter, 1:ParticleSize), adios2_mode_sync, ierr)
               !write(*,*) array1(turnCounter, 1:5)
               
               selStart(1) = pos2;
               call adios2_set_selection(var_handle, 3, selStart, selCount, ierr)
               call adios2_get(a2_Reader, var_handle, array2(turnCounter, 1:ParticleSize), adios2_mode_sync, ierr)
               !write(*,*) array2(turnCounter, 1:5)
            endif

            turnCounter = turnCounter + 1;
         endif
         call adios2_end_step(a2_Reader, ierr)
      enddo

      !write (*,*) turnCounter, TurnSize, ParticleSize

      if (turnCounter == TurnSize+1) then
         call footprint(array1, array2,tune1,tune2, TurnSize, ParticleSize);      
         do i = 1, ParticleSize
          write(4,*)tune1(i),tune2(i)
        enddo
     else if (turnCounter == 1) then
        write (*,*) "turn starts is too large. No actions taken"
        write (4,*) ""
     else
        write (*,*) "turn count is over the limit. no actions taken"
        write (4,*) ""
      endif

      deallocate(tune1);
      deallocate(tune2);

      deallocate(array1)
      deallocate(array2)

      
    end SUBROUTINE tunefoot_run
    

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
