MODULE SENSEI_TUNEFOOT
  use mpi
  use sensei_beam_r_mod

  public
  integer :: m_rank, m_size

  CONTAINS
    !----------------------------------------------------------------------------!
    SUBROUTINE tunefoot_run_xml(comm, ierr)      
    !----------------------------------------------------------------------------!          
      integer, intent(in) :: comm 
      integer, intent(out) :: ierr
      
      integer*8 :: readInNumTurns, readInNptls;
      integer*8 :: numTurns, ptlPerTurn; !! actual turns & ptls read
      integer   :: nTurns_I4, ptl_I4
      integer :: ptls, i, bunchID;

      real*8, allocatable, dimension(:,:) :: data1, data2
      real*8, allocatable, dimension(:)   :: d1,    d2

      real*8, allocatable, dimension(:)   :: tune1, tune2

      type(SenseiBeamReader) :: senseiR  

      call sensei_r_init(senseiR, comm, char(0), ierr);
      call sensei_r_get_input(senseiR, readInNumTurns, readInNptls, bunchID);
      
      call isNotPowerOfTwo(readInNumTurns, ierr);
      if (ierr > 0) then
         write(*,*) "Error: num  of turns needs to be power of 2."
      else
         allocate(d1(readInNumTurns * readInNptls))
         allocate(d2(readInNumTurns * readInNptls))
         
         call sensei_r_get_data(senseiR, d1, d2, ierr);     
         
         call sensei_r_get_params(senseiR, numTurns, ptlPerTurn);

         if (numTurns <  readInNumTurns)  then
            write (*,*) "Not enough turns. "
            deallocate(d1)
            deallocate(d2)
            return;
         endif
         write(*,*) "matched"
         
         nTurns_I4 = numTurns;
         ptl_I4 = ptlPerTurn;
         !ptlPerTurn = sensei_r_get_ptls(senseiR);
         !numTurns = sensei_r_get_nturns(senseiR);
         
         allocate(data1(numTurns, ptlPerTurn));
         allocate(data2(numTurns, ptlPerTurn));
         
         do i=1, numTurns
            data1(i, :) = d1(1+(i-1)*ptlPerTurn:i*ptlPerTurn)
            data2(i, :) = d2(1+(i-1)*ptlPerTurn:i*ptlPerTurn)
         enddo

         deallocate(d1); 
         deallocate(d2); 

      
         write(*,*) "checking: turns/npts::", numTurns, ptlPerTurn
         
         ptls = ptlPerTurn
         if ((ptls > 0)) then
            !write(*,*) data1(1:6,1)
            write (*,*) "Here are values at first 6 TURNS for the FIRST particle on the two designed attributes:";
            write(*,'(6e15.5)') data1(1:6,1)
            write(*,'(6e15.5)') data2(1:6,1)
            
            allocate(tune1(ptls)); 
            allocate(tune2(ptls)); 
            
            call footprint(data1, data2,tune1,tune2, nTurns_I4, ptl_I4);      
            
            do i = 1, ptls
               write(4,*)tune1(i),tune2(i)
            enddo
            
            deallocate(tune1)
            deallocate(tune2)
            
         else
            if (m_rank == 0) then
               write(*,*) "no ptls retrieved. bye!"
            endif
         endif
         
         deallocate(data1)
         deallocate(data2)
      endif  ! isNotPowerOf2
      call sensei_r_close(senseiR, ierr)

    end SUBROUTINE tunefoot_run_xml


    recursive subroutine isNotPowerOfTwo(n, ierr)
      integer*8, intent(in) :: n
      integer, intent(out) :: ierr
      ierr = 0;
      if (n .eq. 2) then
         return
      endif
      if (mod(n,2) == 1) then
         ierr = 1
      else
         call isNotPowerOfTwo(n/2, ierr);
      endif
    end SUBROUTINE isNotPowerOfTwo
    
    !----------------------------------------------------------------------------!
    !find the working points (tunes) for npt points with M turns.
    subroutine footprint(xtr,ytr,tunex,tuney,M,npt) 
    !----------------------------------------------------------------------------!
      implicit none
      integer :: M,npt     
      real*8, dimension(M,npt), intent(in) :: xtr,ytr
      real*8, dimension(npt), intent(out)  :: tunex,tuney
      !real*8 data1(M)
      !real power(M/2)
      real,   allocatable, dimension(:) ::  power
      real*8, allocatable, dimension(:) ::  data1
      
      !double complex output(M)
      double complex, allocatable, dimension(:) ::output
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
      allocate(data1(M))
      allocate(power(M/2))
      allocate(output(M))

      xfrqmin = xfrqmin1
      xfrqmax = xfrqmax1

      do ilp = 1, Nloop
        do i = 1, M
          data1(i) = xtr(i,ilp) 
        enddo


        call realft(data1,M,isign)
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

      deallocate(output)
      deallocate(power)
      deallocate(data1)
      deallocate(btf)
      

      end subroutine footprint

    !----------------------------------------------------------------------------!
    SUBROUTINE realft(data,n,isign)
    !----------------------------------------------------------------------------!
      INTEGER, intent(in):: isign,n
      REAL*8, intent(inout) :: data(n)
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

END MODULE SENSEI_TUNEFOOT
