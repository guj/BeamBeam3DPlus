MODULE ADIOS2_UTIL
  !  this will cause conflict with Commun.f90 class on MPI definition
  !  use mpi
    use adios2
    public

    type(adios2_adios)      :: a2_handle

    type(adios2_io)         :: a2_io_file
    type(adios2_engine)     :: a2_engine_file

    type(adios2_io)         :: a2_io
    type(adios2_engine)     :: a2_engine

    type(adios2_io)         :: a2_phase_io
    type(adios2_engine)     :: a2_phaseout_engine

    type(adios2_variable)   :: var_handle   !sampled (ptFrac)
    type(adios2_variable)   :: phvar_handle !checkpoint, phaseout 
    type(adios2_attribute)  :: attr_NumTurns
    type(adios2_attribute)  :: attr_TurnGap
    type(adios2_attribute)  :: attr_Turn

    logical                 :: adios2_initialized
    integer*8               :: pNbunch
    integer*8, dimension(3) :: vTotal, vStart, vCount    ! for var_handle
    integer*8, dimension(3) :: phTotal, phStart, phCount ! for phvar_handle

    CONTAINS

    !----------------------------------------------------------------------------!
    SUBROUTINE adios2_io_init (comm, ierr)
    !----------------------------------------------------------------------------!
        implicit none

        integer, intent(out)            :: ierr
        integer, intent(in)             :: comm
        character*(20)                  :: fname
        character*(20)                  :: currVarName
        integer                         :: i
        character(len=:), allocatable :: engineType

        ierr = 0

        ! Init adios2
        !call adios2_init_config (a2_handle, "adios2_config.xml", comm, &
        call adios2_init (a2_handle, "adios2_config.xml", comm, &
             adios2_debug_mode_on, ierr)
        
        ! Init IO object
        call adios2_declare_io (a2_io, a2_handle, "SimulationOutput", ierr)

        ! if using SST, then provide a file backup
        call adios2_io_engine_type(engineType, a2_io, ierr)

        if ((engineType == 'BPFile') .or. (engineType == 'HDF5')) then
           ! no need to look for file backup
        else
           !call adios2_declare_io(a2_io_file, a2_handle, "SimulationOutputFile", ierr) 
           !if (a2_io_file %valid .eqv. .true.) then 
           !   call adios2_open (a2_engine_file, a2_io_file, "beam3d.bp", adios2_mode_write, comm, ierr)
           !endif
        endif

        call adios2_declare_io (a2_phase_io, a2_handle, "PhaseOutput", ierr)

        fname = "beam3d.bp"
        ! Open file
        call adios2_open (a2_engine, a2_io, fname, adios2_mode_write, &
             comm, ierr)
        call adios2_open (a2_phaseout_engine, a2_phase_io, "phaseout", adios2_mode_write, &
             comm, ierr)
      end subroutine adios2_io_init



    !----------------------------------------------------------------------------!
    SUBROUTINE adios2_turn_start(ierr)
    !----------------------------------------------------------------------------!
      implicit none
      integer, intent(out)            :: ierr
      call adios2_begin_step (a2_engine, ierr)
      if (a2_io_file %valid .eqv. .true.) then 
         call adios2_begin_step (a2_engine_file, ierr)
      endif
    end SUBROUTINE adios2_turn_start


    !----------------------------------------------------------------------------!
    SUBROUTINE adios2_turn_end(ierr)
    !----------------------------------------------------------------------------!
      implicit none
      integer, intent(out)            :: ierr
      call adios2_end_step (a2_engine, ierr)
      if (a2_io_file %valid .eqv. .true.) then 
         call adios2_end_step(a2_engine_file, ierr)
      endif
    end SUBROUTINE adios2_turn_end



    SUBROUTINE adios2_setup_ptl(nbunchs, global, start, count, outRate, nturns, ierr)
        implicit none

        integer, intent(out)            :: ierr
        integer, intent(in)             :: nbunchs
        integer, intent(in)             :: global, start, count, outRate, nturns
        integer                         :: i

        ierr = 0

        ! Define variables
        vStart(1) = 0; vStart(2) = start;  vStart(3) = 0;
        vCount(1) = 6; vCount(2) = count;  vCount(3) = 1;
        vTotal(1) = 6; vTotal(2) = global; vTotal(3) = nbunchs;
        pNbunch = nbunchs

        ! not constant dim. (fixed shape/start/count)
        call adios2_define_variable(var_handle, a2_io, "particles", adios2_type_dp, 3, &
             vTotal, vStart, vCount, .false., ierr)
        
        ! var attr
        call adios2_define_attribute(attr_TurnGap, a2_io, 'ptRate', &
             outRate, var_handle%name, '/', ierr)

        ! global attr
        call adios2_define_attribute(attr_NumTurns, a2_io, 'nturns', &
             nturns, ierr)

    end subroutine adios2_setup_ptl


    !----------------------------------------------------------------------------!
    SUBROUTINE adios2_writePtlData(currBunch, iturn, data)
    !----------------------------------------------------------------------------!
        !--------------------------------------------------------------------
        ! ERROR INDICATORS AND WARNINGS
        !
        !--------------------------------------------------------------------
        ! MPI library
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        ! Declare variables
        integer, intent(in)                 :: currBunch, iturn
        integer                             :: ierr
        !real(kind=8), intent(in)            :: xcoords(:), ycoords(:), zcoords(:)
        double precision, dimension(vCount(1),vCount(2)), intent(in) :: data
        ! create character array with full filename
        ! write out using 2DECOMP&FFT MPI-IO routines

        ierr = 0

        vStart(3) = currBunch-1;
        !write (*,*) vStart(1), vStart(2), vStart(3)
        
        call adios2_set_selection(var_handle, 3, vStart, vCount, ierr)
         
        ! attributes are not step aware
        !call adios2_define_attribute(attr_Turn, a2_io, 'iturn', &
        !     iturn, var_handle%name,  '/', ierr)

        call adios2_put (a2_engine, var_handle, data, adios2_mode_sync, ierr)
        if (a2_io_file %valid .eqv. .true.) then 
           call adios2_put(a2_engine_file, var_handle, data, adios2_mode_sync, ierr)
        endif

    END SUBROUTINE adios2_writePtlData


    !----------------------------------------------------------------------------!
    SUBROUTINE adios2_setup_phaseout(nbunchs, global, start, count, ierr)
    !----------------------------------------------------------------------------!
        implicit none

        integer, intent(out)            :: ierr
        integer, intent(in)             :: nbunchs
        integer, intent(in)             :: global, start, count

        ierr = 0

        ! Define variables
        
        phStart(1) = 0; phStart(2) = start;  phStart(3) = 0;
        phCount(1) = 6; phCount(2) = count;  phCount(3) = 1;
        phTotal(1) = 6; phTotal(2) = global; phTotal(3) = nbunchs;
        pNbunch = nbunchs
        ! not constant dim. (fixed shape/start/count)
        call adios2_define_variable(phvar_handle, a2_phase_io, "phaseOutParticles", adios2_type_dp, 3, &
             phTotal, phStart, phCount, .false., ierr)
    end subroutine adios2_setup_phaseout


    !----------------------------------------------------------------------------!
    SUBROUTINE adios2_phaseout(data, currbunch, itic, iseedinit, close2g)
    !----------------------------------------------------------------------------!
      implicit none
      integer, intent(in) ::  itic, currbunch
      integer*8, intent(in) :: iseedinit
      double precision, dimension(:,:) :: data
      double precision, dimension(2) :: close2g
      integer  :: ierr
      ! no timesteps, only one copy of checkpoint
      phStart(3) = currBunch-1;
      call adios2_set_selection(phvar_handle, 3, phStart, phCount, ierr)             
      call adios2_put (a2_phaseout_engine, phvar_handle, data, adios2_mode_sync, ierr)

      !call adios2_phaseoutPtl(idbunch, iturn, itic, iseedinit, close2g, Bpts);
      !write(*,*)Nplocal,iturn,itic,iseedinit,close2g(1),close2g(2)
      !write(9)Bpts(1:6,1:Nplocal)

    end SUBROUTINE adios2_phaseout


    !----------------------------------------------------------------------------!
    SUBROUTINE adios2_io_finalize(ierr)
    !----------------------------------------------------------------------------!
        implicit none

        integer, intent(out) :: ierr

        ierr = 0
        call adios2_close    (a2_phaseout_engine, ierr)
        call adios2_close    (a2_engine, ierr)

        if (a2_io_file %valid .eqv. .true.) then
           call adios2_close    (a2_engine_file, ierr)
        endif
        !call adios2_close    (a2_engine, ierr)

        call adios2_finalize (a2_handle, ierr)

    END SUBROUTINE adios2_io_finalize

END MODULE ADIOS2_UTIL

