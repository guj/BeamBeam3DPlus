module sensei_beam_r_mod
    !! the c pointer definitions
    use sensei_def_mod 

    implicit none

    interface reader_init
        module procedure sensei_r_init
        module procedure sensei_r_close
        module procedure sensei_r_get_data
        module procedure sensei_r_get_params
    end interface


contains
    subroutine sensei_r_init(sensei, comm, name, ierr)
        type(SenseiBeamReader), intent(out) :: sensei
        integer, intent(in) :: comm
        integer, intent(out) :: ierr
        character*(*), intent(in) :: name
        call SenseiHandleReader_createFile_f2c(sensei%f2c, comm, TRIM(ADJUSTL(name))//char(0), ierr)
    end subroutine


    subroutine sensei_r_close(sensei, ierr)
        type(SenseiBeamReader), intent(in) :: sensei
        integer, intent(out) :: ierr
        call SenseiHandleReader_remove_f2c (sensei%f2c, ierr)
    end subroutine


    subroutine sensei_r_get_params(sensei, nTurns, nptls)
      type(SenseiBeamReader), intent(in) :: sensei
      integer*8, intent(out) :: nTurns
      integer*8, intent(out) :: nptls
      integer :: ierr
      call SenseiHandleReader_get_num_turns_f2c(sensei%f2c, nTurns, ierr)
      call SenseiHandleReader_get_num_particles_f2c(sensei%f2c, nptls, ierr)
    end subroutine sensei_r_get_params


    subroutine sensei_r_get_input(sensei, nTurns, nptls, bunchID)
      type(SenseiBeamReader), intent(in) :: sensei
      integer*8, intent(out) :: nTurns
      integer*8, intent(out) :: nptls
      integer*4, intent(out) :: bunchID
      integer :: ierr
      call SenseiHandleReader_get_input_f2c(sensei%f2c, nTurns, nptls, bunchID);
    end subroutine sensei_r_get_input


    subroutine sensei_r_get_actual(sensei, nTurns, nptls, ierr)
      type(SenseiBeamReader), intent(in) :: sensei
      integer*8, intent(out) :: nTurns
      integer*8, intent(out) :: nptls
      integer, intent(out) :: ierr
      call SenseiHandleReader_get_num_turns_f2c(sensei%f2c, nTurns, ierr)
      call SenseiHandleReader_get_num_particles_f2c(sensei%f2c, nptls, ierr)
    end subroutine sensei_r_get_actual

    subroutine sensei_r_get_data(sensei, data1, data2, ierr)
      type(SenseiBeamReader), intent(in) :: sensei
      real(kind=8), dimension (:), intent(in) :: data1, data2
      integer, intent(out) :: ierr
      
      call SenseiHandleReader_get_data_f2c(sensei%f2c, data1, data2, ierr)
    end subroutine sensei_r_get_Data
     
  end module sensei_beam_r_mod


