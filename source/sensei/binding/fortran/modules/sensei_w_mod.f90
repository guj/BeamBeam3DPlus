module sensei_beam_w_mod
    !! the c pointer definitions
    use sensei_def_mod 

    implicit none

    interface writer_init
        module procedure sensei_w_init
        module procedure sensei_w_close
        module procedure sensei_w_setupPtl
    end interface


contains

    subroutine sensei_w_init(sensei, comm, name, ierr)
        type(SenseiBeamWriter), intent(out) :: sensei
        integer, intent(in) :: comm
        integer, intent(out) :: ierr
        character*(*), intent(in) :: name
        !call SenseiHandleWriter_createFile_f2c (sensei%f2c, comm, name, ierr)
        call SenseiHandleWriter_createFile_f2c(sensei%f2c, comm, TRIM(ADJUSTL(name))//char(0), ierr)
    end subroutine


    subroutine sensei_w_close(sensei, ierr)
        type(SenseiBeamWriter), intent(in) :: sensei
        integer, intent(out) :: ierr
        call SenseiHandleWriter_remove_f2c (sensei%f2c, ierr)
    end subroutine

    subroutine sensei_w_setupPtl(sensei, nbunchs, start, count, outRate, ierr)
        type(SenseiBeamWriter), intent(in) :: sensei
        integer, intent(out) :: ierr
        integer, intent(in) :: nbunchs
        integer*8, intent(in) :: start, count
        integer, intent(in):: outRate
        call SenseiHandleWriter_setupPtl_f2c(sensei%f2c, nbunchs, start, count, outRate,  ierr)
    end subroutine sensei_w_setupPtl

    subroutine sensei_w_turn_start(sensei, turnID, ierr)
      type(SenseiBeamWriter), intent(in) :: sensei
      integer, intent(out) :: ierr
      integer, intent(in) :: turnID
      call SenseiHandleWriter_turnStart_f2c(sensei%f2c, turnID, ierr)
    end subroutine sensei_w_turn_start

    subroutine sensei_w_turn_end(sensei, ierr)
      type(SenseiBeamWriter), intent(in) :: sensei
      integer, intent(out) :: ierr
      call SenseiHandleWriter_turnEnd_f2c(sensei%f2c, ierr)
    end subroutine sensei_w_turn_end

    subroutine sensei_writeAttr(sensei, bunchID, attrName, data, ierr)
      type(SenseiBeamWriter), intent(in) :: sensei
      integer, intent(in) :: bunchID
      character*(*), intent(in) :: attrName
      integer, intent(out) :: ierr
      real(kind=8), dimension(:), intent(in):: data
      call SenseiHandleWriter_attr_f2c(sensei%f2c, bunchID, TRIM(ADJUSTL(attrName))//char(0), data, ierr)
    end subroutine sensei_writeAttr
     
end module sensei_beam_w_mod


