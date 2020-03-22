module sensei_def_mod
    implicit none
 
    !! C API handle
    type SenseiBeamWriter
        integer(kind=8):: f2c = 0_8
        logical :: valid = .false.
    end type

    type SenseiBeamReader
        integer(kind=8):: f2c = 0_8
        logical :: valid = .false.
    end type

end module sensei_def_mod

