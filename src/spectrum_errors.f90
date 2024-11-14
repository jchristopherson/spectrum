module spectrum_errors
    use iso_fortran_env
    implicit none

    integer(int32), parameter :: SPCTRM_MEMORY_ERROR = 10000
    integer(int32), parameter :: SPCTRM_INVALID_INPUT_ERROR = 10001
    integer(int32), parameter :: SPCTRM_ARRAY_SIZE_MISMATCH_ERROR = 10002
    integer(int32), parameter :: SPCTRM_ARRAY_SIZE_ERROR = 10003
    integer(int32), parameter :: SPCTRM_SINGULAR_MATRIX_ERROR = 10004
end module