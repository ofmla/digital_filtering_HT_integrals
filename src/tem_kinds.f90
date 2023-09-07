!*******************************************************************************
!>
!  Numeric kinds for TEM.

module tem_kinds

    use iso_fortran_env, only: real32, real64  ! precision kinds

    private

    integer,parameter,public :: dp = real64  
    integer,parameter,public :: sp = real32  

end module tem_kinds
!*******************************************************************************

