! Pre-defined kind parameters and types. Based on Paolo Armani's NN
! digonalization code.

module feb_type
    implicit none

    integer, parameter :: standard_error_unit = 0
    integer, parameter :: standard_input_unit = 5
    integer, parameter :: standard_output_unit = 6
    INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
    INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
    INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
    INTEGER, PARAMETER :: SP = KIND(1.0)
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
    INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
    INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
    INTEGER, PARAMETER :: LGT = KIND(.true.)
    REAL(SP), PARAMETER :: PI_SP     = 3.141592653589793238462643383279502884197_sp
    REAL(SP), PARAMETER :: PIO2_SP   = 1.57079632679489661923132169163975144209858_sp
    REAL(SP), PARAMETER :: TWOPI_SP  = 6.283185307179586476925286766559005768394_sp
    REAL(SP), PARAMETER :: SQRT2_SP  = 1.41421356237309504880168872420969807856967_sp
    REAL(SP), PARAMETER :: EULER_SP  = 0.5772156649015328606065120900824024310422_sp
    REAL(DP), PARAMETER :: PI_DP     = 3.141592653589793238462643383279502884197_dp
    REAL(DP), PARAMETER :: PIO2_DP   = 1.57079632679489661923132169163975144209858_dp
    REAL(DP), PARAMETER :: TWOPI_DP  = 6.283185307179586476925286766559005768394_dp
    ! i, the unit imaginary number
    COMPLEX(DPC), PARAMETER :: IMGUNT_DP = (0.0_DP, 1.0_DP)
    COMPLEX(SPC), PARAMETER :: IMGUNT_SP = (0.0_SP, 1.0_SP)

end module feb_type
