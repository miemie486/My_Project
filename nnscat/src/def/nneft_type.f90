! Bingwei Long
! nneft dialect for types
! changes to be made here when playing with different precisions

module nneft_type

  use feb_type
  implicit none

  integer, parameter  :: NER = DP
  integer, parameter  :: NEC = DPC
  ! integer, parameter  :: NER = SP
  ! integer, parameter  :: NEC = SPC


  real(NER), parameter :: PI_NE      = 3.141592653589793238462643383279502884197_NER
  real(NER), parameter :: PIO2_NE    = 1.57079632679489661923132169163975144209858_NER
! OEZ = 180.0
  real(NER), parameter :: OEZOPI_NE       = 180.0_NER / PI_NE
  real(NER), parameter :: PIOOEZ_NE       = PI_NE / 180.0_NER
  real(NER), parameter :: NINETYOPI_NE    = 90.0_NER / PI_NE
  real(NER), parameter :: PIONINETY_NE    = PI_NE / 90.0_NER

  complex(NEC), parameter :: IMGUNT_NE = (0.0_NER, 1.0_NER)
end module nneft_type
