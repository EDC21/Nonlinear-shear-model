!DIR$ FREEFORM

!> Global variables for Fatigue CZM
module Fatigue_Globals

  integer, parameter :: maxCohElemNum = 3000000
  integer, parameter :: maxNodeNum = 3000000

  double precision :: GV4IP(maxCohElemNum, 32)
  double precision :: GV4Elem(maxCohElemNum, 8)
  double precision :: ElemDimens(maxCohElemNum, 4)

  double precision :: NodeCoords(maxNodeNum, 3)

  double precision :: NNW(500,3)
  double precision :: NNB(500)
  double precision :: NNC(500)
  double precision :: NNMinMax(4,2)
  double precision :: NNdSRGthP(5,5)

  integer :: CohElem(maxCohElemNum,9)
  integer :: IPEdgNbrElem(maxCohElemNum,9)
  integer :: IP1NbrElem(maxCohElemNum,12)
  integer :: IP2NbrElem(maxCohElemNum,12)
  integer :: IP3NbrElem(maxCohElemNum,12)
  integer :: IP4NbrElem(maxCohElemNum,12)
  integer :: CohElemNum
  integer :: NNN

end module Fatigue_Globals
!DIR$ NOFREEFORM
!DIR$ FIXEDFORMLINESIZE:132