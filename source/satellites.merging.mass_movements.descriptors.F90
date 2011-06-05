!% Contains a module which defines descriptors for satellite merger mass movements.

module Satellite_Merging_Mass_Movements_Descriptors
  !% Defines descriptors for satellite merger mass movements.
  public

  integer, parameter :: doesNotMove    =0
  integer, parameter :: movesToDisk    =1
  integer, parameter :: movesToSpheroid=2
  
  ! Stored mass movement descriptors for the current merging event.
  integer            :: thisMergerGasMovesTo,thisMergerStarsMoveTo,thisHostGasMovesTo,thisHostStarsMoveTo
  !$omp threadprivate(thisMergerGasMovesTo,thisMergerStarsMoveTo,thisHostGasMovesTo,thisHostStarsMoveTo)

end module Satellite_Merging_Mass_Movements_Descriptors
