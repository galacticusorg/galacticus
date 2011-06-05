!% Contains a module which stores properties of merger remnants related to their size.

module Satellite_Merging_Remnant_Sizes_Properties
  !% Stores properties of merger remnants related to their size.
  public
  
  ! Value indicating that there was no change in the remnant spheroid.
  double precision, parameter :: remnantNoChangeValue=-1.0d0

  double precision :: remnantRadius,remnantCircularVelocity,remnantSpecificAngularMomentum
  !$omp threadprivate(remnantRadius,remnantCircularVelocity,remnantSpecificAngularMomentum)
  
end module Satellite_Merging_Remnant_Sizes_Properties
