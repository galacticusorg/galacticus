!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!% Contains a module which computes scaling factors to convert between model and observed cosmologies.

module Galacticus_Output_Analyses_Cosmology_Scalings
  !% Computes scaling factors to convert between model and observed cosmologies.
  private
  public :: Cosmology_Conversion_Factors
  
contains
  
  subroutine Cosmology_Conversion_Factors(redshift,cosmologyFunctionsModel,cosmologyFunctionsObserved,cosmologyScalingMass, &
       & cosmologyScalingMassFunction,cosmologyScalingSize,cosmologyConversionMass,cosmologyConversionMassFunction,cosmologyConversionSize)
    !% Compute conversion factors for mass, size, and volume to adjust model galaxy properties to the assumed cosmology of the
    !% observations.
    use Cosmology_Functions
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    double precision                         , intent(in   )           :: redshift
    class           (cosmologyFunctionsClass), intent(inout)           :: cosmologyFunctionsModel, cosmologyFunctionsObserved
    type            (varying_string         ), intent(in   ), optional :: cosmologyScalingMass   , cosmologyScalingMassFunction   , cosmologyScalingSize
    double precision                         , intent(  out), optional :: cosmologyConversionMass, cosmologyConversionMassFunction, cosmologyConversionSize
    double precision                                                   :: timeModel              , timeObserved

    ! Check optional arguments for consistency.
    if (present(cosmologyScalingMass        ).neqv.present(cosmologyConversionMass        )) call Galacticus_Error_Report('Cosmology_Conversion_Factors','"cosmologyScalingMass" and "cosmologyConversionMass" must be either both present or both not-present'                )
    if (present(cosmologyScalingSize        ).neqv.present(cosmologyConversionSize        )) call Galacticus_Error_Report('Cosmology_Conversion_Factors','"cosmologyScalingSize" and "cosmologyConversionSize" must be either both present or both not-present'                )
    if (present(cosmologyScalingMassFunction).neqv.present(cosmologyConversionMassFunction)) call Galacticus_Error_Report('Cosmology_Conversion_Factors','"cosmologyScalingMassFunction" and "cosmologyConversionMassFunction" must be either both present or both not-present')

    ! Test redshift.
    if (redshift <= 0.0d0) then
       ! At zero redshifts the conversions are all identity conversions.
       if (present(cosmologyConversionSize        )) cosmologyConversionSize        =1.0d0
       if (present(cosmologyConversionMass        )) cosmologyConversionMass        =1.0d0
       if (present(cosmologyConversionMassFunction)) cosmologyConversionMassFunction=1.0d0
    else
       ! Find the cosmic times in model and observed cosmologies.
       timeModel   =cosmologyFunctionsModel   %cosmicTime(cosmologyFunctionsModel   %expansionFactorFromRedshift(redshift))
       timeObserved=cosmologyFunctionsObserved%cosmicTime(cosmologyFunctionsObserved%expansionFactorFromRedshift(redshift))
       ! Compute scaling factor for masses.
       if (present(cosmologyScalingMass)) then
          select case (char(cosmologyScalingMass))
          case ('none'      )
             ! Nothing to do in this case, observational data does not scale with cosmology.
             cosmologyConversionMass=1.0d0
          case ('luminosity')
             ! Observational data scales as luminosity.
             cosmologyConversionMass=(cosmologyFunctionsObserved%distanceLuminosity(timeObserved)/cosmologyFunctionsModel%distanceLuminosity(timeModel))**2
          case default
             call Galacticus_Error_Report('Cosmology_Conversion_Factors','unrecognized cosmology scaling')
          end select
       end if
       ! Compute scaling factor for sizes.
       if (present(cosmologyScalingSize)) then
          select case (char(cosmologyScalingSize))
          case ('none'      )
             ! Nothing to do in this case, observational data does not scale with cosmology.
             cosmologyConversionSize=1.0d0
          case ('angular')
             ! Observational data scales as angular size.
             cosmologyConversionSize=cosmologyFunctionsObserved%distanceAngular(timeObserved)/cosmologyFunctionsModel%distanceAngular(timeModel)
          case default
             call Galacticus_Error_Report('Cosmology_Conversion_Factors','unrecognized cosmology scaling')
          end select
       end if
       ! Compute scaling factor for mass functions.
       if (present(cosmologyScalingMassFunction)) then
          select case (char(cosmologyScalingMassFunction))
          case ('none'      )
             ! Nothing to do in this case, observational data does not scale with cosmology.
             cosmologyConversionMassFunction=1.0d0
          case ('inverseComovingVolume')
             ! Observational data scales as the inverse of comoving volume determined from angular distance and redshift interval.
             cosmologyConversionMassFunction= (cosmologyFunctionsModel%distanceComoving      (timeModel)/cosmologyFunctionsObserved%distanceComoving      (timeObserved))**2 &
                  &                          /(cosmologyFunctionsModel%hubbleParameterEpochal(timeModel)/cosmologyFunctionsObserved%hubbleParameterEpochal(timeObserved))
          case default
             call Galacticus_Error_Report('Cosmology_Conversion_Factors','unrecognized cosmology scaling')
          end select
       end if
    end if
    return
  end subroutine Cosmology_Conversion_Factors
  
end module Galacticus_Output_Analyses_Cosmology_Scalings
