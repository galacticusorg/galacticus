!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!!{RST
Contains custom functions for the standard spheroid component.
!!}

function Node_Component_Spheroid_Standard_Mass_Distribution(self,componentType,massType,weightBy,weightIndex) result(massDistribution_)
  !!{RST
  Return the mass distribution for the standard spheroid component.
  !!}
  use :: Mass_Distributions                         , only : massDistributionClass                            , massDistributionMatches_
  use :: Node_Component_Spheroid_Standard_Data      , only : massDistributionStellar_                         , massDistributionGas_    , kinematicDistribution_      , scalerStellarPool      , &
       &                                                     scalerGasPool
  use :: Node_Components_Galactic_Mass_Distributions, only : Node_Component_Mass_Distribution_Scaled_Spherical
  use :: Galactic_Structure_Options                 , only : componentTypeSpheroid                            , massTypeStellar         , massTypeGaseous             , enumerationWeightByType, &
       &                                                     weightByMass                                     , weightByLuminosity      , enumerationComponentTypeType, enumerationMassTypeType
  implicit none
  class           (massDistributionClass          ), pointer                 :: massDistribution_
  class           (nodeComponentSpheroidStandard  ), intent(inout)           :: self
  type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
  type            (enumerationMassTypeType        ), intent(in   ), optional :: massType
  type            (enumerationWeightByType        ), intent(in   ), optional :: weightBy
  integer                                          , intent(in   ), optional :: weightIndex
  type            (stellarLuminosities            ), save                    :: luminosities
  !$omp threadprivate(luminosities)
  double precision                                                           :: massStellar      , massGas     , &
       &                                                                        radiusScale
  logical                                                                    :: includeGas       , includeStars
  !![
  <optionalArgument name="weightBy" defaultsTo="weightByMass" />
  !!]
  
  ! Determine which components of the spheroid to include.
  includeGas  =massDistributionMatches_(componentTypeSpheroid,massTypeGaseous,componentType,massType) .and.  weightBy_ == weightByMass
  includeStars=massDistributionMatches_(componentTypeSpheroid,massTypeStellar,componentType,massType) .and. (weightBy_ == weightByMass .or. weightBy_ == weightByLuminosity)  
  ! Get properties of the mass distribution and ensure they are physical.
  if      (weightBy_ == weightByMass      ) then
     massStellar       =  max (0.0d0,self        %massStellar        (           ))
     massGas           =  max (0.0d0,self        %massGas            (           ))
  else if (weightBy_ == weightByLuminosity) then
     luminosities      =             self        %luminositiesStellar(           )
     massStellar       =  max (0.0d0,luminosities%luminosity         (weightIndex))
     massGas           =       0.0d0
  else
     massDistribution_ => null()
     return
  end if
  ! Build the spheroid mass distribution.
  radiusScale       =  self%radius()
  massDistribution_ => Node_Component_Mass_Distribution_Scaled_Spherical(includeStars,includeGas,massStellar,massGas,radiusScale,massDistributionStellar_,massDistributionGas_,kinematicDistribution_,scalerStellarPool,scalerGasPool)
  return
end function Node_Component_Spheroid_Standard_Mass_Distribution

double precision function Node_Component_Spheroid_Standard_Mass_Baryonic(self) result(massBaryonic)
  !!{RST
  Return the baryonic mass for the standard spheroid component.
  !!}
  implicit none
  class(nodeComponentSpheroidStandard), intent(inout) :: self

  massBaryonic=+max(0.0d0,self%massStellar()) &
       &       +max(0.0d0,self%massGas    ())
  return
end function Node_Component_Spheroid_Standard_Mass_Baryonic

double precision function Node_Component_Spheroid_Standard_Half_Mass_Radius(self) result(radiusHalfMass)
  !!{RST
  Return the half-mass radius of the standard spheroid.
  !!}
  use :: Error                                , only : Error_Report
  use :: Mass_Distributions                   , only : massDistributionSpherical
  use :: Node_Component_Spheroid_Standard_Data, only : massDistributionStellar_
  implicit none
  class(nodeComponentSpheroidStandard), intent(inout) :: self


  select type (massDistributionStellar_)
  class is (massDistributionSpherical)
     radiusHalfMass=+massDistributionStellar_%radiusHalfMass() &
          &         *self                    %radius        ()
  class default
     radiusHalfMass=0.0d0
     call Error_Report('spheroid mass distribution is not spherically-symmetric'//{introspection:location})
  end select
  return
end function Node_Component_Spheroid_Standard_Half_Mass_Radius
