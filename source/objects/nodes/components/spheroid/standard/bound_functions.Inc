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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains custom functions for the standard spheroid component.
!!}

function Node_Component_Spheroid_Standard_Mass_Distribution(self,componentType,massType,weightBy,weightIndex) result(massDistribution_)
  !!{RST
  Return the mass distribution for the standard spheroid component.
  !!}
  use :: Mass_Distributions                         , only : massDistributionClass                            , massDistributionMatches_
  use :: Node_Component_Spheroid_Standard_Data      , only : massDistributionStellar_                         , massDistributionGas_                                , kinematicDistribution_      , scalerStellarPool        , &
       &                                                     scalerGasPool                                    , galacticComponentTruncation_                        , tidalTruncationID           , galacticComponentHeating_, &
       &                                                     kinematicDistributionGas_
  use :: Node_Components_Galactic_Mass_Distributions, only : Node_Component_Mass_Distribution_Scaled_Spherical, Node_Component_Mass_Distribution_Truncated_Spherical
  use :: Galactic_Structure_Options                 , only : componentTypeSpheroid                            , massTypeStellar                                     , massTypeGaseous             , enumerationWeightByType  , &
       &                                                     weightByMass                                     , weightByLuminosity                                  , enumerationComponentTypeType, enumerationMassTypeType  , &
       &                                                     componentTypeDarkMatterOnly
  use :: Error                                      , only : Error_Report
  implicit none
  class           (massDistributionClass          ), pointer                 :: massDistribution_, massDistributionDarkMatter_
  class           (nodeComponentDarkMatterProfile ), pointer                 :: darkMatterProfile
  class           (nodeComponentSpheroidStandard  ), intent(inout)           :: self
  type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
  type            (enumerationMassTypeType        ), intent(in   ), optional :: massType
  type            (enumerationWeightByType        ), intent(in   ), optional :: weightBy
  integer                                          , intent(in   ), optional :: weightIndex
  type            (stellarLuminosities            ), save                    :: luminosities
  !$omp threadprivate(luminosities)
  double precision                                                           :: massStellar      , massGas                    , &
       &                                                                        radiusScale      , fractionMassRetained
  logical                                                                    :: includeGas       , includeStars               , &
       &                                                                        isHeated
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
  ! Build the spheroid mass distribution, truncated if the spheroid has been tidally stripped, and heated if it has been tidally
  ! heated.
  radiusScale         =  self%radius()
  fractionMassRetained=  1.0d0
  if (galacticComponentTruncation_%isActive()) then
     if (self%floatRank0MetaPropertyIsCreated(tidalTruncationID)) fractionMassRetained=exp(-self%floatRank0MetaPropertyGet(tidalTruncationID))
  end if
  ! Determine whether the spheroid is heated. Heating is applied only to satellites.
  isHeated=(includeStars .or. includeGas) .and. galacticComponentHeating_%isActive()
  if (isHeated) isHeated=self%hostNode%isSatellite()
  if (isHeated) then
     ! Heat the spheroid consistently with the heated dark-matter-only mass distribution of its host halo.
     darkMatterProfile           => self             %hostNode%darkMatterProfile(                                         )
     massDistributionDarkMatter_ => darkMatterProfile         %massDistribution (componentType=componentTypeDarkMatterOnly)
     if (.not.associated(massDistributionDarkMatter_)) call Error_Report('heating of the spheroid requires a dark matter mass distribution'//{introspection:location})
     massDistribution_ => Node_Component_Mass_Distribution_Truncated_Spherical(includeStars,includeGas,massStellar,massGas,radiusScale,massDistributionStellar_,massDistributionGas_,kinematicDistribution_,kinematicDistributionGas_,galacticComponentTruncation_,fractionMassRetained,galacticComponentHeating_,massDistributionDarkMatter_)
     !![
     <objectDestructor name="massDistributionDarkMatter_"/>
     !!]
  else if (fractionMassRetained < 1.0d0) then
     massDistribution_ => Node_Component_Mass_Distribution_Truncated_Spherical(includeStars,includeGas,massStellar,massGas,radiusScale,massDistributionStellar_,massDistributionGas_,kinematicDistribution_,kinematicDistributionGas_,galacticComponentTruncation_,fractionMassRetained)
  else
     massDistribution_ => Node_Component_Mass_Distribution_Scaled_Spherical   (includeStars,includeGas,massStellar,massGas,radiusScale,massDistributionStellar_,massDistributionGas_,kinematicDistribution_,kinematicDistributionGas_,scalerStellarPool,scalerGasPool                  )
  end if
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
  Return the half-mass radius of the standard spheroid, accounting for any tidal truncation and heating of its mass distribution.
  !!}
  use :: Error                                , only : Error_Report
  use :: Mass_Distributions                   , only : massDistributionSpherical, massDistributionClass
  use :: Node_Component_Spheroid_Standard_Data, only : massDistributionStellar_ , galacticComponentTruncation_, tidalTruncationID, galacticComponentHeating_
  use :: Galactic_Structure_Options           , only : massTypeStellar          , massTypeGaseous
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout) :: self
  class           (massDistributionSpherical    ), pointer       :: massDistributionTruncated
  class           (massDistributionClass        ), pointer       :: massDistributionHeated
  double precision                                               :: fractionMassRetained

  ! If the spheroid is heated, find the half-mass radius from its heated mass distribution. Stars and gas are heated identically,
  ! so either may be used.
  if (galacticComponentHeating_%isActive()) then
     if (self%hostNode%isSatellite()) then
        if      (self%massStellar() > 0.0d0) then
           massDistributionHeated => self%massDistribution(massType=massTypeStellar)
        else if (self%massGas    () > 0.0d0) then
           massDistributionHeated => self%massDistribution(massType=massTypeGaseous)
        else
           massDistributionHeated => null()
        end if
        if (associated(massDistributionHeated)) then
           select type (massDistributionHeated)
           class is (massDistributionSpherical)
              radiusHalfMass=massDistributionHeated%radiusHalfMass()
           class default
              radiusHalfMass=0.0d0
              call Error_Report('spheroid mass distribution is not spherically-symmetric'//{introspection:location})
           end select
           !![
           <objectDestructor name="massDistributionHeated"/>
           !!]
           return
        end if
     end if
  end if
  fractionMassRetained=1.0d0
  if (galacticComponentTruncation_%isActive()) then
     if (self%floatRank0MetaPropertyIsCreated(tidalTruncationID)) fractionMassRetained=exp(-self%floatRank0MetaPropertyGet(tidalTruncationID))
  end if
  select type (massDistributionStellar_)
  class is (massDistributionSpherical)
     if (fractionMassRetained < 1.0d0) then
        massDistributionTruncated => galacticComponentTruncation_%massDistribution(massDistributionStellar_,fractionMassRetained)
        radiusHalfMass=+massDistributionTruncated%radiusHalfMass() &
             &         *self                     %radius        ()
        !![
        <objectDestructor name="massDistributionTruncated"/>
        !!]
     else
        radiusHalfMass=+massDistributionStellar_ %radiusHalfMass() &
             &         *self                     %radius        ()
     end if
  class default
     radiusHalfMass=0.0d0
     call Error_Report('spheroid mass distribution is not spherically-symmetric'//{introspection:location})
  end select
  return
end function Node_Component_Spheroid_Standard_Half_Mass_Radius
