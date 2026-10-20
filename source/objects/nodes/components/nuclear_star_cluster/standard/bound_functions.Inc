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

  !+    Contributions to this file made by: Matías Liempi, Claude.

!!{RST
Contains custom functions for the standard nuclear star cluster component.
!!}

double precision function Node_Component_NSC_Standard_Radius(self) result(radius)
   !!{RST
   Return the radius of the standard nuclear star cluster. If the mass distribution of the nuclear star cluster is tidally
   truncated, this is the scale radius of its untruncated mass distribution, and so is computed from the mass of that
   untruncated distribution.
   !!}
   use :: Node_Component_NSC_Standard_Data, only : radiusNormalization, galacticComponentTruncation_, tidalTruncationID
   implicit none
   class           (nodeComponentNSCStandard), intent(inout) :: self
   double precision                                          :: massNuclearStarCluster
   double precision                          , parameter     :: massNormalization     = 1.0d6 ! M⊙

   massNuclearStarCluster=+self%massGas    () &
        &                 +self%massStellar()
   if (galacticComponentTruncation_%isActive()) then
      if (self%floatRank0MetaPropertyIsCreated(tidalTruncationID)) massNuclearStarCluster=+massNuclearStarCluster                                    &
           &                                                                              *exp(+self%floatRank0MetaPropertyGet(tidalTruncationID))
   end if
   ! The radius is computed using the radius model from F. Antonini, E. Barausse & J. Silk (2015; https://ui.adsabs.harvard.edu/abs/2015ApJ...812...72A).
   if (massNuclearStarCluster > 0.0d0) then
      radius=+radiusNormalization          &
           & *sqrt(                        &
           &       +massNuclearStarCluster &
           &       /massNormalization      &
           &      )
   else 
      radius=+0.0d0
   end if 
   return
end function Node_Component_NSC_Standard_Radius

double precision function Node_Component_NSC_Standard_Velocity(self)
   !!{RST
   Return the circular velocity of the standard nuclear star cluster.
   !!}
   use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
   implicit none
   class(nodeComponentNSCStandard), intent(inout) :: self
   double precision                               :: radiusNuclearStarCluster, massNuclearStarCluster

   radiusNuclearStarCluster=+self%radius     ()
   massNuclearStarCluster  =+self%massGas    () &
        &                   +self%massStellar()
   ! The velocity is returned in km/s units.
   if     (                                  &
        &   massNuclearStarCluster   > 0.0d0 &
        &  .and.                             &
        &   radiusNuclearStarCluster > 0.0d0 &
        & ) then
      Node_Component_NSC_Standard_Velocity=+sqrt(                                &
           &                                     +gravitationalConstant_internal &
           &                                     *massNuclearStarCluster         &
           &                                     /radiusNuclearStarCluster       &
           &                                    )
   else 
      Node_Component_NSC_Standard_Velocity=+0.0d0
   end if 
   return
end function Node_Component_NSC_Standard_Velocity

function Node_Component_NSC_Standard_Mass_Distribution(self,componentType,massType,weightBy,weightIndex) result(massDistribution_)
  !!{RST
  Return the mass distribution for the standard nuclear star cluster component.
  !!}
  use :: Mass_Distributions                         , only : massDistributionClass                            , massDistributionMatches_
  use :: Node_Component_NSC_Standard_Data           , only : massDistributionStellar_                         , massDistributionGas_                                , kinematicDistribution_      , scalerStellarPool        , &
       &                                                     scalerGasPool                                    , galacticComponentTruncation_                        , tidalTruncationID           , galacticComponentHeating_, &
       &                                                     kinematicDistributionGas_
  use :: Node_Components_Galactic_Mass_Distributions, only : Node_Component_Mass_Distribution_Scaled_Spherical, Node_Component_Mass_Distribution_Truncated_Spherical
  use :: Galactic_Structure_Options                 , only : componentTypeNuclearStarCluster                  , massTypeStellar                                     , massTypeGaseous             , enumerationWeightByType  , &
       &                                                     weightByMass                                     , weightByLuminosity                                  , enumerationComponentTypeType, enumerationMassTypeType  , &
       &                                                     componentTypeDarkMatterOnly
  use :: Error                                      , only : Error_Report
  implicit none
  class           (massDistributionClass          ), pointer                 :: massDistribution_, massDistributionDarkMatter_
  class           (nodeComponentDarkMatterProfile ), pointer                 :: darkMatterProfile
  class           (nodeComponentNSCStandard       ), intent(inout)           :: self
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
  
  ! Determine which components of the nuclear star cluster to include.
  includeGas  =massDistributionMatches_(componentTypeNuclearStarCluster,massTypeGaseous,componentType,massType) .and.  weightBy_ == weightByMass
  includeStars=massDistributionMatches_(componentTypeNuclearStarCluster,massTypeStellar,componentType,massType) .and. (weightBy_ == weightByMass .or. weightBy_ == weightByLuminosity)  
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
  ! Build the nuclear star cluster mass distribution, truncated if the nuclear star cluster has been tidally stripped, and heated if it has been tidally
  ! heated.
  radiusScale         =  self%radius()
  fractionMassRetained=  1.0d0
  if (galacticComponentTruncation_%isActive()) then
     if (self%floatRank0MetaPropertyIsCreated(tidalTruncationID)) fractionMassRetained=exp(-self%floatRank0MetaPropertyGet(tidalTruncationID))
  end if
  ! Determine whether the nuclear star cluster is heated. Heating is applied only to satellites.
  isHeated=(includeStars .or. includeGas) .and. galacticComponentHeating_%isActive()
  if (isHeated) isHeated=self%hostNode%isSatellite()
  if (isHeated) then
     ! Heat the nuclear star cluster consistently with the heated dark-matter-only mass distribution of its host halo.
     darkMatterProfile           => self             %hostNode%darkMatterProfile(                                         )
     massDistributionDarkMatter_ => darkMatterProfile         %massDistribution (componentType=componentTypeDarkMatterOnly)
     if (.not.associated(massDistributionDarkMatter_)) call Error_Report('heating of the nuclear star cluster requires a dark matter mass distribution'//{introspection:location})
     massDistribution_ => Node_Component_Mass_Distribution_Truncated_Spherical(includeStars,includeGas,massStellar,massGas,radiusScale,massDistributionStellar_,massDistributionGas_,kinematicDistribution_,kinematicDistributionGas_,galacticComponentTruncation_,fractionMassRetained,galacticComponentHeating_,massDistributionDarkMatter_)
     !![
     <objectDestructor name="massDistributionDarkMatter_"/>
     !!]
  else if (fractionMassRetained < 1.0d0) then
     massDistribution_ => Node_Component_Mass_Distribution_Truncated_Spherical(includeStars,includeGas,massStellar,massGas,radiusScale,massDistributionStellar_,massDistributionGas_,kinematicDistribution_,kinematicDistributionGas_,galacticComponentTruncation_,fractionMassRetained)
  else
     massDistribution_ => Node_Component_Mass_Distribution_Scaled_Spherical   (includeStars,includeGas,massStellar,massGas,radiusScale,massDistributionStellar_,massDistributionGas_,kinematicDistribution_,kinematicDistributionGas_,scalerStellarPool,scalerGasPool                      )
  end if
  return
end function Node_Component_NSC_Standard_Mass_Distribution

double precision function Node_Component_NSC_Standard_Mass_Baryonic(self) result(massBaryonic)
  !!{RST
  Return the baryonic mass for the standard :term:`NSC` component.
  !!}
  implicit none
  class(nodeComponentNSCStandard), intent(inout) :: self

  massBaryonic=+max(0.0d0,self%massStellar()) &
       &       +max(0.0d0,self%massGas    ())
  return
end function Node_Component_NSC_Standard_Mass_Baryonic

double precision function Node_Component_NSC_Standard_Angular_Momentum(self) result(angularMomentum)
   !!{RST
   Return the specific angular momentum of the standard :term:`NSC` at a given radius.
   !!}
   implicit none
   class           (nodeComponentNSCStandard), intent(inout) :: self
   double precision                                          :: massDynamical, velocity, &
        &                                                       radius
   
   ! Get the velocity and the dynamical mass of the nuclear star cluster.
   radius       =+self%radius     ()
   velocity     =+self%velocity   ()
   massDynamical=+self%massGas    () &
        &        +self%massStellar()
   if     (                       &
        &   massDynamical > 0.0d0 &
        &  .and.                  &
        &   radius        > 0.0d0 &
        & ) then
      angularMomentum=+massDynamical &
           &          *velocity      &
           &          *radius
   else
       angularMomentum=+0.0d0
   end if
   return
end function Node_Component_NSC_Standard_Angular_Momentum

double precision function Node_Component_NSC_Standard_Half_Mass_Radius(self) result(radiusHalfMass)
  !!{RST
  Return the half-mass radius of the standard nuclear star cluster, accounting for any tidal truncation and heating of
  its mass distribution.
  !!}
  use :: Error                           , only : Error_Report
  use :: Mass_Distributions              , only : massDistributionSpherical, massDistributionClass
  use :: Node_Component_NSC_Standard_Data, only : massDistributionStellar_ , galacticComponentTruncation_, tidalTruncationID, galacticComponentHeating_
  use :: Galactic_Structure_Options      , only : massTypeStellar          , massTypeGaseous
  implicit none
  class           (nodeComponentNSCStandard ), intent(inout) :: self
  class           (massDistributionSpherical), pointer       :: massDistributionTruncated
  class           (massDistributionClass    ), pointer       :: massDistributionHeated
  double precision                                           :: fractionMassRetained

  ! If the nuclear star cluster is heated, find the half-mass radius from its heated mass distribution. Stars and gas are heated identically,
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
              call Error_Report('nuclear star cluster mass distribution is not spherically-symmetric'//{introspection:location})
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
     radiusHalfMass=+0.0d0
     call Error_Report('nuclear star cluster mass distribution is not spherically-symmetric'//{introspection:location})
  end select
  return
end function Node_Component_NSC_Standard_Half_Mass_Radius
