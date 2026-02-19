!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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

  !+    Contributions to this file made by: Matías Liempi

!!{
Contains custom functions for the standard nuclear star cluster component.
!!}

double precision function Node_Component_NSC_Standard_Radius(self) result(radius)
   !!{
   Return the radius of the standard nuclear star cluster.
   !!}
   use :: Node_Component_NSC_Standard_Data, only : radiusNormalization
   implicit none
   class           (nodeComponentNSCStandard), intent(inout) :: self
   double precision                                          :: massNuclearStarCluster
   double precision                          , parameter     :: massNormalization     = 1.0d6 ! M⊙

   massNuclearStarCluster=+self%massGas    () &
        &                 +self%massStellar()
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
   !!{
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
  !!{
  Return the mass distribution for the standard nuclear star cluster component.
  !!}
  use :: Mass_Distributions              , only : massDistributionClass          , massDistributionSphericalScaler, massDistributionComposite   , massDistributionList   , &
        &                                         massDistributionSpherical      , massDistributionMatches_
  use :: Node_Component_NSC_Standard_Data, only : massDistributionStellar_       , massDistributionGas_           , kinematicDistribution_
  use :: Galactic_Structure_Options      , only : componentTypeNuclearStarCluster, massTypeStellar                , massTypeGaseous             , enumerationWeightByType, &
       &                                          weightByMass                   , weightByLuminosity             , enumerationComponentTypeType, enumerationMassTypeType
  implicit none
  class           (massDistributionClass          ), pointer                 :: massDistribution_
  class           (nodeComponentNSCStandard       ), intent(inout)           :: self
  type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
  type            (enumerationMassTypeType        ), intent(in   ), optional :: massType
  type            (enumerationWeightByType        ), intent(in   ), optional :: weightBy
  integer                                          , intent(in   ), optional :: weightIndex
  type            (massDistributionSphericalScaler), pointer                 :: massDistributionStellar   , massDistributionGas
  type            (massDistributionComposite      ), pointer                 :: massDistributionTotal
  type            (massDistributionList           ), pointer                 :: massDistributionComponents
  type            (stellarLuminosities            ), save                    :: luminosities
  !$omp threadprivate(luminosities)
  double precision                                                           :: massStellar               , massGas            , &
       &                                                                        radiusScale
  logical                                                                    :: includeGas                , includeStars
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
  ! Determine which components to build.
  radiusScale=self%radius()
  if (radiusScale <= 0.0d0 .or. .not.(includeGas .or. includeStars)) then
     ! Nuclear star cluster has non-positive size, or no components matched. Return a null distribution.
     massDistribution_ => null()
  else
     ! Build the individual distributions.
     massDistributionStellar => null()
     massDistributionGas     => null()
     if (includeStars) then
        allocate(massDistributionStellar)
        select type (massDistributionStellar_)
        class is (massDistributionSpherical)
           !![
           <referenceConstruct object="massDistributionStellar" constructor="massDistributionSphericalScaler(factorScalingLength=radiusScale,factorScalingMass=massStellar,massDistribution_=massDistributionStellar_)"/>
           !!]
        end select
        call massDistributionStellar%setKinematicsDistribution(kinematicDistribution_)
     end if
     if (includeGas  ) then
        allocate(massDistributionGas  )
        select type (massDistributionGas_    )
        class is (massDistributionSpherical)
           !![
           <referenceConstruct object="massDistributionGas"     constructor="massDistributionSphericalScaler(factorScalingLength=radiusScale,factorScalingMass=massGas    ,massDistribution_=massDistributionGas_    )"/>
           !!]
        end select
        call massDistributionGas    %setKinematicsDistribution(kinematicDistribution_)
     end if
     ! Combine the distributions as necessary.
     if      (includeStars .and. includeGas) then
        ! Wrap the dimensionless mass distribution inside scaler classes to allow us to re-scale it to any nuclear star cluster system, and then composite those.
        allocate(massDistributionTotal          )
        allocate(massDistributionComponents     )
        allocate(massDistributionComponents%next)
        massDistributionComponents     %massDistribution_ => massDistributionStellar
        massDistributionComponents%next%massDistribution_ => massDistributionGas
        !![
        <referenceConstruct object="massDistributionTotal" constructor="massDistributionComposite(massDistributionComponents)"/>
        <objectDestructor name="massDistributionStellar"/>
        <objectDestructor name="massDistributionGas"    />
        !!]
        nullify(massDistributionComponents)
        ! Return a pointer to the nuclear star cluster mass distribution.
        massDistribution_ => massDistributionTotal
     else if (includeStars                 ) then
        ! Return just the stellar component.
        massDistribution_ => massDistributionStellar
     else if (                   includeGas) then
        ! Return just the gas component.
        massDistribution_ => massDistributionGas
     end if
  end if
  return
end function Node_Component_NSC_Standard_Mass_Distribution

double precision function Node_Component_NSC_Standard_Mass_Baryonic(self) result(massBaryonic)
  !!{
  Return the baryonic mass for the standard \gls{nsc} component.
  !!}
  implicit none
  class(nodeComponentNSCStandard), intent(inout) :: self

  massBaryonic=+max(0.0d0,self%massStellar()) &
       &       +max(0.0d0,self%massGas    ())
  return
end function Node_Component_NSC_Standard_Mass_Baryonic

double precision function Node_Component_NSC_Standard_Angular_Momentum(self) result(angularMomentum)
   !!{
   Return the specific angular momentum of the standard \gls{nsc} at a given radius.
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
  !!{
  Return the half-mass radius of the standard nuclear star cluster.
  !!}
  use :: Error                           , only : Error_Report
  use :: Mass_Distributions              , only : massDistributionSpherical
  use :: Node_Component_NSC_Standard_Data, only : massDistributionStellar_
  implicit none
  class(nodeComponentNSCStandard), intent(inout) :: self

 select type (massDistributionStellar_)
  class is (massDistributionSpherical)
     radiusHalfMass=+massDistributionStellar_%radiusHalfMass() &
          &         *self                    %radius        ()
  class default
     radiusHalfMass=+0.0d0
     call Error_Report('nuclear star cluster mass distribution is not spherically-symmetric'//{introspection:location})
  end select
  return
end function Node_Component_NSC_Standard_Half_Mass_Radius
