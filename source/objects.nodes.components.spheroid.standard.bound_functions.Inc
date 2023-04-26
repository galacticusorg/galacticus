!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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

!!{
Contains custom functions for the standard spheroid component.
!!}

function Node_Component_Spheroid_Standard_Mass_Distribution(self,weightBy,weightIndex) result(massDistribution_)
  !!{
  Return the mass distribution for the standard spheroid component.
  !!}
  use :: Mass_Distributions                   , only : massDistributionClass    , massDistributionSphericalScaler, massDistributionComposite, massDistributionList   , &
        &                                              massDistributionSpherical
  use :: Node_Component_Spheroid_Standard_Data, only : massDistributionStellar_ , massDistributionGas_           , kinematicDistribution_
  use :: Galactic_Structure_Options           , only : componentTypeSpheroid    , massTypeStellar                , massTypeGaseous          , enumerationWeightByType, &
       &                                               weightByMass             , weightByLuminosity
  implicit none
  class           (massDistributionClass          ), pointer                 :: massDistribution_
  class           (nodeComponentSpheroidStandard  ), intent(inout)           :: self
  type            (enumerationWeightByType        ), intent(in   ), optional :: weightBy
  integer                                          , intent(in   ), optional :: weightIndex
  type            (massDistributionSphericalScaler), pointer                 :: massDistributionStellar   , massDistributionGas
  type            (massDistributionComposite      ), pointer                 :: massDistributionTotal
  type            (massDistributionList           ), pointer                 :: massDistributionComponents
  type            (stellarLuminosities            ), save                    :: luminosities
  !$omp threadprivate(luminosities)
  double precision                                                           :: massStellar               , massGas            , &
       &                                                                        radiusScale
  !![
  <optionalArgument name="weightBy" defaultsTo="weightByMass" />
  !!]
  
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
  radiusScale=self%radius()
  if (radiusScale <= 0.0d0) then
     ! Spheroid has non-positive size. Set an arbitrary, positive size, but also set masses to zero. This results in an empty mass
     ! distribution.
     radiusScale=1.0d0
     massStellar=0.0d0
     massGas    =0.0d0
  end if
  ! Wrap the dimensionless mass distribution inside scaler classes to allow us to re-scale it to any spheroid system, and then composite those.
  allocate(massDistributionStellar        )
  allocate(massDistributionGas            )
  allocate(massDistributionTotal          )
  allocate(massDistributionComponents     )
  allocate(massDistributionComponents%next)
  massDistributionComponents     %massDistribution_ => massDistributionStellar
  massDistributionComponents%next%massDistribution_ => massDistributionGas
  select type (massDistributionStellar_)
  class is (massDistributionSpherical)
     !![
     <referenceConstruct object="massDistributionStellar" constructor="massDistributionSphericalScaler(factorScalingLength=radiusScale,factorScalingMass=massStellar,massDistribution_=massDistributionStellar_  )"/>
     !!]
  end select
  select type (massDistributionGas_    )
  class is (massDistributionSpherical)
     !![
     <referenceConstruct object="massDistributionGas"     constructor="massDistributionSphericalScaler(factorScalingLength=radiusScale,factorScalingMass=massGas    ,massDistribution_=massDistributionGas_      )"/>
     !!]
  end select
  call massDistributionStellar%setKinematicsDistribution(kinematicDistribution_)
  call massDistributionGas    %setKinematicsDistribution(kinematicDistribution_)
  !![
  <referenceConstruct    object="massDistributionTotal"   constructor="massDistributionComposite      (                                                                                massDistributionComponents)"/>
  <objectDestructor name="massDistributionStellar"/>
  <objectDestructor name="massDistributionGas"    />
  !!]
  nullify(massDistributionComponents)
  ! Return a pointer to the spheroid mass distribution.
  massDistribution_ => massDistributionTotal
  return
end function Node_Component_Spheroid_Standard_Mass_Distribution

double precision function Node_Component_Spheroid_Standard_Half_Mass_Radius(self) result(radiusHalfMass)
  !!{
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
