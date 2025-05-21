!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Implements an N-body dark matter halo mass error class which
implements a model for errors in spherical overdensity halo finders.
!!}

  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScale , darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMO, darkMatterProfileDMOClass

  !![
  <nbodyHaloMassError name="nbodyHaloMassErrorSOHaloFinder">
   <description>An N-body dark matter halo mass error class which implements a model for errors in spherical overdensity halo finders.</description>
  </nbodyHaloMassError>
  !!]
  type, extends(nbodyHaloMassErrorClass) :: nbodyHaloMassErrorSOHaloFinder
     !!{
     An N-body halo mass error class which implements a model for errors in spherical
     overdensity halo finders.
     !!}
     private
     double precision                                     :: massParticle
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_  => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
   contains
     final     ::                    soHaloFinderDestructor
     procedure :: errorFractional => soHaloFinderErrorFractional
     procedure :: correlation     => soHaloFinderCorrelation
  end type nbodyHaloMassErrorSOHaloFinder

  interface nbodyHaloMassErrorSOHaloFinder
     !!{
     Constructors for the \refClass{nbodyHaloMassErrorSOHaloFinder} N-body halo mass error class.
     !!}
     module procedure soHaloFinderParameters
     module procedure soHaloFinderInternal
  end interface nbodyHaloMassErrorSOHaloFinder

contains

  function soHaloFinderParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nbodyHaloMassErrorSOHaloFinder} N-body halo mass error class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyHaloMassErrorSOHaloFinder)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass     ), pointer       :: darkMatterProfileDMO_
    double precision                                                :: massParticle

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>massParticle</name>
      <source>parameters</source>
      <description>Mass of particle in the simulation to which the spherical overdensity algorithm was applied.</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=nbodyHaloMassErrorSOHaloFinder(darkMatterHaloScale_,darkMatterProfileDMO_,massParticle)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function soHaloFinderParameters

  function soHaloFinderInternal(darkMatterHaloScale_,darkMatterProfileDMO_,massParticle) result(self)
    !!{
    Internal constructor for the \refClass{nbodyHaloMassErrorSOHaloFinder} N-body halo mass error class.
    !!}
    implicit none
    type            (nbodyHaloMassErrorSOHaloFinder)                        :: self
    class           (darkMatterHaloScaleClass      ), target, intent(in   ) :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass     ), target, intent(in   ) :: darkMatterProfileDMO_
    double precision                                        , intent(in   ) :: massParticle
    !![
    <constructorAssign variables="*darkMatterHaloScale_, *darkMatterProfileDMO_, massParticle"/>
    !!]

    return
  end function soHaloFinderInternal

  subroutine soHaloFinderDestructor(self)
    !!{
    Destructor for the \refClass{nbodyHaloMassErrorSOHaloFinder} N-body halo mass error class.
    !!}
    implicit none
    type(nbodyHaloMassErrorSOHaloFinder), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine soHaloFinderDestructor

  double precision function soHaloFinderErrorFractional(self,node)
    !!{
    Return the fractional error on the mass of an N-body halo in the power-law error model.
    !!}
    use :: Coordinates             , only : coordinateSpherical  , assignment(=)
    use :: Galacticus_Nodes        , only : nodeComponentBasic   , treeNode
    use :: Mass_Distributions      , only : massDistributionClass
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (nbodyHaloMassErrorSOHaloFinder), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    class           (nodeComponentBasic            ), pointer       :: basic
    class           (massDistributionClass         ), pointer       :: massDistribution_
    double precision                                , parameter     :: errorConstant                =0.014d0
    type            (coordinateSpherical           )                :: coordinates
    double precision                                                :: radiusHalo                           , densityOuterRadius, &
         &                                                             densityRatioInternalToSurface        , particleCount     , &
         &                                                             errorFractionalFixedSphere

    ! Get the basic component of the node.
    basic                         =>  node                     %basic       (               )
    ! Determine number of particles in the halo.
    particleCount                 =  +basic%mass        () &
         &                           /self %massParticle
    ! Fractional error in mass within fixed sphere, assuming Poisson statistics (which should be valid for a halo which contains a
    ! fraction of all particles in the simulation that is much less than unity).
    errorFractionalFixedSphere    =  +1.0d0               &
         &                           /sqrt(particleCount)
    ! Get the outer radius of the halo.
    radiusHalo                    =  +self             %darkMatterHaloScale_ %radiusVirial(node       )
    ! Get the density at the edge of the halo.
    massDistribution_             =>  self             %darkMatterProfileDMO_%get         (node       )
    coordinates                   =  [radiusHalo,0.0d0,0.0d0]
    densityOuterRadius            =  +massDistribution_                      %density     (coordinates)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Find the ratio of the mean interior density in the halo to the density at the halo outer radius.
    densityRatioInternalToSurface =  +3.0d0                 &
         &                           *basic%mass()          &
         &                           /4.0d0                 &
         &                           /Pi                    &
         &                           /radiusHalo        **3 &
         &                           /densityOuterRadius
    ! Compute the total error accounting for the change in size of the spherical region defining the halo.
    soHaloFinderErrorFractional   =  +sqrt(                                      &
         &                                 +errorFractionalFixedSphere       **2 &
         &                                 *(                                    &
         &                                   +1.0d0                              &
         &                                   +1.0d0                              &
         &                                   /(                                  &
         &                                     +densityRatioInternalToSurface    &
         &                                     -1.0d0                            &
         &                                    )                              **2 &
         &                                  )                                    &
         &                                 +errorConstant                    **2 &
         &                                )
    return
  end function soHaloFinderErrorFractional

  double precision function soHaloFinderCorrelation(self,node1,node2)
    !!{
    Return the correlation of the masses of a pair of N-body halos.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(nbodyHaloMassErrorSOHaloFinder), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node1 , node2
    class(nodeComponentBasic            ), pointer       :: basic1, basic2
    !$GLC attributes unused :: self

    basic1 => node1%basic()
    basic2 => node2%basic()
    if     (                                &
         &   basic1%mass() == basic2%mass() &
         &  .and.                           &
         &   basic1%time() == basic2%time() &
         & ) then
       soHaloFinderCorrelation=1.0d0
    else
       soHaloFinderCorrelation=0.0d0
    end if
    return
  end function soHaloFinderCorrelation

