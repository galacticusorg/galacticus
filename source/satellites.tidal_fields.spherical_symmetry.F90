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

  !!{
  Implements a model of the tidal field acting on a satellite assuming spherical symmetry in the host.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <satelliteTidalField name="satelliteTidalFieldSphericalSymmetry">
   <description>
    A satellite tidal field class which assumes a spherically-symmetric host halo, and computes the tidal field accordingly
    using:
    \begin{equation}
     \mathcal{F} = f_\mathrm{boost} \left[ {\mathrm{G} M_\mathrm{host}(&lt;r_\mathrm{p}) \over r_\mathrm{p}^3} - 4 \pi \mathrm{G}
     \rho_\mathrm{host}(r_\mathrm{p}) + \omega_\mathrm{p}^2 \right],
    \end{equation}
    where $r_\mathrm{p}$ is the pericentric radius. $M_\mathrm{host}(&lt;r)$ is the mass of the host halo enclosed within a sphere
    of radius $r$, $\rho_\mathrm{host}(r)$ is the host density at radius $r$, and $\omega_\mathrm{p}$ is the orbital angular
    velocity at pericenter. The term $f_\mathrm{boost}=${\normalfont \ttfamily [factorBoost]} scales the overall tidal field. Note
    that the centrifugal term, $\omega_\mathrm{p}^2$, is included only if the {\normalfont \ttfamily
    includeCentrifugalAcceleration} argument is set to true. The tidal field is evaluated at the current orbital position of the
    satellite by default, but can be evaluated at the orbital pericenter if the {\normalfont \ttfamily atPericenter} argument is
    set to true.
   </description>
  </satelliteTidalField>
  !!]
  type, extends(satelliteTidalFieldClass) :: satelliteTidalFieldSphericalSymmetry
     !!{
     Implementation of a satellite tidal friction class which assumes spherical symmetry.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: factorBoost
   contains
     !![
     <methods>
       <method method="factors" description="Compute factors needed for tidal tensor calculation."/>
     </methods>
     !!]
     final     ::                        sphericalSymmetryDestructor
     procedure :: tidalTensor         => sphericalSymmetryTidalTensor
     procedure :: tidalTensorRadial   => sphericalSymmetryTidalTensorRadial
     procedure :: tidalTensorDominant => sphericalSymmetryTidalTensorDominant
     procedure :: factors             => sphericalSymmetryFactors
  end type satelliteTidalFieldSphericalSymmetry

  interface satelliteTidalFieldSphericalSymmetry
     !!{
     Constructors for the sphericalSymmetry satellite tidal field class.
     !!}
     module procedure sphericalSymmetryConstructorParameters
     module procedure sphericalSymmetryConstructorInternal
  end interface satelliteTidalFieldSphericalSymmetry

contains

  function sphericalSymmetryConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{satelliteTidalFieldSphericalSymmetry} satellite tidal field class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteTidalFieldSphericalSymmetry)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass            ), pointer       :: darkMatterHaloScale_
    double precision                                                      :: factorBoost

    !![
    <inputParameter>
      <name>factorBoost</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The factor by which to boost satellite tidal fields in the \mono{sphericalSymmetry} tidal field class.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=satelliteTidalFieldSphericalSymmetry(factorBoost,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function sphericalSymmetryConstructorParameters

  function sphericalSymmetryConstructorInternal(factorBoost,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{satelliteTidalFieldSphericalSymmetry} satellite tidal field class.
    !!}
    implicit none
    type            (satelliteTidalFieldSphericalSymmetry)                        :: self
    class           (darkMatterHaloScaleClass            ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                      , intent(in   )         :: factorBoost
    !![
    <constructorAssign variables="factorBoost, *darkMatterHaloScale_"/>
    !!]

    return
  end function sphericalSymmetryConstructorInternal

  subroutine sphericalSymmetryDestructor(self)
    !!{
    Destructor for the \refClass{satelliteTidalFieldSphericalSymmetry} satellite tidal field class.
    !!}
    implicit none
    type(satelliteTidalFieldSphericalSymmetry), intent(inout) :: self
    
    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine sphericalSymmetryDestructor

  subroutine sphericalSymmetryFactors(self,node,nodeHost,atPericenter,massEnclosedHost,densityHost,radiusOrbital,velocityOrbital)
    !!{
    Compute relevant quantities for spherically symmetric tidal field calculations.
    !!}
    use :: Coordinates       , only : coordinateCylindrical                           , assignment(=)
    use :: Galacticus_Nodes  , only : nodeComponentSatellite                          , treeNode
    use :: Kepler_Orbits     , only : keplerOrbit
    use :: Mass_Distributions, only : massDistributionClass
    use :: Satellite_Orbits  , only : Satellite_Orbit_Extremum_Phase_Space_Coordinates, extremumPericenter
    use :: Vectors           , only : Vector_Magnitude
    implicit none
    class           (satelliteTidalFieldSphericalSymmetry), intent(inout)                   :: self
    type            (treeNode                            ), intent(inout)                   :: node
    type            (treeNode                            ), intent(inout), optional, target :: nodeHost
    logical                                               , intent(in   )                   :: atPericenter
    double precision                                      , intent(  out)                   :: massEnclosedHost  , densityHost    , &
         &                                                                                     radiusOrbital     , velocityOrbital
    type            (treeNode                            ), pointer                         :: nodeHost_
    class           (nodeComponentSatellite              ), pointer                         :: satellite
    class           (massDistributionClass               ), pointer                         :: massDistribution_
    type            (keplerOrbit                         )                                  :: orbit
    type            (coordinateCylindrical               )                                  :: coordinatesOrbital

    ! Find the host node.
    if (present(nodeHost)) then
       nodeHost_ => nodeHost
    else
       nodeHost_ => node    %parent
    end if
    ! Get the satellite component.
    satellite => node%satellite()
    ! Get the position and velocity at which the tidal field should be evaluated.
    if (atPericenter) then
       ! Get the orbital radius and velocity at pericenter.
       orbit=satellite%virialOrbit()
       call Satellite_Orbit_Extremum_Phase_Space_Coordinates(nodeHost_,orbit,extremumPericenter,radiusOrbital,velocityOrbital,self%darkMatterHaloScale_)
    else
       ! Use the current orbital position/velocity.
       radiusOrbital  =Vector_Magnitude(satellite%position())
       velocityOrbital=Vector_Magnitude(satellite%velocity())
    end if
    ! Find the mass and density of the host halo at pericenter.
    coordinatesOrbital =  [radiusOrbital,0.0d0,0.0d0]
    massDistribution_  => nodeHost_        %massDistribution    (                  )
    densityHost        =  massDistribution_%density             (coordinatesOrbital)
    massEnclosedHost   =  massDistribution_%massEnclosedBySphere(     radiusOrbital)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end subroutine sphericalSymmetryFactors
  
  function sphericalSymmetryTidalTensor(self,node,nodeHost,atPericenter,includeCentrifugalAcceleration) result(tidalTensor)
    !!{
    Return the tidal tensor for satellite halos assuming spherical symmetry of the host.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Tensors                         , only : tensorNullR2D3Sym             , tensorIdentityR2D3Sym, assignment(=), operator(*)
    use :: Vectors                         , only : Vector_Outer_Product
    implicit none
    type            (tensorRank2Dimension3Symmetric      )                                  :: tidalTensor
    class           (satelliteTidalFieldSphericalSymmetry), intent(inout)                   :: self
    type            (treeNode                            ), intent(inout)                   :: node
    type            (treeNode                            ), intent(inout), optional, target :: nodeHost
    logical                                               , intent(in   ), optional         :: atPericenter     , includeCentrifugalAcceleration
    double precision                                      , dimension(3)                    :: positionCartesian
    type            (tensorRank2Dimension3Symmetric      )                                  :: positionTensor   , accelerationTensor
    double precision                                                                        :: densityHost      , massEnclosedHost              , &
         &                                                                                     radiusOrbital    , velocityOrbital
    !![
    <optionalArgument name="atPericenter"                   defaultsTo=".false."/>
    <optionalArgument name="includeCentrifugalAcceleration" defaultsTo=".false."/>
    !!]
    
    ! For isolated halos, always return zero tidal field.
    if (.not.node%isSatellite().and..not.present(nodeHost)) then
       tidalTensor=tensorNullR2D3Sym
       return
    end if
    ! Get required factors.
    call self%factors(node,nodeHost,atPericenter_,massEnclosedHost,densityHost,radiusOrbital,velocityOrbital)
    ! Construct the tidal tensor.
    positionCartesian=[radiusOrbital,0.0d0,0.0d0]
    positionTensor   =Vector_Outer_Product(positionCartesian,symmetrize=.true.)
    tidalTensor      =+gravitationalConstant_internal                                       &
         &            *(                                                                    &
         &              -(massEnclosedHost         /radiusOrbital**3)*tensorIdentityR2D3Sym &
         &              +(massEnclosedHost*3.0d0   /radiusOrbital**5)*positionTensor        &
         &              -(densityHost     *4.0d0*Pi/radiusOrbital**2)*positionTensor        &
         &             )
    ! Add centrifugal term if requested.
    if (includeCentrifugalAcceleration_) then
       accelerationTensor= tensorRank2Dimension3Symmetric(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0)
       tidalTensor       =+tidalTensor        &
            &             +accelerationTensor &
            &             *velocityOrbital**2 &
            &             /radiusOrbital  **2
    end if
    ! Boost the tidal field.
    tidalTensor=+     tidalTensor &
         &      *self%factorBoost
    return
  end function sphericalSymmetryTidalTensor

  double precision function sphericalSymmetryTidalTensorRadial(self,node,nodeHost,atPericenter,includeCentrifugalAcceleration) result(tidalTensorRadial)
    !!{
    Return the radial part of the tidal tensor for satellite halos assuming spherical symmetry of the host.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (satelliteTidalFieldSphericalSymmetry), intent(inout)                   :: self
    type            (treeNode                            ), intent(inout)                   :: node
    type            (treeNode                            ), intent(inout), optional, target :: nodeHost
    logical                                               , intent(in   ), optional         :: atPericenter , includeCentrifugalAcceleration
    double precision                                                                        :: densityHost  , massEnclosedHost              , &
         &                                                                                     radiusOrbital, velocityOrbital
    !![
    <optionalArgument name="atPericenter"                   defaultsTo=".false."/>
    <optionalArgument name="includeCentrifugalAcceleration" defaultsTo=".false."/>
    !!]
    
    ! For isolated halos, always return zero tidal field.
    if (.not.node%isSatellite().and..not.present(nodeHost)) then
       tidalTensorRadial=+0.0d0
       return
    end if
    ! Get required factors.
    call self%factors(node,nodeHost,atPericenter_,massEnclosedHost,densityHost,radiusOrbital,velocityOrbital)
    ! Compute the tidal field.
    tidalTensorRadial       =+2.0d0   *gravitationalConstant_internal*massEnclosedHost/radiusOrbital**3 &
         &                   -4.0d0*Pi*gravitationalConstant_internal*densityHost
    ! Add centrifugal term if requested.
    if (includeCentrifugalAcceleration_)         &
         & tidalTensorRadial=+tidalTensorRadial  &
         &                   +velocityOrbital**2 &
         &                   /radiusOrbital  **2
    ! Boost the tidal field.
    tidalTensorRadial=+      tidalTensorRadial &
         &            *self%factorBoost
    return
  end function sphericalSymmetryTidalTensorRadial

  double precision function sphericalSymmetryTidalTensorDominant(self,node,nodeHost,atPericenter,includeCentrifugalAcceleration) result(tidalTensorDominant)
    !!{
    Return the dominant eigenvalue of the tidal tensor for satellite halos assuming spherical symmetry of the host (in this case equal to the radial
    component, computed via the \mono{tidalTensorRadial} method).
    !!}
    class  (satelliteTidalFieldSphericalSymmetry), intent(inout)                   :: self
    type   (treeNode                            ), intent(inout)                   :: node
    type   (treeNode                            ), intent(inout), optional, target :: nodeHost
    logical                                      , intent(in   ), optional         :: atPericenter, includeCentrifugalAcceleration

    ! The dominant term is just equal to the radial term.
    tidalTensorDominant=self%tidalTensorRadial(node,nodeHost,atPericenter,includeCentrifugalAcceleration)
    return
  end function sphericalSymmetryTidalTensorDominant
  
