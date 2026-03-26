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
  Implements a model of the tidal field acting on a satellite for arbitrary geometry in the host.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <satelliteTidalField name="satelliteTidalFieldStandard">
   <description>
    A satellite tidal field class that computes the tidal field in arbitrary geometry. Note that the centrifugal term,
    $\omega_\mathrm{p}^2$, is included only if the {\normalfont \ttfamily includeCentrifugalAcceleration} argument is set to
    true. The tidal field is evaluated at the current orbital position of the satellite by default, but can be evaluated at the
    orbital pericenter if the {\normalfont \ttfamily atPericenter} argument is set to true.
   </description>
  </satelliteTidalField>
  !!]
  type, extends(satelliteTidalFieldClass) :: satelliteTidalFieldStandard
     !!{
     Implementation of a satellite tidal field class in arbitrary geometry.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     !![
     <methods>
       <method method="factors"        description="Compute factors needed for tidal tensor calculation."/>
       <method method="tidalTensorGet" description="Get the tidal tensor."                               />
     </methods>
     !!]
     final     ::                        standardDestructor
     procedure :: tidalTensor         => standardTidalTensor
     procedure :: tidalTensorRadial   => standardTidalTensorRadial
     procedure :: tidalTensorDominant => standardTidalTensorDominant
     procedure :: tidalTensorGet      => standardTidalTensorGet
     procedure :: factors             => standardFactors
  end type satelliteTidalFieldStandard

  interface satelliteTidalFieldStandard
     !!{
     Constructors for the standard satellite tidal field class.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface satelliteTidalFieldStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{satelliteTidalFieldStandard} satellite tidal field class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (satelliteTidalFieldStandard)                :: self
    type (inputParameters            ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass   ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=satelliteTidalFieldStandard(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{satelliteTidalFieldStandard} satellite tidal field class.
    !!}
    implicit none
    type (satelliteTidalFieldStandard)                        :: self
    class(darkMatterHaloScaleClass   ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !!{
    Destructor for the \refClass{satelliteTidalFieldStandard} satellite tidal field class.
    !!}
    implicit none
    type(satelliteTidalFieldStandard), intent(inout) :: self
    
    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine standardDestructor

  subroutine standardFactors(self,node,nodeHost,atPericenter,coordinatesOrbital,coordinatesOrbitalVelocity,radiusOrbital,velocityOrbital)
    !!{
    Compute relevant quantities for tidal field calculations.
    !!}
    use :: Coordinates     , only : coordinateCartesian                             , assignment(=)
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    use :: Kepler_Orbits   , only : keplerOrbit
    use :: Satellite_Orbits, only : Satellite_Orbit_Extremum_Phase_Space_Coordinates, extremumPericenter
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class           (satelliteTidalFieldStandard), intent(inout)                   :: self
    type            (treeNode                   ), intent(inout)                   :: node
    type            (treeNode                   ), intent(inout), optional, target :: nodeHost
    logical                                      , intent(in   )                   :: atPericenter
    type            (coordinateCartesian        ), intent(  out)                   :: coordinatesOrbital
    type            (coordinateCartesian        ), intent(  out), optional         :: coordinatesOrbitalVelocity
    double precision                             , intent(  out), optional         :: radiusOrbital             , velocityOrbital
    type            (treeNode                   ), pointer                         :: nodeHost_
    class           (nodeComponentSatellite     ), pointer                         :: satellite
    double precision                                                               :: radiusOrbital_            , velocityOrbital_
    type            (keplerOrbit                )                                  :: orbit

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
       call Satellite_Orbit_Extremum_Phase_Space_Coordinates(nodeHost_,orbit,extremumPericenter,radiusOrbital_,velocityOrbital_,self%darkMatterHaloScale_)
       coordinatesOrbital=[radiusOrbital_,0.0d0,0.0d0]
       if (present(radiusOrbital             )) radiusOrbital=  radiusOrbital_
       if (present(velocityOrbital           )) velocityOrbital=velocityOrbital_
       if (present(coordinatesOrbitalVelocity)) coordinatesOrbitalVelocity=[0.0d0,velocityOrbital_,0.0d0]
    else
       ! Use the current orbital position/velocity.
       coordinatesOrbital=satellite%position()
       if (present(radiusOrbital             )) radiusOrbital             =Vector_Magnitude(satellite%position())
       if (present(velocityOrbital           )) velocityOrbital           =Vector_Magnitude(satellite%velocity())
       if (present(coordinatesOrbitalVelocity)) coordinatesOrbitalVelocity=                 satellite%velocity()
    end if
    return
  end subroutine standardFactors
  
  function standardTidalTensorGet(self,node,nodeHost,atPericenter,includeCentrifugalAcceleration,isSphericallySymmetric) result(tidalTensor)
    !!{
    Return the tidal tensor for satellite halos in arbitrary host geometry.
    !!}
    use :: Coordinates       , only : coordinateCartesian  , assignment(=)
    use :: Mass_Distributions, only : massDistributionClass
    use :: Tensors           , only : tensorNullR2D3Sym    , assignment(=), operator(*)
    implicit none
    type            (tensorRank2Dimension3Symmetric)                                  :: tidalTensor
    class           (satelliteTidalFieldStandard   ), intent(inout)                   :: self
    type            (treeNode                      ), intent(inout)                   :: node
    type            (treeNode                      ), intent(inout), optional, target :: nodeHost
    logical                                         , intent(in   ), optional         :: atPericenter           , includeCentrifugalAcceleration
    logical                                         , intent(  out), optional         :: isSphericallySymmetric
    class           (massDistributionClass         ), pointer                         :: massDistribution_
    type            (treeNode                      ), pointer                         :: nodeHost_
    double precision                                , dimension(3)                    :: velocitiesOrbital
    type            (coordinateCartesian           )                                  :: coordinatesOrbital     , coordinatesOrbitalVelocity
    type            (tensorRank2Dimension3Symmetric)                                  :: accelerationTensor
    double precision                                                                  :: radiusOrbital          , velocityOrbital
    logical                                                                           :: isSphericallySymmetric_
    !![
    <optionalArgument name="atPericenter"                   defaultsTo=".false."/>
    <optionalArgument name="includeCentrifugalAcceleration" defaultsTo=".false."/>
    !!]
    
    ! For isolated halos, always return zero tidal field.
    if (.not.node%isSatellite().and..not.present(nodeHost)) then
       tidalTensor=tensorNullR2D3Sym
       ! A null tidal field is spherically symmetric.
       if (present(isSphericallySymmetric)) isSphericallySymmetric=.true.
       return
    end if
    ! Get required factors.
    call self%factors(node,nodeHost,atPericenter_,coordinatesOrbital,coordinatesOrbitalVelocity,radiusOrbital,velocityOrbital)
    ! Construct the tidal tensor.
    if (present(nodeHost)) then
       nodeHost_ => nodeHost
    else
       nodeHost_ => node    %parent
    end if
    massDistribution_       =>  nodeHost_        %massDistribution      ()
    isSphericallySymmetric_ =   massDistribution_%isSphericallySymmetric()
    if (present(isSphericallySymmetric)) isSphericallySymmetric=isSphericallySymmetric_
    if (present(isSphericallySymmetric) .and. isSphericallySymmetric_) then
       ! Evaluate along the x-axis if the spherically-symmetric version is requested.
       coordinatesOrbital=[radiusOrbital,0.0d0,0.0d0]
       tidalTensor       =+massDistribution_%tidalTensor(coordinatesOrbital)
    else
       ! Evaluate at the orbital coordinates.
       tidalTensor       =+massDistribution_%tidalTensor(coordinatesOrbital)
    end if
    if (.not.isSphericallySymmetric_ .and. atPericenter_) &
         & call Error_Report('can not compute tidal tensor at pericenter in non-spherically symmetric potentials'//{introspection:location})
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Add centrifugal term if requested.
    if (includeCentrifugalAcceleration_) then
       ! Construct the acceleration tensor.
       if (present(isSphericallySymmetric) .and. isSphericallySymmetric_) then
          ! In the spherically symmetric case this has the following form.
          accelerationTensor=+tensorRank2Dimension3Symmetric(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0)    &
               &             *velocityOrbital                                                    **2 &
               &             /radiusOrbital                                                      **2
       else
          ! Use the fully general form.
          velocitiesOrbital=coordinatesOrbitalVelocity
          accelerationTensor=tensorRank2Dimension3Symmetric(                                                                                      &
               &                                            +velocitiesOrbital(2)*velocitiesOrbital(2)+velocitiesOrbital(3)*velocitiesOrbital(3), & ! xx
               &                                            -velocitiesOrbital(1)*velocitiesOrbital(2)                                          , & ! xy
               &                                            -velocitiesOrbital(1)*velocitiesOrbital(3)                                          , & ! xz
               &                                            +velocitiesOrbital(1)*velocitiesOrbital(1)+velocitiesOrbital(3)*velocitiesOrbital(3), & ! yy
               &                                            -velocitiesOrbital(2)*velocitiesOrbital(3)                                          , & ! yz
               &                                            +velocitiesOrbital(1)*velocitiesOrbital(1)+velocitiesOrbital(2)*velocitiesOrbital(2)  & ! zz
               &                                           )                                                                                      &
               &             /radiusOrbital**2
       end if
       tidalTensor=+tidalTensor        &
            &      +accelerationTensor
    end if
    return
  end function standardTidalTensorGet

  function standardTidalTensor(self,node,nodeHost,atPericenter,includeCentrifugalAcceleration) result(tidalTensor)
    !!{
    Return the tidal tensor for satellite halos in arbitrary host geometry.
    !!}
    implicit none
    type   (tensorRank2Dimension3Symmetric)                                  :: tidalTensor
    class  (satelliteTidalFieldStandard   ), intent(inout)                   :: self
    type   (treeNode                      ), intent(inout)                   :: node
    type   (treeNode                      ), intent(inout), optional, target :: nodeHost
    logical                                , intent(in   ), optional         :: atPericenter, includeCentrifugalAcceleration
   
    tidalTensor=self%tidalTensorGet(node,nodeHost,atPericenter,includeCentrifugalAcceleration)
    return
  end function standardTidalTensor

  double precision function standardTidalTensorRadial(self,node,nodeHost,atPericenter,includeCentrifugalAcceleration) result(tidalTensorRadial)
    !!{
    Return the radial part of the tidal tensor for satellite halos in arbitrary host geometry.
    !!}
    use :: Coordinates, only : coordinateCartesian, assignment(=)
    implicit none
    class           (satelliteTidalFieldStandard   ), intent(inout)                   :: self
    type            (treeNode                      ), intent(inout)                   :: node
    type            (treeNode                      ), intent(inout), optional, target :: nodeHost
    logical                                         , intent(in   ), optional         :: atPericenter          , includeCentrifugalAcceleration
    double precision                                , dimension(3)                    :: positionOrbital
    type            (coordinateCartesian           )                                  :: coordinatesOrbital
    double precision                                                                  :: radiusOrbital
    type            (tensorRank2Dimension3Symmetric)                                  :: tidalTensor
    logical                                                                           :: isSphericallySymmetric
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
    call self%factors(node,nodeHost,atPericenter_,coordinatesOrbital,radiusOrbital=radiusOrbital)
    positionOrbital=coordinatesOrbital
    ! Get the tidal tensor.
    tidalTensor=self%tidalTensorGet(node,nodeHost,atPericenter,includeCentrifugalAcceleration,isSphericallySymmetric=isSphericallySymmetric)
    ! Extract the radial component.
    if (isSphericallySymmetric) then
       ! Spherically symmetric case.
       tidalTensorRadial=+tidalTensor%element      (0,0            )
    else
       ! Fully general case.
       tidalTensorRadial=+tidalTensor%vectorProject(positionOrbital)    &
            &            /                            radiusOrbital **2
    end if
    return
  end function standardTidalTensorRadial

  double precision function standardTidalTensorDominant(self,node,nodeHost,atPericenter,includeCentrifugalAcceleration) result(tidalTensorDominant)
    !!{
    Return the dominant tidal field from the tidal tensor for satellite halos in arbitrary host geometry. To allow for
    non-spherical mass distributions, we proceed as follows to determine the tidal field:
    
    Let $\boldsymbol{\mathsf{G}}$ be the gravitational tidal tensor evaluated at the position of the satellite. Consider a unit vector,
    $\boldsymbol{\hat{x}}$ in the satellite. The tidal field along this vector is $\boldsymbol{\mathsf{G}} \boldsymbol{\hat{x}}$. The
    radial component of the tidal field in this direction is then $\boldsymbol{\hat{x}} \boldsymbol{\mathsf{G}}
    \boldsymbol{\hat{x}}$. We want to find the maximum of the tidal field over all possible directions (i.e. all possible unit
    vectors).
    
    We can write our unit vector as
    \begin{equation}
     \boldsymbol{\hat{x}} = \sum_{i=1}^3 a_i \boldsymbol{\hat{e}}_i,
    \end{equation}
    where the $\boldsymbol{\hat{e}}_i$ are the eigenvectors of $\boldsymbol{\mathsf{G}}$ and $\sum_{i=1}^3 a_i^2 = 1$. Then since, by
    definition, $\boldsymbol{\mathsf{G}} \boldsymbol{\hat{e}}_i = \lambda_i \boldsymbol{\hat{e}}_i$, where $\lambda_i$ are the
    eigenvalues of $\boldsymbol{\mathsf{G}}$, we have that
    \begin{equation}
     \boldsymbol{\hat{x}} \boldsymbol{\mathsf{G}} \boldsymbol{\hat{x}}= \sum_{i=1}^3 a_i^2 \lambda_i.
    \end{equation}
    The sum on the right hand side of the above is a weighted average of eigenvalues. Any weighted average is maximized by
    setting the weight of the largest value to $1$, and all other weights to $0$. Therefore, our tidal field is maximized along
    the direction corresponding the eigenvector of $\boldsymbol{\mathsf{G}}$ with the largest eigenvalue. (Note that we want the
    largest positive eigenvalue, not the largest absolute eigenvalue as we're interested in stretching tidal fields, not
    compressive ones.)
    !!}
    use :: Coordinates   , only : assignment(=), coordinateCartesian
    use :: Linear_Algebra, only : assignment(=), matrix                        , vector
    use :: Tensors       , only : assignment(=), tensorRank2Dimension3Symmetric
    implicit none
    class           (satelliteTidalFieldStandard   ), intent(inout)                   :: self
    type            (treeNode                      ), intent(inout)                   :: node
    type            (treeNode                      ), intent(inout), optional, target :: nodeHost
    logical                                         , intent(in   ), optional         :: atPericenter                   , includeCentrifugalAcceleration
    double precision                                , dimension(3  )                  :: tidalTensorEigenValueComponents
    double precision                                , dimension(3,3)                  :: tidalTensorComponents
    type            (tensorRank2Dimension3Symmetric)                                  :: tidalTensor
    type            (matrix                        )                                  :: tidalTensorMatrix              , tidalTensorEigenVectors
    type            (vector                        )                                  :: tidalTensorEigenValues
    logical                                                                           :: isSphericallySymmetric
    !![
    <optionalArgument name="atPericenter"                   defaultsTo=".false."/>
    <optionalArgument name="includeCentrifugalAcceleration" defaultsTo=".false."/>
    !!]

    ! For isolated halos, always return zero tidal field.
    if (.not.node%isSatellite().and..not.present(nodeHost)) then
       tidalTensorDominant=+0.0d0
       return
    end if
    ! Get the tidal tensor.
    tidalTensor=self%tidalTensorGet(node,nodeHost,atPericenter,includeCentrifugalAcceleration,isSphericallySymmetric=isSphericallySymmetric)
    ! Find the maximum of the tidal field over all directions in the satellite. To compute this we multiply a unit vector by
    ! the tidal tensor, which gives the vector tidal acceleration for displacements along that direction. We then take the dot
    ! product with that same unit vector to get the acceleration in that direction. This acceleration is maximized for a vector
    ! coinciding with the eigenvector of the tidal tensor with the largest eigenvalue. Since the acceleration in the direction
    ! of the eigenvector is exactly the eigenvalue corresponding to that eigenvector, we simply take the largest eigenvalue.
    ! For spherical mass distributions this reduces to:
    !
    ! +2GM(r)r⁻³ - 4πGρ(r)
    if (isSphericallySymmetric) then
       ! For spherically-symmetric mass distributions we can avoid the expense of solving for the eigenvalues. We simply
       ! evaluate the tidal tensor at [r,0,0] (since the distribution is spherically-symmetric we can evaluate at any
       ! position on the sphere), and take the 0,0 element of the tensor which will be the largest (and only) positive
       ! eigenvector.
       tidalTensorDominant            =tidalTensor%element(0,0)
    else
       ! Fully general calculation.
       tidalTensorComponents          =tidalTensor
       tidalTensorMatrix              =tidalTensorComponents
       call tidalTensorMatrix%eigenSystem(tidalTensorEigenVectors,tidalTensorEigenValues)
       tidalTensorEigenValueComponents=tidalTensorEigenValues
       tidalTensorDominant            =maxval(tidalTensorEigenValueComponents)
    end if
    return
  end function standardTidalTensorDominant
  
