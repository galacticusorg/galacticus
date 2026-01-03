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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson, Xiaolong Du.

  !!{
  Implementation of a satellite tidal radius class which follows the method of \cite{king_structure_1962}.
  !!}

  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Kind_Numbers           , only : kind_int8

  !![
  <satelliteTidalStrippingRadius name="satelliteTidalStrippingRadiusKing1962">
   <description>
    A satellite tidal radius class which uses the method of \cite{king_structure_1962}:
    \begin{equation}
    r_\mathrm{tidal}=\left(\frac{GM_\mathrm{sat}}{\gamma_\mathrm{c} \omega^2-d^2\Phi/dr^2}\right)^{1/3},
    \end{equation}
    where $\omega$ is the orbital angular velocity of the satellite, $\Phi(r)$ is the gravitational potential due to the host,
    and $\gamma_\mathrm{c}=${\normalfont \ttfamily [efficiencyCentrifugal]} is the a model parameter that controls the efficiency of
    centrifugal force. The calculation is based on the dark matter only density profile of the satellite---no accounting is
    made for the baryonic components.
   </description>
  </satelliteTidalStrippingRadius>
  !!]
  type, extends(satelliteTidalStrippingRadiusClass) :: satelliteTidalStrippingRadiusKing1962
     !!{
     Implementation of a satellite tidal radius class which follows the method of \cite{king_structure_1962}.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_  => null()
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_  => null()
     double precision                                    :: efficiencyCentrifugal          , expandMultiplier, &
          &                                                 fractionDarkMatter
     logical                                             :: applyPreInfall
   contains
     final     ::           king1962Destructor
     procedure :: radius => king1962Radius
  end type satelliteTidalStrippingRadiusKing1962

  interface satelliteTidalStrippingRadiusKing1962
     !!{
     Constructors for the \refClass{satelliteTidalStrippingRadiusKing1962} satellite tidal stripping class.
     !!}
     module procedure king1962ConstructorParameters
     module procedure king1962ConstructorInternal
  end interface satelliteTidalStrippingRadiusKing1962

contains

  function king1962ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{satelliteTidalStrippingRadiusKing1962} satellite tidal stripping class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteTidalStrippingRadiusKing1962)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (cosmologyParametersClass             ), pointer       :: cosmologyParameters_
    class           (darkMatterHaloScaleClass             ), pointer       :: darkMatterHaloScale_
    double precision                                                       :: efficiencyCentrifugal
    logical                                                                :: applyPreInfall

    !![
    <inputParameter>
      <name>efficiencyCentrifugal</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>Efficiency of the centrifugal force. For zero value, the centrifugal force is neglected in computing the tidal radius.</description>
    </inputParameter>
    <inputParameter>
      <name>applyPreInfall</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, tidal radii are computed pre-infall.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=satelliteTidalStrippingRadiusKing1962(efficiencyCentrifugal,applyPreInfall,cosmologyParameters_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function king1962ConstructorParameters

  function king1962ConstructorInternal(efficiencyCentrifugal,applyPreInfall,cosmologyParameters_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{satelliteTidalStrippingRadiusKing1962} satellite tidal stripping class.
    !!}
    implicit none
    type            (satelliteTidalStrippingRadiusKing1962)                        :: self
    class           (cosmologyParametersClass             ), intent(in   ), target :: cosmologyParameters_
    class           (darkMatterHaloScaleClass             ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                       , intent(in   )         :: efficiencyCentrifugal
    logical                                                , intent(in   )         :: applyPreInfall
    double precision                                       , parameter             :: toleranceAbsolute   =0.0d0, toleranceRelative=1.0d-3
    !![
    <constructorAssign variables="efficiencyCentrifugal, applyPreInfall, *cosmologyParameters_, *darkMatterHaloScale_"/>
    !!]

    self%fractionDarkMatter=+(                                         & 
         &                    +self%cosmologyParameters_%OmegaMatter() &
         &                    -self%cosmologyParameters_%OmegaBaryon() &
         &                   )                                         &
         &                  /  self%cosmologyParameters_%OmegaMatter()
    return
  end function king1962ConstructorInternal

  subroutine king1962Destructor(self)
    !!{
    Destructor for the \refClass{satelliteTidalStrippingRadiusKing1962} satellite tidal stripping class.
    !!}
    implicit none
    type(satelliteTidalStrippingRadiusKing1962), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine king1962Destructor

  double precision function king1962Radius(self,node)
    !!{
    Return the tidal radius using the formulation of \cite{king_structure_1962}. To allow for non-spherical mass distributions,
    we proceed as follows to determine the tidal field:
    
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
    use :: Coordinates                     , only : assignment(=)            , coordinateCartesian
    use :: Galactic_Structure_Options      , only : massTypeDark
    use :: Galacticus_Nodes                , only : nodeComponentSatellite   , nodeComponentBasic            , treeNode
    use :: Linear_Algebra                  , only : assignment(=)            , matrix                        , vector
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Constants_Astronomical, only : gigaYear                 , gravitationalConstant_internal, megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Tensors                         , only : assignment(=)            , tensorRank2Dimension3Symmetric
    use :: Vectors                         , only : Vector_Magnitude         , Vector_Product
    implicit none
    class           (satelliteTidalStrippingRadiusKing1962), intent(inout), target :: self
    type            (treeNode                             ), intent(inout), target :: node
    type            (treeNode                             ), pointer               :: nodeHost
    class           (nodeComponentBasic                   ), pointer               :: basic                                  , basicHost
    class           (nodeComponentSatellite               ), pointer               :: satellite
    class           (massDistributionClass                ), pointer               :: massDistribution_                      , massDistributionDark
    double precision                                       , dimension(3  )        :: position                               , velocity                    , &
         &                                                                            tidalTensorEigenValueComponents
    double precision                                       , dimension(3,3)        :: tidalTensorComponents
    double precision                                       , parameter             :: radiusZero                      =0.0d+0
    double precision                                       , parameter             :: radiusTidalTinyFraction         =1.0d-6
    double precision                                                               :: massSatellite                          , frequencyAngular            , &
         &                                                                            radius                                 , tidalFieldRadial            , &
         &                                                                            radiusGuess                            , densityTidal                , &
         &                                                                            tidalPull                              , tidalTensorEigenValueMaximum, &
         &                                                                            radiusDownwardLimit
    type            (tensorRank2Dimension3Symmetric       )                        :: tidalTensor
    type            (matrix                               )                        :: tidalTensorMatrix                      , tidalTensorEigenVectors
    type            (vector                               )                        :: tidalTensorEigenValues
    type            (coordinateCartesian                  )                        :: coordinates

    ! Find the host node.
    if (node%isOnMainBranch().or.(.not.self%applyPreInfall.and..not.node%isSatellite())) then
       ! For nodes on the main branch, always return the virial radius.
       king1962Radius=self%darkMatterHaloScale_%radiusVirial(node)
       return
    else if (node%isSatellite()) then
       ! For satellites, use the node with which they will merge.
       nodeHost => node%mergesWith()
    else
       ! Walk up the branch to find the node with which this node will merge.
       nodeHost => node
       do while (nodeHost%isPrimaryProgenitor())
          nodeHost => nodeHost%parent
       end do
       nodeHost  => nodeHost%parent%firstChild
       ! Follow that node back through descendants until we find the node at the corresponding time.
       basic     => node    %basic()
       basicHost => nodeHost%basic()
       do while (basicHost%time() > basic%time())
          if (.not.associated(nodeHost%firstChild)) exit
          nodeHost  => nodeHost%firstChild
          basicHost => nodeHost%basic     ()
       end do
    end if
    ! Get required quantities from satellite and host nodes.
    satellite         => node     %satellite (        )
    massSatellite     =  satellite%boundMass (        )
    position          =  satellite%position  (        )
    velocity          =  satellite%velocity  (        )
    radius            =  Vector_Magnitude    (position)
    ! Compute the orbital period.
    frequencyAngular  = +Vector_Magnitude(Vector_Product(position,velocity)) &
         &              /radius**2                                           &
         &              *kilo                                                &
         &              *gigaYear                                            &
         &              /megaParsec
    ! Find the maximum of the tidal field over all directions in the satellite. To compute this we multiply a unit vector by
    ! the tidal tensor, which gives the vector tidal acceleration for displacements along that direction. We then take the dot
    ! product with that same unit vector to get the acceleration in that direction. This acceleration is maximized for a vector
    ! coinciding with the eigenvector of the tidal tensor with the largest eigenvalue. Since the acceleration in the direction
    ! of the eigenvector is exactly the eigenvalue corresponding to that eigenvector, we simply take the largest eigenvalue.
    ! For spherical mass distributions this reduces to:
    !
    ! -2GM(r)r⁻³ - 4πGρ(r)
    if (associated(nodeHost)) then
       massDistribution_ => nodeHost%massDistribution()
       if (massDistribution_%isSphericallySymmetric()) then
          ! For spherically-symmetric mass distributions we can avoid the expense of solving for the eigenvalues. We simply
          ! evaluate the tidal tensor at [r,0,0] (since the distribution is spherically-symmetric we can evaluate at any
          ! position on the sphere), and take the 0,0 element of the tensor which will be the largest (and only) positive
          ! eigenvector.
          coordinates                    =[radius,0.0d0,0.0d0]
          tidalTensor                    =massDistribution_%tidalTensor(coordinates)
          tidalTensorEigenValueMaximum   =tidalTensor%element(0,0)
       else
          coordinates                    =position
          tidalTensor                    =massDistribution_%tidalTensor(coordinates)
          tidalTensorComponents          =tidalTensor
          tidalTensorMatrix              =tidalTensorComponents
          call tidalTensorMatrix%eigenSystem(tidalTensorEigenVectors,tidalTensorEigenValues)
          tidalTensorEigenValueComponents=tidalTensorEigenValues
          tidalTensorEigenValueMaximum   =maxval(tidalTensorEigenValueComponents)
       end if       
       tidalFieldRadial=-tidalTensorEigenValueMaximum &
            &           *(                            &
            &             +kilo                       &
            &             *gigaYear                   &
            &             /megaParsec                 &
            &            )**2
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    else
       tidalFieldRadial                =+0.0d0
    end if
    ! If the tidal force is stretching (not compressing), compute the tidal radius.
    massDistribution_ => node%massDistribution()
    tidalPull         =  self%efficiencyCentrifugal*frequencyAngular**2-tidalFieldRadial
    if     (                                                             &
         &   tidalPull                                          >  0.0d0 &
         &  .and.                                                        &
         &   massSatellite                                      >  0.0d0 &
         &  .and.                                                        &
         &   massDistribution_%massEnclosedBySphere(radiusZero) >= 0.0d0 &
         & ) then
       ! Find the tidal density.
       densityTidal=+tidalPull                         &
            &       /(kilo*gigaYear/megaParsec)    **2 &
            &       /gravitationalConstant_internal    &
            &       *3.0d0                             &
            &       /4.0d0                             &
            &       /Pi
       ! Solve for the radius enclosing this density.
       radiusGuess       =  +sqrt(                                              &
               &                  +gravitationalConstant_internal               &
               &                  *massSatellite                                &
               &                  /self%darkMatterHaloScale_%radiusVirial(node) &
               &                  /tidalPull                                    &
               &                  *(kilo*gigaYear/megaParsec)**2                &
               &                 )
       radiusDownwardLimit=radiusTidalTinyFraction*radiusGuess
       if   (                                                                 &
          &  + 3.0                                                            &
          &   /4.0d0                                                          &
          &   /Pi                                                             &
          &   *massDistribution_%massEnclosedBySphere(radiusDownwardLimit)    &
          &   /                                       radiusDownwardLimit **3 &
          &  >                                                                &
          &    densityTidal                                                   &
          & ) then
          king1962Radius    =  massDistribution_%radiusEnclosingDensity(densityTidal,radiusGuess)
       else
          king1962Radius    =  0.0d0
       end if
    else
       ! If the bound mass of the satellite exceeds the original mass (which can happen during failed ODE steps), simply return
       ! the virial radius. Otherwise, solve for the radius enclosing the current bound mass.
       massDistributionDark => node%massDistribution(massType=massTypeDark)
       if (massSatellite > massDistributionDark%massEnclosedBySphere(self%darkMatterHaloScale_%radiusVirial(node))) then
          king1962Radius=self                %darkMatterHaloScale_%radiusVirial       (node                                 )
       else
          king1962Radius=massDistributionDark                     %radiusEnclosingMass(massSatellite*self%fractionDarkMatter)
       end if
       !![
       <objectDestructor name="massDistributionDark"/>
       !!]
    end if
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function king1962Radius
