!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

  !% Implementation of a satellite tidal radius class which follows the method of \cite{king_structure_1962}.

  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Kind_Numbers           , only : kind_int8
  use :: Root_Finder            , only : rootFinder

  !# <satelliteTidalStrippingRadius name="satelliteTidalStrippingRadiusKing1962">
  !#  <description>
  !#   A satellite tidal radius class which uses the method of \cite{king_structure_1962}:
  !#   \begin{equation}
  !#   r_\mathrm{tidal}=\left(\frac{GM_\mathrm{sat}}{\omega^2-d^2\Phi/dr^2}\right)^{1/3},
  !#   \end{equation}
  !#   where $\omega$ is the orbital angular velocity of the satellite, and $\Phi(r)$ is the gravitational potential due to the
  !#   host. The calculation is based on the dark matter only density profile of the satellite---no accounting is made for the
  !#   baryonic components.
  !#  </description>
  !# </satelliteTidalStrippingRadius>
  type, extends(satelliteTidalStrippingRadiusClass) :: satelliteTidalStrippingRadiusKing1962
     !% Implementation of a satellite tidal radius class which follows the method of \cite{king_structure_1962}.
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: radiusTidalPrevious           , expandMultiplier, &
          &                                                 fractionDarkMatter
     integer         (kind_int8               )          :: lastUniqueID
     type            (rootFinder              )          :: finder
   contains
     !# <methods>
     !#   <method description="Reset memoized calculations." method="calculationReset" />
     !# </methods>
     final     ::                     king1962Destructor
     procedure :: autoHook         => king1962AutoHook
     procedure :: calculationReset => king1962CalculationReset
     procedure :: radius           => king1962Radius
  end type satelliteTidalStrippingRadiusKing1962

  interface satelliteTidalStrippingRadiusKing1962
     !% Constructors for the {\normalfont \ttfamily king1962} satellite tidal stripping class.
     module procedure king1962ConstructorParameters
     module procedure king1962ConstructorInternal
  end interface satelliteTidalStrippingRadiusKing1962

  ! Module-scope objects used for root finding.
  type            (treeNode), pointer :: king1962Node
  double precision                    :: king1962TidalPull
  !$omp threadprivate(king1962Node,king1962TidalPull)

contains

  function king1962ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily king1962} satellite tidaql stripping class which builds the object from a parameter set.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (satelliteTidalStrippingRadiusKing1962)                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(cosmologyParametersClass             ), pointer       :: cosmologyParameters_
    class(darkMatterHaloScaleClass             ), pointer       :: darkMatterHaloScale_

    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=satelliteTidalStrippingRadiusKing1962(cosmologyParameters_,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="darkMatterHaloScale_"/>
    return
  end function king1962ConstructorParameters

  function king1962ConstructorInternal(cosmologyParameters_,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily king1962} satellite tidal stripping class.
    implicit none
    type            (satelliteTidalStrippingRadiusKing1962)                        :: self
    class           (cosmologyParametersClass             ), intent(in   ), target :: cosmologyParameters_
    class           (darkMatterHaloScaleClass             ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                       , parameter             :: toleranceAbsolute   =0.0d0, toleranceRelative=1.0d-3
    !# <constructorAssign variables="*cosmologyParameters_, *darkMatterHaloScale_"/>

    self%fractionDarkMatter=+(                                         & 
         &                    +self%cosmologyParameters_%OmegaMatter() &
         &                    -self%cosmologyParameters_%OmegaBaryon() &
         &                   )                                         &
         &                  /  self%cosmologyParameters_%OmegaMatter()
    self%expandMultiplier=2.0d0
    self%finder          =rootFinder(                                             &
         &                           rootFunction     =king1962TidalRadiusSolver, &
         &                           toleranceAbsolute=toleranceAbsolute        , &
         &                           toleranceRelative=toleranceRelative          &
         &                          )
    return
  end function king1962ConstructorInternal

  subroutine king1962AutoHook(self)
    !% Attach to the calculation reset event.
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(satelliteTidalStrippingRadiusKing1962), intent(inout) :: self

    call calculationResetEvent%attach(self,king1962CalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine king1962AutoHook

  subroutine king1962Destructor(self)
    !% Destructor for the {\normalfont \ttfamily king1962} satellite tidal stripping class.
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(satelliteTidalStrippingRadiusKing1962), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    call calculationResetEvent%detach(self,king1962CalculationReset)
    return
  end subroutine king1962Destructor

  double precision function king1962Radius(self,node)
    !% Return the tidal radius using the formulation of \cite{king_structure_1962}. To allow for non-spherical mass distributions,
    !% we proceed as follows to determine the tidal field:
    !%
    !% Let $\bm{\mathsf{G}}$ be the gravitational tidal tensor evaluated at the position of the satellite. Consider a unit vector,
    !% $\boldsymbol{\hat{x}}$ in the satellite. The tidal field along this vector is $\bm{\mathsf{G}} \boldsymbol{\hat{x}}$. The
    !% radial component of the tidal field in this direction is then $\boldsymbol{\hat{x}} \bm{\mathsf{G}}
    !% \boldsymbol{\hat{x}}$. We want to find the maximum of the tidal field over all possible directions (i.e. all possible unit
    !% vectors).
    !%
    !% We can write our unit vector as
    !% \begin{equation}
    !%  \boldsymbol{\hat{x}} = \sum_{i=1}^3 a_i \boldsymbol{\hat{e}}_i,
    !% \end{equation}
    !% where the $\boldsymbol{\hat{e}}_i$ are the eigenvectors of $\bm{\mathsf{G}}$ and $\sum_{i=1}^3 a_i^2 = 1$. Then since, by
    !% definition, $\bm{\mathsf{G}} \boldsymbol{\hat{e}}_i = \lambda_i \boldsymbol{\hat{e}}_i$, where $\lambda_i$ are the
    !% eigenvalues of $\bm{\mathsf{G}}$, we have that
    !% \begin{equation}
    !%  \boldsymbol{\hat{x}} \bm{\mathsf{G}} \boldsymbol{\hat{x}}= \sum_{i=1}^3 a_i^2 \lambda_i.
    !% \end{equation}
    !% The sum on the right hand side of the above is a weighted average of eigenvalues. Any weighted averge is maximized by
    !% setting the weight of the largest value to $1$, and all other weights to $0$. Therefore, our tidal field is maximized along
    !% the direction corresponding the eigenvector of $\bm{\mathsf{G}}$ with the largest eigenvalue. (Note that we want the
    !% largest positive eigenvalue, not the largest absolute eigenvalue as we're interested in stretching tidal fields, not
    !% compressive ones.)
    use :: Galacticus_Error                  , only : Galacticus_Error_Report         , errorStatusSuccess
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass, Galactic_Structure_Radius_Enclosing_Mass
    use :: Galactic_Structure_Options        , only : coordinateSystemCartesian       , massTypeDark
    use :: Galactic_Structure_Tidal_Tensors  , only : Galactic_Structure_Tidal_Tensor
    use :: Galacticus_Nodes                  , only : nodeComponentSatellite          , nodeComponentBasic                      , treeNode
    use :: Linear_Algebra                    , only : vector                          , matrix                                  , assignment(=)
    use :: Numerical_Constants_Astronomical  , only : gigaYear                        , megaParsec                              , gravitationalConstantGalacticus
    use :: Numerical_Constants_Math          , only : Pi
    use :: Numerical_Constants_Prefixes      , only : kilo
    use :: Root_Finder                       , only : rangeExpandMultiplicative       , rangeExpandSignExpectNegative           , rangeExpandSignExpectPositive
    use :: Tensors                           , only : assignment(=)                   , tensorRank2Dimension3Symmetric
    use :: Vectors                           , only : Vector_Magnitude                , Vector_Product
    implicit none
    class           (satelliteTidalStrippingRadiusKing1962), intent(inout)         :: self
    type            (treeNode                             ), intent(inout), target :: node
    type            (treeNode                             ), pointer               :: nodeHost
    class           (nodeComponentBasic                   ), pointer               :: basic
    class           (nodeComponentSatellite               ), pointer               :: satellite
    double precision                                       , dimension(3  )        :: position                               , velocity                        , &
         &                                                                            tidalTensorEigenValueComponents
    double precision                                       , dimension(3,3)        :: tidalTensorComponents                  , tidalTensorEigenVectorComponents
    double precision                                       , parameter             :: radiusZero                      =0.0d0
    double precision                                       , parameter             :: radiusTidalTinyFraction         =1.0d-6
    integer                                                                        :: status
    double precision                                                               :: massSatellite                          , frequencyAngular                , &
         &                                                                            radius                                 , tidalFieldRadial                , &
         &                                                                            radiusLimitDownward
    type            (tensorRank2Dimension3Symmetric       )                        :: tidalTensor
    type            (matrix                               )                        :: tidalTensorMatrix                      , tidalTensorEigenVectors
    type            (vector                               )                        :: tidalTensorEigenValues

    ! Test for a satellite node.
    if (.not.node%isSatellite()) then
       king1962Radius=self%darkMatterHaloScale_%virialRadius(node)
       return
    end if
    ! Get required quantities from satellite and host nodes.
    nodeHost          => node     %mergesWith(        )
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
    ! Find the maximum of the tidal field over all directions in the satellite. To compute this we multiply the a unit vector by
    ! the tidal tensor, which gives the vector tidal acceleration for displacements along that direction. We then take the dot
    ! product with that same unit vector to get the magnitude of this acceleration in the that direction. This magnitude is
    ! maximized for a vector coinciding with the eigenvector of the tidal tensor with the largest eigenvalue. For spherical mass
    ! distributions this reduces to:
    !
    ! -2GM(r)r⁻³ - 4πGρ(r)
    tidalTensor                     = Galactic_Structure_Tidal_Tensor(nodeHost,position)
    tidalTensorComponents           = tidalTensor
    tidalTensorMatrix               = tidalTensorComponents
    call tidalTensorMatrix%eigenSystem(tidalTensorEigenVectors,tidalTensorEigenValues)
    tidalTensorEigenValueComponents = tidalTensorEigenValues
    tidalTensorEigenVectorComponents= tidalTensorEigenVectors
    tidalFieldRadial                =-abs(tidalTensor%vectorProject(tidalTensorEigenVectorComponents(maxloc(tidalTensorEigenValueComponents,dim=1),:))) &
         &                           *(                                                                                                                 &
         &                             +kilo                                                                                                            &
         &                             *gigaYear                                                                                                        &
         &                             /megaParsec                                                                                                      &
         &                            )**2
    ! If the tidal force is stretching (not compressing), compute the tidal radius.
    if     (                                                                      &
         &   frequencyAngular**2                               > tidalFieldRadial &
         &  .and.                                                                 &
         &   massSatellite                                     >  0.0d0           &
         &  .and.                                                                 &
         &   Galactic_Structure_Enclosed_Mass(node,radiusZero) >= 0.0d0           &
         & ) then
       ! Check if node differs from previous one for which we performed calculations.
       if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
       ! Initial estimate of the tidal radius.
       king1962TidalPull=frequencyAngular**2-tidalFieldRadial
       if (self%radiusTidalPrevious <= 0.0d0) then
          self%radiusTidalPrevious=+(                                 &
               &                     +gravitationalConstantGalacticus &
               &                     *massSatellite                   &
               &                     /king1962TidalPull               &
               &                     *(kilo*gigaYear/megaParsec)**2   &
               &                    )**(1.0d0/3.0d0)
          self%expandMultiplier   =+2.0d0
       end if
       ! Find the tidal radius in the dark matter profile.
       radiusLimitDownward=+radiusTidalTinyFraction  &
            &              *self%radiusTidalPrevious
       call self%finder%rangeExpand(                                                              &
            &                       rangeExpandUpward            =+1.0d0*self%expandMultiplier  , &
            &                       rangeExpandDownward          =+1.0d0/self%expandMultiplier  , &
            &                       rangeExpandDownwardSignExpect= rangeExpandSignExpectNegative, &
            &                       rangeExpandUpwardSignExpect  = rangeExpandSignExpectPositive, &
            &                       rangeDownwardLimit           = radiusLimitDownward          , &
            &                       rangeExpandType              = rangeExpandMultiplicative      &
            &                      )
       king1962Node => node
       ! Find the tidal radius, using the previous result as an initial guess.
       self%radiusTidalPrevious=self%finder%find(rootGuess=self%radiusTidalPrevious,status=status)
       if (status == errorStatusSuccess) then
          self%expandMultiplier   =1.2d0
       else if (king1962TidalRadiusSolver(radiusLimitDownward) > 0.0d0) then 
          ! Complete stripping.
          self%radiusTidalPrevious=0.0d0
       else
          ! Find the tidal radius, using the previous result as an initial guess.
          call Galacticus_Error_Report('unable to find tidal radius'//{introspection:location})
       end if
       king1962Radius=self%radiusTidalPrevious
    else
       ! If the bound mass of the satellite exceeds the original mass (which can happen during failed ODE steps), simply return
       ! the virial radius. Otherwise, solve for the radius enclosing the current bound mass.
       basic => node%basic()
       if (massSatellite > basic%mass()) then
          king1962Radius=self%darkMatterHaloScale_%virialRadius(node)
       else
          king1962Radius=Galactic_Structure_Radius_Enclosing_Mass(node,massSatellite*self%fractionDarkMatter,massType=massTypeDark)
       end if
    end if
    return
  end function king1962Radius

  double precision function king1962TidalRadiusSolver(radius)
    !% Root function used to find the tidal radius within a subhalo.
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Numerical_Constants_Astronomical  , only : gravitationalConstantGalacticus , gigaYear, megaParsec  
    use :: Numerical_Constants_Prefixes      , only : kilo
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: enclosedMass

    ! Get the satellite component.
    enclosedMass             =+Galactic_Structure_Enclosed_Mass(king1962Node,radius)
    king1962TidalRadiusSolver=+king1962TidalPull                  &
         &                    -gravitationalConstantGalacticus    &
         &                    *enclosedMass                       &
         &                    /radius                         **3 &
         &                    *(                                  &
         &                      +kilo                             &
         &                      *gigaYear                         &
         &                      /megaParsec                       &
         &                     )                              **2
    return
  end function king1962TidalRadiusSolver

  subroutine king1962CalculationReset(self,node)
    !% Reset the stored tidal radii.
    implicit none
    class(satelliteTidalStrippingRadiusKing1962), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node

    self%radiusTidalPrevious=-1.0d0
    self%lastUniqueID       =node%uniqueID()
    return
  end subroutine king1962CalculationReset
