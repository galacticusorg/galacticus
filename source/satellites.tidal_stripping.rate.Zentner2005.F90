!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implementation of a satellite tidal stripping class which follows the model of \cite{zentner_physics_2005}.

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Kind_Numbers           , only : kind_int8

  !# <satelliteTidalStripping name="satelliteTidalStrippingZentner2005">
  !#  <description>
  !#   A satellite tidal stripping class which uses the formalism of \cite{zentner_physics_2005} to compute the mass loss rate
  !#   $\dot{M}_\mathrm{sat}$:
  !#   \begin{equation}
  !#   \dot{M}_\mathrm{sat}=-\alpha \frac{M_\mathrm{sat}(>r_\mathrm{tidal})}{T_\mathrm{orb}},
  !#   \end{equation}
  !#   where $\alpha=${\normalfont \ttfamily [efficiency]},
  !#   \begin{equation}
  !#   T_\mathrm{orb} = {1 \over \hbox{max}(\omega/2\pi,v_\mathrm{r}/r)},
  !#   \end{equation}
  !#   where $\omega$ is the angular velocity of the satellite, $v_\mathrm{r}$ is the radial velocity, $r$ is the orbital radius,
  !#   and $r_\mathrm{tidal}$ is the tidal radius of the satellite, given by the \cite{king_structure_1962} formula:
  !#   \begin{equation}
  !#   r_\mathrm{tidal}=\left(\frac{GM_\mathrm{sat}}{\omega^2-d^2\Phi/dr^2}\right)^{1/3},
  !#   \end{equation}
  !#   where $\omega$ is the orbital angular velocity of the satellite, and $\Phi(r)$ is the gravitational potential due to the
  !#   host.
  !#  </description>
  !# </satelliteTidalStripping>
  type, extends(satelliteTidalStrippingClass) :: satelliteTidalStrippingZentner2005
     !% Implementation of a satellite tidal stripping class which follows the model of \cite{zentner_physics_2005}.
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: efficiency                    , expandMultiplier, &
          &                                                 radiusTidalPrevious
     integer         (kind_int8               )          :: lastUniqueID
   contains
     !# <methods>
     !#   <method description="Reset memoized calculations." method="calculationReset" />
     !# </methods>
     final     ::                     zentner2005Destructor
     procedure :: autoHook         => zentner2005AutoHook
     procedure :: calculationReset => zentner2005CalculationReset
     procedure :: massLossRate     => zentner2005MassLossRate
  end type satelliteTidalStrippingZentner2005

  interface satelliteTidalStrippingZentner2005
     !% Constructors for the {\normalfont \ttfamily zentner2005} satellite tidal stripping class.
     module procedure zentner2005ConstructorParameters
     module procedure zentner2005ConstructorInternal
  end interface satelliteTidalStrippingZentner2005

  ! Module-scope objects used for root finding.
  type            (treeNode), pointer :: zentner2005Node
  double precision                    :: zentner2005TidalPull
  !$omp threadprivate(zentner2005Node,zentner2005TidalPull)

contains

  function zentner2005ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily zentner2005} satellite tidaql stripping class which builds the object from a parameter set.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteTidalStrippingZentner2005)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass          ), pointer       :: darkMatterHaloScale_
    double precision                                                    :: efficiency

    !# <inputParameter>
    !#   <name>efficiency</name>
    !#   <defaultValue>2.5d0</defaultValue>
    !#   <description>The dimensionless rate coefficient apeparing in the \cite{zentner_physics_2005} expression for the tidal mass loss rate from subhalos.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=satelliteTidalStrippingZentner2005(efficiency,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_"/>
    return
  end function zentner2005ConstructorParameters

  function zentner2005ConstructorInternal(efficiency,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily zentner2005} satellite tidal stripping class.
    implicit none
    type            (satelliteTidalStrippingZentner2005)                        :: self
    class           (darkMatterHaloScaleClass          ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                    , intent(in)            :: efficiency
    !# <constructorAssign variables="efficiency, *darkMatterHaloScale_"/>

    self%expandMultiplier=2.0d0
    return
  end function zentner2005ConstructorInternal

  subroutine zentner2005AutoHook(self)
    !% Attach to the calculation reset event.
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(satelliteTidalStrippingZentner2005), intent(inout) :: self

    call calculationResetEvent%attach(self,zentner2005CalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine zentner2005AutoHook

  subroutine zentner2005Destructor(self)
    !% Destructor for the {\normalfont \ttfamily zentner2005} satellite tidal stripping class.
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(satelliteTidalStrippingZentner2005), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    call calculationResetEvent%detach(self,zentner2005CalculationReset)
    return
  end subroutine zentner2005Destructor

  double precision function zentner2005MassLossRate(self,node)
    !% Return a mass loss rate for satellites due to tidal stripping using the formulation of \cite{zentner_physics_2005}. To
    !% allow for non-spherical mass distributions, we proceed as follows to determine the tidal field:
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
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Galactic_Structure_Options        , only : coordinateSystemCartesian
    use :: Galactic_Structure_Tidal_Tensors  , only : Galactic_Structure_Tidal_Tensor
    use :: Galacticus_Nodes                  , only : nodeComponentSatellite          , treeNode
    use :: Linear_Algebra                    , only : vector                          , matrix                        , assignment(=)
    use :: Numerical_Constants_Astronomical  , only : gigaYear                        , megaParsec                    , gravitationalConstantGalacticus
    use :: Numerical_Constants_Math          , only : Pi
    use :: Numerical_Constants_Prefixes      , only : kilo
    use :: Root_Finder                       , only : rangeExpandMultiplicative       , rangeExpandSignExpectNegative , rangeExpandSignExpectPositive   , rootFinder
    use :: Tensors                           , only : assignment(=)                   , tensorRank2Dimension3Symmetric
    use :: Vectors                           , only : Vector_Magnitude                , Vector_Product
    implicit none
    class           (satelliteTidalStrippingZentner2005), intent(inout)         :: self
    type            (treeNode                          ), intent(inout), target :: node
    type            (treeNode                          ), pointer               :: nodeHost
    class           (nodeComponentSatellite            ), pointer               :: satellite
    double precision                                    , dimension(3  )        :: position                               , velocity, &
         &                                                                         tidalTensorEigenValueComponents
    double precision                                    , dimension(3,3)        :: tidalTensorComponents                  , tidalTensorEigenVectorComponents
    double precision                                    , parameter             :: toleranceAbsolute               =0.0d0 , toleranceRelative               =1.0d-3
    double precision                                    , parameter             :: radiusZero                      =0.0d0
    double precision                                    , parameter             :: radiusTidalTinyFraction         =1.0d-6
    type            (rootFinder                        ), save                  :: finder
    !$omp threadprivate(finder)
    double precision                                                            :: massSatellite                          , frequencyAngular                       , &
         &                                                                         periodOrbital                          , radius                                 , &
         &                                                                         tidalFieldRadial                       , radiusTidal                            , &
         &                                                                         massOuterSatellite                     , frequencyRadial
    type            (tensorRank2Dimension3Symmetric    )                        :: tidalTensor
    type            (matrix                            )                        :: tidalTensorMatrix                      , tidalTensorEigenVectors
    type            (vector                            )                        :: tidalTensorEigenValues

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
    frequencyRadial   = +abs             (   Dot_Product(position,velocity)) &
         &              /radius**2                                           &
         &              *kilo                                                &
         &              *gigaYear                                            &
         &              /megaParsec
    ! Find the orbital period. We use the larger of the angular and radial frequencies to avoid numerical problems for purely
    ! radial or purely circular orbits.
    periodOrbital     = 2.0d0*Pi/max(frequencyAngular,frequencyRadial)
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
       zentner2005TidalPull=  frequencyAngular**2-tidalFieldRadial
       if (self%radiusTidalPrevious <= 0.0d0) then
          self%radiusTidalPrevious=+(                                 &
               &                     +gravitationalConstantGalacticus &
               &                     *massSatellite                   &
               &                     /zentner2005TidalPull            &
               &                     *(kilo*gigaYear/megaParsec)**2   &
               &                    )**(1.0d0/3.0d0)
          self%expandMultiplier   =+2.0d0
       end if
       ! Check if tidal radius will lie outside of current boundary.
       if (Galactic_Structure_Enclosed_Mass(node,self%radiusTidalPrevious) >= massSatellite) then
          ! Tidal radius lies outside current boundary, so no additional stripping will occur.
          zentner2005MassLossRate=0.0d0
       else
          ! Find the tidal radius in the dark matter profile.
          if (.not.finder%isInitialized()) then
             call finder%rootFunction(zentner2005TidalRadiusSolver)
             call finder%tolerance   (toleranceAbsolute,toleranceRelative)
          end if
        call finder%rangeExpand(                                                             &
               &                rangeExpandUpward            =1.0d0*self%expandMultiplier  , &
               &                rangeExpandDownward          =1.0d0/self%expandMultiplier  , &
               &                rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                rangeExpandType              =rangeExpandMultiplicative      &
               &               )
          zentner2005Node => node
          ! Check for extremes.
          if (zentner2005TidalRadiusSolver(radiusTidalTinyFraction*self%radiusTidalPrevious) >  0.0d0) then
             ! Complete stripping.
             radiusTidal       =0.0d0
             massOuterSatellite=massSatellite
          else
             ! Find the tidal radius, using the previous result as an initial guess.
             radiusTidal       =finder%find(rootGuess=self%radiusTidalPrevious)
             massOuterSatellite=max(                                                                  &
                  &                 massSatellite-Galactic_Structure_Enclosed_Mass(node,radiusTidal), &
                  &                 0.0d0                                                             &
                  &                )
             self%radiusTidalPrevious  =radiusTidal
             self%expandMultiplier=1.2d0
          end if
          zentner2005MassLossRate=-self%efficiency    &
               &                  *massOuterSatellite &
               &                  /periodOrbital
       end if
    else
       zentner2005MassLossRate=0.0d0
    end if
    return
  end function zentner2005MassLossRate

  double precision function zentner2005TidalRadiusSolver(radius)
    !% Root function used to find the tidal radius within a subhalo.
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Numerical_Constants_Astronomical  , only : gigaYear                        , megaParsec
    use :: Numerical_Constants_Astronomical      , only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Prefixes      , only : kilo
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: enclosedMass

    ! Get the satellite component.
    enclosedMass                =+Galactic_Structure_Enclosed_Mass(zentner2005Node,radius)
    zentner2005TidalRadiusSolver=+zentner2005TidalPull               &
         &                       -gravitationalConstantGalacticus    &
         &                       *enclosedMass                       &
         &                       /radius                         **3 &
         &                       *(                                  &
         &                         +kilo                             &
         &                         *gigaYear                         &
         &                         /megaParsec                       &
         &                        )                              **2
    return
  end function zentner2005TidalRadiusSolver

  subroutine zentner2005CalculationReset(self,node)
    !% Reset the stored tidal radii.
    implicit none
    class(satelliteTidalStrippingZentner2005), intent(inout) :: self
    type (treeNode                          ), intent(inout) :: node

    self%radiusTidalPrevious=-1.0d0
    self%lastUniqueID       =node%uniqueID()
    return
  end subroutine zentner2005CalculationReset
