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
  Implements a node operator class that updates the properties of prompt cusps following the model of
  \cite{delos_cusp-halo_2025}.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <nodeOperator name="nodeOperatorDarkMatterProfilePromptCuspsUpdate">
   <description>
    A node operator class that updates the properties of prompt cusps following the model of \cite{delos_cusp-halo_2025}.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDarkMatterProfilePromptCuspsUpdate
     !!{     
     A node operator class that updates the properties of prompt cusps following the model of \cite{delos_cusp-halo_2025}.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     integer                                  :: promptCuspMassID              , promptCuspAmplitudeID, &
          &                                      promptCuspNFWYID
   contains
     final     ::                                        darkMatterProfilePromptCuspsUpdateDestructor
     procedure :: differentialEvolutionSolveAnalytics => darkMatterProfilePromptCuspsUpdateSolveAnalytics
  end type nodeOperatorDarkMatterProfilePromptCuspsUpdate
  
  interface nodeOperatorDarkMatterProfilePromptCuspsUpdate
     !!{
     Constructors for the \refClass{nodeOperatorDarkMatterProfilePromptCuspsUpdate} node operator class.
     !!}
     module procedure darkMatterProfilePromptCuspsUpdateConstructorParameters
     module procedure darkMatterProfilePromptCuspsUpdateConstructorInternal
  end interface nodeOperatorDarkMatterProfilePromptCuspsUpdate

  ! Submodule-scope variables used in root-finding.
  double precision            :: massHalo_        , amplitudeCusp_, &
       &                         concentration_   , radiusScale_
  !$omp threadprivate(massHalo_,concentration_,radiusScale_,amplitudeCusp_)

  ! Maximum allowed value of the y-parameter in the cusp-NFW profile. Values of 1 or greater are not valid. We limit here to a
  ! value close to 1.
  double precision, parameter :: yMaximum =0.999d0
  
contains
  
  function darkMatterProfilePromptCuspsUpdateConstructorParameters(parameters) result(self)
    !!{    
    Constructor for the \refClass{nodeOperatorDarkMatterProfilePromptCuspsUpdate} node operator class which takes a parameter set as
    input.    
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorDarkMatterProfilePromptCuspsUpdate)                :: self
    type (inputParameters                               ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass                      ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodeOperatorDarkMatterProfilePromptCuspsUpdate(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function darkMatterProfilePromptCuspsUpdateConstructorParameters

  function darkMatterProfilePromptCuspsUpdateConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDarkMatterProfilePromptCuspsUpdate} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorDarkMatterProfilePromptCuspsUpdate)                        :: self
    class(darkMatterHaloScaleClass                      ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    !![
    <addMetaProperty component="darkMatterProfile" name="promptCuspAmplitude" id="self%promptCuspAmplitudeID" isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspMass"      id="self%promptCuspMassID"      isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspNFWY"      id="self%promptCuspNFWYID"      isEvolvable="no" isCreator="no"/>
    !!]
    return
  end function darkMatterProfilePromptCuspsUpdateConstructorInternal

  subroutine darkMatterProfilePromptCuspsUpdateDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorDarkMatterProfilePromptCuspsUpdate} dark matter halo profile scale radius class.
    !!}
    implicit none
    type(nodeOperatorDarkMatterProfilePromptCuspsUpdate), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine darkMatterProfilePromptCuspsUpdateDestructor

  subroutine darkMatterProfilePromptCuspsUpdateSolveAnalytics(self,node,time)
    !!{
    Compute the value of the $y$-parameter in the prompt cusp.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic, nodeComponentDarkMatterProfile
    use :: Numerical_Constants_Math, only : Pi
    use :: Root_Finder             , only : rootFinder        , rangeExpandMultiplicative     , rangeExpandSignExpectPositive, rangeExpandSignExpectNegative
    implicit none
    class           (nodeOperatorDarkMatterProfilePromptCuspsUpdate), intent(inout) :: self
    type            (treeNode                                      ), intent(inout) :: node
    double precision                                                , intent(in   ) :: time
    class           (nodeComponentBasic                            ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile                ), pointer       :: darkMatterProfile
    type            (rootFinder                                    ), save          :: finder
    logical                                                         , save          :: finderInitialized=.false.
    !$omp threadprivate(finder,finderInitialized)
    double precision                                                                :: densityScale             , y

    ! Initialize the root finder.
    if (.not.finderInitialized) then
       finder=rootFinder(                                                             &
            &            rootFunction                 =densityNormalizationRoot     , &
            &            toleranceRelative            =1.0d-6                       , &
            &            rangeExpandUpward            =2.0d+0                       , &
            &            rangeExpandDownward          =0.5d+0                       , &
            &            rangeExpandType              =rangeExpandMultiplicative    , &
            &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &            rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive  &
            &           )
      finderInitialized=.true.
    end if
    ! Solve for the density normalization.
    basic               =>  node                                  %basic                    (                          )
    darkMatterProfile   =>  node                                  %darkMatterProfile        (                          )
    radiusScale_        =   darkMatterProfile                     %scale                    (                          )
    concentration_      =  +self             %darkMatterHaloScale_%radiusVirial             (node                      ) &
         &                 /                                       radiusScale_
    massHalo_           =   basic                                 %mass                     (                          )
    amplitudeCusp_      =  +darkMatterProfile                     %floatRank0MetaPropertyGet(self%promptCuspAmplitudeID)
    y                   =  +darkMatterProfile                     %floatRank0MetaPropertyGet(self%promptCuspNFWYID     )
    densityScale        =  +massHalo_                                                                                      &
         &                 /radiusScale_**3                                                                                &
         &                 /4.0d0                                                                                          &
         &                 /Pi                                                                                             &
         &                 /(                                                                                              &
         &                   + 2.0d0                       *asinh(sqrt(concentration_            )/                y     ) &
         &                   -(2.0d0-y**2)/sqrt(1.0d0-y**2)*atanh(sqrt(concentration_*(1.0d0-y**2)/(concentration_+y**2))) &
         &                   -sqrt(concentration_*(concentration_+y**2))/(1.0d0+concentration_)                            &
         &                  )
    densityScale        =   finder%find(rootGuess=densityScale)
    ! Compute the cusp y parameter.
    y                 =min(                      &
         &                 +amplitudeCusp_       &
         &                 /densityScale         &
         &                 /radiusScale_**1.5d0, &
         &                 +yMaximum             &
         &                )
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWYID,y)
    return
  end subroutine darkMatterProfilePromptCuspsUpdateSolveAnalytics

  double precision function densityNormalizationRoot(densityScale)
    !!{
    Root function used in finding the density normalization for cusp-NFW density profiles.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: densityScale
    double precision                :: y

    y                 =min(                      &
         &                 +amplitudeCusp_       &
         &                 /densityScale         &
         &                 /radiusScale_**1.5d0, &
         &                 +yMaximum             &
         &                )
    densityNormalizationRoot=+massHalo_                                                                                      &
         &                   /radiusScale_**3                                                                                &
         &                   /4.0d0                                                                                          &
         &                   /Pi                                                                                             &
         &                   /(                                                                                              &
         &                     + 2.0d0                       *asinh(sqrt(concentration_            )/                y     ) &
         &                     -(2.0d0-y**2)/sqrt(1.0d0-y**2)*atanh(sqrt(concentration_*(1.0d0-y**2)/(concentration_+y**2))) &
         &                     -sqrt(concentration_*(concentration_+y**2))/(1.0d0+concentration_)                            &
         &                    )                                                                                              &
         &                   -densityScale
    return
  end function densityNormalizationRoot
  
