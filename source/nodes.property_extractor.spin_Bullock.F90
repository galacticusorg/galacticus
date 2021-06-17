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

!!{
Contains a module which implements a node property extractor class for the \cite{bullock_profiles_2001} definition of spin parameter.
!!}

  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScale , darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMO, darkMatterProfileDMOClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorSpinBullock">
   <description>A node property extractor class for the \cite{bullock_profiles_2001} definition of spin parameter.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorSpinBullock
     !!{
     A property extractor class for the \cite{bullock_profiles_2001} definition of spin parameter..
     !!}
     private
     class  (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_
     class  (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_
     logical                                     :: vectorSpinAvailable
     integer                                     :: elementCount_
   contains
     final     ::                 spinBullockDestructor
     procedure :: elementCount => spinBullockElementCount
     procedure :: extract      => spinBullockExtract
     procedure :: names        => spinBullockNames
     procedure :: descriptions => spinBullockDescriptions
     procedure :: unitsInSI    => spinBullockUnitsInSI
     procedure :: type         => spinBullockType
  end type nodePropertyExtractorSpinBullock

  interface nodePropertyExtractorSpinBullock
     !!{
     Constructors for the ``spinBullock'' output analysis class.
     !!}
     module procedure spinBullockConstructorParameters
     module procedure spinBullockConstructorInternal
  end interface nodePropertyExtractorSpinBullock

contains

  function spinBullockConstructorParameters(parameters) result(self)
    !!{
    Constructor for the y output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorSpinBullock)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass        ), pointer       :: darkMatterHaloScale_
    class(darkMatterProfileDMOClass       ), pointer       :: darkMatterProfileDMO_

    !![
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=nodePropertyExtractorSpinBullock(darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function spinBullockConstructorParameters

  function spinBullockConstructorInternal(darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily spinBullock} output analysis property extractor class.
    !!}
    use :: Galacticus_Nodes, only : defaultSpinComponent
    implicit none
    type (nodePropertyExtractorSpinBullock)                        :: self
    class(darkMatterHaloScaleClass        ), intent(in   ), target :: darkMatterHaloScale_
    class(darkMatterProfileDMOClass       ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="*darkMatterHaloScale_, *darkMatterProfileDMO_"/>
    !!]

    self%vectorSpinAvailable=defaultSpinComponent%spinVectorIsGettable()
    if (self%vectorSpinAvailable) then
       self%elementCount_=4
    else
       self%elementCount_=1
    end if
    return
  end function spinBullockConstructorInternal

  subroutine spinBullockDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily spinBullock} output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorSpinBullock), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine spinBullockDestructor

  integer function spinBullockElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily spinBullock} property extractor.
    !!}
    implicit none
    class           (nodePropertyExtractorSpinBullock), intent(inout) :: self
    double precision                                  , intent(in   ) :: time
    !$GLC attributes unused :: time

    spinBullockElementCount=self%elementCount_
    return
  end function spinBullockElementCount

  function spinBullockExtract(self,node,time,instance)
    !!{
    Implement extraction of the spin parameter under the \cite{bullock_profiles_2001} defintiion.
    !!}
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum
    use :: Galacticus_Nodes      , only : nodeComponentBasic               , nodeComponentSpin, treeNode
    implicit none
    double precision                                   , dimension(:) , allocatable :: spinBullockExtract
    class           (nodePropertyExtractorSpinBullock ), intent(inout), target      :: self
    type            (treeNode                         ), intent(inout), target      :: node
    double precision                                   , intent(in   )              :: time
    type            (multiCounter                     ), intent(inout), optional    :: instance
    class           (nodeComponentBasic               ), pointer                    :: basic
    class           (nodeComponentSpin                ), pointer                    :: spin
    double precision                                   , dimension(3)               :: spinVectorUnit
    double precision                                                                :: spinBullock
    !$GLC attributes unused :: time, instance

    allocate(spinBullockExtract(self%elementCount_))
    basic                 =>  node                     %basic         (                               )
    spinBullock           =  +Dark_Matter_Halo_Angular_Momentum       (node,self%darkMatterProfileDMO_) &
         &                   /sqrt(2.0d0)                                                               &
         &                   /basic                    %mass          (                               ) &
         &                   /self%darkMatterHaloScale_%virialRadius  (node                           ) &
         &                   /self%darkMatterHaloScale_%virialVelocity(node                           )
    spinBullockExtract(1) =   spinBullock
    if (self%vectorSpinAvailable) then
       spin                    =>  node          %spin      ()
       spinVectorUnit          =  +spin          %spinVector() &
            &                     /spin          %spin      ()
       spinBullockExtract(2:4) =  +spinBullock                 &
            &                     *spinVectorUnit
    end if
    return
  end function spinBullockExtract

  function spinBullockNames(self,time)
    !!{
    Return the name of the {\normalfont \ttfamily spinBullock} property.
    !!}
    implicit none
    type            (varying_string                  ), dimension(:) , allocatable :: spinBullockNames
    class           (nodePropertyExtractorSpinBullock), intent(inout)              :: self
    double precision                                  , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(spinBullockNames(self%elementCount_))
    spinBullockNames(1)         = var_str('spinBullock' )
    if (self%vectorSpinAvailable)                                                                          &
         & spinBullockNames(2:4)=[var_str('spinBullockX'),var_str('spinBullockY'),var_str('spinBullockZ')]
    return
  end function spinBullockNames

  function spinBullockDescriptions(self,time)
    !!{
    Return a description of the {\normalfont \ttfamily spinBullock} property.
    !!}
    implicit none
    type            (varying_string                  ), dimension(:) , allocatable :: spinBullockDescriptions
    class           (nodePropertyExtractorSpinBullock), intent(inout)              :: self
    double precision                                  , intent(in   )              :: time
    !$GLC attributes unused :: time

    allocate(spinBullockDescriptions(self%elementCount_))
    spinBullockDescriptions       (1  )= var_str('Spin parameter of the halo under the Bullock et al. (2001) definition [].'                   )
    if (self%vectorSpinAvailable)                                                                                                                 &
         & spinBullockDescriptions(2:4)=[                                                                                                         &
         &                               var_str('x-component of the spin parameter of the halo under the Bullock et al. (2001) definition [].'), &
         &                               var_str('y-component of the spin parameter of the halo under the Bullock et al. (2001) definition [].'), &
         &                               var_str('z-component of the spin parameter of the halo under the Bullock et al. (2001) definition [].')  &
         &                              ]
    return
  end function spinBullockDescriptions

  function spinBullockUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily spinBullock} property in the SI system.
    !!}
    implicit none
    double precision                                  , allocatable  , dimension(:) :: spinBullockUnitsInSI
    class           (nodePropertyExtractorSpinBullock), intent(inout)               :: self
    double precision                                  , intent(in   )               :: time
    !$GLC attributes unused :: self, time

    allocate(spinBullockUnitsInSI(self%elementCount_))
    spinBullockUnitsInSI=0.0d0
    return
  end function spinBullockUnitsInSI

  integer function spinBullockType(self)
    !!{
    Return the type of the {\normalfont \ttfamily spinBullock} property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorSpinBullock), intent(inout) :: self
    !$GLC attributes unused :: self

    spinBullockType=outputAnalysisPropertyTypeLinear
    return
  end function spinBullockType
