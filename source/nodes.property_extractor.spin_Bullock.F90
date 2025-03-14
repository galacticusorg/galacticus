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
Implements a node property extractor class for the \cite{bullock_profiles_2001} definition of spin parameter.
!!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale, darkMatterHaloScaleClass

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
     class  (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     logical                                    :: vectorSpinAvailable
     integer                                    :: elementCount_
   contains
     final     ::                 spinBullockDestructor
     procedure :: elementCount => spinBullockElementCount
     procedure :: extract      => spinBullockExtract
     procedure :: names        => spinBullockNames
     procedure :: descriptions => spinBullockDescriptions
     procedure :: unitsInSI    => spinBullockUnitsInSI
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

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodePropertyExtractorSpinBullock(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function spinBullockConstructorParameters

  function spinBullockConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily spinBullock} output analysis property extractor class.
    !!}
    use :: Galacticus_Nodes, only : defaultSpinComponent
    implicit none
    type (nodePropertyExtractorSpinBullock)                        :: self
    class(darkMatterHaloScaleClass        ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    self%vectorSpinAvailable=defaultSpinComponent%angularMomentumVectorIsGettable()
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
    <objectDestructor name="self%darkMatterHaloScale_"/>
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
    Implement extraction of the spin parameter under the \cite{bullock_profiles_2001} definition.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSpin, treeNode
    implicit none
    double precision                                  , dimension(:) , allocatable :: spinBullockExtract
    class           (nodePropertyExtractorSpinBullock), intent(inout), target      :: self
    type            (treeNode                        ), intent(inout), target      :: node
    double precision                                  , intent(in   )              :: time
    type            (multiCounter                    ), intent(inout), optional    :: instance
    class           (nodeComponentBasic              ), pointer                    :: basic
    class           (nodeComponentSpin               ), pointer                    :: spin
    double precision                                  , dimension(3)               :: spinVectorUnit
    double precision                                                               :: spinBullock
    !$GLC attributes unused :: time, instance

    allocate(spinBullockExtract(self%elementCount_))
    basic       =>  node                     %basic          (    )
    spin        =>  node                     %spin           (    )
    spinBullock =  +spin                     %angularMomentum(    ) &
         &         /sqrt(2.0d0)                                     &
         &         /basic                    %mass           (    ) &
         &         /self%darkMatterHaloScale_%radiusVirial   (node) &
         &         /self%darkMatterHaloScale_%velocityVirial (node)
    spinBullockExtract(1) =   spinBullock
    if (self%vectorSpinAvailable) then
       spinVectorUnit          =  +spin%angularMomentumVector() &
            &                     /spin%angularMomentum      ()
       spinBullockExtract(2:4) =  +spinBullock                  &
            &                     *spinVectorUnit
    end if
    return
  end function spinBullockExtract

  subroutine spinBullockNames(self,time,names)
    !!{
    Return the name of the {\normalfont \ttfamily spinBullock} property.
    !!}
    implicit none
    class           (nodePropertyExtractorSpinBullock), intent(inout)                             :: self
    double precision                                  , intent(in   )                             :: time
    type            (varying_string                  ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(self%elementCount_))
    names    (1)=var_str('spinBullock' )
    if (self%vectorSpinAvailable) then
        names(2)=var_str('spinBullockX')
        names(3)=var_str('spinBullockY')
        names(4)=var_str('spinBullockZ')
     end if
    return
  end subroutine spinBullockNames

  subroutine spinBullockDescriptions(self,time,descriptions)
    !!{
    Return a description of the {\normalfont \ttfamily spinBullock} property.
    !!}
    implicit none
    class           (nodePropertyExtractorSpinBullock), intent(inout)                             :: self
    double precision                                  , intent(in   )                             :: time
    type            (varying_string                  ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(self%elementCount_))
    descriptions    (1)=var_str('Spin parameter of the halo under the Bullock et al. (2001) definition [].'                   )
    if (self%vectorSpinAvailable) then
        descriptions(2)=var_str('x-component of the spin parameter of the halo under the Bullock et al. (2001) definition [].')
        descriptions(3)=var_str('y-component of the spin parameter of the halo under the Bullock et al. (2001) definition [].')
        descriptions(4)=var_str('z-component of the spin parameter of the halo under the Bullock et al. (2001) definition [].')
     end if
    return
  end subroutine spinBullockDescriptions

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

