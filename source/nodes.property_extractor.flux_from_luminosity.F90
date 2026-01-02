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
Implements a class which extracts fluxes from luminosities.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorFluxFromLuminosity">
   <description>A property extractor that converts luminosities to fluxes.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorFluxFromLuminosity
     !!{
     A flux-from-luminosity property extractor class.
     !!}
     private
     class(cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_    => null()
     class(nodePropertyExtractorClass), pointer :: nodePropertyExtractor_ => null()
   contains
     final     ::                 fluxFromLuminosityDestructor
     procedure :: elementCount => fluxFromLuminosityElementCount
     procedure :: extract      => fluxFromLuminosityExtract
     procedure :: quantity     => fluxFromLuminosityQuantity
     procedure :: names        => fluxFromLuminosityNames
     procedure :: descriptions => fluxFromLuminosityDescriptions
     procedure :: unitsInSI    => fluxFromLuminosityUnitsInSI
  end type nodePropertyExtractorFluxFromLuminosity

  interface nodePropertyExtractorFluxFromLuminosity
     !!{
     Constructors for the \refClass{nodePropertyExtractorFluxFromLuminosity} output analysis class.
     !!}
     module procedure fluxFromLuminosityConstructorParameters
     module procedure fluxFromLuminosityConstructorInternal
  end interface nodePropertyExtractorFluxFromLuminosity

contains

  function fluxFromLuminosityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorFluxFromLuminosity} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorFluxFromLuminosity)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(nodePropertyExtractorClass             ), pointer       :: nodePropertyExtractor_
    class(cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    !!]
    self=nodePropertyExtractorFluxFromLuminosity(nodePropertyExtractor_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"/>
    <objectDestructor name="cosmologyFunctions_"   />
    !!]
    return
  end function fluxFromLuminosityConstructorParameters

  function fluxFromLuminosityConstructorInternal(nodePropertyExtractor_,cosmologyFunctions_) result(self)
    implicit none
    type (nodePropertyExtractorFluxFromLuminosity)                        :: self
    class(nodePropertyExtractorClass             ), intent(in   ), target :: nodePropertyExtractor_
    class(cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*nodePropertyExtractor_, *cosmologyFunctions_"/>
    !!]

    return
  end function fluxFromLuminosityConstructorInternal
  
  subroutine fluxFromLuminosityDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorFluxFromLuminosity} output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorFluxFromLuminosity), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"/>
    <objectDestructor name="self%cosmologyFunctions_"   />
    !!]
    return
  end subroutine fluxFromLuminosityDestructor

  integer function fluxFromLuminosityElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily fluxFromLuminosity} property extractors.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (nodePropertyExtractorFluxFromLuminosity), intent(inout) :: self
    double precision                                         , intent(in   ) :: time

    select type (nodePropertyExtractor_=> self%nodePropertyExtractor_)
    class is (nodePropertyExtractorTuple)
       fluxFromLuminosityElementCount=nodePropertyExtractor_%elementCount(time)
    class default
       fluxFromLuminosityElementCount=0
       call Error_Report('unsupported class'//{introspection:location})
    end select
    return
  end function fluxFromLuminosityElementCount

  function fluxFromLuminosityExtract(self,node,time,instance) result(fluxes)
    !!{
    Extract fluxes from luminosities.
    !!}
    use :: Numerical_Constants_Units       , only : ergs
    use :: Numerical_Constants_Prefixes    , only : centi
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Error                           , only : Error_Report
    implicit none
    double precision                                         , allocatable   , dimension(:) :: fluxes
    double precision                                         , allocatable   , dimension(:) :: unitsInSI
    class           (nodePropertyExtractorFluxFromLuminosity), intent(inout) , target       :: self
    type            (treeNode                               ), intent(inout) , target       :: node
    double precision                                         , intent(in   )                :: time
    type            (multiCounter                           ), intent(inout) , optional     :: instance

    select type (nodePropertyExtractor_=> self%nodePropertyExtractor_)
    class is (nodePropertyExtractorTuple)
       fluxes   =nodePropertyExtractor_%extract  (node,time,instance)
       unitsInSI=nodePropertyExtractor_%unitsInSI(     time         )
       if (self%cosmologyFunctions_%expansionFactor(time) < 1.0d0) then
          ! Epoch is z>0 - compute the flux.
          fluxes=+fluxes                                               &
               & *unitsInSI                                            &
               & *centi**2                                             &
               & /ergs                                                 &
               & /4.0d0                                                &
               & /Pi                                                   &
               & /self%cosmologyFunctions_%distanceLuminosity(time)**2 &
               & /megaParsec**2
       else
          ! Epoch is z=0 or later - no flux is defined.
          fluxes=huge(0.0d0)
       end if
    class default
       allocate(fluxes(0))
       call Error_Report('unsupported class'//{introspection:location})
    end select
    return    
  end function fluxFromLuminosityExtract

  function fluxFromLuminosityQuantity(self)
    !!{
    Return the class of the flux property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyQuantityUnknown
    implicit none
    type (enumerationOutputAnalysisPropertyQuantityType)                :: fluxFromLuminosityQuantity
    class(nodePropertyExtractorFluxFromLuminosity      ), intent(inout) :: self
    !$GLC attributes unused :: self

    fluxFromLuminosityQuantity=outputAnalysisPropertyQuantityUnknown
    return
  end function fluxFromLuminosityQuantity

  subroutine fluxFromLuminosityNames(self,time,names)
    !!{
    Return the name of the fluxFromLuminosity property.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (nodePropertyExtractorFluxFromLuminosity), intent(inout)                            :: self
    double precision                                         , intent(in   )                            :: time
    type            (varying_string                         ), intent(inout), dimension(:), allocatable :: names
    integer                                                                                             :: i
    
    select type (nodePropertyExtractor_=> self%nodePropertyExtractor_)
    class is (nodePropertyExtractorTuple)
       call nodePropertyExtractor_%names(time,names)
       do i=1,size(names)
          if (extract(names(i),1,10) == "luminosity") then
             names(i)="flux"//extract(names(i),11)
          else
             call Error_Report('can not construct property name'//{introspection:location})
          end if
       end do
    class default
       allocate(names(0))
       call Error_Report('unsupported class'//{introspection:location})
    end select    
    return
  end subroutine fluxFromLuminosityNames

  subroutine fluxFromLuminosityDescriptions(self,time,descriptions)
    !!{
    Return a description of the fluxFromLuminosity property.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (nodePropertyExtractorFluxFromLuminosity), intent(inout)                            :: self
    double precision                                         , intent(in   )                            :: time
    type            (varying_string                         ), intent(inout), dimension(:), allocatable :: descriptions
    integer                                                                                             :: i

    select type (nodePropertyExtractor_=> self%nodePropertyExtractor_)
    class is (nodePropertyExtractorTuple)
       call nodePropertyExtractor_%descriptions(time,descriptions)
       do i=1,size(descriptions)
          descriptions(i)="Flux (in ergs s¯¹ cm¯²) derived from: "//descriptions(i)
       end do
    class default
       allocate(descriptions(0))
       call Error_Report('unsupported class'//{introspection:location})
    end select
    return
  end subroutine fluxFromLuminosityDescriptions

  function fluxFromLuminosityUnitsInSI(self,time) result(unitsInSI)
    !!{
    Return the units of the fluxFromLuminosity property in the SI system.
    !!}
    use :: Numerical_Constants_Prefixes, only : centi
    use :: Numerical_Constants_Units   , only : ergs
    implicit none
    double precision                                         , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorFluxFromLuminosity), intent(inout)              :: self
    double precision                                         , intent(in   )              :: time

    allocate(unitsInSI(self%elementCount(time)))
    unitsInSI=ergs/centi**2
    return
  end function fluxFromLuminosityUnitsInSI
