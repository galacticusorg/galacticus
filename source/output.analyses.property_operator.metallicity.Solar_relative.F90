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
  Implements a property operator class which converts a metallicity, assumed to be a mass ratio of a
  given element to hydrogen, to $[\mathrm{N}/\mathrm{H}]$ form.
  !!}

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorMetallicitySolarRelative">
   <description>A property operator class which converts a metallicity, assumed to be a mass ratio of a given element to hydrogen, to $[\mathrm{N}/\mathrm{H}]$ form.</description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorMetallicitySolarRelative
     !!{
     A metallicity property operator class which converts a metallicity, assumed to be a mass ratio of a given element to
     hydrogen, to $[\mathrm{N}/\mathrm{H}]$ form.
     !!}
     private
     integer :: atomicNumberElement, atomicNumberHydrogen, &
          &     abundancesIndex
   contains
     procedure :: operate => metallicitySolarRelativeOperate
  end type outputAnalysisPropertyOperatorMetallicitySolarRelative

  interface outputAnalysisPropertyOperatorMetallicitySolarRelative
     !!{
     Constructors for the \refClass{outputAnalysisPropertyOperatorMetallicitySolarRelative} output analysis class.
     !!}
     module procedure metallicitySolarRelativeConstructorParameters
     module procedure metallicitySolarRelativeConstructorInternal
  end interface outputAnalysisPropertyOperatorMetallicitySolarRelative

contains

  function metallicitySolarRelativeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisPropertyOperatorMetallicitySolarRelative} output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (outputAnalysisPropertyOperatorMetallicitySolarRelative)                :: self
    type   (inputParameters                                       ), intent(inout) :: parameters
    integer                                                                        :: atomicNumberElement

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>atomicNumberElement</name>
      <source>parameters</source>
      <description>The atomic number of the element used to define metallicity.</description>
    </inputParameter>
    !!]
    self=outputAnalysisPropertyOperatorMetallicitySolarRelative(atomicNumberElement)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function metallicitySolarRelativeConstructorParameters

  function metallicitySolarRelativeConstructorInternal(atomicNumberElement) result (self)
    !!{
    Internal constructor for the \refClass{outputAnalysisPropertyOperatorMetallicitySolarRelative} output analysis distribution operator class.
    !!}
    use :: Atomic_Data, only : Abundance_Pattern_Lookup
    implicit none
    type   (outputAnalysisPropertyOperatorMetallicitySolarRelative)                :: self
    integer                                                        , intent(in   ) :: atomicNumberElement
    !![
    <constructorAssign variables="atomicNumberElement"/>
    !!]

    self%atomicNumberHydrogen=1
    self%abundancesIndex     =Abundance_Pattern_Lookup(abundanceName='solar')
    return
  end function metallicitySolarRelativeConstructorInternal

  double precision function metallicitySolarRelativeOperate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an metallicity output analysis property operator.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: Atomic_Data  , only : Atomic_Abundance
    implicit none
    class           (outputAnalysisPropertyOperatorMetallicitySolarRelative), intent(inout)           :: self
    double precision                                                        , intent(in   )           :: propertyValue
    type            (treeNode                                              ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType             ), intent(inout), optional :: propertyType
    integer         (c_size_t                                              ), intent(in   ), optional :: outputIndex
    double precision                                                                                  :: metallicityRelativeToSolar, massFractionSolar
    !$GLC attributes unused :: propertyType, outputIndex, node

    massFractionSolar         =+Atomic_Abundance(self%abundancesIndex,atomicNumber=self%atomicNumberElement ) &
         &                     /Atomic_Abundance(self%abundancesIndex,atomicNumber=self%atomicNumberHydrogen)
    metallicityRelativeToSolar=+propertyValue     &
         &                     /massFractionSolar
    if (metallicityRelativeToSolar > 0.0d0) then
       metallicitySolarRelativeOperate=+log10(metallicityRelativeToSolar)
    else
       metallicitySolarRelativeOperate=-huge (+0.0d0                    )
    end if
    return
  end function metallicitySolarRelativeOperate
