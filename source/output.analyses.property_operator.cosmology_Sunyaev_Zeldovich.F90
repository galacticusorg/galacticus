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
Implements a thermal Sunyaev-Zeldovich cosmological scaling corrector analysis property operator class.
!!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass
  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Output_Times        , only : outputTimesClass

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorCosmologySZ">
   <description>
    An output analysis property operator class which scales thermal Sunyaev-Zeldovich properties to $z=0$ using the expected
    cosmological scalings. Specifically, the property is multiplied by a factor of $E^{-2/3}(t)$ where $H(t) = E(t) H_0$ is the
    epoch-dependent Hubble parameter \citep{planck_collaboration_planck_2013}.
  </description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorCosmologySZ
     !!{
     A cosmological angular distance corrector analysis property operator class.
     !!}
     private
     class           (cosmologyFunctionsClass ), pointer                   :: cosmologyFunctions_  => null()
     class           (cosmologyParametersClass), pointer                   :: cosmologyParameters_ => null()
     class           (outputTimesClass        ), pointer                   :: outputTimes_         => null()
     double precision                          , allocatable, dimension(:) :: correctionFactor
   contains
     final     ::            csmlgySZDestructor
     procedure :: operate => csmlgySZOperate
  end type outputAnalysisPropertyOperatorCosmologySZ

  interface outputAnalysisPropertyOperatorCosmologySZ
     !!{
     Constructors for the \refClass{outputAnalysisPropertyOperatorCosmologySZ} output analysis class.
     !!}
     module procedure csmlgySZConstructorParameters
     module procedure csmlgySZConstructorInternal
  end interface outputAnalysisPropertyOperatorCosmologySZ

contains

  function csmlgySZConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisPropertyOperatorCosmologySZ} output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (outputAnalysisPropertyOperatorCosmologySZ)                :: self
    type (inputParameters                          ), intent(inout) :: parameters
    class(cosmologyParametersClass                 ), pointer       :: cosmologyParameters_
    class(cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_
    class(outputTimesClass                         ), pointer       :: outputTimes_

    !![
    <objectBuilder class="outputTimes"         name="outputTimes_"         source="parameters"/>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    ! Construct the object.
    self=outputAnalysisPropertyOperatorCosmologySZ(cosmologyParameters_,cosmologyFunctions_,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"        />
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function csmlgySZConstructorParameters

  function csmlgySZConstructorInternal(cosmologyParameters_,cosmologyFunctions_,outputTimes_) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisPropertyOperatorCosmologySZ} output analysis property operator class.
    !!}
    use, intrinsic :: ISO_C_Binding    , only : c_size_t
    implicit none
    type   (outputAnalysisPropertyOperatorCosmologySZ)                        :: self
    class  (cosmologyParametersClass                 ), intent(in   ), target :: cosmologyParameters_
    class  (cosmologyFunctionsClass                  ), intent(in   ), target :: cosmologyFunctions_
    class  (outputTimesClass                         ), intent(in   ), target :: outputTimes_
    integer(c_size_t                                 )                        :: outputIndex
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologyFunctions_, *outputTimes_"/>
    !!]

    allocate(self%correctionFactor(self%outputTimes_%count()))
    do outputIndex=1,self%outputTimes_%count()
       self%correctionFactor(outputIndex)=(                                                                                       &
            &                              +self%cosmologyFunctions_ %HubbleParameterEpochal(self%outputTimes_%time(outputIndex)) &
            &                              /self%cosmologyParameters_%HubbleConstant        (                                   ) &
            &                             )**(-2.0d0/3.0d0)
    end do
    return
  end function csmlgySZConstructorInternal

  subroutine csmlgySZDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisPropertyOperatorCosmologySZ} output analysis property operator class.
    !!}
    implicit none
    type(outputAnalysisPropertyOperatorCosmologySZ), intent(inout) :: self

    deallocate(self%correctionFactor)
    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%outputTimes_"        />
    !!]
    return
  end subroutine csmlgySZDestructor

  double precision function csmlgySZOperate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an csmlgySZ output analysis property operator.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (outputAnalysisPropertyOperatorCosmologySZ), intent(inout)           :: self
    double precision                                           , intent(in   )           :: propertyValue
    type            (treeNode                                 ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType), intent(inout), optional :: propertyType
    integer         (c_size_t                                 ), intent(in   ), optional :: outputIndex
    !$GLC attributes unused :: propertyType, node

    ! Validate.
    if (.not.present(outputIndex)) call Error_Report('ouputIndex is required'//{introspection:location})
    ! Apply the correction.
    csmlgySZOperate=+propertyValue                      &
         &          *self%correctionFactor(outputIndex)
    return
  end function csmlgySZOperate
