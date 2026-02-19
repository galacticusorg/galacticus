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
Implements a cosmological angular distance corrector analysis property operator class.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Output_Times       , only : outputTimesClass

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorCsmlgyAnglrDstnc">
   <description>
   An output analysis property operator class which corrects properties for the difference in cosmological angular diameter
   distance between true and assumed (i.e. in the observational analysis) cosmologies. Typically the observational data will have
   been analyzed assuming some specific set of cosmological parameters which will differ from that in the current
   model. Therefore, the size of a galaxy must be adjusted to match what would be inferred if they were assessed using the same
   cosmological parameters as were used for the observational data. Typically, this will mean that sizes are scaled in proportion
   to $D^\prime_\mathrm{A}(z)/D_\mathrm{A}(z)$, where $D_\mathrm{A}(z)$ and $D^\prime_\mathrm{A}(z)$ are the luminosity distances
   to redshift $z$ in the true and assumed cosmologies respectively.
  </description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorCsmlgyAnglrDstnc
     !!{
     A cosmological angular distance corrector analysis property operator class.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctionsModel => null(), cosmologyFunctionsData => null()
     class           (outputTimesClass       ), pointer                   :: outputTimes_            => null()
     double precision                         , allocatable, dimension(:) :: correctionFactor
   contains
     final     ::            csmlgyAngularDistanceDestructor
     procedure :: operate => csmlgyAngularDistanceOperate
  end type outputAnalysisPropertyOperatorCsmlgyAnglrDstnc

  interface outputAnalysisPropertyOperatorCsmlgyAnglrDstnc
     !!{
     Constructors for the \refClass{outputAnalysisPropertyOperatorCsmlgyAnglrDstnc} output analysis class.
     !!}
     module procedure csmlgyAngularDistanceConstructorParameters
     module procedure csmlgyAngularDistanceConstructorInternal
  end interface outputAnalysisPropertyOperatorCsmlgyAnglrDstnc

contains

  function csmlgyAngularDistanceConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisPropertyOperatorCsmlgyAnglrDstnc} output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (outputAnalysisPropertyOperatorCsmlgyAnglrDstnc)                :: self
    type (inputParameters                               ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                       ), pointer       :: cosmologyFunctionsModel, cosmologyFunctionsData
    class(outputTimesClass                              ), pointer       :: outputTimes_
    type (inputParameters                               )                :: dataAnalysisParameters

    ! Check and read parameters.
    dataAnalysisParameters=parameters%subParameters('dataAnalysis',requirePresent=.false.,requireValue=.false.)
    !![
    <objectBuilder class="outputTimes"        name="outputTimes_"            source="parameters"            />
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctionsModel" source="parameters"            />
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctionsData"  source="dataAnalysisParameters"/>
    !!]
    ! Construct the object.
    self=outputAnalysisPropertyOperatorCsmlgyAnglrDstnc(cosmologyFunctionsModel,cosmologyFunctionsData,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"           />
    <objectDestructor name="cosmologyFunctionsModel"/>
    <objectDestructor name="cosmologyFunctionsData" />
    !!]
    return
  end function csmlgyAngularDistanceConstructorParameters

  function csmlgyAngularDistanceConstructorInternal(cosmologyFunctionsModel,cosmologyFunctionsData,outputTimes_) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisPropertyOperatorCsmlgyAnglrDstnc} output analysis property operator class.
    !!}
    use            :: Error            , only : Error_Report
    use, intrinsic :: ISO_C_Binding    , only : c_size_t
    implicit none
    type            (outputAnalysisPropertyOperatorCsmlgyAnglrDstnc)                        :: self
    class           (cosmologyFunctionsClass                       ), intent(in   ), target :: cosmologyFunctionsModel       , cosmologyFunctionsData
    class           (outputTimesClass                              ), intent(in   ), target :: outputTimes_
    double precision                                                , parameter             :: distanceSmall          =1.0d-6
    integer         (c_size_t                                      )                        :: outputIndex
    double precision                                                                        :: redshift                      , timeData              , &
         &                                                                                     timeModel                     , distanceData          , &
         &                                                                                     distanceModel
    !![
    <constructorAssign variables="*cosmologyFunctionsModel, *cosmologyFunctionsData, *outputTimes_"/>
    !!]

    allocate(self%correctionFactor(self%outputTimes_%count()))
    do outputIndex=1,self%outputTimes_%count()
       ! Get current redshift.
       redshift        =self%outputTimes_%redshift(outputIndex)
       ! Find corresponding cosmic times in both the data and model cosmological models.
       timeData        =self%cosmologyFunctionsData %cosmicTime                 (          &
            &           self%cosmologyFunctionsData %expansionFactorFromRedshift (         &
            &                                                                     redshift &
            &                                                                    )         &
            &                                                                   )
       timeModel       =self%cosmologyFunctionsModel%cosmicTime                 (          &
            &           self%cosmologyFunctionsModel%expansionFactorFromRedshift (         &
            &                                                                     redshift &
            &                                                                    )         &
            &                                                                   )
       ! Find angular distance in both data and model cosmological models.
       distanceData    =self%cosmologyFunctionsData %distanceAngular(timeData )
       distanceModel   =self%cosmologyFunctionsModel%distanceAngular(timeModel)
       ! Compute the correction factor - the assumption here is that the property was derived from an angular size. Therefore, we
       ! first find the true physical size observed in the model by multiplying by the model angular distance, and then convert
       ! back to the property of interest under the assumptions of the data analysis by dividing by the angular distance in the
       ! data analysis cosmological model.  Handle zero distance (i.e. present day outputs - and in fact we test for very small
       ! rather than zero distance to catch rounding errors) as a special case.
       if (distanceModel < distanceSmall) then
          if (distanceData > distanceSmall) call Error_Report('angular distance in model cosmology is zero, but non-zero for data cosmology - correction is undefined'//{introspection:location})
          self%correctionFactor(outputIndex)= +1.0d0
       else
          self%correctionFactor(outputIndex)=+distanceModel &
               &                             /distanceData
       end if
    end do
    return
  end function csmlgyAngularDistanceConstructorInternal

  subroutine csmlgyAngularDistanceDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisPropertyOperatorCsmlgyAnglrDstnc} output analysis property operator class.
    !!}
    implicit none
    type(outputAnalysisPropertyOperatorCsmlgyAnglrDstnc), intent(inout) :: self

    if (allocated(self%correctionFactor)) deallocate(self%correctionFactor)
    !![
    <objectDestructor name="self%cosmologyFunctionsModel"/>
    <objectDestructor name="self%cosmologyFunctionsData" />
    <objectDestructor name="self%outputTimes_"           />
    !!]
    return
  end subroutine csmlgyAngularDistanceDestructor

  double precision function csmlgyAngularDistanceOperate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an csmlgyAngularDistance output analysis property operator.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (outputAnalysisPropertyOperatorCsmlgyAnglrDstnc), intent(inout)           :: self
    double precision                                                , intent(in   )           :: propertyValue
    type            (treeNode                                      ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType     ), intent(inout), optional :: propertyType
    integer         (c_size_t                                      ), intent(in   ), optional :: outputIndex
    !$GLC attributes unused :: propertyType, node

    ! Validate.
    if (.not.present(outputIndex)) call Error_Report('ouputIndex is required'//{introspection:location})
    ! Apply the correction.
    csmlgyAngularDistanceOperate=+propertyValue                      &
         &                       *self%correctionFactor(outputIndex)
    return
  end function csmlgyAngularDistanceOperate
