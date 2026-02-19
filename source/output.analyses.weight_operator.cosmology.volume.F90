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
Implements a cosmological volume corrector analysis weight operator class.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Geometry_Surveys   , only : surveyGeometryClass

  !![
  <outputAnalysisWeightOperator name="outputAnalysisWeightOperatorCsmlgyVolume">
   <description>
    An output analysis weight operator class which corrects weights for the difference in cosmological volume between true and
    assumed (i.e. in the observational analysis) cosmologies. Typically the observational data will have been analyzed assuming
    some specific set of cosmological parameters which will differ from that in the current model. Therefore, the comoving volume
    occupied by a population of galaxies must be adjusted to match what would be inferred if they were assessed using the same
    cosmological parameters as were used for the observational data. Typically, this will mean that weights are scaled in
    proportion to $V_\mathrm{max} / V^\prime_\mathrm{max}$, where $V_\mathrm{max}$ and $V^\prime_\mathrm{max}$ are the maximum
    volumes within which the galaxy would have been detected in the true and assumed cosmologies respectively.
   </description>
  </outputAnalysisWeightOperator>
  !!]
  type, extends(outputAnalysisWeightOperatorClass) :: outputAnalysisWeightOperatorCsmlgyVolume
     !!{
     A cosmological volume corrector analysis weight operator class.
     !!}
     private
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctionsModel => null(), cosmologyFunctionsData => null()
     class(surveyGeometryClass    ), pointer :: surveyGeometry_         => null()
   contains
     final     ::            csmlgyVolumeDestructor
     procedure :: operate => csmlgyVolumeOperate
  end type outputAnalysisWeightOperatorCsmlgyVolume

  interface outputAnalysisWeightOperatorCsmlgyVolume
     !!{
     Constructors for the \refClass{outputAnalysisWeightOperatorCsmlgyVolume} output analysis class.
     !!}
     module procedure csmlgyVolumeConstructorParameters
     module procedure csmlgyVolumeConstructorInternal
  end interface outputAnalysisWeightOperatorCsmlgyVolume

contains

  function csmlgyVolumeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisWeightOperatorCsmlgyVolume} output analysis weight operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (outputAnalysisWeightOperatorCsmlgyVolume)                :: self
    type   (inputParameters                         ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass                 ), pointer       :: cosmologyFunctionsModel, cosmologyFunctionsData
    class  (surveyGeometryClass                     ), pointer       :: surveyGeometry_
    type   (inputParameters                         )                :: dataAnalysisParameters

    ! Check and read parameters.
    dataAnalysisParameters=parameters%subParameters('dataAnalysis',requirePresent=.false.,requireValue=.false.)
    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctionsModel" source="parameters"            />
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctionsData"  source="dataAnalysisParameters"/>
    <objectBuilder class="surveyGeometry"     name="surveyGeometry_"         source="dataAnalysisParameters"/>
    !!]
    ! Construct the object.
    self=outputAnalysisWeightOperatorCsmlgyVolume(cosmologyFunctionsModel,cosmologyFunctionsData,surveyGeometry_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctionsModel"/>
    <objectDestructor name="cosmologyFunctionsData" />
    <objectDestructor name="surveyGeometry_"        />
    !!]
    return
  end function csmlgyVolumeConstructorParameters

  function csmlgyVolumeConstructorInternal(cosmologyFunctionsModel,cosmologyFunctionsData,surveyGeometry_) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisWeightOperatorCsmlgyVolume} output analysis weight operator class.
    !!}
    implicit none
    type   (outputAnalysisWeightOperatorCsmlgyVolume)                        :: self
    class  (cosmologyFunctionsClass                 ), intent(in   ), target :: cosmologyFunctionsModel, cosmologyFunctionsData
    class  (surveyGeometryClass                     ), intent(in   ), target :: surveyGeometry_
    !![
    <constructorAssign variables="*cosmologyFunctionsModel, *cosmologyFunctionsData, *surveyGeometry_"/>
    !!]

    return
  end function csmlgyVolumeConstructorInternal

  subroutine csmlgyVolumeDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisWeightOperatorCsmlgyVolume} output analysis weight operator class.
    !!}
    implicit none
    type(outputAnalysisWeightOperatorCsmlgyVolume), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctionsModel"/>
    <objectDestructor name="self%cosmologyFunctionsData" />
    <objectDestructor name="self%surveyGeometry_"        />
    !!]
    return
  end subroutine csmlgyVolumeDestructor

  double precision function csmlgyVolumeOperate(self,weightValue,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
    !!{
    Implement an csmlgyVolume output analysis weight operator.
    !!}
    use            :: Error                  , only : Error_Report
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Output_Analyses_Options, only : outputAnalysisPropertyQuantityLuminosity, outputAnalysisPropertyQuantityStarFormationRate,outputAnalysisPropertyQuantityMass, outputAnalysisPropertyTypeLinear, &
          &                                           outputAnalysisPropertyTypeMagnitude     , outputAnalysisPropertyTypeLog10
    implicit none
    class           (outputAnalysisWeightOperatorCsmlgyVolume     ), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    double precision                                               , intent(in   ) :: propertyValue         , propertyValueIntrinsic, &
         &                                                                            weightValue
    type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   ) :: propertyType
    type            (enumerationOutputAnalysisPropertyQuantityType), intent(in   ) :: propertyQuantity
    integer         (c_size_t                                     ), intent(in   ) :: outputIndex
    double precision                                                               :: distanceModelMinimum  , distanceDataMinimum   , &
         &                                                                            distanceModelMaximum  , distanceDataMaximum   , &
         &                                                                            expansionFactorMinimum, expansionFactorMaximum, &
         &                                                                            volumeData            , volumeModel           , &
         &                                                                            correctionFactor
    integer                                                                           field
    !$GLC attributes unused :: outputIndex, propertyValue, propertyType, node

    ! Compute the correction factor - the assumption here is that the volume density was derived from a 1/Vₘₐₓ type approach. To
    ! correct for the distance in model and data cosmological models, we therefore first multiply the weight by Vₘₐₓ for the model
    ! cosmology (which gives the absolute number of galaxies per survey), then divide by the Vₘₐₓ for the data cosmology.
    volumeData =0.0d0
    volumeModel=0.0d0
    do field=1,self%surveyGeometry_%fieldCount()
       select case (propertyQuantity%ID)
       case (outputAnalysisPropertyQuantityMass      %ID)
          select case (propertyType%ID)
          case (outputAnalysisPropertyTypeLinear   %ID)
             distanceDataMinimum   =+self%surveyGeometry_       %distanceMinimum       (                                           &
                  &                                                                      mass             =propertyValue         , &
                  &                                                                      field            =field                   &
                  &                                                                    )
             distanceDataMaximum   =+self%surveyGeometry_       %distanceMaximum       (                                           &
                  &                                                                      mass             =propertyValue         , &
                  &                                                                      field            =field                   &
                  &                                                                    )
          case (outputAnalysisPropertyTypeLog10    %ID)
             distanceDataMinimum   =+self%surveyGeometry_       %distanceMinimum       (                                           &
                  &                                                                      mass             =propertyValueIntrinsic, &
                  &                                                                      field            =field                   &
                  &                                                                    )
             distanceDataMaximum   =+self%surveyGeometry_       %distanceMaximum       (                                           &
                  &                                                                      mass             =propertyValueIntrinsic, &
                  &                                                                      field            =field                   &
                  &                                                                    )
          case default
             call Error_Report('unsupported property type'//{introspection:location})
          end select
       case (outputAnalysisPropertyQuantityLuminosity%ID)
          select case (propertyType%ID)
          case (outputAnalysisPropertyTypeLinear   %ID)
             distanceDataMinimum   =+self%surveyGeometry_       %distanceMinimum       (                                           &
                  &                                                                      luminosity       =propertyValue         , &
                  &                                                                      field            =field                   &
                  &                                                                    )
             distanceDataMaximum   =+self%surveyGeometry_       %distanceMaximum       (                                           &
                  &                                                                      luminosity       =propertyValue         , &
                  &                                                                      field            =field                   &
                  &                                                                    )
          case (outputAnalysisPropertyTypeLog10    %ID)
             distanceDataMinimum   =+self%surveyGeometry_       %distanceMinimum       (                                           &
                  &                                                                      luminosity       =propertyValueIntrinsic, &
                  &                                                                      field            =field                   &
                  &                                                                    )
             distanceDataMaximum   =+self%surveyGeometry_       %distanceMaximum       (                                           &
                  &                                                                      luminosity       =propertyValueIntrinsic, &
                  &                                                                      field            =field                   &
                  &                                                                    )
          case (outputAnalysisPropertyTypeMagnitude    %ID)
             distanceDataMinimum   =+self%surveyGeometry_       %distanceMinimum       (                                           &
                  &                                                                      magnitudeAbsolute=propertyValue         , &
                  &                                                                      field            =field                   &
                  &                                                                    )
             distanceDataMaximum   =+self%surveyGeometry_       %distanceMaximum       (                                           &
                  &                                                                      magnitudeAbsolute=propertyValue         , &
                  &                                                                      field            =field                   &
                  &                                                                    )
          case default
             call Error_Report('unsupported property type'//{introspection:location})
          end select
       case (outputAnalysisPropertyQuantityStarFormationRate%ID)
          select case (propertyType%ID)
          case (outputAnalysisPropertyTypeLinear   %ID)
             distanceDataMinimum   =+self%surveyGeometry_       %distanceMinimum       (                                           &
                  &                                                                      starFormationRate=propertyValue         , &
                  &                                                                      field            =field                   &
                  &                                                                    )
             distanceDataMaximum   =+self%surveyGeometry_       %distanceMaximum       (                                           &
                  &                                                                      starFormationRate=propertyValue         , &
                  &                                                                      field            =field                   &
                  &                                                                    )
          case (outputAnalysisPropertyTypeLog10    %ID)
             distanceDataMinimum   =+self%surveyGeometry_       %distanceMinimum       (                                           &
                  &                                                                      starFormationRate=propertyValueIntrinsic, &
                  &                                                                      field            =field                   &
                  &                                                                    )
             distanceDataMaximum   =+self%surveyGeometry_       %distanceMaximum       (                                           &
                  &                                                                      starFormationRate=propertyValueIntrinsic, &
                  &                                                                      field            =field                   &
                  &                                                                    )
          case default
             call Error_Report('unsupported property type'//{introspection:location})
          end select
       case default
          call Error_Report('unsupported property class'//{introspection:location})
       end select
       expansionFactorMinimum=+self%cosmologyFunctionsData%expansionFactor       (                         &
            &                  self%cosmologyFunctionsData%timeAtDistanceComoving (                        &
            &                                                                      distanceDataMinimum     &
            &                                                                     )                        &
            &                                                                    )
       expansionFactorMaximum=+self%cosmologyFunctionsData%expansionFactor       (                         &
            &                  self%cosmologyFunctionsData%timeAtDistanceComoving (                        &
            &                                                                      distanceDataMaximum     &
            &                                                                     )                        &
            &                                                                    )
       distanceModelMinimum  =+self%cosmologyFunctionsModel%distanceComoving     (                         &
            &                  self%cosmologyFunctionsModel%cosmicTime            (                        &
            &                                                                      expansionFactorMinimum  &
            &                                                                     )                        &
            &                                                                    )
       distanceModelMaximum  =+self%cosmologyFunctionsModel%distanceComoving     (                         &
            &                  self%cosmologyFunctionsModel%cosmicTime            (                        &
            &                                                                      expansionFactorMaximum  &
            &                                                                     )                        &
            &                                                                    )
       volumeData            =+volumeData                                                                  &
            &                 +self%surveyGeometry_        %solidAngle           (                         &
            &                                                                      field                   &
            &                                                                    )                         &
            &                 *(                                                                           &
            &                   +distanceDataMaximum **3                                                   &
            &                   -distanceDataMinimum **3                                                   &
            &                  )
       volumeModel           =+volumeModel                                                                 &
            &                 +self%surveyGeometry_        %solidAngle           (                         &
            &                                                                      field                   &
            &                                                                    )                         &
            &                 *(                                                                           &
            &                   +distanceModelMaximum**3                                                   &
            &                   -distanceModelMinimum**3                                                   &
            &                  )
    end do
    if (volumeData > 0.0d0) then
       correctionFactor=+volumeModel      &
            &           /volumeData
    else if (volumeModel <= 0.0d0) then
       ! Both model and data have zero volume - the galaxy is not within the survey limits, so the factor applied here is arbitrary.
       correctionFactor=1.0d0
    else
       correctionFactor=1.0d0
       call Error_Report('model volume is non-zero, but data volume is zero'//{introspection:location})
    end if
    ! Multiply by the correction factor.
    csmlgyVolumeOperate=+weightValue      &
         &              *correctionFactor
    return
  end function csmlgyVolumeOperate
