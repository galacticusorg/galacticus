!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which implements a cosmological volume corrector analysis weight operator class.

  use Cosmology_Functions
  use Geometry_Surveys
  
  !# <outputAnalysisWeightOperator name="outputAnalysisWeightOperatorCsmlgyVolume" defaultThreadPrivate="yes">
  !#  <description>A cosmological volume corrector analysis weight operator class.</description>
  !# </outputAnalysisWeightOperator>
  type, extends(outputAnalysisWeightOperatorClass) :: outputAnalysisWeightOperatorCsmlgyVolume
     !% A cosmological volume corrector analysis weight operator class.
     private
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctionsModel, cosmologyFunctionsData
     class(surveyGeometryClass    ), pointer :: surveyGeometry_
   contains
     final     ::            csmlgyVolumeDestructor
     procedure :: operate => csmlgyVolumeOperate
  end type outputAnalysisWeightOperatorCsmlgyVolume

  interface outputAnalysisWeightOperatorCsmlgyVolume
     !% Constructors for the ``csmlgyVolume'' output analysis class.
     module procedure csmlgyVolumeConstructorParameters
     module procedure csmlgyVolumeConstructorInternal
  end interface outputAnalysisWeightOperatorCsmlgyVolume

contains

  function csmlgyVolumeConstructorParameters(parameters) result(self)
    !% Constructor for the ``csmlgyVolume'' output analysis weight operator class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type   (outputAnalysisWeightOperatorCsmlgyVolume)                :: self
    type   (inputParameters                         ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass                 ), pointer       :: cosmologyFunctionsModel, cosmologyFunctionsData
    class  (surveyGeometryClass                     ), pointer       :: surveyGeometry_
    type   (inputParameters                         )                :: dataAnalysisParameters
    
    ! Check and read parameters.
    dataAnalysisParameters=parameters%subParameters('dataAnalysis',requirePresent=.false.,requireValue=.false.)
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctionsModel" source="parameters"            />
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctionsData"  source="dataAnalysisParameters"/>
    !# <objectBuilder class="surveyGeometry"     name="surveyGeometry_"         source="dataAnalysisParameters"/>
    ! Construct the object.
    self=outputAnalysisWeightOperatorCsmlgyVolume(cosmologyFunctionsModel,cosmologyFunctionsData,surveyGeometry_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function csmlgyVolumeConstructorParameters

  function csmlgyVolumeConstructorInternal(cosmologyFunctionsModel,cosmologyFunctionsData,surveyGeometry_) result(self)
    !% Internal constructor for the ``csmlgyVolume'' output analysis weight operator class.
    implicit none
    type   (outputAnalysisWeightOperatorCsmlgyVolume)                        :: self
    class  (cosmologyFunctionsClass                 ), intent(in   ), target :: cosmologyFunctionsModel, cosmologyFunctionsData
    class  (surveyGeometryClass                     ), intent(in   ), target :: surveyGeometry_
    !# <constructorAssign variables="*cosmologyFunctionsModel, *cosmologyFunctionsData, *surveyGeometry_"/>

    return
  end function csmlgyVolumeConstructorInternal

  subroutine csmlgyVolumeDestructor(self)
    !% Destructor for the ``csmlgyVolume'' output analysis weight operator class.
    implicit none
    type(outputAnalysisWeightOperatorCsmlgyVolume), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctionsModel"/>
    !# <objectDestructor name="self%cosmologyFunctionsData" />
    !# <objectDestructor name="self%surveyGeometry_"        />
    return
  end subroutine csmlgyVolumeDestructor
  
  double precision function csmlgyVolumeOperate(self,weightValue,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
    !% Implement an csmlgyVolume output analysis weight operator.
    use, intrinsic :: ISO_C_Binding
    use            :: Output_Times
    use            :: Galacticus_Error
    use            :: Output_Analyses_Options
    implicit none
    class           (outputAnalysisWeightOperatorCsmlgyVolume), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    double precision                                          , intent(in   ) :: propertyValue         , propertyValueIntrinsic, &
         &                                                                       weightValue
    integer                                                   , intent(in   ) :: propertyType          , propertyQuantity
    integer         (c_size_t                                ), intent(in   ) :: outputIndex
    double precision                                                          :: distanceModelMinimum  , distanceDataMinimum   , &
         &                                                                       distanceModelMaximum  , distanceDataMaximum   , &
         &                                                                       expansionFactorMinimum, expansionFactorMaximum, &
         &                                                                       volumeData            , volumeModel           , &
         &                                                                       correctionFactor
    integer                                                                      field
    !GCC$ attributes unused :: outputIndex, propertyValue, propertyType, node

    ! Compute the correction factor - the assumption here is that the volume density was derived from a 1/Vₘₐₓ type approach. To
    ! correct for the distance in model and data cosmological models, we therefore first multiply the weight by Vₘₐₓ for the model
    ! cosmology (which gives the absolute number of galaxies per survey), then divide by the Vₘₐₓ for the data cosmology.    
    volumeData =0.0d0
    volumeModel=0.0d0
    do field=1,self%surveyGeometry_%fieldCount()
       select case (propertyQuantity)
       case (outputAnalysisPropertyQuantityMass      )
          select case (propertyType)
          case (outputAnalysisPropertyTypeLinear   )
             distanceDataMinimum   =+self%surveyGeometry_       %distanceMinimum       (                                           &
                  &                                                                      mass             =propertyValue         , &
                  &                                                                      field            =field                   &
                  &                                                                    )
             distanceDataMaximum   =+self%surveyGeometry_       %distanceMaximum       (                                           &
                  &                                                                      mass             =propertyValue         , &
                  &                                                                      field            =field                   &
                  &                                                                    )
          case (outputAnalysisPropertyTypeLog10    )
             distanceDataMinimum   =+self%surveyGeometry_       %distanceMinimum       (                                           &
                  &                                                                      mass             =propertyValueIntrinsic, &
                  &                                                                      field            =field                   &
                  &                                                                    )
             distanceDataMaximum   =+self%surveyGeometry_       %distanceMaximum       (                                           &
                  &                                                                      mass             =propertyValueIntrinsic, &
                  &                                                                      field            =field                   &
                  &                                                                    )
          case default
             call Galacticus_Error_Report('unsupported property type'//{introspection:location})
          end select
       case (outputAnalysisPropertyQuantityLuminosity)
          select case (propertyType)
          case (outputAnalysisPropertyTypeLinear   )
             distanceDataMinimum   =+self%surveyGeometry_       %distanceMinimum       (                                           &
                  &                                                                      luminosity       =propertyValue         , &
                  &                                                                      field            =field                   &
                  &                                                                    )
             distanceDataMaximum   =+self%surveyGeometry_       %distanceMaximum       (                                           &
                  &                                                                      luminosity       =propertyValue         , &
                  &                                                                      field            =field                   &
                  &                                                                    )
          case (outputAnalysisPropertyTypeLog10    )
             distanceDataMinimum   =+self%surveyGeometry_       %distanceMinimum       (                                           &
                  &                                                                      luminosity       =propertyValueIntrinsic, &
                  &                                                                      field            =field                   &
                  &                                                                    )
             distanceDataMaximum   =+self%surveyGeometry_       %distanceMaximum       (                                           &
                  &                                                                      luminosity       =propertyValueIntrinsic, &
                  &                                                                      field            =field                   &
                  &                                                                    )
          case (outputAnalysisPropertyTypeMagnitude    )
             distanceDataMinimum   =+self%surveyGeometry_       %distanceMinimum       (                                           &
                  &                                                                      magnitudeAbsolute=propertyValue         , &
                  &                                                                      field            =field                   &
                  &                                                                    )
             distanceDataMaximum   =+self%surveyGeometry_       %distanceMaximum       (                                           &
                  &                                                                      magnitudeAbsolute=propertyValue         , &
                  &                                                                      field            =field                   &
                  &                                                                    )
          case default
             call Galacticus_Error_Report('unsupported property type'//{introspection:location})
          end select
       case default
          call Galacticus_Error_Report('unsupported property class'//{introspection:location})
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
    correctionFactor   =+volumeModel      &
         &              /volumeData
    ! Multiply by the correction factor.
    csmlgyVolumeOperate=+weightValue      &
         &              *correctionFactor
    return
  end function csmlgyVolumeOperate
