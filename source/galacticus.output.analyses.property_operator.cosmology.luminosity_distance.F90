!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements a cosmological luminosity distance corrector analysis property operator class.

  use, intrinsic :: ISO_C_Binding
  use            :: Cosmology_Functions
  
  !# <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc" defaultThreadPrivate="yes">
  !#  <description>A cosmological luminosity distance corrector analysis property operator class.</description>
  !# </outputAnalysisPropertyOperator>
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc
     !% A cosmological luminosity distance corrector analysis property operator class.
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctionsModel, cosmologyFunctionsData
     integer         (c_size_t               )          :: outputIndexPrevious
     double precision                                   :: correctionFactor
   contains
     final     ::            csmlgyLuminosityDistanceDestructor
     procedure :: operate => csmlgyLuminosityDistanceOperate
  end type outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc

  interface outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc
     !% Constructors for the ``csmlgyLuminosityDistance'' output analysis class.
     module procedure csmlgyLuminosityDistanceConstructorParameters
     module procedure csmlgyLuminosityDistanceConstructorInternal
  end interface outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc

contains

  function csmlgyLuminosityDistanceConstructorParameters(parameters) result(self)
    !% Constructor for the ``csmlgyLuminosityDistance'' output analysis property operator class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type   (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc)                :: self
    type   (inputParameters                                ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass                        ), pointer       :: cosmologyFunctionsModel, cosmologyFunctionsData
    type   (inputParameters                                )                :: dataAnalysisParameters
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)
    dataAnalysisParameters=parameters%subParameters('dataAnalysis',requirePresent=.false.,requireValue=.false.)
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctionsModel" source="parameters"            />
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctionsData"  source="dataAnalysisParameters"/>
    ! Construct the object.
    self=outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc(cosmologyFunctionsModel,cosmologyFunctionsData)
    return
  end function csmlgyLuminosityDistanceConstructorParameters

  function csmlgyLuminosityDistanceConstructorInternal(cosmologyFunctionsModel,cosmologyFunctionsData) result(self)
    !% Internal constructor for the ``randomErrorPolynomial'' output analysis property operator class.
    implicit none
    type   (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc)                        :: self
    class  (cosmologyFunctionsClass                        ), intent(in   ), target :: cosmologyFunctionsModel, cosmologyFunctionsData
    !# <constructorAssign variables="*cosmologyFunctionsModel, *cosmologyFunctionsData"/>

    self%outputIndexPrevious=-1_c_size_t
    self%correctionFactor   =+0.0d0
    return
  end function csmlgyLuminosityDistanceConstructorInternal

  subroutine csmlgyLuminosityDistanceDestructor(self)
    !% Destructorfor the ``randomErrorPolynomial'' output analysis property operator class.
    implicit none
    type(outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctionsModel"/>
    !# <objectDestructor name="self%cosmologyFunctionsData" />
    return
  end subroutine csmlgyLuminosityDistanceDestructor
  
  double precision function csmlgyLuminosityDistanceOperate(self,propertyValue,propertyType,outputIndex)
    !% Implement an csmlgyLuminosityDistance output analysis property operator.
    use Galacticus_Output_Times
    use Galacticus_Error
    implicit none
    class           (outputAnalysisPropertyOperatorCsmlgyLmnstyDstnc), intent(inout)           :: self
    double precision                                                 , intent(in   )           :: propertyValue
    integer                                                          , intent(inout), optional :: propertyType
    integer         (c_size_t                                       ), intent(in   ), optional :: outputIndex
    double precision                                                                           :: distanceSmall=1.0d-6
    double precision                                                                           :: redshift            , timeData    , &
         &                                                                                        timeModel           , distanceData, &
         &                                                                                        distanceModel
    !GCC$ attributes unused :: propertyType

    ! Validate.
    if (.not.present(outputIndex)) call Galacticus_Error_Report('csmlgyLuminosityDistanceOperate','ouputIndex is required')
    ! Recompute correction factor if necessary.
    if (self%outputIndexPrevious /= outputIndex) then
       ! Store output index for which correction factor was computed.
       self%outputIndexPrevious=outputIndex
       ! Get current redshift.
       redshift        =Galacticus_Output_Redshift(outputIndex)
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
       ! Find luminosity distance in both data and model cosmological models.
       distanceData    =self%cosmologyFunctionsData %distanceLuminosity(timeData )
       distanceModel   =self%cosmologyFunctionsModel%distanceLuminosity(timeModel)
       ! Compute the correction factor - the assumption here is that the property was derived from a flux. Therefore, we first find
       ! the true flux observed in the model by dividing by the model luminosity distance squared, and then convert back to the
       ! property of interest under the assumptions of the data analysis by multiplying by the luminosity distance squared in the
       ! data analysis cosmological model.    
       ! Handle zero distance (i.e. present day outputs - and in fact we test for very small rather than zero distance to catch
       ! rounding errors) as a special case.
       if (distanceModel < distanceSmall) then
          if (distanceData > distanceSmall) call Galacticus_Error_Report('csmlgyLuminosityDistanceOperate','luminosity distance in model cosmology is zero, but non-zero for data cosmology - correction is undefined')
          self%correctionFactor= +1.0d0
       else
          self%correctionFactor=(               &
               &                 +distanceData  &
               &                 /distanceModel &
               &                )**2
       end if
    end if
    ! Multiply by the correction factor.
    csmlgyLuminosityDistanceOperate=+propertyValue         &
         &                          *self%correctionFactor
    return
  end function csmlgyLuminosityDistanceOperate
