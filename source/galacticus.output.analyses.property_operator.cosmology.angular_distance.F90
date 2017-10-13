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

!% Contains a module which implements a cosmological angular distance corrector analysis property operator class.

  use, intrinsic :: ISO_C_Binding
  use            :: Cosmology_Functions
  
  !# <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorCsmlgyAnglrDstnc" defaultThreadPrivate="yes">
  !#  <description>A cosmological angular distance corrector analysis property operator class.</description>
  !# </outputAnalysisPropertyOperator>
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorCsmlgyAnglrDstnc
     !% A cosmological angular distance corrector analysis property operator class.
     private
     class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctionsModel, cosmologyFunctionsData
     double precision                         , allocatable, dimension(:) :: correctionFactor
   contains
     final     ::            csmlgyAngularDistanceDestructor
     procedure :: operate => csmlgyAngularDistanceOperate
  end type outputAnalysisPropertyOperatorCsmlgyAnglrDstnc

  interface outputAnalysisPropertyOperatorCsmlgyAnglrDstnc
     !% Constructors for the ``csmlgyAngularDistance'' output analysis class.
     module procedure csmlgyAngularDistanceConstructorParameters
     module procedure csmlgyAngularDistanceConstructorInternal
  end interface outputAnalysisPropertyOperatorCsmlgyAnglrDstnc

contains

  function csmlgyAngularDistanceConstructorParameters(parameters) result(self)
    !% Constructor for the ``csmlgyAngularDistance'' output analysis property operator class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type (outputAnalysisPropertyOperatorCsmlgyAnglrDstnc)                :: self
    type (inputParameters                               ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                       ), pointer       :: cosmologyFunctionsModel, cosmologyFunctionsData
    type (inputParameters                               )                :: dataAnalysisParameters
    
    ! Check and read parameters.
    dataAnalysisParameters=parameters%subParameters('dataAnalysis',requirePresent=.false.,requireValue=.false.)
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctionsModel" source="parameters"            />
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctionsData"  source="dataAnalysisParameters"/>
    ! Construct the object.
    self=outputAnalysisPropertyOperatorCsmlgyAnglrDstnc(cosmologyFunctionsModel,cosmologyFunctionsData)
    !# <inputParametersValidate source="parameters"/>
    return
  end function csmlgyAngularDistanceConstructorParameters

  function csmlgyAngularDistanceConstructorInternal(cosmologyFunctionsModel,cosmologyFunctionsData) result(self)
    !% Internal constructor for the ``randomErrorPolynomial'' output analysis property operator class.
    use, intrinsic :: ISO_C_Binding
    use               Galacticus_Output_Times
    use               Memory_Management
    use               Galacticus_Error
    implicit none
    type            (outputAnalysisPropertyOperatorCsmlgyAnglrDstnc)                        :: self
    class           (cosmologyFunctionsClass                       ), intent(in   ), target :: cosmologyFunctionsModel       , cosmologyFunctionsData
    double precision                                                , parameter             :: distanceSmall          =1.0d-6
    integer         (c_size_t                                      )                        :: outputIndex
    double precision                                                                        :: redshift                      , timeData              , &
         &                                                                                     timeModel                     , distanceData          , &
         &                                                                                     distanceModel
    !# <constructorAssign variables="*cosmologyFunctionsModel, *cosmologyFunctionsData"/>

    call allocateArray(self%correctionFactor,[Galacticus_Output_Time_Count()])
    do outputIndex=1,Galacticus_Output_Time_Count()
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
       ! Find angular distance in both data and model cosmological models.
       distanceData    =self%cosmologyFunctionsData %distanceAngular(timeData )
       distanceModel   =self%cosmologyFunctionsModel%distanceAngular(timeModel)
       ! Compute the correction factor - the assumption here is that the property was derived from an angular size. Therefore, we
       ! first find the true physical size observed in the model by multiplying by the model angular distance, and then convert
       ! back to the property of interest under the assumptions of the data analysis by dividing by the angular distance in the
       ! data analysis cosmological model.  Handle zero distance (i.e. present day outputs - and in fact we test for very small
       ! rather than zero distance to catch rounding errors) as a special case.
       if (distanceModel < distanceSmall) then
          if (distanceData > distanceSmall) call Galacticus_Error_Report('angular distance in model cosmology is zero, but non-zero for data cosmology - correction is undefined'//{introspection:location})
          self%correctionFactor(outputIndex)= +1.0d0
       else
          self%correctionFactor(outputIndex)=+distanceModel &
               &                             /distanceData
       end if
    end do
    return
  end function csmlgyAngularDistanceConstructorInternal

  subroutine csmlgyAngularDistanceDestructor(self)
    !% Destructor for the ``csmlgyAnglrDstnc'' output analysis property operator class.
    use               Memory_Management
    implicit none
    type(outputAnalysisPropertyOperatorCsmlgyAnglrDstnc), intent(inout) :: self

    call deallocateArray(self%correctionFactor)
    !# <objectDestructor name="self%cosmologyFunctionsModel"/>
    !# <objectDestructor name="self%cosmologyFunctionsData" />
    return
  end subroutine csmlgyAngularDistanceDestructor
  
  double precision function csmlgyAngularDistanceOperate(self,propertyValue,node,propertyType,outputIndex)
    !% Implement an csmlgyAngularDistance output analysis property operator.
    use Galacticus_Error
    implicit none
    class           (outputAnalysisPropertyOperatorCsmlgyAnglrDstnc), intent(inout)           :: self
    double precision                                                , intent(in   )           :: propertyValue
    type            (treeNode                                      ), intent(inout), optional :: node
    integer                                                         , intent(inout), optional :: propertyType
    integer         (c_size_t                                      ), intent(in   ), optional :: outputIndex
    !GCC$ attributes unused :: propertyType, node

    ! Validate.
    if (.not.present(outputIndex)) call Galacticus_Error_Report('ouputIndex is required'//{introspection:location})
    ! Apply the correction.
    csmlgyAngularDistanceOperate=+propertyValue                      &
         &                       *self%correctionFactor(outputIndex)
    return
  end function csmlgyAngularDistanceOperate
