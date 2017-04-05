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

!% Contains a module which provides a collection of utilities useful for on-the-fly analyses.

module Output_Analysis_Utilities
  !% Provides a collection of utilities useful for on-the-fly analyses.
  implicit none
  private
  public :: Output_Analysis_Output_Weight_Survey_Volume

contains

  function Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,cosmologyFunctions_,massLimit,allowSingleEpoch) result (outputWeight)
    !% Compute output weights corresponding to the cosmological volumes associated with the given survey.
    use, intrinsic :: ISO_C_Binding
    use            :: Galacticus_Output_Times
    use            :: Geometry_Surveys
    use            :: Cosmology_Functions
    use            :: Galacticus_Error
    implicit none
    double precision                         , dimension(:) , allocatable :: outputWeight
    class           (surveyGeometryClass    ), intent(inout)              :: surveyGeometry_
    class           (cosmologyFunctionsClass), intent(inout)              :: cosmologyFunctions_
    double precision                         , intent(in   )              :: massLimit
    logical                                  , intent(in   ), optional    :: allowSingleEpoch
    double precision                         , parameter                  :: timeTolerance      =1.0d-6
    integer         (c_size_t               )                             :: iOutput
    integer                                                               :: iField
    double precision                                                      :: timeMinimum               , timeMaximum    , &
         &                                                                   distanceMinimum           , distanceMaximum
    !# <optionalArgument name="allowSingleEpoch" defaultsTo=".false." />
    
    allocate(outputWeight(Galacticus_Output_Time_Count()))
    outputWeight=0.0d0
    if (Galacticus_Output_Time_Count() == 1_c_size_t) then
       ! Handle cases where we have just a single epoch output.
       if (allowSingleEpoch_) then
          ! Iterate over all fields.
          do iField=1,surveyGeometry_%fieldCount()
             ! Test whether the output epoch lies within the range of comoving distances for this field.
             timeMinimum=cosmologyFunctions_%timeAtDistanceComoving(surveyGeometry_%distanceMaximum(massLimit,iField))
             timeMaximum=cosmologyFunctions_%timeAtDistanceComoving(surveyGeometry_%distanceMinimum(massLimit,iField))
             if     (                                                                         &
                  &   Galacticus_Output_Time(1_c_size_t) >= timeMinimum*(1.0d0-timeTolerance) &
                  &  .and.                                                                    &
                  &   Galacticus_Output_Time(1_c_size_t) <= timeMaximum*(1.0d0+timeTolerance) &
                  & ) outputWeight=1.0d0
          end do
          if (all(outputWeight == 0.0d0)) call Galacticus_Error_Report('Output_Analysis_Output_Weight_Survey_Volume','zero weights')
       else
          call Galacticus_Error_Report('Output_Analysis_Output_Weight_Survey_Volume','single epoch output not permitted')
       end if
    else
       do iOutput=1,Galacticus_Output_Time_Count()
          do iField=1,surveyGeometry_%fieldCount()
             if (iOutput == Galacticus_Output_Time_Count()) then
                timeMaximum=     Galacticus_Output_Time(iOutput)
             else
                timeMaximum=sqrt(Galacticus_Output_Time(iOutput)*Galacticus_Output_Time(iOutput+1))
             end if
             if (iOutput ==                              1) then
                timeMinimum=     Galacticus_Output_Time(iOutput)
             else
                timeMinimum=sqrt(Galacticus_Output_Time(iOutput)*Galacticus_Output_Time(iOutput-1))
             end if
             distanceMinimum=max(                                                          &
                  &              cosmologyFunctions_%distanceComoving(timeMaximum       ), &
                  &              surveyGeometry_    %distanceMinimum (massLimit  ,iField)  &
                  &             )
             distanceMaximum=min(                                                          &
                  &              cosmologyFunctions_%distanceComoving(timeMinimum       ), &
                  &              surveyGeometry_    %distanceMaximum (massLimit  ,iField)  &
                  &             )
             outputWeight                      (iOutput)  &
                  & =outputWeight              (iOutput)  &
                  & +surveyGeometry_%solidAngle(iField )  &
                  & /3.0d0                                &
                  & *max(                                 &
                  &      +0.0d0                         , &
                  &      +distanceMaximum**3              &
                  &      -distanceMinimum**3              &
                  &    )
          end do
       end do
       where(outputWeight < 0.0d0)
          outputWeight=0.0d0
       end where
       if (any(outputWeight > 0.0d0)) then
          outputWeight=+    outputWeight  &
               &       /sum(outputWeight)
       else
          call Galacticus_Error_Report('Output_Analysis_Output_Weight_Survey_Volume','zero weights')
       end if
    end if
    return
  end function Output_Analysis_Output_Weight_Survey_Volume
  
end module Output_Analysis_Utilities
