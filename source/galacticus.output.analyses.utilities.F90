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

!% Contains a module which provides a collection of utilities useful for on-the-fly analyses.

module Output_Analysis_Utilities
  !% Provides a collection of utilities useful for on-the-fly analyses.
  implicit none
  private
  public :: Output_Analysis_Output_Weight_Survey_Volume

contains

  function Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,cosmologyFunctions_,outputTimes_,massLimit,magnitudeAbsoluteLimit,luminosity,allowSingleEpoch) result (outputWeight)
    !% Compute output weights corresponding to the cosmological volumes associated with the given survey.
    use, intrinsic :: ISO_C_Binding
    use            :: Output_Times
    use            :: Geometry_Surveys
    use            :: Cosmology_Functions
    use            :: Galacticus_Error
    implicit none
    double precision                         , dimension(:) , allocatable :: outputWeight
    class           (surveyGeometryClass    ), intent(inout)              :: surveyGeometry_
    class           (cosmologyFunctionsClass), intent(inout)              :: cosmologyFunctions_
    class           (outputTimesClass       ), intent(inout)              :: outputTimes_
    double precision                         , intent(in   ), optional    :: massLimit                 , magnitudeAbsoluteLimit, &
         &                                                                   luminosity
    logical                                  , intent(in   ), optional    :: allowSingleEpoch
    double precision                         , parameter                  :: timeTolerance      =1.0d-6
    integer         (c_size_t               )                             :: iOutput
    integer                                                               :: iField
    double precision                                                      :: timeMinimum               , timeMaximum           , &
         &                                                                   distanceMinimum           , distanceMaximum       , &
         &                                                                   redshiftMinimum           , redshiftMaximum       , &
         &                                                                   time
    !# <optionalArgument name="allowSingleEpoch" defaultsTo=".false." />
    
    allocate(outputWeight(outputTimes_%count()))
    outputWeight=0.0d0
    if (outputTimes_%count() == 1_c_size_t) then
       ! Handle cases where we have just a single epoch output.
       if (allowSingleEpoch_) then
          ! Iterate over all fields.
          do iField=1,surveyGeometry_%fieldCount()
             ! Test whether the output epoch lies within the range of comoving distances for this field.
             time       =cosmologyFunctions_%cosmicTime            (cosmologyFunctions_%expansionFactorFromRedshift(outputTimes_%redshift(1_c_size_t)                                                         ))
             timeMinimum=cosmologyFunctions_%timeAtDistanceComoving(surveyGeometry_    %distanceMaximum            (mass=massLimit,magnitudeAbsolute=magnitudeAbsoluteLimit,luminosity=luminosity,field=iField))
             timeMaximum=cosmologyFunctions_%timeAtDistanceComoving(surveyGeometry_    %distanceMinimum            (mass=massLimit,magnitudeAbsolute=magnitudeAbsoluteLimit,luminosity=luminosity,field=iField))
             if     (                                           &
                  &   time >= timeMinimum*(1.0d0-timeTolerance) &
                  &  .and.                                      &
                  &   time <= timeMaximum*(1.0d0+timeTolerance) &
                  & ) outputWeight=1.0d0
          end do
          if (all(outputWeight == 0.0d0)) call Galacticus_Error_Report('zero weights'//{introspection:location})
       else
          call Galacticus_Error_Report('single epoch output not permitted'//{introspection:location})
       end if
    else
       do iOutput=1,outputTimes_%count()
          do iField=1,surveyGeometry_%fieldCount()
             if (iOutput == outputTimes_%count()) then
                redshiftMinimum=     outputTimes_%redshift(iOutput)
             else
                redshiftMinimum=sqrt(outputTimes_%redshift(iOutput)*outputTimes_%redshift(iOutput+1))
             end if
             if (iOutput ==                              1) then
                redshiftMaximum=     outputTimes_%redshift(iOutput)
             else
                redshiftMaximum=sqrt(outputTimes_%redshift(iOutput)*outputTimes_%redshift(iOutput-1))
             end if
             timeMinimum    =    cosmologyFunctions_%cosmicTime      (cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
             timeMaximum    =    cosmologyFunctions_%cosmicTime      (cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum))
             distanceMinimum=max(                                                                                                                                    &
                  &              cosmologyFunctions_%distanceComoving(     timeMaximum                                                                            ), &
                  &              surveyGeometry_    %distanceMinimum (mass=massLimit  ,magnitudeAbsolute=magnitudeAbsoluteLimit,luminosity=luminosity,field=iField)  &
                  &             )
             distanceMaximum=min(                                                                                                                                    &
                  &              cosmologyFunctions_%distanceComoving(     timeMinimum                                                                            ), &
                  &              surveyGeometry_    %distanceMaximum (mass=massLimit  ,magnitudeAbsolute=magnitudeAbsoluteLimit,luminosity=luminosity,field=iField)  &
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
          call Galacticus_Error_Report('zero weights'//{introspection:location})
       end if
    end if
    return
  end function Output_Analysis_Output_Weight_Survey_Volume
  
end module Output_Analysis_Utilities
