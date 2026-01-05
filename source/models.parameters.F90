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
Contains a module that implements a class of parameter mapping functions.
!!}

module Model_Parameters
  !!{
  Implements a class of unary operators.
  !!}
  private
  public :: modelParameterList, modelParameterListLogPrior

  !![
  <functionClass>
   <name>modelParameter</name>
   <descriptiveName>Model Parameters</descriptiveName>
   <description>Class providing model parameters.</description>
   <default>active</default>
   <method name="name">
     <type>type(varying_string)</type>
     <pass>yes</pass>
     <description>Return the name of this parameter.</description>
   </method>
   <method name="logPrior">
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: x</argument>
     <description>Return the log-prior for this parameter.</description>
   </method>
   <method name="priorSample">
     <type>double precision</type>
     <pass>yes</pass>
     <description>Sample from the parameter's prior.</description>
   </method>
   <method name="priorInvert">
     <type>double precision</type>
     <argument>double precision, intent(in   ) :: f</argument>
     <pass>yes</pass>
     <description>Invert the prior, returning the parameter value given the cumulative probability.</description>
   </method>
   <method name="priorMinimum">
     <type>double precision</type>
     <pass>yes</pass>
     <description>Return the minimum non-zero value of the prior for this parameter.</description>
   </method>
   <method name="priorMaximum">
     <type>double precision</type>
     <pass>yes</pass>
     <description>Return the maximum non-zero value of the prior for this parameter.</description>
   </method>
   <method name="randomPerturbation">
     <type>double precision</type>
     <pass>yes</pass>
     <description>Return a random perturbation for this parameter.</description>
   </method>
   <method name="map">
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: x</argument>
     <description>Map the parameter value.</description>
   </method>
   <method name="unmap">
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: x</argument>
     <description>Unmap the parameter value.</description>
   </method>
   <method name="mapJacobian">
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: x</argument>
     <description>Compute the Jacobian of the map at the given parameter value.</description>
   </method>
  </functionClass>
  !!]

  type :: modelParameterList
     !!{
     Class used to construct lists of model parameters.
     !!}
     class(modelParameterClass), public, pointer :: modelParameter_ => null()
  end type modelParameterList

contains

  double precision function modelParameterListLogPrior(modelParameterList_,posteriorSampleState_)
    !!{
    Compute the log-prior of a list of parameters.
    !!}
    use :: Models_Likelihoods_Constants, only : logImpossible
    use :: Posterior_Sampling_State    , only : posteriorSampleStateClass
    implicit none
    class           (modelParameterList       ), intent(in   ), dimension(:                        ) :: modelParameterList_
    class           (posteriorSampleStateClass), intent(inout)                                       :: posteriorSampleState_
    double precision                                          , dimension(size(modelParameterList_)) :: stateVector
    integer                                                                                          :: i
    double precision                                                                                 :: logPrior             , valueUnmapped

    stateVector               =posteriorSampleState_%get()
    modelParameterListLogPrior=0.0d0
    do i=1,size(modelParameterList_)
       valueUnmapped=modelParameterList_(i)%modelParameter_%unmap(stateVector(i))
       logPrior=modelParameterList_(i)%modelParameter_%logPrior   (valueUnmapped)
       if (logPrior <= logImpossible) then
          modelParameterListLogPrior=+logImpossible
          exit
       else
          modelParameterListLogPrior=+modelParameterListLogPrior                                                &
               &                     +logPrior                                                                  &
               &                     -log(                                                                      &
               &                          abs(                                                                  &
               &                              modelParameterList_(i)%modelParameter_%mapJacobian(valueUnmapped) &
               &                             )                                                                  &
               &                         )
       end if
    end do
    return
  end function modelParameterListLogPrior

end module Model_Parameters
