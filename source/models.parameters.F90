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

!!{RST
Contains a module that implements a class of parameter mapping functions.
!!}

module Model_Parameters
  !!{RST
  Implements a class of unary operators.
  !!}
  private
  public :: modelParameterList, modelParameterListLogPrior

  !![
  <functionClass docformat="rst">
   <name>modelParameter</name>
   <descriptiveName>Model Parameters</descriptiveName>
   <description>
   Class providing model parameters for Bayesian inference---the individual free parameters of a Galacticus model that are explored during parameter estimation (e.g.\ via MCMC). Each parameter has a name, a prior distribution (with ``logPrior``, ``priorSample``, ``priorInvert``, ``priorMinimum``, and ``priorMaximum`` methods), and a bijective mapping to an unconstrained real line for efficient sampling (via ``map``/``unmap``). Implementations include active parameters that vary during inference, and fixed parameters held at constant values.
   </description>
   <default>active</default>
   <method name="name">
     <description>
     Return the name of this parameter as it appears in the Galacticus parameter file and in output metadata, used to identify the parameter when applying posterior sampler updates.
     </description>
     <type>type(varying_string)</type>
     <pass>yes</pass>
   </method>
   <method name="logPrior">
     <description>
     Return the natural logarithm of the prior probability density evaluated at the physical parameter value ``x``, used in computing the log-posterior during Bayesian inference.
     </description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: x</argument>
   </method>
   <method name="priorSample">
     <description>
     Draw a random sample from the prior distribution for this parameter, returning a physical parameter value; used to initialize the posterior sampler or generate prior predictive samples.
     </description>
     <type>double precision</type>
     <pass>yes</pass>
   </method>
   <method name="priorInvert">
     <description>
     Invert the prior, returning the parameter value given the cumulative probability.
     </description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: f</argument>
   </method>
   <method name="priorMinimum">
     <description>
     Return the minimum non-zero value of the prior for this parameter.
     </description>
     <type>double precision</type>
     <pass>yes</pass>
   </method>
   <method name="priorMaximum">
     <description>
     Return the maximum non-zero value of the prior for this parameter.
     </description>
     <type>double precision</type>
     <pass>yes</pass>
   </method>
   <method name="randomPerturbation">
     <description>
     Return a random perturbation for this parameter.
     </description>
     <type>double precision</type>
     <pass>yes</pass>
   </method>
   <method name="map">
     <description>
     Apply the bijective mapping to transform the physical parameter value ``x`` onto the unconstrained real line used by the posterior sampler (e.g.\ logarithm for positive-definite parameters).
     </description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: x</argument>
   </method>
   <method name="unmap">
     <description>
     Apply the inverse bijective mapping to transform the sampler's unconstrained variable ``x`` back to the physical parameter value used to evaluate the model.
     </description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: x</argument>
   </method>
   <method name="mapJacobian">
     <description>
     Compute the Jacobian of the map at the given parameter value.
     </description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: x</argument>
   </method>
  </functionClass>
  !!]

  type :: modelParameterList
     !!{RST
     Class used to construct lists of model parameters.
     !!}
     class(modelParameterClass), public, pointer :: modelParameter_ => null()
   contains
     !![
     <methods>
       <method method="assignment(=)" description="Assign postprocessor list objects."/>
     </methods>
     !!]
     final     ::                  modelParameterListDestructor
     procedure ::                  modelParameterListAssign
     generic   :: assignment(=) => modelParameterListAssign
  end type modelParameterList

contains

  double precision function modelParameterListLogPrior(modelParameterList_,posteriorSampleState_)
    !!{RST
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

  subroutine modelParameterListDestructor(self)
    !!{RST
    Destructor for elements of model parameter lists.
    !!}
    implicit none
    type(modelParameterList), intent(inout) :: self

    !![
    <objectDestructor name="self%modelParameter_"/>
    !!]
    return
  end subroutine modelParameterListDestructor

  recursive subroutine modelParameterListAssign(self,from)
    !!{RST
    Perform assignment for the ``modelParameterList`` class.
    !!}
    implicit none
    class(modelParameterList), intent(  out) :: self
    class(modelParameterList), intent(in   ) :: from

    nullify(self%modelParameter_)
    if (associated(from%modelParameter_)) then
       self%modelParameter_ => from%modelParameter_
       !![
       <referenceCountIncrement owner="self" object="modelParameter_"/>
       !!]
    end if
    return
  end subroutine modelParameterListAssign

end module Model_Parameters
