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
  Implementation of an active model parameter class.
  !!}

  use :: Math_Operators_Unary    , only : operatorUnaryClass
  use :: Statistics_Distributions, only : distributionFunction1DClass

  !![
  <modelParameter name="modelParameterActive">
   <description>An active model parameter class.</description>
  </modelParameter>
  !!]
  type, extends(modelParameterClass) :: modelParameterActive
     !!{
     Implementation of an active model parameter class.
     !!}
     private
     class  (distributionFunction1DClass), pointer :: prior  => null(), perturber => null()
     class  (operatorUnaryClass         ), pointer :: mapper => null()
     type   (varying_string             )          :: name_
     logical                                       :: slow
   contains
     !![
     <methods>
       <method method="isSlow" description="Return true if changes in this paper may lead to slow likelihood evaluation."/>
     </methods>
     !!]
     final     ::                       activeDestructor
     procedure :: name               => activeName
     procedure :: logPrior           => activeLogPrior
     procedure :: priorSample        => activePriorSample
     procedure :: priorInvert        => activePriorInvert
     procedure :: priorMinimum       => activePriorMinimum
     procedure :: priorMaximum       => activePriorMaximum
     procedure :: randomPerturbation => activeRandomPerturbation
     procedure :: map                => activeMap
     procedure :: unmap              => activeUnmap
     procedure :: mapJacobian        => activeMapJacobian
     procedure :: isSlow             => activeIsSlow
  end type modelParameterActive

  interface modelParameterActive
     !!{
     Constructors for the \refClass{modelParameterActive} 1D distribution function class.
     !!}
     module procedure activeConstructorParameters
     module procedure activeConstructorInternal
  end interface modelParameterActive

contains

  function activeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{modelParameterActive} model parameter class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (modelParameterActive       )                :: self
    type   (inputParameters            ), intent(inout) :: parameters
    class  (distributionFunction1DClass), pointer       :: prior     , perturber
    class  (operatorUnaryClass         ), pointer       :: mapper
    type   (varying_string             )                :: name
    logical                                             :: slow
    
    !![
    <inputParameter>
      <name>name</name>
      <description>The name of the parameter.</description>
      <defaultValue>var_str('')</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>slow</name>
      <description>If true, changes in the parameter are considered to result in slow likelihood evaluations.</description>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="distributionFunction1D" parameterName="distributionFunction1DPrior"     name="prior"     source="parameters"/>
    <objectBuilder class="distributionFunction1D" parameterName="distributionFunction1DPerturber" name="perturber" source="parameters"/>
    <objectBuilder class="operatorUnary"          parameterName="operatorUnaryMapper"             name="mapper"    source="parameters"/>
    !!]
    self=modelParameterActive(name,slow,prior,perturber,mapper)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="prior"    />
    <objectDestructor name="perturber"/>
    <objectDestructor name="mapper"   />
    !!]
    return
  end function activeConstructorParameters

  function activeConstructorInternal(name_,slow,prior,perturber,mapper) result(self)
    !!{
    Internal constructor for the \refClass{modelParameterActive} model parameter class.
    !!}
    implicit none
    type   (modelParameterActive       )                        :: self
    class  (distributionFunction1DClass), intent(in   ), target :: prior , perturber
    class  (operatorUnaryClass         ), intent(in   ), target :: mapper
    type   (varying_string             ), intent(in   )         :: name_
    logical                             , intent(in   )         :: slow
    !![
    <constructorAssign variables="name_, slow, *prior, *perturber, *mapper"/>
    !!]

    return
  end function activeConstructorInternal

  subroutine activeDestructor(self)
    !!{
    Destructor for \refClass{modelParameterActive} model parameter class.
    !!}
    implicit none
    type(modelParameterActive), intent(inout) :: self

    !![
    <objectDestructor name="self%prior"    />
    <objectDestructor name="self%perturber"/>
    <objectDestructor name="self%mapper"   />
    !!]
    return
  end subroutine activeDestructor

  function activeName(self)
    !!{
    Return the name of this parameter.
    !!}
    implicit none
    type (varying_string      )                :: activeName
    class(modelParameterActive), intent(inout) :: self

    activeName=self%name_
    return
  end function activeName

  double precision function activeLogPrior(self,x)
    !!{
    Return the log-prior on this parameter.
    !!}
    use :: Models_Likelihoods_Constants, only : logImpossible
    implicit none
    class           (modelParameterActive), intent(inout) :: self
    double precision                      , intent(in   ) :: x

    activeLogPrior=self%prior%density(x)
    if (activeLogPrior > 0.0d0) then
       activeLogPrior=log(activeLogPrior)
    else
       activeLogPrior=logImpossible
    end if
    return
  end function activeLogPrior

  double precision function activePriorSample(self)
    !!{
    Sample from the of this parameter.
    !!}
    implicit none
    class(modelParameterActive), intent(inout) :: self

    activePriorSample=self%prior%sample()
    return
  end function activePriorSample

  double precision function activePriorInvert(self,f)
    !!{
    Invert the prior, returning the parameter value given the cumulative probability.
    !!}
    implicit none
    class           (modelParameterActive), intent(inout) :: self
    double precision                      , intent(in   ) :: f

    activePriorInvert=self%prior%inverse(f)
    return
  end function activePriorInvert

  double precision function activePriorMinimum(self)
    !!{
    Return the minimum value for which the prior is non-zero.
    !!}
    implicit none
    class(modelParameterActive), intent(inout) :: self

    activePriorMinimum=self%prior%minimum()
    return
  end function activePriorMinimum

  double precision function activePriorMaximum(self)
    !!{
    Return the maximum value for which the prior is non-zero.
    !!}
    implicit none
    class(modelParameterActive), intent(inout) :: self

    activePriorMaximum=self%prior%maximum()
    return
  end function activePriorMaximum

  double precision function activeRandomPerturbation(self)
    !!{
    Return a random perturbation to this parameter.
    !!}
    implicit none
    class(modelParameterActive), intent(inout) :: self

    activeRandomPerturbation=self%perturber%sample()
    return
  end function activeRandomPerturbation

  double precision function activeMap(self,x)
    !!{
    Map this parameter.
    !!}
    implicit none
    class           (modelParameterActive), intent(inout) :: self
    double precision                      , intent(in   ) :: x

    activeMap=self%mapper%operate(x)
    return
  end function activeMap

  double precision function activeUnmap(self,x)
    !!{
    Unmap this parameter.
    !!}
    implicit none
    class           (modelParameterActive), intent(inout) :: self
    double precision                      , intent(in   ) :: x

    activeUnmap=self%mapper%unoperate(x)
    return
  end function activeUnmap

  double precision function activeMapJacobian(self,x)
    !!{
    Compute the Jacobian of the map for this parameter.
    !!}
    implicit none
    class           (modelParameterActive), intent(inout) :: self
    double precision                      , intent(in   ) :: x

    activeMapJacobian=self%mapper%jacobian(x)
    return
  end function activeMapJacobian

  logical function activeIsSlow(self)
    !!{
    Return true if changes in this parameter may result in slow likelihood evaluation.
    !!}
    implicit none
    class(modelParameterActive), intent(inout) :: self

    activeIsSlow=self%slow
    return
  end function activeIsSlow

