!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Implementation of a posterior sampling likelihood class which simply heats another likelihood to a specified temperature
  !!}

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodHeated">
   <description>
     The likelihood of the provided class is heated to a temperature $T=${\normalfont \ttfamily [temperature]}, such that
     \begin{equation}
     \log \mathcal{L} \rightarrow T^{-1} \log \mathcal{L}.
     \end{equation}
   </description>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodHeated
     !!{
     Implementation of a posterior sampling likelihood class which simply heats another likelihood to a specified temperature.
     !!}
     private
     class           (posteriorSampleLikelihoodClass), pointer :: posteriorSampleLikelihood_ => null()
     double precision                                          :: temperature
   contains
     final     ::                    heatedDestructor
     procedure :: evaluate        => heatedEvaluate
     procedure :: functionChanged => heatedFunctionChanged
  end type posteriorSampleLikelihoodHeated

  interface posteriorSampleLikelihoodHeated
     !!{
     Constructors for the {\normalfont \ttfamily heated} posterior sampling likelihood class.
     !!}
     module procedure heatedConstructorParameters
     module procedure heatedConstructorInternal
  end interface posteriorSampleLikelihoodHeated

contains

  function heatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily heated} posterior sampling likelihood class which builds the object
    from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleLikelihoodHeated)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (posteriorSampleLikelihoodClass ), pointer       :: posteriorSampleLikelihood_
    double precision                                                 :: temperature

    !![
    <inputParameter>
      <name>temperature</name>
      <description>The (dimensionless) temperature to which to heat the provided likelihood.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="posteriorSampleLikelihood" name="posteriorSampleLikelihood_" source="parameters"/>
    !!]
    self=posteriorSampleLikelihoodHeated(temperature,posteriorSampleLikelihood_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="posteriorSampleLikelihood_"/>
    !!]
    return
  end function heatedConstructorParameters

  function heatedConstructorInternal(temperature,posteriorSampleLikelihood_) result(self)
    !!{
    Constructor for ``heated'' posterior sampling likelihood class.
    !!}
    implicit none
    type            (posteriorSampleLikelihoodHeated)                        :: self
    class           (posteriorSampleLikelihoodClass ), intent(in   ), target :: posteriorSampleLikelihood_
    double precision                                 , intent(in   )         :: temperature
    !![
    <constructorAssign variables="temperature, *posteriorSampleLikelihood_"/>
    !!]
    return
  end function heatedConstructorInternal

  subroutine heatedDestructor(self)
    !!{
    Destructor for ``heated posterior sampling likelihood class.
    !!}
    implicit none
    type(posteriorSampleLikelihoodHeated), intent(inout) :: self

    !![
    <objectDestructor name="self%posteriorSampleLikelihood_"/>
    !!]
    return
  end subroutine heatedDestructor

  double precision function heatedEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{
    Return the log-likelihood for a ``heated'' likelihood function.
    !!}
    use :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    class           (posteriorSampleLikelihoodHeated), intent(inout)               :: self
    class           (posteriorSampleStateClass      ), intent(inout)               :: simulationState
    type            (modelParameterList             ), intent(in   ), dimension(:) :: modelParametersActive_, modelParametersInactive_
    class           (posteriorSampleConvergenceClass), intent(inout)               :: simulationConvergence
    double precision                                 , intent(in   )               :: temperature           , logLikelihoodCurrent    , &
         &                                                                            logPriorCurrent       , logPriorProposed
    real                                             , intent(inout)               :: timeEvaluate
    double precision                                 , intent(  out), optional     :: logLikelihoodVariance
    logical                                          , intent(inout), optional     :: forceAcceptance

    heatedEvaluate=self%posteriorSampleLikelihood_%evaluate(simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    if (heatedEvaluate > logImprobable) then
       heatedEvaluate             =heatedEvaluate       /self%temperature
       if (present(logLikelihoodVariance))                                   &
            &logLikelihoodVariance=logLikelihoodVariance/self%temperature**2
    end if
    return
  end function heatedEvaluate

  subroutine heatedFunctionChanged(self)
    !!{
    Respond to possible changes in the likelihood function.
    !!}
    implicit none
    class(posteriorSampleLikelihoodHeated), intent(inout) :: self

    call self%posteriorSampleLikelihood_%functionChanged()
    return
  end subroutine heatedFunctionChanged
