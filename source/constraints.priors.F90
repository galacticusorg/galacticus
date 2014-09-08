!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements priors for use when constraining \glc.

module Constraints_Priors
  !% Implements priors for use when constraining \glc.
  use Statistics_Distributions
  private
  public :: priorsEvaluateLog, priorsSample

  type, public :: prior
     !% A class used to describe priors for use when constraining \glc.
     integer                        :: variableIndex
     logical                        :: autoDistribution
     class  (distribution), pointer :: priorDistribution     
   contains
     !# <workaround type="gfortran" PR="58471 58470" url="http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58471 http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58470">
     !# final     ::               priorDestructor
     !# </workaround>
     !@ <objectMethods>
     !@   <object>prior</object>
     !@   <objectMethod>
     !@     <method>sample</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone state\arginout</arguments>
     !@     <description>Sample from the prior, populating the given {\tt state}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>logDensity</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless class(state)\textgreater} simulationState\argin</arguments>
     !@     <description>Return the log of the probability density of the prior at the given {\tt simulationState}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>invert</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero p\argin</arguments>
     !@     <description>Return the value of the random variable corresponding to the given cumulative probability.</description>
     !@   </objectMethod>
     !@ </objectMethods>
    procedure :: sample     => priorSample
    procedure :: logDensity => priorLogDensity
    procedure :: invert     => priorInvert
  end type prior
  
  interface prior
     !% Interface to constructors for the {\tt prior} class.
     module procedure priorConstructor
     module procedure priorConstructorXML
  end interface prior

contains

  function priorConstructor(priorDistribution,variableIndex)
    !% Construct a prior object.
    implicit none
    type   (prior       )                        :: priorConstructor
    class  (distribution), intent(in   ), target :: priorDistribution
    integer              , intent(in   )         :: variableIndex
    
    priorConstructor%variableIndex     =  variableIndex
    priorConstructor%priorDistribution => priorDistribution
    priorConstructor%autoDistribution  =  .false.
    return
  end function priorConstructor

  function priorConstructorXML(definition,variableIndex)
    !% Construct a prior object.
    use FoX_DOM
    use IO_XML
    implicit none
    type   (prior)                         :: priorConstructorXML
    type   (node ), pointer, intent(in   ) :: definition
    integer                , intent(in   ) :: variableIndex
    type   (node ), pointer                :: distributionDefinition

    ! Assign the variable index to which this prior applies.
    priorConstructorXML%variableIndex     =  variableIndex
    ! Construct the distribution for this prior.
    distributionDefinition                => XML_Get_First_Element_By_Tag_Name(definition,"distribution")
    priorConstructorXML%priorDistribution => distributionNew(distributionDefinition)
    priorConstructorXML%autoDistribution  =  .true.
    return
  end function priorConstructorXML
  
  elemental subroutine priorDestructor(self)
    !% Destroy a prior object.
    implicit none
    type(prior), intent(inout) :: self
 
    if (self%autoDistribution) deallocate(self%priorDistribution)
    return
  end subroutine priorDestructor

  subroutine priorSample(self,state)
    !% Sample from the prior, populating the given {\tt state}.
    implicit none
    class           (prior), intent(inout)               :: self
    double precision       , intent(inout), dimension(:) :: state

    select type (priorDistribution => self%priorDistribution)
    class is (distribution1D)
       state(self%variableIndex)=priorDistribution%sample()
    end select
    return
  end subroutine priorSample

  subroutine priorsSample(priors,priorState)
    !% Sample from all priors populating the given {\tt priorState}.
    use Constraints_State
    implicit none
    type            (prior), intent(inout), dimension(:) :: priors
    class           (state), intent(inout)               :: priorState
    double precision       , allocatable  , dimension(:) :: stateVector
    integer                                              :: i

    stateVector=priorState%get()
    do i=1,size(priors)
       call priors(i)%sample(stateVector)
    end do
    call priorState%update(stateVector,.false.,.false.)
    deallocate(stateVector)
    return
  end subroutine priorsSample

  double precision function priorLogDensity(self,simulationState)
    !% Return the logarithm of the density of the prior at the given {\tt simulationState}.
    use Constraints_Constants
    use Constraints_State
    implicit none
    class           (prior), intent(in   )             :: self
    class           (state), intent(in   )             :: simulationState
    double precision       , dimension(:), allocatable :: stateVector
    double precision                                   :: x              , density

    select type (priorDistribution => self%priorDistribution)
    class is (distribution1D)
       stateVector=simulationState%get()
       x          =stateVector(self%variableIndex)
       density    =priorDistribution%density(x)
       if (density <= 0.0d0) then
          priorLogDensity=logImpossible
       else
          priorLogDensity=log(density)
       end  if
    end select
    return
  end function priorLogDensity

  double precision function priorsEvaluateLog(priors,simulationState)
    !% Evaluate the logarithm of the prior for the given {\tt state}.
    use Constraints_Constants
    use Constraints_State
    implicit none
    class           (state), intent(in   )               :: simulationState
    type            (prior), intent(in   ), dimension(:) :: priors
    double precision                                     :: thisLogPrior
    integer                                              :: i
    
    priorsEvaluateLog=0.0d0
    do i=1,size(priors)
       thisLogPrior=priors(i)%logDensity(simulationState)
       if (thisLogPrior <= logImpossible) then
          priorsEvaluateLog=logImpossible
          exit
       end if
       priorsEvaluateLog=priorsEvaluateLog+thisLogPrior
    end do
    return
  end function priorsEvaluateLog

  double precision function priorInvert(self,p)
    !% Return the value of the prior given the cumulative probability, {\tt p}.
    use Galacticus_Error
    implicit none
    class           (prior), intent(in   )             :: self
    double precision       , intent(in   ) :: p

    select type (priorDistribution => self%priorDistribution)
    class is (distribution1D)
       priorInvert=priorDistribution%inverse(p)
    class default
       call Galacticus_Error_Report('priorInvert','unable to invert prior')
    end select
    return
  end function priorInvert

end module Constraints_Priors
