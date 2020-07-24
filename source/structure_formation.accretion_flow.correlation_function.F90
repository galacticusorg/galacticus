!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% An accretion flow class which models the accretion flow using the 2-halo correlation function.

  use :: Correlation_Functions_Two_Point, only : correlationFunctionTwoPointClass
  use :: Cosmology_Functions            , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Biases        , only : darkMatterHaloBiasClass

  !# <accretionFlows name="accretionFlowsCorrelationFunction">
  !#  <description>An accretion flow class which models the accretion flow using the 2-halo correlation function.</description>
  !# </accretionFlows>
  type, extends(accretionFlowsClass) :: accretionFlowsCorrelationFunction
     !% An accretion flow class which models the accretion flow using the 2-halo correlation function.
     private
     class(correlationFunctionTwoPointClass), pointer :: correlationFunctionTwoPoint_ => null()
     class(cosmologyFunctionsClass         ), pointer :: cosmologyFunctions_          => null()
     class(darkMatterHaloBiasClass         ), pointer :: darkMatterHaloBias_          => null()
   contains
     final     ::             correlationFunctionDestructor
     procedure :: density  => correlationFunctionDensity
     procedure :: velocity => correlationFunctionVelocity
  end type accretionFlowsCorrelationFunction

  interface accretionFlowsCorrelationFunction
     !% Constructors for the {\normalfont \ttfamily correlationFunction} accretion flows class.
     module procedure correlationFunctionConstructorParameters
     module procedure correlationFunctionConstructorInternal
  end interface accretionFlowsCorrelationFunction
  
contains

  function correlationFunctionConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily correlationFunction} accretion flow class that takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (accretionFlowsCorrelationFunction)                :: self
    type (inputParameters                  ), intent(inout) :: parameters
    class(correlationFunctionTwoPointClass ), pointer       :: correlationFunctionTwoPoint_
    class(cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class(darkMatterHaloBiasClass          ), pointer       :: darkMatterHaloBias_

    !# <objectBuilder class="correlationFunctionTwoPoint" name="correlationFunctionTwoPoint_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"          name="cosmologyFunctions_"          source="parameters"/>
    !# <objectBuilder class="darkMatterHaloBias"          name="darkMatterHaloBias_"          source="parameters"/>
    self=accretionFlowsCorrelationFunction(cosmologyFunctions_,correlationFunctionTwoPoint_,darkMatterHaloBias_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"         />
    !# <objectDestructor name="correlationFunctionTwoPoint_"/>
    !# <objectDestructor name="darkMatterHaloBias_"         />
    return
  end function correlationFunctionConstructorParameters

  function correlationFunctionConstructorInternal(cosmologyFunctions_,correlationFunctionTwoPoint_,darkMatterHaloBias_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily correlationFunction} accretion flows class.
    implicit none
    type (accretionFlowsCorrelationFunction)                        :: self
    class(correlationFunctionTwoPointClass ), intent(in   ), target :: correlationFunctionTwoPoint_
    class(cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
    class(darkMatterHaloBiasClass          ), intent(in   ), target :: darkMatterHaloBias_
    !# <constructorAssign variables="*cosmologyFunctions_, *correlationFunctionTwoPoint_, *darkMatterHaloBias_"/>

    return
  end function correlationFunctionConstructorInternal

  subroutine correlationFunctionDestructor(self)
    !% Destructor for the {\normalfont \ttfamily correlationFunction} accretion flows class.
    implicit none
    type(accretionFlowsCorrelationFunction), intent(inout) :: self
    
    !# <objectDestructor name="self%correlationFunctionTwoPoint_"/>
    !# <objectDestructor name="self%cosmologyFunctions_"         />
    !# <objectDestructor name="self%darkMatterHaloBias_"         />
    return
  end subroutine correlationFunctionDestructor
  
  double precision function correlationFunctionDensity(self,node,radius)
    !% Compute the density of the accretion flow at the given radius.
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (accretionFlowsCorrelationFunction), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    double precision                                   , intent(in   ) :: radius
    class           (nodeComponentBasic               ), pointer       :: basic
    double precision                                                   :: time
    
    basic                      =>  node %basic()
    time                       =   basic%time ()
    correlationFunctionDensity =  +(                                                                          &
         &                          +1.0d0                                                                    &
         &                          +self%darkMatterHaloBias_         %bias                (node,radius     ) &
         &                          *self%correlationFunctionTwoPoint_%correlation         (     radius,time) &
         &                         )                                                                          &
         &                        *  self%cosmologyFunctions_         %matterDensityEpochal(            time)
    return
  end function correlationFunctionDensity

  double precision function correlationFunctionVelocity(self,node,radius)
    !% Compute the velocity of the accretion flow at the given radius.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (accretionFlowsCorrelationFunction), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    double precision                                   , intent(in   ) :: radius
    !$GLC attributes unused :: self, node, radius
    
    correlationFunctionVelocity=0.0d0
    call Galacticus_Error_Report('velocity is currently unsupported'//{introspection:location})
    return
  end function correlationFunctionVelocity
