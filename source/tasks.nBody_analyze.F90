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

  use :: NBody_Importers, only : nBodyImporter, nBodyImporterClass
  use :: NBody_Operators, only : nBodyOperator, nBodyOperatorClass

  !# <task name="taskNBodyAnalyze">
  !#  <description>A task which analyzes N-body simulation data.</description>
  !# </task>
  type, extends(taskClass) :: taskNBodyAnalyze
     !% Implementation of a task which analyzes N-body simulation data.
     private
     class  (nBodyImporterClass), pointer :: nBodyImporter_
     class  (nBodyOperatorClass), pointer :: nBodyOperator_
   contains
     final     ::                       nbodyAnalyzeDestructor
     procedure :: perform            => nbodyAnalyzePerform
     procedure :: requiresOutputFile => nbodyAnalyzeRequiresOutputFile
  end type taskNBodyAnalyze

  interface taskNBodyAnalyze
     !% Constructors for the {\normalfont \ttfamily nbodyAnalyze} task.
     module procedure nbodyAnalyzeConstructorParameters
     module procedure nbodyAnalyzeConstructorInternal
  end interface taskNBodyAnalyze

contains

  function nbodyAnalyzeConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily nbodyAnalyze} task class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (taskNBodyAnalyze  )                :: self
    type (inputParameters   ), intent(inout) :: parameters
    class(nBodyImporterClass), pointer       :: nBodyImporter_
    class(nBodyOperatorClass), pointer       :: nBodyOperator_

    !# <objectBuilder class="nBodyImporter" name="nBodyImporter_" source="parameters"/>
    !# <objectBuilder class="nBodyOperator" name="nBodyOperator_" source="parameters"/>
    self=taskNBodyAnalyze(nBodyImporter_,nBodyOperator_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="nBodyImporter_"/>
    !# <objectDestructor name="nBodyOperator_"/>
    return
  end function nbodyAnalyzeConstructorParameters

  function nbodyAnalyzeConstructorInternal(nBodyImporter_,nBodyOperator_) result(self)
    !% Constructor for the {\normalfont \ttfamily nbodyAnalyze} task class which takes a parameter set as input.
    implicit none
    type (taskNBodyAnalyze  )                          :: self
    class(nBodyImporterClass), intent(in   ), target   :: nBodyImporter_
    class(nBodyOperatorClass), intent(in   ), target   :: nBodyOperator_
    !# <constructorAssign variables="*nBodyImporter_, *nBodyOperator_"/>

    return
  end function nbodyAnalyzeConstructorInternal

  subroutine nbodyAnalyzeDestructor(self)
    !% Destructor for the {\normalfont \ttfamily nbodyAnalyze} task class.
    implicit none
    type(taskNBodyAnalyze), intent(inout) :: self

    !# <objectDestructor name="self%nBodyOperator_"/>
    !# <objectDestructor name="self%nBodyImporter_"/>
    return
  end subroutine nbodyAnalyzeDestructor

  subroutine nbodyAnalyzePerform(self,status)
    !% Compute and output the halo mass function.
    use :: Galacticus_Display   , only : Galacticus_Display_Indent, Galacticus_Display_Unindent
    use :: Galacticus_Error     , only : errorStatusSuccess
    use :: NBody_Simulation_Data, only : nBodyData
    implicit none
    class  (taskNBodyAnalyze), intent(inout), target       :: self
    integer                  , intent(  out), optional     :: status
    type   (nBodyData       ), allocatable  , dimension(:) :: simulations
    integer                                                :: i
    
    call Galacticus_Display_Indent('Begin task: N-body analyze')
    ! Import N-body data.
    call self%nBodyImporter_%import (simulations)
    ! Operate on the N-body data.
    call self%nBodyOperator_%operate(simulations)
    ! Close the analysis group.
    do i=1,size(simulations)
       if (simulations(i)%analysis%isOpen()) call simulations(i)%analysis%close()
    end do
    ! Done.
    if (present(status)) status=errorStatusSuccess
    call Galacticus_Display_Unindent('Done task: N-body analyze' )
    return
  end subroutine nbodyAnalyzePerform

  logical function nbodyAnalyzeRequiresOutputFile(self)
    !% Specifies that this task does not requires the main output file.
    implicit none
    class(taskNBodyAnalyze), intent(inout) :: self
    !$GLC attributes unused :: self

    nbodyAnalyzeRequiresOutputFile=.false.
    return
  end function nbodyAnalyzeRequiresOutputFile
