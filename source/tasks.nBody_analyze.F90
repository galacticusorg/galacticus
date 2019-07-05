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

  use NBody_Importers, only : nBodyImporterClass, nBodyImporter
  use NBody_Operators, only : nBodyOperatorClass, nBodyOperator

  !# <task name="taskNBodyAnalyze">
  !#  <description>A task which analyzes N-body simulation data.</description>
  !# </task>
  type, extends(taskClass) :: taskNBodyAnalyze
     !% Implementation of a task which analyzes N-body simulation data.
     private
     class  (nBodyImporterClass), pointer :: nBodyImporter_
     class  (nBodyOperatorClass), pointer :: nBodyOperator_
     type   (varying_string    )          :: nBodyFileName , nbodyFileNamePrevious
     logical                              :: usePrevious
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
    use Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (taskNBodyAnalyze  )                :: self
    type (inputParameters   ), intent(inout) :: parameters
    class(nBodyImporterClass), pointer       :: nBodyImporter_
    class(nBodyOperatorClass), pointer       :: nBodyOperator_
    type (varying_string    )                :: nBodyFileName , nbodyFileNamePrevious

    !# <inputParameter>
    !#   <name>nBodyFileName</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The name of the file containing the N-body data.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    if (parameters%isPresent('nbodyFileNamePrevious')) then
       !# <inputParameter>
       !#   <name>nbodyFileNamePrevious</name>
       !#   <cardinality>1</cardinality>
       !#   <description>The name of the file containing the N-body data.</description>
       !#   <source>parameters</source>
       !#   <type>string</type>
       !# </inputParameter>
    end if
    !# <objectBuilder class="nBodyImporter" name="nBodyImporter_" source="parameters"/>
    !# <objectBuilder class="nBodyOperator" name="nBodyOperator_" source="parameters"/>
    !# <conditionalCall>
    !#  <call>self=taskNBodyAnalyze(nBodyImporter_,nBodyOperator_,nBodyFileName{conditions})</call>
    !#  <argument name="nbodyFileNamePrevious" value="nbodyFileNamePrevious" parameterPresent="parameters"/>
    !# </conditionalCall>
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="nBodyImporter_"/>
    !# <objectDestructor name="nBodyOperator_"/>
    return
  end function nbodyAnalyzeConstructorParameters

  function nbodyAnalyzeConstructorInternal(nBodyImporter_,nBodyOperator_,nBodyFileName,nbodyFileNamePrevious) result(self)
    !% Constructor for the {\normalfont \ttfamily nbodyAnalyze} task class which takes a parameter set as input.
    implicit none
    type (taskNBodyAnalyze  )                          :: self
    class(nBodyImporterClass), intent(in   ), target   :: nBodyImporter_
    class(nBodyOperatorClass), intent(in   ), target   :: nBodyOperator_
    type (varying_string    ), intent(in   )           :: nBodyFileName
    type (varying_string    ), intent(in   ), optional :: nbodyFileNamePrevious
    !# <constructorAssign variables="*nBodyImporter_, *nBodyOperator_, nBodyFileName, nbodyFileNamePrevious"/>

    self%usePrevious=present(nbodyFileNamePrevious)
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
    use NBody_Simulation_Data, only : nBodyData
    use Galacticus_Error     , only : errorStatusSuccess
    use Galacticus_Display   , only : Galacticus_Display_Indent, Galacticus_Display_Unindent
    implicit none
    class  (taskNBodyAnalyze), intent(inout), target   :: self
    integer                  , intent(  out), optional :: status
    type   (nBodyData       )                          :: simulation

    call Galacticus_Display_Indent('Begin task: N-body analyze')
    if (.not.self%usePrevious) then
       simulation=self%nBodyImporter_%import(char(self%nbodyFileName)                                 )
    else
       simulation=self%nBodyImporter_%import(char(self%nbodyFileName),char(self%nbodyFileNamePrevious))
    end if
    ! Operate on the N-body data.
    call self%nBodyOperator_%operate(simulation)
    ! Close the analysis group.
    call simulation%analysis%close()
    ! Done.
    if (present(status)) status=errorStatusSuccess
    call Galacticus_Display_Unindent('Done task: N-body analyze' )
    return
  end subroutine nbodyAnalyzePerform

  logical function nbodyAnalyzeRequiresOutputFile(self)
    !% Specifies that this task does not requires the main output file.
    implicit none
    class(taskNBodyAnalyze), intent(inout) :: self    
    !GCC$ attributes unused :: self

    nbodyAnalyzeRequiresOutputFile=.false.
    return
  end function nbodyAnalyzeRequiresOutputFile
