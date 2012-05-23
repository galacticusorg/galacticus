!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which defines the structure used for describing chemical abundances in \glc.

module Chemical_Abundances_Structure
  !% Defines the structure used for describing chemical abundances in \glc.
  use ISO_Varying_String
  implicit none
  private
  public :: chemicalAbundancesStructure, Chemicals_Names, Chemicals_Index, Chemicals_Property_Count
  
  type chemicalAbundancesStructure
     !% The structure used for describing chemical abundances in \glc.
     private
     double precision, allocatable, dimension(:) :: chemicalValue
   contains
     ! Pack/unpack methods.
     !@ <objectMethods>
     !@   <object>chemicalAbundancesStructure</object>
     !@   <objectMethod>
     !@     <method>pack</method>
     !@     <description>Packs chemical abundance data into a 1-D array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>unpack(array)</method>
     !@     <description>Unpacks chemical abundance data from a 1-D array.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                 :: pack         => Chemical_Abundances_Pack
     procedure                 :: unpack       => Chemical_Abundances_Unpack
     ! Abundance methods.
     !@ <objectMethods>
     !@   <object>chemicalAbundancesStructure</object>
     !@   <objectMethod>
     !@     <method>abundance</method>
     !@     <description>Returns the abundance of a chemical given its index.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>abundanceSet</method>
     !@     <description>Sets the abundance of a chemical given its index.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>reset</method>
     !@     <description>Resets abundances to zero.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>numberToMass</method>
     !@     <description>Converts from abundances by number to abundances by mass.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>massToNumber</method>
     !@     <description>Converts from abundances by mass to abundances by number.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>multiply</method>
     !@     <description>Multiplies abundances by given factor.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>divide</method>
     !@     <description>Divides abundances by given factor.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                 :: abundance         => Chemicals_Abundances
     procedure                 :: abundanceSet      => Chemicals_Abundances_Set
     procedure                 :: reset             => Chemicals_Abundances_Reset
     procedure                 :: numberToMass      => Chemicals_Number_To_Mass
     procedure                 :: massToNumber      => Chemicals_Mass_To_Number
     procedure                 :: multiply          => Chemicals_Abundances_Multiply
     procedure                 :: divide            => Chemicals_Abundances_Divide
  end type chemicalAbundancesStructure

  ! Count of the number of elements being tracked.
  integer                                         :: chemicalsCount=0
  integer                                         :: propertyCount

  ! Names of chemicals to track.
  type(varying_string), allocatable, dimension(:) :: chemicalsToTrack

  ! Indices of chemicals as used in the Chemical_Structures module.
  integer,              allocatable, dimension(:) :: chemicalsIndices

  ! Net charge and mass (in atomic units) of chemicals.
  double precision,     allocatable, dimension(:) :: chemicalsCharges,chemicalsMasses

  ! Flag indicating if this module has been initialized.
  logical                                         :: chemicalAbundancesInitialized=.false.

contains

  subroutine Chemical_Abundances_Initialize
    !% Initialize the {\tt chemicalAbundanceStructure} object module. Determines which chemicals are to be tracked.
    use Input_Parameters
    use Memory_Management
    use Chemical_Structures
    implicit none
    integer                 :: iChemical
    type(chemicalStructure) :: thisChemical

    ! Check if this module has been initialized already.    
    !$omp critical (Chemical_Abundances_Module_Initialize)
    if (.not.chemicalAbundancesInitialized) then

       ! Determine how many elements we are required to track.
       chemicalsCount=Get_Input_Parameter_Array_Size('chemicalsToTrack')
       ! Number of properties to track is the same as the number of chemicals.
       propertyCount=chemicalsCount
       ! If tracking chemicals, read names of which ones to track.
       if (chemicalsCount > 0) then
          allocate(chemicalsToTrack(chemicalsCount))
          call Alloc_Array(chemicalsIndices,[chemicalsCount])
          call Alloc_Array(chemicalsCharges,[chemicalsCount])
          call Alloc_Array(chemicalsMasses ,[chemicalsCount])
          !@ <inputParameter>
          !@   <name>chemicalsToTrack</name>
          !@   <defaultValue></defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The names of the chemicals to be tracked.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1..*</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('chemicalsToTrack',chemicalsToTrack)
          ! Validate the input names by looking them up in the list of chemical names.
          do iChemical=1,chemicalsCount
             chemicalsIndices(iChemical)=Chemical_Database_Get_Index(char(chemicalsToTrack(iChemical)))
             call thisChemical%retrieve(char(chemicalsToTrack(iChemical)))
             chemicalsCharges(iChemical)=dble(thisChemical%charge())
             chemicalsMasses (iChemical)=     thisChemical%mass  ()
          end do
       end if
       ! Flag that this module is now initialized.
       chemicalAbundancesInitialized=.true.
    end if
    !$omp end critical (Chemical_Abundances_Module_Initialize)

    return
  end subroutine Chemical_Abundances_Initialize

  integer function Chemicals_Property_Count()
    !% Return the number of properties required to track chemicals. This is equal to the number of chemicals tracked, {\tt
    !% chemicalsCount}.
    implicit none

    ! Ensure module is initialized.
    call Chemical_Abundances_Initialize

    Chemicals_Property_Count=propertyCount
    return
  end function Chemicals_Property_Count

  function Chemicals_Names(index)
    !% Return a name for the specified entry in the chemicals structure.
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    type(varying_string)             :: Chemicals_Names
    integer,              intent(in) :: index

    ! Ensure module is initialized.
    call Chemical_Abundances_Initialize

    if (index >= 1 .and. index <= chemicalsCount) then
       Chemicals_Names=trim(chemicalsToTrack(index))
    else
       call Galacticus_Error_Report('Chemicals_Names','index out of range')
    end if
    return
  end function Chemicals_Names

  integer function Chemicals_Index(chemicalName)
    !% Returns the index of a chemical in the chemical abundances structure given the {\tt chemicalName}.
    implicit none
    character(len=*), intent(in) :: chemicalName
    integer                      :: iChemical

    Chemicals_Index=-1 ! Indicates chemical not found.
    do iChemical=1,chemicalsCount
       if (chemicalsToTrack(iChemical) == trim(chemicalName)) then
          Chemicals_Index=iChemical
          return
       end if
    end do
    return
  end function Chemicals_Index
  
  subroutine Chemicals_Abundances_Multiply(theseAbundances,scaleFactor)
    !% Multiply chemical abundances in {\tt theseAbundances} by a scalar {\tt scaleFactor}.
    implicit none
    class(chemicalAbundancesStructure), intent(inout) :: theseAbundances
    double precision,                   intent(in)    :: scaleFactor
    
    ! Ensure the target structure is allocated.
    call Chemical_Abundances_Allocate_Values(theseAbundances)

    ! Do the multiplication.
    select type (theseAbundances)
    type is (chemicalAbundancesStructure)
       theseAbundances%chemicalValue=theseAbundances%chemicalValue*scaleFactor
    end select
    return
  end subroutine Chemicals_Abundances_Multiply

  subroutine Chemicals_Abundances_Divide(theseAbundances,scaleFactor)
    !% Divide chemical abundances in {\tt theseAbundances} by a scalar {\tt scaleFactor}.
    implicit none
    class(chemicalAbundancesStructure), intent(inout) :: theseAbundances
    double precision,                   intent(in)    :: scaleFactor
    
    ! Ensure the target structure is allocated.
    call Chemical_Abundances_Allocate_Values(theseAbundances)

    ! Do the division.
    select type (theseAbundances)
    type is (chemicalAbundancesStructure)
       theseAbundances%chemicalValue=theseAbundances%chemicalValue/scaleFactor
    end select
    return
  end subroutine Chemicals_Abundances_Divide

  double precision function Chemicals_Abundances(chemicals,moleculeIndex)
    !% Returns the abundance of a molecule in the chemical abundances structure given the {\tt moleculeIndex}.
    implicit none
    class(chemicalAbundancesStructure), intent(in) :: chemicals
    integer,                             intent(in) :: moleculeIndex

    Chemicals_Abundances=chemicals%chemicalValue(moleculeIndex)
    return
  end function Chemicals_Abundances

  subroutine Chemicals_Number_To_Mass(chemicals,chemicalsByMass)
    !% Multiply all chemical species by their mass in units of the atomic mass. This converts abundances by number into abundances by mass.
    implicit none
    class(chemicalAbundancesStructure), intent(in)    :: chemicals
    type(chemicalAbundancesStructure),  intent(inout) :: chemicalsByMass

    ! Ensure values array exists.
    call Chemical_Abundances_Allocate_Values(chemicalsByMass)
    
    ! If the input chemical abundances structure is uninitialized, just return zero abundances.
    if (.not.allocated(chemicals%chemicalValue)) then
       call Chemicals_Abundances_Reset(chemicalsByMass)
       return
    end if

    ! Scale by the masses of the chemicals.
    chemicalsByMass%chemicalValue=chemicals%chemicalValue*chemicalsMasses
    return
  end subroutine Chemicals_Number_To_Mass

  subroutine Chemicals_Mass_To_Number(chemicals,chemicalsByNumber)
    !% Divide all chemical species by their mass in units of the atomic mass. This converts abundances by mass into abundances by number.
    implicit none
    class(chemicalAbundancesStructure), intent(in)    :: chemicals
    type(chemicalAbundancesStructure),  intent(inout) :: chemicalsByNumber

    ! Ensure values array exists.
    call Chemical_Abundances_Allocate_Values(chemicalsByNumber)
    
    ! If the input chemical abundances structure is uninitialized, just return zero abundances.
    if (.not.allocated(chemicals%chemicalValue)) then
       call Chemicals_Abundances_Reset(chemicalsByNumber)
       return
    end if

    ! Scale by the masses of the chemicals.
    chemicalsByNumber%chemicalValue=chemicals%chemicalValue/chemicalsMasses
    return
  end subroutine Chemicals_Mass_To_Number

  subroutine Chemicals_Abundances_Set(chemicals,moleculeIndex,abundance)
    !% Sets the abundance of a molecule in the chemical abundances structure given the {\tt moleculeIndex}.
    implicit none
    class(chemicalAbundancesStructure), intent(inout) :: chemicals
    integer,                             intent(in)    :: moleculeIndex
    double precision,                    intent(in)    :: abundance
    
    ! Ensure values array exists.
    call Chemical_Abundances_Allocate_Values(chemicals)
    
    select type (chemicals)
    type is (chemicalAbundancesStructure)
       chemicals%chemicalValue(moleculeIndex)=abundance
    end select
    return
  end subroutine Chemicals_Abundances_Set

  subroutine Chemicals_Abundances_Reset(chemicals)
    !% Resets all chemical abundances to zero.
    implicit none
    class(chemicalAbundancesStructure), intent(inout) :: chemicals

    ! Ensure values array exists.
    call Chemical_Abundances_Allocate_Values(chemicals)

    ! Reset to zero.
    select type (chemicals)
    type is (chemicalAbundancesStructure)
       chemicals%chemicalValue=0.0d0
    end select
    return
  end subroutine Chemicals_Abundances_Reset

  subroutine Chemical_Abundances_Allocate_Values(chemicals)
    !% Ensure that the {\tt chemicalValue} array in an {\tt chemicalsStructure} is allocated.
    use Memory_Management
    implicit none
    class(chemicalAbundancesStructure), intent(inout) :: chemicals

    select type (chemicals)
    type is (chemicalAbundancesStructure)
       if (.not.allocated(chemicals%chemicalValue)) then
          allocate(chemicals%chemicalValue(chemicalsCount))
          call Memory_Usage_Record(sizeof(chemicals%chemicalValue))
       end if
    end select
    return
  end subroutine Chemical_Abundances_Allocate_Values

  subroutine Chemical_Abundances_Pack(chemicals,chemicalAbundancesArray)
    !% Pack abundances from an array into an abundances structure.
    implicit none
    class(chemicalAbundancesStructure), intent(inout)            :: chemicals
    double precision,                   intent(in), dimension(:) :: chemicalAbundancesArray

    ! Ensure module is initialized.
    call Chemical_Abundances_Initialize

    ! Ensure values array exists.
    call Chemical_Abundances_Allocate_Values(chemicals)
    ! Extract chemical values from array.
    select type (chemicals)
    type is (chemicalAbundancesStructure)
       chemicals%chemicalValue=chemicalAbundancesArray
    end select
    return
  end subroutine Chemical_Abundances_Pack

  subroutine Chemical_Abundances_Unpack(chemicals,chemicalAbundancesArray)
    !% Pack abundances from an array into an abundances structure.
    implicit none
    double precision,                   intent(out), dimension(:) :: chemicalAbundancesArray(:)
    class(chemicalAbundancesStructure), intent(in)                :: chemicals

    ! Ensure module is initialized.
    call Chemical_Abundances_Initialize

    ! Place elemental values into arrays.
    if (allocated(chemicals%chemicalValue)) then
       chemicalAbundancesArray=chemicals%chemicalValue
    else
       chemicalAbundancesArray=0.0d0
    end if
    return
  end subroutine Chemical_Abundances_Unpack

end module Chemical_Abundances_Structure
