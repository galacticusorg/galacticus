!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which defines the structure used for describing molecular abundances in \glc.

module Molecular_Abundances_Structure
  !% Defines the structure used for describing molecular abundances in \glc.
  use ISO_Varying_String
  private
  public :: molecularAbundancesStructure, Molecules_Names, Molecules_Index, Molecules_Property_Count
  ! <gfortran 4.6> The following method should only be callable as a type-bound procedure, but currently this does not seem to be
  ! recognized when called on a function result.
  public :: Molecules_Abundances_Reset
  
  type molecularAbundancesStructure
     !% The structure used for describing molecular abundances in \glc.
    ! private
     double precision, allocatable, dimension(:) :: molecularValue
   contains
     ! Pack/unpack methods.
     !@ <objectMethods>
     !@   <object>molecularAbundancesStructure</object>
     !@   <objectMethod>
     !@     <method>pack</method>
     !@     <description>Packs molecular abundance data into a 1-D array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>unpack(array)</method>
     !@     <description>Unpacks molecular abundance data from a 1-D array.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                 :: pack         => Molecular_Abundances_Pack
     procedure                 :: unpack       => Molecular_Abundances_Unpack
     ! Abundance methods.
     !@ <objectMethods>
     !@   <object>molecularAbundancesStructure</object>
     !@   <objectMethod>
     !@     <method>abundance</method>
     !@     <description>Returns the abundance of a molecule given its index.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>abundanceSet</method>
     !@     <description>Sets the abundance of a molecule given its index.</description>
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
     procedure                 :: abundance         => Molecules_Abundances
     procedure                 :: abundanceSet      => Molecules_Abundances_Set
     procedure                 :: reset             => Molecules_Abundances_Reset
     procedure                 :: numberToMass      => Molecules_Number_To_Mass
     procedure                 :: massToNumber      => Molecules_Mass_To_Number
     procedure                 :: multiply          => Molecules_Abundances_Multiply
     procedure                 :: divide            => Molecules_Abundances_Divide
  end type molecularAbundancesStructure

  ! Count of the number of elements being tracked.
  integer                                         :: moleculesCount=0
  integer                                         :: propertyCount

  ! Names of molecules to track.
  type(varying_string), allocatable, dimension(:) :: moleculesToTrack

  ! Indices of molecules as used in the Molecular_Structures module.
  integer,              allocatable, dimension(:) :: moleculesIndices

  ! Net charge and mass (in atomic units) of molecules.
  double precision,     allocatable, dimension(:) :: moleculesCharges,moleculesMasses

  ! Flag indicating if this module has been initialized.
  logical                                         :: molecularAbundancesInitialized=.false.

contains

  subroutine Molecular_Abundances_Initialize
    !% Initialize the {\tt molecularAbundanceStructure} object module. Determines which molecules are to be tracked.
    use Input_Parameters
    use Memory_Management
    use Molecular_Structures
    implicit none
    integer                  :: iMolecule
    type(molecularStructure) :: thisMolecule

    ! Check if this module has been initialized already.    
    !$omp critical (Molecular_Abundances_Module_Initalize)
    if (.not.molecularAbundancesInitialized) then

       ! Determine how many elements we are required to track.
       moleculesCount=Get_Input_Parameter_Array_Size('moleculesToTrack')
       ! Number of properties to track is the same as the number of molecules.
       propertyCount=moleculesCount
       ! If tracking molecules, read names of which ones to track.
       if (moleculesCount > 0) then
          allocate(moleculesToTrack(moleculesCount))
          call Alloc_Array(moleculesIndices,[moleculesCount])
          call Alloc_Array(moleculesCharges,[moleculesCount])
          call Alloc_Array(moleculesMasses ,[moleculesCount])
          !@ <inputParameter>
          !@   <name>moleculesToTrack</name>
          !@   <defaultValue></defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The names of the molecules to be tracked.
          !@   </description>
          !@ </inputParameter>
          call Get_Input_Parameter('moleculesToTrack',moleculesToTrack)
          ! Validate the input names by looking them up in the list of molecular names.
          do iMolecule=1,moleculesCount
             moleculesIndices(iMolecule)=Molecular_Database_Get_Index(char(moleculesToTrack(iMolecule)))
             call thisMolecule%retrieve(char(moleculesToTrack(iMolecule)))
             moleculesCharges(iMolecule)=dble(thisMolecule%charge())
             moleculesMasses (iMolecule)=     thisMolecule%mass  ()
          end do
       end if
       ! Flag that this module is now initialized.
       molecularAbundancesInitialized=.true.
    end if
    !$omp end critical (Molecular_Abundances_Module_Initalize)

    return
  end subroutine Molecular_Abundances_Initialize

  integer function Molecules_Property_Count()
    !% Return the number of properties required to track molecules. This is equal to the number of molecules tracked, {\tt
    !% moleculesCount}.
    implicit none

    ! Ensure module is initialized.
    call Molecular_Abundances_Initialize

    Molecules_Property_Count=propertyCount
    return
  end function Molecules_Property_Count

  function Molecules_Names(index)
    !% Return a name for the specified entry in the molecules structure.
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    type(varying_string)             :: Molecules_Names
    integer,              intent(in) :: index

    ! Ensure module is initialized.
    call Molecular_Abundances_Initialize

    if (index >= 1 .and. index <= moleculesCount) then
       Molecules_Names=trim(moleculesToTrack(index))
    else
       call Galacticus_Error_Report('Molecules_Names','index out of range')
    end if
    return
  end function Molecules_Names

  integer function Molecules_Index(moleculeName)
    !% Returns the index of a molecule in the molecular abundances structure given the {\tt moleculeName}.
    implicit none
    character(len=*), intent(in) :: moleculeName
    integer                      :: iMolecule

    Molecules_Index=-1 ! Indicates molecule not found.
    do iMolecule=1,moleculesCount
       if (moleculesToTrack(iMolecule) == trim(moleculeName)) then
          Molecules_Index=iMolecule
          return
       end if
    end do
    return
  end function Molecules_Index
  
  subroutine Molecules_Abundances_Multiply(theseAbundances,scaleFactor)
    !% Multiply molecular abundances in {\tt theseAbundances} by a scalar {\tt scaleFactor}.
    implicit none
#ifdef GCC45
    class(molecularAbundancesStructure), intent(inout) :: theseAbundances
#else
    type(molecularAbundancesStructure),  intent(inout) :: theseAbundances
#endif
    double precision,                    intent(in)    :: scaleFactor
    
    ! Ensure the target structure is allocated.
    call Molecular_Abundances_Allocate_Values(theseAbundances)

    ! Do the multiplication.
#ifdef GCC45
    select type (theseAbundances)
    type is (molecularAbundancesStructure)
#endif
       theseAbundances%molecularValue=theseAbundances%molecularValue*scaleFactor
#ifdef GCC45
    end select
#endif
    return
  end subroutine Molecules_Abundances_Multiply

  subroutine Molecules_Abundances_Divide(theseAbundances,scaleFactor)
    !% Divide molecular abundances in {\tt theseAbundances} by a scalar {\tt scaleFactor}.
    implicit none
#ifdef GCC45
    class(molecularAbundancesStructure), intent(inout) :: theseAbundances
#else
    type(molecularAbundancesStructure),  intent(inout) :: theseAbundances
#endif
    double precision,                    intent(in)    :: scaleFactor
    
    ! Ensure the target structure is allocated.
    call Molecular_Abundances_Allocate_Values(theseAbundances)

    ! Do the division.
#ifdef GCC45
    select type (theseAbundances)
    type is (molecularAbundancesStructure)
#endif
       theseAbundances%molecularValue=theseAbundances%molecularValue/scaleFactor
#ifdef GCC45
    end select
#endif
    return
  end subroutine Molecules_Abundances_Divide

  double precision function Molecules_Abundances(molecules,moleculeIndex)
    !% Returns the abundance of a molecule in the molecular abundances structure given the {\tt moleculeIndex}.
    implicit none
#ifdef GCC45
    class(molecularAbundancesStructure), intent(in) :: molecules
#else
    type(molecularAbundancesStructure),  intent(in) :: molecules
#endif
    integer,                             intent(in) :: moleculeIndex

    Molecules_Abundances=molecules%molecularValue(moleculeIndex)
    return
  end function Molecules_Abundances

  subroutine Molecules_Number_To_Mass(molecules,moleculesByMass)
    !% Multiply all molecular species by their mass in units of the atomic mass. This converts abundances by number into abundances by mass.
    implicit none
#ifdef GCC45
    class(molecularAbundancesStructure), intent(in)    :: molecules
#else
    type(molecularAbundancesStructure),  intent(in)    :: molecules 
#endif
    type(molecularAbundancesStructure),  intent(inout) :: moleculesByMass

    ! Ensure values array exists.
    call Molecular_Abundances_Allocate_Values(moleculesByMass)
    
    ! If the input molecular abundances structure is uninitialized, just return zero abundances.
    if (.not.allocated(molecules%molecularValue)) then
       call Molecules_Abundances_Reset(moleculesByMass)
       return
    end if

    ! Scale by the masses of the molecules.
    moleculesByMass%molecularValue=molecules%molecularValue*moleculesMasses
    return
  end subroutine Molecules_Number_To_Mass

  subroutine Molecules_Mass_To_Number(molecules,moleculesByNumber)
    !% Divide all molecular species by their mass in units of the atomic mass. This converts abundances by mass into abundances by number.
    implicit none
#ifdef GCC45
    class(molecularAbundancesStructure), intent(in)    :: molecules
#else
    type(molecularAbundancesStructure),  intent(in)    :: molecules 
#endif
    type(molecularAbundancesStructure),  intent(inout) :: moleculesByNumber

    ! Ensure values array exists.
    call Molecular_Abundances_Allocate_Values(moleculesByNumber)
    
    ! If the input molecular abundances structure is uninitialized, just return zero abundances.
    if (.not.allocated(molecules%molecularValue)) then
       call Molecules_Abundances_Reset(moleculesByNumber)
       return
    end if

    ! Scale by the masses of the molecules.
    moleculesByNumber%molecularValue=molecules%molecularValue/moleculesMasses
    return
  end subroutine Molecules_Mass_To_Number

  subroutine Molecules_Abundances_Set(molecules,moleculeIndex,abundance)
    !% Sets the abundance of a molecule in the molecular abundances structure given the {\tt moleculeIndex}.
    implicit none
#ifdef GCC45
    class(molecularAbundancesStructure), intent(inout) :: molecules
#else
    type(molecularAbundancesStructure),  intent(inout) :: molecules
#endif
    integer,                             intent(in)    :: moleculeIndex
    double precision,                    intent(in)    :: abundance
    
    ! Ensure values array exists.
    call Molecular_Abundances_Allocate_Values(molecules)
    
#ifdef GCC45
    select type (molecules)
    type is (molecularAbundancesStructure)
#endif
       molecules%molecularValue(moleculeIndex)=abundance
#ifdef GCC45
    end select
#endif
    return
  end subroutine Molecules_Abundances_Set

  subroutine Molecules_Abundances_Reset(molecules)
    !% Resets all molecular abundances to zero.
    implicit none
#ifdef GCC45
    class(molecularAbundancesStructure), intent(inout) :: molecules
#else
    type(molecularAbundancesStructure),  intent(inout) :: molecules
#endif

    ! Ensure values array exists.
    call Molecular_Abundances_Allocate_Values(molecules)
    
    ! Reset to zero.
#ifdef GCC45
    select type (molecules)
    type is (molecularAbundancesStructure)
#endif
       molecules%molecularValue=0.0d0
#ifdef GCC45
    end select
#endif
    return
  end subroutine Molecules_Abundances_Reset

  subroutine Molecular_Abundances_Allocate_Values(molecules)
    !% Ensure that the {\tt molecularValue} array in an {\tt moleculesStructure} is allocated.
    use Memory_Management
    implicit none
#ifdef GCC45
    class(molecularAbundancesStructure), intent(inout) :: molecules
#else
    type(molecularAbundancesStructure),  intent(inout) :: molecules
#endif

#ifdef GCC45
    select type (molecules)
    type is (molecularAbundancesStructure)
#endif
       if (.not.allocated(molecules%molecularValue)) then
          allocate(molecules%molecularValue(moleculesCount))
          call Memory_Usage_Record(sizeof(molecules%molecularValue))
       end if
#ifdef GCC45
    end select
#endif
    return
  end subroutine Molecular_Abundances_Allocate_Values

  subroutine Molecular_Abundances_Pack(molecules,molecularAbundancesArray)
    !% Pack abundances from an array into an abundances structure.
    implicit none
#ifdef GCC45
    class(molecularAbundancesStructure), intent(inout)            :: molecules
#else
    type(molecularAbundancesStructure),  intent(inout)            :: molecules
#endif
    double precision,                    intent(in), dimension(:) :: molecularAbundancesArray

    ! Ensure module is initialized.
    call Molecular_Abundances_Initialize

    ! Ensure values array exists.
    call Molecular_Abundances_Allocate_Values(molecules)
    ! Extract molecular values from array.
#ifdef GCC45
    select type (molecules)
    type is (molecularAbundancesStructure)
#endif    
       molecules%molecularValue=molecularAbundancesArray
#ifdef GCC45
    end select
#endif    
    return
  end subroutine Molecular_Abundances_Pack

  subroutine Molecular_Abundances_Unpack(molecules,molecularAbundancesArray)
    !% Pack abundances from an array into an abundances structure.
    implicit none
    double precision,                    intent(out), dimension(:) :: molecularAbundancesArray(:)
#ifdef GCC45
    class(molecularAbundancesStructure), intent(in)                :: molecules
#else
    type(molecularAbundancesStructure),  intent(in)                :: molecules
#endif

    ! Ensure module is initialized.
    call Molecular_Abundances_Initialize

    ! Place elemental values into arrays.
    if (allocated(molecules%molecularValue)) then
       molecularAbundancesArray=molecules%molecularValue
    else
       molecularAbundancesArray=0.0d0
    end if
    return
  end subroutine Molecular_Abundances_Unpack

end module Molecular_Abundances_Structure
