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

!% Contains a module which filters output on stellar mass.

module Galacticus_Merger_Tree_Output_Filter_Luminosities
  !% Filters output for lightcone geometry.
  implicit none
  private
  public :: Galacticus_Merger_Tree_Output_Filter_Luminosity,Galacticus_Merger_Tree_Output_Filter_Luminosity_Initialize

  ! Flags indicating if the module has been initialized and if this filter is active.
  logical                                     :: luminosityFilterInitialized                =.false.
  logical                                     :: luminosityFilterActive

  ! The absolute magnitude thresholds for this filter.
  double precision, allocatable, dimension(:) :: luminosityFilterAbsoluteMagnitudeThresholds

contains

  !# <mergerTreeOutputFilterInitialize>
  !#   <unitName>Galacticus_Merger_Tree_Output_Filter_Luminosity_Initialize</unitName>
  !# </mergerTreeOutputFilterInitialize>
  subroutine Galacticus_Merger_Tree_Output_Filter_Luminosity_Initialize(filterNames)
    !% Initializes the stellar mass filter module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    use Memory_Management
    use Stellar_Luminosities_Structure
    implicit none
    type(varying_string), dimension(:), intent(in   ) :: filterNames

    ! Initialize the filter if necessary.
    if (.not.luminosityFilterInitialized) then
       ! Determine if this filter has been selected.
       luminosityFilterActive=any(filterNames == "luminosity")
       ! If this filter is active, read the minimum stellar mass.
       if (luminosityFilterActive) then
          ! Check that the list of thresholds has the correct size.
          if (Get_Input_Parameter_Array_Size('luminosityFilterAbsoluteMagnitudeThresholds') /= unitStellarLuminosities%luminosityCount())                      &
               & call  Galacticus_Error_Report(                                                                                                                &
               &                               'Galacticus_Merger_Tree_Output_Filter_Luminosity_Initialize'                                                  , &
               &                               'luminosityFilterAbsoluteMagnitudeThresholds input arrays must have same dimension as other luminosity arrays'  &
               &                              )
          call Alloc_Array(luminosityFilterAbsoluteMagnitudeThresholds,[unitStellarLuminosities%luminosityCount()])

          ! Get the minimum stellar mass for output
          !@ <inputParameter>
          !@   <name>luminosityFilterAbsoluteMagnitudeThresholds</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The maximum absolute magnitudes (in the AB system) of a galaxy to pass the {\tt luminosity} output filter.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..*</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('luminosityFilterAbsoluteMagnitudeThresholds',luminosityFilterAbsoluteMagnitudeThresholds)
       end if
       ! Flag that this filter is now initialized.
       luminosityFilterInitialized=.true.
    end if
    return
  end subroutine Galacticus_Merger_Tree_Output_Filter_Luminosity_Initialize

  !# <mergerTreeOutputFilter>
  !#   <unitName>Galacticus_Merger_Tree_Output_Filter_Luminosity</unitName>
  !# </mergerTreeOutputFilter>
  subroutine Galacticus_Merger_Tree_Output_Filter_Luminosity(thisNode,doOutput)
    !% Determines whether {\tt thisNode} has sufficient stellar mass to be output.
    use Galacticus_Nodes
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    logical                             , intent(inout)          :: doOutput
    class           (nodeComponentBasic)               , pointer :: thisBasicComponent
    integer                                                      :: iLuminosity
    double precision                                             :: abMagnitude       , luminosity, time

    ! Return immediately if this filter is not active.
    if (.not.luminosityFilterActive) return

    ! Get the basic component.
    thisBasicComponent => thisNode%basic()

    ! Get the time for this node.
    time=thisBasicComponent%time()

    ! Loop over all luminosities.
    do iLuminosity=1,unitStellarLuminosities%luminosityCount()

       ! Only check those luminosities which are being output at this output time.
       if (unitStellarLuminosities%isOutput(iLuminosity,time)) then

          ! Get the total stellar luminosity of the galaxy.
          luminosity=Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeStellar,weightBy=weightByLuminosity,weightIndex=iLuminosity)

          ! Test only if the luminosity is greater than zero.
          if (luminosity > 0.0d0) then

             ! Convert to absolute magnitude.
             abMagnitude=-2.5d0*log10(luminosity)

             ! Filter out the galaxy if it is below the stellar mass threshold.
             if (abMagnitude > luminosityFilterAbsoluteMagnitudeThresholds(iLuminosity)) then
                doOutput=.false.
                return
             end if
          else
             doOutput=.false.
             return
          end if
       end if
    end do
    return
  end subroutine Galacticus_Merger_Tree_Output_Filter_Luminosity

end module Galacticus_Merger_Tree_Output_Filter_Luminosities
