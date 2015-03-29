!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  double precision, allocatable, dimension(:) :: luminosityFilterAbsoluteMagnitudeThresholdMinima
  double precision, allocatable, dimension(:) :: luminosityFilterAbsoluteMagnitudeThresholdMaxima

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
          if     (                                                                                                                                                  &
               &   (                                                                                                                                                &
               &     Input_Parameter_Is_Present    ('luminosityFilterAbsoluteMagnitudeThresholdMinima')                                                             &
               &    .and.                                                                                                                                           &
               &     Get_Input_Parameter_Array_Size('luminosityFilterAbsoluteMagnitudeThresholdMinima') /= unitStellarLuminosities%luminosityCount(unmapped=.true.) &
               &   )                                                                                                                                                &
               &  .or.                                                                                                                                              &
               &   (                                                                                                                                                &
               &     Input_Parameter_Is_Present    ('luminosityFilterAbsoluteMagnitudeThresholdMaxima')                                                             &
               &    .and.                                                                                                                                           &
               &     Get_Input_Parameter_Array_Size('luminosityFilterAbsoluteMagnitudeThresholdMaxima') /= unitStellarLuminosities%luminosityCount(unmapped=.true.) &
               &   )                                                                                                                                                &
               & )                                                                                                                                                  &
               & call  Galacticus_Error_Report(                                                                                                                     &
               &                                 'Galacticus_Merger_Tree_Output_Filter_Luminosity_Initialize'             ,                                         &
               &                                 'luminosityFilterAbsoluteMagnitudeThreshold(Minima|Maxima) input arrays '                                          &
               &                               //'must have same dimension as other luminosity arrays'                                                              &
               &                              )
          call Alloc_Array(luminosityFilterAbsoluteMagnitudeThresholdMinima,[unitStellarLuminosities%luminosityCount(unmapped=.true.)])
          call Alloc_Array(luminosityFilterAbsoluteMagnitudeThresholdMaxima,[unitStellarLuminosities%luminosityCount(unmapped=.true.)])
          ! Get the magnitude limits.
          !@ <inputParameter>
          !@   <name>luminosityFilterAbsoluteMagnitudeThresholdMinima</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The minimum absolute magnitudes (in the AB system) of a galaxy to pass the {\tt luminosity} output filter.
          !@   </description>
          !@   <defaultValue>$-\infty$</defaultValue>
          !@   <type>real</type>
          !@   <cardinality>0..*</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('luminosityFilterAbsoluteMagnitudeThresholdMinima',luminosityFilterAbsoluteMagnitudeThresholdMinima,defaultValue=spread(-HUGE(0.0d0),1,unitStellarLuminosities%luminosityCount(unmapped=.true.)))
          !@ <inputParameter>
          !@   <name>luminosityFilterAbsoluteMagnitudeThresholdMaxima</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The maximum absolute magnitudes (in the AB system) of a galaxy to pass the {\normalfont \ttfamily luminosity} output filter.
          !@   </description>
          !@   <defaultValue>$+\infty$</defaultValue>
          !@   <type>real</type>
          !@   <cardinality>0..*</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('luminosityFilterAbsoluteMagnitudeThresholdMaxima',luminosityFilterAbsoluteMagnitudeThresholdMaxima,defaultValue=spread(+HUGE(0.0d0),1,unitStellarLuminosities%luminosityCount(unmapped=.true.)))
          ! Map magnitude limits onto the expanded filter set.
          call Stellar_Luminosities_Parameter_Map(luminosityFilterAbsoluteMagnitudeThresholdMinima)
          call Stellar_Luminosities_Parameter_Map(luminosityFilterAbsoluteMagnitudeThresholdMaxima)
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
    !% Determines whether {\normalfont \ttfamily thisNode} has sufficient stellar mass to be output.
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
             if     (                                                                             &
                  &   abMagnitude > luminosityFilterAbsoluteMagnitudeThresholdMaxima(iLuminosity) &
                  &  .or.                                                                         &
                  &   abMagnitude < luminosityFilterAbsoluteMagnitudeThresholdMinima(iLuminosity) &
                  & ) then
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
