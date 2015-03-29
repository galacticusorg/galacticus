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

!% Contains a module which filters output on flux (apparent magnitude).

module Galacticus_Merger_Tree_Output_Filter_Fluxes
  !% Filters output for lightcone geometry.
  implicit none
  private
  public :: Galacticus_Merger_Tree_Output_Filter_Flux,Galacticus_Merger_Tree_Output_Filter_Flux_Initialize

  ! Flags indicating if the module has been initialized and if this filter is active.
  logical                                     :: fluxFilterInitialized                =.false.
  logical                                     :: fluxFilterActive

  ! The absolute magnitude thresholds for this filter.
  double precision, allocatable, dimension(:) :: fluxFilterAbsoluteMagnitudeThresholdMinima
  double precision, allocatable, dimension(:) :: fluxFilterAbsoluteMagnitudeThresholdMaxima

contains

  !# <mergerTreeOutputFilterInitialize>
  !#   <unitName>Galacticus_Merger_Tree_Output_Filter_Flux_Initialize</unitName>
  !# </mergerTreeOutputFilterInitialize>
  subroutine Galacticus_Merger_Tree_Output_Filter_Flux_Initialize(filterNames)
    !% Initializes the stellar mass filter module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    use Memory_Management
    use Stellar_Luminosities_Structure
    implicit none
    type(varying_string), dimension(:), intent(in   ) :: filterNames

    ! Initialize the filter if necessary.
    if (.not.fluxFilterInitialized) then
       ! Determine if this filter has been selected.
       fluxFilterActive=any(filterNames == "flux")
       ! If this filter is active, read the minimum stellar mass.
       if (fluxFilterActive) then
          ! Check that the list of thresholds has the correct size.
          if     (                                                                                                                                            &
               &   (                                                                                                                                          &
               &     Input_Parameter_Is_Present    ('fluxFilterAbsoluteMagnitudeThresholdMinima')                                                             &
               &    .and.                                                                                                                                     &
               &     Get_Input_Parameter_Array_Size('fluxFilterAbsoluteMagnitudeThresholdMinima') /= unitStellarLuminosities%luminosityCount(unmapped=.true.) &
               &   )                                                                                                                                          &
               &  .or.                                                                                                                                        &
               &   (                                                                                                                                          &
               &     Input_Parameter_Is_Present    ('fluxFilterAbsoluteMagnitudeThresholdMaxima')                                                             &
               &    .and.                                                                                                                                     &
               &     Get_Input_Parameter_Array_Size('fluxFilterAbsoluteMagnitudeThresholdMaxima') /= unitStellarLuminosities%luminosityCount(unmapped=.true.) &
               &   )                                                                                                                                          &
               & )                                                                                                                                            &
               & call  Galacticus_Error_Report(                                                                                                               &
               &                                 'Galacticus_Merger_Tree_Output_Filter_Flux_Initialize'             ,                                         &
               &                                 'fluxFilterAbsoluteMagnitudeThreshold(Minima|Maxima) input arrays '                                          &
               &                               //'must have same dimension as other luminosity arrays'                                                        &
               &                              )
          call Alloc_Array(fluxFilterAbsoluteMagnitudeThresholdMinima,[unitStellarLuminosities%luminosityCount(unmapped=.true.)])
          call Alloc_Array(fluxFilterAbsoluteMagnitudeThresholdMaxima,[unitStellarLuminosities%luminosityCount(unmapped=.true.)])
          ! Get the magnitude limits.
          !@ <inputParameter>
          !@   <name>fluxFilterAbsoluteMagnitudeThresholdMinima</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The minimum absolute magnitudes (in the AB system) of a galaxy to pass the {\tt flux} output filter.
          !@   </description>
          !@   <defaultValue>$-\infty$</defaultValue>
          !@   <type>real</type>
          !@   <cardinality>0..*</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('fluxFilterAbsoluteMagnitudeThresholdMinima',fluxFilterAbsoluteMagnitudeThresholdMinima,defaultValue=spread(-HUGE(0.0d0),1,unitStellarLuminosities%luminosityCount(unmapped=.true.)))
          !@ <inputParameter>
          !@   <name>fluxFilterAbsoluteMagnitudeThresholdMaxima</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The maximum absolute magnitudes (in the AB system) of a galaxy to pass the {\normalfont \ttfamily flux} output filter.
          !@   </description>
          !@   <defaultValue>$+\infty$</defaultValue>
          !@   <type>real</type>
          !@   <cardinality>0..*</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('fluxFilterAbsoluteMagnitudeThresholdMaxima',fluxFilterAbsoluteMagnitudeThresholdMaxima,defaultValue=spread(+HUGE(0.0d0),1,unitStellarLuminosities%luminosityCount(unmapped=.true.)))
          ! Map magnitude limits onto the expanded filter set.
          call Stellar_Luminosities_Parameter_Map(fluxFilterAbsoluteMagnitudeThresholdMinima)
          call Stellar_Luminosities_Parameter_Map(fluxFilterAbsoluteMagnitudeThresholdMaxima)
       end if
       ! Flag that this filter is now initialized.
       fluxFilterInitialized=.true.
    end if
    return
  end subroutine Galacticus_Merger_Tree_Output_Filter_Flux_Initialize

  !# <mergerTreeOutputFilter>
  !#   <unitName>Galacticus_Merger_Tree_Output_Filter_Flux</unitName>
  !# </mergerTreeOutputFilter>
  subroutine Galacticus_Merger_Tree_Output_Filter_Flux(node,doOutput)
    !% Determines whether {\normalfont \ttfamily node} has sufficient flux to be output.
    use Galacticus_Nodes
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Stellar_Luminosities_Structure
    use Cosmology_Functions
    implicit none
    type            (treeNode               ), intent(inout), pointer :: node
    logical                                  , intent(inout)          :: doOutput
    class           (nodeComponentBasic     )               , pointer :: basic
    class           (cosmologyFunctionsClass)               , pointer :: cosmologyFunctions_
    double precision                         , parameter              :: expansionFactorTolerance=1.0d-6
    integer                                                           :: iLuminosity
    double precision                                                  :: abMagnitude                    , luminosity     , &
         &                                                               time                           , distanceModulus, &
         &                                                               expansionFactor

    ! Return immediately if this filter is not active.
    if (.not.fluxFilterActive) return
    ! Get the basic component.
    basic => node%basic()
    ! Get the time for this node.
    time=basic%time()
    ! Find the expansion factor.
    expansionFactor=cosmologyFunctions_%expansionFactor(time)
    ! If present day, do not apply this filter.
    if (expansionFactor > 1.0d0-expansionFactorTolerance) return
    ! Compute the distance modulus.
    distanceModulus=25.0d0+5.0d0*log10(cosmologyFunctions_%distanceLuminosity(time))+2.5d0*log10(expansionFactor)
    ! Loop over all luminosities.
    do iLuminosity=1,unitStellarLuminosities%luminosityCount()
       ! Only check those luminosities which are being output at this output time.
       if (unitStellarLuminosities%isOutput(iLuminosity,time)) then
          ! Get the total stellar luminosity of the galaxy.
          luminosity=Galactic_Structure_Enclosed_Mass(node,massType=massTypeStellar,weightBy=weightByLuminosity,weightIndex=iLuminosity)
          ! Test only if the luminosity is greater than zero.
          if (luminosity > 0.0d0) then
             ! Convert to absolute magnitude.
             abMagnitude=-2.5d0*log10(luminosity)+distanceModulus
             ! Filter out the galaxy if it is below the stellar mass threshold.
             if     (                                                                       &
                  &   abMagnitude > fluxFilterAbsoluteMagnitudeThresholdMaxima(iLuminosity) &
                  &  .or.                                                                   &
                  &   abMagnitude < fluxFilterAbsoluteMagnitudeThresholdMinima(iLuminosity) &
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
  end subroutine Galacticus_Merger_Tree_Output_Filter_Flux

end module Galacticus_Merger_Tree_Output_Filter_Fluxes
