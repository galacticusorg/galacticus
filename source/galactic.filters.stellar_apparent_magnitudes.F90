!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a galactic low-pass (i.e. bright-pass) filter for stellar apparent magnitudes.

  use Cosmology_Functions
  
  !# <galacticFilter name="galacticFilterStellarApparentMagnitudes" defaultThreadPrivate="yes">
  !#  <description>
  !#  A galactic low-pass (i.e. bright-pass) filter for stellar apparent magnitudes. Galaxies with apparent magnitude in each
  !#  band, $i$, less than or equal to a fixed threshold, $m_{0,i}=${\normalfont \ttfamily [apparentMagnitudeThreshold]}.
  !#  </description>
  !# </galacticFilter>
  type, extends(galacticFilterClass) :: galacticFilterStellarApparentMagnitudes
     !% A galactic low-pass (i.e. bright pass) filter class for stellar apparent magnitudes.
     private
     class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctions_
     double precision                         , allocatable, dimension(:) :: apparentMagnitudeThreshold
   contains
     procedure :: passes => stellarApparentMagnitudesPasses
  end type galacticFilterStellarApparentMagnitudes

  interface galacticFilterStellarApparentMagnitudes
     !% Constructors for the ``stellarApparentMagnitudes'' galactic filter class.
     module procedure stellarApparentMagnitudesConstructorParameters
     module procedure stellarApparentMagnitudesConstructorInternal
  end interface galacticFilterStellarApparentMagnitudes

contains

  function stellarApparentMagnitudesConstructorParameters(parameters)
    !% Constructor for the ``stellarApparentMagnitudes'' galactic filter class which takes a parameter set as input.
    use Input_Parameters2
    use Stellar_Luminosities_Structure
    use Memory_Management
    use Galacticus_Error
    implicit none
    type            (galacticFilterStellarApparentMagnitudes)                              :: stellarApparentMagnitudesConstructorParameters
    type            (inputParameters                        ), intent(inout)               :: parameters
    double precision                                         , allocatable  , dimension(:) :: apparentMagnitudeThreshold
    class           (cosmologyFunctionsClass                ), pointer                     :: cosmologyFunctions_
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    if (parameters%count('apparentMagnitudeThreshold') /= unitStellarLuminosities%luminosityCount(unmapped=.true.)) &
         & call  Galacticus_Error_Report(                                                                                                  &
         &                                 'stellarApparentMagnitudesConstructorParameters'                                              , &
         &                                 '[apparentMagnitudeThreshold] input array must have same dimension as other luminosity arrays'  &
         &                              )
    call allocateArray(stellarApparentMagnitudesConstructorParameters%apparentMagnitudeThreshold,[unitStellarLuminosities%luminosityCount(unmapped=.true.)])
    !# <inputParameter>
    !#   <name>apparentMagnitudeThreshold</name>
    !#   <source>parameters</source>
    !#   <description>The parameter $m_0$ appearing in the stellar apparent magnitude threshold for the stellar apparent magnitude galactic filter class.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>          
    stellarApparentMagnitudesConstructorParameters=galacticFilterStellarApparentMagnitudes(apparentMagnitudeThreshold,cosmologyFunctions_)
   return
  end function stellarApparentMagnitudesConstructorParameters

  function stellarApparentMagnitudesConstructorInternal(apparentMagnitudeThreshold,cosmologyFunctions_)
    !% Internal constructor for the ``stellarApparentMagnitudes'' galactic filter class.
    use Stellar_Luminosities_Structure
    use Galacticus_Error
    implicit none
    type            (galacticFilterStellarApparentMagnitudes)                              :: stellarApparentMagnitudesConstructorInternal
    double precision                                         , intent(in   ), dimension(:) :: apparentMagnitudeThreshold
    class           (cosmologyFunctionsClass                ), intent(in   ), target       :: cosmologyFunctions_
    !# <constructorAssign variables="apparentMagnitudeThreshold, *cosmologyFunctions_"/>

    if (size(apparentMagnitudeThreshold) /= unitStellarLuminosities%luminosityCount(unmapped=.true.))                                    &
         & call  Galacticus_Error_Report(                                                                                                &
         &                               'stellarApparentMagnitudesConstructorInternal'                                                , &
         &                               '[apparentMagnitudeThreshold] input array must have same dimension as other luminosity arrays'  &
         &                              )    
    ! Map magnitude limits onto the expanded filter set.
    call Stellar_Luminosities_Parameter_Map(stellarApparentMagnitudesConstructorInternal%apparentMagnitudeThreshold)
    return
  end function stellarApparentMagnitudesConstructorInternal

  logical function stellarApparentMagnitudesPasses(self,node)
    !% Implement a stellar apparent magnitude low-pass galactic filter.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Stellar_Luminosities_Structure
    implicit none
    class           (galacticFilterStellarApparentMagnitudes), intent(inout) :: self
    type            (treeNode                               ), intent(inout) :: node
    class           (nodeComponentBasic                     ), pointer       :: basic
    double precision                                         , parameter     :: expansionFactorTolerance=1.0d-6
    double precision                                                         :: time                           , luminosity     , &
         &                                                                      abMagnitude                    , expansionFactor, &
         &                                                                      distanceModulus
    integer                                                                  :: iLuminosity
    
    ! Get the basic component.
    basic => node%basic()
    ! Get the time for this node.
    time  =  basic%time()
    ! Find the expansion factor.
    expansionFactor=self%cosmologyFunctions_%expansionFactor(time)
    ! If present day, do not apply this filter.
    if (expansionFactor > 1.0d0-expansionFactorTolerance) then
       stellarApparentMagnitudesPasses=.true.
    else
       ! Compute the distance modulus.
       distanceModulus=25.0d0+5.0d0*log10(self%cosmologyFunctions_%distanceLuminosity(time))+2.5d0*log10(expansionFactor)
       ! Loop over all luminosities.
       do iLuminosity=1,unitStellarLuminosities%luminosityCount()
          ! Only check those luminosities which are being output at this output time.
          if (unitStellarLuminosities%isOutput(iLuminosity,time)) then
             ! Get the total stellar luminosity of the galaxy.
             luminosity=Galactic_Structure_Enclosed_Mass(node,massType=massTypeStellar,weightBy=weightByLuminosity,weightIndex=iLuminosity)
             ! Test only if the luminosity is greater than zero.
             if (luminosity > 0.0d0) then
                ! Convert to apparent magnitude.
                abMagnitude=-2.5d0*log10(luminosity)+distanceModulus
                ! Filter out the galaxy if it is below the threshold.
                if (abMagnitude < self%apparentMagnitudeThreshold(iLuminosity)) then
                   stellarApparentMagnitudesPasses=.false.
                   return
                end if
             else
                stellarApparentMagnitudesPasses=.false.
                return
             end if
          end if
       end do
       stellarApparentMagnitudesPasses=.true.
    end if
    return
  end function stellarApparentMagnitudesPasses
