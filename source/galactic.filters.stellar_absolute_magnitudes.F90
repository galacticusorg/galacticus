!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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

!!{
Contains a module which implements a galactic low-pass (i.e. bright-pass) filter for stellar absolute magnitudes.
!!}
  
  use :: Galactic_Structure, only : galacticStructureClass

  !![
  <galacticFilter name="galacticFilterStellarAbsoluteMagnitudes">
   <description>
   A galactic low-pass (i.e. bright-pass) filter for stellar absolute magnitudes. Galaxies with absolute magnitude in each
   band, $i$, less than or equal to a fixed threshold, $M_{0,i}=${\normalfont \ttfamily [absoluteMagnitudeThreshold]}.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterStellarAbsoluteMagnitudes
     !!{
     A galactic low-pass (i.e. bright pass) filter class for stellar absolute magnitudes.
     !!}
     private
     class           (galacticStructureClass), pointer                   :: galacticStructure_         => null()
     double precision                        , allocatable, dimension(:) :: absoluteMagnitudeThreshold
   contains
     final     ::           stellarAbsoluteMagnitudesDestructor
     procedure :: passes => stellarAbsoluteMagnitudesPasses
  end type galacticFilterStellarAbsoluteMagnitudes

  interface galacticFilterStellarAbsoluteMagnitudes
     !!{
     Constructors for the ``stellarAbsoluteMagnitudes'' galactic filter class.
     !!}
     module procedure stellarAbsoluteMagnitudesConstructorParameters
     module procedure stellarAbsoluteMagnitudesConstructorInternal
  end interface galacticFilterStellarAbsoluteMagnitudes

contains

  function stellarAbsoluteMagnitudesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``stellarAbsoluteMagnitudes'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Error                         , only : Error_Report
    use :: Input_Parameters              , only : inputParameter         , inputParameters
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    type            (galacticFilterStellarAbsoluteMagnitudes)                              :: self
    type            (inputParameters                        ), intent(inout)               :: parameters
    double precision                                         , allocatable  , dimension(:) :: absoluteMagnitudeThreshold
    class           (galacticStructureClass                 ), pointer                     :: galacticStructure_

    ! Check and read parameters.
    if (parameters%count('absoluteMagnitudeThreshold') /= unitStellarLuminosities%luminosityCount(unmapped=.true.))            &
         & call  Error_Report(                                                                                                 &
         &                    '[absoluteMagnitudeThreshold] input array must have same dimension as other luminosity arrays'// &
         &                    {introspection:location}                                                                         &
         &                   )
    allocate(absoluteMagnitudeThreshold(unitStellarLuminosities%luminosityCount(unmapped=.true.)))
    !![
    <inputParameter>
      <name>absoluteMagnitudeThreshold</name>
      <source>parameters</source>
      <description>The parameter $M_0$ appearing in the stellar absolute magnitude threshold for the stellar absolute magnitude galactic filter class.</description>
    </inputParameter>
    <objectBuilder class="galacticStructure" name="galacticStructure_" source="parameters"/>
    !!]
    self=galacticFilterStellarAbsoluteMagnitudes(absoluteMagnitudeThreshold,galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticStructure_"/>
    !!]
    return
  end function stellarAbsoluteMagnitudesConstructorParameters

  function stellarAbsoluteMagnitudesConstructorInternal(absoluteMagnitudeThreshold,galacticStructure_) result(self)
    !!{
    Internal constructor for the ``stellarAbsoluteMagnitudes'' galactic filter class.
    !!}
    use :: Error                         , only : Error_Report
    use :: Stellar_Luminosities_Structure, only : Stellar_Luminosities_Parameter_Map, unitStellarLuminosities
    implicit none
    type            (galacticFilterStellarAbsoluteMagnitudes)                              :: self
    class           (galacticStructureClass                 ), intent(in   ), target       :: galacticStructure_
    double precision                                         , intent(in   ), dimension(:) :: absoluteMagnitudeThreshold
    !![
    <constructorAssign variables="absoluteMagnitudeThreshold, *galacticStructure_"/>
    !!]

    if (size(absoluteMagnitudeThreshold) /= unitStellarLuminosities%luminosityCount(unmapped=.true.))                          &
         & call  Error_Report(                                                                                                 &
         &                    '[absoluteMagnitudeThreshold] input array must have same dimension as other luminosity arrays'// &
         &                    {introspection:location}                                                                         &
         &                   )
    ! Map magnitude limits onto the expanded filter set.
    call Stellar_Luminosities_Parameter_Map(self%absoluteMagnitudeThreshold)
    return
  end function stellarAbsoluteMagnitudesConstructorInternal

  subroutine stellarAbsoluteMagnitudesDestructor(self)
    !!{
    Destructor for the ``stellarAbsoluteMagnitudes'' galactic filter class.
    !!}
    implicit none
    type(galacticFilterStellarAbsoluteMagnitudes), intent(inout) :: self
    
    !![
    <objectDestructor name="self%galacticStructure_"/>
    !!]
    return
  end subroutine stellarAbsoluteMagnitudesDestructor

  logical function stellarAbsoluteMagnitudesPasses(self,node)
    !!{
    Implement a stellar absolute magnitude low-pass galactic filter.
    !!}
    use :: Galactic_Structure_Options    , only : massTypeStellar        , weightByLuminosity
    use :: Galacticus_Nodes              , only : nodeComponentBasic     , treeNode
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    class           (galacticFilterStellarAbsoluteMagnitudes), intent(inout)         :: self
    type            (treeNode                               ), intent(inout), target :: node
    class           (nodeComponentBasic                     ), pointer               :: basic
    double precision                                                                 :: time       , luminosity, &
         &                                                                              abMagnitude
    integer                                                                          :: iLuminosity

    ! Get the basic component.
    basic => node%basic()
    ! Get the time for this node.
    time  =  basic%time()
    ! Loop over all luminosities.
    do iLuminosity=1,unitStellarLuminosities%luminosityCount()
       ! Only check those luminosities which are being output at this output time.
       if (unitStellarLuminosities%isOutput(iLuminosity,time)) then
          ! Get the total stellar luminosity of the galaxy.
          luminosity=self%galacticStructure_%massEnclosed(node,massType=massTypeStellar,weightBy=weightByLuminosity,weightIndex=iLuminosity)
          ! Test only if the luminosity is greater than zero.
          if (luminosity > 0.0d0) then
             ! Convert to absolute magnitude.
             abMagnitude=-2.5d0*log10(luminosity)
             ! Filter out the galaxy if it is below the threshold.
             if (abMagnitude < self%absoluteMagnitudeThreshold(iLuminosity)) then
                stellarAbsoluteMagnitudesPasses=.false.
                return
             end if
          else
             stellarAbsoluteMagnitudesPasses=.false.
             return
          end if
       end if
    end do
    stellarAbsoluteMagnitudesPasses=.true.
    return
  end function stellarAbsoluteMagnitudesPasses
