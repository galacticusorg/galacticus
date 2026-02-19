!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implements a Population III supernovae class based on \cite{heger_nucleosynthetic_2002}.
  !!}
  
  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Numerical_Interpolation, only : interpolator
  use            :: Stellar_Astrophysics   , only : stellarAstrophysics, stellarAstrophysicsClass

  !![
  <supernovaePopulationIII name="supernovaePopulationIIIHegerWoosley2002">
   <description>
    A Population III supernovae class that computes the energies of pair instability supernovae from the results of
    \cite{heger_nucleosynthetic_2002}.
   </description>
  </supernovaePopulationIII>
  !!]
  type, extends(supernovaePopulationIIIClass) :: supernovaePopulationIIIHegerWoosley2002
     !!{
     A Population III supernovae class based on \cite{heger_nucleosynthetic_2002}
     !!}
     private
     class           (stellarAstrophysicsClass), pointer                   :: stellarAstrophysics_ => null()
     integer         (c_size_t                )                            :: countTable
     double precision                          , allocatable, dimension(:) :: energy                        , massHeliumCore
     type            (interpolator            )                            :: interpolator
   contains
     final     ::                     hegerWoosley2002Destructor
     procedure :: energyCumulative => hegerWoosley2002EnergyCumulative
  end type supernovaePopulationIIIHegerWoosley2002

  interface supernovaePopulationIIIHegerWoosley2002
     !!{
     Constructors for the \refClass{supernovaePopulationIIIHegerWoosley2002} Population III supernovae class.
     !!}
     module procedure hegerWoosley2002ConstructorParameters
     module procedure hegerWoosley2002ConstructorInternal
  end interface supernovaePopulationIIIHegerWoosley2002

contains

  function hegerWoosley2002ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{supernovaePopulationIIIHegerWoosley2002} Population III supernovae class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (supernovaePopulationIIIHegerWoosley2002)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(stellarAstrophysicsClass               ), pointer       :: stellarAstrophysics_

    !![
    <objectBuilder class="stellarAstrophysics" name="stellarAstrophysics_" source="parameters"/>
    !!]
    self=supernovaePopulationIIIHegerWoosley2002(stellarAstrophysics_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="stellarAstrophysics_"/>
    !!]
    return
  end function hegerWoosley2002ConstructorParameters

  function hegerWoosley2002ConstructorInternal(stellarAstrophysics_) result(self)
    !!{
    Internal constructor for the \refClass{supernovaePopulationIIIHegerWoosley2002} Population III supernovae class.
    !!}
    use :: FoX_dom                         , only : destroy       , node                          , parseFile
    use :: Error                           , only : Error_Report
    use :: Input_Paths                     , only : inputPath     , pathTypeDataStatic
    use :: IO_XML                          , only : XML_Array_Read, XML_Count_Elements_By_Tag_Name, XML_Get_First_Element_By_Tag_Name
    use :: Numerical_Constants_Astronomical, only : massSolar
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Constants_Units       , only : ergs
    implicit none
    type   (supernovaePopulationIIIHegerWoosley2002)                         :: self
    class  (stellarAstrophysicsClass               ), intent(in   ), target  :: stellarAstrophysics_
    type   (node                                   )               , pointer :: doc                 , energyElement, &
         &                                                                      massElement
    integer                                                                  :: ioErr
    !![
    <constructorAssign variables="*stellarAstrophysics_"/>
    !!]

    ! Read in pair instability supernova energies.
    !$omp critical (FoX_DOM_Access)
    ! Open the XML file containing yields.
    doc => parseFile(char(inputPath(pathTypeDataStatic))//'stellarAstrophysics/Supernovae_Pair_Instability_Heger_Woosley_1992.xml',iostat=ioErr)
    if (ioErr /= 0) call Error_Report('Unable to parse supernovae file'//{introspection:location})
    ! Get the mass and energy elements.
    massElement   => XML_Get_First_Element_By_Tag_Name(doc,"heliumCoreMass" )
    energyElement => XML_Get_First_Element_By_Tag_Name(doc,"supernovaEnergy")
    ! Read the arrays.
    call XML_Array_Read(massElement  ,"data",self%massHeliumCore)
    call XML_Array_Read(energyElement,"data",self%energy        )
    self%countTable=XML_Count_Elements_By_Tag_Name(massElement,"data")
    ! Convert energies to M☉ (km/s)².
    self%energy=self%energy*(1.0d51*ergs/massSolar/kilo**2)
    ! Destroy the document.
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    self%interpolator=interpolator(self%massHeliumCore,self%energy)
    return
  end function hegerWoosley2002ConstructorInternal

  subroutine hegerWoosley2002Destructor(self)
    !!{
    Destructor for the \refClass{supernovaePopulationIIIHegerWoosley2002} Population III supernovae class.
    !!}
    implicit none
    type(supernovaePopulationIIIHegerWoosley2002), intent(inout) :: self

    !![
    <objectDestructor name="self%stellarAstrophysics_"/>
    !!]
    return
  end subroutine hegerWoosley2002Destructor

  double precision function hegerWoosley2002EnergyCumulative(self,initialMass,age,metallicity)
    !!{
    Compute the cumulative energy input from Population III star pair instability supernovae using the results of
    \cite{heger_nucleosynthetic_2002}.
    !!}
    implicit none
    class           (supernovaePopulationIIIHegerWoosley2002), intent(inout) :: self
    double precision                                         , intent(in   ) :: age        , initialMass   , &
         &                                                                      metallicity
    double precision                                                         :: lifetime   , massHeliumCore

    ! Get the lifetime of a star of this initial mass and metallicity.
    lifetime=self%stellarAstrophysics_%lifetime(initialMass,metallicity)
    ! Check if star has reached end of life.
    if (lifetime <= age) then
       ! Star has reached end of life. Compute core helium mass using simple scaling given by Heger & Woosley.
       massHeliumCore=(13.0d0/24.0d0)*(initialMass-20.0d0)
       ! Check if this is within the range tabulated.
       if     (                                                        &
            &   massHeliumCore >= self%massHeliumCore(1              ) &
            &  .and.                                                   &
            &   massHeliumCore <= self%massHeliumCore(self%countTable) &
            & ) then
          hegerWoosley2002EnergyCumulative=self%interpolator%interpolate(massHeliumCore)
       else
          hegerWoosley2002EnergyCumulative=0.0d0
       end if
    else
       ! Star has not gone supernova yet.
       hegerWoosley2002EnergyCumulative=0.0d0
    end if
    return
  end function hegerWoosley2002EnergyCumulative
