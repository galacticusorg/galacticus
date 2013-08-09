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

!% Contains a module which implements calculations related to Population III supernovae.

module Supernovae_Population_III_HegerWoosley
  !% Implements calculations related to Population III supernovae.
  implicit none
  private
  public :: Supernovae_Population_III_HegerWoosley_Initialize

  ! Variables holding the supernovae data tables.
  integer                                     :: supernovaeTableCount
  double precision, allocatable, dimension(:) :: supernovaeTableEnergy, supernovaeTableHeliumCoreMass

contains

  !# <supernovaePopIIIMethod>
  !#  <unitName>Supernovae_Population_III_HegerWoosley_Initialize</unitName>
  !# </supernovaePopIIIMethod>
  subroutine Supernovae_Population_III_HegerWoosley_Initialize(supernovaePopIIIMethod,SNePopIII_Cumulative_Energy_Get)
    !% Initialize the ``Heger-Woosley2002'' Population III supernovae module.
    use Numerical_Constants_Astronomical
    use ISO_Varying_String
    use Galacticus_Error
    use FoX_dom
    use IO_XML
    implicit none
    type     (varying_string  ), intent(in   )          :: supernovaePopIIIMethod
    procedure(double precision), intent(inout), pointer :: SNePopIII_Cumulative_Energy_Get
    type     (Node            )               , pointer :: doc                            , energyElement, &
         &                                                 massElement
    integer                                             :: ioErr

    if (supernovaePopIIIMethod == 'Heger-Woosley2002') then
       ! Set up pointers to our procedures.
       SNePopIII_Cumulative_Energy_Get => SNePopIII_Cumulative_Energy_HegerWoosley

       ! Read in pair instability supernova energies.
       !$omp critical (FoX_DOM_Access)
       ! Open the XML file containing yields.
       doc => parseFile('data/stellarAstrophysics/Supernovae_Pair_Instability_Heger_Woosley_1992.xml',iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Supernovae_Population_III_HegerWoosley_Initialize','Unable to parse supernovae file')

       ! Get the mass and energy elements.
       massElement    => XML_Get_First_Element_By_Tag_Name(doc,"heliumCoreMass" )
       energyElement  => XML_Get_First_Element_By_Tag_Name(doc,"supernovaEnergy")

       ! Read the arrays.
       call XML_Array_Read(massElement  ,"data",supernovaeTableHeliumCoreMass)
       call XML_Array_Read(energyElement,"data",supernovaeTableEnergy        )
       supernovaeTableCount=XML_Array_Length(massElement,"data")

       ! Convert energies to MSolar (km/s)^2.
       supernovaeTableEnergy=supernovaeTableEnergy*(1.0d51*ergs/massSolar/kilo**2)

       ! Destroy the document.
       call destroy(doc)

       !$omp end critical (FoX_DOM_Access)
    end if
    return
  end subroutine Supernovae_Population_III_HegerWoosley_Initialize

  double precision function SNePopIII_Cumulative_Energy_HegerWoosley(initialMass,age,metallicity)
    !% Compute the cumulative energy input from Population III star pair instability supernovae using the results of
    !% \cite{heger_nucleosynthetic_2002}.
    use Stellar_Astrophysics
    use Numerical_Interpolation
    use FGSL
    implicit none
    double precision                   , intent(in   ) :: age                            , initialMass   , &
         &                                                metallicity
    double precision                                   :: lifetime                       , massHeliumCore
    type            (fgsl_interp      ), save          :: interpolationObject
    type            (fgsl_interp_accel), save          :: interpolationAccelerator
    logical                            , save          :: interpolationReset      =.true.
    !$omp threadprivate(interpolationObject,interpolationAccelerator,interpolationReset)
    ! Get the lifetime of a star of this initial mass and metallicity.
    lifetime=Star_Lifetime(initialMass,metallicity)

    ! Check if star has reached end of life.
    if (lifetime <= age) then
       ! Star has reached end of life. Compute core helium mass using simple scaling given by Heger & Woosley.
       massHeliumCore=(13.0d0/24.0d0)*(initialMass-20.0d0)
       ! Check if this is within the range tabulated.
       if (massHeliumCore >= supernovaeTableHeliumCoreMass(1) .and. massHeliumCore <= supernovaeTableHeliumCoreMass(supernovaeTableCount)) then
          SNePopIII_Cumulative_Energy_HegerWoosley=Interpolate(supernovaeTableCount,supernovaeTableHeliumCoreMass&
               &,supernovaeTableEnergy,interpolationObject,interpolationAccelerator,massHeliumCore ,reset=interpolationReset)
       else
          SNePopIII_Cumulative_Energy_HegerWoosley=0.0d0
       end if
    else
       ! Star has not gone supernova yet.
       SNePopIII_Cumulative_Energy_HegerWoosley=0.0d0
    end if

    return
  end function SNePopIII_Cumulative_Energy_HegerWoosley

end module Supernovae_Population_III_HegerWoosley
