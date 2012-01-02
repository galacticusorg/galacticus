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


!% Contains a module which implements calculations related to Population III supernovae.

module Supernovae_Population_III_HegerWoosley
  !% Implements calculations related to Population III supernovae.
  implicit none
  private
  public :: Supernovae_Population_III_HegerWoosley_Initialize
  
  ! Variables holding the supernovae data tables.
  integer                                     :: supernovaeTableCount
  double precision, allocatable, dimension(:) :: supernovaeTableHeliumCoreMass,supernovaeTableEnergy

contains

  !# <supernovaePopIIIMethod>
  !#  <unitName>Supernovae_Population_III_HegerWoosley_Initialize</unitName>
  !# </supernovaePopIIIMethod>
  subroutine Supernovae_Population_III_HegerWoosley_Initialize(supernovaePopIIIMethod,SNePopIII_Cumulative_Energy_Get)
    !% Initialize the ``Heger-Woosley2002'' Population III supernovae module.
    use Numerical_Constants_Units
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use ISO_Varying_String
    use Galacticus_Error
    use FoX_dom
    use Memory_Management
    implicit none
    type(varying_string),                 intent(in)    :: supernovaePopIIIMethod
    procedure(double precision), pointer, intent(inout) :: SNePopIII_Cumulative_Energy_Get
    type(Node),                  pointer                :: doc,massElement,energyElement,thisDatum
    type(NodeList),              pointer                :: massList,energyList,massDataList,energyDataList
    integer                                             :: ioErr,iSupernovae

    if (supernovaePopIIIMethod == 'Heger-Woosley2002') then
       ! Set up pointers to our procedures.
       SNePopIII_Cumulative_Energy_Get => SNePopIII_Cumulative_Energy_HegerWoosley

       ! Read in pair instability supernova energies.
       !$omp critical (FoX_DOM_Access)
       ! Open the XML file containing yields.
       doc => parseFile('data/Supernovae_Pair_Instability_Heger_Woosley_1992.xml',iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Supernovae_Population_III_HegerWoosley_Initialize','Unable to parse supernovae file')

       ! Get the mass and energy elements.
       massList       => getElementsByTagname(doc,"heliumCoreMass")
       massElement    => item(massList,0)
       massDataList   => getElementsByTagname(massElement,"data")
       energyList     => getElementsByTagname(doc,"supernovaEnergy")
       energyElement  => item(massList,0)
       energyDataList => getElementsByTagname(energyElement,"data")

       ! Count how many elements are present and allocate arrays.
       supernovaeTableCount=getLength(massDataList)
       call Alloc_Array(supernovaeTableHeliumCoreMass,[supernovaeTableCount])
       call Alloc_Array(supernovaeTableEnergy        ,[supernovaeTableCount])

       ! Loop through isotopes and compute the net metal yield.
       do iSupernovae=0,getLength(massDataList)-1
          thisDatum => item(massDataList  ,iSupernovae)
          call extractDataContent(thisDatum,supernovaeTableHeliumCoreMass(iSupernovae+1))
          thisDatum => item(energyDataList,iSupernovae)
          call extractDataContent(thisDatum,supernovaeTableEnergy        (iSupernovae+1))
       end do

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
    double precision,        intent(in) :: initialMass,age,metallicity
    double precision                    :: lifetime,massHeliumCore
    type(fgsl_interp),       save       :: interpolationObject
    type(fgsl_interp_accel), save       :: interpolationAccelerator
    logical,                 save       :: interpolationReset=.true.
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
