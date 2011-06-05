!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements calculations related to Type Ia supernovae.

module Supernovae_Type_Ia_Nagashima
  !% Implements calculations related to Type Ia supernovae.
  private
  public :: Supernovae_Type_Ia_Nagashima_Initialize
  
  ! Parameters of the distribution of binaries from Nagashima et al. (2005; MNRAS; 358; 1427; eqn. 17).
  double precision :: binaryMassMinimum  =3.00d0, binaryMassMaximum=12.0d0
  double precision :: typeIaNormalization=0.07d0, gamma            = 2.0d0

  ! Total yield of metals from Type Ia supernova.
  double precision                            :: totalYield
  double precision, allocatable, dimension(:) :: elementYield

contains

  !# <supernovaeIaMethod>
  !#  <unitName>Supernovae_Type_Ia_Nagashima_Initialize</unitName>
  !# </supernovaeIaMethod>
  subroutine Supernovae_Type_Ia_Nagashima_Initialize(supernovaeIaMethod,SNeIa_Cumulative_Number_Get,SNeIa_Cumulative_Yield_Get)
    !% Initialize the ``Nagashima'' Type Ia supernovae module.
    use Numerical_Constants_Units
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use ISO_Varying_String
    use Galacticus_Error
    use FoX_dom
    use Atomic_Data
    use Memory_Management
    implicit none
    type(varying_string),                 intent(in)    :: supernovaeIaMethod
    procedure(double precision), pointer, intent(inout) :: SNeIa_Cumulative_Number_Get,SNeIa_Cumulative_Yield_Get
    type(Node),                  pointer                :: doc,thisIsotope,thisYield,thisAtom
    type(NodeList),              pointer                :: isotopesList,propertyList
    integer                                             :: iIsotope,ioErr,atomicNumber,atomicIndex
    double precision                                    :: isotopeYield

    if (supernovaeIaMethod == 'Nagashima') then
       ! Set up pointers to our procedures.
       SNeIa_Cumulative_Number_Get => SNeIa_Cumulative_Number_Nagashima
       SNeIa_Cumulative_Yield_Get  => SNeIa_Cumulative_Yield_Nagashima

       ! Allocate an array to store individual element yields.
       call Alloc_Array(elementYield,[Atomic_Data_Atoms_Count()])
       elementYield=0.0d0

       ! Read in Type Ia yields.
       !$omp critical (FoX_DOM_Access)
       ! Open the XML file containing yields.
       doc => parseFile('data/Supernovae_Type_Ia_Yields.xml',iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Supernovae_Type_Ia_Nagashima_Initialize','Unable to parse yields file')

       ! Get a list of all isotopes.
       isotopesList => getElementsByTagname(doc,"isotope")

       ! Loop through isotopes and compute the net metal yield.
       do iIsotope=0,getLength(isotopesList)-1
          thisIsotope  => item(isotopesList,iIsotope)
          propertyList => getElementsByTagname(thisIsotope,"yield")
          if (getLength(propertyList) /= 1) call Galacticus_Error_Report('Supernovae_Type_Ia_Nagashima_Initialize' &
               & ,'isotope must have precisely one yield')
          thisYield => item(propertyList,0)
          call extractDataContent(thisYield,isotopeYield)
          totalYield=totalYield+isotopeYield
          propertyList => getElementsByTagname(thisIsotope,"atomicNumber")
          if (getLength(propertyList) /= 1) call Galacticus_Error_Report('Supernovae_Type_Ia_Nagashima_Initialize' &
               & ,'isotope must have precisely one atomic number')
          thisAtom => item(propertyList,0)
          call extractDataContent(thisAtom,atomicNumber)
          atomicIndex=Atom_Lookup(atomicNumber=atomicNumber)
          elementYield(atomicIndex)=elementYield(atomicIndex)+isotopeYield
       end do

       ! Destroy the document.
       call destroy(doc)

       !$omp end critical (FoX_DOM_Access)
    end if
    return
  end subroutine Supernovae_Type_Ia_Nagashima_Initialize

  double precision function SNeIa_Cumulative_Number_Nagashima(initialMass,age,metallicity)
    !% Compute the cumulative number of Type Ia supernovae originating per unit mass of stars that form with given {\tt
    !% initialMass} and {\tt metallicity} after a time {\tt age}. The calculation is based on that of \cite{nagashima_metal_2005}. The
    !% number returned here assumes a distribution of binary mass ratios and so only makes sense once it is integrated over an initial
    !% mass function.
    use Stellar_Astrophysics
    implicit none
    double precision, intent(in) :: initialMass,age,metallicity
    double precision             :: dyingStarMass,muMinimum

    ! Check if initial mass is within the range of binary masses that lead to Type Ia supernovae.
    if (initialMass > binaryMassMinimum .and. initialMass < binaryMassMaximum) then
       
       ! Get the initial mass of a star which is just dying at this age.
       dyingStarMass=Star_Initial_Mass(age,metallicity)
       
       ! Compute the cumulative number of Type Ia supernovae originating from stars of this mass.
       muMinimum=max(dyingStarMass/initialMass,(1.0d0-binaryMassMaximum/2.0d0/initialMass))
       if (muMinimum < 0.5d0) then
          SNeIa_Cumulative_Number_Nagashima=typeIaNormalization*(1.0d0-(2.0d0*muMinimum)**(1.0d0+gamma))
       else
          SNeIa_Cumulative_Number_Nagashima=0.0d0      
       end if
      
    else
       ! Mass is not in range - assume that no Type Ia SNe are produced.
       SNeIa_Cumulative_Number_Nagashima=0.0d0
    end if
    return
  end function SNeIa_Cumulative_Number_Nagashima

  double precision function SNeIa_Cumulative_Yield_Nagashima(initialMass,age,metallicity,atomIndex)
    !% Compute the cumulative yield from Type Ia supernovae originating per unit mass of stars that form with given {\tt
    !% initialMass} and {\tt metallicity} after a time {\tt age}. The calculation is based on the Type Ia rate calculation of
    !% \cite{nagashima_metal_2005} and the Type Ia yields from \cite{nomoto_nucleosynthesis_1997}. The number returned here
    !% assumes a distribution of binary mass ratios and so only makes sense once it is integrated over an initial mass function.
    use Stellar_Astrophysics
    implicit none
    double precision, intent(in)           :: initialMass,age,metallicity
    integer,          intent(in), optional :: atomIndex
    double precision                       :: yield

    if (present(atomIndex)) then
       ! Return yield for requested atomic index.
       yield=elementYield(atomIndex)
    else
       ! No atomic index given, therefore return total metal yield.
       yield=totalYield
    end if
    SNeIa_Cumulative_Yield_Nagashima=SNeIa_Cumulative_Number_Nagashima(initialMass,age,metallicity)*yield
    return
  end function SNeIa_Cumulative_Yield_Nagashima

end module Supernovae_Type_Ia_Nagashima
