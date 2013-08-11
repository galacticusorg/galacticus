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

!% Contains a module which implements calculations related to Type Ia supernovae.

module Supernovae_Type_Ia_Nagashima
  !% Implements calculations related to Type Ia supernovae.
  implicit none
  private
  public :: Supernovae_Type_Ia_Nagashima_Initialize

  ! Parameters of the distribution of binaries from Nagashima et al. (2005; MNRAS; 358; 1427; eqn. 17).
  double precision                            :: binaryMassMaximum=12.0d0, binaryMassMinimum  =3.00d0
  double precision                            :: gamma            =2.0d0 , typeIaNormalization=0.07d0

  ! Total yield of metals from Type Ia supernova.
  double precision                            :: totalYield
  double precision, allocatable, dimension(:) :: elementYield

contains

  !# <supernovaeIaMethod>
  !#  <unitName>Supernovae_Type_Ia_Nagashima_Initialize</unitName>
  !# </supernovaeIaMethod>
  subroutine Supernovae_Type_Ia_Nagashima_Initialize(supernovaeIaMethod,SNeIa_Cumulative_Number_Get,SNeIa_Cumulative_Yield_Get)
    !% Initialize the ``Nagashima'' Type Ia supernovae module.
    use ISO_Varying_String
    use Galacticus_Error
    use FoX_dom
    use IO_XML
    use Atomic_Data
    use Memory_Management
    use Galacticus_Input_Paths
    implicit none
    type            (varying_string                   ), intent(in   )          :: supernovaeIaMethod
    procedure       (SNeIa_Cumulative_Number_Nagashima), intent(inout), pointer :: SNeIa_Cumulative_Number_Get
    procedure       (SNeIa_Cumulative_Yield_Nagashima ), intent(inout), pointer :: SNeIa_Cumulative_Yield_Get
    type            (Node                             )               , pointer :: doc                        , thisAtom    , &
         &                                                                         thisIsotope                , thisYield
    type            (NodeList                         )               , pointer :: isotopesList
    integer                                                                     :: atomicIndex                , atomicNumber, &
         &                                                                         iIsotope                   , ioErr
    double precision                                                            :: isotopeYield

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
       doc => parseFile(char(Galacticus_Input_Path())//'data/stellarAstrophysics/Supernovae_Type_Ia_Yields.xml',iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Supernovae_Type_Ia_Nagashima_Initialize','Unable to parse yields file')

       ! Get a list of all isotopes.
       isotopesList => getElementsByTagname(doc,"isotope")

       ! Loop through isotopes and compute the net metal yield.
       do iIsotope=0,getLength(isotopesList)-1
          thisIsotope  => item(isotopesList,iIsotope)
          if (XML_Array_Length(thisIsotope,"yield") /= 1) call Galacticus_Error_Report('Supernovae_Type_Ia_Nagashima_Initialize' &
               & ,'isotope must have precisely one yield')
          thisYield => XML_Get_First_Element_By_Tag_Name(thisIsotope,"yield")
          call extractDataContent(thisYield,isotopeYield)
          totalYield=totalYield+isotopeYield
          if (XML_Array_Length(thisIsotope,"atomicNumber") /= 1) call Galacticus_Error_Report('Supernovae_Type_Ia_Nagashima_Initialize' &
               & ,'isotope must have precisely one atomic number')
          thisAtom => XML_Get_First_Element_By_Tag_Name(thisIsotope,"atomicNumber")
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
    double precision, intent(in   ) :: age          , initialMass, metallicity
    double precision                :: dyingStarMass, muMinimum

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
    implicit none
    double precision, intent(in   )           :: age      , initialMass, metallicity
    integer         , intent(in   ), optional :: atomIndex
    double precision                          :: yield

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
