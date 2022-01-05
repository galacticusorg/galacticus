!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  Implements a supernovae type Ia class based on \cite{nagashima_metal_2005}.
  !!}

  use :: Stellar_Astrophysics, only : stellarAstrophysicsClass

  !![
  <supernovaeTypeIa name="supernovaeTypeIaNagashima2005">
   <description>
    A supernovae type Ia class which uses the prescriptions from \cite{nagashima_metal_2005} to compute the numbers and yields
    of Type Ia supernovae.
   </description>
  </supernovaeTypeIa>
  !!]
  type, extends(supernovaeTypeIaClass) :: supernovaeTypeIaNagashima2005
     !!{
     A supernovae type Ia class based on \cite{nagashima_metal_2005}.
     !!}
     private
     class           (stellarAstrophysicsClass), pointer                   :: stellarAstrophysics_ => null()
     double precision                                                      :: totalYield
     double precision                          , allocatable, dimension(:) :: elementYield
   contains
     final     ::           nagashima2005Destructor
     procedure :: number => nagashima2005Number
     procedure :: yield  => nagashima2005Yield
  end type supernovaeTypeIaNagashima2005

  interface supernovaeTypeIaNagashima2005
     !!{
     Constructors for the {\normalfont \ttfamily nagashima2005} supernovae type Ia class.
     !!}
     module procedure nagashima2005ConstructorParameters
     module procedure nagashima2005ConstructorInternal
  end interface supernovaeTypeIaNagashima2005

contains

  function nagashima2005ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily nagashima2005} supernovae type Ia class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (supernovaeTypeIaNagashima2005)                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(stellarAstrophysicsClass     ), pointer       :: stellarAstrophysics_

    !![
    <objectBuilder class="stellarAstrophysics" name="stellarAstrophysics_" source="parameters"/>
    !!]
    self=supernovaeTypeIaNagashima2005(stellarAstrophysics_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="stellarAstrophysics_"/>
    !!]
    return
  end function nagashima2005ConstructorParameters

  function nagashima2005ConstructorInternal(stellarAstrophysics_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily nagashima2005} supernovae type Ia class.
    !!}
    use :: Atomic_Data      , only : Atom_Lookup                   , Atomic_Data_Atoms_Count
    use :: FoX_dom          , only : destroy                       , node
    use :: Galacticus_Error , only : Galacticus_Error_Report
    use :: Galacticus_Paths , only : galacticusPath                , pathTypeDataStatic
    use :: IO_XML           , only : XML_Count_Elements_By_Tag_Name, XML_Get_First_Element_By_Tag_Name                        , XML_Get_Elements_By_Tag_Name, xmlNodeList, &
         &                           XML_Parse                     , extractDataContent                => extractDataContentTS
    use :: Memory_Management, only : allocateArray
    implicit none
    type            (supernovaeTypeIaNagashima2005)                              :: self
    class           (stellarAstrophysicsClass     ), intent(in   ), target       :: stellarAstrophysics_
    type            (node                         ), pointer                     :: doc                 , atom        , &
         &                                                                          isotope             , yield
    type            (xmlNodeList                  ), allocatable  , dimension(:) :: isotopesList
    integer                                                                      :: atomicIndex         , atomicNumber, &
         &                                                                          iIsotope            , ioErr
    double precision                                                             :: isotopeYield
    !![
    <constructorAssign variables="*stellarAstrophysics_"/>
    !!]

    ! Allocate an array to store individual element yields.
    call allocateArray(self%elementYield,[Atomic_Data_Atoms_Count()])
    self%elementYield=0.0d0
    self%totalYield  =0.0d0
    ! Read in Type Ia yields.
    !$omp critical (FoX_DOM_Access)
    ! Open the XML file containing yields.
    doc => XML_Parse(char(galacticusPath(pathTypeDataStatic))//'stellarAstrophysics/Supernovae_Type_Ia_Yields.xml',iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('Unable to parse yields file'//{introspection:location})
    ! Get a list of all isotopes.
    call XML_Get_Elements_By_Tag_Name(doc,"isotope",isotopesList)
    ! Loop through isotopes and compute the net metal yield.
    do iIsotope=0,size(isotopesList)-1
       isotope  => isotopesList(iIsotope)%element
       if (XML_Count_Elements_By_Tag_Name(isotope,"yield") /= 1) call Galacticus_Error_Report('isotope must have precisely one yield'//{introspection:location})
       yield => XML_Get_First_Element_By_Tag_Name(isotope,"yield")
       call extractDataContent(yield,isotopeYield)
       self%totalYield=self%totalYield+isotopeYield
       if (XML_Count_Elements_By_Tag_Name(isotope,"atomicNumber") /= 1) call Galacticus_Error_Report('isotope must have precisely one atomic number'//{introspection:location})
       atom => XML_Get_First_Element_By_Tag_Name(isotope,"atomicNumber")
       call extractDataContent(atom,atomicNumber)
       atomicIndex=Atom_Lookup(atomicNumber=atomicNumber)
       self%elementYield(atomicIndex)=self%elementYield(atomicIndex)+isotopeYield
    end do
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    return
  end function nagashima2005ConstructorInternal

  subroutine nagashima2005Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily nagashima2005} supernovae type Ia class.
    !!}
    implicit none
    type(supernovaeTypeIaNagashima2005), intent(inout) :: self

    !![
    <objectDestructor name="self%stellarAstrophysics_"/>
    !!]
    return
  end subroutine nagashima2005Destructor

  double precision function nagashima2005Number(self,initialMass,age,metallicity)
    !!{
    Compute the cumulative number of Type Ia supernovae originating per unit mass of stars that form with given {\normalfont \ttfamily
    initialMass} and {\normalfont \ttfamily metallicity} after a time {\normalfont \ttfamily age}. The calculation is based on that of \cite{nagashima_metal_2005}. The
    number returned here assumes a distribution of binary mass ratios and so only makes sense once it is integrated over an initial
    mass function.
    !!}
    implicit none
    class           (supernovaeTypeIaNagashima2005), intent(inout) :: self
    double precision                               , intent(in   ) :: age          , initialMass, metallicity
    double precision                                               :: dyingStarMass, muMinimum
    ! Parameters of the distribution of binaries from Nagashima et al. (2005; MNRAS; 358; 1427; eqn. 17).
    double precision                               , parameter     :: binaryMassMaximum=12.0d0, binaryMassMinimum  =3.00d0
    double precision                               , parameter     :: gamma            = 2.0d0, typeIaNormalization=0.07d0

    ! Check if initial mass is within the range of binary masses that lead to Type Ia supernovae.
    if (initialMass > binaryMassMinimum .and. initialMass < binaryMassMaximum) then
       ! Get the initial mass of a star which is just dying at this age.
       dyingStarMass=self%stellarAstrophysics_%massInitial(age,metallicity)
       ! Compute the cumulative number of Type Ia supernovae originating from stars of this mass.
       muMinimum=max(dyingStarMass/initialMass,(1.0d0-binaryMassMaximum/2.0d0/initialMass))
       if (muMinimum < 0.5d0) then
          nagashima2005Number=typeIaNormalization*(1.0d0-(2.0d0*muMinimum)**(1.0d0+gamma))
       else
          nagashima2005Number=0.0d0
       end if
    else
       ! Mass is not in range - assume that no Type Ia SNe are produced.
       nagashima2005Number=0.0d0
    end if
    return
  end function nagashima2005Number

  double precision function nagashima2005Yield(self,initialMass,age,metallicity,atomIndex)
    !!{
    Compute the cumulative yield from Type Ia supernovae originating per unit mass of stars that form with given {\normalfont
    \ttfamily initialMass} and {\normalfont \ttfamily metallicity} after a time {\normalfont \ttfamily age}. The calculation is
    based on the Type Ia rate calculation of \cite{nagashima_metal_2005} and the Type Ia yields from
    \cite{nomoto_nucleosynthesis_1997}. The number returned here assumes a distribution of binary mass ratios and so only makes
    sense once it is integrated over an initial mass function.
    !!}
    implicit none
    class           (supernovaeTypeIaNagashima2005), intent(inout)           :: self
    double precision                               , intent(in   )           :: age      , initialMass, metallicity
    integer                                        , intent(in   ), optional :: atomIndex
    double precision                                                         :: yield

    if (present(atomIndex)) then
       ! Return yield for requested atomic index.
       yield=self%elementYield(atomIndex)
    else
       ! No atomic index given, therefore return total metal yield.
       yield=self%totalYield
    end if
    nagashima2005Yield=self%number(initialMass,age,metallicity)*yield
    return
  end function nagashima2005Yield
