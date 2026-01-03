!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  Implements a supernovae type Ia class in which the yield is independent of progenitor.
  !!}
  
  !![
  <supernovaeTypeIa name="supernovaeTypeIaFixedYield" abstract="yes">
   <description>
    A supernovae type Ia class in which the yield is independent of progenitor.
   </description>
  </supernovaeTypeIa>
  !!]
  type, abstract, extends(supernovaeTypeIaClass) :: supernovaeTypeIaFixedYield
     !!{
     A supernovae type Ia class in which the yield is independent of progenitor.
     !!}
     private
     double precision                            :: totalYield
     double precision, allocatable, dimension(:) :: elementYield
     logical                                     :: initialized
   contains
     !![
     <methods>
       <method method="initialize" description="Initialize yield data."/>
     </methods>
     !!]
     procedure :: yield            => fixedYieldYield
     procedure :: initialize       => fixedYieldInitialize
  end type supernovaeTypeIaFixedYield

contains

  subroutine fixedYieldInitialize(self)
    !!{
    Read data for the {\normalfont \ttfamily fixedYield} supernovae type Ia class.
    !!}
    use :: Atomic_Data       , only : Atom_Lookup                   , Atomic_Data_Atoms_Count
    use :: FoX_dom           , only : destroy                       , node                             , extractDataContent
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath                     , pathTypeDataStatic
    use :: IO_XML            , only : XML_Count_Elements_By_Tag_Name, XML_Get_First_Element_By_Tag_Name, XML_Get_Elements_By_Tag_Name, xmlNodeList, &
         &                            XML_Parse
    use :: ISO_Varying_String, only : varying_string
    implicit none
    class           (supernovaeTypeIaFixedYield), intent(inout)               :: self
    type            (node                      ), pointer                     :: doc         , atom        , &
         &                                                                       isotope     , yield
    type            (xmlNodeList               ), allocatable  , dimension(:) :: isotopesList
    integer                                                                   :: atomicIndex , atomicNumber, &
         &                                                                       iIsotope    , ioErr
    double precision                                                          :: isotopeYield
    type            (varying_string            )                              :: fileName
    
    if (self%initialized) return
    ! Allocate an array to store individual element yields.
    allocate(self%elementYield(Atomic_Data_Atoms_Count()))
    self%elementYield=0.0d0
    self%totalYield  =0.0d0
    ! Read in Type Ia yields.
    fileName=char(inputPath(pathTypeDataStatic))//'stellarAstrophysics/Supernovae_Type_Ia_Yields.xml'
    !$omp critical (FoX_DOM_Access)
    doc => XML_Parse(char(fileName),iostat=ioErr)
    if (ioErr /= 0) call Error_Report('Unable to parse yields file'//{introspection:location})
    ! Get a list of all isotopes.
    call XML_Get_Elements_By_Tag_Name(doc,"isotope",isotopesList)
    ! Loop through isotopes and compute the net metal yield.
    do iIsotope=0,size(isotopesList)-1
       isotope  => isotopesList(iIsotope)%element
       if (XML_Count_Elements_By_Tag_Name(isotope,"yield") /= 1) call Error_Report('isotope must have precisely one yield'//{introspection:location})
       yield => XML_Get_First_Element_By_Tag_Name(isotope,"yield")
       call extractDataContent(yield,isotopeYield)
       self%totalYield=self%totalYield+isotopeYield
       if (XML_Count_Elements_By_Tag_Name(isotope,"atomicNumber") /= 1) call Error_Report('isotope must have precisely one atomic number'//{introspection:location})
       atom => XML_Get_First_Element_By_Tag_Name(isotope,"atomicNumber")
       call extractDataContent(atom,atomicNumber)
       atomicIndex=Atom_Lookup(atomicNumber=atomicNumber)
       self%elementYield(atomicIndex)=self%elementYield(atomicIndex)+isotopeYield
    end do
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    self%initialized=.true.
    return
  end subroutine fixedYieldInitialize

  double precision function fixedYieldYield(self,initialMassFunction_,initialMass,age,metallicity,atomIndex) result(yield)
    !!{
    Compute the cumulative yield from Type Ia supernovae originating per unit interval of secondary star mass with given
    {\normalfont \ttfamily initialMass} and {\normalfont \ttfamily metallicity} after a time {\normalfont \ttfamily age}. The
    calculation is based on that of \cite{nagashima_metal_2005} with Type Ia yields from \cite{nomoto_nucleosynthesis_1997}. The
    number returned here assumes a distribution of binary mass ratios and so only makes sense once it is integrated over an
    initial mass function.
    !!}
    implicit none
    class           (supernovaeTypeIaFixedYield), intent(inout)           :: self
    class           (initialMassFunctionClass  ), intent(inout)           :: initialMassFunction_
    double precision                            , intent(in   )           :: age                 , initialMass, &
         &                                                                   metallicity
    integer                                     , intent(in   ), optional :: atomIndex

    call self%initialize()
    if (present(atomIndex)) then
       ! Return yield for requested atomic index.
       yield=self%elementYield(atomIndex)
    else
       ! No atomic index given, therefore return total metal yield.
       yield=self%totalYield
    end if
    yield=self%number(initialMassFunction_,initialMass,age,metallicity)*yield
    return
  end function fixedYieldYield
