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

  !!{RST
  Implements a supernovae type Ia class in which the yield is independent of progenitor.
  !!}
  
  !![
  <supernovaeTypeIa name="supernovaeTypeIaFixedYield" abstract="yes" docformat="rst">
   <description>
   A supernovae type Ia class in which the yield is independent of progenitor mass. The yields are read from an XML file, which
   may optionally tabulate them at several metallicities, in which case they are interpolated linearly in metallicity (and held
   constant beyond the tabulated range). The file should have the following structure:

   .. code-block:: none

       &lt;supernovaeYields&gt;
        &lt;isotope&gt;
          &lt;name&gt;56Fe&lt;/name&gt;
          &lt;atomicMass&gt;56&lt;/atomicMass&gt;
          &lt;atomicNumber&gt;26&lt;/atomicNumber&gt;
          &lt;yield&gt;6.26e-1&lt;/yield&gt;
        &lt;/isotope&gt;
        .
        .
        .
       &lt;/supernovaeYields&gt;

   Only metals should be included---the total metal yield is formed by summing over every isotope present in the file. To make the
   yields depend on metallicity, wrap each set of isotopes in a ``yieldsMetallicity`` element carrying the metallicity at which
   that set applies:

   .. code-block:: none

       &lt;supernovaeYields&gt;
        &lt;yieldsMetallicity&gt;
          &lt;metallicity&gt;1.49e-4&lt;/metallicity&gt;
          &lt;isotope&gt;
            .
            .
          &lt;/isotope&gt;
        &lt;/yieldsMetallicity&gt;
        &lt;yieldsMetallicity&gt;
          &lt;metallicity&gt;1.49e-2&lt;/metallicity&gt;
          .
          .
        &lt;/yieldsMetallicity&gt;
       &lt;/supernovaeYields&gt;
   </description>
  </supernovaeTypeIa>
  !!]
  type, abstract, extends(supernovaeTypeIaClass) :: supernovaeTypeIaFixedYield
     !!{RST
     A supernovae type Ia class in which the yield is independent of progenitor mass.
     !!}
     private
     type            (varying_string)              :: fileName
     double precision, allocatable, dimension(:  ) :: totalYield  , metallicity
     double precision, allocatable, dimension(:,:) :: elementYield
     logical                                       :: initialized , metallicityDependent
   contains
     !![
     <methods docformat="rst">
       <method method="initialize" description="Initialize yield data."/>
     </methods>
     !!]
     procedure :: yield      => fixedYieldYield
     procedure :: initialize => fixedYieldInitialize
  end type supernovaeTypeIaFixedYield

contains

  subroutine fixedYieldInitialize(self)
    !!{RST
    Read data for the ``fixedYield`` supernovae type Ia class.
    !!}
    use            :: Atomic_Data       , only : Atom_Lookup                   , Atomic_Data_Atoms_Count
    use            :: FoX_dom           , only : destroy                       , node                             , extractDataContent
    use            :: Error             , only : Error_Report
    use            :: IO_XML            , only : XML_Count_Elements_By_Tag_Name, XML_Get_First_Element_By_Tag_Name, XML_Get_Elements_By_Tag_Name, xmlNodeList, &
         &                                       XML_Parse
    use            :: Sorting           , only : sortIndex
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    implicit none
    class           (supernovaeTypeIaFixedYield), intent(inout)               :: self
    type            (node                      ), pointer                     :: doc          , atom             , &
         &                                                                       isotope      , yield            , &
         &                                                                       metallicity  , yieldsContainer
    type            (xmlNodeList               ), allocatable  , dimension(:) :: isotopesList , metallicitiesList
    integer                                                                   :: atomicIndex  , atomicNumber     , &
         &                                                                       iIsotope     , ioErr            , &
         &                                                                       iMetallicity , countMetallicity
    double precision                                                          :: isotopeYield
    double precision                            , allocatable, dimension(:  ) :: metallicities
    double precision                            , allocatable, dimension(:,:) :: elementYields
    double precision                            , allocatable, dimension(:  ) :: totalYields
    integer         (c_size_t                  ), allocatable, dimension(:  ) :: order

    if (self%initialized) return
    ! Read in Type Ia yields.
    !$omp critical (FoX_DOM_Access)
    doc => XML_Parse(self%fileName,iostat=ioErr)
    if (ioErr /= 0) call Error_Report('Unable to parse yields file: "'//char(self%fileName)//'"'//{introspection:location})
    ! Determine whether yields are tabulated as a function of metallicity. If no `yieldsMetallicity` elements are
    ! present the file carries a single, metallicity-independent set of yields.
    call XML_Get_Elements_By_Tag_Name(doc,"yieldsMetallicity",metallicitiesList)
    self%metallicityDependent=    size(metallicitiesList) > 0
    countMetallicity         =max(size(metallicitiesList),1)
    allocate(metallicities(                          countMetallicity))
    allocate(totalYields  (                          countMetallicity))
    allocate(elementYields(Atomic_Data_Atoms_Count(),countMetallicity))
    elementYields=0.0d0
    totalYields  =0.0d0
    metallicities=0.0d0
    do iMetallicity=1,countMetallicity
       if (self%metallicityDependent) then
          yieldsContainer => metallicitiesList(iMetallicity-1)%element
          if (XML_Count_Elements_By_Tag_Name(yieldsContainer,"metallicity") /= 1) &
               & call Error_Report('yieldsMetallicity must have precisely one metallicity'//{introspection:location})
          metallicity => XML_Get_First_Element_By_Tag_Name(yieldsContainer,"metallicity")
          call extractDataContent(metallicity,metallicities(iMetallicity))
       else
          yieldsContainer => doc
       end if
       ! Get a list of all isotopes for this metallicity.
       call XML_Get_Elements_By_Tag_Name(yieldsContainer,"isotope",isotopesList)
       if (size(isotopesList) == 0) call Error_Report('no isotopes found in yields file'//{introspection:location})
       ! Loop through isotopes and compute the net metal yield.
       do iIsotope=0,size(isotopesList)-1
          isotope  => isotopesList(iIsotope)%element
          if (XML_Count_Elements_By_Tag_Name(isotope,"yield") /= 1) call Error_Report('isotope must have precisely one yield'//{introspection:location})
          yield => XML_Get_First_Element_By_Tag_Name(isotope,"yield")
          call extractDataContent(yield,isotopeYield)
          totalYields(iMetallicity)=totalYields(iMetallicity)+isotopeYield
          if (XML_Count_Elements_By_Tag_Name(isotope,"atomicNumber") /= 1) call Error_Report('isotope must have precisely one atomic number'//{introspection:location})
          atom => XML_Get_First_Element_By_Tag_Name(isotope,"atomicNumber")
          call extractDataContent(atom,atomicNumber)
          atomicIndex=Atom_Lookup(atomicNumber=atomicNumber)
          elementYields(atomicIndex,iMetallicity)=elementYields(atomicIndex,iMetallicity)+isotopeYield
       end do
    end do
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    ! Sort into order of increasing metallicity, as the interpolation assumes this.
    order=sortIndex(metallicities)
    allocate(self%metallicity (                          countMetallicity))
    allocate(self%totalYield  (                          countMetallicity))
    allocate(self%elementYield(Atomic_Data_Atoms_Count(),countMetallicity))
    do iMetallicity=1,countMetallicity
       self%metallicity (  iMetallicity)=metallicities(  order(iMetallicity))
       self%totalYield  (  iMetallicity)=totalYields  (  order(iMetallicity))
       self%elementYield(:,iMetallicity)=elementYields(:,order(iMetallicity))
    end do
    if (self%metallicityDependent .and. countMetallicity > 1) then
       do iMetallicity=2,countMetallicity
          if (self%metallicity(iMetallicity) <= self%metallicity(iMetallicity-1)) &
               & call Error_Report('yields file contains repeated metallicities'//{introspection:location})
       end do
    end if
    self%initialized=.true.
    return
  end subroutine fixedYieldInitialize

  double precision function fixedYieldInterpolate(metallicities,yields,metallicity) result(yieldInterpolated)
    !!{RST
    Interpolate a tabulated yield linearly in metallicity, holding it constant beyond the tabulated range.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: metallicities, yields
    double precision, intent(in   )               :: metallicity
    integer                                       :: iMetallicity
    double precision                              :: fraction

    ! With a single tabulated metallicity---or a file which does not tabulate metallicity at all---the yield is
    ! independent of metallicity.
    if (size(yields) == 1) then
       yieldInterpolated=yields(1)
       return
    end if
    ! Hold the yield fixed beyond the tabulated range rather than extrapolating, which for yields varying
    ! steeply with metallicity could otherwise produce negative values.
    if      (metallicity <= metallicities(1           )) then
       yieldInterpolated=yields(1           )
    else if (metallicity >= metallicities(size(yields))) then
       yieldInterpolated=yields(size(yields))
    else
       do iMetallicity=1,size(yields)-1
          if (metallicity <= metallicities(iMetallicity+1)) then
             fraction         =+(metallicity                     -metallicities(iMetallicity)) &
                  &            /(metallicities(iMetallicity+1)   -metallicities(iMetallicity))
             yieldInterpolated=+yields(iMetallicity  )*(1.0d0-fraction) &
                  &            +yields(iMetallicity+1)*       fraction
             return
          end if
       end do
       yieldInterpolated=yields(size(yields))
    end if
    return
  end function fixedYieldInterpolate

  double precision function fixedYieldYield(self,initialMassFunction_,initialMass,age,metallicity,atomIndex) result(yield)
    !!{RST
    Compute the cumulative yield from Type Ia supernovae originating per unit interval of secondary star mass with given ``initialMass`` and ``metallicity`` after a time ``age``. Yields are read from the file named by the ``fileName`` parameter---by default those of :cite:t:`nomoto_nucleosynthesis_1997`---and are interpolated in metallicity if that file tabulates them at more than one metallicity. The number returned here assumes a distribution of binary mass ratios and so only makes sense once it is integrated over an initial mass function.
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
       yield=fixedYieldInterpolate(self%metallicity,self%elementYield(atomIndex,:),metallicity)
    else
       ! No atomic index given, therefore return total metal yield.
       yield=fixedYieldInterpolate(self%metallicity,self%totalYield               ,metallicity)
    end if
    yield=self%number(initialMassFunction_,initialMass,age,metallicity)*yield
    return
  end function fixedYieldYield
