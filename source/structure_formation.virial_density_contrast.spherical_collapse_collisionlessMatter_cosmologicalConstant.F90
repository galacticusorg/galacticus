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
  An implementation of dark matter halo virial density contrasts based on spherical collapse in a matter plus cosmological constant universe.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctions                              , cosmologyFunctionsClass
  use :: Spherical_Collapse_Solvers, only : sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt
  use :: Tables                    , only : table1D

  !![
  <virialDensityContrast name="virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt">
   <description>
    A class implementing dark matter halo virial density contrasts based on the spherical collapse model in a universe which
    contains collisionless matter and a cosmological constant (see, for example, \citealt{percival_cosmological_2005}).
   </description>
  </virialDensityContrast>
  !!]
  type, extends(virialDensityContrastClass) :: virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
     !!{
     A dark matter halo virial density contrast class based on spherical collapse in a matter plus cosmological constant universe.
     !!}
     private
     logical                                                                         :: tableInitialized        =  .false., turnaroundInitialized=.false.
     double precision                                                                :: tableTimeMinimum                  , tableTimeMaximum             , &
          &                                                                             turnaroundTimeMinimum             , turnaroundTimeMaximum
     logical                                                                         :: tableStore
     class           (table1D                                         ), allocatable :: deltaVirial                       , turnaround
     class           (cosmologyFunctionsClass                         ), pointer     :: cosmologyFunctions_      => null()
     class           (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt), pointer     :: sphericalCollapseSolver_ => null()
   contains
     !![
     <methods>
       <method description="Tabulate spherical collapse virial density contrast." method="retabulate" />
     </methods>
     !!]
     final     ::                                sphericalCollapseClsnlssMttrCsmlgclCnstntDestructor
     procedure :: densityContrast             => sphericalCollapseClsnlssMttrCsmlgclCnstntDensityContrast
     procedure :: densityContrastRateOfChange => sphericalCollapseClsnlssMttrCsmlgclCnstntDensityContrastRtChng
     procedure :: turnAroundOverVirialRadii   => sphericalCollapseClsnlssMttrCsmlgclCnstntTrnrndVrlRd
     procedure :: retabulate                  => sphericalCollapseClsnlssMttrCsmlgclCnstntRetabulate
  end type virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt

  interface virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
     !!{
     Constructors for the \refClass{virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt} dark matter halo virial density contrast class.
     !!}
     module procedure sphericalCollapseClsnlssMttrCsmlgclCnstntConstructorParameters
     module procedure sphericalCollapseClsnlssMttrCsmlgclCnstntConstructorInternal
  end interface virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt

contains

  function sphericalCollapseClsnlssMttrCsmlgclCnstntConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt} dark matter halo virial density contrast class that takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)                :: self
    type   (inputParameters                                               ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass                                       ), pointer       :: cosmologyFunctions_
    logical                                                                                :: tableStore

    !![
    <inputParameter>
      <name>tableStore</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, store/restore the tabulated solution to/from file when possible.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(tableStore,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function sphericalCollapseClsnlssMttrCsmlgclCnstntConstructorParameters

  function sphericalCollapseClsnlssMttrCsmlgclCnstntConstructorInternal(tableStore,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt} dark matter halo virial density contrast class.
    !!}
    implicit none
    type   (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)                        :: self
    class  (cosmologyFunctionsClass                                       ), intent(in   ), target :: cosmologyFunctions_
    logical                                                                , intent(in   )         :: tableStore
    !![
    <constructorAssign variables="tableStore, *cosmologyFunctions_"/>
    !!]

    self%tableInitialized     =.false.
    self%turnaroundInitialized=.false.
    allocate(sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt :: self%sphericalCollapseSolver_)
    select type (sphericalCollapseSolver_ => self%sphericalCollapseSolver_)
    type is (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt)
       !![
       <referenceConstruct isResult="yes" owner="self" nameAssociated="sphericalCollapseSolver_" object="sphericalCollapseSolver_" constructor="sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt(self%cosmologyFunctions_)"/>
       !!]
    end select
    return
  end function sphericalCollapseClsnlssMttrCsmlgclCnstntConstructorInternal

  subroutine sphericalCollapseClsnlssMttrCsmlgclCnstntDestructor(self)
    !!{
    Destructor for the \refClass{virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt} dark matter halo virial density contrast class.
    !!}
    implicit none
    type (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), intent(inout) :: self

    if (self%tableInitialized) then
       call self%deltaVirial%destroy()
       deallocate(self%deltaVirial)
    end if
    if (self%turnaroundInitialized) then
       call self%turnaround%destroy()
       deallocate(self%turnaround)
    end if
    !![
    <objectDestructor name="self%cosmologyFunctions_"     />
    <objectDestructor name="self%sphericalCollapseSolver_"/>
    !!]
    return
  end subroutine sphericalCollapseClsnlssMttrCsmlgclCnstntDestructor

  subroutine sphericalCollapseClsnlssMttrCsmlgclCnstntRetabulate(self,time)
    !!{
    Recompute the look-up tables for virial density contrast.
    !!}
    implicit none
    class           (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), intent(inout) :: self
    double precision                                                                , intent(in   ) :: time
    logical                                                                                         :: remakeTable

    ! Check if we need to recompute our table.
    if (self%tableInitialized) then
       remakeTable=(time < self%tableTimeMinimum .or. time > self%tableTimeMaximum)
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call self%sphericalCollapseSolver_%virialDensityContrast(time,self%tableStore,self%deltaVirial)
       self%tableInitialized=.true.
       self%tableTimeMinimum=self%deltaVirial%x(+1)
       self%tableTimeMaximum=self%deltaVirial%x(-1)
    end if
    return
  end subroutine sphericalCollapseClsnlssMttrCsmlgclCnstntRetabulate

  double precision function sphericalCollapseClsnlssMttrCsmlgclCnstntDensityContrast(self,mass,time,expansionFactor,collapsing)
    !!{
    Return the virial density contrast at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), intent(inout)           :: self
    double precision                                                                , intent(in   )           :: mass
    double precision                                                                , intent(in   ), optional :: time            , expansionFactor
    logical                                                                         , intent(in   ), optional :: collapsing
    logical                                                                                                   :: collapsingActual
    double precision                                                                                          :: timeActual
    !$GLC attributes unused :: mass

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Error_Report('only one argument can be specified'//{introspection:location})
       else
          timeActual=time
       end if
    else
       if (present(expansionFactor)) then
          if (present(collapsing)) then
             collapsingActual=collapsing
          else
             collapsingActual=.false.
          end if
          timeActual=self%cosmologyFunctions_%cosmicTime(expansionFactor,collapsingActual)
       else
          call Error_Report('at least one argument must be given'//{introspection:location})
       end if
    end if
    ! Remake the table if necessary.
    call self%retabulate(timeActual)
    ! Interpolate to get the density contrast.
    sphericalCollapseClsnlssMttrCsmlgclCnstntDensityContrast=self%deltaVirial%interpolate(timeActual)
    return
  end function sphericalCollapseClsnlssMttrCsmlgclCnstntDensityContrast

  double precision function sphericalCollapseClsnlssMttrCsmlgclCnstntDensityContrastRtChng(self,mass,time,expansionFactor,collapsing)
    !!{
    Return the virial density contrast at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), intent(inout)           :: self
    double precision                                                                , intent(in   )           :: mass
    double precision                                                                , intent(in   ), optional :: time            , expansionFactor
    logical                                                                         , intent(in   ), optional :: collapsing
    logical                                                                                                   :: collapsingActual
    double precision                                                                                          :: timeActual
    !$GLC attributes unused :: mass

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Error_Report('only one argument can be specified'//{introspection:location})
       else
          timeActual=time
       end if
    else
       if (present(expansionFactor)) then
          if (present(collapsing)) then
             collapsingActual=collapsing
          else
             collapsingActual=.false.
          end if
          timeActual=self%cosmologyFunctions_%cosmicTime(expansionFactor,collapsingActual)
       else
          call Error_Report('at least one argument must be given'//{introspection:location})
       end if
    end if
    ! Remake the table if necessary.
    call self%retabulate(timeActual)
    ! Interpolate to get the expansion factor.
    sphericalCollapseClsnlssMttrCsmlgclCnstntDensityContrastRtChng=self%deltaVirial%interpolateGradient(timeActual)
   return
  end function sphericalCollapseClsnlssMttrCsmlgclCnstntDensityContrastRtChng

  double precision function sphericalCollapseClsnlssMttrCsmlgclCnstntTrnrndVrlRd(self,mass,time,expansionFactor,collapsing)
    !!{
    Return the ratio of turnaround and virial radii at the given epoch, based spherical collapse in a matter plus cosmological
    constant universe.
    !!}
    implicit none
    class           (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), intent(inout)           :: self
    double precision                                                                 , intent(in   )           :: mass
    double precision                                                                 , intent(in   ), optional :: time       , expansionFactor
    logical                                                                          , intent(in   ), optional :: collapsing
    logical                                                                                                    :: remakeTable
    double precision                                                                                           :: time_
    !$GLC attributes unused :: mass

    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Check if we need to recompute our table.
    if (self%turnaroundInitialized) then
       remakeTable=(time_ < self%turnaroundTimeMinimum .or. time_ > self%turnaroundTimeMaximum)
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call self%sphericalCollapseSolver_%radiusTurnaround(time_,self%tableStore,self%turnaround)
       self%turnaroundInitialized=.true.
       self%turnaroundTimeMinimum=self%turnaround%x(+1)
       self%turnaroundTimeMaximum=self%turnaround%x(-1)
    end if
    ! Interpolate to get the ratio.
    sphericalCollapseClsnlssMttrCsmlgclCnstntTrnrndVrlRd=self%turnaround%interpolate(time_)
    return
  end function sphericalCollapseClsnlssMttrCsmlgclCnstntTrnrndVrlRd
