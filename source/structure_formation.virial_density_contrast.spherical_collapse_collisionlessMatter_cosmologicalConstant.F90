!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% An implementation of dark matter halo virial density contrasts based on spherical collapse in a matter plus cosmological constant universe.

  use Tables                    , only : table1D
  use Cosmology_Functions       , only : cosmologyFunctionsClass                         , cosmologyFunctions
  use Spherical_Collapse_Solvers, only : sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt

  !# <virialDensityContrast name="virialDensityContrastSphericalCollapseCllsnlssMttrCsmlgclCnstnt">
  !#  <description>Dark matter halo virial density contrasts based on the spherical collapse in a matter plus cosmological constant universe.</description>
  !#  <deepCopy>
  !#   <functionClass variables="sphericalCollapseSolver_"/>
  !#  </deepCopy>
  !#  <stateStorable>
  !#   <functionClass variables="sphericalCollapseSolver_"/>
  !#  </stateStorable>
  !# </virialDensityContrast>
  type, extends(virialDensityContrastClass) :: virialDensityContrastSphericalCollapseCllsnlssMttrCsmlgclCnstnt
     !% A dark matter halo virial density contrast class based on spherical collapse in a matter plus cosmological constant universe.
     private
     logical                                                                         :: tableInitialized        =  .false., turnaroundInitialized=.false.
     double precision                                                                :: tableTimeMinimum                  , tableTimeMaximum             , &
          &                                                                             turnaroundTimeMinimum             , turnaroundTimeMaximum
     logical                                                                         :: tableStore
     class           (table1D                                         ), allocatable :: deltaVirial                       , turnaround
     class           (cosmologyFunctionsClass                         ), pointer     :: cosmologyFunctions_      => null()
     class           (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt), pointer     :: sphericalCollapseSolver_ => null()
   contains
     !@ <objectMethods>
     !@   <object>virialDensityContrastSphericalCollapseCllsnlssMttrCsmlgclCnstnt</object>
     !@   <objectMethod>
     !@     <method>retabulate</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ time\argin</arguments>
     !@     <description>Tabulate spherical collapse virial density contrast.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                sphericalCollapseCllsnlssMttrCsmlgclCnstntDestructor
     procedure :: densityContrast             => sphericalCollapseCllsnlssMttrCsmlgclCnstntDensityContrast
     procedure :: densityContrastRateOfChange => sphericalCollapseCllsnlssMttrCsmlgclCnstntDensityContrastRtChng
     procedure :: turnAroundOverVirialRadii   => sphericalCollapseCllsnlssMttrCsmlgclCnstntTrnrndVrlRd
     procedure :: retabulate                  => sphericalCollapseCllsnlssMttrCsmlgclCnstntRetabulate
  end type virialDensityContrastSphericalCollapseCllsnlssMttrCsmlgclCnstnt

  interface virialDensityContrastSphericalCollapseCllsnlssMttrCsmlgclCnstnt
     !% Constructors for the {\normalfont \ttfamily sphericalCollapseCllsnlssMttrCsmlgclCnstnt} dark matter halo virial density contrast class.
     module procedure sphericalCollapseCllsnlssMttrCsmlgclCnstntConstructorParameters
     module procedure sphericalCollapseCllsnlssMttrCsmlgclCnstntConstructorInternal
  end interface virialDensityContrastSphericalCollapseCllsnlssMttrCsmlgclCnstnt

contains

  function sphericalCollapseCllsnlssMttrCsmlgclCnstntConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily sphericalCollapseCllsnlssMttrCsmlgclCnstnt} dark matter halo virial density contrast class that takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (virialDensityContrastSphericalCollapseCllsnlssMttrCsmlgclCnstnt)                :: self
    type   (inputParameters                                                ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass                                        ), pointer       :: cosmologyFunctions_
    logical                                                                                 :: tableStore

    !# <inputParameter>
    !#   <name>tableStore</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>If true, store/restore the tabulated solution to/from file when possible.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    self=virialDensityContrastSphericalCollapseCllsnlssMttrCsmlgclCnstnt(tableStore,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"/>
    return
  end function sphericalCollapseCllsnlssMttrCsmlgclCnstntConstructorParameters

  function sphericalCollapseCllsnlssMttrCsmlgclCnstntConstructorInternal(tableStore,cosmologyFunctions_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily sphericalCollapseCllsnlssMttrCsmlgclCnstnt} dark matter halo virial density contrast class.
    implicit none
    type   (virialDensityContrastSphericalCollapseCllsnlssMttrCsmlgclCnstnt)                        :: self
    class  (cosmologyFunctionsClass                                        ), intent(in   ), target :: cosmologyFunctions_
    logical                                                                 , intent(in   )         :: tableStore
    !# <constructorAssign variables="tableStore, *cosmologyFunctions_"/>

    self%tableInitialized     =.false.
    self%turnaroundInitialized=.false.
    allocate(sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt :: self%sphericalCollapseSolver_)
    select type (sphericalCollapseSolver_ => self%sphericalCollapseSolver_)
    type is (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt)
       !# <referenceConstruct isResult="yes" object="sphericalCollapseSolver_" constructor="sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt(self%cosmologyFunctions_)"/>
    end select
    return
  end function sphericalCollapseCllsnlssMttrCsmlgclCnstntConstructorInternal

  subroutine sphericalCollapseCllsnlssMttrCsmlgclCnstntDestructor(self)
    !% Destructor for the {\normalfont \ttfamily sphericalCollapseCllsnlssMttrCsmlgclCnstnt} dark matter halo virial density contrast class.
    implicit none
    type (virialDensityContrastSphericalCollapseCllsnlssMttrCsmlgclCnstnt), intent(inout) :: self

    if (self%tableInitialized) then
       call self%deltaVirial%destroy()
       deallocate(self%deltaVirial)
    end if
    if (self%turnaroundInitialized) then
       call self%turnaround%destroy()
       deallocate(self%turnaround)
    end if
    !# <objectDestructor name="self%cosmologyFunctions_"     />
    !# <objectDestructor name="self%sphericalCollapseSolver_"/>
    return
  end subroutine sphericalCollapseCllsnlssMttrCsmlgclCnstntDestructor

  subroutine sphericalCollapseCllsnlssMttrCsmlgclCnstntRetabulate(self,time)
    !% Recompute the look-up tables for virial density contrast.
    implicit none
    class           (virialDensityContrastSphericalCollapseCllsnlssMttrCsmlgclCnstnt), intent(inout) :: self
    double precision                                                                 , intent(in   ) :: time
    logical                                                                                          :: remakeTable

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
  end subroutine sphericalCollapseCllsnlssMttrCsmlgclCnstntRetabulate

  double precision function sphericalCollapseCllsnlssMttrCsmlgclCnstntDensityContrast(self,mass,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (virialDensityContrastSphericalCollapseCllsnlssMttrCsmlgclCnstnt), intent(inout)           :: self
    double precision                                                                 , intent(in   )           :: mass
    double precision                                                                 , intent(in   ), optional :: time            , expansionFactor
    logical                                                                          , intent(in   ), optional :: collapsing
    logical                                                                                                    :: collapsingActual
    double precision                                                                                           :: timeActual
    !GCC$ attributes unused :: mass

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Galacticus_Error_Report('only one argument can be specified'//{introspection:location})
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
          call Galacticus_Error_Report('at least one argument must be given'//{introspection:location})
       end if
    end if
    ! Remake the table if necessary.
    call self%retabulate(timeActual)
    ! Interpolate to get the expansion factor.
    sphericalCollapseCllsnlssMttrCsmlgclCnstntDensityContrast=self%deltaVirial%interpolate(timeActual)
    return
  end function sphericalCollapseCllsnlssMttrCsmlgclCnstntDensityContrast

  double precision function sphericalCollapseCllsnlssMttrCsmlgclCnstntDensityContrastRtChng(self,mass,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (virialDensityContrastSphericalCollapseCllsnlssMttrCsmlgclCnstnt), intent(inout)           :: self
    double precision                                                                 , intent(in   )           :: mass
    double precision                                                                 , intent(in   ), optional :: time            , expansionFactor
    logical                                                                          , intent(in   ), optional :: collapsing
    logical                                                                                                    :: collapsingActual
    double precision                                                                                           :: timeActual
    !GCC$ attributes unused :: mass

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Galacticus_Error_Report('only one argument can be specified'//{introspection:location})
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
          call Galacticus_Error_Report('at least one argument must be given'//{introspection:location})
       end if
    end if
    ! Remake the table if necessary.
    call self%retabulate(timeActual)
    ! Interpolate to get the expansion factor.
    sphericalCollapseCllsnlssMttrCsmlgclCnstntDensityContrastRtChng=self%deltaVirial%interpolateGradient(timeActual)
   return
  end function sphericalCollapseCllsnlssMttrCsmlgclCnstntDensityContrastRtChng

  double precision function sphericalCollapseCllsnlssMttrCsmlgclCnstntTrnrndVrlRd(self,mass,time,expansionFactor,collapsing)
    !% Return the ratio of turnaround and virial radii at the given epoch, based spherical collapse in a matter plus cosmological
    !% constant universe.
    implicit none
    class           (virialDensityContrastSphericalCollapseCllsnlssMttrCsmlgclCnstnt), intent(inout)           :: self
    double precision                                                                 , intent(in   )           :: mass
    double precision                                                                 , intent(in   ), optional :: time       , expansionFactor
    logical                                                                          , intent(in   ), optional :: collapsing
    logical                                                                                                    :: remakeTable
    double precision                                                                                           :: time_
    !GCC$ attributes unused :: mass

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
    sphericalCollapseCllsnlssMttrCsmlgclCnstntTrnrndVrlRd=self%turnaround%interpolate(time_)
    return
  end function sphericalCollapseCllsnlssMttrCsmlgclCnstntTrnrndVrlRd
