!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  use Tables
  use Cosmology_Functions

  !# <virialDensityContrast name="virialDensityContrastSphericalCollapseMatterLambda">
  !#  <description>Dark matter halo virial density contrasts based on the spherical collapse in a matter plus cosmological constant universe.</description>
  !# </virialDensityContrast>
  type, extends(virialDensityContrastClass) :: virialDensityContrastSphericalCollapseMatterLambda
     !% A dark matter halo virial density contrast class based on spherical collapse in a matter plus cosmological constant universe.
     private
     logical                                                :: tableInitialized    =  .false.
     double precision                                       :: tableTimeMinimum              , tableTimeMaximum
     class           (table1D                ), allocatable :: deltaVirial
     class           (cosmologyFunctionsClass), pointer     :: cosmologyFunctions_ => null()
   contains
     !@ <objectMethods>
     !@   <object>virialDensityContrastSphericalCollapseMatterLambda</object>
     !@   <objectMethod>
     !@     <method>retabulate</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ time\argin</arguments>
     !@     <description>Tabulate spherical collapse virial density contrast.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                sphericalCollapseMatterLambdaDestructor
     procedure :: densityContrast             => sphericalCollapseMatterLambdaDensityContrast
     procedure :: densityContrastRateOfChange => sphericalCollapseMatterLambdaDensityContrastRateOfChange
     procedure :: turnAroundOverVirialRadii   => sphericalCollapseMatterLambdaTurnAroundOverVirialRadii
     procedure :: retabulate                  => sphericalCollapseMatterLambdaRetabulate
  end type virialDensityContrastSphericalCollapseMatterLambda

  interface virialDensityContrastSphericalCollapseMatterLambda
     !% Constructors for the {\normalfont \ttfamily sphericalCollapseMatterLambda} dark matter halo virial density contrast class.
     module procedure sphericalCollapseMatterLambdaConstructorParameters
     module procedure sphericalCollapseMatterLambdaConstructorInternal
  end interface virialDensityContrastSphericalCollapseMatterLambda

contains

  function sphericalCollapseMatterLambdaConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily sphericalCollapseMatterLambda} dark matter halo virial density contrast class that takes a parameter set as input.
    use Input_Parameters
    implicit none
    type (virialDensityContrastSphericalCollapseMatterLambda)                :: self
    type (inputParameters                                   ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                           ), pointer       :: cosmologyFunctions_
    
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    self=virialDensityContrastSphericalCollapseMatterLambda(cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function sphericalCollapseMatterLambdaConstructorParameters

  function sphericalCollapseMatterLambdaConstructorInternal(cosmologyFunctions_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily sphericalCollapseMatterLambda} dark matter halo virial density contrast class.
    implicit none
    type (virialDensityContrastSphericalCollapseMatterLambda)                       :: self
    class(cosmologyFunctionsClass                          ), intent(in   ), target :: cosmologyFunctions_
    !# <constructorAssign variables="*cosmologyFunctions_"/>

    self%tableInitialized=.false.
    return
  end function sphericalCollapseMatterLambdaConstructorInternal

  subroutine sphericalCollapseMatterLambdaDestructor(self)
    !% Destructor for the {\normalfont \ttfamily sphericalCollapseMatterLambda} dark matter halo virial density contrast class.
    implicit none
    type (virialDensityContrastSphericalCollapseMatterLambda), intent(inout) :: self
    
    if (self%tableInitialized) then
       call self%deltaVirial%destroy()
       deallocate(self%deltaVirial)
    end if
    !# <objectDestructor name="self%cosmologyFunctions_" />
    return
  end subroutine sphericalCollapseMatterLambdaDestructor

  subroutine sphericalCollapseMatterLambdaRetabulate(self,time)
    !% Recompute the look-up tables for virial density contrast.
    use Spherical_Collapse_Matter_Lambda
    implicit none
    class           (virialDensityContrastSphericalCollapseMatterLambda), intent(inout) :: self
    double precision                                                    , intent(in   ) :: time
    logical                                                                             :: remakeTable

    ! Check if we need to recompute our table.
    if (self%tableInitialized) then
       remakeTable=(time < self%tableTimeMinimum .or. time > self%tableTimeMaximum)
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Spherical_Collape_Matter_Lambda_Delta_Virial_Tabulate(time,self%deltaVirial,self%cosmologyFunctions_)
       self%tableInitialized=.true.
       self%tableTimeMinimum=self%deltaVirial%x(+1)
       self%tableTimeMaximum=self%deltaVirial%x(-1)
    end if
    return
  end subroutine sphericalCollapseMatterLambdaRetabulate

  double precision function sphericalCollapseMatterLambdaDensityContrast(self,mass,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    use Galacticus_Error
    implicit none
    class           (virialDensityContrastSphericalCollapseMatterLambda), intent(inout)           :: self
    double precision                                                    , intent(in   )           :: mass
    double precision                                                    , intent(in   ), optional :: time            , expansionFactor
    logical                                                             , intent(in   ), optional :: collapsing
    logical                                                                                       :: collapsingActual
    double precision                                                                              :: timeActual
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
    sphericalCollapseMatterLambdaDensityContrast=self%deltaVirial%interpolate(timeActual)
    return
  end function sphericalCollapseMatterLambdaDensityContrast

  double precision function sphericalCollapseMatterLambdaDensityContrastRateOfChange(self,mass,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    use Galacticus_Error
    implicit none
    class           (virialDensityContrastSphericalCollapseMatterLambda), intent(inout)           :: self
    double precision                                                    , intent(in   )           :: mass
    double precision                                                    , intent(in   ), optional :: time            , expansionFactor
    logical                                                             , intent(in   ), optional :: collapsing
    logical                                                                                       :: collapsingActual
    double precision                                                                              :: timeActual
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
    sphericalCollapseMatterLambdaDensityContrastRateOfChange=self%deltaVirial%interpolateGradient(timeActual)
   return
  end function sphericalCollapseMatterLambdaDensityContrastRateOfChange

  double precision function sphericalCollapseMatterLambdaTurnAroundOverVirialRadii(self,time,expansionFactor,collapsing)
    !% Return the ratio of turnaround and virial radii at the given epoch, based spherical collapse in a matter plus cosmological
    !% constant universe.
    implicit none
    class           (virialDensityContrastSphericalCollapseMatterLambda), intent(inout)           :: self
    double precision                                                    , intent(in   ), optional :: time      , expansionFactor
    logical                                                             , intent(in   ), optional :: collapsing
    !GCC$ attributes unused :: self, time, expansionFactor, collapsing
    
    ! In simple cosmological constant dark energy universes, this ratio is always precisely 2 (e.g. Percival 2005;
    ! http://adsabs.harvard.edu/abs/2005A%26A...443..819P)
    sphericalCollapseMatterLambdaTurnAroundOverVirialRadii=2.0d0
    return
  end function sphericalCollapseMatterLambdaTurnAroundOverVirialRadii
