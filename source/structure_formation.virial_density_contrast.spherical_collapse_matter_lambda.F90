!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !# <virialDensityContrast name="virialDensityContrastSphericalCollapseMatterLambda">
  !#  <description>Dark matter halo virial density contrasts based on the spherical collapse in a matter plus cosmological constant universe.</description>
  !# </virialDensityContrast>
  use Tables

  type, extends(virialDensityContrastClass) :: virialDensityContrastSphericalCollapseMatterLambda
     !% A dark matter halo virial density contrast class based on spherical collapse in a matter plus cosmological constant universe.
     private
     logical                                :: tableInitialized
     double precision                       :: tableTimeMinimum, tableTimeMaximum
     class           (table1D), allocatable :: deltaVirial
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
     procedure :: stateStore                  => sphericalCollapseMatterLambdaStateStore
     procedure :: stateRestore                => sphericalCollapseMatterLambdaStateRestore
     procedure :: densityContrast             => sphericalCollapseMatterLambdaDensityContrast
     procedure :: densityContrastRateOfChange => sphericalCollapseMatterLambdaDensityContrastRateOfChange
     procedure :: retabulate                  => sphericalCollapseMatterLambdaRetabulate
  end type virialDensityContrastSphericalCollapseMatterLambda

  interface virialDensityContrastSphericalCollapseMatterLambda
     !% Constructors for the {\tt sphericalCollapseMatterLambda} dark matter halo virial density contrast class.
     module procedure sphericalCollapseMatterLambdaDefaultConstructor
  end interface virialDensityContrastSphericalCollapseMatterLambda

contains

  function sphericalCollapseMatterLambdaDefaultConstructor()
    !% Default constructor for the {\tt sphericalCollapseMatterLambda} dark matter halo virial density contrast class.
    use Input_Parameters
    implicit none
    type (virialDensityContrastSphericalCollapseMatterLambda), target  :: sphericalCollapseMatterLambdaDefaultConstructor
    
    sphericalCollapseMatterLambdaDefaultConstructor%tableInitialized=.false.
    return
  end function sphericalCollapseMatterLambdaDefaultConstructor

  subroutine sphericalCollapseMatterLambdaDestructor(self)
    !% Destructor for the {\tt sphericalCollapseMatterLambda} dark matter halo virial density contrast class.
    implicit none
    type (virialDensityContrastSphericalCollapseMatterLambda), intent(inout) :: self
    
    if (self%tableInitialized) then
       call self%deltaVirial%destroy()
       deallocate(self%deltaVirial)
    end if
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
       call Spherical_Collape_Matter_Lambda_Delta_Virial_Tabulate(time,self%deltaVirial)
       self%tableInitialized=.true.
       self%tableTimeMinimum=self%deltaVirial%x(+1)
       self%tableTimeMaximum=self%deltaVirial%x(-1)
    end if
    return
  end subroutine sphericalCollapseMatterLambdaRetabulate

  double precision function sphericalCollapseMatterLambdaDensityContrast(self,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    use Galacticus_Error
    use Cosmology_Functions
    implicit none
    class           (virialDensityContrastSphericalCollapseMatterLambda), intent(inout)           :: self
    double precision                                                    , intent(in   ), optional :: time               , expansionFactor
    logical                                                             , intent(in   ), optional :: collapsing
    class           (cosmologyFunctionsClass                           ), pointer                 :: cosmologyFunctions_
    logical                                                                                       :: collapsingActual
    double precision                                                                              :: timeActual

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Galacticus_Error_Report('sphericalCollapseMatterLambdaDensityContrast','only one argument can be specified')
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
          ! Get the default cosmology functions object.
          cosmologyFunctions_ => cosmologyFunctions()
          timeActual=cosmologyFunctions_%cosmicTime(expansionFactor,collapsingActual)
       else
          call Galacticus_Error_Report('sphericalCollapseMatterLambdaDensityContrast','at least one argument must be given')
       end if
    end if
    ! Remake the table if necessary.
    call self%retabulate(timeActual)
    ! Interpolate to get the expansion factor.
    sphericalCollapseMatterLambdaDensityContrast=self%deltaVirial%interpolate(timeActual)
    return
  end function sphericalCollapseMatterLambdaDensityContrast

  double precision function sphericalCollapseMatterLambdaDensityContrastRateOfChange(self,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    use Galacticus_Error
    use Cosmology_Functions
    implicit none
    class           (virialDensityContrastSphericalCollapseMatterLambda), intent(inout)           :: self
    double precision                                                    , intent(in   ), optional :: time      , expansionFactor
    logical                                                             , intent(in   ), optional :: collapsing
    class           (cosmologyFunctionsClass                           ), pointer                 :: cosmologyFunctions_
    logical                                                                                       :: collapsingActual
    double precision                                                                              :: timeActual

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Galacticus_Error_Report('sphericalCollapseMatterLambdaDensityContrastRateOfChange','only one argument can be specified')
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
          ! Get the default cosmology functions object.
          cosmologyFunctions_ => cosmologyFunctions()
          timeActual=cosmologyFunctions_%cosmicTime(expansionFactor,collapsingActual)
       else
          call Galacticus_Error_Report('sphericalCollapseMatterLambdaDensityContrastRateOfChange','at least one argument must be given')
       end if
    end if
    ! Remake the table if necessary.
    call self%retabulate(timeActual)
    ! Interpolate to get the expansion factor.
    sphericalCollapseMatterLambdaDensityContrastRateOfChange=self%deltaVirial%interpolateGradient(timeActual)
   return
  end function sphericalCollapseMatterLambdaDensityContrastRateOfChange

  subroutine sphericalCollapseMatterLambdaStateStore(self,stateFile,fgslStateFile)
    !% Reset the tabulation if state is to be retrieved. This will force tables to be rebuilt.
    use FGSL
    implicit none
    class  (virialDensityContrastSphericalCollapseMatterLambda), intent(inout) :: self
    integer                                                    , intent(in   ) :: stateFile
    type   (fgsl_file                                         ), intent(in   ) :: fgslStateFile
    
    write (stateFile) self%tableTimeMinimum,self%tableTimeMaximum
    return
  end subroutine sphericalCollapseMatterLambdaStateStore
  
  subroutine sphericalCollapseMatterLambdaStateRestore(self,stateFile,fgslStateFile)
    !% Reset the tabulation if state is to be retrieved. This will force tables to be rebuilt.
    use FGSL
    implicit none
    class  (virialDensityContrastSphericalCollapseMatterLambda), intent(inout) :: self
    integer                                                    , intent(in   ) :: stateFile
    type   (fgsl_file                                         ), intent(in   ) :: fgslStateFile

    read (stateFile) self%tableTimeMinimum,self%tableTimeMaximum
    self%tableInitialized=.false.
    return
  end subroutine sphericalCollapseMatterLambdaStateRestore
