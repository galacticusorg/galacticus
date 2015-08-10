!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of critical overdensity for collapse based on spherical collapse in a
  !% matter plus cosmological constant universe.

  !# <criticalOverdensity name="criticalOverdensitySphericalCollapseMatterLambda" defaultThreadPrivate="yes">
  !#  <description>Critical overdensity for collapse based on the spherical collapse in a matter plus cosmological constant universe (see, for example, \citealt{percival_cosmological_2005}).</description>
  !# </criticalOverdensity>
  use Tables

  type, extends(criticalOverdensityClass) :: criticalOverdensitySphericalCollapseMatterLambda
     !% A critical overdensity class based on spherical collapse in a matter plus cosmological constant universe.
     private
     logical                                :: tableInitialized
     double precision                       :: tableTimeMinimum   , tableTimeMaximum
     class           (table1D), allocatable :: overdensityCritical
   contains
     !@ <objectMethods>
     !@   <object>criticalOverdensitySphericalCollapseMatterLambda</object>
     !@   <objectMethod>
     !@     <method>retabulate</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ time\argin</arguments>
     !@     <description>Tabulate spherical collapse critical overdensity.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                 sphericalCollapseMatterLambdaDestructor
     procedure :: stateStore   => sphericalCollapseMatterLambdaStateStore
     procedure :: stateRestore => sphericalCollapseMatterLambdaStateRestore
     procedure :: value        => sphericalCollapseMatterLambdaValue
     procedure :: gradientTime => sphericalCollapseMatterLambdaGradientTime
     procedure :: gradientMass => sphericalCollapseMatterLambdaGradientMass
     procedure :: retabulate   => sphericalCollapseMatterLambdaRetabulate
  end type criticalOverdensitySphericalCollapseMatterLambda

  interface criticalOverdensitySphericalCollapseMatterLambda
     !% Constructors for the {\normalfont \ttfamily sphericalCollapseMatterLambda} critical overdensity for collapse class.
     module procedure sphericalCollapseMatterLambdaConstructorParameters
     module procedure sphericalCollapseMatterLambdaConstructorInternal
  end interface criticalOverdensitySphericalCollapseMatterLambda

contains

  function sphericalCollapseMatterLambdaConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily sphericalCollapseMatterLambda} critical overdensity class
    !% which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(criticalOverdensitySphericalCollapseMatterLambda)                :: sphericalCollapseMatterLambdaConstructorParameters
    type(inputParameters                                 ), intent(in   ) :: parameters

    sphericalCollapseMatterLambdaConstructorParameters%tableInitialized=.false.
    return
  end function sphericalCollapseMatterLambdaConstructorParameters

  function sphericalCollapseMatterLambdaConstructorInternal()
    !% Internal constructor for the {\normalfont \ttfamily sphericalCollapseMatterLambda} critical overdensity class.
    implicit none
    type(criticalOverdensitySphericalCollapseMatterLambda) :: sphericalCollapseMatterLambdaConstructorInternal

    sphericalCollapseMatterLambdaConstructorInternal%tableInitialized=.false.
    return
  end function sphericalCollapseMatterLambdaConstructorInternal

  subroutine sphericalCollapseMatterLambdaDestructor(self)
    !% Destructor for the {\normalfont \ttfamily sphericalCollapseMatterLambda} critical overdensity for collapse class.
    implicit none
    type (criticalOverdensitySphericalCollapseMatterLambda), intent(inout) :: self
    
    if (self%tableInitialized) then
       call self%overdensityCritical%destroy()
       deallocate(self%overdensityCritical)
    end if
    return
  end subroutine sphericalCollapseMatterLambdaDestructor

  subroutine sphericalCollapseMatterLambdaRetabulate(self,time)
    !% Recompute the look-up tables for critical overdensity for collapse.
    use Spherical_Collapse_Matter_Lambda
    implicit none
    class           (criticalOverdensitySphericalCollapseMatterLambda), intent(inout) :: self
    double precision                                                  , intent(in   ) :: time
    logical                                                                           :: remakeTable

    ! Check if we need to recompute our table.
    if (self%tableInitialized) then
       remakeTable=(time < self%tableTimeMinimum .or. time > self%tableTimeMaximum)
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Spherical_Collapse_Matter_Lambda_Critical_Overdensity_Tabulate(time,self%overdensityCritical)
       self%tableInitialized=.true.
       self%tableTimeMinimum=self%overdensityCritical%x(+1)
       self%tableTimeMaximum=self%overdensityCritical%x(-1)
    end if
    return
  end subroutine sphericalCollapseMatterLambdaRetabulate

  double precision function sphericalCollapseMatterLambdaValue(self,time,expansionFactor,collapsing,mass)
    !% Return the critical overdensity at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    use Galacticus_Error
    use Cosmology_Functions
    implicit none
    class           (criticalOverdensitySphericalCollapseMatterLambda), intent(inout)           :: self
    double precision                                                  , intent(in   ), optional :: time               , expansionFactor, &
         &                                                                                         mass
    logical                                                           , intent(in   ), optional :: collapsing
    class           (cosmologyFunctionsClass                         ), pointer                 :: cosmologyFunctions_
    logical                                                                                     :: collapsingActual
    double precision                                                                            :: time_

    ! Determine cosmological time.
    cosmologyFunctions_ => cosmologyFunctions()
    call cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Interpolate to get the expansion factor.
    sphericalCollapseMatterLambdaValue=self%overdensityCritical%interpolate(time_)
    return
  end function sphericalCollapseMatterLambdaValue

  double precision function sphericalCollapseMatterLambdaGradientTime(self,time,expansionFactor,collapsing,mass)
    !% Return the time derivative of the critical overdensity at the given epoch, based spherical collapse in a matter plus
    !% cosmological constant universe.
    use Cosmology_Functions
    implicit none
    class           (criticalOverdensitySphericalCollapseMatterLambda), intent(inout)           :: self
    double precision                                                  , intent(in   ), optional :: time               , expansionFactor, &
         &                                                                                         mass
    logical                                                           , intent(in   ), optional :: collapsing
    class           (cosmologyFunctionsClass                         ), pointer                 :: cosmologyFunctions_
    logical                                                                                     :: collapsingActual
    double precision                                                                            :: time_

    ! Determine cosmological time.
    cosmologyFunctions_ => cosmologyFunctions()
    call cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Interpolate to get the expansion factor.
    sphericalCollapseMatterLambdaGradientTime=self%overdensityCritical%interpolateGradient(time_)
    return
  end function sphericalCollapseMatterLambdaGradientTime

  double precision function sphericalCollapseMatterLambdaGradientMass(self,time,expansionFactor,collapsing,mass)
    !% Return the gradient with respect to mass of critical overdensity at the given time and mass.
    use Linear_Growth
    use Cosmology_Functions
    implicit none
    class           (criticalOverdensitySphericalCollapseMatterLambda), intent(inout)           :: self
    double precision                                                  , intent(in   ), optional :: time      , expansionFactor
    logical                                                           , intent(in   ), optional :: collapsing
    double precision                                                  , intent(in   ), optional :: mass

    sphericalCollapseMatterLambdaGradientMass=0.0d0
    return
  end function sphericalCollapseMatterLambdaGradientMass

  subroutine sphericalCollapseMatterLambdaStateStore(self,stateFile,fgslStateFile)
    !% Reset the tabulation if state is to be retrieved. This will force tables to be rebuilt.
    use FGSL
    implicit none
    class  (criticalOverdensitySphericalCollapseMatterLambda), intent(inout) :: self
    integer                                                  , intent(in   ) :: stateFile
    type   (fgsl_file                                       ), intent(in   ) :: fgslStateFile
    
    write (stateFile) self%tableTimeMinimum,self%tableTimeMaximum
    return
  end subroutine sphericalCollapseMatterLambdaStateStore
  
  subroutine sphericalCollapseMatterLambdaStateRestore(self,stateFile,fgslStateFile)
    !% Reset the tabulation if state is to be retrieved. This will force tables to be rebuilt.
    use FGSL
    implicit none
    class  (criticalOverdensitySphericalCollapseMatterLambda), intent(inout) :: self
    integer                                                  , intent(in   ) :: stateFile
    type   (fgsl_file                                       ), intent(in   ) :: fgslStateFile

    read (stateFile) self%tableTimeMinimum,self%tableTimeMaximum
    self%tableInitialized=.false.
    call self%retabulate(sqrt(self%tableTimeMinimum*self%tableTimeMaximum))
    return
  end subroutine sphericalCollapseMatterLambdaStateRestore
