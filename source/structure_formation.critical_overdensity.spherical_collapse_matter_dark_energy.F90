!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

  !% An implementation of critical overdensity for collapse based on spherical collapse in a
  !% matter plus dark energy universe.

  !# <criticalOverdensity name="criticalOverdensitySphericalCollapseMatterDE" defaultThreadPrivate="yes">
  !#  <description>Critical overdensity for collapse based on the spherical collapse in a matter plus dark energy universe.</description>
  !# </criticalOverdensity>
  type, extends(criticalOverdensitySphericalCollapseMatterLambda) :: criticalOverdensitySphericalCollapseMatterDE
     !% A dark matter halo virial density contrast class based on spherical collapse in a matter plus dark energy universe.
     private
   contains
     final     ::                 sphericalCollapseMatterDEDestructor
     procedure :: retabulate   => sphericalCollapseMatterDERetabulate
  end type criticalOverdensitySphericalCollapseMatterDE

  interface criticalOverdensitySphericalCollapseMatterDE
     !% Constructors for the {\normalfont \ttfamily sphericalCollapseMatterDE} critical overdensity for collapse class.
     module procedure sphericalCollapseMatterDEConstructorParameters
     module procedure sphericalCollapseMatterDEConstructorInternal
  end interface criticalOverdensitySphericalCollapseMatterDE

contains

  function sphericalCollapseMatterDEConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily sphericalCollapseMatterDE} critical overdensity class
    !% which takes a parameter set as input.
    use Input_Parameters2
    use Dark_Matter_Particles
    use Galacticus_Error
    implicit none
    type (criticalOverdensitySphericalCollapseMatterDE)                :: sphericalCollapseMatterDEConstructorParameters
    type (inputParameters                             ), intent(inout) :: parameters
    class(darkMatterParticleClass                     ), pointer       :: darkMatterParticle_

    !# <objectBuilder class="linearGrowth"             name="sphericalCollapseMatterDEConstructorParameters%linearGrowth_"             source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"       name="sphericalCollapseMatterDEConstructorParameters%cosmologyFunctions_"       source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="sphericalCollapseMatterDEConstructorParameters%cosmologicalMassVariance_" source="parameters"/>
    !# <objectBuilder class="darkMatterParticle"       name="darkMatterParticle_"                                                      source="parameters"/>
    sphericalCollapseMatterDEConstructorParameters%tableInitialized=.false.
    select type (darkMatterParticle_)
    class is (darkMatterParticleCDM)
       ! Cold dark matter particle - this is as expected.
    class default
       call Galacticus_Error_Report('sphericalCollapseMatterDEConstructorParameters','critical overdensity expects a cold dark matter particle')
    end select
    return
  end function sphericalCollapseMatterDEConstructorParameters

  function sphericalCollapseMatterDEConstructorInternal(linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_)
    !% Internal constructor for the {\normalfont \ttfamily sphericalCollapseMatterDE} critical overdensity class.
    use Dark_Matter_Particles
    use Galacticus_Error
    implicit none
    type (criticalOverdensitySphericalCollapseMatterDE)                        :: sphericalCollapseMatterDEConstructorInternal
    class(cosmologyFunctionsClass                     ), target, intent(in   ) :: cosmologyFunctions_    
    class(linearGrowthClass                           ), target, intent(in   ) :: linearGrowth_    
    class(cosmologicalMassVarianceClass               ), target, intent(in   ) :: cosmologicalMassVariance_
    class(darkMatterParticleClass                     )        , intent(in   ) :: darkMatterParticle_

    sphericalCollapseMatterDEConstructorInternal%tableInitialized          =  .false.
    sphericalCollapseMatterDEConstructorInternal%cosmologyFunctions_       => cosmologyFunctions_
    sphericalCollapseMatterDEConstructorInternal%linearGrowth_             => linearGrowth_
    sphericalCollapseMatterDEConstructorInternal%cosmologicalMassVariance_ => cosmologicalMassVariance_
    ! Require that the dark matter be cold dark matter.
    select type (darkMatterParticle_)
    class is (darkMatterParticleCDM)
       ! Cold dark matter particle - this is as expected.
    class default
       call Galacticus_Error_Report('sphericalCollapseMatterDEConstructorInternal','critical overdensity expects a cold dark matter particle')
    end select
    return
  end function sphericalCollapseMatterDEConstructorInternal

  subroutine sphericalCollapseMatterDEDestructor(self)
    !% Destructor for the {\normalfont \ttfamily sphericalCollapseMatterDE} critical overdensity for collapse class.
    implicit none
    type (criticalOverdensitySphericalCollapseMatterDE), intent(inout) :: self
    
    !# <objectDestructor name="self%cosmologyFunctions_"/>
    !# <objectDestructor name="self%linearGrowth_"      />
    if (self%tableInitialized) then
       call self%overdensityCritical%destroy()
       deallocate(self%overdensityCritical)
    end if
    return
  end subroutine sphericalCollapseMatterDEDestructor

  subroutine sphericalCollapseMatterDERetabulate(self,time)
    !% Recompute the look-up tables for critical overdensity for collapse.
    use Spherical_Collapse_Matter_Dark_Energy
    implicit none
    class           (criticalOverdensitySphericalCollapseMatterDE), intent(inout) :: self
    double precision                                              , intent(in   ) :: time
    logical                                                                       :: remakeTable

    ! Check if we need to recompute our table.
    if (self%tableInitialized) then
       remakeTable=(time < self%tableTimeMinimum .or. time > self%tableTimeMaximum)
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Spherical_Collapse_Dark_Energy_Critical_Overdensity_Tabulate(time,self%overdensityCritical,self%linearGrowth_,self%cosmologyFunctions_)
       self%tableInitialized=.true.
       self%tableTimeMinimum=self%overdensityCritical%x(+1)
       self%tableTimeMaximum=self%overdensityCritical%x(-1)
    end if
    return
  end subroutine sphericalCollapseMatterDERetabulate
