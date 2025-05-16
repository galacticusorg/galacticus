!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Implements an environmental critical overdensity class.
!!}

  !![
  <criticalOverdensity name="criticalOverdensityEnvironmental">
   <description>The critical overdensity is given by some other critical overdensity class multiplied some environment-dependent factor.</description>
  </criticalOverdensity>
  !!]
  type, extends(criticalOverdensityClass) :: criticalOverdensityEnvironmental
     !!{
     A critical overdensity class in which critical overdensity is given by some other critical overdensity class multiplied some environment-dependent factor.
     !!}
     private
     class           (criticalOverdensityClass), pointer :: criticalOverdensity_ => null()
     class           (haloEnvironmentClass    ), pointer :: haloEnvironment_     => null()
     double precision                                    :: a                             , massEnvironment
    contains
     final     ::                    environmentalDestructor
     procedure :: value           => environmentalValue
     procedure :: gradientTime    => environmentalGradientTime
     procedure :: gradientMass    => environmentalGradientMass
     procedure :: isMassDependent => environmentalIsMassDependent
     procedure :: isNodeDependent => environmentalIsNodeDependent
  end type criticalOverdensityEnvironmental

  interface criticalOverdensityEnvironmental
     !!{
     Constructors for the \refClass{criticalOverdensityEnvironmental} critical overdensity class.
     !!}
     module procedure environmentalConstructorParameters
     module procedure environmentalConstructorInternal
  end interface criticalOverdensityEnvironmental

contains

  function environmentalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{criticalOverdensityEnvironmental} critical overdensity class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (criticalOverdensityEnvironmental)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (criticalOverdensityClass        ), pointer       :: criticalOverdensity_
    class           (haloEnvironmentClass            ), pointer       :: haloEnvironment_
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass   ), pointer       :: cosmologicalMassVariance_
    class           (linearGrowthClass               ), pointer       :: linearGrowth_
    double precision                                                  :: a

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>a</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>Parameter controlling environmental dependence of critical overdensity.</description>
    </inputParameter>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="haloEnvironment"          name="haloEnvironment_"          source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=criticalOverdensityEnvironmental(a,criticalOverdensity_,haloEnvironment_,cosmologyFunctions_,linearGrowth_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="haloEnvironment_"         />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function environmentalConstructorParameters

  function environmentalConstructorInternal(a,criticalOverdensity_,haloEnvironment_,cosmologyFunctions_,linearGrowth_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the \refClass{criticalOverdensityEnvironmental} critical overdensity class.
    !!}
    implicit none
    type            (criticalOverdensityEnvironmental)                        :: self
    class           (criticalOverdensityClass        ), target, intent(in   ) :: criticalOverdensity_
    class           (haloEnvironmentClass            ), target, intent(in   ) :: haloEnvironment_
    class           (cosmologyFunctionsClass         ), target, intent(in   ) :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass   ), target, intent(in   ) :: cosmologicalMassVariance_
    class           (linearGrowthClass               ), target, intent(in   ) :: linearGrowth_
    double precision                                          , intent(in   ) :: a
    !![
    <constructorAssign variables="a, *criticalOverdensity_, *haloEnvironment_, *cosmologyFunctions_, *linearGrowth_, *cosmologicalMassVariance_"/>
    !!]

    self%massEnvironment=self%haloEnvironment_%environmentMass()
    return
  end function environmentalConstructorInternal

  subroutine environmentalDestructor(self)
    !!{
    Destructor for the \refClass{criticalOverdensityEnvironmental} critical overdensity class.
    !!}
    implicit none
    type(criticalOverdensityEnvironmental), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%haloEnvironment_"         />
    <objectDestructor name="self%linearGrowth_"            />
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine environmentalDestructor

  double precision function environmentalValue(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the critical overdensity for collapse at the given time and mass.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (criticalOverdensityEnvironmental), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: time      , expansionFactor
    logical                                           , intent(in   ), optional :: collapsing
    double precision                                  , intent(in   ), optional :: mass
    type            (treeNode                        ), intent(inout), optional :: node
    class           (nodeComponentBasic              ), pointer                 :: basic

    if (.not.present(node)) call Error_Report('"node" must be provided to give access to environment'//{introspection:location})
    ! Get the critical overdensity at zero environmental overdensity and scale by some power of the linear growth factor.
    environmentalValue   =  +  self%criticalOverdensity_         %value            (time,expansionFactor,collapsing,mass,node                  )
    basic                =>    node%hostTree            %nodeBase%basic            (                                                           )
    if (basic%mass() < self%massEnvironment) then
       environmentalValue=  +environmentalValue                                                                                                  &
            &               *  self%linearGrowth_                %value            (time,expansionFactor,collapsing                            ) &
            &               **(                                                                                                                  &
            &                  +self%a                                                                                                           &
            &                  *self%haloEnvironment_            %overdensityLinear(                                     node,presentDay=.true.) &
            &                 )
    end if
    return
  end function environmentalValue

  double precision function environmentalGradientTime(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to time of critical overdensity at the given time and mass.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use :: Error           , only : Error_Report
    implicit none
    class           (criticalOverdensityEnvironmental), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: time      , expansionFactor
    logical                                           , intent(in   ), optional :: collapsing
    double precision                                  , intent(in   ), optional :: mass
    type            (treeNode                        ), intent(inout), optional :: node
    class           (nodeComponentBasic              ), pointer                 :: basic
    !![
    <optionalArgument name="expansionFactor" defaultsTo="self%cosmologyFunctions_%expansionFactor(time)" />
    !!]

    if (.not.present(node)) call Error_Report('"node" must be provided to give access to environment'//{introspection:location})
    basic => node%hostTree%nodeBase%basic()
    if (basic%mass() < self%massEnvironment) then
       environmentalGradientTime=+   self%criticalOverdensity_%gradientTime                        (time,expansionFactor ,collapsing,mass,node                  ) &
            &                    *   self%linearGrowth_       %value                               (time,expansionFactor ,collapsing                            ) &
            &                    **(                                                                                                                              &
            &                       +self%a                                                                                                                       &
            &                       *self%haloEnvironment_    %overdensityLinear                   (                                      node,presentDay=.true.) &
            &                      )                                                                                                                              &
            &                    +   self%criticalOverdensity_%value                               (time,expansionFactor ,collapsing,mass,node                  ) &
            &                    *   self%a                                                                                                                       &
            &                    *   self%haloEnvironment_    %overdensityLinear                   (                                      node,presentDay=.true.) &
            &                    *   self%linearGrowth_       %value                               (time,expansionFactor ,collapsing                            ) &
            &                    **(                                                                                                                              &
            &                       +self%a                                                                                                                       &
            &                       *self%haloEnvironment_    %overdensityLinear                   (                                      node,presentDay=.true.) &
            &                      )                                                                                                                              &
            &                    *   self%linearGrowth_       %logarithmicDerivativeExpansionFactor(time,expansionFactor ,collapsing                            ) &
            &                    *   self%cosmologyFunctions_ %expansionRate                       (     expansionFactor_                                       )
    else
       environmentalGradientTime=+   self%criticalOverdensity_%gradientTime                        (time,expansionFactor ,collapsing,mass,node                  )
    end if
    return
  end function environmentalGradientTime

  double precision function environmentalGradientMass(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to mass of critical overdensity at the given time and mass.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use :: Error           , only : Error_Report
    implicit none
    class           (criticalOverdensityEnvironmental), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: time      , expansionFactor
    logical                                           , intent(in   ), optional :: collapsing
    double precision                                  , intent(in   ), optional :: mass
    type            (treeNode                        ), intent(inout), optional :: node
    class           (nodeComponentBasic              ), pointer                 :: basic

    if (.not.present(node)) call Error_Report('"node" must be provided to give access to environment'//{introspection:location})
    basic => node%hostTree%nodeBase%basic()
    if (basic%mass() < self%massEnvironment) then
       environmentalGradientMass=  +self%criticalOverdensity_%gradientMass     (time,expansionFactor,collapsing,mass,node                  ) &
            &                    *  self%linearGrowth_       %value            (time,expansionFactor,collapsing                            ) &
            &                    **(                                                                                                         &
            &                       +self%a                                                                                                  &
            &                       *self%haloEnvironment_   %overdensityLinear(                                     node,presentDay=.true.) &
            &                      )
    else
       environmentalGradientMass=  +self%criticalOverdensity_%gradientMass     (time,expansionFactor,collapsing,mass,node                  )
    end if
    return
  end function environmentalGradientMass

  logical function environmentalIsMassDependent(self)
    !!{
    Return whether the critical overdensity is mass dependent.
    !!}
    implicit none
    class(criticalOverdensityEnvironmental), intent(inout) :: self

    environmentalIsMassDependent=self%criticalOverdensity_%isMassDependent()
    return
  end function environmentalIsMassDependent

  logical function environmentalIsNodeDependent(self)
    !!{
    Return whether the critical overdensity is node dependent.
    !!}
    implicit none
    class(criticalOverdensityEnvironmental), intent(inout) :: self
    !$GLC attributes unused :: self

    environmentalIsNodeDependent=.true.
    return
  end function environmentalIsNodeDependent
