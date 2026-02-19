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
  An implementation of dark matter halo profile shapes  using the \cite{klypin_multidark_2014} algorithm.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !![
  <enumeration>
   <name>klypin2015Sample</name>
   <description>Enumeration of sample types for the {\normalfont \ttfamily klypin2015} dark matter profile shape parameter class.</description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <entry label="all"    />
   <entry label="relaxed"/>
  </enumeration>
  !!]

  !![
  <darkMatterProfileShape name="darkMatterProfileShapeKlypin2015">
   <description>Dark matter halo shape parameters are computed using the algorithm of \cite{klypin_multidark_2014}.</description>
  </darkMatterProfileShape>
  !!]
  type, extends(darkMatterProfileShapeClass) :: darkMatterProfileShapeKlypin2015
     !!{
     A dark matter halo profile shape parameter class implementing the algorithm of \cite{klypin_multidark_2014}.
     !!}
     private
     class(criticalOverdensityClass       ), pointer :: criticalOverdensity_      => null()
     class(cosmologicalMassVarianceClass  ), pointer :: cosmologicalMassVariance_ => null()
     type (enumerationKlypin2015SampleType)          :: sample
   contains
     final     ::          klypin2015Destructor
     procedure :: shape => klypin2015Shape
  end type darkMatterProfileShapeKlypin2015

  interface darkMatterProfileShapeKlypin2015
     !!{
     Constructors for the \refClass{darkMatterProfileShapeKlypin2015} dark matter halo profile shape parameter class.
     !!}
     module procedure klypin2015ConstructorParameters
     module procedure klypin2015ConstructorInternal
  end interface darkMatterProfileShapeKlypin2015

contains

  function klypin2015ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily klypin2015} dark matter halo profile
    shape class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileShapeKlypin2015)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(criticalOverdensityClass        ), pointer       :: criticalOverdensity_
    class(cosmologicalMassVarianceClass   ), pointer       :: cosmologicalMassVariance_
    type (varying_string                  )                :: sample

    !![
    <inputParameter>
      <name>sample</name>
      <defaultValue>var_str('all')</defaultValue>
      <description>The sample to use for the halo shape parameter algorithm of \cite{klypin_multidark_2014}.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=darkMatterProfileShapeKlypin2015(enumerationKlypin2015SampleEncode(char(sample),includesPrefix=.false.),criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function klypin2015ConstructorParameters

  function klypin2015ConstructorInternal(sample,criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileShapeKlypin2015} dark matter halo profile
    shape class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (darkMatterProfileShapeKlypin2015)                        :: self
    type (enumerationKlypin2015SampleType ), intent(in   )         :: sample
    class(criticalOverdensityClass        ), intent(in   ), target :: criticalOverdensity_
    class(cosmologicalMassVarianceClass   ), intent(in   ), target :: cosmologicalMassVariance_
    !![
    <constructorAssign variables="sample, *criticalOverdensity_, *cosmologicalMassVariance_"/>
    !!]

    if (.not.enumerationKlypin2015SampleIsValid(sample)) call Error_Report('invalid sample'//{introspection:location})
    return
  end function klypin2015ConstructorInternal

  subroutine klypin2015Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileShapeKlypin2015} dark matter halo profile shape class.
    !!}
    implicit none
    type(darkMatterProfileShapeKlypin2015), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine klypin2015Destructor

  double precision function klypin2015Shape(self,node)
    !!{
    Return the Einasto profile shape parameter of the dark matter halo profile of {\normalfont \ttfamily node} using the
    \cite{klypin_multidark_2014} algorithm.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterProfileShapeKlypin2015), intent(inout) :: self
    type            (treeNode                        ), intent(inout) :: node
    class           (nodeComponentBasic              ), pointer       :: basic
    double precision                                                  :: nu

    ! Get the basic component.
    basic => node%basic()
    ! Compute the shape parameter.
    nu     =+self%criticalOverdensity_     %value       (time=basic%time(),mass=basic%mass()) &
         &  /self%cosmologicalMassVariance_%rootVariance(time=basic%time(),mass=basic%mass())
    select case (self%sample%ID)
    case (klypin2015SampleAll    %ID)
       klypin2015Shape=0.115d0+0.0165d0*nu**2
    case (klypin2015SampleRelaxed%ID)
       klypin2015Shape=0.115d0+0.0140d0*nu**2
    case default
       klypin2015Shape=0.0d0
       call Error_Report('unknown sample'//{introspection:location})
    end select
    return
  end function klypin2015Shape
