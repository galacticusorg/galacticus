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
  An implementation of dark matter halo profile shapes  using the \cite{gao_redshift_2008} algorithm.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !![
  <darkMatterProfileShape name="darkMatterProfileShapeGao2008">
   <description>
    A dark matter profile shape class in which the shape parameter for Einasto density profiles\index{Einasto
    profile}\index{density profile!Einasto} is computed using a fitting function from \cite{gao_redshift_2008}:
    \begin{equation}
    \alpha = \left\{ \begin{array}{ll} 0.155 + 0.0095\nu^2 &amp; \hbox{ if } \nu &lt; 3.907 \\ 0.3 &amp; \hbox{ if } \nu \ge 3.907,
    \end{array} \right.
    \end{equation}
    where $\nu=\delta_\mathrm{c}(t)/\sigma(M)$ is the peak height of the halo. The truncation at $\alpha = 0.3$ is included
    since \cite{gao_redshift_2008}'s fits do not probe this region and extremely large values of $\alpha$ are numerically
    troublesome.
   </description>
  </darkMatterProfileShape>
  !!]
  type, extends(darkMatterProfileShapeClass) :: darkMatterProfileShapeGao2008
     !!{
     A dark matter halo profile shape parameter class implementing the algorithm of \cite{gao_redshift_2008}.
     !!}
     private
     class(criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
   contains
     final     ::          gao2008Destructor
     procedure :: shape => gao2008Shape
  end type darkMatterProfileShapeGao2008

  interface darkMatterProfileShapeGao2008
     !!{
     Constructors for the \refClass{darkMatterProfileShapeGao2008} dark matter halo profile shape class.
     !!}
     module procedure gao2008ConstructorParameters
     module procedure gao2008ConstructorInternal
  end interface darkMatterProfileShapeGao2008

contains

  function gao2008ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily gao2008} dark matter halo profile
    shape class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileShapeGao2008)                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(criticalOverdensityClass     ), pointer       :: criticalOverdensity_
    class(cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_

    !![
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=darkMatterProfileShapeGao2008(criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function gao2008ConstructorParameters

  function gao2008ConstructorInternal(criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileShapeGao2008} dark matter halo profile
    shape class.
    !!}
    implicit none
    type (darkMatterProfileShapeGao2008)                        :: self
    class(criticalOverdensityClass     ), intent(in   ), target :: criticalOverdensity_
    class(cosmologicalMassVarianceClass), intent(in   ), target :: cosmologicalMassVariance_
    !![
    <constructorAssign variables=" *criticalOverdensity_, *cosmologicalMassVariance_"/>
    !!]

    return
  end function gao2008ConstructorInternal

  subroutine gao2008Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileShapeGao2008} dark matter halo profile shape class.
    !!}
    implicit none
    type(darkMatterProfileShapeGao2008), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine gao2008Destructor

  double precision function gao2008Shape(self,node)
    !!{
    Return the Einasto profile shape parameter of the dark matter halo profile of {\normalfont \ttfamily node} using the
    \cite{gao_redshift_2008} algorithm. More specifically, the parameter is given by:
    \begin{equation}
    \alpha = \left\{ \begin{array}{ll} 0.155 + 0.0095\nu^2 & \hbox{ if } \nu < 3.907 \\ 0.3 & \hbox{ if } \nu \ge 3.907, \end{array} \right.
    \end{equation}
    where $\nu=\delta_\mathrm{c}(t)/\sigma(M)$ is the peak height of the halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterProfileShapeGao2008), intent(inout)  :: self
    type            (treeNode                     ), intent(inout), :: node
    double precision                               , parameter      :: nuMaximum=3.907d0
    class           (nodeComponentBasic           ), pointer        :: basic
    double precision                                                :: nu

    ! Get the basic component.
    basic => node%basic()
    ! Compute the shape parameter.
    nu     =+self%criticalOverdensity_     %value       (time=basic%time(),mass=basic%mass()) &
         &  /self%cosmologicalMassVariance_%rootVariance(time=basic%time(),mass=basic%mass())
    if (nu < nuMaximum) then
       gao2008Shape=0.155d0+0.0095d0*nu**2
    else
       gao2008Shape=0.300d0
    end if
    return
  end function gao2008Shape
