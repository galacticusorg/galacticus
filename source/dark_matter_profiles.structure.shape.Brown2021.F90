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
  An implementation of dark matter halo profile shapes  using the \cite{brown_towards_2022} algorithm.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !![
  <darkMatterProfileShape name="darkMatterProfileShapeBrown2021">
   <description>
    A dark matter profile shape class in which the shape parameter for Einasto density profiles\index{Einasto
    profile}\index{density profile!Einasto} is computed using a fitting function from \cite[][eqn. 21]{brown_towards_2022}:
    \begin{equation}
    \alpha = \left\{ \begin{array}{ll} 8.52 \times 10^{-4} \nu_\alpha^4 + 0.166 &amp; \hbox{ if } \nu_\alpha &lt; 3.541 \\ 0.3 &amp; \hbox{ if } \nu_\alpha \ge 3.541, \end{array} \right.
    \end{equation}
    where $\nu_\alpha=\delta_\mathrm{c}(t)/\sigma_\alpha(M)$ is the peak height of the halo. The truncation at $\alpha = 0.3$ is included
    since \cite{brown_towards_2022}'s fits do not probe this region and extremely large values of $\alpha$ are numerically
    troublesome.
     
     This implementation accepts any \refClass{cosmologicalMassVarianceClass} object for use in computing for computing
     $\sigma_\mathrm{c}(M)$. \emph{However}, \cite{brown_towards_2022} recommend using $\sigma_\mathrm{c}(M)$ computed using a
     generalized top-hat window function (\refClass{powerSpectrumWindowFunctionTopHatGeneralized}) with $\mu_\mathrm{g}=0.01$.
   </description>
  </darkMatterProfileShape>
  !!]
  type, extends(darkMatterProfileShapeClass) :: darkMatterProfileShapeBrown2021
     !!{
     A dark matter halo profile shape parameter class implementing the algorithm of \cite{brown_towards_2022}.
     !!}
     private
     class(criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
   contains
     final     ::          brown2021Destructor
     procedure :: shape => brown2021Shape
  end type darkMatterProfileShapeBrown2021

  interface darkMatterProfileShapeBrown2021
     !!{
     Constructors for the \refClass{darkMatterProfileShapeBrown2021} dark matter halo profile shape class.
     !!}
     module procedure brown2021ConstructorParameters
     module procedure brown2021ConstructorInternal
  end interface darkMatterProfileShapeBrown2021

contains

  function brown2021ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily brown2021} dark matter halo profile
    shape class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileShapeBrown2021)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(criticalOverdensityClass       ), pointer       :: criticalOverdensity_
    class(cosmologicalMassVarianceClass  ), pointer       :: cosmologicalMassVariance_

    !![
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=darkMatterProfileShapeBrown2021(criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function brown2021ConstructorParameters

  function brown2021ConstructorInternal(criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileShapeBrown2021} dark matter halo profile
    shape class.
    !!}
    implicit none
    type (darkMatterProfileShapeBrown2021)                        :: self
    class(criticalOverdensityClass       ), intent(in   ), target :: criticalOverdensity_
    class(cosmologicalMassVarianceClass  ), intent(in   ), target :: cosmologicalMassVariance_
    !![
    <constructorAssign variables=" *criticalOverdensity_, *cosmologicalMassVariance_"/>
    !!]

    return
  end function brown2021ConstructorInternal

  subroutine brown2021Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileShapeBrown2021} dark matter halo profile shape class.
    !!}
    implicit none
    type(darkMatterProfileShapeBrown2021), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine brown2021Destructor

  double precision function brown2021Shape(self,node)
    !!{
    Return the Einasto profile shape parameter of the dark matter halo profile of {\normalfont \ttfamily node} using the
    \cite[][eqn. 21]{brown_towards_2022} algorithm. Specifically, the parameter is given by:
    \begin{equation}
     \alpha = \left\{ \begin{array}{ll} 8.52 \times 10^{-4} \nu_\alpha^4 + 0.166 & \hbox{ if } \nu_\alpha < 3.541 \\ 0.3 & \hbox{ if } \nu_\alpha \ge 3.541, \end{array} \right.
    \end{equation}
    where $\nu_\alpha=\delta_\mathrm{c}(t)/\sigma_\alpha(M)$ is the peak height of the halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterProfileShapeBrown2021), intent(inout)  :: self
    type            (treeNode                       ), intent(inout), :: node
    double precision                                 , parameter      :: peakHeightMaximum=3.541d0
    class           (nodeComponentBasic             ), pointer        :: basic
    double precision                                                  :: peakHeight

    ! Get the basic component.
    basic => node%basic()
    ! Compute the shape parameter.
    peakHeight=+self%criticalOverdensity_     %value       (time=basic%time(),mass=basic%mass()) &
         &     /self%cosmologicalMassVariance_%rootVariance(time=basic%time(),mass=basic%mass())
    if (peakHeight < peakHeightMaximum) then
       brown2021Shape=+8.520d-4      &
            &         *peakHeight**4 &
            &         +0.166d+0
    else
       brown2021Shape=+0.300d+0
    end if
    return
  end function brown2021Shape
