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
  Implementation of a simple cooling rate class.
  !!}

  !![
  <coolingRate name="coolingRateSimple">
   <description>
    A cooling rate class in which the cooling rate equals the mass of hot gas divided by a fixed timescale. Specifically,
    \begin{equation}
    \dot{M}_\mathrm{cool} = M_\mathrm{hot}/\tau_\mathrm{cool} ,
    \end{equation}
    where $\tau_\mathrm{cool}=${\normalfont \ttfamily [timeScale]}.
   </description>
  </coolingRate>
  !!]
  type, extends(coolingRateClass) :: coolingRateSimple
     !!{
     Implementation of cooling rate class in which the cooling rate equals the mass of hot gas divided by a fixed timescale.
     !!}
     private
     double precision :: timeScale
   contains
     procedure :: rate => simpleRate
  end type coolingRateSimple

  interface coolingRateSimple
     !!{
     Constructors for the simple cooling rate class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface coolingRateSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the simple cooling rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (coolingRateSimple)                :: self
    type            (inputParameters  ), intent(inout) :: parameters
    double precision                                   :: timeScale

    !![
    <inputParameter>
      <name>timeScale</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The timescale (in Gyr) for cooling in the simple cooling rate model.</description>
    </inputParameter>
    !!]
    self=coolingRateSimple(timeScale)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(timeScale) result(self)
    !!{
    Internal constructor for the simple cooling rate class.
    !!}
    implicit none
    type            (coolingRateSimple)                :: self
    double precision                   , intent(in   ) :: timeScale

    !![
    <constructorAssign variables="timeScale"/>
    !!]
    return
  end function simpleConstructorInternal

  double precision function simpleRate(self,node)
    !!{
    Returns the cooling rate (in $M_\odot$ Gyr$^{-1}$) in the hot atmosphere for a model in which this rate is always simple.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, treeNode
    implicit none
    class(coolingRateSimple   ), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    class(nodeComponentHotHalo), pointer       :: hotHalo

    hotHalo    =>  node   %hotHalo  ()
    simpleRate =  +hotHalo%mass     () &
         &        /self   %timeScale
    return
  end function simpleRate

