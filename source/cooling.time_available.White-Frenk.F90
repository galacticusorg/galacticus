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
  Implementation of the \cite{white_galaxy_1991} time available for cooling class.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <coolingTimeAvailable name="coolingTimeAvailableWhiteFrenk1991">
   <description>
    A time available for cooling class which implements the algorithm of \cite{white_galaxy_1991}. The time available for
    cooling is equal to
    \begin{equation}
     t_\mathrm{available} = \exp\left[ f \ln t_\mathrm{Universe} + (1-f)\ln t_\mathrm{dynamical} \right],
    \end{equation}
    where $f=${\normalfont \ttfamily [ageFactor]} is an interpolating factor, $t_\mathrm{Universe}$ is the age of the Universe
    and $t_\mathrm{dynamical}$ is the dynamical time in the halo. The original \cite{white_galaxy_1991} algorithm corresponds
    to $f=1$.
   </description>
  </coolingTimeAvailable>
  !!]
  type, extends(coolingTimeAvailableClass) :: coolingTimeAvailableWhiteFrenk1991
     !!{
     Implementation of a time available for cooling class which implements the algorithm of \cite{white_galaxy_1991}.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: ageFactor
   contains
     final     ::                              whiteFrenk1991Destructor
     procedure :: timeAvailable             => whiteFrenk1991TimeAvailable
     procedure :: timeAvailableIncreaseRate => whiteFrenk1991TimeAvailableIncreaseRate
  end type coolingTimeAvailableWhiteFrenk1991

  interface coolingTimeAvailableWhiteFrenk1991
     !!{
     Constructors for the \cite{white_galaxy_1991} time available for cooling class.
     !!}
     module procedure whiteFrenk1991ConstructorParameters
     module procedure whiteFrenk1991ConstructorInternal
  end interface coolingTimeAvailableWhiteFrenk1991

contains

  function whiteFrenk1991ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{white_galaxy_1991} time available for cooling class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (coolingTimeAvailableWhiteFrenk1991)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass          ), pointer       :: darkMatterHaloScale_
    double precision                                                    :: ageFactor

    !![
    <inputParameter>
     <name>ageFactor</name>
     <source>parameters</source>
     <defaultValue>0.0d0</defaultValue>
     <description>Interpolates (geometrically) between the age of the Universe and the halo dynamical time for the time available for cooling in the \cite{white_galaxy_1991} method.</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=coolingTimeAvailableWhiteFrenk1991(ageFactor,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function whiteFrenk1991ConstructorParameters

  function whiteFrenk1991ConstructorInternal(ageFactor,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \cite{white_galaxy_1991} cooling rate class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (coolingTimeAvailableWhiteFrenk1991)                        :: self
    double precision                                    , intent(in   )         :: ageFactor
    class           (darkMatterHaloScaleClass          ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="ageFactor, *darkMatterHaloScale_"/>
    !!]

    if     (                   &
         &   ageFactor < 0.0d0 &
         &  .or.               &
         &   ageFactor > 1.0d0 &
         & ) call Error_Report('0 ≤ coolingTimeAvailableAgeFactor ≤ 1 is required'//{introspection:location})
    return
  end function whiteFrenk1991ConstructorInternal

  subroutine whiteFrenk1991Destructor(self)
    !!{
    Destructor for the \cite{white_galaxy_1991} cooling rate class.
    !!}
    implicit none
    type(coolingTimeAvailableWhiteFrenk1991), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine whiteFrenk1991Destructor

  double precision function whiteFrenk1991TimeAvailable(self,node)
    !!{
    Returns the time available for cooling (in units of Gyr) in the hot atmosphere for the \cite{white_galaxy_1991} model.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (coolingTimeAvailableWhiteFrenk1991), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    class           (nodeComponentBasic                ), pointer       :: basic
    double precision                                                    :: timeNow, timeDynamical

    ! Return the appropriate time.
    if (self%ageFactor == 1.0d0) then
       ! Time available equals the age of the Universe, which is just the time for this node.
       basic                       => node %basic()
       whiteFrenk1991TimeAvailable =  basic%time ()
    else if (self%ageFactor == 0.0d0) then
       ! Time available equals the halo dynamical time.
       whiteFrenk1991TimeAvailable =  self%darkMatterHaloScale_%timescaleDynamical(node)
    else
       ! Time is interpolated between age of Universe and dynamical time. Do the interpolation.
       basic                       =>  node %basic                                  (    )
       timeDynamical               =   self %darkMatterHaloScale_%timescaleDynamical(node)
       timeNow                     =   basic                     %time              (    )
       whiteFrenk1991TimeAvailable =  +  timeDynamical   &
            &                         *(                 &
            &                           +timeNow         &
            &                           /timeDynamical   &
            &                          )**self%ageFactor
    end if
    return
  end function whiteFrenk1991TimeAvailable

  double precision function whiteFrenk1991TimeAvailableIncreaseRate(self,node)
    !!{
    Compute the rate of increase of the time available for cooling using the \cite{white_galaxy_1991} method. We return a rate
    of 1, even though technically it can depend on halo properties.
    !!}
    implicit none
    class(coolingTimeAvailableWhiteFrenk1991), intent(inout) :: self
    type (treeNode                          ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    ! Simply return unit rate.
    whiteFrenk1991TimeAvailableIncreaseRate=1.0d0
    return
  end function whiteFrenk1991TimeAvailableIncreaseRate
