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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Implements a galactic filter on lightcone geometry.
!!}

  use :: Geometry_Lightcones, only : geometryLightconeClass

  !![
  <galacticFilter name="galacticFilterLightcone" docformat="rst">
   <description>
   Passes nodes that fall within a lightcone geometry as defined by a :galacticus-class:`geometryLightconeClass` object, enabling mock catalog construction by selecting only galaxies observable within the survey volume.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterLightcone
     !!{RST
     A galactic filter class on lightcone geometry.
     !!}
     private
     class(geometryLightconeClass), pointer :: geometryLightcone_ => null()
   contains
     final     ::           lightconeDestructor
     procedure :: passes => lightconePasses
  end type galacticFilterLightcone

  interface galacticFilterLightcone
     !!{RST
     Constructors for the :galacticus-class:`galacticFilterLightcone` galactic filter class.
     !!}
     module procedure lightconeConstructorParameters
     module procedure lightconeConstructorInternal
  end interface galacticFilterLightcone

contains

  function lightconeConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`galacticFilterLightcone` galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (galacticFilterLightcone)                :: self
    type (inputParameters        ), intent(inout) :: parameters
    class(geometryLightconeClass ), pointer       :: geometryLightcone_

    !![
    <objectBuilder class="geometryLightcone" name="geometryLightcone_" source="parameters"/>
    !!]
    self=galacticFilterLightcone(geometryLightcone_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="geometryLightcone_"/>
    !!]
    return
  end function lightconeConstructorParameters

  function lightconeConstructorInternal(geometryLightcone_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`galacticFilterLightcone` galactic filter class.
    !!}
    use :: Display            , only : displayGreen         , displayMagenta, displayMessage, displayReset
    use :: Geometry_Lightcones, only : geometryLightconeNull
    use :: ISO_Varying_String , only : varying_string       , assignment(=) , operator(//)
    implicit none
    type (galacticFilterLightcone)                        :: self
    class(geometryLightconeClass ), intent(in   ), target :: geometryLightcone_
    type (varying_string         )                        :: message
    !![
    <constructorAssign variables="*geometryLightcone_"/>
    !!]

    ! Warn if the lightcone geometry is null. This is a valid choice - a null lightcone has zero solid angle, so nothing lies
    ! within it - but it makes this filter reject every galaxy, which is unlikely to be what was intended.
    select type (geometryLightcone_)
    class is (geometryLightconeNull)
       message=displayMagenta()//"WARNING:"//displayReset()//" the `galacticFilterLightcone` class is in use with a null lightcone geometry"        //char(10)// &
            &  displayGreen()//"    HELP:"//displayReset()//" a null lightcone has zero solid angle, so no galaxy ever lies within it and this"     //char(10)// &
            &                                               "          filter will therefore reject every galaxy. If that is not what you intended,"//char(10)// &
            &                                               "          set a non-null `geometryLightcone`."
       call displayMessage(message)
    end select
    return
  end function lightconeConstructorInternal

  subroutine lightconeDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`galacticFilterLightcone` galactic filter class.
    !!}
    implicit none
    type(galacticFilterLightcone), intent(inout) :: self

    !![
    <objectDestructor name="self%geometryLightcone_"/>
    !!]
    return
  end subroutine lightconeDestructor

  logical function lightconePasses(self,node)
    !!{RST
    Implement a lightcone geometry galactic filter.
    !!}
    implicit none
    class(galacticFilterLightcone), intent(inout)         :: self
    type (treeNode               ), intent(inout), target :: node

    lightconePasses=self%geometryLightcone_%isInLightcone(node,atPresentEpoch=.true.)
    return
  end function lightconePasses
