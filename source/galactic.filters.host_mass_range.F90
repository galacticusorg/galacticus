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
Implements a galactic filter which passes nodes with host halo basic mass within a specified range.
!!}

  !![
  <galacticFilter name="galacticFilterHostMassRange">
   <description>Passes nodes with host halo basic mass, $M_\mathrm{host}$, in the range {\normalfont \ttfamily [massMinimum]}$\le M_\mathrm{host}&lt;${\normalfont \ttfamily [massMaximum]}.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterHostMassRange
     !!{
     A galactic filter which passes nodes with host halo basic mass within a specified range.
     !!}
     private
     double precision :: massMinimum , massMaximum
     logical          :: useFinalHost
   contains
     procedure :: passes => hostMassRangePasses
  end type galacticFilterHostMassRange

  interface galacticFilterHostMassRange
     !!{
     Constructors for the \refClass{galacticFilterHostMassRange} galactic filter class.
     !!}
     module procedure hostMassRangeConstructorParameters
     module procedure hostMassRangeConstructorInternal
  end interface galacticFilterHostMassRange

contains

  function hostMassRangeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterHostMassRange} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterHostMassRange)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    double precision                                             :: massMinimum , massMaximum
    logical                                                      :: useFinalHost

    !![
    <inputParameter>
      <name>massMinimum</name>
      <source>parameters</source>
      <description>The minimum mass of host halo to pass.</description>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <source>parameters</source>
      <description>The maximum mass of host halo to pass.</description>
    </inputParameter>
    <inputParameter>
      <name>useFinalHost</name>
      <source>parameters</source>
      <description>If true, the final host (i.e. the isolated host halo in the subhalo hierarchy) is used for filtering, otherwise the immediate host is used.</description>
    </inputParameter>
    !!]
    self=galacticFilterHostMassRange(massMinimum,massMaximum,useFinalHost)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function hostMassRangeConstructorParameters

  function hostMassRangeConstructorInternal(massMinimum,massMaximum,useFinalHost) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterHostMassRange} galactic filter class.
    !!}
    implicit none
    type            (galacticFilterHostMassRange)                :: self
    double precision                             , intent(in   ) :: massMinimum , massMaximum
    logical                                      , intent(in   ) :: useFinalHost
    !![
    <constructorAssign variables="massMinimum, massMaximum, useFinalHost"/>
    !!]

    return
  end function hostMassRangeConstructorInternal

  logical function hostMassRangePasses(self,node)
    !!{
    Implement a filter which passes nodes with host halo basic mass in a specified range.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(galacticFilterHostMassRange), intent(inout)          :: self
    type (treeNode                   ), intent(inout), target  :: node
    type (treeNode                   )               , pointer :: nodeHost
    class(nodeComponentBasic         )               , pointer :: basic

    if (node%isSatellite()) then
       nodeHost => node%parent
       if (self%useFinalHost) then
          do while (associated(nodeHost%parent))
             nodeHost => nodeHost%parent
          end do
       end if
    else
       nodeHost => node
    end if
    basic => nodeHost%basic()
    hostMassRangePasses= basic%mass() >= self %massMinimum &
         &              .and.                              &
         &               basic%mass() <  self %massMaximum
    return
  end function hostMassRangePasses
