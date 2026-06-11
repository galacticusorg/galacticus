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

!!{RST
Implements a galactic high-pass filter for total ISM mass.
!!}

  !![
  <galacticFilter name="galacticFilterISMMass" docformat="rst">
   <description>
   A galactic high-pass filter for ISM mass. Galaxies with a combined disk plus spheroid ISM mass greater than or equal to a fixed threshold, :math:`M_\mathrm{ISM,0}=`\ ``[massThreshold]``.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterISMMass
     !!{RST
     A galactic high-pass filter class for ISM mass.
     !!}
     private
     double precision :: massThreshold
   contains
     procedure :: passes => ismMassPasses
  end type galacticFilterISMMass

  interface galacticFilterISMMass
     !!{RST
     Constructors for the ``galacticFilterISMMass`` galactic filter class.
     !!}
     module procedure ismMassConstructorParameters
     module procedure ismMassConstructorInternal
  end interface galacticFilterISMMass

contains

  function ismMassConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``galacticFilterISMMass`` galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterISMMass)                :: self
    type            (inputParameters      ), intent(inout) :: parameters
    double precision                                       :: massThreshold
    
    !![
    <inputParameter docformat="rst">
      <name>massThreshold</name>
      <source>parameters</source>
      <description>
      The parameter :math:`M_0` (in units of :math:`\mathrm{M}_\odot`) appearing in the ISM mass threshold for the ISM mass galactic filter class.
      </description>
    </inputParameter>
    !!]
    self=galacticFilterISMMass(massThreshold)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function ismMassConstructorParameters

  function ismMassConstructorInternal(massThreshold) result(self)
    !!{RST
    Internal constructor for the ``galacticFilterISMMass`` galactic filter class.
    !!}
    implicit none
    type            (galacticFilterISMMass)                :: self
    double precision                       , intent(in   ) :: massThreshold
    !![
    <constructorAssign variables="massThreshold"/>
    !!]
    
    return
  end function ismMassConstructorInternal

  logical function ismMassPasses(self,node)
    !!{RST
    Implement an ismMass-pass galactic filter.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, treeNode
    implicit none
    class           (galacticFilterISMMass), intent(inout)         :: self
    type            (treeNode             ), intent(inout), target :: node
    class           (nodeComponentDisk    ), pointer               :: disk
    class           (nodeComponentSpheroid), pointer               :: spheroid
    double precision                                               :: ismMass

    disk          => node    %disk    ()
    spheroid      => node    %spheroid()
    ismMass       = +disk    %massGas () &
         &          +spheroid%massGas ()
    ismMassPasses =  ismMass             &
         &          >=                   &
         &           self%massThreshold
    return
  end function ismMassPasses
