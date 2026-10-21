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
  Implements a node operator class that performs tidal heating in spheroids.
  !!}
  
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Satellites_Tidal_Fields, only : satelliteTidalFieldClass

  !![
  <nodeOperator name="nodeOperatorTidalHeatingSpheroids" docformat="rst">
   <description>
   Computes tidal heating of satellite galaxy spheroids by the tidal field of the host halo. The heating is applied as a torque which increases the (pseudo-)angular momentum of the spheroid, and so causes it to expand:

   .. math::

      \dot{J}_\mathrm{spheroid} = |\mathcal{T}| ( M_\mathrm{spheroid, gas} + M_\mathrm{spheroid, stars} ) r_\mathrm{spheroid}^2,

   where :math:`\mathcal{T}` is the radial tidal tensor evaluated at pericenter (see :galacticus-class:`satelliteTidalFieldClass`), and :math:`r_\mathrm{spheroid}` is the scale radius of the spheroid. This node operator can not be used together with heating of the spheroid by a :galacticus-class:`galacticComponentHeatingClass` (other than :galacticus-class:`galacticComponentHeatingNone`). No heating is applied if the spheroid has non-positive mass or radius, if its radius already exceeds the virial radius of the halo, or if its angular momentum already exceeds :math:`(M_\mathrm{spheroid, gas} + M_\mathrm{spheroid, stars}) r_\mathrm{virial} V_\mathrm{virial}`.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorTidalHeatingSpheroids
     !!{RST
     A node operator class that performs tidal heating in spheroids.
     !!}
     private
     class(satelliteTidalFieldClass), pointer :: satelliteTidalField_ => null()
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::                          tidalHeatingSpheroidsDestructor
     procedure :: differentialEvolution => tidalHeatingSpheroidsDifferentialEvolution
  end type nodeOperatorTidalHeatingSpheroids

  interface nodeOperatorTidalHeatingSpheroids
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorTidalHeatingSpheroids` node operator class.
     !!}
     module procedure tidalHeatingSpheroidsConstructorParameters
     module procedure tidalHeatingSpheroidsConstructorInternal
  end interface nodeOperatorTidalHeatingSpheroids

contains
  
  function tidalHeatingSpheroidsConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorTidalHeatingSpheroids` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorTidalHeatingSpheroids)                :: self
    type (inputParameters                  ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass         ), pointer       :: darkMatterHaloScale_
    class(satelliteTidalFieldClass         ), pointer       :: satelliteTidalField_

    call tidalHeatingSpheroidsExclusivityCheck(parameters)
    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="satelliteTidalField" name="satelliteTidalField_" source="parameters"/>
    !!]
    self=nodeOperatorTidalHeatingSpheroids(darkMatterHaloScale_,satelliteTidalField_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="satelliteTidalField_"/>
    !!]
    return
  end function tidalHeatingSpheroidsConstructorParameters

  subroutine tidalHeatingSpheroidsExclusivityCheck(parameters)
    !!{RST
    Check that heating of spheroids by a :galacticus-class:`galacticComponentHeatingClass` is not also in use. Heating by this
    node operator and by such a class would double count the heating.
    !!}
    use :: Error             , only : Error_Report
    use :: Input_Parameters  , only : inputParameters
    use :: ISO_Varying_String, only : varying_string , operator(/=)   , assignment(=)
    implicit none
    type(inputParameters), intent(inout) :: parameters
    type(inputParameters), pointer       :: parametersRoot
    type(inputParameters)                :: parametersSpheroid
    type(varying_string )                :: galacticComponentHeating

    ! Find the root parameters. The spheroid component finds its heating class in its own parameters, or, failing that, in the root
    ! parameters, so check both in that order.
    parametersRoot => parameters%parent
    if (.not.associated(parametersRoot)) return
    do while (associated(parametersRoot%parent))
       parametersRoot => parametersRoot%parent
    end do
    galacticComponentHeating='none'
    if (parametersRoot%isPresent('componentSpheroid')) then
       parametersSpheroid=parametersRoot%subParameters('componentSpheroid')
       if (parametersSpheroid%isPresent('galacticComponentHeating')) then
          call parametersSpheroid%value('galacticComponentHeating',galacticComponentHeating)
       else if (parametersRoot%isPresent('galacticComponentHeating')) then
          call parametersRoot    %value('galacticComponentHeating',galacticComponentHeating)
       end if
    else if (parametersRoot%isPresent('galacticComponentHeating')) then
       call parametersRoot%value('galacticComponentHeating',galacticComponentHeating)
    end if
    if (galacticComponentHeating /= 'none') call Error_Report('this node operator can not be used together with heating of the spheroid by a `galacticComponentHeating` class'//{introspection:location})
    return
  end subroutine tidalHeatingSpheroidsExclusivityCheck

  function tidalHeatingSpheroidsConstructorInternal(darkMatterHaloScale_,satelliteTidalField_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorTidalHeatingSpheroids` node operator class.
    !!}
    implicit none
    type (nodeOperatorTidalHeatingSpheroids)                        :: self
    class(darkMatterHaloScaleClass         ), intent(in   ), target :: darkMatterHaloScale_
    class(satelliteTidalFieldClass         ), intent(in   ), target :: satelliteTidalField_
    !![
    <constructorAssign variables="*darkMatterHaloScale_, *satelliteTidalField_"/>
    !!]

    return
  end function tidalHeatingSpheroidsConstructorInternal

  subroutine tidalHeatingSpheroidsDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorTidalHeatingSpheroids` node operator class.
    !!}
    implicit none
    type(nodeOperatorTidalHeatingSpheroids), intent(inout) :: self
    
    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%satelliteTidalField_"/>
    !!]
    return
  end subroutine tidalHeatingSpheroidsDestructor
  
  subroutine tidalHeatingSpheroidsDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{RST
    Perform star formation in a spheroid.
    !!}
    use :: Galacticus_Nodes, only : propertyInactive, nodeComponentSpheroid
    implicit none
    class           (nodeOperatorTidalHeatingSpheroids), intent(inout), target  :: self
    type            (treeNode                         ), intent(inout), target  :: node
    logical                                            , intent(inout)          :: interrupt
    procedure       (interruptTask                    ), intent(inout), pointer :: functionInterrupt
    integer                                            , intent(in   )          :: propertyType
    class           (nodeComponentSpheroid            )               , pointer :: spheroid
    double precision                                                            :: tidalField       , tidalTorque
    
    ! Do nothing during inactive property solving.
    if     (propertyInactive(propertyType)      ) return
    ! Return if the node is not a satellite.
    if     (.not.node%isSatellite()             ) return
    ! Return if the spheroid is unphysical.
    spheroid => node%spheroid()
    if     (    spheroid%angularMomentum() < 0.0d0 &
         & .or. spheroid%radius         () < 0.0d0 &
         & .or. spheroid%massGas        () < 0.0d0 &
         & .or. spheroid%massStellar    () < 0.0d0 &
         & ) return
    ! Return if the spheroid has excessive angular momentum or radius.
    if     (                                                                                                                                                                      &
         &   spheroid%angularMomentum() > (spheroid%massGas()+spheroid%massStellar())*self%darkMatterHaloScale_%radiusVirial(node)*self%darkMatterHaloScale_%velocityVirial(node) &
         &  .or.                                                                                                                                                                  &
         &   spheroid%radius         () >                                             self%darkMatterHaloScale_%radiusVirial(node)                                                &
         & ) return
    tidalField =self%satelliteTidalField_%tidalTensorRadial(node,atPericenter=.true.,includeCentrifugalAcceleration=.false.)
    tidalTorque=+abs(tidalField)              &
         &      *(                           &
         &        +spheroid%massGas    ()    &
         &        +spheroid%massStellar()    &
         &       )                           &
         &      *  spheroid%radius     ()**2
    call spheroid%angularMomentumRate(tidalTorque)
    return
  end subroutine tidalHeatingSpheroidsDifferentialEvolution

