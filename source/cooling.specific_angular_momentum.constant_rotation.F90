!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which calculates the specific angular momentum of cooling gas assuming a constant rotation velocity as a
!% function of radius.

module Cooling_Specific_Angular_Momenta_Constant_Rotation
  !% Calculates the specific angular momentum of cooling gas assuming a constant rotation velocity as a function of radius.
  use, intrinsic :: ISO_C_Binding
  use Kind_Numbers
  implicit none
  private
  public :: Cooling_Specific_AM_Constant_Rotation_Initialize, Cooling_Specific_AM_Constant_Rotation_Reset

  ! Record of unique ID of node which we last computed results for.
  integer         (kind=kind_int8)            :: lastUniqueID                          =-1
  !$omp threadprivate(lastUniqueID)
  ! Record of whether or not specific angular momentum of cooling has has already been computed for this node.
  logical                                     :: coolingSpecificAngularMomentumComputed=.false.
  !$omp threadprivate(coolingSpecificAngularMomentumComputed)
  ! Stored values of specific angular momentum of cooling gas.
  double precision                            :: coolingSpecificAngularMomentumStored
  !$omp threadprivate(coolingSpecificAngularMomentumStored)
  ! Parameters controlling the calculation.
  integer                         , parameter :: profileDarkMatter                     =0
  integer                         , parameter :: profileHotGas                         =1
  integer                                     :: meanSpecificAngularMomentumFrom               , rotationNormalizationFrom
  logical                                     :: coolingAngularMomentumUseInteriorMean

contains

  !# <coolingSpecificAngularMomentumMethod>
  !#  <unitName>Cooling_Specific_AM_Constant_Rotation_Initialize</unitName>
  !# </coolingSpecificAngularMomentumMethod>
  subroutine Cooling_Specific_AM_Constant_Rotation_Initialize(coolingSpecificAngularMomentumMethod,Cooling_Specific_Angular_Momentum_Get)
    !% Initializes the ``constant rotation'' specific angular momentum of cooling gas module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                                     ), intent(in   )          :: coolingSpecificAngularMomentumMethod
    procedure(Cooling_Specific_Angular_Momentum_Constant_Rotation), intent(inout), pointer :: Cooling_Specific_Angular_Momentum_Get
    type     (varying_string                                     )                         :: inputOption

    if (coolingSpecificAngularMomentumMethod == 'constantRotation') then
       Cooling_Specific_Angular_Momentum_Get => Cooling_Specific_Angular_Momentum_Constant_Rotation

       !@ <inputParameter>
       !@   <name>coolingMeanAngularMomentumFrom</name>
       !@   <defaultValue>hotGas</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The component (``{\tt hotGas}'' or ``{\tt darkMatter}'') from which the mean specific angular momentum should be computed for
       !@     calculations of cooling gas specific angular momentum.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingMeanAngularMomentumFrom',inputOption,defaultValue='hotGas')
       select case (char(inputOption))
       case ("darkMatter")
          meanSpecificAngularMomentumFrom=profileDarkMatter
       case ("hotGas")
          meanSpecificAngularMomentumFrom=profileHotGas
       end select

       !@ <inputParameter>
       !@   <name>coolingRotationVelocityFrom</name>
       !@   <defaultValue>hotGas</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The component (``{\tt hotGas}'' or ``{\tt darkMatter}'') from which the constant rotation speed should be computed for
       !@     calculations of cooling gas specific angular momentum.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingRotationVelocityFrom',inputOption,defaultValue='hotGas')
       select case (char(inputOption))
       case ("darkMatter")
          rotationNormalizationFrom=profileDarkMatter
       case ("hotGas")
          rotationNormalizationFrom=profileHotGas
       end select

       !@ <inputParameter>
       !@   <name>coolingAngularMomentumUseInteriorMean</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether to use the specific angular momentum at the cooling radius, or the mean specific angular momentum interior to that radius.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingAngularMomentumUseInteriorMean',coolingAngularMomentumUseInteriorMean,defaultValue=.false.)

    end if
    return
  end subroutine Cooling_Specific_AM_Constant_Rotation_Initialize

  !# <calculationResetTask>
  !# <unitName>Cooling_Specific_AM_Constant_Rotation_Reset</unitName>
  !# </calculationResetTask>
  subroutine Cooling_Specific_AM_Constant_Rotation_Reset(thisNode)
    !% Reset the specific angular momentum of cooling gas calculation.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    coolingSpecificAngularMomentumComputed=.false.
    lastUniqueID                          =thisNode%uniqueID()
    return
  end subroutine Cooling_Specific_AM_Constant_Rotation_Reset

  double precision function Cooling_Specific_Angular_Momentum_Constant_Rotation(thisNode,radius)
    !% Return the specific angular momentum of cooling gas in the constant rotation model.
    use Galacticus_Nodes
    use Cooling_Infall_Radii
    use Dark_Matter_Profiles
    use Hot_Halo_Density_Profile
    use Numerical_Constants_Physical
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    double precision                      , intent(in   )          :: radius
    class           (nodeComponentBasic  )               , pointer :: thisBasicComponent
    class           (nodeComponentSpin   )               , pointer :: thisSpinComponent
    class           (nodeComponentHotHalo)               , pointer :: thisHotHaloComponent
    double precision                                               :: meanSpecificAngularMomentum, rotationNormalization

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Cooling_Specific_AM_Constant_Rotation_Reset(thisNode)

    ! Check if specific angular momentum of cooling gas is already computed.
    if (.not.coolingSpecificAngularMomentumComputed) then
       ! Flag that cooling radius is now computed.
       coolingSpecificAngularMomentumComputed=.true.

       ! Compute the mean specific angular momentum.
       select case (meanSpecificAngularMomentumFrom)
       case (profileDarkMatter)
          ! Compute mean specific angular momentum of the dark matter halo from the spin parameter, mass and energy of the halo.
          thisBasicComponent => thisNode%basic()
          thisSpinComponent  => thisNode%spin ()
          meanSpecificAngularMomentum= gravitationalConstantGalacticus                   &
               &                      *thisSpinComponent %spin()                         &
               &                      *thisBasicComponent%mass()**1.5d0                  &
               &                      /sqrt(abs(Dark_Matter_Profile_Energy(thisNode)))
       case (profileHotGas    )
          ! Compute mean specific angular momentum from the hot halo component.
          thisHotHaloComponent => thisNode%hotHalo()
          meanSpecificAngularMomentum= thisHotHaloComponent%angularMomentum() &
               &                      /thisHotHaloComponent%mass           ()
       end select

       ! Compute the rotation normalization.
       select case (rotationNormalizationFrom      )
       case (profileDarkMatter)
          rotationNormalization=Dark_Matter_Profile_Rotation_Normalization(thisNode)
       case (profileHotGas    )
          rotationNormalization=Hot_Halo_Profile_Rotation_Normalization   (thisNode)
       end select

       ! Compute the specific angular momentum of the cooling gas.
       coolingSpecificAngularMomentumStored= rotationNormalization       &
            &                               *meanSpecificAngularMomentum
    end if

    ! Check that the radius is positive.
    if (radius > 0.0d0) then
       ! Return the computed value.
       if (coolingAngularMomentumUseInteriorMean) then
          ! Find the specific angular momentum interior to the specified radius.
          Cooling_Specific_Angular_Momentum_Constant_Rotation= coolingSpecificAngularMomentumStored                  &
               &                                              *Hot_Halo_Profile_Radial_Moment(thisNode,3.0d0,radius) &
               &                                              /Hot_Halo_Profile_Radial_Moment(thisNode,2.0d0,radius)
       else
          ! Find the specific angular momentum at the specified radius.
          Cooling_Specific_Angular_Momentum_Constant_Rotation= coolingSpecificAngularMomentumStored                  &
               &                                              *                                              radius
       end if
    else
       ! Radius is non-positive - return zero.
       Cooling_Specific_Angular_Momentum_Constant_Rotation=0.0d0
    end if
    return
  end function Cooling_Specific_Angular_Momentum_Constant_Rotation

end module Cooling_Specific_Angular_Momenta_Constant_Rotation
