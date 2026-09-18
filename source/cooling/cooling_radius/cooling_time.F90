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

  !+    Contributions to this file made by: Claude.

  !!{RST
  Implements an abstract cooling radius class for models in which the cooling radius is found by comparing the cooling time with the time available for cooling.
  !!}

  use :: Cooling_Times          , only : coolingTimeClass
  use :: Cooling_Times_Available, only : coolingTimeAvailableClass
  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Kind_Numbers           , only : kind_int8
  use :: Radiation_Fields       , only : radiationFieldCosmicMicrowaveBackground

  !![
  <coolingRadius name="coolingRadiusCoolingTime" abstract="yes" docformat="rst">
   <description>
   An abstract cooling radius class for models in which the cooling radius is the radius at which the cooling time of gas in the hot atmosphere (see :galacticus-class:`coolingTime`) equals the time available for cooling (see :galacticus-class:`coolingTimeAvailable`). The cooling radius and its growth rate are memoized for the most recent node, and reset by the calculation reset event. This is an abstract class---the cooling radius and its growth rate must be provided by a concrete class.
   </description>
   <deepCopy>
    <functionClass variables="radiation"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="radiation"/>
   </stateStorable>
  </coolingRadius>
  !!]
  type, abstract, extends(coolingRadiusClass) :: coolingRadiusCoolingTime
     !!{RST
     An abstract cooling radius class for models in which the cooling radius is the radius at which the cooling time equals the time available for cooling.
     !!}
     private
     class           (cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_    => null()
     class           (coolingTimeAvailableClass              ), pointer :: coolingTimeAvailable_  => null()
     class           (coolingTimeClass                       ), pointer :: coolingTime_           => null()
     type            (radiationFieldCosmicMicrowaveBackground), pointer :: radiation              => null()
     integer         (kind=kind_int8                         )          :: lastUniqueID           =  -1_kind_int8
     integer                                                            :: abundancesCount                       , chemicalsCount
     ! Stored values of cooling radius.
     logical                                                            :: radiusComputed                        , radiusGrowthRateComputed
     double precision                                                   :: radiusGrowthRateStored                , radiusStored
   contains
     !![
     <methods docformat="rst">
       <method description="Reset memoized calculations." method="calculationReset" />
       <method description="Initialize the state shared by all cooling time-based cooling radius classes. Must be called by the constructor of each concrete class, after its cosmology functions object has been assigned." method="initialize" />
     </methods>
     !!]
     procedure :: autoHook         => coolingTimeAutoHook
     procedure :: calculationReset => coolingTimeCalculationReset
     procedure :: initialize       => coolingTimeInitialize
  end type coolingRadiusCoolingTime

contains

  subroutine coolingTimeInitialize(self)
    !!{RST
    Initialize the state shared by all cooling time-based cooling radius classes.
    !!}
    use :: Abundances_Structure         , only : Abundances_Property_Count
    use :: Array_Utilities              , only : operator(.intersection.)
    use :: Chemical_Abundances_Structure, only : Chemicals_Property_Count
    use :: Error                        , only : Component_List           , Error_Report
    use :: Galacticus_Nodes             , only : defaultHotHaloComponent
    implicit none
    class(coolingRadiusCoolingTime), intent(inout) :: self

    ! Initial state of stored solutions.
    self%radiusComputed          =.false.
    self%radiusGrowthRateComputed=.false.
    ! Get a count of the number of abundances and chemicals properties.
    self%abundancesCount=Abundances_Property_Count()
    self%chemicalsCount =Chemicals_Property_Count ()
    ! Initialize radiation field.
    allocate(self%radiation)
    !![
    <referenceConstruct owner="self" object="radiation" constructor="radiationFieldCosmicMicrowaveBackground(self%cosmologyFunctions_)"/>
    !!]
    ! Check that required components are gettable.
    if     (                                                                                                             &
         &  .not.(                                                                                                       &
         &         defaultHotHaloComponent%       massIsGettable() .and.                                                 &
         &         defaultHotHaloComponent% abundancesIsGettable() .and.                                                 &
         &         defaultHotHaloComponent%outerRadiusIsGettable() .and.                                                 &
         &        (defaultHotHaloComponent%  chemicalsIsGettable() .or.  self%chemicalsCount == 0)                       &
         &       )                                                                                                       &
         & ) call Error_Report                                                                                           &
         & (                                                                                                             &
         &  'This method requires that the "mass", "abundances", "outerRadius", and "chemicals" '//                      &
         &  '(if any chemicals are being used) properties of the hot halo are gettable.'         //                      &
         &  Component_List(                                                                                              &
         &                 'hotHalo'                                                                                  ,  &
         &                  defaultHotHaloComponent%massAttributeMatch       (requireGettable=.true.                 )   &
         &                 .intersection.                                                                                &
         &                  defaultHotHaloComponent%abundancesAttributeMatch (requireGettable=.true.                 )   &
         &                 .intersection.                                                                                &
         &                  defaultHotHaloComponent%outerRadiusAttributeMatch(requireGettable=.true.                 )   &
         &                 .intersection.                                                                                &
         &                  defaultHotHaloComponent%chemicalsAttributeMatch  (requireGettable=self%chemicalsCount > 0)   &
         &                )                                                                                           // &
         &  {introspection:location}                                                                                     &
         & )
    return
  end subroutine coolingTimeInitialize

  subroutine coolingTimeAutoHook(self)
    !!{RST
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(coolingRadiusCoolingTime), intent(inout) :: self

    call calculationResetEvent%attach(self,coolingTimeCalculationReset,openMPThreadBindingAllLevels,label='coolingRadiusCoolingTime')
    return
  end subroutine coolingTimeAutoHook

  subroutine coolingTimeCalculationReset(self,node,uniqueID)
    !!{RST
    Reset the cooling radius calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (coolingRadiusCoolingTime), intent(inout) :: self
    type   (treeNode                ), intent(inout) :: node
    integer(kind_int8               ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%radiusComputed          =.false.
    self%radiusGrowthRateComputed=.false.
    self%lastUniqueID            =uniqueID
    return
  end subroutine coolingTimeCalculationReset
