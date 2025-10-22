!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a module which implements a dark matter profile method that provides no properties (but does provide a mass distribution factory).
!!}

module Node_Component_Dark_Matter_Profile_Scale_Free
  !!{
  Implements a dark matter profile method that provides no properties (but does provide a mass distribution factory).
  !!}
  use :: Dark_Matter_Profiles    , only : darkMatterProfileClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Mass_Distributions      , only : massDistributionClass
  implicit none
  private
  public :: Node_Component_Dark_Matter_Profile_Scale_Free_Thread_Init, Node_Component_Dark_Matter_Profile_Scale_Free_Thread_Uninit, &
       &    Node_Component_Dark_Matter_Profile_Scale_Free_State_Store, Node_Component_Dark_Matter_Profile_Scale_Free_State_Restore, &
       &    Node_Component_Dark_Matter_Profile_Scale_Free_Init

  !![
  <component>
   <class>darkMatterProfile</class>
   <name>scaleFree</name>
   <isDefault>false</isDefault>
   <bindings>
     <binding method="massDistribution" bindsTo="component" isDeferred="true" >
      <interface>
       <type>class(massDistributionClass), pointer</type>
       <rank>0</rank>
       <module>Galactic_Structure_Options, only : enumerationWeightByType, enumerationComponentTypeType, enumerationMassTypeType</module>
       <module>Mass_Distributions        , only : massDistributionClass                                                         </module>
       <self pass="true" intent="inout" />
       <argument>type   (enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
       <argument>type   (enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
       <argument>type   (enumerationWeightByType     ), intent(in   ), optional :: weightBy     </argument>
       <argument>integer                              , intent(in   ), optional :: weightIndex  </argument>
      </interface>
     </binding>
   </bindings>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterProfileClass   ), pointer :: darkMatterProfile_
  class(darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_
  !$omp threadprivate(darkMatterProfile_,darkMatterProfileDMO_)
  
  ! Procedure pointers to mass distribution functions.
  procedure(Node_Component_Dark_Matter_Profile_Scale_Free_Mass_Dist), pointer :: Node_Component_Dark_Matter_Profile_Scale_Free_Mass_Dist_

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_Free_Init</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_Free_Init(parameters)
    !!{
    Initializes the scale dark matter profile component.
    !!}
    use :: Input_Parameters, only : inputParameters
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent, nodeComponentDarkMatterProfileScaleFree
    type(inputParameters                        ), intent(inout) :: parameters
    type(nodeComponentDarkMatterProfileScaleFree)                :: darkMatterProfile
    !$GLC attributes unused :: parameters

    !$omp critical (Node_Component_Dark_Matter_Profile_Init)
    if (defaultDarkMatterProfileComponent%scaleFreeIsActive()) then
       Node_Component_Dark_Matter_Profile_Scale_Free_Mass_Dist_ => Node_Component_Dark_Matter_Profile_Scale_Free_Mass_Dist
       call darkMatterProfile%massDistributionFunction(Node_Component_Dark_Matter_Profile_Scale_Free_Mass_Dist_)
    end if
    !$omp end critical (Node_Component_Dark_Matter_Profile_Init)
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Free_Init
    
  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_Free_Thread_Init</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_Free_Thread_Init(parameters)
    !!{
    Initializes the tree node scale dark matter profile module.
    !!}
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    use :: Input_Parameters, only : inputParameter                   , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters

    if (defaultDarkMatterProfileComponent%scaleFreeIsActive()) then
       !![
       <objectBuilder class="darkMatterProfile"    name="darkMatterProfile_"    source="parameters"/>
       <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
       !!]
     end if
     return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Free_Thread_Init

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_Free_Thread_Uninit</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_Free_Thread_Uninit()
    !!{
    Uninitializes the tree node scale dark matter profile module.
    !!}
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    implicit none

    if (defaultDarkMatterProfileComponent%scaleFreeIsActive()) then
       !![
       <objectDestructor name="darkMatterProfile_"   />
       <objectDestructor name="darkMatterProfileDMO_"/>
       !!]
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Free_Thread_Uninit

  !![
  <stateStoreTask>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_Free_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_Free_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentDarkMatterProfile -> scaleFree',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="darkMatterProfile_ darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Free_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_Free_State_Restore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_Free_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentDarkProfile -> scaleFree',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="darkMatterProfile_ darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Free_State_Restore

  function Node_Component_Dark_Matter_Profile_Scale_Free_Mass_Dist(self,componentType,massType,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the mass distribution associated with the dark matter profile.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentDarkMatterProfileScaleFree
    use :: Galactic_Structure_Options, only : enumerationWeightByType                , enumerationComponentTypeType, enumerationMassTypeType, componentTypeAll, &
         &                                    componentTypeDarkHalo                  , componentTypeDarkMatterOnly , massTypeAll            , massTypeDark
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class  (massDistributionClass                  ), pointer                 :: massDistribution_
    class  (nodeComponentDarkMatterProfileScaleFree), intent(inout)           :: self
    type   (enumerationComponentTypeType           ), intent(in   ), optional :: componentType
    type   (enumerationMassTypeType                ), intent(in   ), optional :: massType
    type   (enumerationWeightByType                ), intent(in   ), optional :: weightBy
    integer                                         , intent(in   ), optional :: weightIndex
    !![
    <optionalArgument name="componentType" defaultsTo="componentTypeAll"/>
    <optionalArgument name="massType"      defaultsTo="massTypeAll"     />
    !!]

    massDistribution_ => null()
    if         (                                              &
         &       massType_      == massTypeAll                &
         &      .or.                                          &
         &       massType_      == massTypeDark               &
         &     ) then
       if      (                                              &
            &    componentType_ == componentTypeAll           &
            &   .or.                                          &
            &    componentType_ == componentTypeDarkHalo      &
            &  ) then
          massDistribution_ => darkMatterProfile_   %get(self%hostNode,weightBy,weightIndex)
       else if (                                              &
            &   componentType_ == componentTypeDarkMatterOnly &
            &  ) then
          massDistribution_ => darkMatterProfileDMO_%get(self%hostNode,weightBy,weightIndex)
          call massDistribution_%setTypes(componentType=componentTypeDarkMatterOnly)
       end if
    end if
    return
  end function Node_Component_Dark_Matter_Profile_Scale_Free_Mass_Dist

end module Node_Component_Dark_Matter_Profile_Scale_Free
