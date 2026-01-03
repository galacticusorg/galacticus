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
  Implements a node operator class that computes the accretion rate onto a halo corresponding to the ``Bertschinger''
  mass. Typically this is the mass corresponding to a spherical top-hat collapse scenario, although in practice this class can be
  provided with any \refClass{virialDensityContrastClass}.
  !!}

  !![
  <nodeOperator name="nodeOperatorBertschingerMass">
   <description>
     A node operator class that computes the accretion rate onto a halo corresponding to the ``Bertschinger'' mass. Typically
     this is the mass corresponding to a spherical top-hat collapse scenario, although in practice this class can be provided with
     any \refClass{virialDensityContrastClass}.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorBertschingerMass
     !!{
     A node operator class that computes the accretion rate onto a halo corresponding to the ``Bertschinger'' mass. Typically
     this is the mass corresponding to a spherical top-hat collapse scenario, although in practice this class can be provided with
     any \refClass{virialDensityContrastClass}.
     !!}
     private
     class  (cosmologyParametersClass  ), pointer :: cosmologyParameters_        => null()
     class  (cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_         => null()
     class  (virialDensityContrastClass), pointer :: virialDensityContrast_      => null()
     integer                                      :: massBertschingerID                   , massBertschingerTargetID, &
          &                                          accretionRateBertschingerID
   contains
     final     ::                                        bertschingerMassDestructor
     procedure :: nodeTreeInitialize                  => bertschingerMassNodeTreeInitialize
     procedure :: nodeInitialize                      => bertschingerMassNodeInitialize
     procedure :: nodePromote                         => bertschingerMassNodePromote
     procedure :: nodesMerge                          => bertschingerMassNodesMerge
     procedure :: differentialEvolutionAnalytics      => bertschingerMassDifferentialEvolutionAnalytics
     procedure :: differentialEvolutionSolveAnalytics => bertschingerMassDifferentialEvolutionSolveAnalytics
  end type nodeOperatorBertschingerMass
  
  interface nodeOperatorBertschingerMass
     !!{
     Constructors for the \refClass{nodeOperatorBertschingerMass} node operator class.
     !!}
     module procedure bertschingerMassConstructorParameters
     module procedure bertschingerMassConstructorInternal
  end interface nodeOperatorBertschingerMass
  
contains
  
  function bertschingerMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorBertschingerMass} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorBertschingerMass)                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(cosmologyParametersClass    ), pointer       :: cosmologyParameters_
    class(cosmologyFunctionsClass     ), pointer       :: cosmologyFunctions_
    class(virialDensityContrastClass  ), pointer       :: virialDensityContrast_

    !![
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self=nodeOperatorBertschingerMass(cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function bertschingerMassConstructorParameters

  function bertschingerMassConstructorInternal(cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorBertschingerMass} node operator class.
    !!}
    implicit none
    type (nodeOperatorBertschingerMass)                        :: self
    class(cosmologyParametersClass    ), intent(in   ), target :: cosmologyParameters_
    class(cosmologyFunctionsClass     ), intent(in   ), target :: cosmologyFunctions_
    class(virialDensityContrastClass  ), intent(in   ), target :: virialDensityContrast_
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologyFunctions_, *virialDensityContrast_"/>
    !!]
    
    !![
    <addMetaProperty component="basic" name="massBertschinger"          id="self%massBertschingerID"          isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="basic" name="massBertschingerTarget"    id="self%massBertschingerTargetID"    isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="basic" name="accretionRateBertschinger" id="self%accretionRateBertschingerID" isEvolvable="no"  isCreator="yes"/>
    !!]
    return
  end function bertschingerMassConstructorInternal

  subroutine bertschingerMassDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorBertschingerMass} node operator class.
    !!}
    implicit none
    type(nodeOperatorBertschingerMass), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]
    return
  end subroutine bertschingerMassDestructor

  subroutine bertschingerMassNodeTreeInitialize(self,node)
    !!{
    Initialize the Bertschinger mass of all nodes in the tree.    
    !!}
    implicit none
    class(nodeOperatorBertschingerMass), intent(inout), target  :: self
    type (treeNode                    ), intent(inout), target  :: node
    
    call self%nodeInitialize(node)
    return
  end subroutine bertschingerMassNodeTreeInitialize
  
  recursive subroutine bertschingerMassNodeInitialize(self,node)
    !!{
    Compute the rate of growth of the ``\gls{dmou}'' mass of a halo assuming a constant growth rate.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic
    implicit none
    class           (nodeOperatorBertschingerMass), intent(inout), target  :: self
    type            (treeNode                    ), intent(inout), target  :: node
    type            (treeNode                    )               , pointer :: nodeChild          , nodeParent
    class           (nodeComponentBasic          )               , pointer :: basicChild         , basicParent   , &
         &                                                                    basic
    double precision                                                       :: timeDelta          , massUnresolved, &
         &                                                                    massTotalProgenitor
    !$GLC attributes unused :: self

    basic => node%basic()
    ! Set the Bertschinger mass of the node.
    call basic%floatRank0MetaPropertySet(                                                                                                                                 &
         &                               self%massBertschingerID                                                                                                        , &
         &                               Dark_Matter_Profile_Mass_Definition(                                                                                             &
         &                                                                                          node                                                                , &
         &                                                                                          self%virialDensityContrast_%densityContrast(                          &
         &                                                                                                                                      basic%mass            (), &
         &                                                                                                                                      basic%timeLastIsolated()  &
         &                                                                                                                                     )                        , &
         &                                                                   cosmologyParameters_  =self%cosmologyParameters_                                           , &
         &                                                                   cosmologyFunctions_   =self%cosmologyFunctions_                                            , &
         &                                                                   virialDensityContrast_=self%virialDensityContrast_                                           &
         &                                                                  )                                                                                             &
         &                              )
    ! Determine node status.
    if (node%isSatellite()) then
       ! Node is a satellite - we assume no accretion.
       call basic%floatRank0MetaPropertySet(self%accretionRateBertschingerID,0.0d0                                                   )
       call basic%floatRank0MetaPropertySet(self%   massBertschingerTargetID,basic%floatRank0MetaPropertyGet(self%massBertschingerID))
    else if (.not.associated(node%parent)) then
       ! For parent-less nodes (i.e. the root node of the tree), the rate is set equal to that of the
       ! progenitor, if it has one.
       nodeChild => node%firstChild
       if (associated(nodeChild)) then
          ! Get the basic component of the child node.
          basicChild => nodeChild%basic()
          ! Get the growth rate of the child.
          call basic%floatRank0MetaPropertySet(self%accretionRateBertschingerID,basicChild%floatRank0MetaPropertyGet(self%accretionRateBertschingerID))
          call basic%floatRank0MetaPropertySet(self%   massBertschingerTargetID,basic     %floatRank0MetaPropertyGet(self%         massBertschingerID))
       else
          ! Parentless node has no child - set a zero growth rate.
          call basic%floatRank0MetaPropertySet(self%accretionRateBertschingerID,0.0d0                                                                 )
          call basic%floatRank0MetaPropertySet(self%   massBertschingerTargetID,basic     %floatRank0MetaPropertyGet(self%         massBertschingerID))
       end if
    else
       ! Get the parent node.
       nodeParent => node%parent
       ! Get the basic component of the parent node.
       basicParent => nodeParent%basic()
       ! Compute the unresolved mass.
       massUnresolved=nodeMassUnresolved(nodeParent)
       if (massUnresolved > 0.0d0) then
          ! Positive mass growth - assume this occurs entirely in the main progenitor.
          if (node%isPrimaryProgenitor()) then
             ! Main progenitor - compute required growth rate.
             timeDelta=basicParent%time()-basic%time()
             if (timeDelta > 0.0d0)                                                                                                                                          &
                  & call basic%floatRank0MetaPropertySet(self%accretionRateBertschingerID,                                                         massUnresolved/timeDelta)
             call        basic%floatRank0MetaPropertySet(self%   massBertschingerTargetID,basic%floatRank0MetaPropertyGet(self%massBertschingerID)+massUnresolved          )
          else
             ! Non-main progenitor - assume zero growth rate.
             call        basic%floatRank0MetaPropertySet(self%accretionRateBertschingerID,                                                         0.0d0                   )
             call        basic%floatRank0MetaPropertySet(self%   massBertschingerTargetID,basic%floatRank0MetaPropertyGet(self%massBertschingerID)                         )
          end if
       else
          ! Negative mass growth - assume all progenitors lose mass at proportionally equal rates.
          ! Compute the total mass in progenitors.
          massTotalProgenitor=+Dark_Matter_Profile_Mass_Definition(                                                                                                   &
               &                                                                          nodeParent                                                                , &
               &                                                                          self%virialDensityContrast_%densityContrast(                                &
               &                                                                                                                      basicParent%mass            (), &
               &                                                                                                                      basicParent%timeLastIsolated()  &
               &                                                                                                                     )                              , &
               &                                                   cosmologyParameters_  =self%cosmologyParameters_                                                 , &
               &                                                   cosmologyFunctions_   =self%cosmologyFunctions_                                                  , &
               &                                                   virialDensityContrast_=self%virialDensityContrast_                                                 &
               &                                                  )                                                                                                   &
               &              -massUnresolved
          ! Compute the time available for accretion.
          timeDelta=basicParent%time()-basic%time()
          ! Compute mass growth rate.
          if (timeDelta > 0.0d0)                                                                                          &
               & call basic%floatRank0MetaPropertySet(                                                                    &
               &                                                                       self%accretionRateBertschingerID , &
               &                                      +basic%floatRank0MetaPropertyGet(self%         massBertschingerID)  &
               &                                      *massUnresolved                                                     &
               &                                      /massTotalProgenitor                                                &
               &                                      /timeDelta                                                          &
               &                                     )
          call        basic%floatRank0MetaPropertySet(                                                                    &
               &                                                                       self%   massBertschingerTargetID , &
               &                                      +basic%floatRank0MetaPropertyGet(self%         massBertschingerID)  &
               &                                      +basic%floatRank0MetaPropertyGet(self%         massBertschingerID)  &
               &                                      *massUnresolved                                                     &
               &                                      /massTotalProgenitor                                                &
               &                                     )
       end if
    end if
    return

  contains
    
    recursive double precision function nodeMassUnresolved(node)
      !!{
      Return the unresolved mass for {\normalfont \ttfamily node}.
      !!}
      use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
      implicit none
      type (treeNode          ), intent(inout), pointer :: node
      type (treeNode          )               , pointer :: nodeChild
      class(nodeComponentBasic)               , pointer :: basicChild, basic

      ! Get the basic component.
      basic => node%basic()
      ! Initialize the unresolved mass to the mass of the current node's basic component.
      if (basic%floatRank0MetaPropertyGet(self%massBertschingerID) == 0.0d0) call self%nodeInitialize(node)
      nodeMassUnresolved=basic%floatRank0MetaPropertyGet(self%massBertschingerID)
      ! Remove the mass of all child nodes.
      nodeChild => node%firstChild
      do while (associated(nodeChild))
         basicChild         =>  nodeChild %basic                    (                       )
         if (basic%floatRank0MetaPropertyGet(self%massBertschingerID) == 0.0d0) call self%nodeInitialize(nodeChild)
         nodeMassUnresolved =  +           nodeMassUnresolved                                 &
              &                -basicChild%floatRank0MetaPropertyGet(self%massBertschingerID)
         nodeChild          =>  nodeChild %sibling
      end do
      return
    end function nodeMassUnresolved

  end subroutine bertschingerMassNodeInitialize

  subroutine bertschingerMassDifferentialEvolutionAnalytics(self,node)
    !!{
    Mark analytically-solvable properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorBertschingerMass), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    class(nodeComponentBasic          ), pointer       :: basic

    basic => node%basic()
    call basic%floatRank0MetaPropertyAnalytic(self%massBertschingerID)
    return
  end subroutine bertschingerMassDifferentialEvolutionAnalytics

  subroutine bertschingerMassDifferentialEvolutionSolveAnalytics(self,node,time)
    !!{
    Evolve ``\gls{dmou}'' mass at a constant rate, to achieve linear interpolation in time.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (nodeOperatorBertschingerMass), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: time
    class           (nodeComponentBasic          ), pointer       :: basic            , basicParent
    double precision                                              :: massTarget       , timeTarget , &
         &                                                           massRateAccretion
    
    basic             => node %basic                    (                                )
    massRateAccretion =  basic%floatRank0MetaPropertyGet(self%accretionRateBertschingerID)
    if (massRateAccretion == 0.0d0) return
    basicParent => node       %parent%basic                    (                             )
    massTarget  =  basic             %floatRank0MetaPropertyGet(self%massBertschingerTargetID)
    timeTarget  =  basicParent       %time                     (                             )
    ! The mass is assumed to grow linearly with time.
    call basic%floatRank0MetaPropertySet(self%massBertschingerID,massTarget+massRateAccretion*(time-timeTarget))
    return
  end subroutine bertschingerMassDifferentialEvolutionSolveAnalytics

  subroutine bertschingerMassNodePromote(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent.
    !!}
    use :: Error             , only : Error_Report
    use :: Galacticus_Nodes  , only : nodeComponentBasic
    use :: ISO_Varying_String, only : var_str           , varying_string, operator(//)
    use :: String_Handling   , only : operator(//)
    implicit none
    class    (nodeOperatorBertschingerMass), intent(inout) :: self
    type     (treeNode                    ), intent(inout) :: node
    type     (treeNode                    ), pointer       :: nodeParent
    class    (nodeComponentBasic          ), pointer       :: basicParent, basic
    type     (varying_string              )                :: message
    character(len=12                      )                :: label
    !$GLC attributes unused :: self

    nodeParent  => node      %parent
    basic       => node      %basic ()
    basicParent => nodeParent%basic ()
    ! Ensure the two halos exist at the same time.
    if (basic%time() /= basicParent%time()) then
       message=var_str("node [")//node%index()//"] has not been evolved to its parent ["//nodeParent%index()//"]"//char(10)
       write (label,'(f12.6)') basic%time()
       message=message//"    node is at time: "//label//" Gyr"//char(10)
       write (label,'(f12.6)') basicParent%time()
       message=message//"  parent is at time: "//label//" Gyr"
       call Error_Report(message//{introspection:location})
    end if
    ! Adjust the mass, and accretion rate to that of the parent node.
    call basic%floatRank0MetaPropertySet(self%         massBertschingerID,basicParent%floatRank0MetaPropertyGet(self%         massBertschingerID))
    call basic%floatRank0MetaPropertySet(self%   massBertschingerTargetID,basicParent%floatRank0MetaPropertyGet(self%   massBertschingerTargetID))
    call basic%floatRank0MetaPropertySet(self%accretionRateBertschingerID,basicParent%floatRank0MetaPropertyGet(self%accretionRateBertschingerID))
    return
  end subroutine bertschingerMassNodePromote

  subroutine bertschingerMassNodesMerge(self,node)
    !!{
    Act on a merger between nodes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorBertschingerMass), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    class(nodeComponentBasic          ), pointer       :: basic

    ! Shut down mass accretion onto the halo now that it is a satellite.
    basic => node%basic()
    call basic%floatRank0MetaPropertySet(self%accretionRateBertschingerID,0.0d0)
    return
  end subroutine bertschingerMassNodesMerge
