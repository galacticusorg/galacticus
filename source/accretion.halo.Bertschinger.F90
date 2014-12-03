!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of accretion from the \gls{igm} onto halos using simple truncation to
  !% mimic the effects of reionization, and the Bertschinger mass to define available mass.

  !# <accretionHalo name="accretionHaloBertschinger">
  !#  <description>Accretion onto halos using simple truncation to mimic the effects of reionization, and the Bertschinger mass to define available mass.</description>
  !# </accretionHalo>

  type, extends(accretionHaloSimple) :: accretionHaloBertschinger
     !% A halo accretion class using simple truncation to mimic the effects of reionization, and the Bertschinger mass to define
     !% available mass.
     private
   contains
     procedure :: velocityScale      => bertschingerVelocityScale
     procedure :: accretionRateTotal => bertschingerAccretionRateTotal
     procedure :: massTotal          => bertschingerMassTotal
  end type accretionHaloBertschinger

  interface accretionHaloBertschinger
     !% Constructors for the {\tt bertschinger} halo accretion class.
     module procedure bertschingerConstructor
     module procedure bertschingerDefaultConstructor
  end interface accretionHaloBertschinger

  interface assignment(=)
     module procedure bertschingerFromSimple
  end interface assignment(=)

  ! Initialization state.
  logical :: bertschingerInitialized=.false.

contains
  
  subroutine bertschingerFromSimple(bertschinger,simple)
    !% Assign a {\tt simple} halo accretion object to a {\tt bertschinger} halo accretion object.
    implicit none
    type(accretionHaloBertschinger), intent(inout) :: bertschinger
    type(accretionHaloSimple      ), intent(in   ) :: simple
    
    bertschinger%reionizationSuppressionTime    =simple%reionizationSuppressionTime
    bertschinger%reionizationSuppressionVelocity=simple%reionizationSuppressionVelocity
    bertschinger%negativeAccretionAllowed       =simple%negativeAccretionAllowed
    bertschinger%accreteNewGrowthOnly           =simple%accreteNewGrowthOnly
    bertschinger%radiation                      =simple%radiation
    return
  end subroutine bertschingerFromSimple
  
  subroutine bertschingerInitialize()
    !% Initialize the {\tt bertschinger} halo accretion class.
    use Array_Utilities
    use Galacticus_Nodes
    use Galacticus_Error
    implicit none
    
    if (.not.bertschingerInitialized) then
       !$omp critical(accretionHaloBertschingerInitialize)
       if (.not.bertschingerInitialized) then
          if     (                                                                                                                     &
               &  .not.(                                                                                                               &
               &         defaultBasicComponent%         massBertschingerIsGettable()                                                   &
               &        .and.                                                                                                          &
               &         defaultBasicComponent%accretionRateBertschingerIsGettable()                                                   &
               &       )                                                                                                               &
               & )                                                                                                                     &
               & call Galacticus_Error_Report                                                                                          &
               &   (                                                                                                                   &
               &    'bertschingerInitialize'                                                                                         , &
               &    'the "massBertschinger" and "accretionRateBertschinger" properties of the basic component must be gettable.'    // &
               &    Galacticus_Component_List(                                                                                         &
               &                              'basic'                                                                                , &
               &                                defaultBasicComponent%         massBertschingerAttributeMatch(requireGettable=.true.)  &
               &                               .intersection.                                                                          &
               &                                defaultBasicComponent%accretionRateBertschingerAttributeMatch(requireGettable=.true.)  &
               &                             )                                                                                         &
               &   )
          bertschingerInitialized=.true.
       end if
       !$omp end critical(accretionHaloBertschingerInitialize)
    end if
    return
  end subroutine bertschingerInitialize

  function bertschingerDefaultConstructor()
    !% Default constructor for the {\tt bertschinger} halo accretion class.
    implicit none
    type(accretionHaloBertschinger), target :: bertschingerDefaultConstructor

    call bertschingerInitialize()
    bertschingerDefaultConstructor=accretionHaloSimple()
    return
  end function bertschingerDefaultConstructor
       
  function bertschingerConstructor(reionizationSuppressionTime,reionizationSuppressionVelocity,negativeAccretionAllowed,accreteNewGrowthOnly)
    !% Default constructor for the {\tt bertschinger} halo accretion class.
    implicit none
    type            (accretionHaloBertschinger), target        :: bertschingerConstructor
    double precision                           , intent(in   ) :: reionizationSuppressionTime, reionizationSuppressionVelocity
    logical                                    , intent(in   ) :: negativeAccretionAllowed   , accreteNewGrowthOnly

    call bertschingerInitialize()
    bertschingerConstructor=accretionHaloSimple(reionizationSuppressionTime,reionizationSuppressionVelocity,negativeAccretionAllowed,accreteNewGrowthOnly)
    return
  end function bertschingerConstructor

  double precision function bertschingerVelocityScale(self,node)
    !% Returns the velocity scale to use for {\tt node}. Use the virial velocity.
    use Galacticus_Nodes
    use Dark_Matter_Profiles
    implicit none
    class(accretionHaloBertschinger), intent(inout)          :: self
    type (treeNode                 ), intent(inout), pointer :: node
    class(darkMatterProfileClass   )               , pointer :: darkMatterProfile_

    darkMatterProfile_        => darkMatterProfile                         (    )
    bertschingerVelocityScale =  darkMatterProfile_%circularVelocityMaximum(node)
    return
  end function bertschingerVelocityScale

  double precision function bertschingerAccretionRateTotal(self,node)
    !% Returns the velocity scale to use for {\tt node}. Use the virial velocity.
    use Galacticus_Nodes
    implicit none
    class(accretionHaloBertschinger), intent(inout)          :: self
    type (treeNode                 ), intent(inout), pointer :: node
    class(nodeComponentBasic       )               , pointer :: basic

    basic                          => node %basic                    ()
    bertschingerAccretionRateTotal =  basic%accretionRateBertschinger()
    return
  end function bertschingerAccretionRateTotal

  double precision function bertschingerMassTotal(self,node)
    !% Returns the velocity scale to use for {\tt node}. Use the virial velocity.
    use Galacticus_Nodes
    implicit none
    class(accretionHaloBertschinger), intent(inout)          :: self
    type (treeNode                 ), intent(inout), pointer :: node
    class(nodeComponentBasic       )               , pointer :: basic

    basic                 => node %basic           ()
    bertschingerMassTotal =  basic%massBertschinger()
    return
  end function bertschingerMassTotal
