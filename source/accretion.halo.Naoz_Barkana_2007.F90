!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!+ Contributions to this file made by: Daniel McAndrew.

  !% An implementation of accretion from the \gls{igm} onto halos using filtering mass of the \gls{igm}
  !% calculated from an equation from \cite{naoz_formation_2007}.

  !# <accretionHalo name="accretionHaloNaozBarkana2007">
  !#  <description>Accretion onto halos using filtering mass of the \gls{igm} calculated from an equation from \cite{naoz_formation_2007}.</description>
  !# </accretionHalo>

  type, extends(accretionHaloSimple) :: accretionHaloNaozBarkana2007
     !% A halo accretion class using filtering mass of the \gls{igm} calculated from an equation from \cite{naoz_formation_2007}.
     private
   contains
     procedure :: branchHasBaryons => naozBarkana2007BranchHasBaryons
     procedure :: failedFraction   => naozBarkana2007FailedFraction
  end type accretionHaloNaozBarkana2007

  interface accretionHaloNaozBarkana2007
     !% Constructors for the {\normalfont \ttfamily naozBarkana2007} halo accretion class.
     module procedure naozBarkana2007Constructor
     module procedure naozBarkana2007DefaultConstructor
  end interface accretionHaloNaozBarkana2007

  interface assignment(=)
     module procedure naozBarkana2007FromSimple
  end interface assignment(=)

contains

  subroutine naozBarkana2007FromSimple(naozBarkana2007,simple)
    !% Assign a {\normalfont \ttfamily simple} halo accretion object to a {\normalfont \ttfamily naozBarkana2007} halo accretion object.
    implicit none
    type(accretionHaloNaozBarkana2007), intent(inout) :: naozBarkana2007
    type(accretionHaloSimple         ), intent(in   ) :: simple
    
    naozBarkana2007%reionizationSuppressionTime    =simple%reionizationSuppressionTime
    naozBarkana2007%reionizationSuppressionVelocity=simple%reionizationSuppressionVelocity
    naozBarkana2007%negativeAccretionAllowed       =simple%negativeAccretionAllowed
    naozBarkana2007%accreteNewGrowthOnly           =simple%accreteNewGrowthOnly
    naozBarkana2007%radiation                      =simple%radiation
    return
  end subroutine naozBarkana2007FromSimple
  
  function naozBarkana2007DefaultConstructor()
    !% Default constructor for the {\normalfont \ttfamily naozBarkana2007} halo accretion class.
    implicit none
    type(accretionHaloNaozBarkana2007), target :: naozBarkana2007DefaultConstructor

    naozBarkana2007DefaultConstructor=accretionHaloSimple()
    return
  end function naozBarkana2007DefaultConstructor
       
  function naozBarkana2007Constructor(reionizationSuppressionTime,reionizationSuppressionVelocity,negativeAccretionAllowed,accreteNewGrowthOnly)
    !% Default constructor for the {\normalfont \ttfamily naozBarkana2007} halo accretion class.
    implicit none
    type            (accretionHaloNaozBarkana2007), target        :: naozBarkana2007Constructor
    double precision                              , intent(in   ) :: reionizationSuppressionTime, reionizationSuppressionVelocity
    logical                                                       :: negativeAccretionAllowed   , accreteNewGrowthOnly
    
    naozBarkana2007Constructor=accretionHaloSimple(reionizationSuppressionTime,reionizationSuppressionVelocity,negativeAccretionAllowed,accreteNewGrowthOnly)
    return
  end function naozBarkana2007Constructor

  logical function naozBarkana2007BranchHasBaryons(self,node)
    !% Returns true if this branch can accrete any baryons.
    use Galacticus_Nodes
    implicit none
    class(accretionHaloNaozBarkana2007), intent(inout)          :: self
    type (treeNode                    ), intent(inout), pointer :: node
    
    naozBarkana2007BranchHasBaryons=.true.
    return
  end function naozBarkana2007BranchHasBaryons

  double precision function naozBarkana2007FailedFraction(self,node)
    !% Returns the velocity scale to use for {\normalfont \ttfamily node}. Use the virial velocity.
    use Galacticus_Nodes
    use Intergalactic_Medium_State
    use Cosmology_Parameters
    use Galacticus_Error
    implicit none
    class           (accretionHaloNaozBarkana2007 ), intent(inout)          :: self
    type            (treeNode                     ), intent(inout), pointer :: node
    class           (nodeComponentBasic           )               , pointer :: basic
    class           (intergalacticMediumStateClass)               , pointer :: igmState_
    class           (cosmologyParametersClass     )               , pointer :: cosmologyParameters_
    double precision                                                        :: massFiltering       , accretionFraction

    igmState_ => intergalacticMediumState()
    select type (igmState_)
    class is (intergalacticMediumStateInternal)
       basic         => node     %basic        (            )
       massFiltering =  igmState_%filteringMass(basic%time())
    class default
       call Galacticus_Error_Report('naozBarkana2007FailedFraction','requires [intergalacticMediumStateMethod]=internal')
    end select
    cosmologyParameters_ => cosmologyParameters()
    accretionFraction=+1.0d0                              &
         &            /(                                  &
         &              +1.0d0                            &
         &              +(                                &
         &                +2.0d0**(1.0d0/3.0d0)           &
         &                -1.0d0                          &
         &               )                                &
         &              *8.0d0                            &
         &              *massFiltering                    &
         &              /basic%mass()                     &
         &             )**3
    naozBarkana2007FailedFraction=1.0d0-accretionFraction
    return
  end function naozBarkana2007FailedFraction
