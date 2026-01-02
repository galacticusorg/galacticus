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
Contains a module which implements calculations of dark matter halo formation times.
!!}

module Dark_Matter_Halo_Formation_Times
  !!{
  Implements calculations of dark matter halo formation times.
  !!}
  implicit none
  private
  public :: Dark_Matter_Halo_Formation_Time

contains

  double precision function Dark_Matter_Halo_Formation_Time(node,formationMassFraction,darkMatterHaloMassAccretionHistory_,nodeFormation)
    !!{
    Returns the time at which the main branch progenitor of {\normalfont \ttfamily node} first had a mass equal to {\normalfont \ttfamily
    formationMassFraction} of the current mass.
    !!}
    use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistoryClass
    use :: Error                                    , only : Error_Report
    use :: Galacticus_Nodes                         , only : nodeComponentBasic                     , treeNode
    implicit none
    type            (treeNode                               ), intent(inout), target            :: node
    double precision                                         , intent(in   )                    :: formationMassFraction
    class           (darkMatterHaloMassAccretionHistoryClass), intent(inout)                    :: darkMatterHaloMassAccretionHistory_
    type            (treeNode                               ), intent(inout), pointer, optional :: nodeFormation
    type            (treeNode                               ), pointer                          :: formationNode                      , workNode
    class           (nodeComponentBasic                     ), pointer                          :: basicParent                        , basic   , &
         &                                                                                         basicWork
    double precision                                                                            :: massNode                           , massWork, &
         &                                                                                         timeWork
    
    ! Validate the formation mass fraction.
    if     (                                                                                            &
         &   formationMassFraction <= 0.0d0                                                             &
         &  .or.                                                                                        &
         &   formationMassFraction >  1.0d0                                                             &
         & ) call Error_Report('`formationMassFraction` ∈ [0,1) is required'//{introspection:location}) 
    ! Get the mass of the starting node.
    basic    => node %basic()
    massNode =  basic%mass ()
    if (present(nodeFormation)) then
       ! We have an initial guess for the formation node. If necessary, walk back up the tree to find the appropriate starting
       ! node.
       workNode => nodeFormation
       do while (associated(workNode))
          basicWork     => workNode%basic()
          if (basicWork%mass() > formationMassFraction*massNode) exit
          workNode => workNode%parent
       end do
    else
       workNode => node
    end if  
    formationNode => null()
    do while (associated(workNode))
       formationNode => workNode
       basicWork     => workNode%basic()
       if (basicWork%mass() <= formationMassFraction*massNode) exit
       workNode => workNode%firstChild
    end do
    if (.not.associated(workNode)) then
       ! Find the formation time based on the mass accretion history.
       Dark_Matter_Halo_Formation_Time=darkMatterHaloMassAccretionHistory_%time(formationNode,formationMassFraction*massNode)
    else
       ! Interpolate to get the exact time of formation.
       basicParent => workNode %parent%basic()
       massWork    =  basicWork       %mass ()
       timeWork    =  basicWork       %time ()
       Dark_Matter_Halo_Formation_Time=                                 timeWork  &
            &                          +(basicParent%time()            -timeWork) &
            &                          *(formationMassFraction*massNode-massWork) &
            &                          /(basicParent%mass()            -massWork)
       if (present(nodeFormation)) nodeFormation => workNode
    end if
    return
  end function Dark_Matter_Halo_Formation_Time

end module Dark_Matter_Halo_Formation_Times
