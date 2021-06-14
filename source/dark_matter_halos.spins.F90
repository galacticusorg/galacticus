!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which implements calculations of dark matter halo angular momentum.

module Dark_Matter_Halo_Spins
  !% Implements calculations of dark matter halo angular momentum.
  implicit none
  private
  public :: Dark_Matter_Halo_Angular_Momentum, Dark_Matter_Halo_Spin, Dark_Matter_Halo_Angular_Momentum_Growth_Rate

  ! Record of whether the module has been initialized.
  logical :: propertiesAsserted=.false.

contains

  subroutine assertPropertiesGettable()
    !% Assert that properties required for spin calculations are gettable.
    use :: Galacticus_Error  , only : Galacticus_Component_List, Galacticus_Error_Report
    use :: Galacticus_Nodes  , only : defaultBasicComponent    , defaultSpinComponent
    use :: ISO_Varying_String, only : operator(//)
    implicit none

    if (.not.propertiesAsserted) then
       !$omp critical(darkMatterHaloSpinsAssertions)
       if (.not.propertiesAsserted) then
          ! Ensure that the spin property is available.
          if (.not.defaultSpinComponent%spinIsGettable())                                                           &
               & call Galacticus_Error_Report                                                                       &
               &      (                                                                                             &
               &       'spin property of spin component must be gettable.'//                                        &
               &       Galacticus_Component_List(                                                                   &
               &                                 'spin'                                                          ,  &
               &                                 defaultSpinComponent %spinAttributeMatch(requireGettable=.true.)   &
               &                                )                                                                // &
               &       {introspection:location}                                                                     &
               &      )
          if (.not.defaultBasicComponent%massIsGettable())                                                          &
               & call Galacticus_Error_Report                                                                       &
               &      (                                                                                             &
               &       'mass property of basic component must be gettable.'//                                       &
               &       Galacticus_Component_List(                                                                   &
               &                                 'basic'                                                         ,  &
               &                                 defaultBasicComponent%massAttributeMatch(requireGettable=.true.)   &
               &                                )                                                                // &
               &       {introspection:location}                                                                     &
               &      )
          ! Record that the module is now initialized.
          propertiesAsserted=.true.
       end if
       !$omp end critical(darkMatterHaloSpinsAssertions)
    end if
    return
  end subroutine assertPropertiesGettable

  double precision function Dark_Matter_Halo_Angular_Momentum(node,darkMatterProfileDMO_)
    !% Returns the total anuglar momentum of {\normalfont \ttfamily node} based on its mass, energy and spin parameter.
    use :: Dark_Matter_Profiles_DMO    , only : darkMatterProfileDMOClass
    use :: Galacticus_Nodes            , only : nodeComponentBasic             , nodeComponentSpin, treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    type (treeNode                 ), intent(inout) :: node
    class(darkMatterProfileDMOClass), intent(inout) :: darkMatterProfileDMO_
    class(nodeComponentBasic       ), pointer       :: basic
    class(nodeComponentSpin        ), pointer       :: spin

    call assertPropertiesGettable()
    basic => node%basic(                 )
    spin  => node%spin (autoCreate=.true.)
    Dark_Matter_Halo_Angular_Momentum=+gravitationalConstantGalacticus                      &
         &                            *         spin                 %spin  (    )          &
         &                            *         basic                %mass  (    )**2.5d0   &
         &                            /sqrt(abs(darkMatterProfileDMO_%energy(node)       ))
    return
  end function Dark_Matter_Halo_Angular_Momentum

  double precision function Dark_Matter_Halo_Angular_Momentum_Growth_Rate(node,darkMatterProfileDMO_)
    !% Returns the rate of change of the total anuglar momentum of {\normalfont \ttfamily node} based on its mass, energy and spin parameter.
    use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
    use :: Galacticus_Error        , only : Galacticus_Error_Report
    use :: Galacticus_Nodes        , only : nodeComponentBasic       , nodeComponentSpin, treeNode
    implicit none
    type            (treeNode                 ), intent(inout) :: node
    class           (darkMatterProfileDMOClass), intent(inout) :: darkMatterProfileDMO_
    class           (nodeComponentBasic       ), pointer       :: basic
    class           (nodeComponentSpin        ), pointer       :: spin
    double precision                                           :: rateFractional

    call assertPropertiesGettable()
    basic                                         =>  node%basic(                 )
    spin                                          =>  node%spin (autoCreate=.true.)
    Dark_Matter_Halo_Angular_Momentum_Growth_Rate =  +Dark_Matter_Halo_Angular_Momentum(node,darkMatterProfileDMO_)
    if (Dark_Matter_Halo_Angular_Momentum_Growth_Rate == 0.0d0) return
    rateFractional=0.0d0
    if     (spin                  %spin            (    ) >  0.0d0) then
       rateFractional=+rateFractional&
            &         +spin                 %spinGrowthRate   (    ) &
            &         /spin                 %spin             (    )
    else if (spin                 %spinGrowthRate  (    ) /= 0.0d0) then
       call Galacticus_Error_Report('spin is zero, but growth rate is non-zero'  //{introspection:location})
    end if
    if     (basic                 %mass            (    ) >  0.0d0) then
       rateFractional=+rateFractional        &
            &         +2.5d0                 &
            &         *basic                %accretionRate   (    ) &
            &         /basic                %mass            (    )
    else if (basic                %accretionRate   (    ) /= 0.0d0) then
       call Galacticus_Error_Report('mass is zero, but growth rate is non-zero'  //{introspection:location})
    end if
    if      (darkMatterProfileDMO_%energy          (node) /= 0.0d0) then
       rateFractional=+rateFractional                               &
            &         -0.5d0                                        &
            &         *darkMatterProfileDMO_%energyGrowthRate(node) &
            &         /darkMatterProfileDMO_%energy          (node)
    else if (darkMatterProfileDMO_%energyGrowthRate(node) /= 0.0d0) then
       call Galacticus_Error_Report('energy is zero, but growth rate is non-zero'//{introspection:location})
    end if
    Dark_Matter_Halo_Angular_Momentum_Growth_Rate=+Dark_Matter_Halo_Angular_Momentum_Growth_Rate &
         &                                        *rateFractional
    return
  end function Dark_Matter_Halo_Angular_Momentum_Growth_Rate

  double precision function Dark_Matter_Halo_Spin(node,angularMomentum,darkMatterProfileDMO_)
    !% Returns the spin of {\normalfont \ttfamily node} given its angular momentum.
    use :: Dark_Matter_Profiles_DMO    , only : darkMatterProfileDMOClass
    use :: Galacticus_Nodes            , only : nodeComponentBasic             , treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    type            (treeNode                 ), intent(inout) :: node
    double precision                           , intent(in   ) :: angularMomentum
    class           (darkMatterProfileDMOClass), intent(inout) :: darkMatterProfileDMO_
    class           (nodeComponentBasic       ), pointer       :: basic

    call assertPropertiesGettable()
    basic                 =>  node             %basic()
    Dark_Matter_Halo_Spin =  +angularMomentum                                &
         &                   /gravitationalConstantGalacticus                &
         &                   *abs(darkMatterProfileDMO_%energy(node))**0.5d0 &
         &                   /    basic                %mass  (    ) **2.5d0
    return
  end function Dark_Matter_Halo_Spin

end module Dark_Matter_Halo_Spins
