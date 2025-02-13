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
Contains a module which implements calculations of dark matter halo angular momentum.
!!}

module Dark_Matter_Halo_Spins
  !!{
  Implements calculations of dark matter halo angular momentum.
  !!}
  implicit none
  private
  public :: Dark_Matter_Halo_Angular_Momentum_Scale

  ! Record of whether the module has been initialized.
  logical :: propertiesAsserted=.false.

contains

  subroutine assertPropertiesGettable()
    !!{
    Assert that properties required for spin calculations are gettable.
    !!}
    use :: Error             , only : Component_List       , Error_Report
    use :: Galacticus_Nodes  , only : defaultBasicComponent
    use :: ISO_Varying_String, only : operator(//)
    implicit none

    if (.not.propertiesAsserted) then
       !$omp critical(darkMatterHaloSpinsAssertions)
       if (.not.propertiesAsserted) then
          if (.not.defaultBasicComponent%massIsGettable())                                               &
               & call Error_Report                                                                       &
               &      (                                                                                  &
               &       'mass property of basic component must be gettable.'//                            &
               &       Component_List(                                                                   &
               &                      'basic'                                                         ,  &
               &                      defaultBasicComponent%massAttributeMatch(requireGettable=.true.)   &
               &                     )                                                                // &
               &       {introspection:location}                                                          &
               &      )
          ! Record that the module is now initialized.
          propertiesAsserted=.true.
       end if
       !$omp end critical(darkMatterHaloSpinsAssertions)
    end if
    return
  end subroutine assertPropertiesGettable

  double precision function Dark_Matter_Halo_Angular_Momentum_Scale(node,darkMatterHaloScale_,darkMatterProfileDMO_,useBullockDefinition) result(angularMomentumScale)
    !!{
    Returns the characteristic angular momentum scale of {\normalfont \ttfamily node} (as used in spin definitions) based on its mass, and energy.
    !!}
    use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleClass
    use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMOClass
    use :: Galacticus_Nodes                , only : nodeComponentBasic            , treeNode
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Galactic_Structure_Options      , only : componentTypeDarkMatterOnly   , massTypeDark
    implicit none
    type   (treeNode                 ), intent(inout)           :: node
    class  (darkMatterHaloScaleClass ), intent(inout)           :: darkMatterHaloScale_
    class  (darkMatterProfileDMOClass), intent(inout), optional :: darkMatterProfileDMO_
    logical                           , intent(in   ), optional :: useBullockDefinition
    class  (nodeComponentBasic       ), pointer                 :: basic
    class  (massDistributionClass    ), pointer                 :: massDistribution_
   !![
    <optionalArgument name="useBullockDefinition" defaultsTo=".false." />
    !!]
    
    call assertPropertiesGettable()
    basic => node%basic()
    if (useBullockDefinition_) then
       ! Use the halo angular momentum scale used in the Bullock et al. (2001; http://adsabs.harvard.edu/abs/2001ApJ...555..240B)
       ! definition of halo spin.
       angularMomentumScale=+sqrt(2.0d0)                               &
            &               *basic               %mass          (    ) &
            &               *darkMatterHaloScale_%velocityVirial(node) &
            &               *darkMatterHaloScale_%radiusVirial  (node)
    else
       ! Use the halo angular momentum scale used in the Peebles (1971; http://adsabs.harvard.edu/abs/1971A%26A....11..377P)
       ! definition of halo spin.
       if (present(darkMatterProfileDMO_)) then
          massDistribution_    =>  darkMatterProfileDMO_%get             (node                                    )
       else
          massDistribution_    =>  node                 %massDistribution(componentTypeDarkMatterOnly,massTypeDark)
       end if
       angularMomentumScale =  +gravitationalConstant_internal                                                                        &
            &                  *         basic            %mass  (                                                         )  **2.5d0 &
            &                  /sqrt(abs(massDistribution_%energy(darkMatterHaloScale_%radiusVirial(node),massDistribution_)))
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    end if
    return
  end function Dark_Matter_Halo_Angular_Momentum_Scale

end module Dark_Matter_Halo_Spins
