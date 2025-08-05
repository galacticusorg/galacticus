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
  Implements a dark matter profile scale radius class using a Weiner process random walk around the mean expectation.
  !!}
  
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Galacticus_Nodes        , only : nodeComponentDarkMatterProfile
  
  !![
  <darkMatterProfileScaleRadius name="darkMatterProfileScaleRadiusRandomWalk">
    <description>
      A dark matter profile scale radius class that assigns dark matter profile scale radii using a Weiner process random walk around the mean
      expectation, and then interpolates linearly between child and parent nodes. For primary progenitor nodes $\dot{r}_\mathrm{s}
      = (r_{\mathrm{s},i+1}-r_{\mathrm{s},i})/(t_{i+1}-t_i)$, where $r_{\mathrm{s},i}$ is the scale radius of the dark matter
      profile of the node in the initialized tree, $r_{\mathrm{s},i+1}$ is the spin of its parent node, and $t_i$ and $t_{i+1}$
      are the corresponding times. For non-primary progenitors the rate of change is set to zero, i.e. $\dot{r}_\mathrm{s}=0$.

      The energy of each halo is set to the energy expected from the mean scale radius for halos of the mass and epoch, plus a
      perturbation given by:
      \begin{equation}
      \Delta E_\mathrm{i}(t_2) = \Delta E_\mathrm{i}(t_1) + \left[ \sigma^2 \left\{ E_\mathrm{v}^2(t_2) - E_\mathrm{v}^2(t_1) \right\} \right]^{1/2} N(0,1).
      \end{equation}
      where $J_\mathrm{v}(t) = M_\mathrm{v}(t) V_\mathrm{v}(t) R_\mathrm{v}(t)$ is the characteristic virial angular momentum,
      $M_\mathrm{v}(t)$, $V_\mathrm{v}(t)$, and $R_\mathrm{v}(t)$ are the virial mass, velocity, and radius respectively,
      $\sigma^2$ represents the variance in angular momentum per unit increase in $J_\mathrm{v}^2$, and $N(0,1)$ is a random
      variable distributed as a standard normal.
    </description>
  </darkMatterProfileScaleRadius>
  !!]
  type, extends(darkMatterProfileScaleRadiusClass) :: darkMatterProfileScaleRadiusRandomWalk
     !!{
     A dark matter profile scale radius class that assigns dark matter profile scale radii using a Weiner process random walk
     around the mean expectation.
     !!}
     private
     class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_ => null()
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_         => null()
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_          => null()
     double precision                                             :: energyVarianceSpecific
   contains
     final     ::           darkMatterProfileScaleRandomWalkDestructor
     procedure :: radius => darkMatterProfileScaleRandomWalkRadius
  end type darkMatterProfileScaleRadiusRandomWalk
  
  interface darkMatterProfileScaleRadiusRandomWalk
     !!{
     Constructors for the \refClass{darkMatterProfileScaleRadiusRandomWalk} node operator class.
     !!}
     module procedure darkMatterProfileScaleRandomWalkConstructorParameters
     module procedure darkMatterProfileScaleRandomWalkConstructorInternal
  end interface darkMatterProfileScaleRadiusRandomWalk

  ! Sub-module-scope variables used in root finding.
  double precision                                                  :: energyTarget_     , radiusVirial_
  class           (darkMatterProfileScaleRadiusRandomWalk), pointer :: self_
  class           (nodeComponentDarkMatterProfile        ), pointer :: darkMatterProfile_
  type            (treeNode                              ), pointer :: node_
  !$omp threadprivate(energyTarget_,self_,node_,darkMatterProfile_,radiusVirial_)

contains
  
  function darkMatterProfileScaleRandomWalkConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileScaleRadiusRandomWalk} dark matter profile scale radius class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (darkMatterProfileScaleRadiusRandomWalk)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (darkMatterProfileScaleRadiusClass     ), pointer       :: darkMatterProfileScaleRadius_
    class           (darkMatterProfileDMOClass             ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass              ), pointer       :: darkMatterHaloScale_
    double precision                                                        :: energyVarianceSpecific

    !![
    <inputParameter>
      <name>energyVarianceSpecific</name>
      <description>The variance in the difference in the energy of a halo per unit energy growth.</description>
      <source>parameters</source>
      <defaultValue>0.006d0</defaultValue>
    </inputParameter>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    !!]
    self=darkMatterProfileScaleRadiusRandomWalk(energyVarianceSpecific,darkMatterProfileScaleRadius_,darkMatterProfileDMO_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    <objectDestructor name="darkMatterProfileDMO_"        />
    <objectDestructor name="darkMatterHaloScale_"         />
    !!]
    return
  end function darkMatterProfileScaleRandomWalkConstructorParameters

  function darkMatterProfileScaleRandomWalkConstructorInternal(energyVarianceSpecific,darkMatterProfileScaleRadius_,darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileScaleRadiusRandomWalk} dark matter profile scale radius class.
    !!}
    implicit none
    type            (darkMatterProfileScaleRadiusRandomWalk)                        :: self
    class           (darkMatterProfileScaleRadiusClass     ), intent(in   ), target :: darkMatterProfileScaleRadius_
    class           (darkMatterProfileDMOClass             ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass              ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                        , intent(in   )         :: energyVarianceSpecific
    !![
    <constructorAssign variables="energyVarianceSpecific, *darkMatterProfileScaleRadius_, *darkMatterProfileDMO_, *darkMatterHaloScale_"/>
    !!]

    return
  end function darkMatterProfileScaleRandomWalkConstructorInternal

  subroutine darkMatterProfileScaleRandomWalkDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileScaleRadiusRandomWalk} dark matter halo profile scale radius class.
    !!}
    implicit none
    type(darkMatterProfileScaleRadiusRandomWalk), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    <objectDestructor name="self%darkMatterProfileDMO_"        />
    <objectDestructor name="self%darkMatterHaloScale_"         />
    !!]
    return
  end subroutine darkMatterProfileScaleRandomWalkDestructor

  double precision function darkMatterProfileScaleRandomWalkRadius(self,node) result(radiusScale)
    !!{
    Initialize dark matter profile scale radii.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic
    use :: Root_Finder                     , only : rootFinder                    , rangeExpandMultiplicative, rangeExpandSignExpectPositive, rangeExpandSignExpectNegative    
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Mass_Distributions              , only : massDistributionClass
    implicit none
    class           (darkMatterProfileScaleRadiusRandomWalk), intent(inout), target  :: self
    type            (treeNode                              ), intent(inout), target  :: node
    class           (nodeComponentDarkMatterProfile        )               , pointer :: darkMatterProfileChild
    class           (nodeComponentBasic                    )               , pointer :: basic                         , basicChild
    class           (massDistributionClass                 )               , pointer :: massDistribution_             , massDistributionChild_
    type            (rootFinder                            ), save                   :: finder
    logical                                                 , save                   :: finderConstructed     =.false.
    !$omp threadprivate(finder,finderConstructed)
    double precision                                                                 :: energyPerturbation            , radiusScaleOriginal   , &
         &                                                                              energyScaleChild              , energyScale           , &
         &                                                                              radiusVirialChild

    ! Set the scale radius to that given by the lower level scale radius class.
    radiusScale=self%darkMatterProfileScaleRadius_%radius(node)
    ! For nodes with a child, evaluate the perturbation to this radius.
    if (associated(node%firstChild)) then
       basic                  =>  node                                        %basic            (                                        )
       basicChild             =>  node                  %firstChild           %basic            (                                        )
       darkMatterProfileChild =>  node                  %firstChild           %darkMatterProfile(                                        )
       radiusVirial_          =   self                  %darkMatterHaloScale_ %radiusVirial     (node                                    )
       radiusVirialChild      =   self                  %darkMatterHaloScale_ %radiusVirial     (node%firstChild                         )
       energyScale            =  +gravitationalConstant_internal                                                                              &
            &                    *basic                                       %mass             (                                        )**2 &
            &                    /self                  %darkMatterHaloScale_ %radiusVirial     (node                                    )
       energyScaleChild       =  +gravitationalConstant_internal                                                                              &
            &                    *basicChild                                  %mass             (                                        )**2 &
            &                    /radiusVirialChild       
       massDistributionChild_ =>  self                  %darkMatterProfileDMO_%get              (node%firstChild                         )
       energyPerturbation     =  +massDistributionChild_                      %energy           (radiusVirialChild,massDistributionChild_)
       !![
       <objectDestructor name="massDistributionChild_"/>
       !!]
       radiusScaleOriginal    =  +darkMatterProfileChild                      %scale            (                                        )
       call darkMatterProfileChild%scaleSet(self%darkMatterProfileScaleRadius_%radius           (node%firstChild                         ))
       massDistributionChild_ =>  self                  %darkMatterProfileDMO_%get              (node%firstChild                         )
       energyPerturbation     =  +energyPerturbation                                                                                          &
            &                    -massDistributionChild_                      %energy           (radiusVirialChild,massDistributionChild_)
       !![
       <objectDestructor name="massDistributionChild_"/>
       !!]
       call darkMatterProfileChild%scaleSet(                                   radiusScaleOriginal                                        )
       if (energyScale > energyScaleChild) then
          energyPerturbation =+energyPerturbation                                          &
               &              +sqrt(                                                       &
               &                    +self%energyVarianceSpecific                           &
               &                    *(                                                     &
               &                      +energyScale     **2                                 &
               &                      -energyScaleChild**2                                 &
               &                     )                                                     &
               &                   )                                                       &
               &              *node%hostTree%randomNumberGenerator_%standardNormalSample()
       end if
       ! Solve for the new scale radius that gives us this energy.
       if (.not.finderConstructed) then 
          finder=rootFinder(                                                             &
               &            rootFunction                 =energyRoot                   , &
               &            toleranceAbsolute            =1.0d-6                       , &
               &            toleranceRelative            =1.0d-6                       , &
               &            rangeExpandUpward            =2.0d+0                       , &
               &            rangeExpandDownward          =0.5d+0                       , &
               &            rangeExpandType              =rangeExpandMultiplicative    , &
               &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &            rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
               &           )
          finderConstructed=.true.
       end if
       self_               =>  self
       node_               =>  node
       darkMatterProfile_  =>  node                                   %darkMatterProfile (autoCreate=.true.                    )
       radiusScaleOriginal =   darkMatterProfile_                     %scale             (                                     )
       call darkMatterProfile_%scaleSet(radiusScale        )
       massDistribution_   =>  self             %darkMatterProfileDMO_%get               (node                                 )
       energyTarget_       =  +massDistribution_                      %energy            (radiusVirial_,massDistribution_      ) &
            &                 +                                        energyPerturbation
       radiusScale         =   finder                                 %find              (rootGuess =darkMatterProfile_%scale())
       call darkMatterProfile_%scaleSet(radiusScaleOriginal)
       !![
       <objectDestructor name="massDistribution_"     />
       !!]
    end if
    return
  end function darkMatterProfileScaleRandomWalkRadius

  double precision function energyRoot(radiusScale)
    !!{
    Root function used in finding the scale radius corresponding to a given halo energy.
    !!}
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    double precision                       , intent(in   ) :: radiusScale
    class           (massDistributionClass), pointer       :: massDistribution_

    call darkMatterProfile_%scaleSet(radiusScale)
    massDistribution_ =>  self_            %darkMatterProfileDMO_%get          (node_                          ) 
    energyRoot        =  +massDistribution_                      %energy       (radiusVirial_,massDistribution_) &
         &               -                                        energyTarget_
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function energyRoot
