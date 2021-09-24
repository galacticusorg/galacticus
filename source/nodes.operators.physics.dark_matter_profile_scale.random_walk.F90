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

  !!{
  Implements a node operator class that assigns dark matter profile scale radii using a Weiner process random walk around the mean
  expectation, and then to be interpolated linearly between child and parent nodes.
  !!}
  
  use :: Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadiusClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass
  use :: Galacticus_Nodes          , only : nodeComponentDarkMatterProfile
  
  !![
  <nodeOperator name="nodeOperatorDarkMatterProfileScaleRandomWalk">
    <description>
      A node operator class that assigns dark matter profile scale radii using a Weiner process random walk around the mean
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
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDarkMatterProfileScaleRandomWalk
     !!{
     A node operator class that assigns dark matter profile scale radii using a Weiner process random walk around the mean
     expectation, and then to be interpolated linearly between child and parent nodes.
     !!}
     private
     class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_ => null()
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_         => null()
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_          => null()
     double precision                                             :: energyVarianceSpecific
   contains
     final     ::                          darkMatterProfileScaleRandomWalkConstructorDestructor
     procedure :: nodeTreeInitialize    => darkMatterProfileScaleRandomWalkNodeTreeInitialize
     procedure :: nodeInitialize        => darkMatterProfileScaleRandomWalkNodeInitialize
     procedure :: differentialEvolution => darkMatterProfileScaleRandomWalkDifferentialEvolution
  end type nodeOperatorDarkMatterProfileScaleRandomWalk
  
  interface nodeOperatorDarkMatterProfileScaleRandomWalk
     !!{
     Constructors for the {\normalfont \ttfamily darkMatterProfileScaleRandomWalk} node operator class.
     !!}
     module procedure darkMatterProfileScaleRandomWalkConstructorParameters
     module procedure darkMatterProfileScaleRandomWalkConstructorInternal
  end interface nodeOperatorDarkMatterProfileScaleRandomWalk

  ! Sub-module-scope variables used in root finding.
  double precision                                                        :: energyTarget_
  class           (nodeOperatorDarkMatterProfileScaleRandomWalk), pointer :: self_
  class           (nodeComponentDarkMatterProfile              ), pointer :: darkMatterProfileWork_
  type            (treeNode                                    ), pointer :: nodeWork_
  !$omp threadprivate(energyTarget_,self_,nodeWork_,darkMatterProfileWork_)

contains
  
  function darkMatterProfileScaleRandomWalkConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily darkMatterProfileScaleRandomWalk} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorDarkMatterProfileScaleRandomWalk)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    class           (darkMatterProfileScaleRadiusClass           ), pointer       :: darkMatterProfileScaleRadius_
    class           (darkMatterProfileDMOClass                   ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass                    ), pointer       :: darkMatterHaloScale_
    double precision                                                              :: energyVarianceSpecific

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
    self=nodeOperatorDarkMatterProfileScaleRandomWalk(energyVarianceSpecific,darkMatterProfileScaleRadius_,darkMatterProfileDMO_,darkMatterHaloScale_)
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
    Constructor for the {\normalfont \ttfamily darkMatterProfileScaleRandomWalk} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorDarkMatterProfileScaleRandomWalk)                        :: self
    class           (darkMatterProfileScaleRadiusClass           ), intent(in   ), target :: darkMatterProfileScaleRadius_
    class           (darkMatterProfileDMOClass                   ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass                    ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                              , intent(in   )         :: energyVarianceSpecific
    !![
    <constructorAssign variables="energyVarianceSpecific, *darkMatterProfileScaleRadius_, *darkMatterProfileDMO_, *darkMatterHaloScale_"/>
    !!]

    return
  end function darkMatterProfileScaleRandomWalkConstructorInternal

  subroutine darkMatterProfileScaleRandomWalkConstructorDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily darkMatterProfileScaleRandomWalk} dark matter halo profile scale radius class.
    !!}
    implicit none
    type(nodeOperatorDarkMatterProfileScaleRandomWalk), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    <objectDestructor name="self%darkMatterProfileDMO_"        />
    <objectDestructor name="self%darkMatterHaloScale_"         />
    !!]
    return
  end subroutine darkMatterProfileScaleRandomWalkConstructorDestructor

  subroutine darkMatterProfileScaleRandomWalkNodeTreeInitialize(self,node)
    !!{
    Initialize dark matter profile scale radii.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic
    use :: Merger_Tree_Walkers             , only : mergerTreeWalkerAllNodes
    use :: Root_Finder                     , only : rootFinder                     , rangeExpandMultiplicative, rangeExpandSignExpectPositive, rangeExpandSignExpectNegative    
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
   implicit none
    class           (nodeOperatorDarkMatterProfileScaleRandomWalk), intent(inout), target  :: self
    type            (treeNode                                    ), intent(inout), target  :: node
    type            (treeNode                                    )               , pointer :: nodeWork
    class           (nodeComponentDarkMatterProfile              )               , pointer :: darkMatterProfile  , darkMatterProfileWork
    class           (nodeComponentBasic                          )               , pointer :: basicWork
    type            (mergerTreeWalkerAllNodes                    )                         :: treeWalker
    type            (rootFinder                                  )                         :: finder
    double precision                                                                       :: radiusScale        , radiusScaleOriginal  , &
         &                                                                                    energyScalePrevious, energyScaleCurrent   , &
         &                                                                                    energyPerturbation

    ! Initialize the scale radius, if necessary.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    if (darkMatterProfile%scale() < 0.0d0) then
       ! Perform our own depth-first tree walk to set scales in all nodes of the tree. This is necessary as we require access
       ! to the parent scale to set scale growth rates, but must initialize scales in a strictly depth-first manner as some
       ! algorithms rely on knowing the progenitor structure of the tree to compute scale radii.
       treeWalker=mergerTreeWalkerAllNodes(node%hostTree,spanForest=.true.)
       do while (treeWalker%next(nodeWork))
          ! Get the scale radius - this will initialize the radius if necessary.
          darkMatterProfileWork => nodeWork%darkMatterProfile(autoCreate=.true.)
          call darkMatterProfileWork%scaleSet(self%darkMatterProfileScaleRadius_%radius(nodeWork))
       end do
       ! Repeat the tree walk, and for each leaf node walk up the branch applying random perturbations in energy.
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
       self_      => self
       treeWalker =  mergerTreeWalkerAllNodes(node%hostTree,spanForest=.true.)
       do while (treeWalker%next(nodeWork))
          ! Skip non-leaf nodes.
          if (associated(nodeWork%firstChild)) cycle
          ! Walk up the branch, applying the Weiner process.
          energyPerturbation =0.0d0
          energyScalePrevious=0.0d0
          do while (associated(nodeWork))
             basicWork          =>  nodeWork                       %basic       (        )
             energyScaleCurrent =  +gravitationalConstantGalacticus                           &
                  &                *basicWork                      %mass        (        )**2 &
                  &                /self     %darkMatterHaloScale_ %virialRadius(nodeWork)
             if (energyScaleCurrent > energyScalePrevious) then
                energyPerturbation =+energyPerturbation                                              &
                     &              +sqrt(                                                           &
                     &                    +self%energyVarianceSpecific                               &
                     &                    *(                                                         &
                     &                      +energyScaleCurrent **2                                  &
                     &                      -energyScalePrevious**2                                  &
                     &                     )                                                         &
                     &                   )                                                           &
                     &              *nodeWork%hostTree%randomNumberGenerator_%standardNormalSample()
                energyScalePrevious=+energyScaleCurrent
             end if
             ! Solve for the new scale radius that gives us this energy.
             nodeWork_              =>  nodeWork
             darkMatterProfileWork  =>  nodeWork                                   %darkMatterProfile (                                       )
             darkMatterProfileWork_ =>  darkMatterProfileWork
             radiusScaleOriginal    =   darkMatterProfileWork                      %scale             (                                       )
             energyTarget_          =  +self                 %darkMatterProfileDMO_%energy            (          nodeWork                     ) &
                  &                    +                                            energyPerturbation
             radiusScale            =   finder                                     %find              (rootGuess=darkMatterProfileWork%scale())
             call darkMatterProfileWork%scaleSet(radiusScale)
             ! Move to the next node on the branch.
             if (nodeWork%isPrimaryProgenitor()) then
                nodeWork => nodeWork%parent
             else
                nodeWork => null()
             end if
          end do
       end do
    end if
    return
  end subroutine darkMatterProfileScaleRandomWalkNodeTreeInitialize

  double precision function energyRoot(radiusScale)
    !!{
    Root function used in finding the scale radius corresponding to a given halo energy.
    !!}
    implicit none
    double precision, intent(in   ) :: radiusScale

    call darkMatterProfileWork_%scaleSet(radiusScale)
    energyRoot=+self_%darkMatterProfileDMO_%energy       (nodeWork_) &
         &     -                            energyTarget_
    return
  end function energyRoot
  
  subroutine darkMatterProfileScaleRandomWalkNodeInitialize(self,node)
    !!{
    Compute the rate of growth of dark matter profile scale radius assuming a constant growth rate.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorDarkMatterProfileScaleRandomWalk), intent(inout)          :: self
    type            (treeNode                                    ), intent(inout), target  :: node
    class           (nodeComponentBasic                          )               , pointer :: basic            , basicParent
    class           (nodeComponentDarkMatterProfile              )               , pointer :: darkMatterProfile, darkMatterProfileParent
    double precision                                                                       :: timeInterval

    ! Set the growth rate for the scale radius.
    call self%nodeTreeInitialize(node)
    darkMatterProfile => node%darkMatterProfile()
    if (node%isPrimaryProgenitor()) then
       ! Node is the primary progenitor, so compute the scale radius growth rate.
       basic        =>  node              %basic()
       basicParent  =>  node       %parent%basic()
       timeInterval =  +basicParent       %time () &
            &          -basic             %time ()
       if (timeInterval > 0.0d0) then
          darkMatterProfileParent => node%parent%darkMatterProfile()
          call darkMatterProfile%scaleGrowthRateSet(                                         &
               &                                    (                                        &
               &                                     +darkMatterProfileParent%scale       () &
               &                                     -darkMatterProfile      %scale       () &
               &                                    )                                        &
               &                                    /                         timeInterval   &
               &                                   )
       else
          ! Time interval is non-positive - assume zero growth rate.
          call darkMatterProfile%scaleGrowthRateSet(0.0d0)
       end if
    else
       ! Node is a non-primary progenitor - assume zero growth rate.
       call    darkMatterProfile%scaleGrowthRateSet(0.0d0)
    end if
    return
  end subroutine darkMatterProfileScaleRandomWalkNodeInitialize
  
  subroutine darkMatterProfileScaleRandomWalkDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Evolve dark matter profile scale radius at a constant rate, to achieve linear interpolation in time.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, propertyTypeInactive
    implicit none
    class    (nodeOperatorDarkMatterProfileScaleRandomWalk), intent(inout), target  :: self
    type     (treeNode                                    ), intent(inout)          :: node
    logical                                                , intent(inout)          :: interrupt
    procedure(interruptTask                               ), intent(inout), pointer :: functionInterrupt
    integer                                                , intent(in   )          :: propertyType
    class    (nodeComponentDarkMatterProfile              )               , pointer :: darkMatterProfile
    !$GLC attributes unused :: interrupt, functionInterrupt
    
    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    call darkMatterProfile%scaleRate(darkMatterProfile%scaleGrowthRate())
    return
  end subroutine darkMatterProfileScaleRandomWalkDifferentialEvolution
