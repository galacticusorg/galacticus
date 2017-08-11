!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

  !% An implementation of dark matter halo profile concentrations using the \cite{ludlow_mass-concentration-redshift_2016}
  !% algorithm.

  use Cosmology_Functions
  use Cosmology_Parameters
  use Dark_Matter_Halo_Scales
  use Dark_Matter_Profiles
  
  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationLudlow2016">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{ludlow_mass-concentration-redshift_2016}.</description>
  !# </darkMatterProfileConcentration>
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationLudlow2016
     !% A dark matter halo profile concentration class implementing the algorithm of
     !% \cite{ludlow_mass-concentration-redshift_2016}.
     private
     class           (cosmologyFunctionsClass            ), pointer :: cosmologyFunctions_             => null()
     class           (cosmologyParametersClass           ), pointer :: cosmologyParameters_            => null()
     class           (darkMatterProfileConcentrationClass), pointer :: darkMatterProfileConcentration_ => null()
     class           (darkMatterHaloScaleClass           ), pointer :: darkMatterHaloScale_            => null()
     class           (darkMatterProfileClass             ), pointer :: darkMatterProfile_              => null()
     double precision                                               :: C                                        , f              , &
          &                                                            timeFormationSeekDelta                   , densityContrast
   contains
     final     ::                                ludlow2016Destructor
     procedure :: concentration               => ludlow2016Concentration
     procedure :: densityContrastDefinition   => ludlow2016DensityContrastDefinition
     procedure :: darkMatterProfileDefinition => ludlow2016DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationLudlow2016

  interface darkMatterProfileConcentrationLudlow2016
     !% Constructors for the {\normalfont \ttfamily ludlow2016} dark matter halo profile concentration class.
     module procedure ludlow2016ConstructorParameters
     module procedure ludlow2016ConstructorInternal
  end interface darkMatterProfileConcentrationLudlow2016

  ! Module-scope variables used in root finding.
  class           (darkMatterProfileConcentrationLudlow2016), pointer :: ludlow2016Self
  type            (treeNode                                ), pointer :: ludlow2016Node
  double precision                                                    :: ludlow2016MassHaloCharacteristic, ludlow2016MassLimit
  !$omp threadprivate(ludlow2016Self,ludlow2016Node,ludlow2016MassHaloCharacteristic,ludlow2016MassLimit)
  
contains

  function ludlow2016ConstructorParameters(parameters) result(self)
    !% Default constructor for the {\normalfont \ttfamily ludlow2016} dark matter halo profile concentration class.
    implicit none
    type            (darkMatterProfileConcentrationLudlow2016)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                 ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass                ), pointer       :: cosmologyParameters_     
    class           (darkMatterProfileConcentrationClass     ), pointer       :: darkMatterProfileConcentration_
    class           (darkMatterHaloScaleClass                ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileClass                  ), pointer       :: darkMatterProfile_
    double precision                                                          :: C                              , f, &
         &                                                                       timeFormationSeekDelta

    ! Check and read parameters.
    !# <inputParameter>
    !#   <name>C</name>
    !#   <source>parameters</source>
    !#   <defaultValue>400.0d0</defaultValue>
    !#   <description>The parameter $C$ appearing in the halo concentration algorithm of \cite{ludlow_mass-concentration-redshift_2016}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>f</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.02d0</defaultValue>
    !#   <description>The parameter $f$ appearing in the halo concentration algorithm of \cite{ludlow_mass-concentration-redshift_2016}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>timeFormationSeekDelta</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The parameter $\Delta \log t$ by which the logarithm of the trial formation time is incremented when stepping through the formation history of a node to find the formation time. If set to zero (or a negative value) the cumulative mass histories of nodes are assumed to be monotonic functions of time, and the formation time is instead found by a root finding algorithm,</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"             name="cosmologyFunctions_"             source="parameters"/>
    !# <objectBuilder class="cosmologyParameters"            name="cosmologyParameters_"            source="parameters"/>
    !# <objectBuilder class="darkMatterProfileConcentration" name="darkMatterProfileConcentration_" source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"            name="darkMatterHaloScale_"            source="parameters"/>
    !# <objectBuilder class="darkMatterProfile"              name="darkMatterProfile_"              source="parameters"/>
    self=darkMatterProfileConcentrationLudlow2016(C,f,timeFormationSeekDelta,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileConcentration_,darkMatterHaloScale_,darkMatterProfile_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function ludlow2016ConstructorParameters
  
  function ludlow2016ConstructorInternal(C,f,timeFormationSeekDelta,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileConcentration_,darkMatterHaloScale_,darkMatterProfile_) result(self)
    !% Constructor for the {\normalfont \ttfamily ludlow2016} dark matter halo profile concentration class.
    use Galacticus_Error
    implicit none
    type            (darkMatterProfileConcentrationLudlow2016)                        :: self
    double precision                                          , intent(in   )         :: C                              , f, &
         &                                                                               timeFormationSeekDelta
    class           (cosmologyFunctionsClass                 ), intent(in   ), target :: cosmologyFunctions_
    class           (cosmologyParametersClass                ), intent(in   ), target :: cosmologyParameters_     
    class           (darkMatterProfileConcentrationClass     ), intent(in   ), target :: darkMatterProfileConcentration_
    class           (darkMatterHaloScaleClass                ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileClass                  ), intent(in   ), target :: darkMatterProfile_
    !# <constructorAssign variables="C, f, timeFormationSeekDelta, *cosmologyFunctions_, *cosmologyParameters_, *darkMatterProfileConcentration_, *darkMatterHaloScale_, *darkMatterProfile_"/>

    ! Find the density contrast as used to define masses by Ludlow et al. (2016).
    self%densityContrast=200.0d0/self%cosmologyParameters_%omegaMatter()
    return
  end function ludlow2016ConstructorInternal

  subroutine ludlow2016Destructor(self)
    !% Destructor for the {\normalfont \ttfamily ludlow2016} dark matter halo profile concentration class.
    implicit none
    type(darkMatterProfileConcentrationLudlow2016), intent(inout) :: self
    
    !# <objectDestructor name="self%cosmologyFunctions_"            />
    !# <objectDestructor name="self%cosmologyParameters_"           />
    !# <objectDestructor name="self%darkMatterProfileConcentration_"/>
    !# <objectDestructor name="self%darkMatterHaloScale_"           />
    !# <objectDestructor name="self%darkMatterProfile_"             />
    return
  end subroutine ludlow2016Destructor

  double precision function ludlow2016Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the
    !% \cite{ludlow_mass-concentration-redshift_2016} algorithm.
    use Root_Finder
    use Numerical_Constants_Math
    use Numerical_Comparison
    use Galacticus_Error
    use Galacticus_Display
    use Dark_Matter_Profile_Mass_Definitions
    implicit none
    class           (darkMatterProfileConcentrationLudlow2016), intent(inout), target  :: self
    type            (treeNode                                ), intent(inout), pointer :: node
    type            (treeNode                                )               , pointer :: nodeBranch
    class           (nodeComponentBasic                      )               , pointer :: basic                       , basicBranch
    class           (nodeComponentDarkMatterProfile          )               , pointer :: darkMatterProfile_
    double precision                                          , parameter              :: toleranceAbsolute     =0.0d0, toleranceRelative     =1.0d-4
    integer                                                   , parameter              :: iterationCountMaximum =100
    type            (rootFinder                              ), save                   :: finder
    !$omp threadprivate(finder) 
    double precision                                                                   :: massHaloCharacteristic      , timeBranchEarliest           , &
         &                                                                                densityMeanScaleRadius      , timeFormation                , &
         &                                                                                radiusScale                 , radiusScalePrevious          , &
         &                                                                                timeFormationTrial
    integer                                                                            :: iterationCount
    
    ! Get the dark matter profile component of the node.
    darkMatterProfile_ => node%darkMatterProfile()
    ! Set an initial guess to the scale radius using the fall-back concentration method.
    radiusScalePrevious=+self%darkMatterHaloScale_           %virialRadius     (node) &
         &              /self%darkMatterProfileConcentration_%concentrationMean(node)
    call darkMatterProfile_%scaleSet(radiusScalePrevious)
    ! For halos with no progenitors, simply keep the fall-back result.
    if (.not.associated(node%firstChild)) then
       radiusScale=radiusScalePrevious
    else
       ! Begin iteratively seeking a solution for the scale radius.
       iterationCount=0
       do while (iterationCount < iterationCountMaximum)
          iterationCount=iterationCount+1
          ! Compute the characteristic halo mass, M₋₂.
          basic                            =>  node                    %basic       (                               )
          massHaloCharacteristic           =  +self %darkMatterProfile_%enclosedMass(node,darkMatterProfile_%scale())
          ludlow2016Self                   =>  self
          ludlow2016Node                   =>  node
          ludlow2016MassHaloCharacteristic =  +massHaloCharacteristic
          ludlow2016MassLimit              =  +self%f                                                                                                                   &
               &                              *Dark_Matter_Profile_Mass_Definition(                                                                                     &
               &                                                                      node                                                                            , &
               &                                                                   +  self                    %densityContrast                                          &
               &                                                                   *(                                                                                   &
               &                                                                     +self%cosmologyFunctions_%hubbleParameterEpochal(time           =basic%time())     &
               &                                                                     /self%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=      1.0d0 )     &
               &                                                                    )                                                                              **2  &
               &                                                                   *  self%cosmologyFunctions_%expansionFactor       (                basic%time())**3  &
               &                                                                  )
          ! Find the earliest time in the branch.
          nodeBranch         => node
          timeBranchEarliest =  huge(0.0d0)
          do while (associated(nodeBranch))
             basicBranch        =>                        nodeBranch %basic     (    )
             timeBranchEarliest =  min(timeBranchEarliest,basicBranch%time      (    ))
             nodeBranch         =>                        nodeBranch %walkBranch(node)
          end do
          ! Test if the formation time is before the earliest time in the branch.
          if (ludlow2016FormationTimeRoot(timeBranchEarliest) > 0.0d0) then
             ! The characteristic halo mass is never resolved this branch - fall though to an alternative concentration calculation.
             radiusScale        =+self%darkMatterHaloScale_           %virialRadius (node) &
                  &              /self%darkMatterProfileConcentration_%concentration(node)
             ! Force convergence by setting the previous scale radius to that which we just set - since we're choosing a (possibly
             ! random) concentration from a distribution we no longer need to iterate to find a solution.
             radiusScalePrevious=+radiusScale
          else
             ! Find the time at which the mass in progenitors equals this characteristic halo mass.
             if (.not.finder%isInitialized()) then
                call finder%rootFunction(ludlow2016FormationTimeRoot        )
                call finder%tolerance   (toleranceAbsolute,toleranceRelative)
             end if
             timeFormation=finder%find(rootRange=[timeBranchEarliest,basic%time()])
             ! If requested, check for possible earlier formation times by simply stepping through trial times and finding the
             ! earliest at which the required mass threshold is reached. This is used for cases where the cumulative mass history
             ! is not monotonic.
             if (self%timeFormationSeekDelta > 0.0d0) then
                timeFormationTrial=timeBranchEarliest
                do while (timeFormationTrial < timeFormation)
                   if (ludlow2016FormationTimeRoot(timeFormationTrial) > 0.0d0) then
                      timeFormation=timeFormationTrial
                      exit
                   end if
                   timeFormationTrial=+timeFormationTrial                           &
                        &             *exp(self%timeFormationSeekDelta)
                end do
             end if
             ! Find the mean density with the scale radius, ⟨ρ₋₂⟩.
             densityMeanScaleRadius=+  self%C                                                                          &
                  &                 *  self%cosmologyParameters_%densityCritical       (                             ) &
                  &                 *(                                                                                 &
                  &                   +self%cosmologyFunctions_ %hubbleParameterEpochal(time           =timeFormation) &
                  &                   /self%cosmologyFunctions_ %hubbleParameterEpochal(expansionFactor=1.0d0        ) &
                  &                  )**2
             ! Find the corresponding scale radius.
             radiusScale=(                        &
                  &       +3.0d0                  &
                  &       *massHaloCharacteristic &
                  &       /densityMeanScaleRadius &
                  &       /4.0d0                  &
                  &       /Pi                     &
                  &      )**(1.0d0/3.0d0)
          end if
          ! Test for convergence.
          if (Values_Agree(radiusScale,radiusScalePrevious,relTol=1.0d-3)) exit
          ! Convergence was not attained - record current results and perform another iteration.
          call darkMatterProfile_%scaleSet(radiusScale)
          radiusScalePrevious=radiusScale
       end do
       ! Note that we do not check if convergence was actually reached. Due to the nature of the algorithm it is possible that the
       ! function for which we are seeking the root is discontinuous. As the scale radius changes, so does M₀ and therefore the
       ! mass limit f·M₀. As a result some progenitors will discontinuously pass/fail the mass limit check, making the root
       ! function discontinuous. Oscillating solutions can't be avoid in such situations, so we simply take the final iteration as
       ! our best estimate of the scale radius.
    end if
    ! Find the corresponding concentration.
    ludlow2016Concentration=+self%darkMatterHaloScale_%virialRadius(node) &
         &                  /                          radiusScale
    return
  end function ludlow2016Concentration

  double precision function ludlow2016FormationTimeRoot(timeFormation)
    !% Function used to find the formation time of a halo in the {\normalfont \ttfamily ludlow2016} concentration algorithm.
    use Dark_Matter_Profile_Mass_Definitions
    implicit none
    double precision                    , intent(in   ) :: timeFormation
    type            (treeNode          ), pointer       :: nodeBranch   , nodeChild        , &
         &                                                 nodeSibling
    class           (nodeComponentBasic), pointer       :: basicBranch  , basicChild       , &
         &                                                 basicSibling
    double precision                                    :: massBranch   , massAccretionRate, &
         &                                                 massSiblings , massProgenitor
    
    nodeBranch => ludlow2016Node
    massBranch =  0.0d0
    do while (associated(nodeBranch))
       basicBranch => nodeBranch%basic()
       if (associated(nodeBranch%firstChild).and.basicBranch%time() >= timeFormation) then
          nodeChild => nodeBranch%firstChild
          do while (associated(nodeChild))
             basicChild => nodeChild%basic()
             if (basicChild%time() < timeFormation) then
                if (nodeChild%isPrimaryProgenitor()) then
                   ! Interpolate in mass for primary progenitors.
                   massSiblings =  0.0d0
                   nodeSibling  => nodeChild
                   do while (associated(nodeSibling))
                      basicSibling =>  nodeSibling %basic  ()
                      massSiblings =  +massSiblings                                                                                                                              &
                           &          +Dark_Matter_Profile_Mass_Definition(                                                                                                      &
                           &                                                  nodeSibling                                                                                      , &
                           &                                               +  ludlow2016Self                    %densityContrast                                                 &
                           &                                               *(                                                                                                    &
                           &                                                 +ludlow2016Self%cosmologyFunctions_%hubbleParameterEpochal(time           =basicSibling%time())     &
                           &                                                 /ludlow2016Self%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=             1.0d0 )     &
                           &                                                )                                                                                               **2  &
                           &                                               *  ludlow2016Self%cosmologyFunctions_%expansionFactor       (                basicSibling%time())**3  &
                           &                                              )
                      nodeSibling  =>  nodeSibling %sibling
                   end do
                   massAccretionRate=+(                                                                                                                                          &
                        &              +Dark_Matter_Profile_Mass_Definition(                                                                                                     &
                        &                                                      nodeBranch                                                                                      , &
                        &                                                   +  ludlow2016Self                    %densityContrast                                                &
                        &                                                   *(                                                                                                   &
                        &                                                     +ludlow2016Self%cosmologyFunctions_%hubbleParameterEpochal(time           =basicBranch%time())     &
                        &                                                     /ludlow2016Self%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=            1.0d0 )     &
                        &                                                    )                                                                                              **2  &
                        &                                                   *  ludlow2016Self%cosmologyFunctions_%expansionFactor       (                basicBranch%time())**3  &
                        &                                                  )                                                                                                     &
                        &              -massSiblings                                                                                                                             &
                        &             )                                                                                                                                          &
                        &            /(                                                                                                                                          &
                        &              +basicBranch%time()                                                                                                                       &
                        &              -basicChild %time()                                                                                                                       &
                        &             )
                   massProgenitor   =+Dark_Matter_Profile_Mass_Definition(                                                                                                    &
                        &                                                    nodeChild                                                                                      , &
                        &                                                 +  ludlow2016Self                    %densityContrast                                               &
                        &                                                 *(                                                                                                  &
                        &                                                   +ludlow2016Self%cosmologyFunctions_%hubbleParameterEpochal(time           =basicChild%time())     &
                        &                                                   /ludlow2016Self%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=           1.0d0 )     &
                        &                                                  )                                                                                             **2  &
                        &                                                 *  ludlow2016Self%cosmologyFunctions_%expansionFactor       (                basicChild%time())**3  &
                        &                                                )                                                                                                    &
                        &            +massAccretionRate                                                                                                                       &
                        &            *(                                                                                                                                       &
                        &              +timeFormation                                                                                                                         &
                        &              -basicChild   %time()                                                                                                                  &
                        &             )
                   if (massProgenitor >= ludlow2016MassLimit) &
                        & massBranch=+massBranch              &
                        &            +massProgenitor
                else
                   ! No interpolation in mass for non-primary progenitors.
                   massProgenitor=Dark_Matter_Profile_Mass_Definition(                                                                                                    &
                        &                                                nodeChild                                                                                      , &
                        &                                             +  ludlow2016Self                    %densityContrast                                               &
                        &                                             *(                                                                                                  &
                        &                                               +ludlow2016Self%cosmologyFunctions_%hubbleParameterEpochal(time           =basicChild%time())     &
                        &                                               /ludlow2016Self%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=           1.0d0 )     &
                        &                                              )                                                                                             **2  &
                        &                                             *  ludlow2016Self%cosmologyFunctions_%expansionFactor       (                basicChild%time())**3  &
                        &                                            )
                   if (massProgenitor >= ludlow2016MassLimit) &
                        & massBranch=+massBranch              &
                        &            +massProgenitor 
                end if
             end if
             nodeChild => nodeChild%sibling
          end do
       else if (.not.associated(nodeBranch%firstChild).and.basicBranch%time() == timeFormation) then
          massProgenitor=Dark_Matter_Profile_Mass_Definition(                                                                                                     & 
               &                                                nodeBranch                                                                                      , &
               &                                             +  ludlow2016Self                    %densityContrast                                                &
               &                                             *(                                                                                                   &
               &                                               +ludlow2016Self%cosmologyFunctions_%hubbleParameterEpochal(time           =basicBranch%time())     &
               &                                               /ludlow2016Self%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=            1.0d0 )     &
               &                                              )                                                                                              **2  &
               &                                             *  ludlow2016Self%cosmologyFunctions_%expansionFactor       (                basicBranch%time())**3  &
               &                                            )
          if (massProgenitor >= ludlow2016MassLimit) &
               & massBranch=+massBranch              &
               &            +massProgenitor
       end if
       nodeBranch => nodeBranch%walkBranch(ludlow2016Node)
    end do
    ludlow2016FormationTimeRoot=+massBranch                       &
         &                      -ludlow2016MassHaloCharacteristic
   return
  end function ludlow2016FormationTimeRoot
  
  function ludlow2016DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of concentration in the
    !% \cite{ludlow_mass-concentration-redshift_2016} algorithm. While \cite{ludlow_mass-concentration-redshift_2016} used $\Delta
    !% = 200 \rho_{\rm crit}$ to define halos, their model actually predicts the scale radius, $r_{-2}$, rather than the
    !% concentration. Therefore, here we report that the \cite{ludlow_mass-concentration-redshift_2016} concentrations are defined
    !% using the model's own virial density contrast definition --- this ensures that the predicted scale radii are applied
    !% directly to model halos.
    use Virial_Density_Contrast
    implicit none
    class(virialDensityContrastClass              ), pointer       :: ludlow2016DensityContrastDefinition
    class(darkMatterProfileConcentrationLudlow2016), intent(inout) :: self
    !GCC$ attributes unused :: self
    
    ludlow2016DensityContrastDefinition => virialDensityContrast()
    return
  end function ludlow2016DensityContrastDefinition

  function ludlow2016DarkMatterProfileDefinition(self)
    !% Return a dark matter density profile object defining that used in the definition of concentration in the
    !% \cite{ludlow_mass-concentration-redshift_2016} algorithm.
    use Dark_Matter_Halo_Scales
    implicit none
    class  (darkMatterProfileClass                            ), pointer                     :: ludlow2016DarkMatterProfileDefinition
    class  (darkMatterProfileConcentrationLudlow2016          ), intent(inout)               :: self
    class  (darkMatterHaloScaleVirialDensityContrastDefinition), pointer                     :: darkMatterHaloScaleDefinition
    class  (darkMatterProfileClass                            ), allocatable  , target, save :: densityProfileDefinition
    logical                                                                           , save :: densityProfileDefinitionInitialized=.false.
    !$omp threadprivate(densityProfileDefinition,densityProfileDefinitionInitialized)
    
    if (.not.densityProfileDefinitionInitialized) then
       allocate(darkMatterProfileEinasto                           :: densityProfileDefinition     )
       allocate(darkMatterHaloScaleVirialDensityContrastDefinition :: darkMatterHaloScaleDefinition)
       select type (densityProfileDefinition)
       type is (darkMatterProfileEinasto)
          select type (darkMatterHaloScaleDefinition)
          type is (darkMatterHaloScaleVirialDensityContrastDefinition)
             darkMatterHaloScaleDefinition=darkMatterHaloScaleVirialDensityContrastDefinition(self%densityContrastDefinition())
             densityProfileDefinition     =darkMatterProfileEinasto                          (darkMatterHaloScaleDefinition   )
             call densityProfileDefinition%makeIndestructible()
          end select
       end select
       densityProfileDefinitionInitialized=.true.
    end if
    ludlow2016DarkMatterProfileDefinition => densityProfileDefinition
    return
  end function ludlow2016DarkMatterProfileDefinition
