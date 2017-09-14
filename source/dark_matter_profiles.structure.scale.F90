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

!% Contains a module which implements calculations of dark matter profile scale radii from concentrations.

module Dark_Matter_Profile_Scales
  !% Implements calculations of dark matter profile scale radii from concentrations.
  use Galacticus_Nodes
  use Dark_Matter_Profiles
  use Dark_Matter_Profiles_Concentration
  use Dark_Matter_Halo_Scales
  use Cosmology_Functions
  use Virial_Density_Contrast
  private
  public :: Dark_Matter_Profile_Scale
  
  ! Initialization state.
  logical :: moduleInitialized=.false.

  ! Scale calculation option.
  logical :: darkMatterProfileScaleCorrectForConcentrationDefinition

contains

  double precision function Dark_Matter_Profile_Scale(node,meanConcentration,concentrationMethod)
    !% Compute the scale radius of the dark matter profile of {\normalfont \ttfamily node}.
    use Root_Finder
    use Galacticus_Calculations_Resets
    use Numerical_Constants_Math
    use Input_Parameters
    use Numerical_Comparison
    implicit none    
    type            (treeNode                                          ), pointer    , intent(inout)           :: node
    class           (darkMatterProfileConcentrationClass               ), target     , intent(inout), optional :: concentrationMethod
    logical                                                                          , intent(in   ), optional :: meanConcentration
    class           (virialDensityContrastClass                        ), pointer                              :: virialDensityContrastDefinition
    class           (virialDensityContrastClass                        ), pointer    , save                    :: virialDensityContrastDefinitionPrevious => null()
    !$omp threadprivate(virialDensityContrastDefinitionPrevious)
    class           (darkMatterProfileClass                            ), pointer                              :: darkMatterProfile_
    class           (darkMatterHaloScaleClass                          ), pointer                              :: darkMatterHaloScale_
    class           (nodeComponentBasic                                ), pointer                              :: basic
    type            (treeNode                                          ), pointer                              :: workNode
    class           (darkMatterProfileConcentrationClass               ), pointer                              :: darkMatterProfileConcentrationDefinition
    type            (darkMatterHaloScaleVirialDensityContrastDefinition), allocatable, save                    :: darkMatterHaloScaleDefinition
    !$omp threadprivate(darkMatterHaloScaleDefinition)
    class           (darkMatterProfileClass                            ), pointer                              :: darkMatterProfileDefinition
    class           (cosmologyFunctionsClass                           ), pointer                              :: cosmologyFunctions_
    class           (nodeComponentBasic                                ), pointer                              :: workBasic
    class           (nodeComponentDarkMatterProfile                    ), pointer                              :: workDarkMatterProfile
    class           (virialDensityContrastClass                        ), pointer                              :: virialDensityContrast_
    double precision                                                                 , save                    :: massRatioPrevious                       =  2.0d0
    !$omp threadprivate(massRatioPrevious)
    double precision                                                    , parameter                            :: massRatioBuffer                         =  1.1d0 , massRatioShrink=0.99d0
    type            (rootFinder                                        )                                       :: finder
    double precision                                                                                           :: mass                                             , massDefinition        , &
         &                                                                                                        concentration                                    , massRatio
    double precision                                                                                           :: concentrationOriginal

    ! Initialize as necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Dark_Matter_Profile_Scale_Initialize)
       if (.not.moduleInitialized) then
          ! Get parameter controlling scale calculation method.
          !# <inputParameter>
          !#   <name>darkMatterProfileScaleCorrectForConcentrationDefinition</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>If true, then when computing dark matter profile scale radii using concentrations, any difference between the current definition of halo scales
          !#     (i.e. typically virial density contrast definitions) and density profiles and those assumed in measuring the concentrations will be taken into account.
          !#     If false, the concentration is applied blindly.</description>
          !#   <source>globalParameters</source>
          !#   <type>string</type>
          !# </inputParameter>
          ! Record that module is initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical(Dark_Matter_Profile_Scale_Initialize)
    end if
    ! Get required objects.
    darkMatterHaloScale_                        => darkMatterHaloScale           ()
    if (present(concentrationMethod)) then
       darkMatterProfileConcentrationDefinition => concentrationMethod
    else
       darkMatterProfileConcentrationDefinition => darkMatterProfileConcentration()
    end if
    ! Find the original concentration.
    if (present(meanConcentration).and.meanConcentration) then
       concentrationOriginal=darkMatterProfileConcentrationDefinition%concentrationMean(node)
    else
       concentrationOriginal=darkMatterProfileConcentrationDefinition%concentration    (node)
    end if
    ! Determine if concentration must be corrected.
    if (darkMatterProfileScaleCorrectForConcentrationDefinition) then
       ! Get objects for the concentration definition.
       virialDensityContrast_          => virialDensityContrast                                             ()
       virialDensityContrastDefinition => darkMatterProfileConcentrationDefinition%densityContrastDefinition()
       ! Get the basic component of the supplied node and extract its mass.
       basic => node %basic()
       mass  =  basic%mass ()
       ! If there is no difference between the alt and non-alt virial density contrasts, then no correction need be made.
       if     (                                                                                                 &
            &  Values_Differ(                                                                                   &
            &                       virialDensityContrast_         %densityContrast(basic%mass(),basic%time()), &
            &                       virialDensityContrastDefinition%densityContrast(basic%mass(),basic%time()), &
            &                relTol=1.0d-6                                                                      &
            &               )                                                                                   &
            & ) then
          ! Get other required objects.
          cosmologyFunctions_         => cosmologyFunctions                                                  ()
          darkMatterProfile_          => darkMatterProfile                                                   ()
          darkMatterProfileDefinition => darkMatterProfileConcentrationDefinition%darkMatterProfileDefinition()
          ! Construct a new dark matter halo scale definition object only if the virial density contrast definition differs from that
          ! previously used. Otherwise, just re-use the saved dark matter halo scale definition.
          if (.not.associated(virialDensityContrastDefinition,virialDensityContrastDefinitionPrevious)) then
             if (allocated(darkMatterHaloScaleDefinition)) deallocate(darkMatterHaloScaleDefinition)
             allocate(darkMatterHaloScaleDefinition)
             darkMatterHaloScaleDefinition           =  darkMatterHaloScaleVirialDensityContrastDefinition(virialDensityContrastDefinition)
             virialDensityContrastDefinitionPrevious => virialDensityContrastDefinition
          end if
          ! Create a node and set the mass and time.
          workNode              => treeNode                  (                 )
          workBasic             => workNode%basic            (autoCreate=.true.)
          workDarkMatterProfile => workNode%darkMatterProfile(autoCreate=.true.)
          call workBasic            %timeSet            (basic%time())
          call workBasic            %timeLastIsolatedSet(basic%time())
          call workDarkMatterProfile%scaleIsLimitedSet  (.false.     )
          ! The finder is initialized each time as it is allocated on the stack - this allows this function to be called recursively.
          call finder               %tolerance          (                                                             &
               &                                         toleranceRelative            =1.0d-3                         &
               &                                        )
          call finder               %rangeExpand        (                                                             &
               &                                         rangeExpandUpward            =1.0d0*massRatioPrevious      , &
               &                                         rangeExpandDownward          =1.0d0/massRatioPrevious      , &
               &                                         rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                                         rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                                         rangeExpandType              =rangeExpandMultiplicative      &
               &                                        )
          call finder               %rootFunction       (                                                             &
               &                                                                       massRootFunction               &
               &                                        )
          massDefinition=finder%find(rootGuess=mass)
          ! Find the ratio of the recovered mass under the given definition to the input mass, defined to be always greater than
          ! unity. This will be used as the basis of the range expansion for the next solution.
          if (massDefinition > mass) then
             massRatio=1.0d0*(massDefinition/mass)
          else
             massRatio=1.0d0/(massDefinition/mass)
          end if
          ! Increase this mass ratio by a small factor.
          massRatio=massRatioBuffer*massRatio
          ! If this new mass ratio exceeds our previous mass ratio, update the previous mass ratio for use in the next
          ! solution. Otherwise, shrink the previous mass ratio by a small amount.
          if (massRatio > massRatioPrevious) then
             massRatioPrevious=massRatio
          else
             massRatioPrevious=massRatioShrink*massRatioPrevious
          end if
          ! Update the work node properties and computed concentration.
          call workBasic%massSet(massDefinition)
          call Galacticus_Calculations_Reset                 (workNode)
          call darkMatterProfileDefinition  %calculationReset(workNode)
          ! Find the concentration.
          if (present(meanConcentration).and.meanConcentration) then
             ! We are simply using the mean concentration-mass relation here.
             concentration=+darkMatterProfileConcentrationDefinition%concentrationMean(workNode)
          else
             ! In this case we need to allow for possible scatter in the concentration mass relation. Therefore, we take the original
             ! concentration (which may include some scatter away from the mean relation) and scale it by the ratio of the mean
             ! concentrations for the corrected and original nodes.
             concentration=+concentrationOriginal                                                &
                  &        *darkMatterProfileConcentrationDefinition%concentrationMean(workNode) &
                  &        /darkMatterProfileConcentrationDefinition%concentrationMean(    node)
          end if
          Dark_Matter_Profile_Scale=+darkMatterHaloScaleDefinition%virialRadius(workNode) &
               &                    /concentration
          call workNode%destroy()
          deallocate(workNode)
          ! Destroy objects as necessary.
          if (darkMatterProfileDefinition%isFinalizable()) deallocate(darkMatterProfileDefinition)
       else
          Dark_Matter_Profile_Scale=+darkMatterHaloScale_%virialRadius(node) &
               &                    /concentrationOriginal
       end if
    else
       Dark_Matter_Profile_Scale=+darkMatterHaloScale_%virialRadius(node) &
            &                    /concentrationOriginal
    end if
    ! Nullify the concentration definition so that it isn't automatically finalized.
    darkMatterProfileConcentrationDefinition => null()
    return

  contains
    
    double precision function massRootFunction(massDefinitionTrial)
      !% Root function used to find the mass of a halo corresponding to the definition used for a particular concentration class.
      use Galacticus_Calculations_Resets
      implicit none
      double precision, intent(in   ) :: massDefinitionTrial
      double precision                :: radiusOuterDefinition, concentrationDefinition, &
           &                             radiusCore           , massOuter              , &
           &                             radiusOuter          , densityOuter

      ! Set the mass of the worker node.
      call workBasic%massSet(massDefinitionTrial)
      call Galacticus_Calculations_Reset                 (workNode)
      call darkMatterHaloScaleDefinition%calculationReset(workNode)
      ! Get outer radius for this trial definition mass.
      radiusOuterDefinition=darkMatterHaloScaleDefinition%virialRadius(workNode)
      ! Get concentration for this a trial definition mass.
      if (present(meanConcentration).and.meanConcentration) then
         ! We are simply using the mean concentration-mass relation here.
         concentrationDefinition=darkMatterProfileConcentrationDefinition%concentrationMean(workNode)
      else
         ! In this case we need to allow for possible scatter in the concentration mass relation. Therefore, we take the original
         ! concentration (which may include some scatter away from the mean relation) and scale it by the ratio of the mean
         ! concentrations for the corrected and original nodes.
         concentrationDefinition=+concentrationOriginal                                               &
              &                 *darkMatterProfileConcentrationDefinition%concentrationMean(workNode) &
              &                 /darkMatterProfileConcentrationDefinition%concentrationMean(node    )
      end if
      ! Get core radius.      
      radiusCore             =radiusOuterDefinition/concentrationDefinition
      call workDarkMatterProfile%scaleSet(radiusCore)
      call Galacticus_Calculations_Reset                 (workNode)
      call darkMatterProfileDefinition  %calculationReset(workNode)
      ! Find the non-alt density.
      densityOuter=+cosmologyFunctions_   %matterDensityEpochal(                 workBasic%time()) &
           &       *virialDensityContrast_%densityContrast     (workBasic%mass(),workBasic%time())      
      ! Solve for radius which encloses required non-alt density.
      radiusOuter=darkMatterProfileDefinition%radiusEnclosingDensity(workNode,densityOuter)
      ! Get the mass within this radius.
      massOuter  =darkMatterProfileDefinition%enclosedMass          (workNode,radiusOuter)
      ! Return root function.
      massRootFunction=massOuter-mass
      return
    end function massRootFunction

  end function Dark_Matter_Profile_Scale

end module Dark_Matter_Profile_Scales
