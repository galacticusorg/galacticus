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
  private
  public :: Dark_Matter_Profile_Scale

  ! Initialization state.
  logical :: moduleInitialized=.false.

  ! Scale calculation option.
  logical :: darkMatterProfileScaleCorrectForConcentrationDefinition

contains

  double precision function Dark_Matter_Profile_Scale(node,meanConcentration,concentrationMethod)
    !% Compute the scale radius of the dark matter profile of {\normalfont \ttfamily node}.
    use Galacticus_Nodes
    use Dark_Matter_Profiles
    use Dark_Matter_Profiles_Concentration
    use Root_Finder
    use Virial_Density_Contrast
    use Dark_Matter_Halo_Scales
    use Cosmology_Parameters
    use Cosmology_Functions
    use Galacticus_Calculations_Resets
    use Numerical_Constants_Math
    use Input_Parameters
    implicit none    
    type            (treeNode                                          ), pointer, intent(inout)           :: node
    class           (darkMatterProfileConcentrationClass               ), target , intent(inout), optional :: concentrationMethod
    logical                                                                      , intent(in   ), optional :: meanConcentration
    class           (darkMatterProfileConcentrationClass               ), pointer                          :: darkMatterProfileConcentrationDefinition
    class           (virialDensityContrastClass                        ), pointer                          :: virialDensityContrastDefinition         , virialDensityContrast_
    class           (darkMatterProfileClass                            ), pointer                          :: darkMatterProfileDefinition             , darkMatterProfile_
    class           (cosmologyParametersClass                          ), pointer                          :: cosmologyParameters_
    class           (cosmologyFunctionsClass                           ), pointer                          :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                          ), pointer                          :: darkMatterHaloScale_
    type            (treeNode                                          ), pointer                          :: workNode
    class           (nodeComponentBasic                                ), pointer                          :: workBasic                               , basic
    class           (nodeComponentDarkMatterProfile                    ), pointer                          :: workDarkMatterProfile
    type            (rootFinder                                        ), save                             :: finder
    !$omp threadprivate(finder)
    type            (darkMatterHaloScaleVirialDensityContrastDefinition)                                   :: darkMatterHaloScaleDefinition
    double precision                                                                                       :: mass                                    , massDefinition        , &
         &                                                                                                    concentration                           , concentrationOriginal

    ! Initialize as necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Dark_Matter_Profile_Scale_Initialize)
       if (.not.moduleInitialized) then
          ! Get parameter controlling scale calculation method.
          !@ <inputParameter>
          !@   <name>darkMatterProfileScaleCorrectForConcentrationDefinition</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     If true, then when computing dark matter profile scale radii using concentrations, any difference between the current definition of halo scales
          !@    (i.e. typically virial density contrast definitions) and density profiles and those assumed in measuring the concentrations will be taken into account.
          !@    If false, the concentration is applied blindly.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('darkMatterProfileScaleCorrectForConcentrationDefinition',darkMatterProfileScaleCorrectForConcentrationDefinition,defaultValue=.false.)
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
       cosmologyParameters_            => cosmologyParameters                                                           (                               )
       cosmologyFunctions_             => cosmologyFunctions                                                            (                               )
       virialDensityContrast_          => virialDensityContrast                                                         (                               )
       darkMatterProfile_              => darkMatterProfile                                                             (                               )
       virialDensityContrastDefinition => darkMatterProfileConcentrationDefinition          %densityContrastDefinition  (                               )
       darkMatterHaloScaleDefinition   =  darkMatterHaloScaleVirialDensityContrastDefinition                            (virialDensityContrastDefinition)
       darkMatterProfileDefinition     => darkMatterProfileConcentrationDefinition          %darkMatterProfileDefinition(                               )
       ! Get the basic component of the supplied node and extract its mass.
       basic => node %basic()
       mass  =  basic%mass ()
       ! Create a node and set the mass and time.
       workNode              => treeNode                  (                 )
       workBasic             => workNode%basic            (autoCreate=.true.)
       workDarkMatterProfile => workNode%darkMatterProfile(autoCreate=.true.)
       call workBasic%timeSet            (basic%time())
       call workBasic%timeLastIsolatedSet(basic%time())
       call finder   %tolerance          (                                               &
            &                             toleranceRelative  =1.0d-3                     &
            &                            )
       call finder   %rangeExpand        (                                               &
            &                             rangeExpandUpward  =2.0d0                    , &
            &                             rangeExpandDownward=0.5d0                    , &
            &                             rangeExpandType    =rangeExpandMultiplicative  &
            &                            )
       call finder   %rootFunction       (                                               &
            &                                                 massRootFunction           &
            &                            )
       massDefinition=finder%find(rootGuess=mass)
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
       deallocate(workNode                   )
       ! Destroy objects as necessary.
       if (darkMatterProfileDefinition%isFinalizable()) deallocate(darkMatterProfileDefinition)
    else
        Dark_Matter_Profile_Scale=+darkMatterHaloScale_%virialRadius(node) &
            &                    /concentrationOriginal
    end if
    ! Nullify the concentration definition so that it isn't automatically finalized.
    darkMatterProfileConcentrationDefinition => null()
    
  contains
    
    double precision function massRootFunction(massDefinitionTrial)
      !% Root function used to find the mass of a halo corresponding to the definition used for a particular concentration class.
      implicit none
      double precision            , intent(in   ) :: massDefinitionTrial
      type            (rootFinder), save          :: radiusFinder
      !$omp threadprivate(radiusFinder)
      double precision                            :: radiusOuterDefinition, concentrationDefinition, &
           &                                         radiusCore           , massOuter              , &
           &                                         radiusOuter
      
      ! Set the mass of the worker node.
      call workBasic%massSet(massDefinitionTrial)
      call Galacticus_Calculations_Reset                 (workNode)
      call darkMatterProfileDefinition  %calculationReset(workNode)
      ! Get outer radius for this trial definition mass.
      radiusOuterDefinition  =darkMatterHaloScaleDefinition           %virialRadius (workNode)
      ! Get concentration for this a trial definition mass.
      if (present(meanConcentration).and.meanConcentration) then
         ! We are simply using the mean concentration-mass relation here.
         concentrationDefinition=darkMatterProfileConcentrationDefinition%concentrationMean(workNode)
      else
         ! In this case we need to allow for possible scatter in the concentration mass relation. Therefore, we take the original
         ! concentration (which may include some scatter away from the mean relation) and scale it by the ratio of the mean
         ! concentrations for the corrected and original nodes.
         concentrationDefinition=+concentrationOriginal                                                &
               &                 *darkMatterProfileConcentrationDefinition%concentrationMean(workNode) &
               &                 /darkMatterProfileConcentrationDefinition%concentrationMean(    node)
      end if
      ! Get core radius.      
      radiusCore             =radiusOuterDefinition/concentrationDefinition
      call workDarkMatterProfile%scaleSet(radiusCore)
      call Galacticus_Calculations_Reset                 (workNode)
      call darkMatterProfileDefinition  %calculationReset(workNode)
      ! Solve for radius which encloses required non-alt density contrast.
      call radiusFinder%tolerance   (                                                 &
           &                         toleranceRelative  =1.0d-3                       &
           &                        )
      call radiusFinder%rangeExpand (                                                 &
           &                         rangeExpandUpward  =2.0d0                      , &
           &                         rangeExpandDownward=0.5d0                      , &
           &                         rangeExpandType    =rangeExpandMultiplicative    &
           &                        )
      call radiusFinder%rootFunction(                                                 &
           &                                             densityContrastRootFunction  &
           &                        )
      radiusOuter=radiusFinder%find(rootGuess=radiusOuterDefinition)
      ! Get the mass within this radius.
      massOuter  =darkMatterProfileDefinition%enclosedMass(workNode,radiusOuter)
      ! Return root function.
      massRootFunction=massOuter-mass
      return
    end function massRootFunction
    
    double precision function densityContrastRootFunction(radiusTrial)
      !% Root function used to find the radius in a dark matter profile which encloses a density contrast equal to the currently
      !% specified density contrast.
      implicit none
      double precision, intent(in   ) :: radiusTrial
      double precision                :: massTrial  , densityContrastTrial
      
      massTrial                  = darkMatterProfileDefinition%enclosedMass(workNode,radiusTrial)
      densityContrastTrial       =+3.0d0                                                      &
           &                      *massTrial                                                  &
           &                      /4.0d0                                                      &
           &                      /Pi                                                         &
           &                      /cosmologyFunctions_%matterDensityEpochal(workBasic%time()) &
           &                      /radiusTrial**3
      densityContrastRootFunction=+densityContrastTrial-virialDensityContrast_%densityContrast(workBasic%mass(),workBasic%time())
      return
    end function densityContrastRootFunction

  end function Dark_Matter_Profile_Scale

end module Dark_Matter_Profile_Scales
