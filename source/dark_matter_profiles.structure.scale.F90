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

  double precision function Dark_Matter_Profile_Scale(node,concentrationMethod)
    !% Compute the scale radius of the dark matter profile of {\tt node}.
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
    class           (darkMatterProfileConcentrationClass               ), pointer                          :: darkMatterProfileConcentrationDefinition
    class           (virialDensityContrastClass                        ), pointer                          :: virialDensityContrastDefinition         , virialDensityContrast_
    class           (darkMatterProfileClass                            ), pointer                          :: darkMatterProfileDefinition             , darkMatterProfile_
    class           (cosmologyParametersClass                          ), pointer                          :: cosmologyParameters_
    class           (cosmologyFunctionsClass                           ), pointer                          :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                          ), pointer                          :: darkMatterHaloScale_
    type            (treeNode                                          ), pointer                          :: workNode
    class           (nodeComponentBasic                                ), pointer                          :: workBasic                               , basic
    class           (nodeComponentDarkMatterProfile                    ), pointer                          :: workDarkMatterProfile
    type            (darkMatterHaloScaleVirialDensityContrastDefinition)                                   :: darkMatterHaloScaleDefinition
    double precision                                                                                       :: mass                                    , massDefinition
    type            (rootFinder                                        )                                   :: finder

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
            &                             toleranceRelative  =1.0d-6                     &
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
       Dark_Matter_Profile_Scale=+darkMatterHaloScaleDefinition           %virialRadius (workNode) &
            &                    /darkMatterProfileConcentrationDefinition%concentration(workNode)
       call workNode%destroy()
       deallocate(workNode)
    else
       Dark_Matter_Profile_Scale= darkMatterHaloScale_                    %virialRadius (node) &
            &                    /darkMatterProfileConcentrationDefinition%concentration(node)
    end if
    ! Nullify the concentration definition so that it isn't automatically finalized.
    darkMatterProfileConcentrationDefinition => null()
    
  contains
    
    double precision function massRootFunction(massDefinitionTrial)
      !% Root function used to find the mass of a halo corresponding to the definition used for a particular concentration class.
      implicit none
      double precision            , intent(in   ) :: massDefinitionTrial
      double precision                            :: radiusOuterDefinition, concentrationDefinition  , &
           &                                         radiusCore           , massOuterDefinition      , &
           &                                         radiusOuter          , densityContrastDefinition, &
           &                                         densityContrast      , massOuter
      type            (rootFinder)                :: radiusFinder
      
      ! Set the mass of the worker node.
      call workBasic%massSet(massDefinitionTrial)
      call Galacticus_Calculations_Reset                 (workNode)
      call darkMatterProfileDefinition  %calculationReset(workNode)
      ! Get outer radius for this trial definition mass.
      radiusOuterDefinition  =darkMatterHaloScaleDefinition           %virialRadius (workNode                      )
      ! Get mass normalization.
      massOuterDefinition    =darkMatterProfileDefinition             %enclosedMass (workNode,radiusOuterDefinition)
      ! Get concentration for this a trial definition mass.
      concentrationDefinition=darkMatterProfileConcentrationDefinition%concentration(workNode                      )
      ! Get core radius.
      radiusCore             =radiusOuterDefinition/concentrationDefinition
      call workDarkMatterProfile%scaleSet(radiusCore)
      call Galacticus_Calculations_Reset                 (workNode)
      call darkMatterProfileDefinition  %calculationReset(workNode)
      ! Solve for radius which encloses required non-alt density contrast.
      call radiusFinder%tolerance   (                                                 &
           &                         toleranceRelative  =1.0d-6                       &
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
      densityContrastTrial       =+3.0d0                                  &
           &                      *massTrial                              &
           &                      /4.0d0                                  &
           &                      /Pi                                     &
           &                      /cosmologyParameters_%OmegaMatter    () &
           &                      /cosmologyParameters_%densityCritical() &
           &                      /radiusTrial**3
      densityContrastRootFunction=+densityContrastTrial-virialDensityContrast_%densityContrast(workBasic%time())
      return
    end function densityContrastRootFunction

  end function Dark_Matter_Profile_Scale

end module Dark_Matter_Profile_Scales
