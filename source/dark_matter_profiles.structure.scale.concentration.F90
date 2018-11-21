!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% An implementation of dark matter halo profile scale radii in which radii are computed from the concentration.
  
  use Dark_Matter_Profiles              , only : darkMatterProfile                                 , darkMatterProfileClass
  use Dark_Matter_Profiles_Concentration, only : darkMatterProfileConcentration                    , darkMatterProfileConcentrationClass
  use Dark_Matter_Halo_Scales           , only : darkMatterHaloScale                               , darkMatterHaloScaleClass           , &
       &                                         darkMatterHaloScaleVirialDensityContrastDefinition
  use Cosmology_Functions               , only : cosmologyFunctions                                , cosmologyFunctionsClass
  use Cosmology_Parameters              , only : cosmologyParameters                               , cosmologyParametersClass
  use Virial_Density_Contrast           , only : virialDensityContrast                             , virialDensityContrastClass
  
  !# <darkMatterProfileScaleRadius name="darkMatterProfileScaleRadiusConcentration">
  !#  <description>Dark matter halo scale radii are computed from the concentration.</description>
  !# </darkMatterProfileScaleRadius>
  type, extends(darkMatterProfileScaleRadiusClass) :: darkMatterProfileScaleRadiusConcentration
     !% A dark matter halo profile scale radius class which computes radii from the concentration.
     private
     class           (cosmologyParametersClass                          ), pointer :: cosmologyParameters_              => null()
     class           (cosmologyFunctionsClass                           ), pointer :: cosmologyFunctions_               => null()
     class           (darkMatterHaloScaleClass                          ), pointer :: darkMatterHaloScale_              => null()
     class           (darkMatterProfileClass                            ), pointer :: darkMatterProfile_                => null(), darkMatterProfileDefinition     => null()
     class           (virialDensityContrastClass                        ), pointer :: virialDensityContrast_            => null(), virialDensityContrastDefinition => null()
     class           (darkMatterProfileConcentrationClass               ), pointer :: darkMatterProfileConcentration_   => null()
     type            (darkMatterHaloScaleVirialDensityContrastDefinition)          :: darkMatterHaloScaleDefinition
     logical                                                                       :: correctForConcentrationDefinition          , useMeanConcentration
     double precision                                                              :: massRatioPrevious
   contains
     final     ::           concentrationDestructor
     procedure :: radius => concentrationRadius
  end type darkMatterProfileScaleRadiusConcentration
  
  interface darkMatterProfileScaleRadiusConcentration
     !% Constructors for the {\normalfont \ttfamily concentration} dark matter halo profile scale radius class.
     module procedure concentrationConstructorParameters
     module procedure concentrationConstructorInternal
  end interface darkMatterProfileScaleRadiusConcentration

contains

  function concentrationConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily concentration} dark matter halo profile scale radius class which takes a
    !% parameter list as input.
    use Input_Parameters
    implicit none
    type   (darkMatterProfileScaleRadiusConcentration)                :: self
    type   (inputParameters                          ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_
    class  (cosmologyParametersClass                 ), pointer       :: cosmologyParameters_
    class  (darkMatterHaloScaleClass                 ), pointer       :: darkMatterHaloScale_
    class  (darkMatterProfileClass                   ), pointer       :: darkMatterProfile_ 
    class  (virialDensityContrastClass               ), pointer       :: virialDensityContrast_
    class  (darkMatterProfileConcentrationClass      ), pointer       :: darkMatterProfileConcentration_
    logical                                                           :: correctForConcentrationDefinition, useMeanConcentration

    !# <inputParameter>
    !#   <name>correctForConcentrationDefinition</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true, then when computing dark matter profile scale radii using concentrations, any difference between the current definition of halo scales
    !#     (i.e. typically virial density contrast definitions) and density profiles and those assumed in measuring the concentrations will be taken into account.
    !#     If false, the concentration is applied blindly.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>useMeanConcentration</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true, then when computing dark matter profile scale radii using concentrations do not account for any possible scatter in the concentration-mass relation.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters"            name="cosmologyParameters_"            source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"             name="cosmologyFunctions_"             source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"            name="darkMatterHaloScale_"            source="parameters"/>
    !# <objectBuilder class="darkMatterProfile"              name="darkMatterProfile_"              source="parameters"/>
    !# <objectBuilder class="virialDensityContrast"          name="virialDensityContrast_"          source="parameters"/>
    !# <objectBuilder class="darkMatterProfileConcentration" name="darkMatterProfileConcentration_" source="parameters"/>
    self=darkMatterProfileScaleRadiusConcentration(correctForConcentrationDefinition,useMeanConcentration,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfile_,virialDensityContrast_,darkMatterProfileConcentration_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function concentrationConstructorParameters

  function concentrationConstructorInternal(correctForConcentrationDefinition,useMeanConcentration,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfile_,virialDensityContrast_,darkMatterProfileConcentration_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily concentration} dark matter halo profile scale radius class.
    implicit none
    type   (darkMatterProfileScaleRadiusConcentration)                        :: self
    class  (cosmologyParametersClass                 ), intent(in   ), target :: cosmologyParameters_
    class  (cosmologyFunctionsClass                  ), intent(in   ), target :: cosmologyFunctions_
    class  (darkMatterHaloScaleClass                 ), intent(in   ), target :: darkMatterHaloScale_
    class  (darkMatterProfileClass                   ), intent(in   ), target :: darkMatterProfile_ 
    class  (virialDensityContrastClass               ), intent(in   ), target :: virialDensityContrast_
    class  (darkMatterProfileConcentrationClass      ), intent(in   ), target :: darkMatterProfileConcentration_
    logical                                           , intent(in   )         :: correctForConcentrationDefinition, useMeanConcentration
    !# <constructorAssign variables="correctForConcentrationDefinition, useMeanConcentration, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterHaloScale_, *darkMatterProfile_, *virialDensityContrast_, *darkMatterProfileConcentration_"/>

    ! Get definitions of virial density contrast and dark matter profile as used by the concentration definition,
    self%virialDensityContrastDefinition => self%darkMatterProfileConcentration_%  densityContrastDefinition()
    self%darkMatterProfileDefinition     => self%darkMatterProfileConcentration_%darkMatterProfileDefinition()
    self%darkMatterHaloScaleDefinition   =  darkMatterHaloScaleVirialDensityContrastDefinition(                                      &
         &                                                                                     self%cosmologyParameters_           , &
         &                                                                                     self%cosmologyFunctions_            , &
         &                                                                                     self%virialDensityContrastDefinition  &
         &                                                                                    )    
    self%massRatioPrevious               =  2.0d0
  return
  end function concentrationConstructorInternal

  subroutine concentrationDestructor(self)
    !% Destructor for the {\normalfont \ttfamily concentration} dark matter halo profile scale radius class.
    implicit none
    type(darkMatterProfileScaleRadiusConcentration), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"           />
    !# <objectDestructor name="self%cosmologyFunctions_"            />
    !# <objectDestructor name="self%darkMatterHaloScale_"           />
    !# <objectDestructor name="self%darkMatterProfile_"             />
    !# <objectDestructor name="self%virialDensityContrast_"         />
    !# <objectDestructor name="self%darkMatterProfileConcentration_"/>
    return
  end subroutine concentrationDestructor
  
  double precision function concentrationRadius(self,node)
    !% Compute the scale radius of the dark matter profile of {\normalfont \ttfamily node}.
    use Galacticus_Nodes              , only : nodeComponentBasic, nodeComponentDarkMatterProfile
    use Root_Finder
    use Galacticus_Calculations_Resets
    use Numerical_Constants_Math
    use Numerical_Comparison
    implicit none
    class           (darkMatterProfileScaleRadiusConcentration), intent(inout) :: self
    type            (treeNode                                 ), intent(inout) :: node
    class           (nodeComponentBasic                       ), pointer       :: basic
    type            (treeNode                                 ), pointer       :: workNode
    class           (nodeComponentBasic                       ), pointer       :: workBasic
    class           (nodeComponentDarkMatterProfile           ), pointer       :: workDarkMatterProfile
    double precision                                           , parameter     :: massRatioBuffer      =1.1d0, massRatioShrink=0.99d0
    type            (rootFinder                               )                :: finder
    double precision                                                           :: mass                       , massDefinition        , &
         &                                                                        concentration              , massRatio
    double precision                                                           :: concentrationOriginal

    ! Find the original concentration.
    if (self%useMeanConcentration) then
       concentrationOriginal=self%darkMatterProfileConcentration_%concentrationMean(node)
    else
       concentrationOriginal=self%darkMatterProfileConcentration_%concentration    (node)
    end if
    ! Determine if concentration must be corrected.
    if (self%correctForConcentrationDefinition) then
       ! Get the basic component of the supplied node and extract its mass.
       basic => node %basic()
       mass  =  basic%mass ()
       ! If there is no difference between the alt and non-alt virial density contrasts, then no correction need be made.
       if     (                                                                                                      &
            &  Values_Differ(                                                                                        &
            &                       self%virialDensityContrast_         %densityContrast(basic%mass(),basic%time()), &
            &                       self%virialDensityContrastDefinition%densityContrast(basic%mass(),basic%time()), &
            &                relTol=1.0d-6                                                                           &
            &               )                                                                                        &
            & ) then
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
               &                                         rangeExpandUpward            =1.0d0*self%massRatioPrevious , &
               &                                         rangeExpandDownward          =1.0d0/self%massRatioPrevious , &
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
          if (massRatio > self%massRatioPrevious) then
             self%massRatioPrevious=massRatio
          else
             self%massRatioPrevious=massRatioShrink*self%massRatioPrevious
          end if
          ! Update the work node properties and computed concentration.
          call workBasic%massSet(massDefinition)
          call Galacticus_Calculations_Reset                 (workNode)
          call self%darkMatterProfileDefinition  %calculationReset(workNode)
          ! Find the concentration.
          if (self%useMeanConcentration) then
             ! We are simply using the mean concentration-mass relation here.
             concentration=+self%darkMatterProfileConcentration_%concentrationMean(workNode)
          else
             ! In this case we need to allow for possible scatter in the concentration mass relation. Therefore, we take the original
             ! concentration (which may include some scatter away from the mean relation) and scale it by the ratio of the mean
             ! concentrations for the corrected and original nodes.
             concentration=+                                     concentrationOriginal           &
                  &        *self%darkMatterProfileConcentration_%concentrationMean    (workNode) &
                  &        /self%darkMatterProfileConcentration_%concentrationMean    (    node)
          end if
          concentrationRadius=+self%darkMatterHaloScaleDefinition%virialRadius (workNode) &
               &              /                                   concentration
          call workNode%destroy()
          deallocate(workNode)
       else
          concentrationRadius=+self%darkMatterHaloScale_%virialRadius         (node) &
               &              /                          concentrationOriginal
       end if
    else
       concentrationRadius=+self%darkMatterHaloScale_%virialRadius         (node) &
            &              /                          concentrationOriginal
    end if
    return

  contains
    
    double precision function massRootFunction(massDefinitionTrial)
      !% Root function used to find the mass of a halo corresponding to the definition used for a particular concentration class.
      implicit none
      double precision, intent(in   ) :: massDefinitionTrial
      double precision                :: radiusOuterDefinition, concentrationDefinition, &
           &                             radiusCore           , massOuter              , &
           &                             radiusOuter          , densityOuter

      ! Set the mass of the worker node.
      call workBasic%massSet(massDefinitionTrial)
      call Galacticus_Calculations_Reset                      (workNode)
      call self%darkMatterHaloScaleDefinition%calculationReset(workNode)
      ! Get outer radius for this trial definition mass.
      radiusOuterDefinition=self%darkMatterHaloScaleDefinition%virialRadius(workNode)
      ! Get concentration for this a trial definition mass.
      if (self%useMeanConcentration) then
         ! We are simply using the mean concentration-mass relation here.
         concentrationDefinition=self%darkMatterProfileConcentration_%concentrationMean(workNode)
      else
         ! In this case we need to allow for possible scatter in the concentration mass relation. Therefore, we take the original
         ! concentration (which may include some scatter away from the mean relation) and scale it by the ratio of the mean
         ! concentrations for the corrected and original nodes.
         concentrationDefinition=+                                     concentrationOriginal           &
              &                  *self%darkMatterProfileConcentration_%concentrationMean    (workNode) &
              &                  /self%darkMatterProfileConcentration_%concentrationMean    (node    )
      end if
      ! Get core radius.      
      radiusCore=radiusOuterDefinition/concentrationDefinition
      call workDarkMatterProfile%scaleSet(radiusCore)
      call Galacticus_Calculations_Reset                 (workNode)
      call self%darkMatterProfileDefinition  %calculationReset(workNode)
      ! Find the non-alt density.
      densityOuter=+self%cosmologyFunctions_   %matterDensityEpochal(                 workBasic%time()) &
           &       *self%virialDensityContrast_%densityContrast     (workBasic%mass(),workBasic%time())      
      ! Solve for radius which encloses required non-alt density.
      radiusOuter=self%darkMatterProfileDefinition%radiusEnclosingDensity(workNode,densityOuter)
      ! Get the mass within this radius.
      massOuter  =self%darkMatterProfileDefinition%enclosedMass          (workNode, radiusOuter)
      ! Return root function.
      massRootFunction=massOuter-mass
      return
    end function massRootFunction

  end function concentrationRadius
