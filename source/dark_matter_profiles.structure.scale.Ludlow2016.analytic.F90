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

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Linear_Growth             , only : linearGrowthClass

  !!{
  An implementation of dark matter halo profile concentrations using the \cite{ludlow_mass-concentration-redshift_2016}
  algorithm with an analytic model for formation time.
  !!}

  !![
  <darkMatterProfileScaleRadius name="darkMatterProfileScaleRadiusLudlow2016Analytic">
    <description>
      Dark matter halo scale radii are computed using the algorithm of \cite{ludlow_mass-concentration-redshift_2016} with an
      analytic model for formation time.
     </description>
  </darkMatterProfileScaleRadius>
  !!]
  type, extends(darkMatterProfileScaleRadiusLudlow2016) :: darkMatterProfileScaleRadiusLudlow2016Analytic
     !!{     
     A dark matter halo profile scale radii class implementing the algorithm of \cite{ludlow_mass-concentration-redshift_2016}
     with an analytic model for formation time.     
     !!}
     private
     class(criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class(linearGrowthClass            ), pointer :: linearGrowth_             => null()
   contains
     final             ::                                 ludlow2016AnalyticDestructor
     procedure, nopass :: formationTimeRoot            => ludlow2016AnalyticFormationTimeRoot
     procedure         :: formationTimeRootFunctionSet => ludlow2016AnalyticFormationTimeRootFunctionSet
     procedure         :: requireBranchHistory         => ludlow2016AnalyticRequireBranchHistory
  end type darkMatterProfileScaleRadiusLudlow2016Analytic

  interface darkMatterProfileScaleRadiusLudlow2016Analytic
     !!{
     Constructors for the \refClass{darkMatterProfileScaleRadiusLudlow2016Analytic} dark matter halo profile scale radius class.
     !!}
     module procedure ludlow2016AnalyticConstructorParameters
     module procedure ludlow2016AnalyticConstructorInternal
  end interface darkMatterProfileScaleRadiusLudlow2016Analytic

contains

  function ludlow2016AnalyticConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily ludlow2016Analytic} dark matter halo profile scale radius class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(darkMatterProfileScaleRadiusLudlow2016Analytic)                :: self
    type(inputParameters                               ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass           ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass          ), pointer       :: cosmologyParameters_
    class           (darkMatterProfileScaleRadiusClass ), pointer       :: darkMatterProfileScaleRadius_
    class           (virialDensityContrastClass        ), pointer       :: virialDensityContrast_
    class           (darkMatterProfileDMOClass         ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass          ), pointer       :: darkMatterHaloScale_
    class           (criticalOverdensityClass          ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass     ), pointer       :: cosmologicalMassVariance_
    class           (linearGrowthClass                 ), pointer       :: linearGrowth_
    double precision                                                    :: C                            , f

    !![
    <inputParameter>
      <name>C</name>
      <source>parameters</source>
      <defaultValue>650.0d0</defaultValue>
      <description>The parameter $C$ appearing in the halo concentration algorithm of \cite[][see their footnote~7]{ludlow_mass-concentration-redshift_2016}.</description>
    </inputParameter>
    <inputParameter>
      <name>f</name>
      <source>parameters</source>
      <defaultValue>0.02d0</defaultValue>
      <description>The parameter $f$ appearing in the halo concentration algorithm of \cite{ludlow_mass-concentration-redshift_2016}.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    <objectBuilder class="cosmologyParameters"          name="cosmologyParameters_"          source="parameters"/>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    <objectBuilder class="virialDensityContrast"        name="virialDensityContrast_"        source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    <objectBuilder class="criticalOverdensity"          name="criticalOverdensity_"          source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"     name="cosmologicalMassVariance_"     source="parameters"/>
    <objectBuilder class="linearGrowth"                 name="linearGrowth_"                 source="parameters"/>
    !!]
    self=darkMatterProfileScaleRadiusLudlow2016Analytic(C,f,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileScaleRadius_,virialDensityContrast_,darkMatterProfileDMO_,darkMatterHaloScale_,criticalOverdensity_,cosmologicalMassVariance_,linearGrowth_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"          />
    <objectDestructor name="cosmologyParameters_"         />
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    <objectDestructor name="virialDensityContrast_"       />
    <objectDestructor name="darkMatterProfileDMO_"        />
    <objectDestructor name="darkMatterHaloScale_"         />
    <objectDestructor name="criticalOverdensity_"         />
    <objectDestructor name="cosmologicalMassVariance_"    />
    <objectDestructor name="linearGrowth_"                />
    !!]
    return
  end function ludlow2016AnalyticConstructorParameters

  function ludlow2016AnalyticConstructorInternal(C,f,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileScaleRadius_,virialDensityContrast_,darkMatterProfileDMO_,darkMatterHaloScale_,criticalOverdensity_,cosmologicalMassVariance_,linearGrowth_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileScaleRadiusLudlow2016Analytic} dark matter halo profile scale radius class.
    !!}
    implicit none
    type            (darkMatterProfileScaleRadiusLudlow2016Analytic)                        :: self
    double precision                                                , intent(in   )         :: C                            , f
    class           (cosmologyFunctionsClass                       ), intent(in   ), target :: cosmologyFunctions_
    class           (cosmologyParametersClass                      ), intent(in   ), target :: cosmologyParameters_
    class           (darkMatterProfileScaleRadiusClass             ), intent(in   ), target :: darkMatterProfileScaleRadius_
    class           (virialDensityContrastClass                    ), intent(in   ), target :: virialDensityContrast_
    class           (darkMatterProfileDMOClass                     ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass                      ), intent(in   ), target :: darkMatterHaloScale_
    class           (criticalOverdensityClass                      ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                 ), intent(in   ), target :: cosmologicalMassVariance_
    class           (linearGrowthClass                             ), intent(in   ), target :: linearGrowth_
    !![
    <constructorAssign variables="C, f, *cosmologyFunctions_, *cosmologyParameters_, *darkMatterProfileScaleRadius_, *virialDensityContrast_, *darkMatterProfileDMO_, *darkMatterHaloScale_, *criticalOverdensity_, *cosmologicalMassVariance_, *linearGrowth_"/>
    !!]

    ! No need to seek formation times as they are analytically computable.
    self%timeFormationSeekDelta        =  -huge(0.0d0)
    ! Find the density contrast as used to define masses by Ludlow et al. (2016).
    self%densityContrast               =  +200.0d0                                 &
         &                                /self%cosmologyParameters_%omegaMatter()
    return
  end function ludlow2016AnalyticConstructorInternal

  subroutine ludlow2016AnalyticDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileScaleRadiusLudlow2016Analytic} dark matter halo profile scale radius class.
    !!}
    implicit none
    type(darkMatterProfileScaleRadiusLudlow2016Analytic), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%linearGrowth_"            />
    !!]
    return
  end subroutine ludlow2016AnalyticDestructor

  subroutine ludlow2016AnalyticFormationTimeRootFunctionSet(self,finder)
    !!{
    Initialize the finder object to compute the relevant formation history.
    !!}
    use :: Root_Finder, only : rootFinder
    implicit none
    class           (darkMatterProfileScaleRadiusLudlow2016Analytic), intent(inout) :: self
    type            (rootFinder                                    ), intent(inout) :: finder
    double precision                                                , parameter     :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-4
    !$GLC attributes unused :: self

    finder=rootFinder(                                                       &
         &            rootFunction     =ludlow2016AnalyticFormationTimeRoot, &
         &            toleranceAbsolute=toleranceAbsolute                  , &
         &            toleranceRelative=toleranceRelative                    &
         &           )
    return
  end subroutine ludlow2016AnalyticFormationTimeRootFunctionSet

  double precision function ludlow2016AnalyticFormationTimeRoot(timeFormation) result(formationTimeRoot)
    !!{
    Function used to find the formation time of a halo in the {\normalfont \ttfamily ludlow2016Analytic} scale radius algorithm.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic
    use :: Error                               , only : Error_Report
    implicit none
    double precision                    , intent(in   ) :: timeFormation
    class           (nodeComponentBasic), pointer       :: basic
    double precision                                    :: massCurrent    , massFormation    , &
         &                                                 timeCurrent    , timePresent      , &
         &                                                 varianceCurrent, varianceFormation

    basic             => states(stateCount)%node%basic()
    massCurrent       =  Dark_Matter_Profile_Mass_Definition(                                                                                   &
         &                                                                          states(stateCount)%node                                   , &
         &                                                                          ludlow2016DensityContrast(states(stateCount),basiC%time()), &
         &                                                   cosmologyParameters_  =states(stateCount)%self%cosmologyParameters_              , &
         &                                                   cosmologyFunctions_   =states(stateCount)%self%cosmologyFunctions_               , &
         &                                                   virialDensityContrast_=states(stateCount)%self%virialDensityContrast_            , &
         &                                                   darkMatterProfileDMO_ =states(stateCount)%self%darkMatterProfileDMO_               &
         &                                                  )
    timeCurrent=basic%time()
    timePresent=states(stateCount)%self%cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)
    select type (self => states(stateCount)%self)
    class is (darkMatterProfileScaleRadiusLudlow2016Analytic)
       varianceCurrent  =+              self%cosmologicalMassVariance_%rootVariance(time=timePresent  ,mass=       massCurrent                             )**2
       varianceFormation=+              self%cosmologicalMassVariance_%rootVariance(time=timePresent  ,mass=self%f*massCurrent                             )**2
       if (varianceFormation <= varianceCurrent) then
          ! No valid solution can be found at this mass. Return a negative root - this will cause the parent class to fall through to the alternative method.
          massFormation =+0.0d0
       else
          massFormation =+massCurrent                                                                                                                           &
               &         *erfc(                                                                                                                                 &
               &               +    (                                                                                                                           &
               &                       +self%criticalOverdensity_     %value       (time=timeFormation,mass=self%f*massCurrent,node=states(stateCount)%node)    &
               &                       /self%linearGrowth_            %value       (time=timeFormation                                                     )    &
               &                       -self%criticalOverdensity_     %value       (time=timeCurrent  ,mass=       massCurrent,node=states(stateCount)%node)    &
               &                       /self%linearGrowth_            %value       (time=timeCurrent                                                       )    &
               &                    )                                                                                                                           &
               &               /sqrt(                                                                                                                           &
               &                     +2.0d0                                                                                                                     &
               &                     *(                                                                                                                         &
               &                       +varianceFormation                                                                                                       &
               &                       -varianceCurrent                                                                                                         &
               &                      )                                                                                                                         &
               &                    )                                                                                                                           &
               &              )
       end if
    class default
       massFormation=-huge(0.0d0)
       call Error_Report('incorrect class'//{introspection:location})
    end select
    formationTimeRoot =  +                   massFormation          &
         &               -states(stateCount)%massHaloCharacteristic
    return
  end function ludlow2016AnalyticFormationTimeRoot

  logical function ludlow2016AnalyticRequireBranchHistory(self) result(requireBranchHistory)
    !!{
    Specify if the branch history is required for the scale radius calculation.
    !!}
    implicit none
    class(darkMatterProfileScaleRadiusLudlow2016Analytic), intent(inout) :: self
    !$GLC attributes unused :: self

    requireBranchHistory=.false.
    return
  end function ludlow2016AnalyticRequireBranchHistory
