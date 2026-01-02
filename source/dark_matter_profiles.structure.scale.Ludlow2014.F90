!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  An implementation of dark matter halo profile concentrations using the \cite{ludlow_mass-concentration-redshift_2014}
  algorithm.
  !!}

  !![
  <darkMatterProfileScaleRadius name="darkMatterProfileScaleRadiusLudlow2014">
   <description>Dark matter halo scale radii are computed using the algorithm of \cite{ludlow_mass-concentration-redshift_2014}.</description>
  </darkMatterProfileScaleRadius>
  !!]
  type, extends(darkMatterProfileScaleRadiusLudlow2016) :: darkMatterProfileScaleRadiusLudlow2014
     !!{
     A dark matter halo profile scale radii class implementing the algorithm of \cite{ludlow_mass-concentration-redshift_2014}.
     !!}
     private
   contains
     procedure, nopass :: formationTimeRoot            => ludlow2014FormationTimeRoot
     procedure         :: formationTimeRootFunctionSet => ludlow2014FormationTimeRootFunctionSet
  end type darkMatterProfileScaleRadiusLudlow2014

  interface darkMatterProfileScaleRadiusLudlow2014
     !!{
     Constructors for the \refClass{darkMatterProfileScaleRadiusLudlow2014} dark matter halo profile concentration class.
     !!}
     module procedure ludlow2014ConstructorParameters
     module procedure ludlow2014ConstructorInternal
  end interface darkMatterProfileScaleRadiusLudlow2014

contains

  function ludlow2014ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily ludlow2014} dark matter halo profile concentration class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(darkMatterProfileScaleRadiusLudlow2014)                :: self
    type(inputParameters                       ), intent(inout) :: parameters

    self%darkMatterProfileScaleRadiusLudlow2016=darkMatterProfileScaleRadiusLudlow2016(parameters)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function ludlow2014ConstructorParameters

  function ludlow2014ConstructorInternal(C,f,timeFormationSeekDelta,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileScaleRadius_,virialDensityContrast_,darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileScaleRadiusLudlow2014} dark matter halo profile concentration class.
    !!}
    implicit none
    type            (darkMatterProfileScaleRadiusLudlow2014)                        :: self
    double precision                                        , intent(in   )         :: C                            , f, &
         &                                                                             timeFormationSeekDelta
    class           (cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
    class           (cosmologyParametersClass              ), intent(in   ), target :: cosmologyParameters_
    class           (darkMatterProfileScaleRadiusClass     ), intent(in   ), target :: darkMatterProfileScaleRadius_
    class           (virialDensityContrastClass            ), intent(in   ), target :: virialDensityContrast_
    class           (darkMatterProfileDMOClass             ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass              ), intent(in   ), target :: darkMatterHaloScale_

    self%darkMatterProfileScaleRadiusLudlow2016=darkMatterProfileScaleRadiusLudlow2016(C,f,timeFormationSeekDelta,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileScaleRadius_,virialDensityContrast_,darkMatterProfileDMO_,darkMatterHaloScale_)
    return
  end function ludlow2014ConstructorInternal

  subroutine ludlow2014FormationTimeRootFunctionSet(self,finder)
    !!{
    Initialize the finder object to compute the relevant formation history.
    !!}
    use :: Root_Finder, only : rootFinder
    implicit none
    class           (darkMatterProfileScaleRadiusLudlow2014), intent(inout) :: self
    type            (rootFinder                            ), intent(inout) :: finder
    double precision                                        , parameter     :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-4
    !$GLC attributes unused :: self

    finder=rootFinder(                                               &
         &            rootFunction     =ludlow2014FormationTimeRoot, &
         &            toleranceAbsolute=toleranceAbsolute          , &
         &            toleranceRelative=toleranceRelative            &
         &           )
    return
  end subroutine ludlow2014FormationTimeRootFunctionSet

  double precision function ludlow2014FormationTimeRoot(timeFormation)
    !!{
    Function used to find the formation time of a halo in the {\normalfont \ttfamily ludlow2014} concentration algorithm.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    double precision                    , intent(in   ) :: timeFormation
    type            (treeNode          ), pointer       :: nodeBranch   , nodeChild        , &
         &                                                 nodeSibling
    class           (nodeComponentBasic), pointer       :: basicBranch  , basicChild       , &
         &                                                 basicSibling
    double precision                                    :: massBranch   , massAccretionRate, &
         &                                                 massSiblings , massProgenitor

    nodeBranch => states(stateCount)%node
    massBranch =  0.0d0
    do while (associated(nodeBranch))
       basicBranch => nodeBranch%basic()
       if (associated(nodeBranch%firstChild).and.basicBranch%time() >= timeFormation) then
          nodeChild => nodeBranch%firstChild
          do while (associated(nodeChild))
             basicChild => nodeChild%basic()
             if (basicChild%time() < timeFormation) then
                ! Interpolate in mass for primary progenitors.
                massSiblings =  0.0d0
                nodeSibling  => nodeChild
                do while (associated(nodeSibling))
                   basicSibling =>  nodeSibling %basic  ()
                   massSiblings =  +massSiblings                                                                                                                                                         &
                        &          +Dark_Matter_Profile_Mass_Definition(                                                                                                                                 &
                        &                                                                         nodeSibling                                                                                          , &
                        &                                                                      +  states(stateCount)%self                    %densityContrast                                            &
                        &                                                                      *(                                                                                                        &
                        &                                                                        +states(stateCount)%cosmologyFunctions_%hubbleParameterEpochal(time           =basicSibling%time())     &
                        &                                                                        /states(stateCount)%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=             1.0d0 )     &
                        &                                                                       )                                                                                                   **2  &
                        &                                                                      *  states(stateCount)%cosmologyFunctions_%expansionFactor       (                basicSibling%time())**3, &
                        &                                               cosmologyParameters_  =states(stateCount)%self%cosmologyParameters_                                                            , &
                        &                                               cosmologyFunctions_   =states(stateCount)%self%cosmologyFunctions_                                                             , &
                        &                                               virialDensityContrast_=states(stateCount)%self%virialDensityContrast_                                                            &
                        &                                              )
                   nodeSibling  =>  nodeSibling %sibling
                end do
                massAccretionRate=+(                                                                                                                                                                     &
                     &              +Dark_Matter_Profile_Mass_Definition(                                                                                                                                &
                     &                                                                             nodeBranch                                                                                          , &
                     &                                                                          +  states(stateCount)%self                    %densityContrast                                           &
                     &                                                                          *(                                                                                                       &
                     &                                                                            +states(stateCount)%cosmologyFunctions_%hubbleParameterEpochal(time           =basicBranch%time())     &
                     &                                                                            /states(stateCount)%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=            1.0d0 )     &
                     &                                                                           )                                                                                                  **2  &
                     &                                                                          *  states(stateCount)%cosmologyFunctions_%expansionFactor       (                basicBranch%time())**3, &
                     &                                                   cosmologyParameters_  =states(stateCount)%self%cosmologyParameters_                                                           , &
                     &                                                   cosmologyFunctions_   =states(stateCount)%self%cosmologyFunctions_                                                            , &
                     &                                                   virialDensityContrast_=states(stateCount)%self%virialDensityContrast_                                                           &
                     &                                                  )                                                                                                                             &
                     &              -massSiblings                                                                                                                                                     &
                     &             )                                                                                                                                                                  &
                     &            /(                                                                                                                                                                  &
                     &              +basicBranch%time()                                                                                                                                               &
                     &              -basicChild %time()                                                                                                                                               &
                     &             )
                massProgenitor   =+Dark_Matter_Profile_Mass_Definition(                                                                                                                               &
                     &                                                                           nodeChild                                                                                          , &
                     &                                                                        +  states(stateCount)%self                    %densityContrast                                          &
                     &                                                                        *(                                                                                                      &
                     &                                                                          +states(stateCount)%cosmologyFunctions_%hubbleParameterEpochal(time           =basicChild%time())     &
                     &                                                                          /states(stateCount)%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=           1.0d0 )     &
                     &                                                                         )                                                                                                 **2  &
                     &                                                                        *  states(stateCount)%cosmologyFunctions_%expansionFactor       (                basicChild%time())**3, &
                     &                                                 cosmologyParameters_  =states(stateCount)%self%cosmologyParameters_                                                          , &
                     &                                                 cosmologyFunctions_   =states(stateCount)%self%cosmologyFunctions_                                                           , &
                     &                                                 virialDensityContrast_=states(stateCount)%self%virialDensityContrast_                                                          &
                     &                                                )                                                                                                                               &
                     &            +massAccretionRate                                                                                                                                                  &
                     &            *(                                                                                                                                                                  &
                     &              +timeFormation                                                                                                                                                    &
                     &              -basicChild   %time()                                                                                                                                             &
                     &             )
                if (massProgenitor >= states(stateCount)%massLimit) &
                     & massBranch=+massBranch                       &
                     &            +massProgenitor
             end if
             nodeChild => nodeChild%sibling
          end do
       else if (.not.associated(nodeBranch%firstChild).and.basicBranch%time() == timeFormation) then
          massProgenitor=Dark_Matter_Profile_Mass_Definition(                                                                                                                                &
               &                                                                       nodeBranch                                                                                          , &
               &                                                                    +  states(stateCount)%self                    %densityContrast                                           &
               &                                                                    *(                                                                                                       &
               &                                                                      +states(stateCount)%cosmologyFunctions_%hubbleParameterEpochal(time           =basicBranch%time())     &
               &                                                                      /states(stateCount)%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=            1.0d0 )     &
               &                                                                     )                                                                                              **2      &
               &                                                                    *  states(stateCount)%cosmologyFunctions_%expansionFactor       (                basicBranch%time())**3, &
               &                                             cosmologyParameters_  =states(stateCount)%self%cosmologyParameters_                                                           , &
               &                                             cosmologyFunctions_   =states(stateCount)%self%cosmologyFunctions_                                                            , &
               &                                             virialDensityContrast_=states(stateCount)%self%virialDensityContrast_                                                           &
               &                                            )
          if (massProgenitor >= states(stateCount)%massLimit) &
               & massBranch=+massBranch                       &
               &            +massProgenitor
       end if
       nodeBranch => nodeBranch%firstChild
    end do
    ludlow2014FormationTimeRoot=+massBranch                                &
         &                      -states(stateCount)%massHaloCharacteristic
    return
  end function ludlow2014FormationTimeRoot

