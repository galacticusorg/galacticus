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

  !+    Contributions to this file made by: Claude.

  !!{RST
  An implementation of dark matter halo profile concentrations using the :cite:t:`ludlow_mass-concentration-redshift_2014` algorithm.
  !!}

  !![
  <darkMatterProfileScaleRadius name="darkMatterProfileScaleRadiusLudlow2014" docformat="rst">
   <description>
   Dark matter halo scale radii are computed using the algorithm of :cite:t:`ludlow_mass-concentration-redshift_2014`. This
   differs from the :galacticus-class:`darkMatterProfileScaleRadiusLudlow2016` algorithm, from which it derives, only in the
   definition of the formation time: following :cite:t:`ludlow_mass-concentration-redshift_2014`, it is the time at which the
   mass of the main progenitor equals the characteristic mass, :math:`M_{-2}`, rather than that at which the collapsed mass
   in all progenitors more massive than :math:`f` times the halo mass does so. The parameter :math:`f` is used only to bound
   the search for the formation time.
   </description>
  </darkMatterProfileScaleRadius>
  !!]
  type, extends(darkMatterProfileScaleRadiusLudlow2016) :: darkMatterProfileScaleRadiusLudlow2014
     !!{RST
     A dark matter halo profile scale radii class implementing the algorithm of :cite:t:`ludlow_mass-concentration-redshift_2014`.
     !!}
     private
   contains
     procedure, nopass :: formationTimeRoot            => ludlow2014FormationTimeRoot
     procedure         :: formationTimeRootFunctionSet => ludlow2014FormationTimeRootFunctionSet
  end type darkMatterProfileScaleRadiusLudlow2014

  interface darkMatterProfileScaleRadiusLudlow2014
     !!{RST
     Constructors for the :galacticus-class:`darkMatterProfileScaleRadiusLudlow2014` dark matter halo profile scale radius class.
     !!}
     module procedure ludlow2014ConstructorParameters
     module procedure ludlow2014ConstructorInternal
  end interface darkMatterProfileScaleRadiusLudlow2014

contains

  function ludlow2014ConstructorParameters(parameters) result(self)
    !!{RST
    Default constructor for the ``ludlow2014`` dark matter halo profile concentration class.
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
    !!{RST
    Constructor for the :galacticus-class:`darkMatterProfileScaleRadiusLudlow2014` dark matter halo profile scale radius class.
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
    !!{RST
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
    !!{RST
    Function used to find the formation time of a halo in the ``ludlow2014`` concentration algorithm. Following
    :cite:t:`ludlow_mass-concentration-redshift_2014`, the formation time is that at which the mass of the main progenitor
    equals the characteristic mass, :math:`M_{-2}`. Within each step of the main branch the mass of the main progenitor is
    interpolated linearly in time, from the mass of the primary progenitor to the mass of the branch halo less those of the
    other progenitors, which merge with it only at the end of the step.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    double precision                    , intent(in   ) :: timeFormation
    type            (treeNode          ), pointer       :: nodeBranch     , nodeChild        , &
         &                                                 nodeSibling
    class           (nodeComponentBasic), pointer       :: basicBranch    , basicChild
    double precision                                    :: massProgenitor , massAccretionRate, &
         &                                                 massMainBranch

    nodeBranch     => states(stateCount)%node
    massMainBranch =  0.0d0
    do while (associated(nodeBranch))
       basicBranch => nodeBranch%basic()
       if (basicBranch%time() == timeFormation) then
          ! The trial formation time coincides with a node on the main branch, which is then the main progenitor.
          massMainBranch=ludlow2014MassDefined(nodeBranch)
          exit
       else if (associated(nodeBranch%firstChild) .and. basicBranch%time() > timeFormation) then
          nodeChild  => nodeBranch%firstChild
          basicChild => nodeChild %basic     ()
          if (basicChild%time() < timeFormation) then
             ! This step of the main branch straddles the trial formation time. The accretion rate is that of the primary
             ! progenitor alone, so the masses of all of its siblings are subtracted from that of the branch.
             massProgenitor    =  ludlow2014MassDefined(nodeChild )
             massAccretionRate =  ludlow2014MassDefined(nodeBranch) &
                  &              -massProgenitor
             nodeSibling       => nodeChild%sibling
             do while (associated(nodeSibling))
                massAccretionRate=+massAccretionRate                  &
                     &            -ludlow2014MassDefined(nodeSibling)
                nodeSibling      =>  nodeSibling%sibling
             end do
             massAccretionRate=+massAccretionRate    &
                  &            /(                    &
                  &              +basicBranch%time() &
                  &              -basicChild %time() &
                  &             )
             massMainBranch   =+massProgenitor               &
                  &            +massAccretionRate            &
                  &            *(                            &
                  &              +           timeFormation   &
                  &              -basicChild%time         () &
                  &             )
             exit
          end if
       end if
       nodeBranch => nodeBranch%firstChild
    end do
    ludlow2014FormationTimeRoot=+massMainBranch                            &
         &                      -states(stateCount)%massHaloCharacteristic
    return

  contains

    double precision function ludlow2014MassDefined(node)
      !!{RST
      Return the mass of the given node in the halo definition used by the ``ludlow2014`` concentration algorithm.
      !!}
      use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
      implicit none
      type (treeNode          ), intent(inout), target :: node
      class(nodeComponentBasic), pointer               :: basic

      basic                 => node%basic()
      ludlow2014MassDefined =  Dark_Matter_Profile_Mass_Definition(                                                                                                                          &
           &                                                                                 node                                                                                          , &
           &                                                                              +  states(stateCount)%self               %densityContrast                                          &
           &                                                                              *(                                                                                                 &
           &                                                                                +states(stateCount)%cosmologyFunctions_%hubbleParameterEpochal(time           =basic%time())     &
           &                                                                                /states(stateCount)%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=      1.0d0 )     &
           &                                                                               )                                                                                            **2  &
           &                                                                              *  states(stateCount)%cosmologyFunctions_%expansionFactor       (                basic%time())**3, &
           &                                                       cosmologyParameters_  =states(stateCount)%self%cosmologyParameters_                                                     , &
           &                                                       cosmologyFunctions_   =states(stateCount)%self%cosmologyFunctions_                                                      , &
           &                                                       virialDensityContrast_=states(stateCount)%self%virialDensityContrast_                                                     &
           &                                                      )
      return
    end function ludlow2014MassDefined

  end function ludlow2014FormationTimeRoot

