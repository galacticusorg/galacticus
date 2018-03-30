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

  !% An implementation of dark matter halo profile concentrations using the \cite{ludlow_mass-concentration-redshift_2014}
  !% algorithm.
  
  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationLudlow2014">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{ludlow_mass-concentration-redshift_2014}.</description>
  !# </darkMatterProfileConcentration>
  type, extends(darkMatterProfileConcentrationLudlow2016) :: darkMatterProfileConcentrationLudlow2014
     !% A dark matter halo profile concentration class implementing the algorithm of
     !% \cite{ludlow_mass-concentration-redshift_2014}.
     private
   contains
     !@ <objectMethods>
     !@   <object>darkMatterProfileConcentrationLudlow2014</object>
     !@   <objectMethod>
     !@     <method>formationTimeRoot</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ timeFormation\argin</arguments>
     !@     <description>Evalute a function which goes to zero at the formation time of the tree.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>formationTimeRootFunctionSet</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(rootFinder)\textgreater} finder\arginout</arguments>
     !@     <description>Initialize a root finder object for use in finding the formation time of the tree.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure, nopass :: formationTimeRoot            => ludlow2014FormationTimeRoot
     procedure         :: formationTimeRootFunctionSet => ludlow2014FormationTimeRootFunctionSet
  end type darkMatterProfileConcentrationLudlow2014

  interface darkMatterProfileConcentrationLudlow2014
     !% Constructors for the {\normalfont \ttfamily ludlow2014} dark matter halo profile concentration class.
     module procedure ludlow2014ConstructorParameters
     module procedure ludlow2014ConstructorInternal
  end interface darkMatterProfileConcentrationLudlow2014

contains

  function ludlow2014ConstructorParameters(parameters) result(self)
    !% Default constructor for the {\normalfont \ttfamily ludlow2014} dark matter halo profile concentration class.
    use Input_Parameters
    implicit none
    type(darkMatterProfileConcentrationLudlow2014)                :: self
    type(inputParameters                         ), intent(inout) :: parameters

    self%darkMatterProfileConcentrationLudlow2016=darkMatterProfileConcentrationLudlow2016(parameters)
    return
  end function ludlow2014ConstructorParameters
  
  function ludlow2014ConstructorInternal(C,f,timeFormationSeekDelta,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileConcentration_,darkMatterHaloScale_,darkMatterProfile_) result(self)
    !% Constructor for the {\normalfont \ttfamily ludlow2014} dark matter halo profile concentration class.
    use Galacticus_Error
    implicit none
    type            (darkMatterProfileConcentrationLudlow2014)                        :: self
    double precision                                          , intent(in   )         :: C                              , f, &
         &                                                                               timeFormationSeekDelta
    class           (cosmologyFunctionsClass                 ), intent(in   ), target :: cosmologyFunctions_
    class           (cosmologyParametersClass                ), intent(in   ), target :: cosmologyParameters_     
    class           (darkMatterProfileConcentrationClass     ), intent(in   ), target :: darkMatterProfileConcentration_
    class           (darkMatterHaloScaleClass                ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileClass                  ), intent(in   ), target :: darkMatterProfile_

    self%darkMatterProfileConcentrationLudlow2016=darkMatterProfileConcentrationLudlow2016(C,f,timeFormationSeekDelta,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileConcentration_,darkMatterHaloScale_,darkMatterProfile_)
    return
  end function ludlow2014ConstructorInternal

  subroutine ludlow2014FormationTimeRootFunctionSet(self,finder)
    !% Initialize the finder object to compute the relevant formation history.
    use Root_Finder
    implicit none
    class           (darkMatterProfileConcentrationLudlow2014), intent(inout) :: self
    type            (rootFinder                              ), intent(inout) :: finder
    double precision                                          , parameter     :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-4
    !GCC$ attributes unused :: self

    if (.not.finder%isInitialized()) then
       call finder%rootFunction(ludlow2014FormationTimeRoot        )
       call finder%tolerance   (toleranceAbsolute,toleranceRelative)
    end if
    return
  end subroutine ludlow2014FormationTimeRootFunctionSet
  
  double precision function ludlow2014FormationTimeRoot(timeFormation)
    !% Function used to find the formation time of a halo in the {\normalfont \ttfamily ludlow2014} concentration algorithm.
    use Dark_Matter_Profile_Mass_Definitions
    implicit none
    double precision                    , intent(in   ) :: timeFormation
    type            (treeNode          ), pointer       :: nodeBranch   , nodeChild        , &
         &                                                 nodeSibling
    class           (nodeComponentBasic), pointer       :: basicBranch  , basicChild       , &
         &                                                 basicSibling
    double precision                                    :: massBranch   , massAccretionRate, &
         &                                                 massSiblings , massProgenitor
    
    nodeBranch => ludlow2016States(ludlow2016StateCount)%node
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
                   massSiblings =  +massSiblings                                                                                                                              &
                        &          +Dark_Matter_Profile_Mass_Definition(                                                                                                      &
                        &                                                  nodeSibling                                                                                      , &
                        &                                               +  ludlow2016States(ludlow2016StateCount)%self                    %densityContrast                                                 &
                        &                                               *(                                                                                                    &
                        &                                                 +ludlow2016States(ludlow2016StateCount)%cosmologyFunctions_%hubbleParameterEpochal(time           =basicSibling%time())     &
                        &                                                 /ludlow2016States(ludlow2016StateCount)%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=             1.0d0 )     &
                        &                                                )                                                                                               **2  &
                        &                                               *  ludlow2016States(ludlow2016StateCount)%cosmologyFunctions_%expansionFactor       (                basicSibling%time())**3  &
                        &                                              )
                   nodeSibling  =>  nodeSibling %sibling
                end do
                massAccretionRate=+(                                                                                                                                          &
                     &              +Dark_Matter_Profile_Mass_Definition(                                                                                                     &
                     &                                                      nodeBranch                                                                                      , &
                     &                                                   +  ludlow2016States(ludlow2016StateCount)%self                    %densityContrast                                                &
                     &                                                   *(                                                                                                   &
                     &                                                     +ludlow2016States(ludlow2016StateCount)%cosmologyFunctions_%hubbleParameterEpochal(time           =basicBranch%time())     &
                     &                                                     /ludlow2016States(ludlow2016StateCount)%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=            1.0d0 )     &
                     &                                                    )                                                                                              **2  &
                     &                                                   *  ludlow2016States(ludlow2016StateCount)%cosmologyFunctions_%expansionFactor       (                basicBranch%time())**3  &
                     &                                                  )                                                                                                     &
                     &              -massSiblings                                                                                                                             &
                     &             )                                                                                                                                          &
                     &            /(                                                                                                                                          &
                     &              +basicBranch%time()                                                                                                                       &
                     &              -basicChild %time()                                                                                                                       &
                     &             )
                massProgenitor   =+Dark_Matter_Profile_Mass_Definition(                                                                                                    &
                     &                                                    nodeChild                                                                                      , &
                     &                                                 +  ludlow2016States(ludlow2016StateCount)%self                    %densityContrast                                               &
                     &                                                 *(                                                                                                  &
                     &                                                   +ludlow2016States(ludlow2016StateCount)%cosmologyFunctions_%hubbleParameterEpochal(time           =basicChild%time())     &
                     &                                                   /ludlow2016States(ludlow2016StateCount)%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=           1.0d0 )     &
                     &                                                  )                                                                                             **2  &
                     &                                                 *  ludlow2016States(ludlow2016StateCount)%cosmologyFunctions_%expansionFactor       (                basicChild%time())**3  &
                     &                                                )                                                                                                    &
                     &            +massAccretionRate                                                                                                                       &
                     &            *(                                                                                                                                       &
                     &              +timeFormation                                                                                                                         &
                     &              -basicChild   %time()                                                                                                                  &
                     &             )
                if (massProgenitor >= ludlow2016States(ludlow2016StateCount)%massLimit) &
                     & massBranch=+massBranch              &
                     &            +massProgenitor
             end if
             nodeChild => nodeChild%sibling
          end do
       else if (.not.associated(nodeBranch%firstChild).and.basicBranch%time() == timeFormation) then
          massProgenitor=Dark_Matter_Profile_Mass_Definition(                                                                                                     & 
               &                                                nodeBranch                                                                                      , &
               &                                             +  ludlow2016States(ludlow2016StateCount)%self                    %densityContrast                                                &
               &                                             *(                                                                                                   &
               &                                               +ludlow2016States(ludlow2016StateCount)%cosmologyFunctions_%hubbleParameterEpochal(time           =basicBranch%time())     &
               &                                               /ludlow2016States(ludlow2016StateCount)%cosmologyFunctions_%hubbleParameterEpochal(expansionFactor=            1.0d0 )     &
               &                                              )                                                                                              **2  &
               &                                             *  ludlow2016States(ludlow2016StateCount)%cosmologyFunctions_%expansionFactor       (                basicBranch%time())**3  &
               &                                            )
          if (massProgenitor >= ludlow2016States(ludlow2016StateCount)%massLimit) &
               & massBranch=+massBranch              &
               &            +massProgenitor
       end if
       nodeBranch => nodeBranch%firstChild
    end do
    ludlow2014FormationTimeRoot=+massBranch                       &
         &                      -ludlow2016States(ludlow2016StateCount)%massHaloCharacteristic
    return
  end function ludlow2014FormationTimeRoot
  
