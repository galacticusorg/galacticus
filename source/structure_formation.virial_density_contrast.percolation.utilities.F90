!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module of utilities needed by the {\normalfont \ttfamily percolation} virial density contrast class.

module Virial_Density_Contrast_Percolation_Utilities
  !% Provides utilities needed by the {\normalfont \ttfamily percolation} virial density contrast class.
  use Galacticus_Nodes                  , only : treeNode                                 , nodeComponentDarkMatterProfile
  use Dark_Matter_Profile_Scales        , only : darkMatterProfileScaleRadius             , darkMatterProfileScaleRadiusClass  , &
       &                                         darkMatterProfileScaleRadiusConcentration
  use Dark_Matter_Profiles_DMO          , only : darkMatterProfileDMO                     , darkMatterProfileDMOClass
  use Cosmology_Parameters              , only : cosmologyParameters                      , cosmologyParametersClass
  use Cosmology_Functions               , only : cosmologyFunctions                       , cosmologyFunctionsClass
  use Dark_Matter_Halo_Scales           , only : darkMatterHaloScale                      , darkMatterHaloScaleClass
  use Dark_Matter_Profiles_Concentration, only : darkMatterProfileConcentration           , darkMatterProfileConcentrationClass
  private
  public :: Virial_Density_Contrast_Percolation_Solver, Virial_Density_Contrast_Percolation_Objects_Constructor, percolationObjects, percolationObjectsDeepCopy

  ! Module-scope variables used in root finding.
  type            (treeNode                                 ), pointer :: workNode
  class           (nodeComponentDarkMatterProfile           ), pointer :: workDarkMatterProfile
  class           (darkMatterProfileDMOClass                ), pointer :: darkMatterProfileDMO_
  type            (darkMatterProfileScaleRadiusConcentration), pointer :: darkMatterProfileScaleRadius_
  double precision                                                     :: boundingDensity              , densityMatterMean, &
       &                                                                  massHalo
  double precision                                           , pointer :: densityContrast
  !$omp threadprivate(workNode,workDarkMatterProfile,boundingDensity,densityMatterMean,massHalo,densityContrast,darkMatterProfileScaleRadius_,darkMatterProfileDMO_)

  type :: percolationObjects
     !% Type used to store pointers to objects
     class(darkMatterProfileDMOClass          ), pointer :: darkMatterProfileDMO_
     class(cosmologyParametersClass           ), pointer :: cosmologyParameters_
     class(cosmologyFunctionsClass            ), pointer :: cosmologyFunctions_
     class(darkMatterHaloScaleClass           ), pointer :: darkMatterHaloScale_
     class(darkMatterProfileConcentrationClass), pointer :: darkMatterProfileConcentration_
   contains
     final :: percolationObjectsDestructor
  end type percolationObjects
  
contains

  !# <functionGlobal>
  !#  <unitName>Virial_Density_Contrast_Percolation_Objects_Constructor</unitName>
  !#  <type>class(*), pointer</type>
  !#  <module>Input_Parameters, only : inputParameter, inputParameters</module>
  !#  <arguments>type(inputParameters), intent(inout), target :: parameters</arguments>
  !# </functionGlobal>
  function Virial_Density_Contrast_Percolation_Objects_Constructor(parameters) result(self)
    !% Construct an instance of the container type for percolation virial density contrast objects from a parameter structure.
    use Input_Parameters, only : inputParameter, inputParameters
    implicit none
    class(*                         ), pointer               :: self
    type (percolationObjects        ), pointer               :: percolationObjects_
    type (inputParameters           ), intent(inout), target :: parameters
    
    allocate(percolationObjects_)
    !# <objectBuilder class="darkMatterProfileDMO"           name="percolationObjects_%darkMatterProfileDMO_"           source="parameters"/>
    !# <objectBuilder class="cosmologyParameters"            name="percolationObjects_%cosmologyParameters_"            source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"             name="percolationObjects_%cosmologyFunctions_"             source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"            name="percolationObjects_%darkMatterHaloScale_"            source="parameters"/>
    !# <objectBuilder class="darkMatterProfileConcentration" name="percolationObjects_%darkMatterProfileConcentration_" source="parameters"/>
    self => percolationObjects_
    nullify(percolationObjects_)
    return
  end function Virial_Density_Contrast_Percolation_Objects_Constructor

  subroutine percolationObjectsDestructor(self)
    !% Destruct an instance of the container type for percolation virial density contrast objects.
    implicit none
    type(percolationObjects), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterProfileDMO_"          />
    !# <objectDestructor name="self%cosmologyParameters_"           />
    !# <objectDestructor name="self%cosmologyFunctions_"            />
    !# <objectDestructor name="self%darkMatterHaloScale_"           />
    !# <objectDestructor name="self%darkMatterProfileConcentration_"/>
    return
  end subroutine percolationObjectsDestructor
  
  !# <functionGlobal>
  !#  <unitName>percolationObjectsDeepCopy</unitName>
  !#  <type>void</type>
  !#  <arguments>class(*), intent(inout) :: self, destination</arguments>
  !# </functionGlobal>
  subroutine percolationObjectsDeepCopy(self,destination)
    !% Perform a deep copy of percolation virial density contrast objects.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(*), intent(inout) :: self, destination

    select type (self)
    class is (percolationObjects)
       select type (destination)
       class is (percolationObjects)
          nullify(destination%darkMatterProfileDMO_          )
          nullify(destination%cosmologyParameters_           )
          nullify(destination%cosmologyFunctions_            )
          nullify(destination%darkMatterHaloScale_           )
          nullify(destination%darkMatterProfileConcentration_)
          allocate(destination%darkMatterProfileDMO_          ,mold=self%darkMatterProfileDMO_          )
          allocate(destination%cosmologyParameters_           ,mold=self%cosmologyParameters_           )
          allocate(destination%cosmologyFunctions_            ,mold=self%cosmologyFunctions_            )
          allocate(destination%darkMatterHaloScale_           ,mold=self%darkMatterHaloScale_           )
          allocate(destination%darkMatterProfileConcentration_,mold=self%darkMatterProfileConcentration_)
          !# <deepCopy source="self%darkMatterProfileDMO_"           destination="destination%darkMatterProfileDMO_"          />
          !# <deepCopy source="self%cosmologyParameters_"            destination="destination%cosmologyParameters_"           />
          !# <deepCopy source="self%cosmologyFunctions_"             destination="destination%cosmologyFunctions_"            />
          !# <deepCopy source="self%darkMatterHaloScale_"            destination="destination%darkMatterHaloScale_"           />
          !# <deepCopy source="self%darkMatterProfileConcentration_" destination="destination%darkMatterProfileConcentration_"/>
       class default
          call Galacticus_Error_Report("destination must be of 'percolationObjects' class"//{introspection:location})
       end select
    class default
       call Galacticus_Error_Report("self must be of 'percolationObjects' class"//{introspection:location})
    end select
    return
  end subroutine percolationObjectsDeepCopy
  
  !# <functionGlobal>
  !#  <unitName>Virial_Density_Contrast_Percolation_Solver</unitName>
  !#  <type>double precision</type>
  !#  <arguments>double precision                            , intent(in   )         :: mass, time, linkingLength</arguments>
  !#  <arguments>double precision                            , intent(in   ), target :: densityContrastCurrent</arguments>
  !#  <arguments>class           (*                         ), intent(in   )         :: percolationObjects_</arguments>
  !#  <arguments>class           (*                         ), intent(inout)         :: virialDensityContrast_</arguments>
  !# </functionGlobal>
  double precision function Virial_Density_Contrast_Percolation_Solver(mass,time,linkingLength,densityContrastCurrent,percolationObjects_,virialDensityContrast_)
    !% Return the virial density contrast at the given epoch, based on the percolation algorithm of \cite{more_overdensity_2011}.
    use Root_Finder                   , only : rootFinder                    , rangeExpandMultiplicative
    use Galacticus_Error              , only : Galacticus_Error_Report
    use Numerical_Constants_Math      , only : Pi
    use Galacticus_Calculations_Resets, only : Galacticus_Calculations_Reset
    use Galacticus_Nodes              , only : nodeComponentBasic
    use Virial_Density_Contrast       , only : virialDensityContrastClass
   implicit none
    double precision                                     , intent(in   )         :: mass                                      , time, &
         &                                                                          linkingLength
    double precision                                     , intent(in   ), target :: densityContrastCurrent
    class           (*                                  ), intent(in   )         :: percolationObjects_
    class           (*                                  ), intent(inout)         :: virialDensityContrast_
    double precision                                     , parameter             :: percolationThreshold           =0.652960d0
    class           (cosmologyParametersClass           ), pointer               :: cosmologyParameters_
    class           (cosmologyFunctionsClass            ), pointer               :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass           ), pointer               :: darkMatterHaloScale_
    class           (darkMatterProfileConcentrationClass), pointer               :: darkMatterProfileConcentration_
    class           (nodeComponentBasic                 ), pointer               :: workBasic
    type            (rootFinder                         )                        :: finder
    double precision                                                             :: radiusHalo

    ! Initialize module-scope variables.
    massHalo        =  mass
    densityContrast => densityContrastCurrent
    ! Extract required objects from the container.
    select type (percolationObjects_)
    type is (percolationObjects)
       darkMatterProfileDMO_           => percolationObjects_%darkMatterProfileDMO_
       cosmologyFunctions_             => percolationObjects_%cosmologyFunctions_
       cosmologyParameters_            => percolationObjects_%cosmologyParameters_
       darkMatterHaloScale_            => percolationObjects_%darkMatterHaloScale_
       darkMatterProfileConcentration_ => percolationObjects_%darkMatterProfileConcentration_
       class default
       call Galacticus_Error_Report('percolationObjects_ must be of "percolationObjects" type'//{introspection:location})
    end select
    ! Build a scale radius object.
    allocate(darkMatterProfileScaleRadius_)
    select type (virialDensityContrast_)
    class is (virialDensityContrastClass)
       !# <referenceConstruct object="darkMatterProfileScaleRadius_">
       !#  <constructor>
       !#   darkMatterProfileScaleRadiusConcentration(                                                                   &amp;
       !#    &amp;                                    correctForConcentrationDefinition=.true.                         , &amp;
       !#    &amp;                                    useMeanConcentration             =.true.                         , &amp;
       !#    &amp;                                    cosmologyParameters_             =cosmologyParameters_           , &amp;
       !#    &amp;                                    cosmologyFunctions_              =cosmologyFunctions_            , &amp;
       !#    &amp;                                    darkMatterHaloScale_             =darkMatterHaloScale_           , &amp;
       !#    &amp;                                    darkMatterProfileDMO_            =darkMatterProfileDMO_          , &amp;
       !#    &amp;                                    virialDensityContrast_           =virialDensityContrast_         , &amp;
       !#    &amp;                                    darkMatterProfileConcentration_  =darkMatterProfileConcentration_  &amp;
       !#    &amp;                                   )
       !#  </constructor>
       !# </referenceConstruct>
    class default
       call Galacticus_Error_Report('virialDensityContrast_ must be of "virialDensityContrastClass" class'//{introspection:location})
    end select
    ! Compute the bounding density, based on percolation theory (eq. 5 of More et al.).
    densityMatterMean=cosmologyFunctions_%matterDensityEpochal(time)
    boundingDensity=+densityMatterMean       &
         &          *percolationThreshold    &
         &          /linkingLength       **3
    ! Create a node and set the mass and time.
    workNode              => treeNode                  (                 )
    workBasic             => workNode%basic            (autoCreate=.true.)
    workDarkMatterProfile => workNode%darkMatterProfile(autoCreate=.true.)
    call workBasic            %massSet            (mass   )
    call workBasic            %timeSet            (time   )
    call workBasic            %timeLastIsolatedSet(time   )
    call workDarkMatterProfile%scaleIsLimitedSet  (.false.)
    call Galacticus_Calculations_Reset(workNode)
    ! Make an initial guess at the halo radius.
    radiusHalo=(mass/4.0d0/Pi/boundingDensity)**(1.0d0/3.0d0)
    ! Find the corresponding halo radius.
    call finder   %tolerance          (                                               &
         &                             toleranceRelative  =1.0d-3                     &
         &                            )
    call finder   %rangeExpand        (                                               &
         &                             rangeExpandUpward  =2.0d0                    , &
         &                             rangeExpandDownward=0.5d0                    , &
         &                             rangeExpandType    =rangeExpandMultiplicative  &
         &                            )
    call finder   %rootFunction       (                                               &
         &                                                 haloRadiusRootFunction     &
         &                            )
    radiusHalo=finder%find(rootGuess=radiusHalo)
    call workNode%destroy()
    deallocate(workNode)
    !# <objectDestructor name="darkMatterProfileScaleRadius_"/>
    ! Compute the corresponding density contrast.    
    Virial_Density_Contrast_Percolation_Solver=+3.0d0                &
         &                                     *mass                 &
         &                                     /4.0d0                &
         &                                     /Pi                   &
         &                                     /radiusHalo       **3 &
         &                                     /densityMatterMean
    return
  end function Virial_Density_Contrast_Percolation_Solver

  double precision function haloRadiusRootFunction(haloRadiusTrial)
    !% Root function used to find the radius of a halo giving the correct bounding density.
    use Dark_Matter_Profiles_Shape    , only : darkMatterProfileShape       , darkMatterProfileShapeClass
    use Numerical_Constants_Math      , only : Pi
    use Galacticus_Calculations_Resets, only : Galacticus_Calculations_Reset
    implicit none
    double precision                             , intent(in   ) :: haloRadiusTrial
    double precision                                             :: scaleRadius            , densityHaloRadius
    class           (darkMatterProfileShapeClass), pointer       :: darkMatterProfileShape_

    ! Construct the current density contrast.
    densityContrast=+3.0d0                &
         &          *massHalo             &
         &          /4.0d0                &
         &          /Pi                   &
         &          /haloRadiusTrial  **3 &
         &          /densityMatterMean
    ! Find scale radius of the halo.
    scaleRadius=darkMatterProfileScaleRadius_%radius(workNode)
    call workDarkMatterProfile%scaleSet(scaleRadius)
    if (workDarkMatterProfile%shapeIsSettable()) then
       darkMatterProfileShape_ => darkMatterProfileShape()
       call workDarkMatterProfile%shapeSet(darkMatterProfileShape_%shape(workNode))
    end if
    call Galacticus_Calculations_Reset(workNode)
    ! Compute density at the halo radius.
    densityHaloRadius=darkMatterProfileDMO_%density(workNode,haloRadiusTrial)
    ! Find difference from target density.
    haloRadiusRootFunction=boundingDensity-densityHaloRadius
    return
  end function haloRadiusRootFunction

end module Virial_Density_Contrast_Percolation_Utilities
