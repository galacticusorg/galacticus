!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  A dark matter halo profile heating class which accounts for heating arising from impulsive outflows.
  !!}
  
  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingImpulsiveOutflow">
   <description>
    A dark matter profile heating model which accounts for heating due to impulsive outflows---i.e. outflows occuring on
    timescales that are small relative to the dynamical time of the halo. The model assumed is that the energy injection is given by
    \begin{equation}
    \dot{\epsilon}(r) = \alpha \frac{\mathrm{G} \dot{M}_\mathrm{outflow}(r)}{r} f\left( \frac{t_\phi}{t_\mathrm{dyn}} \right),
    \end{equation}
    where $\alpha$ is a normalization factor, $t_\phi = M_\mathrm{gas}/\dot{M}_\mathrm{outflow}$ is the timescale for the
    outflow, and $t_\mathrm{dyn} = r_{1/2}/v_{1/2}$ is the dynamical time at the half-mass radius.
    
    In practice, the quantity
    \begin{equation}
    \dot{\epsilon}^\prime = \dot{M}_\mathrm{outflow} f\left( \frac{t_\phi}{t_\mathrm{dyn}} \right),
    \end{equation}
    has been accumulated by the \refClass{nodeOperatorImpulsiveOutflowEnergy} object---radially-dependent factors are then applied here.
   </description>
  </darkMatterProfileHeating>
  !!]
  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingImpulsiveOutflow
     !!{
     A dark matter profile heating class which accounts for heating arising from impulsive outflows.
     !!}
     private
     integer                      :: energyImpulsiveOutflowDiskID          , energyImpulsiveOutflowSpheroidID
     double precision             :: impulsiveEnergyFactor
     class           (*), pointer :: galacticStructure_           => null()
   contains
     final     ::                                    impulsiveOutflowDestructor
     procedure :: specificEnergy                  => impulsiveOutflowSpecificEnergy
     procedure :: specificEnergyGradient          => impulsiveOutflowSpecificEnergyGradient
     procedure :: specificEnergyIsEverywhereZero  => impulsiveOutflowSpecificEnergyIsEverywhereZero
     procedure :: deepCopyReset                   => impulsiveOutflowDeepCopyReset
     procedure :: deepCopy                        => impulsiveOutflowDeepCopy
     procedure :: deepCopyFinalize                => impulsiveOutflowDeepCopyFinalize
  end type darkMatterProfileHeatingImpulsiveOutflow

  interface darkMatterProfileHeatingImpulsiveOutflow
     !!{
     Constructors for the {\normalfont \ttfamily impulsiveOutflow} dark matter profile heating class.
     !!}
     module procedure impulsiveOutflowConstructorParameters
     module procedure impulsiveOutflowConstructorInternal
  end interface darkMatterProfileHeatingImpulsiveOutflow

contains

  function impulsiveOutflowConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily impulsiveOutflow} dark matter profile heating scales class which takes a parameter set as input.
    !!}
    use :: Functions_Global, only : galacticStructureConstruct_, galacticStructureDestruct_
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (darkMatterProfileHeatingImpulsiveOutflow), target        :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (*                                       ), pointer       :: galacticStructure_
    double precision                                                          :: impulsiveEnergyFactor
    !![
    <inputParameter>
      <name>impulsiveEnergyFactor</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The parameter $\alpha$ appearing in the impulsive outflow heating rate.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    call galacticStructureConstruct_(parameters,galacticStructure_)
    self=darkMatterProfileHeatingImpulsiveOutflow(impulsiveEnergyFactor,galacticStructure_)
    !![
    <inputParametersValidate source="parameters" extraAllowedNames="galacticStructure"/>
    !!]
    call galacticStructureDestruct_(galacticStructure_)
    return
  end function impulsiveOutflowConstructorParameters

  function impulsiveOutflowConstructorInternal(impulsiveEnergyFactor,galacticStructure_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily impulsiveOutflow} dark matter profile heating scales class.
    !!}
    implicit none
    type            (darkMatterProfileHeatingImpulsiveOutflow)                        :: self
    class           (*                                       ), intent(in   ), target :: galacticStructure_
    double precision                                          , intent(in   )         :: impulsiveEnergyFactor
    !![
    <constructorAssign variables="impulsiveEnergyFactor, *galacticStructure_"/>
    !!]

    !![
    <addMetaProperty component="darkMatterProfile" name="energyImpulsiveOutflowDisk"     id="self%energyImpulsiveOutflowDiskID"     isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="energyImpulsiveOutflowSpheroid" id="self%energyImpulsiveOutflowSpheroidID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function impulsiveOutflowConstructorInternal

  subroutine impulsiveOutflowDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily impulsiveOutflow} node operator class.
    !!}
    use :: Functions_Global, only : galacticStructureDestruct_
    implicit none
    type(darkMatterProfileHeatingImpulsiveOutflow), intent(inout) :: self

    if (associated(self%galacticStructure_)) call galacticStructureDestruct_(self%galacticStructure_)
    return
  end subroutine impulsiveOutflowDestructor
  
  double precision function impulsiveOutflowSpecificEnergy(self,node,radius,darkMatterProfileDMO_)
    !!{
    Returns the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDarkMatterProfile
    use :: Galactic_Structure_Options      , only : componentTypeDisk             , componentTypeSpheroid, radiusLarge,massTypeDark 
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Functions_Global                , only : galacticStructureMassEnclosed_
    implicit none
    class           (darkMatterProfileHeatingImpulsiveOutflow), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    class           (darkMatterProfileDMOClass               ), intent(inout) :: darkMatterProfileDMO_
    double precision                                          , intent(in   ) :: radius
    class           (nodeComponentDarkMatterProfile          ), pointer       :: darkMatterProfile
    double precision                                                          :: massTotalDisk        , massTotalSpheroid    , &
         &                                                                       fractionMassDisk     , fractionMassSpheroid
    
    darkMatterProfile => node%darkMatterProfile        (                                                                            )
    massTotalDisk     =  galacticStructureMassEnclosed_(self%galacticStructure_,node,radiusLarge,componentType=componentTypeDisk    )
    massTotalSpheroid =  galacticStructureMassEnclosed_(self%galacticStructure_,node,radiusLarge,componentType=componentTypeSpheroid)
    if (massTotalDisk     > 0.0d0) then
       fractionMassDisk    =+galacticStructureMassEnclosed_(self%galacticStructure_,node,radius,componentType=componentTypeDisk    ) &
            &               /massTotalDisk
    else
       fractionMassDisk    =+0.0d0
    end if
    if (massTotalSpheroid > 0.0d0) then
       fractionMassSpheroid=+galacticStructureMassEnclosed_(self%galacticStructure_,node,radius,componentType=componentTypeSpheroid) &
            &               /massTotalSpheroid
    else
       fractionMassSpheroid=+0.0d0
    end if
    impulsiveOutflowSpecificEnergy =  +self                           %impulsiveEnergyFactor                                            &
         &                            *gravitationalConstantGalacticus                                                                  &
         &                            *(                                                                                                &
         &                             +darkMatterProfile             %floatRank0MetaPropertyGet(self%energyImpulsiveOutflowDiskID    ) &
         &                             *fractionMassDisk                                                                                &
         &                             +darkMatterProfile             %floatRank0MetaPropertyGet(self%energyImpulsiveOutflowSpheroidID) &
         &                             *fractionMassSpheroid                                                                            &
         &                            )                                                                                                 &
         &                            /radius
    return
  end function impulsiveOutflowSpecificEnergy

  double precision function impulsiveOutflowSpecificEnergyGradient(self,node,radius,darkMatterProfileDMO_)
    !!{
    Returns the gradient of the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDarkMatterProfile
    use :: Galactic_Structure_Options      , only : componentTypeDisk              , componentTypeSpheroid, radiusLarge, coordinateSystemSpherical
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Math        , only : Pi
    use :: Functions_Global                , only : galacticStructureMassEnclosed_ , galacticStructureDensity_
    implicit none
    class           (darkMatterProfileHeatingImpulsiveOutflow), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    class           (darkMatterProfileDMOClass               ), intent(inout) :: darkMatterProfileDMO_
    double precision                                          , intent(in   ) :: radius
    class           (nodeComponentDarkMatterProfile          ), pointer       :: darkMatterProfile
    double precision                                                          :: massTotalDisk        , massTotalSpheroid      , &
         &                                                                       fractionMassDisk     , fractionMassSpheroid   , &
         &                                                                       fractionDensityDisk  , fractionDensitySpheroid

    darkMatterProfile => node%darkMatterProfile        (                                                                            )
    massTotalDisk     =  galacticStructureMassEnclosed_(self%galacticStructure_,node,radiusLarge,componentType=componentTypeDisk    )
    massTotalSpheroid =  galacticStructureMassEnclosed_(self%galacticStructure_,node,radiusLarge,componentType=componentTypeSpheroid)
    if (massTotalDisk     > 0.0d0) then
       fractionMassDisk       =+galacticStructureMassEnclosed_(self%galacticStructure_,node, radius             ,componentType=componentTypeDisk                                               ) &
            &                  /massTotalDisk
       fractionDensityDisk    =+galacticStructureDensity_     (self%galacticStructure_,node,[radius,0.0d0,0.0d0],componentType=componentTypeDisk    ,coordinateSystem=coordinateSystemSpherical) &
            &                  /massTotalDisk
    else
       fractionMassDisk       =+0.0d0
       fractionDensityDisk    =+0.0d0
    end if
    if (massTotalSpheroid > 0.0d0) then
       fractionMassSpheroid   =+galacticStructureMassEnclosed_(self%galacticStructure_,node, radius             ,componentType=componentTypeSpheroid                                           ) &
            &                  /massTotalSpheroid
       fractionDensitySpheroid=+galacticStructureDensity_     (self%galacticStructure_,node,[radius,0.0d0,0.0d0],componentType=componentTypeSpheroid,coordinateSystem=coordinateSystemSpherical) &
            &                  /massTotalSpheroid
    else
       fractionMassSpheroid   =+0.0d0
       fractionDensitySpheroid=+0.0d0
    end if
    impulsiveOutflowSpecificEnergyGradient =  +self%impulsiveEnergyFactor                                                                          &
         &                                    *gravitationalConstantGalacticus                                                                     &
         &                                    *(                                                                                                   &
         &                                      +(                                                                                                 &
         &                                        +darkMatterProfile             %floatRank0MetaPropertyGet(self%energyImpulsiveOutflowDiskID    ) &
         &                                        *fractionDensityDisk                                                                             &
         &                                        +darkMatterProfile             %floatRank0MetaPropertyGet(self%energyImpulsiveOutflowSpheroidID) &
         &                                        *fractionDensitySpheroid                                                                         &
         &                                       )                                                                                                 &
         &                                      *4.0d0                                                                                             &
         &                                      *Pi                                                                                                &
         &                                      *radius                                                                                            &
         &                                      -(                                                                                                 &
         &                                        +darkMatterProfile             %floatRank0MetaPropertyGet(self%energyImpulsiveOutflowDiskID    ) &
         &                                        *fractionMassDisk                                                                                &
         &                                        +darkMatterProfile             %floatRank0MetaPropertyGet(self%energyImpulsiveOutflowSpheroidID) &
         &                                        *fractionMassSpheroid                                                                            &
         &                                       )                                                                                                 &
         &                                      /radius**2                                                                                         &
         &                                     )
    return
  end function impulsiveOutflowSpecificEnergyGradient

  logical function impulsiveOutflowSpecificEnergyIsEverywhereZero(self,node,darkMatterProfileDMO_)
    !!{
    Returns true if the specific energy is everywhere zero in the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(darkMatterProfileHeatingImpulsiveOutflow), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node
    class(darkMatterProfileDMOClass               ), intent(inout) :: darkMatterProfileDMO_
    class(nodeComponentDarkMatterProfile          ), pointer       :: darkMatterProfile
    !$GLC attributes unused :: darkMatterProfileDMO_

    darkMatterProfile                              =>  node             %darkMatterProfile        (                                     )
    impulsiveOutflowSpecificEnergyIsEverywhereZero =   darkMatterProfile%floatRank0MetaPropertyGet(self%energyImpulsiveOutflowDiskID    ) <= 0.0d0 &
         &                                            .and.                                                                                        &
         &                                             darkMatterProfile%floatRank0MetaPropertyGet(self%energyImpulsiveOutflowSpheroidID) <= 0.0d0

    return
  end function impulsiveOutflowSpecificEnergyIsEverywhereZero

  subroutine impulsiveOutflowDeepCopyReset(self)
    !!{
    Perform a deep copy reset of the object.
    !!}
    use :: Functions_Global, only : galacticStructureDeepCopyReset_
    implicit none
    class(darkMatterProfileHeatingImpulsiveOutflow), intent(inout) :: self
    
    self%copiedSelf => null()
    if (associated(self%galacticStructure_)) call galacticStructureDeepCopyReset_(self%galacticStructure_)
    return
  end subroutine impulsiveOutflowDeepCopyReset
  
  subroutine impulsiveOutflowDeepCopyFinalize(self)
    !!{
    Finalize a deep reset of the object.
    !!}
    use :: Functions_Global, only : galacticStructureDeepCopyFinalize_
    implicit none
    class(darkMatterProfileHeatingImpulsiveOutflow), intent(inout) :: self

   if (associated(self%galacticStructure_)) call galacticStructureDeepCopyFinalize_(self%galacticStructure_)
    return
  end subroutine impulsiveOutflowDeepCopyFinalize
  
  subroutine impulsiveOutflowDeepCopy(self,destination)
    !!{
    Perform a deep copy of the object.
    !!}
    use :: Error           , only : Error_Report
    use :: Functions_Global, only : galacticStructureDeepCopy_
    implicit none
    class(darkMatterProfileHeatingImpulsiveOutflow), intent(inout) :: self
    class(darkMatterProfileHeatingClass           ), intent(inout) :: destination

    call self%darkMatterProfileHeatingClass%deepCopy(destination)
    select type (destination)
    type is (darkMatterProfileHeatingImpulsiveOutflow)
       destination%energyImpulsiveOutflowDiskID    =self%energyImpulsiveOutflowDiskID
       destination%energyImpulsiveOutflowSpheroidID=self%energyImpulsiveOutflowSpheroidID
       destination%impulsiveEnergyFactor           =self%impulsiveEnergyFactor
       nullify(destination%galacticStructure_)
       if (associated(self%galacticStructure_)) then
          allocate(destination%galacticStructure_,mold=self%galacticStructure_)
          call galacticStructureDeepCopy_(self%galacticStructure_,destination%galacticStructure_)
          destination%referenceCount=1          
       end if
    class default
       call Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine impulsiveOutflowDeepCopy
