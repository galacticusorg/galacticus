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

  !+    Contributions to this file made by: Yu Zhao

  !!{
  Implements a node operator class that accumulates the tidal-heating source term from the FDM solitonic core, following the model of \cite{Du_tidal_2018}.
  !!}

  !![
  <nodeOperator name="nodeOperatorTidalMassLossSoliton">
   <description>A node operator class that accumulates the tidal-heating source term from the FDM solitonic core, following the model of \cite{Du_tidal_2018}.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorTidalMassLossSoliton
     !!{
     A node operator class that accumulates the tidal-heating source term from the FDM solitonic core, following the model of \cite{Du_tidal_2018}.
     !!}
     private
     integer   :: massCoreNormalID, massCoreID, densityCoreID
   contains
     procedure :: differentialEvolution       => tidalMassLossSolitonDifferentialEvolution
     procedure :: differentialEvolutionScales => tidalMassLossSolitonDifferentialEvolutionScales
  end type nodeOperatorTidalMassLossSoliton
  
  interface nodeOperatorTidalMassLossSoliton
     !!{
     Constructors for the \refClass{nodeOperatorTidalMassLossSoliton} node operator class.
     !!}
     module procedure tidalMassLossSolitonConstructorParameters
     module procedure tidalMassLossSolitonConstructorInternal
  end interface nodeOperatorTidalMassLossSoliton
  
contains

  function tidalMassLossSolitonConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorTidalMassLossSoliton} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorTidalMassLossSoliton )                :: self
    type (inputParameters                  ), intent(inout) :: parameters

    self = nodeOperatorTidalMassLossSoliton()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function tidalMassLossSolitonConstructorParameters

  function tidalMassLossSolitonConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorTidalMassLossSoliton} node operator class.
    !!}
    implicit none
    type (nodeOperatorTidalMassLossSoliton )                        :: self
    !![
    <addMetaProperty component="darkMatterProfile" name="solitonMassCoreNormal" id="self%massCoreNormalID" isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCore"       id="self%massCoreID"       isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonDensityCore"    id="self%densityCoreID"    isEvolvable="no"  isCreator="no"/>
    !!]
    return
  end function tidalMassLossSolitonConstructorInternal

  subroutine tidalMassLossSolitonDifferentialEvolutionScales(self,node)
    !!{
    Set absolute ODE solver scale for the energy radiated from the hot halo due to cooling following the model of \cite{benson_galaxy_2010-1}.
    !!}
    implicit none
    class           (nodeOperatorTidalMassLossSoliton), intent(inout) :: self
    type            (treeNode                        ), intent(inout) :: node

    return
  end subroutine tidalMassLossSolitonDifferentialEvolutionScales

  subroutine tidalMassLossSolitonDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Accumulates an estimate of the energy radiated from the hot halo due to cooling following the model of \cite{benson_galaxy_2010-1}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDarkMatterProfile, nodeComponentSatellite, treeNode
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Constants_Astronomical, only : gigaYear                      , megaParsec           
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Vectors                         , only : Vector_Magnitude              , Vector_Product
    implicit none
    class           (nodeOperatorTidalMassLossSoliton), intent   (inout), target  :: self    
    type            (treeNode                        ), intent   (inout), target  :: node
    class           (nodeComponentSatellite          ), pointer                   :: satellite
    class           (nodeComponentDarkMatterProfile  ), pointer                   :: darkMatterProfile
    class           (massDistributionClass           ), pointer                   :: massDistribution_
    double precision                                  , dimension(3    )          :: position                      , velocity
    double precision                                  , parameter                 :: frequencyFractionalTiny=1.0d-6
    double precision                                  , parameter                 :: a=5.89794d-5                  , b=-8.72733d-2     , &
         &                                                                           c=1.6774d0                    , gamma=1.5d0
    double precision                                                              :: massSatellite                 , frequencyAngular  , &
         &                                                                           periodOrbital                 , radius            , &
         &                                                                           frequencyOrbital              , frequencyRadial   , &
         &                                                                           massHost                      , densityHost       , &
         &                                                                           densityCore                   , densityRatio      , &
         &                                                                           energyIm                      , massCore
    logical                                           , intent   (inout)          :: interrupt
    procedure       (interruptTask                   ), intent   (inout), pointer :: functionInterrupt
    integer                                           , intent   (in   )          :: propertyType
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType
    
    ! Get required quantities from the satellite node.
    if (.not.node%isSatellite()) return
    satellite          =>  node     %satellite (        )
    massSatellite      =   satellite%boundMass (        )
    position           =   satellite%position  (        )
    velocity           =   satellite%velocity  (        )
    radius             =   Vector_Magnitude    (position)
    ! Compute the orbital frequency. We use the larger of the angular and radial frequencies to avoid numerical problems for purely
    ! radial or purely circular orbits.
    frequencyAngular   =  +Vector_Magnitude(Vector_Product(position,velocity)) &
         &                /radius**2                                           &
         &                *kilo                                                &
         &                *gigaYear                                            &
         &                /megaParsec
    frequencyRadial    =  +abs             (   Dot_Product(position,velocity)) &
         &                /radius**2                                           &
         &                *kilo                                                &
         &                *gigaYear                                            &
         &                /megaParsec
    ! Find the orbital frequency. We use the larger of the angular and radial frequencies to avoid numerical problems for purely
    ! radial or purely circular orbits.
    frequencyOrbital   =max(                  &
         &                  frequencyAngular, &
         &                  frequencyRadial   &
         &                 )
    ! Find the orbital period.
    periodOrbital      =   +2.0d0             &
            &              *Pi                &
            &              /frequencyOrbital
    
    ! Get required quantities from the host node.
    massDistribution_  => node%parent%massDistribution()
    massHost=max(0.0d0,massDistribution_         %massEnclosedBySphere(radius))
    densityHost = 3.0d0*massHost/(4.0d0*Pi*radius**3)

    darkMatterProfile => node%darkMatterProfile()
    densityCore       = darkMatterProfile%floatRank0MetaPropertyGet(self%densityCoreID)
    massCore          = darkMatterProfile%floatRank0MetaPropertyGet(self%massCoreID)

    ! Compute the density ratio between the central density of the soliton and the average density of the host within the orbital radius
    densityRatio = densityCore/densityHost

    if (densityRatio>300.0d0) then
       densityRatio = 300.0d0
    end if

    ! Compute the imaginary part of the energy eigenvalue E, using the fitting formula, Equation (7) of Du et al. (2018; PRD; 97; 3507; https://ui.adsabs.harvard.edu/abs/2018PhRvD..97f3507D/abstract).
    energyIm = -exp(+a*(3.0d0*densityRatio/2.0d0/gamma)**2  &
         &          +b*(3.0d0*densityRatio/2.0d0/gamma)     &
         &          +c                            &
         &         )/periodOrbital
    
    ! Set the soliton core mass evolution rate following Du et al. (2018).
    ! From Eq. (17), the core mass obeys (1/M_c) dM_c/dt = (1/2) Im(E), where Im(E) is given by the fitting formula in Eq. (7).
    ! Here we set dM_c/dt = (1/2) Im(E) * M_c for time integration.
    call darkMatterProfile%floatRank0MetaPropertyRate(                       &
            &                                         self%massCoreNormalID, &
            &                                         +0.5d0                 &
            &                                         *energyIm              &
            &                                         *massCore              &
            &                                 )
    return
  end subroutine tidalMassLossSolitonDifferentialEvolution

