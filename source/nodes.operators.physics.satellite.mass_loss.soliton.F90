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

  !!{
  Implements a node operator class that accumulates the tidal-heating source term from the FDM solitonic core, following the model of \cite{Du_tidal_2018}.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <nodeOperator name="nodeOperatorSatelliteSolitonMassLoss">
   <description>A node operator class that accumulates the tidal-heating source term from the FDM solitonic core, following the model of \cite{Du_tidal_2018}.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteSolitonMassLoss
     !!{
     A node operator class that accumulates the tidal-heating source term from the FDM solitonic core, following the model of \cite{Du_tidal_2018}.
     !!}
     private
     class  (darkMatterHaloScaleClass               ), pointer :: darkMatterHaloScale_ => null()
     integer                                                   :: densityCoreEvolutionID
     integer                                                   :: densityCoreID
   contains
     final     ::                                tidalMassLossSolitonDestructor
     procedure :: differentialEvolution       => tidalMassLossSolitonDifferentialEvolution
     procedure :: differentialEvolutionScales => tidalMassLossSolitonDifferentialEvolutionScales
     procedure :: galaxiesMerge               => tidalMassLossSolitonGalaxiesMerge
  end type nodeOperatorSatelliteSolitonMassLoss
  
  interface nodeOperatorSatelliteSolitonMassLoss
     !!{
     Constructors for the \refClass{nodeOperatorSatelliteSolitonMassLoss} node operator class.
     !!}
     module procedure tidalMassLossSolitonConstructorParameters
     module procedure tidalMassLossSolitonConstructorInternal
  end interface nodeOperatorSatelliteSolitonMassLoss
  
contains

  function tidalMassLossSolitonConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSatelliteSolitonMassLoss} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorSatelliteSolitonMassLoss )                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass             ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
     self=nodeOperatorSatelliteSolitonMassLoss(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function tidalMassLossSolitonConstructorParameters

  function tidalMassLossSolitonConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorSatelliteSolitonMassLoss} node operator class.
    !!}
    implicit none
    type (nodeOperatorSatelliteSolitonMassLoss)                        :: self
    class(darkMatterHaloScaleClass            ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    <addMetaProperty component="basic" name="densityCoreEvolution" id="self%densityCoreEvolutionID" isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="basic" name="densityCore"          id="self%densityCoreID"          isEvolvable="no" isCreator="no"/>
    !!]
    return
  end function tidalMassLossSolitonConstructorInternal

  subroutine tidalMassLossSolitonDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorSatelliteSolitonMassLoss} node operator class.
    !!}
    implicit none
    type(nodeOperatorSatelliteSolitonMassLoss), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine tidalMassLossSolitonDestructor

  subroutine tidalMassLossSolitonDifferentialEvolutionScales(self,node)
    !!{
    Set absolute ODE solver scale for the energy radiated from the hot halo due to cooling following the model of \cite{benson_galaxy_2010-1}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic
    implicit none
    class           (nodeOperatorSatelliteSolitonMassLoss), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , parameter     :: scaleRelative     =+1.0d-1
    class           (nodeComponentBasic                  ), pointer       :: basic

    basic=>node%basic()
    call basic%floatRank0MetaPropertyScale(self%densityCoreEvolutionID, scaleRelative)

    return
  end subroutine tidalMassLossSolitonDifferentialEvolutionScales

  subroutine tidalMassLossSolitonDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Accumulates an estimate of the energy radiated from the hot halo due to cooling following the model of \cite{benson_galaxy_2010-1}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic    , nodeComponentSatellite, treeNode
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Constants_Astronomical, only : gigaYear              , megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Vectors                         , only : Vector_Magnitude      , Vector_Product
    implicit none
    class           (nodeOperatorSatelliteSolitonMassLoss), intent   (inout), target  :: self    
    type            (treeNode                            ), intent   (inout), target  :: node
    logical                                               , intent   (inout)          :: interrupt
    procedure       (interruptTask                       ), intent   (inout), pointer :: functionInterrupt
    integer                                               , intent   (in   )          :: propertyType
    class           (nodeComponentSatellite              ), pointer                   :: satellite
    class           (nodeComponentBasic                  ), pointer                   :: basic
    class           (massDistributionClass               ), pointer                   :: massDistribution_
    double precision                                      , dimension(3    )          :: position                      , velocity
    double precision                                      , parameter                 :: frequencyFractionalTiny=1.0d-6
    double precision                                      , parameter                 :: a=5.89794d-5                  , b=-8.72733d-2   , &
         &                                                                               c=1.6774d0                    , gamma=1.5d0
    double precision                                                                  :: frequencyAngular              , frequencyRadial , &
         &                                                                               periodOrbital                 , radius          , &
         &                                                                               densityCore                   , densityHost     , &
         &                                                                               massHost                      , energyIm        , &
         &                                                                               frequencyOrbital              , densityRatio
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    ! Get required quantities from the satellite node.
    satellite          =>  node     %satellite (        )
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
    frequencyOrbital   =  max(                  &
         &                    frequencyAngular, &
         &                    frequencyRadial   &
         &                   )
    ! Find the orbital period.
    periodOrbital      =  +2.0d0                &
         &                *Pi                   &
         &                /frequencyOrbital

    massDistribution_  => node%massDistribution()
    massHost           =  max(0.0d0,massDistribution_%massEnclosedBySphere(radius))
    densityHost        =  3.0d0*massHost/(4.0d0*Pi*radius**3)

    basic              => node%basic           ()
    densityCore        =  basic%floatRank0MetaPropertyGet(self%densityCoreID)

    ! Compute the density ratio between the central density of the soliton and the average density of the host within the orbital radius
    densityRatio       =  densityCore/densityHost

    if (densityRatio>300.0d0) then
        densityRatio = 300.0d0
    end if

    ! Compute the imaginary part of the energy eigenvalue E, using the fitting formula, Equation (7) of Du et al. (2018; PRD; 97; 3507; https://ui.adsabs.harvard.edu/abs/2018PhRvD..97f3507D/abstract).
    energyIm           =  -exp(+a*(3.0d0*densityRatio/2.0d0/gamma)**2  &
         &                     +b*(3.0d0*densityRatio/2.0d0/gamma)     &
         &                     +c                                      &
         &                    )/periodOrbital
    
    ! Set the rate to be integrated by Galacticus for the soliton core density evolution. From Du et al. (2018), equation (17) can be written as: d(log ρ_c) / dt = 2 Im(E)
    ! Im(E) is given by their fitting formula (eq. 7) as a function of densityRatio = ρ_c / ρ̄_host(<r_orbit) and the orbital period T_orbit.
    ! Here we evaluate the instantaneous rate 2*Im(E), which will be integrated over time to obtain Δlog ρ_c(t) = ∫ [2 Im(E)] dt.
    call basic%floatRank0MetaPropertyRate(                             &
         &                                self%densityCoreEvolutionID, &
         &                                +2.0d0                       &
         &                                *energyIm                    &
         &                               )
    return
  end subroutine tidalMassLossSolitonDifferentialEvolution

  subroutine tidalMassLossSolitonGalaxiesMerge(self,node)
    !!{
    Zero the soliton core density evolution rate for a subhalo that is about to merge.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (nodeOperatorSatelliteSolitonMassLoss), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    class           (nodeComponentBasic                  ), pointer       :: basic

    basic=>node%basic()
    call basic%floatRank0MetaPropertySet(                             &
         &                               self%densityCoreEvolutionID, &
         &                               +0.0d0                       &
         &                              )
    return
  end subroutine tidalMassLossSolitonGalaxiesMerge
