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

  !!{RST
  Implements a node operator class that accumulates the tidal-heating source term from the FDM solitonic core, following the model of :cite:t:`du_tidal_2018`.
  !!}

  !![
  <satelliteTidalStripping name="satelliteTidalStrippingDu2018" docformat="rst">
   <description>
   A node operator class that accumulates the tidal-heating source term from the FDM solitonic core, following the model of :cite:t:`du_tidal_2018`.
   </description>
  </satelliteTidalStripping>
  !!]
  type, extends(satelliteTidalStrippingClass) :: satelliteTidalStrippingDu2018
     !!{RST
     A node operator class that accumulates the tidal-heating source term from the FDM solitonic core, following the model of :cite:t:`du_tidal_2018`.
     !!}
     private
     integer :: massCoreID, densityCoreID
   contains
     procedure :: massLossRate => du2018MassLossRate
  end type satelliteTidalStrippingDu2018

  interface satelliteTidalStrippingDu2018
     !!{RST
     Constructors for the :galacticus-class:`satelliteTidalStrippingDu2018` satellite tidal stripping class.
     !!}
     module procedure du2018ConstructorParameters
     module procedure du2018ConstructorInternal
  end interface satelliteTidalStrippingDu2018

contains

  function du2018ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`satelliteTidalStrippingDu2018` satellite tidal stripping class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(satelliteTidalStrippingDu2018)                :: self
    type(inputParameters              ), intent(inout) :: parameters

    self=satelliteTidalStrippingDu2018()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function du2018ConstructorParameters

  function du2018ConstructorInternal() result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`satelliteTidalStrippingDu2018` satellite tidal stripping class.
    !!}
    implicit none
    type(satelliteTidalStrippingDu2018) :: self
    !![
    <addMetaProperty component="darkMatterProfile" name="solitonMassCore"       id="self%massCoreID"       isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonDensityCore"    id="self%densityCoreID"    isEvolvable="no"  isCreator="no"/>
    !!]
    return
  end function du2018ConstructorInternal

  double precision function du2018MassLossRate(self,node)
    !!{RST
    Set the rate of tidal mass loss from the soliton following the model of :cite:t:`du_tidal_2018`.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDarkMatterProfile, nodeComponentSatellite, treeNode
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Constants_Astronomical, only : gigaYear                      , megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Vectors                         , only : Vector_Magnitude              , Vector_Product
    implicit none
    class           (satelliteTidalStrippingDu2018 ), intent   (inout):: self
    type            (treeNode                      ), intent   (inout):: node
    class           (nodeComponentSatellite        ), pointer         :: satellite
    class           (nodeComponentDarkMatterProfile), pointer         :: darkMatterProfile
    class           (massDistributionClass         ), pointer         :: massDistribution_
    double precision                                , dimension(3    ):: position                          , velocity
    double precision                                , parameter       :: frequencyFractionalTiny=1.00000d-6
    ! Parameters appearing in the fitting function for soliton energy from equation (7) of Du et al. (2018; PRD; 97; 3507;
    ! https://ui.adsabs.harvard.edu/abs/2018PhRvD..97f3507D).
    double precision                                  , parameter        :: energyFitA      =5.89794d-5, energyFitB      =-8.72733d-2, &
         &                                                                  energyFitC      =1.67740d+0, energyFitGamma  =+1.50000d+0
    double precision                                                     :: massSatellite              , frequencyAngular            , &
         &                                                                  periodOrbital              , radius                      , &
         &                                                                  frequencyOrbital           , frequencyRadial             , &
         &                                                                  massHost                   , densityHost                 , &
         &                                                                  densityCore                , densityRatio                , &
         &                                                                  energyIm                   , massCore

    ! Get required quantities from the satellite node.
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
    darkMatterProfile  => node       %darkMatterProfile()
    massDistribution_  => node%parent%massDistribution ()
    massHost           =  max(0.0d0,massDistribution_%massEnclosedBySphere     (     radius       ))
    massCore           =            darkMatterProfile%floatRank0MetaPropertyGet(self%   massCoreID)
    densityCore        =            darkMatterProfile%floatRank0MetaPropertyGet(self%densityCoreID)
    densityHost        =  +3.0d0       &
         &                /4.0d0       &
         &                /Pi          &
         &                *massHost    &
         &                /radius  **3
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Compute the density ratio between the central density of the soliton and the average density of the host within the orbital
    ! radius. Note that equation (7) of Du et al. (2018; PRD; 97; 3507; https://ui.adsabs.harvard.edu/abs/2018PhRvD..97f3507D) is
    ! calibrated for μ < 300. When μ > 300, the core stripping through quantum tunneling is negligible, so Im(E) can either be
    ! fixed at its value at μ = 300 or simply set it to 0. Otherwise, Im(E) from the fitting formula will have an unphysical
    ! increase at very large μ. Here we choose to fix Im(E) at its value at μ = 300.
    densityRatio=min(              &
         &           +densityCore  &
         &           /densityHost, &
         &           +300.0d0      &
         &          )
    ! Compute the imaginary part of the energy eigenvalue E, using the fitting formula, equation (7) of Du et al. (2018; PRD; 97;
    ! 3507; https://ui.adsabs.harvard.edu/abs/2018PhRvD..97f3507D).
    energyIm=-exp(+energyFitA*(3.0d0*densityRatio/2.0d0/energyFitGamma)**2 &
         &        +energyFitB*(3.0d0*densityRatio/2.0d0/energyFitGamma)    &
         &        +energyFitC                                              &
         &       )                                                         &
         &   /periodOrbital
    ! Set the soliton core mass evolution rate following Du et al. (2018). From Eq. (17), the core mass obeys (1/M_c) dM_c/dt =
    ! (1/2) Im(E), where Im(E) is given by the fitting formula in Eq. (7). Here we set dM_c/dt = (1/2) Im(E) M_c for time
    ! integration.

    ! In the soliton-only phase, assume M_bound ≈ 4*M_core. Therefore, evolve the bound mass using dM_bound/dt = 4*dM_c/dt.

    du2018MassLossRate=+4.0d0    &
            &          *0.5d0    &
            &          *energyIm &
            &          *massCore
    return
  end function du2018MassLossRate
