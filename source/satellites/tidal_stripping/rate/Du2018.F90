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
  Implements a satellite tidal stripping class which computes the rate of mass loss from an :term:`FDM` solitonic core, following the model of :cite:t:`du_tidal_2018`.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <satelliteTidalStripping name="satelliteTidalStrippingDu2018" docformat="rst">
   <description>
   A satellite tidal stripping class which computes the rate at which a satellite loses mass from its :term:`FDM` solitonic core,
   following :cite:t:`du_tidal_2018`. Tides excite the core out of its ground state, allowing it to lose mass by tunneling through
   the tidal barrier; the rate is set by the imaginary part of the energy eigenvalue of the perturbed core, for which equation (7)
   of :cite:t:`du_tidal_2018` provides a fitting function in terms of the ratio of the central density of the core to the mean
   density of the host within the orbital radius.

   Note that this class assumes that the satellite has been stripped down to its solitonic core, such that its bound mass is the
   total mass of the soliton. That total is taken to be four times the mass enclosed within the core radius, and it is therefore
   four times the rate of loss of the core mass which is returned here, as this class must report the rate of loss of the *bound*
   mass. The class is consequently not meaningful when applied to a satellite which retains an envelope around its core, and is
   intended for use only through :galacticus-class:`nodeOperatorSatelliteConditionalMassLoss`, which applies it only in that
   state.
   </description>
  </satelliteTidalStripping>
  !!]
  type, extends(satelliteTidalStrippingClass) :: satelliteTidalStrippingDu2018
     !!{RST
     A satellite tidal stripping class which computes the rate of mass loss from an :term:`FDM` solitonic core, following the model of :cite:t:`du_tidal_2018`.
     !!}
     private
     class  (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     integer                                    :: massCoreID                    , densityCoreID
   contains
     final     ::                 du2018Destructor
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
    type (satelliteTidalStrippingDu2018)                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass     ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=satelliteTidalStrippingDu2018(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function du2018ConstructorParameters

  function du2018ConstructorInternal(darkMatterHaloScale_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`satelliteTidalStrippingDu2018` satellite tidal stripping class.
    !!}
    implicit none
    type (satelliteTidalStrippingDu2018)                        :: self
    class(darkMatterHaloScaleClass     ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCore"       id="self%massCoreID"       isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonDensityCore"    id="self%densityCoreID"    isEvolvable="no"  isCreator="no"/>
    !!]
    return
  end function du2018ConstructorInternal

  subroutine du2018Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`satelliteTidalStrippingDu2018` satellite tidal stripping class.
    !!}
    implicit none
    type(satelliteTidalStrippingDu2018), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine du2018Destructor

  double precision function du2018MassLossRate(self,node)
    !!{RST
    Return the rate of tidal mass loss from a satellite which has been stripped to its :term:`FDM` solitonic core, following the model of :cite:t:`du_tidal_2018`.
    !!}
    use :: Galacticus_Nodes                    , only : nodeComponentDarkMatterProfile, nodeComponentSatellite, treeNode
    use :: Mass_Distributions                  , only : massDistributionClass
    use :: Numerical_Constants_Astronomical    , only : gigaYear                      , megaParsec
    use :: Mass_Distribution_Soliton_Schive2014, only : massTotalPerMassCore
    use :: Numerical_Constants_Math            , only : Pi
    use :: Numerical_Constants_Prefixes        , only : kilo
    use :: Vectors                             , only : Vector_Magnitude              , Vector_Product
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
    double precision                                                     :: frequencyAngular           , periodOrbital               , &
         &                                                                  radius                     , &
         &                                                                  frequencyOrbital           , frequencyRadial             , &
         &                                                                  massHost                   , densityHost                 , &
         &                                                                  densityCore                , densityRatio                , &
         &                                                                  energyIm                   , massCore                    , &
         &                                                                  timescaleDynamical

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
    frequencyOrbital   =max(                  &
         &                  frequencyAngular, &
         &                  frequencyRadial   &
         &                 )
    ! Find the orbital period. Guard against a degenerate orbit, for which both the angular and radial frequencies vanish, by
    ! falling back on the dynamical timescale of the host halo.
    timescaleDynamical =   self%darkMatterHaloScale_%timescaleDynamical(node%parent)
    if (frequencyOrbital > frequencyFractionalTiny/timescaleDynamical) then
       periodOrbital   =   +2.0d0             &
            &              *Pi                &
            &              /frequencyOrbital
    else
       periodOrbital   =   +timescaleDynamical
    end if
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
    ! Compute the rate of loss of mass from the core. From equation (17) of Du et al. (2018; PRD; 97; 3507;
    ! https://ui.adsabs.harvard.edu/abs/2018PhRvD..97f3507D) the core mass obeys (1/M_c) dM_c/dt = (1/2) Im(E).
    !
    ! This class must report the rate of loss of the *bound* mass, so scale to that. The satellite is assumed to have been
    ! stripped to its core, such that its bound mass is the total mass of the soliton, and so the two rates differ by the ratio
    ! of the total soliton mass to the core mass.
    du2018MassLossRate=+massTotalPerMassCore &
            &          *0.5d0                &
            &          *energyIm             &
            &          *massCore
    return
  end function du2018MassLossRate
