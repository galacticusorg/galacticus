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

  !!{RST
  Implements the functional forms shared by the :cite:t:`diemer_splashback_2017` and :cite:t:`diemer_splashback_2020`
  splashback radius models. The two models use identical functional forms, and differ only in the values of their fitted
  coefficients---so this abstract class provides all of the machinery, and the concrete classes derived from it supply only
  the coefficient sets.
  !!}

  !+ Implemented by Andrew Benson with assistance from Claude.

  use :: Cosmology_Functions                        , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters                       , only : cosmologyParametersClass
  use :: Cosmological_Density_Field                 , only : cosmologicalMassVarianceClass          , criticalOverdensityClass
  use :: Dark_Matter_Halo_Mass_Accretion_Histories  , only : darkMatterHaloMassAccretionHistoryClass
  use :: Virial_Density_Contrast                    , only : virialDensityContrastClass
  use :: Dark_Matter_Halo_Splashback_Radii_Reference, only : splashbackReferenceScales

  ! Enumeration of the statistic of the distribution of particle apocenters which is to be predicted. The models of
  ! Diemer et al. (2017) and Diemer (2020) are calibrated separately to the mean of that distribution, and to percentiles of
  ! it.
  !![
  <enumeration docformat="rst">
   <name>splashbackDefinition</name>
   <description>
   Specifies which statistic of the distribution of particle apocenters defines the splashback radius.
   </description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="mean"         />
   <entry label="percentile50" />
   <entry label="percentile75" />
   <entry label="percentile87" />
   <entry label="percentile90" />
  </enumeration>
  !!]

  type, public :: diemerParameters
     !!{RST
     Coefficients of the fitting function for :math:`R_\mathrm{sp}/R_\mathrm{200m}` or
     :math:`M_\mathrm{sp}/M_\mathrm{200m}` in the :cite:t:`diemer_splashback_2017` and
     :cite:t:`diemer_splashback_2020` models.
     !!}
     double precision :: a0      , b0       , bOmega , bNu     , &
          &              c0      , cOmega   , cNu    , cOmega2 , &
          &              cNu2    , a0P      , b0P    , bOmegaP , &
          &              bOmegaP2, bNuP     , cOmegaP, cOmegaP2, &
          &              cOmega2P, cOmega2P2, cNuP   , cNu2P
  end type diemerParameters

  type, public :: diemerScatterParameters
     !!{RST
     Coefficients of the fitting function for the 68% scatter in :math:`R_\mathrm{sp}/R_\mathrm{200m}` or
     :math:`M_\mathrm{sp}/M_\mathrm{200m}` in the :cite:t:`diemer_splashback_2017` and
     :cite:t:`diemer_splashback_2020` models.
     !!}
     double precision :: sigma0, sigmaGamma, sigmaNu, sigmaP
  end type diemerScatterParameters

  ! Minimum value permitted for the coefficient "C" in the Diemer fitting functions. Since "C" appears in the denominator of
  ! an exponential it must be bounded away from zero.
  double precision, parameter, public :: diemerCoefficientCMinimum=1.0d-4
  ! Minimum value permitted for the scatter predicted by the Diemer fitting functions.
  double precision, parameter, public :: diemerScatterMinimum     =2.0d-2

  ! Fitting-function coefficients for the Diemer2017 model.
  type(diemerParameters), parameter, public :: diemer2017ParametersRadiusMean      =diemerParameters(                        &
     &                                                                                               a0       = +0.649783d0, &
     &                                                                                               b0       = +0.600362d0, &
     &                                                                                               bOmega   = +0.091996d0, &
     &                                                                                               bNu      = +0.061557d0, &
     &                                                                                               c0       = -0.806288d0, &
     &                                                                                               cOmega   =+17.520522d0, &
     &                                                                                               cNu      = -0.293465d0, &
     &                                                                                               cOmega2  = -9.624342d0, &
     &                                                                                               cNu2     = +0.039196d0, &
     &                                                                                               a0P      = +0.000000d0, &
     &                                                                                               b0P      = +0.000000d0, &
     &                                                                                               bOmegaP  = +0.000000d0, &
     &                                                                                               bOmegaP2 = +0.000000d0, &
     &                                                                                               bNuP     = +0.000000d0, &
     &                                                                                               cOmegaP  = +0.000000d0, &
     &                                                                                               cOmegaP2 = +0.000000d0, &
     &                                                                                               cOmega2P = +0.000000d0, &
     &                                                                                               cOmega2P2= +0.000000d0, &
     &                                                                                               cNuP     = +0.000000d0, &
     &                                                                                               cNu2P    = +0.000000d0  &
     &                                                                                              )
  type(diemerParameters), parameter, public :: diemer2017ParametersRadiusPercentile=diemerParameters(                        &
     &                                                                                               a0       = +0.320332d0, &
     &                                                                                               b0       = +0.267433d0, &
     &                                                                                               bOmega   = +0.113389d0, &
     &                                                                                               bNu      = +0.207989d0, &
     &                                                                                               c0       = -0.959629d0, &
     &                                                                                               cOmega   =+16.245894d0, &
     &                                                                                               cNu      = +0.000000d0, &
     &                                                                                               cOmega2  = -9.497861d0, &
     &                                                                                               cNu2     = -0.018484d0, &
     &                                                                                               a0P      = +0.614807d0, &
     &                                                                                               b0P      = +0.545238d0, &
     &                                                                                               bOmegaP  = +0.000000d0, &
     &                                                                                               bOmegaP2 = +0.000000d0, &
     &                                                                                               bNuP     = -0.223282d0, &
     &                                                                                               cOmegaP  = +0.003941d0, &
     &                                                                                               cOmegaP2 = +8.969094d0, &
     &                                                                                               cOmega2P = -0.000485d0, &
     &                                                                                               cOmega2P2=+10.613168d0, &
     &                                                                                               cNuP     = -0.451066d0, &
     &                                                                                               cNu2P    = +0.088029d0  &
     &                                                                                              )
  type(diemerParameters), parameter, public :: diemer2017ParametersMassMean        =diemerParameters(                       &
     &                                                                                               a0       =+0.679244d0, &
     &                                                                                               b0       =+0.405083d0, &
     &                                                                                               bOmega   =+0.291925d0, &
     &                                                                                               bNu      =+0.000000d0, &
     &                                                                                               c0       =+3.365943d0, &
     &                                                                                               cOmega   =+1.469818d0, &
     &                                                                                               cNu      =-0.075635d0, &
     &                                                                                               cOmega2  =+0.000000d0, &
     &                                                                                               cNu2     =+0.000000d0, &
     &                                                                                               a0P      =+0.000000d0, &
     &                                                                                               b0P      =+0.000000d0, &
     &                                                                                               bOmegaP  =+0.000000d0, &
     &                                                                                               bOmegaP2 =+0.000000d0, &
     &                                                                                               bNuP     =+0.000000d0, &
     &                                                                                               cOmegaP  =+0.000000d0, &
     &                                                                                               cOmegaP2 =+0.000000d0, &
     &                                                                                               cOmega2P =+0.000000d0, &
     &                                                                                               cOmega2P2=+0.000000d0, &
     &                                                                                               cNuP     =+0.000000d0, &
     &                                                                                               cNu2P    =+0.000000d0  &
     &                                                                                              )
  type(diemerParameters), parameter, public :: diemer2017ParametersMassPercentile  =diemerParameters(                       &
     &                                                                                               a0       =+0.264765d0, &
     &                                                                                               b0       =+0.666040d0, &
     &                                                                                               bOmega   =+0.168814d0, &
     &                                                                                               bNu      =+0.000000d0, &
     &                                                                                               c0       =+4.728709d0, &
     &                                                                                               cOmega   =+2.388866d0, &
     &                                                                                               cNu      =-0.084108d0, &
     &                                                                                               cOmega2  =+0.000000d0, &
     &                                                                                               cNu2     =+0.000000d0, &
     &                                                                                               a0P      =+0.843509d0, &
     &                                                                                               b0P      =-0.639169d0, &
     &                                                                                               bOmegaP  =+0.003195d0, &
     &                                                                                               bOmegaP2 =+4.939266d0, &
     &                                                                                               bNuP     =+0.225399d0, &
     &                                                                                               cOmegaP  =-0.705712d0, &
     &                                                                                               cOmegaP2 =-1.241920d0, &
     &                                                                                               cOmega2P =+0.000000d0, &
     &                                                                                               cOmega2P2=+0.000000d0, &
     &                                                                                               cNuP     =-0.391103d0, &
     &                                                                                               cNu2P    =+0.074216d0  &
     &                                                                                              )
  type(diemerScatterParameters), parameter, public :: diemer2017ScatterRadiusMean         =diemerScatterParameters(                        &
     &                                                                                                             sigma0    =+0.052645d0, &
     &                                                                                                             sigmaGamma=+0.003846d0, &
     &                                                                                                             sigmaNu   =-0.012054d0, &
     &                                                                                                             sigmaP    =+0.000000d0  &
     &                                                                                                            )
  type(diemerScatterParameters), parameter, public :: diemer2017ScatterRadiusPercentile   =diemerScatterParameters(                        &
     &                                                                                                             sigma0    =+0.044548d0, &
     &                                                                                                             sigmaGamma=+0.004404d0, &
     &                                                                                                             sigmaNu   =-0.014636d0, &
     &                                                                                                             sigmaP    =+0.022637d0  &
     &                                                                                                            )
  type(diemerScatterParameters), parameter, public :: diemer2017ScatterMassMean           =diemerScatterParameters(                        &
     &                                                                                                             sigma0    =+0.052815d0, &
     &                                                                                                             sigmaGamma=+0.002456d0, &
     &                                                                                                             sigmaNu   =-0.011182d0, &
     &                                                                                                             sigmaP    =+0.000000d0  &
     &                                                                                                            )
  type(diemerScatterParameters), parameter, public :: diemer2017ScatterMassPercentile     =diemerScatterParameters(                        &
     &                                                                                                             sigma0    =+0.027594d0, &
     &                                                                                                             sigmaGamma=+0.002330d0, &
     &                                                                                                             sigmaNu   =-0.012491d0, &
     &                                                                                                             sigmaP    =+0.047344d0  &
     &                                                                                                            )

  ! Fitting-function coefficients for the Diemer2020 model.
  type(diemerParameters), parameter, public :: diemer2020ParametersRadiusMean      =diemerParameters(                        &
     &                                                                                               a0       = +0.659745d0, &
     &                                                                                               b0       = +0.556171d0, &
     &                                                                                               bOmega   = +0.114053d0, &
     &                                                                                               bNu      = +0.069775d0, &
     &                                                                                               c0       = -0.850819d0, &
     &                                                                                               cOmega   =+18.446356d0, &
     &                                                                                               cNu      = -0.333245d0, &
     &                                                                                               cOmega2  =-10.059591d0, &
     &                                                                                               cNu2     = +0.047432d0, &
     &                                                                                               a0P      = +0.000000d0, &
     &                                                                                               b0P      = +0.000000d0, &
     &                                                                                               bOmegaP  = +0.000000d0, &
     &                                                                                               bOmegaP2 = +0.000000d0, &
     &                                                                                               bNuP     = +0.000000d0, &
     &                                                                                               cOmegaP  = +0.000000d0, &
     &                                                                                               cOmegaP2 = +0.000000d0, &
     &                                                                                               cOmega2P = +0.000000d0, &
     &                                                                                               cOmega2P2= +0.000000d0, &
     &                                                                                               cNuP     = +0.000000d0, &
     &                                                                                               cNu2P    = +0.000000d0  &
     &                                                                                              )
  type(diemerParameters), parameter, public :: diemer2020ParametersRadiusPercentile=diemerParameters(                        &
     &                                                                                               a0       = +0.307107d0, &
     &                                                                                               b0       = +0.250844d0, &
     &                                                                                               bOmega   = +0.152731d0, &
     &                                                                                               bNu      = +0.195627d0, &
     &                                                                                               c0       = -1.221448d0, &
     &                                                                                               cOmega   =+17.537448d0, &
     &                                                                                               cNu      = +0.000000d0, &
     &                                                                                               cOmega2  =-10.315817d0, &
     &                                                                                               cNu2     = -0.018938d0, &
     &                                                                                               a0P      = +0.642779d0, &
     &                                                                                               b0P      = +0.507395d0, &
     &                                                                                               bOmegaP  = +0.000000d0, &
     &                                                                                               bOmegaP2 = +0.000000d0, &
     &                                                                                               bNuP     = -0.212788d0, &
     &                                                                                               cOmegaP  = +0.002358d0, &
     &                                                                                               cOmegaP2 = +9.711469d0, &
     &                                                                                               cOmega2P = -0.000550d0, &
     &                                                                                               cOmega2P2=+10.762591d0, &
     &                                                                                               cNuP     = -0.473487d0, &
     &                                                                                               cNu2P    = +0.094019d0  &
     &                                                                                              )
  type(diemerParameters), parameter, public :: diemer2020ParametersMassMean        =diemerParameters(                       &
     &                                                                                               a0       =+0.696229d0, &
     &                                                                                               b0       =+0.373627d0, &
     &                                                                                               bOmega   =+0.300490d0, &
     &                                                                                               bNu      =+0.000000d0, &
     &                                                                                               c0       =+3.344518d0, &
     &                                                                                               cOmega   =+1.371785d0, &
     &                                                                                               cNu      =-0.082529d0, &
     &                                                                                               cOmega2  =+0.000000d0, &
     &                                                                                               cNu2     =+0.000000d0, &
     &                                                                                               a0P      =+0.000000d0, &
     &                                                                                               b0P      =+0.000000d0, &
     &                                                                                               bOmegaP  =+0.000000d0, &
     &                                                                                               bOmegaP2 =+0.000000d0, &
     &                                                                                               bNuP     =+0.000000d0, &
     &                                                                                               cOmegaP  =+0.000000d0, &
     &                                                                                               cOmegaP2 =+0.000000d0, &
     &                                                                                               cOmega2P =+0.000000d0, &
     &                                                                                               cOmega2P2=+0.000000d0, &
     &                                                                                               cNuP     =+0.000000d0, &
     &                                                                                               cNu2P    =+0.000000d0  &
     &                                                                                              )
  type(diemerParameters), parameter, public :: diemer2020ParametersMassPercentile  =diemerParameters(                       &
     &                                                                                               a0       =+0.287428d0, &
     &                                                                                               b0       =+0.661551d0, &
     &                                                                                               bOmega   =+0.132142d0, &
     &                                                                                               bNu      =+0.000000d0, &
     &                                                                                               c0       =+4.591265d0, &
     &                                                                                               cOmega   =+3.092819d0, &
     &                                                                                               cNu      =-0.115458d0, &
     &                                                                                               cOmega2  =+0.000000d0, &
     &                                                                                               cNu2     =+0.000000d0, &
     &                                                                                               a0P      =+0.822820d0, &
     &                                                                                               b0P      =-0.656698d0, &
     &                                                                                               bOmegaP  =+0.003182d0, &
     &                                                                                               bOmegaP2 =+4.953590d0, &
     &                                                                                               bNuP     =+0.288168d0, &
     &                                                                                               cOmegaP  =-0.660871d0, &
     &                                                                                               cOmegaP2 =-1.105000d0, &
     &                                                                                               cOmega2P =+0.000000d0, &
     &                                                                                               cOmega2P2=+0.000000d0, &
     &                                                                                               cNuP     =-0.376065d0, &
     &                                                                                               cNu2P    =+0.078376d0  &
     &                                                                                              )
  type(diemerScatterParameters), parameter, public :: diemer2020ScatterRadiusMean         =diemerScatterParameters(                        &
     &                                                                                                             sigma0    =+0.050135d0, &
     &                                                                                                             sigmaGamma=+0.003548d0, &
     &                                                                                                             sigmaNu   =-0.010766d0, &
     &                                                                                                             sigmaP    =+0.000000d0  &
     &                                                                                                            )
  type(diemerScatterParameters), parameter, public :: diemer2020ScatterRadiusPercentile   =diemerScatterParameters(                        &
     &                                                                                                             sigma0    =+0.041872d0, &
     &                                                                                                             sigmaGamma=+0.004279d0, &
     &                                                                                                             sigmaNu   =-0.014068d0, &
     &                                                                                                             sigmaP    =+0.023470d0  &
     &                                                                                                            )
  type(diemerScatterParameters), parameter, public :: diemer2020ScatterMassMean           =diemerScatterParameters(                        &
     &                                                                                                             sigma0    =+0.045588d0, &
     &                                                                                                             sigmaGamma=+0.001679d0, &
     &                                                                                                             sigmaNu   =-0.007945d0, &
     &                                                                                                             sigmaP    =+0.000000d0  &
     &                                                                                                            )
  type(diemerScatterParameters), parameter, public :: diemer2020ScatterMassPercentile     =diemerScatterParameters(                        &
     &                                                                                                             sigma0    =+0.022368d0, &
     &                                                                                                             sigmaGamma=+0.000916d0, &
     &                                                                                                             sigmaNu   =-0.009101d0, &
     &                                                                                                             sigmaP    =+0.045386d0  &
     &                                                                                                            )

  ! The fitting functions are made public so that they can be validated directly against reference values by the
  ! "tests.splashback_radii.exe" unit test.
  public :: diemerRatio, diemerScatter, diemerPercentileValue

  !![
  <darkMatterHaloSplashbackRadius name="darkMatterHaloSplashbackRadiusDiemer" abstract="yes" docformat="rst">
   <description>
   An abstract splashback radius class providing the functional forms shared by the :cite:t:`diemer_splashback_2017` and
   :cite:t:`diemer_splashback_2020` models.
   </description>
  </darkMatterHaloSplashbackRadius>
  !!]
  type, abstract, extends(darkMatterHaloSplashbackRadiusClass) :: darkMatterHaloSplashbackRadiusDiemer
     !!{RST
     An abstract splashback radius class implementing the :cite:t:`diemer_splashback_2017` and
     :cite:t:`diemer_splashback_2020` functional forms.
     !!}
     private
     class(cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_                 => null()
     class(cosmologyParametersClass               ), pointer :: cosmologyParameters_                => null()
     class(criticalOverdensityClass               ), pointer :: criticalOverdensity_                => null()
     class(cosmologicalMassVarianceClass          ), pointer :: cosmologicalMassVariance_           => null()
     class(virialDensityContrastClass             ), pointer :: virialDensityContrast_              => null()
     class(darkMatterHaloMassAccretionHistoryClass), pointer :: darkMatterHaloMassAccretionHistory_ => null()
     type (enumerationSplashbackDefinitionType    )          :: definition
     type (diemerParameters                       )          :: parametersRadius                            , parametersMass
     type (diemerScatterParameters                )          :: parametersRadiusScatter                     , parametersMassScatter
   contains
     !![
     <methods docformat="rst">
       <method method="propertiesGet" description="Return the properties of the halo needed to evaluate the splashback fitting functions."/>
     </methods>
     !!]
     procedure :: radius        => diemerRadius
     procedure :: radiusRatio   => diemerRadiusRatio
     procedure :: mass          => diemerMass
     procedure :: massRatio     => diemerMassRatio
     procedure :: radiusScatter => diemerRadiusScatter
     procedure :: massScatter   => diemerMassScatter
     procedure :: propertiesGet => diemerPropertiesGet
  end type darkMatterHaloSplashbackRadiusDiemer

contains

  subroutine diemerPropertiesGet(self,node,gamma,peakHeight,omegaMatter,radius200Mean,mass200Mean,massDistribution_)
    !!{RST
    Return the properties of the halo in ``node`` needed to evaluate the :cite:t:`diemer_splashback_2017` and
    :cite:t:`diemer_splashback_2020` splashback fitting functions. Specifically, the mass accretion rate,
    :math:`\Gamma = \mathrm{d}\log M / \mathrm{d}\log a`, the peak height, :math:`\nu`, computed from
    :math:`M_\mathrm{200m}`, the matter density parameter, :math:`\Omega_\mathrm{M}`, and :math:`R_\mathrm{200m}` and
    :math:`M_\mathrm{200m}` themselves.

    Note that :math:`\Gamma` is evaluated here as the *instantaneous* logarithmic growth rate of the halo mass, whereas
    these models were calibrated against :math:`\Gamma_\mathrm{dyn}`---the growth rate averaged over one dynamical time. For
    the :galacticus-class:`darkMatterHaloMassAccretionHistoryDiemer2020` mass accretion history class the two are identical
    by construction, as that class directly parameterizes :math:`\Gamma_\mathrm{dyn}(\nu,z)`.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (darkMatterHaloSplashbackRadiusDiemer), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    double precision                                      , intent(  out)           :: gamma            , peakHeight   , &
         &                                                                             omegaMatter      , radius200Mean, &
         &                                                                             mass200Mean
    class           (massDistributionClass               ), intent(inout), optional :: massDistribution_
    class           (nodeComponentBasic                  ), pointer                 :: basic
    double precision                                                                :: time             , rootVariance

    basic        =>  node                    %basic             (    )
    time         =   basic                   %time              (    )
    omegaMatter  =   self%cosmologyFunctions_%omegaMatterEpochal(time)
    ! Find the mass and radius enclosing a mean density contrast of 200 relative to the mean density of the universe.
    call splashbackReferenceScales(                                                    &
         &                                                node                       , &
         &                         cosmologyParameters_  =self%cosmologyParameters_  , &
         &                         cosmologyFunctions_   =self%cosmologyFunctions_   , &
         &                         virialDensityContrast_=self%virialDensityContrast_, &
         &                         mass                  =     mass200Mean           , &
         &                         radius                =     radius200Mean         , &
         &                         massDistribution_     =     massDistribution_       &
         &                        )
    ! Compute the peak height corresponding to that mass.
    rootVariance = self%cosmologicalMassVariance_%rootVariance(time=time,mass=mass200Mean          )
    if (rootVariance > 0.0d0) then
       peakHeight=+self%criticalOverdensity_     %value       (time=time,mass=mass200Mean,node=node) &
            &     /                                            rootVariance
    else
       peakHeight=+0.0d0
    end if
    ! Compute the logarithmic growth rate of the halo mass, Γ = dlogM/dloga = (dM/dt)/(M H).
    if (basic%mass() > 0.0d0) then
       gamma  =+self%darkMatterHaloMassAccretionHistory_%massAccretionRate(node,time)                                      &
            &  /basic                                   %mass             (         )                                      &
            &  /self%cosmologyFunctions_                %expansionRate    (self%cosmologyFunctions_%expansionFactor(time))
    else
       gamma  =+0.0d0
    end if
    return
  end subroutine diemerPropertiesGet

  double precision function diemerRadiusRatio(self,node,massDistribution_) result(radiusRatio)
    !!{RST
    Return the ratio :math:`R_\mathrm{sp}/R_\mathrm{200m}` for the halo in ``node``.
    !!}
    implicit none
    class           (darkMatterHaloSplashbackRadiusDiemer), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    class           (massDistributionClass               ), intent(inout), optional :: massDistribution_
    double precision                                                                :: gamma            , peakHeight   , &
         &                                                                             omegaMatter      , radius200Mean, &
         &                                                                             mass200Mean

    call self%propertiesGet(node,gamma,peakHeight,omegaMatter,radius200Mean,mass200Mean,massDistribution_)
    radiusRatio=diemerRatio(self%parametersRadius,diemerPercentileValue(self%definition),gamma,peakHeight,omegaMatter)
    return
  end function diemerRadiusRatio

  double precision function diemerMassRatio(self,node,massDistribution_) result(massRatio)
    !!{RST
    Return the ratio :math:`M_\mathrm{sp}/M_\mathrm{200m}` for the halo in ``node``.
    !!}
    implicit none
    class           (darkMatterHaloSplashbackRadiusDiemer), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    class           (massDistributionClass               ), intent(inout), optional :: massDistribution_
    double precision                                                                :: gamma            , peakHeight   , &
         &                                                                             omegaMatter      , radius200Mean, &
         &                                                                             mass200Mean

    call self%propertiesGet(node,gamma,peakHeight,omegaMatter,radius200Mean,mass200Mean,massDistribution_)
    massRatio=diemerRatio(self%parametersMass,diemerPercentileValue(self%definition),gamma,peakHeight,omegaMatter)
    return
  end function diemerMassRatio

  double precision function diemerRadius(self,node,massDistribution_) result(radius)
    !!{RST
    Return the splashback radius (in Mpc) of the halo in ``node``.
    !!}
    implicit none
    class           (darkMatterHaloSplashbackRadiusDiemer), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    class           (massDistributionClass               ), intent(inout), optional :: massDistribution_
    double precision                                                                :: gamma            , peakHeight   , &
         &                                                                             omegaMatter      , radius200Mean, &
         &                                                                             mass200Mean

    call self%propertiesGet(node,gamma,peakHeight,omegaMatter,radius200Mean,mass200Mean,massDistribution_)
    radius=+diemerRatio(self%parametersRadius,diemerPercentileValue(self%definition),gamma,peakHeight,omegaMatter) &
         & *radius200Mean
    return
  end function diemerRadius

  double precision function diemerMass(self,node,massDistribution_) result(mass)
    !!{RST
    Return the splashback mass (in :math:`\mathrm{M}_\odot`) of the halo in ``node``.
    !!}
    implicit none
    class           (darkMatterHaloSplashbackRadiusDiemer), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    class           (massDistributionClass               ), intent(inout), optional :: massDistribution_
    double precision                                                                :: gamma            , peakHeight   , &
         &                                                                             omegaMatter      , radius200Mean, &
         &                                                                             mass200Mean

    call self%propertiesGet(node,gamma,peakHeight,omegaMatter,radius200Mean,mass200Mean,massDistribution_)
    mass=+diemerRatio(self%parametersMass,diemerPercentileValue(self%definition),gamma,peakHeight,omegaMatter) &
         & *mass200Mean
    return
  end function diemerMass

  double precision function diemerRadiusScatter(self,node,massDistribution_) result(scatter)
    !!{RST
    Return the 68% scatter in :math:`\log_{10}(R_\mathrm{sp}/R_\mathrm{200m})` for the halo in ``node``.
    !!}
    implicit none
    class           (darkMatterHaloSplashbackRadiusDiemer), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    class           (massDistributionClass               ), intent(inout), optional :: massDistribution_
    double precision                                                                :: gamma            , peakHeight   , &
         &                                                                             omegaMatter      , radius200Mean, &
         &                                                                             mass200Mean

    call self%propertiesGet(node,gamma,peakHeight,omegaMatter,radius200Mean,mass200Mean,massDistribution_)
    scatter=diemerScatter(self%parametersRadiusScatter,diemerPercentileValue(self%definition),gamma,peakHeight)
    return
  end function diemerRadiusScatter

  double precision function diemerMassScatter(self,node,massDistribution_) result(scatter)
    !!{RST
    Return the 68% scatter in :math:`\log_{10}(M_\mathrm{sp}/M_\mathrm{200m})` for the halo in ``node``.
    !!}
    implicit none
    class           (darkMatterHaloSplashbackRadiusDiemer), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    class           (massDistributionClass               ), intent(inout), optional :: massDistribution_
    double precision                                                                :: gamma            , peakHeight   , &
         &                                                                             omegaMatter      , radius200Mean, &
         &                                                                             mass200Mean

    call self%propertiesGet(node,gamma,peakHeight,omegaMatter,radius200Mean,mass200Mean,massDistribution_)
    scatter=diemerScatter(self%parametersMassScatter,diemerPercentileValue(self%definition),gamma,peakHeight)
    return
  end function diemerMassScatter

  double precision function diemerPercentileValue(definition) result(p)
    !!{RST
    Return the numerical value, :math:`p`, of the percentile corresponding to the given splashback ``definition``. A value of
    :math:`p=-1` indicates the mean of the distribution of particle apocenters, following :cite:t:`diemer_splashback_2020`.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type(enumerationSplashbackDefinitionType), intent(in   ) :: definition

    select case (definition%ID)
    case (splashbackDefinitionMean        %ID)
       p=-1.00d0
    case (splashbackDefinitionPercentile50%ID)
       p=+0.50d0
    case (splashbackDefinitionPercentile75%ID)
       p=+0.75d0
    case (splashbackDefinitionPercentile87%ID)
       p=+0.87d0
    case (splashbackDefinitionPercentile90%ID)
       p=+0.90d0
    case default
       p=+0.00d0
       call Error_Report('unknown splashback definition'//{introspection:location})
    end select
    return
  end function diemerPercentileValue

  double precision function diemerRatio(parameters_,p,gamma,peakHeight,omegaMatter) result(ratio)
    !!{RST
    Evaluate the fitting function

    .. math::

       {X_\mathrm{sp} \over X_\mathrm{200m}} = A + B \exp\left( - {\Gamma \over C} \right),

    where :math:`X` is either radius or mass, and

    .. math::

       A &= A_0 + p A_{0,\mathrm{p}}, \\
       B &= \left( B_0 + B_\Omega \Omega_\mathrm{M} \right) \left( 1 + B_\nu \nu \right), \\
       C &= \left( C_0 + C_\Omega \Omega_\mathrm{M} + C_{\Omega,2} \Omega_\mathrm{M}^2 \right) \left( 1 + C_\nu \nu + C_{\nu,2} \nu^2 \right),

    following :cite:t:`diemer_splashback_2017` and :cite:t:`diemer_splashback_2020`. Here :math:`\Gamma` is the mass
    accretion rate, :math:`\nu` is the peak height computed from :math:`M_\mathrm{200m}`, :math:`\Omega_\mathrm{M}` is
    the matter density parameter at the epoch of interest, and :math:`p` is the percentile of the distribution of particle
    apocenters being predicted (with :math:`p=-1` indicating the mean of that distribution).
    !!}
    implicit none
    type            (diemerParameters), intent(in   ) :: parameters_
    double precision                  , intent(in   ) :: p          , gamma      , &
         &                                               peakHeight , omegaMatter
    double precision                                  :: a          , b          , &
         &                                               c          , b0         , &
         &                                               bOmega     , bNu        , &
         &                                               cOmega     , cOmega2    , &
         &                                               cNu        , cNu2

    ! Evaluate the percentile-dependent coefficients.
    a      =+parameters_%a0      +p*parameters_%a0P
    b0     =+parameters_%b0      +p*parameters_%b0P
    bOmega =+parameters_%bOmega  +  parameters_%bOmegaP *exp(p*parameters_%bOmegaP2 )
    bNu    =+parameters_%bNu     +p*parameters_%bNuP
    cOmega =+parameters_%cOmega  +  parameters_%cOmegaP *exp(p*parameters_%cOmegaP2 )
    cOmega2=+parameters_%cOmega2 +  parameters_%cOmega2P*exp(p*parameters_%cOmega2P2)
    cNu    =+parameters_%cNu     +p*parameters_%cNuP
    cNu2   =+parameters_%cNu2    +p*parameters_%cNu2P
    ! Evaluate the coefficients of the fitting function.
    b      =+(                          &
         &    +b0                       &
         &    +bOmega *omegaMatter      &
         &   )                          &
         &  *(                          &
         &    +1.0d0                    &
         &    +bNu    *peakHeight       &
         &   )
    c      =+(                          &
         &    +parameters_%c0           &
         &    +cOmega *omegaMatter      &
         &    +cOmega2*omegaMatter**2   &
         &   )                          &
         &  *(                          &
         &    +1.0d0                    &
         &    +cNu    *peakHeight       &
         &    +cNu2   *peakHeight   **2 &
         &   )
    ! The coefficient "C" appears in the denominator of an exponential, so must be bounded away from zero.
    c      =max(c,diemerCoefficientCMinimum)
    ratio  =+a                          &
         &  +b                          &
         &  *exp(-gamma/c)
    return
  end function diemerRatio

  double precision function diemerScatter(parameters_,p,gamma,peakHeight) result(scatter)
    !!{RST
    Evaluate the 68% scatter,

    .. math::

       \sigma = \sigma_0 + \sigma_\Gamma \Gamma + \sigma_\nu \nu + \sigma_\mathrm{p} p,

    in :math:`\log_{10}(X_\mathrm{sp}/X_\mathrm{200m})` following :cite:t:`diemer_splashback_2017` and
    :cite:t:`diemer_splashback_2020`. The result is bounded below at :math:`\sigma=0.02`.
    !!}
    implicit none
    type            (diemerScatterParameters), intent(in   ) :: parameters_
    double precision                         , intent(in   ) :: p          , gamma, &
         &                                                      peakHeight

    scatter=+parameters_%sigma0                &
         &  +parameters_%sigmaGamma*gamma      &
         &  +parameters_%sigmaNu   *peakHeight &
         &  +parameters_%sigmaP    *p
    scatter=max(scatter,diemerScatterMinimum)
    return
  end function diemerScatter
