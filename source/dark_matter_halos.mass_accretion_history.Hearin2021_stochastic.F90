!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

  !% An implementation of dark matter halo mass accretion histories using Andrew Hearin's rolling power-law model with stochastic sampling of parameters.

  !# <enumeration>
  !#  <name>mass</name>
  !#  <description>Enumeration of mass limits.</description>
  !#  <decodeFunction>yes</decodeFunction>
  !#  <indexing>1</indexing>
  !#  <entry label="low" />
  !#  <entry label="high"/>
  !# </enumeration>

  !# <enumeration>
  !#  <name>formation</name>
  !#  <description>Enumeration of formation epochs.</description>
  !#  <decodeFunction>yes</decodeFunction>
  !#  <indexing>1</indexing>
  !#  <entry label="early"/>
  !#  <entry label="late" />
  !# </enumeration>

  !# <enumeration>
  !#  <name>parameter</name>
  !#  <description>Enumeration of parameters.</description>
  !#  <decodeFunction>yes</decodeFunction>
  !#  <indexing>1</indexing>
  !#  <entry label="log10PowerLawIndexEarly"/>
  !#  <entry label="log10PowerLawIndexLate" />
  !#  <entry label="log10TimeZero"          />
  !# </enumeration>

  !# <darkMatterHaloMassAccretionHistory name="darkMatterHaloMassAccretionHistoryHearin2021Stochastic">
  !#  <description>Dark matter halo mass accretion histories using Andrew Hearin's rolling power-law model.</description>
  !# </darkMatterHaloMassAccretionHistory>
  type, extends(darkMatterHaloMassAccretionHistoryHearin2021) :: darkMatterHaloMassAccretionHistoryHearin2021Stochastic
     !% A dark matter halo mass accretion history class using Andrew Hearin's rolling power-law model with stochastic sampling of parameters.
     private
     double precision                     :: fractionLateLow     , fractionLateHigh, &
          &                                  timeZeroLogarithmic_
     double precision, dimension(2,2,3  ) :: means
     double precision, dimension(2,2,3,3) :: cholesky
   contains
     !# <methods>
     !#  <method description="Sample parameters for the given node from the distribution function." method="sample"      />
     !#  <method description="Return the fraction of late-forming halos."                           method="fractionLate"/>
     !# </methods>
     procedure :: powerLawIndexEarly_ => hearin2021StochasticPowerLawIndexEarly
     procedure :: powerLawIndexLate_  => hearin2021StochasticPowerLawIndexLate
     procedure :: timeMaximum_        => hearin2021StochasticTimeMaximum
     procedure :: timeZeroLogarithmic => hearin2021StochasticLog10TimeZero
     procedure :: sample              => hearin2021StochasticSample
     procedure :: fractionLate        => hearin2021StochasticFractionLate
  end type darkMatterHaloMassAccretionHistoryHearin2021Stochastic

  interface darkMatterHaloMassAccretionHistoryHearin2021Stochastic
     !% Constructors for the {\normalfont \ttfamily hearin2021Stochastic} dark matter halo mass accretion history class.
     module procedure hearin2021StochasticConstructorParameters
     module procedure hearin2021StochasticConstructorInternal
  end interface darkMatterHaloMassAccretionHistoryHearin2021Stochastic

  ! Set default values for the mean parameter functions.
  double precision, dimension(2,2,3  ), parameter :: meanDefault                =reshape(                                         &
       &                                                                                 [                                        &
       !                                                                                 Each row is the mean of the stated parameter for:
       !                                                                                   early foming/low mass : early forming/high mass : late forming/low mass : late forming/high mass
       !                                                                                  log₁₀(α_early)
       &                                                                                  +0.460d0, +0.680d0, +0.030d0, +0.600d0, &
       !                                                                                  log₁₀(α_late )
       &                                                                                  -0.800d0, +0.760d0, -1.400d0, +0.790d0, &
       !                                                                                  log₁₀(t₀     )
       &                                                                                  -0.250d0, -0.240d0, +0.380d0, +0.750d0  &
       &                                                                                 ]                                      , &
       &                                                                                 [2,2,3]                                  &
       &                                                                                )
  
  ! Set default values for Cholesky matrix elements of parameter functions.
  double precision, dimension(2,2,3,3), parameter :: choleskyDefault            =reshape(                                         &
       &                                                                                 [                                        &
       !                                                                                  Each row is the mean of the stated parameter for:
       !                                                                                   early foming/low mass : early forming/high mass : late forming/low mass : late forming/high mass
       !                                                                                  log₁₀(α_early)-log₁₀(α_early)
       &                                                                                  -0.630d0, -0.750d0, -0.830d0, -1.000d0, &
       !                                                                                  log₁₀(α_early)-log₁₀(α_late )
       &                                                                                  +0.010d0, -0.120d0, -0.330d0, -0.130d0, &
       !                                                                                  log₁₀(α_early)-log₁₀(t₀     )
       &                                                                                  -0.070d0, -0.050d0, -0.350d0, -0.050d0, &
       !                                                                                  log₁₀(α_late )-log₁₀(α_early)
       &                                                                                  +0.010d0, -0.120d0, -0.330d0, -0.130d0, &
       !                                                                                  log₁₀(α_late )-log₁₀(α_late )
       &                                                                                  -0.050d0, -1.200d0, -0.290d0, -0.620d0, &
       !                                                                                  log₁₀(α_late )-log₁₀(t₀     )
       &                                                                                  -0.300d0, -0.050d0, +0.150d0, -0.070d0, &
       !                                                                                  log₁₀(t₀     )-log₁₀(α_early)
       &                                                                                  -0.070d0, -0.050d0, -0.350d0, -0.050d0, &
       !                                                                                  log₁₀(t₀     )-log₁₀(α_late )
       &                                                                                  -0.300d0, -0.050d0, +0.150d0, -0.070d0, &
       !                                                                                  log₁₀(t₀     )-log₁₀(t₀     )
       &                                                                                  -0.810d0, -0.990d0, -0.630d0, -0.700d0  &
       &                                                                                 ]                                      , &
       &                                                                                 [2,2,3,3]                                &
       &                                                                                )

  double precision, dimension(2,3), parameter :: log10MassZeroMean              =reshape(                 &
       &                                                                                 [                &
       !                                                                                  Each row is the zero-point of the mass transition for the mean for:
       !                                                                                   early foming : late forming
       !                                                                                  log₁₀(α_early)
       &                                                                                  13.0d0, 13.0d0, &
       !                                                                                  log10(alpha_late)
       &                                                                                  13.0d0, 14.0d0, &
       !                                                                                  log10(t_0)
       &                                                                                  13.0d0, 13.0d0  &
       &                                                                                 ]              , &
       &                                                                                 [2,3]            &
       &                                                                                )
  
  double precision, dimension(2,3), parameter :: rollOverRateLogMassMean        =reshape(               &
       &                                                                                 [              &
       !                                                                                  Each row is the roll-over rate of the mass transition for the mean for:
       !                                                                                   early foming : late forming
       !                                                                                  log₁₀(α_early)
       &                                                                                  0.5d0, 0.5d0, &
       !                                                                                  log10(alpha_late)
       &                                                                                  0.5d0, 0.5d0, &
       !                                                                                  log10(t_0)
       &                                                                                  0.5d0, 0.5d0  &
       &                                                                                 ]            , &
       &                                                                                 [2,3]          &
       &                                                                                 )

  double precision, dimension(2,3,3), parameter :: log10MassZeroCholesky        =reshape(                   &
       &                                                                                 [                  &
       !                                                                                  Each row is the zero-point of the mass transition for the Cholesky matrix element for:
       !                                                                                   early foming : late forming
       !                                                                                  log₁₀(α_early)-log₁₀(α_early)
       &                                                                                  12.00d0, 13.00d0, &
       !                                                                                  log₁₀(α_early)-log₁₀(α_late )
       &                                                                                  13.00d0, 13.00d0, &
       !                                                                                  log₁₀(α_early)-log₁₀(t₀     )
       &                                                                                  13.00d0, 13.00d0, &
       !                                                                                  log₁₀(α_late )-log₁₀(α_early)
       &                                                                                  13.00d0, 13.00d0, &
       !                                                                                  log₁₀(α_late )-log₁₀(α_late )
       &                                                                                  12.75d0, 14.00d0, &
       !                                                                                  log₁₀(α_late )-log₁₀(t₀     )
       &                                                                                  13.00d0, 13.00d0, &
       !                                                                                  log₁₀(t₀     )-log₁₀(α_early)
       &                                                                                  13.00d0, 13.00d0, &
       !                                                                                  log₁₀(t₀     )-log₁₀(α_late )
       &                                                                                  13.00d0, 13.00d0, &
       !                                                                                  log₁₀(t₀     )-log₁₀(t₀     )
       &                                                                                  13.00d0, 13.00d0  &
       &                                                                                 ]                , &
       &                                                                                 [2,3,3]            &
       &                                                                                )
  
  double precision, dimension(2,3,3), parameter :: rollOverRateLogMassCholesky  =reshape(               &
       &                                                                                 [              &
       !                                                                                  Each row is the roll-over rate of the mass transition for the Cholesky matrix element for:
       !                                                                                   early foming : late forming
       !                                                                                  log₁₀(α_early)-log₁₀(α_early)
       &                                                                                  1.0d0, 0.5d0, &
       !                                                                                  log₁₀(α_early)-log₁₀(α_late )
       &                                                                                  0.5d0, 0.5d0, &
       !                                                                                  log₁₀(α_early)-log₁₀(t₀     )
       &                                                                                  0.5d0, 0.5d0, &
       !                                                                                  log₁₀(α_late )-log₁₀(α_early)
       &                                                                                  0.5d0, 0.5d0, &
       !                                                                                  log₁₀(α_late )-log₁₀(α_late )
       &                                                                                  0.5d0, 1.0d0, &
       !                                                                                  log₁₀(α_late )-log₁₀(t₀     )
       &                                                                                  0.5d0, 0.5d0, &
       !                                                                                  log₁₀(t₀     )-log₁₀(α_early)
       &                                                                                  0.5d0, 0.5d0, &
       !                                                                                  log₁₀(t₀     )-log₁₀(α_late )
       &                                                                                  0.5d0, 0.5d0, &
       !                                                                                  log₁₀(t₀     )-log₁₀(t₀     )
       &                                                                                  0.5d0, 0.5d0  &
       &                                                                                 ]            , &
       &                                                                                 [2,3,3]        &
       &                                                                                )

contains

  function hearin2021StochasticConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily hearin2021Stochastic} dark matter halo mass accretion history class which takes a parameter
    !% set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterHaloMassAccretionHistoryHearin2021Stochastic)                     :: self
    type            (inputParameters                                       ), intent(inout)      :: parameters
    double precision                                                                             :: fractionLateLow, fractionLateHigh, &
         &                                                                                          rateRollOver
    double precision                                                        , dimension(2,2,3  ) :: means
    double precision                                                        , dimension(2,2,3,3) :: cholesky

    ! Population fractions.
    !# <inputParameter>
    !#   <name>fractionLateLow</name>
    !#   <description>The fraction of late-forming halos in the low halo mass limit.</description>
    !#   <source>parameters</source>
    !#   <defaultValue>0.60d0</defaultValue>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>fractionLateHigh</name>
    !#   <description>The fraction of late-forming halos in the high halo mass limit.</description>
    !#   <source>parameters</source>
    !#   <defaultValue>0.62d0</defaultValue>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>rateRollOver</name>
    !#   <description>The roll over rate parameter, $k$.</description>
    !#   <source>parameters</source>
    !#   <defaultValue>3.5d0</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    ! Means.
    !# <inputParameter>
    !#   <name>meanLog10PowerLawIndexEarlyLowEarlyForming</name>
    !#   <description>The mean low-mass limit of $\log_{10}$ early-time power law index for early-forming halos.</description>
    !#   <variable>means(massLow,formationEarly,parameterLog10PowerLawIndexEarly)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>meanDefault(massLow,formationEarly,parameterLog10PowerLawIndexEarly)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>meanLog10PowerLawIndexEarlyHighEarlyForming</name>
    !#   <description>The mean high-mass limit of $\log_{10}$ early-time power law index for early-forming halos.</description>
    !#   <variable>means(massHigh,formationEarly,parameterLog10PowerLawIndexEarly)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>meanDefault(massHigh,formationEarly,parameterLog10PowerLawIndexEarly)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>meanLog10PowerLawIndexEarlyLowLateForming</name>
    !#   <description>The mean low-mass limit of $\log_{10}$ early-time power law index for late-forming halos.</description>
    !#   <variable>means(massLow,formationLate,parameterLog10PowerLawIndexEarly)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>meanDefault(massLow,formationLate,parameterLog10PowerLawIndexEarly)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>meanLog10PowerLawIndexEarlyHighLateForming</name>
    !#   <description>The mean high-mass limit of $\log_{10}$ early-time power law index for late-forming halos.</description>
    !#   <variable>means(massHigh,formationLate,parameterLog10PowerLawIndexEarly)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>meanDefault(massHigh,formationLate,parameterLog10PowerLawIndexEarly)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter> 
    !# <inputParameter>
    !#   <name>meanLog10PowerLawIndexLateLowEarlyForming</name>
    !#   <description>The mean low-mass limit of $\log_{10}$ late-time power law index for early-forming halos.</description>
    !#   <variable>means(massLow,formationEarly,parameterLog10PowerLawIndexLate)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>meanDefault(massLow,formationEarly,parameterLog10PowerLawIndexLate)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>meanLog10PowerLawIndexLateHighEarlyForming</name>
    !#   <description>The mean high-mass limit of $\log_{10}$ late-time power law index for early-forming halos.</description>
    !#   <variable>means(massHigh,formationEarly,parameterLog10PowerLawIndexLate)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>meanDefault(massHigh,formationEarly,parameterLog10PowerLawIndexLate)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>meanLog10PowerLawIndexLateLowLateForming</name>
    !#   <description>The mean low-mass limit of $\log_{10}$ late-time power law index for late-forming halos.</description>
    !#   <variable>means(massLow,formationLate,parameterLog10PowerLawIndexLate)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>meanDefault(massLow,formationLate,parameterLog10PowerLawIndexLate)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>meanLog10PowerLawIndexLateHighLateForming</name>
    !#   <description>The mean high-mass limit of $\log_{10}$ late-time power law index for late-forming halos.</description>
    !#   <variable>means(massHigh,formationLate,parameterLog10PowerLawIndexLate)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>meanDefault(massHigh,formationLate,parameterLog10PowerLawIndexLate)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>meanLog10TimeZeroLowEarlyForming</name>
    !#   <description>The mean low-mass limit of $\log_{10}$ $t_0$ for early-forming halos.</description>
    !#   <variable>means(massLow,formationEarly,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>meanDefault(massLow,formationEarly,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>meanLog10TimeZeroHighEarlyForming</name>
    !#   <description>The mean high-mass limit of $\log_{10}$ $t_0$ for early-forming halos.</description>
    !#   <variable>means(massHigh,formationEarly,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>meanDefault(massHigh,formationEarly,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>meanLog10TimeZeroLowLateForming</name>
    !#   <description>The mean low-mass limit of $\log_{10}$ $t_0$ for late-forming halos.</description>
    !#   <variable>means(massLow,formationLate,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>meanDefault(massLow,formationLate,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>meanLog10TimeZeroHighLateForming</name>
    !#   <description>The mean high-mass limit of $\log_{10}$ $t_0$ for late-forming halos.</description>
    !#   <variable>means(massHigh,formationLate,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>meanDefault(massHigh,formationLate,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    ! Cholesky matrix elements.
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexEarlyLog10PowerLawIndexEarlyLowEarlyForming</name>
    !#   <description>The Cholesky matrix element low-mass limit of $\log_{10}$ early-time power law index for early-forming halos.</description>
    !#   <variable>cholesky(massLow,formationEarly,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexEarly)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massLow,formationEarly,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexEarly)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexEarlyLog10PowerLawIndexEarlyHighEarlyForming</name>
    !#   <description>The Cholesky matrix element high-mass limit of $\log_{10}$ early-time power law index for early-forming halos.</description>
    !#   <variable>cholesky(massHigh,formationEarly,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexEarly)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massHigh,formationEarly,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexEarly)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexEarlyLog10PowerLawIndexEarlyLowLateForming</name>
    !#   <description>The Cholesky matrix element low-mass limit of $\log_{10}$ early-time power law index for late-forming halos.</description>
    !#   <variable>cholesky(massLow,formationLate,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexEarly)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massLow,formationLate,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexEarly)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexEarlyLog10PowerLawIndexEarlyHighLateForming</name>
    !#   <description>The Cholesky matrix element high-mass limit of $\log_{10}$ early-time power law index for late-forming halos.</description>
    !#   <variable>cholesky(massHigh,formationLate,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexEarly)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massHigh,formationLate,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexEarly)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>       
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexEarlyLog10PowerLawIndexLateLowEarlyForming</name>
    !#   <description>The Cholesky matrix element low-mass limit of $\log_{10}$ early-time power law index for early-forming halos.</description>
    !#   <variable>cholesky(massLow,formationEarly,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexLate)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massLow,formationEarly,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexLate)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexEarlyLog10PowerLawIndexLateHighEarlyForming</name>
    !#   <description>The Cholesky matrix element high-mass limit of $\log_{10}$ early-time power law index for early-forming halos.</description>
    !#   <variable>cholesky(massHigh,formationEarly,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexLate)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massHigh,formationEarly,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexLate)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexEarlyLog10PowerLawIndexLateLowLateForming</name>
    !#   <description>The Cholesky matrix element low-mass limit of $\log_{10}$ early-time power law index for late-forming halos.</description>
    !#   <variable>cholesky(massLow,formationLate,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexLate)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massLow,formationLate,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexLate)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexEarlyLog10PowerLawIndexLateHighLateForming</name>
    !#   <description>The Cholesky matrix element high-mass limit of $\log_{10}$ early-time power law index for late-forming halos.</description>
    !#   <variable>cholesky(massHigh,formationLate,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexLate)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massHigh,formationLate,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexLate)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>       
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexEarlyLog10TimeZeroLowEarlyForming</name>
    !#   <description>The Cholesky matrix element low-mass limit of $\log_{10}$ early-time power law index for early-forming halos.</description>
    !#   <variable>cholesky(massLow,formationEarly,parameterLog10PowerLawIndexEarly,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massLow,formationEarly,parameterLog10PowerLawIndexEarly,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexEarlyLog10TimeZeroHighEarlyForming</name>
    !#   <description>The Cholesky matrix element high-mass limit of $\log_{10}$ early-time power law index for early-forming halos.</description>
    !#   <variable>cholesky(massHigh,formationEarly,parameterLog10PowerLawIndexEarly,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massHigh,formationEarly,parameterLog10PowerLawIndexEarly,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexEarlyLog10TimeZeroLowLateForming</name>
    !#   <description>The Cholesky matrix element low-mass limit of $\log_{10}$ early-time power law index for late-forming halos.</description>
    !#   <variable>cholesky(massLow,formationLate,parameterLog10PowerLawIndexEarly,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massLow,formationLate,parameterLog10PowerLawIndexEarly,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexEarlyLog10TimeZeroHighLateForming</name>
    !#   <description>The Cholesky matrix element high-mass limit of $\log_{10}$ early-time power law index for late-forming halos.</description>
    !#   <variable>cholesky(massHigh,formationLate,parameterLog10PowerLawIndexEarly,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massHigh,formationLate,parameterLog10PowerLawIndexEarly,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexLateLog10PowerLawIndexLateLowEarlyForming</name>
    !#   <description>The Cholesky matrix element low-mass limit of $\log_{10}$ early-time power law index for early-forming halos.</description>
    !#   <variable>cholesky(massLow,formationEarly,parameterLog10PowerLawIndexLate,parameterLog10PowerLawIndexLate)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massLow,formationEarly,parameterLog10PowerLawIndexLate,parameterLog10PowerLawIndexLate)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexLateLog10PowerLawIndexLateHighEarlyForming</name>
    !#   <description>The Cholesky matrix element high-mass limit of $\log_{10}$ early-time power law index for early-forming halos.</description>
    !#   <variable>cholesky(massHigh,formationEarly,parameterLog10PowerLawIndexLate,parameterLog10PowerLawIndexLate)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massHigh,formationEarly,parameterLog10PowerLawIndexLate,parameterLog10PowerLawIndexLate)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexLateLog10PowerLawIndexLateLowLateForming</name>
    !#   <description>The Cholesky matrix element low-mass limit of $\log_{10}$ early-time power law index for late-forming halos.</description>
    !#   <variable>cholesky(massLow,formationLate,parameterLog10PowerLawIndexLate,parameterLog10PowerLawIndexLate)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massLow,formationLate,parameterLog10PowerLawIndexLate,parameterLog10PowerLawIndexLate)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexLateLog10PowerLawIndexLateHighLateForming</name>
    !#   <description>The Cholesky matrix element high-mass limit of $\log_{10}$ early-time power law index for late-forming halos.</description>
    !#   <variable>cholesky(massHigh,formationLate,parameterLog10PowerLawIndexLate,parameterLog10PowerLawIndexLate)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massHigh,formationLate,parameterLog10PowerLawIndexLate,parameterLog10PowerLawIndexLate)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexLateLog10TimeZeroLowEarlyForming</name>
    !#   <description>The Cholesky matrix element low-mass limit of $\log_{10}$ early-time power law index for early-forming halos.</description>
    !#   <variable>cholesky(massLow,formationEarly,parameterLog10PowerLawIndexLate,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massLow,formationEarly,parameterLog10PowerLawIndexLate,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexLateLog10TimeZeroHighEarlyForming</name>
    !#   <description>The Cholesky matrix element high-mass limit of $\log_{10}$ early-time power law index for early-forming halos.</description>
    !#   <variable>cholesky(massHigh,formationEarly,parameterLog10PowerLawIndexLate,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massHigh,formationEarly,parameterLog10PowerLawIndexLate,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexLateLog10TimeZeroLowLateForming</name>
    !#   <description>The Cholesky matrix element low-mass limit of $\log_{10}$ early-time power law index for late-forming halos.</description>
    !#   <variable>cholesky(massLow,formationLate,parameterLog10PowerLawIndexLate,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massLow,formationLate,parameterLog10PowerLawIndexLate,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10PowerLawIndexLateLog10TimeZeroHighLateForming</name>
    !#   <description>The Cholesky matrix element high-mass limit of $\log_{10}$ early-time power law index for late-forming halos.</description>
    !#   <variable>cholesky(massHigh,formationLate,parameterLog10PowerLawIndexLate,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massHigh,formationLate,parameterLog10PowerLawIndexLate,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10TimeZeroLog10TimeZeroLowEarlyForming</name>
    !#   <description>The Cholesky matrix element low-mass limit of $\log_{10}$ early-time power law index for early-forming halos.</description>
    !#   <variable>cholesky(massLow,formationEarly,parameterLog10TimeZero,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massLow,formationEarly,parameterLog10TimeZero,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10TimeZeroLog10TimeZeroHighEarlyForming</name>
    !#   <description>The Cholesky matrix element high-mass limit of $\log_{10}$ early-time power law index for early-forming halos.</description>
    !#   <variable>cholesky(massHigh,formationEarly,parameterLog10TimeZero,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massHigh,formationEarly,parameterLog10TimeZero,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10TimeZeroLog10TimeZeroLowLateForming</name>
    !#   <description>The Cholesky matrix element low-mass limit of $\log_{10}$ early-time power law index for late-forming halos.</description>
    !#   <variable>cholesky(massLow,formationLate,parameterLog10TimeZero,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massLow,formationLate,parameterLog10TimeZero,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>choleskyLog10TimeZeroLog10TimeZeroHighLateForming</name>
    !#   <description>The Cholesky matrix element high-mass limit of $\log_{10}$ early-time power law index for late-forming halos.</description>
    !#   <variable>cholesky(massHigh,formationLate,parameterLog10TimeZero,parameterLog10TimeZero)</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>choleskyDefault(massHigh,formationLate,parameterLog10TimeZero,parameterLog10TimeZero)</defaultValue>
    !# <type>real</type>
    !# <cardinality>0..1</cardinality>
    !# </inputParameter>
    ! Symmetrize the Cholesky matrix.
    cholesky(:,:,parameterLog10PowerLawIndexLate,parameterLog10PowerLawIndexEarly)=cholesky(:,:,parameterLog10PowerLawIndexEarly,parameterLog10PowerLawIndexLate)
    cholesky(:,:,parameterLog10TimeZero         ,parameterLog10PowerLawIndexEarly)=cholesky(:,:,parameterLog10PowerLawIndexEarly,parameterLog10TimeZero         )
    cholesky(:,:,parameterLog10TimeZero         ,parameterLog10PowerLawIndexLate )=cholesky(:,:,parameterLog10PowerLawIndexLate ,parameterLog10TimeZero         )
    self=darkMatterHaloMassAccretionHistoryHearin2021Stochastic(rateRollOver,fractionLateLow,fractionLateHigh,means,cholesky)
    !# <inputParametersValidate source="parameters"/>
    return
  end function hearin2021StochasticConstructorParameters

  function hearin2021StochasticConstructorInternal(rateRollOver,fractionLateLow,fractionLateHigh,means,cholesky) result(self)
    !% Internal constructor for the {\normalfont \ttfamily hearin2021Stochastic} dark matter halo mass accretion history class.
    implicit none
    type            (darkMatterHaloMassAccretionHistoryHearin2021Stochastic)                                    :: self
    double precision                                                        , intent(in   )                     :: fractionLateLow, fractionLateHigh, &
          &                                                                                                        rateRollOver
    double precision                                                        , intent(in   ), dimension(2,2,3  ) :: means
    double precision                                                        , intent(in   ), dimension(2,2,3,3) :: cholesky
    !# <constructorAssign variables="rateRollOver, fractionLateLow, fractionLateHigh, means, cholesky"/>
    
    return
  end function hearin2021StochasticConstructorInternal

  subroutine hearin2021StochasticSample(self,node)
    !% Sample parameters from the distribution function for the given {\normalfont \ttfamily node}.
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    use            :: Galacticus_Nodes, only : nodeComponentBasic
    use            :: Linear_Algebra  , only : vector            , matrix, operator(*), assignment(=)
    implicit none
    class           (darkMatterHaloMassAccretionHistoryHearin2021Stochastic), intent(inout)  :: self
    type            (treeNode                                              ), intent(inout)  :: node
    class           (nodeComponentBasic                                    ), pointer        :: basic
    double precision                                                        , dimension(3  ) :: means_      , sample_
    double precision                                                        , dimension(3,3) :: cholesky_
    type            (vector                                                )                 :: means       , sample         , &
         &                                                                                      randoms
    type            (matrix                                                )                 :: cholesky
    double precision                                                                         :: fractionLate, massLogarithmic
    logical                                                                                  :: isLate
    integer                                                                                  :: iFormation  , iParameter     , &
         &                                                                                      jParameter

    ! If parameters were already computed for this tree, simply retrieve them. Otherwise, sample parameters and store them.
    if (node%hostTree%properties%exists('treeMAHHearin2021PowerLawIndexEarly')) then
       self%powerLawIndexEarly  =node%hostTree%properties%value('treeMAHHearin2021PowerLawIndexEarly')
       self%powerLawIndexLate   =node%hostTree%properties%value('treeMAHHearin2021PowerLawIndexLate' )
       self%timeZeroLogarithmic_=node%hostTree%properties%value('treeMAHHearin2021Log10TimeZero'     )
    else
       ! Find late-forming fraction and decide if this tree is late-forming.
       fractionLate=self%fractionLate(node)
       isLate      =node%hostTree%randomNumberGenerator_%uniformSample() <= fractionLate
       ! Evaluate the parameters.
       basic           =>       node %basic()
       massLogarithmic =  log10(basic%mass ())
       if (isLate) then
          iFormation=formationLate
       else
          iFormation=formationEarly
       end if
       ! Compute the means and Cholesky matrix.
       cholesky_=0.0d0
       do iParameter=1,3
          means_        (iParameter         )=self%sigmoid(                                                                             &
               &                                                massLogarithmic                                                       , &
               &                                                log10MassZeroMean          (         iFormation,iParameter           ), &
               &                                                rollOverRateLogMassMean    (         iFormation,iParameter           ), &
               &                                           self%means                      (massLow ,iFormation,iParameter           ), &
               &                                           self%means                      (massHigh,iFormation,iParameter           )  &
               &                                          )
          do jParameter=iParameter,3
             cholesky_(jParameter,iParameter)=self%sigmoid(                                                                             &
                  &                                             massLogarithmic                                                       , &
                  &                                             log10MassZeroCholesky      (         iFormation,iParameter,jParameter), &
                  &                                             rollOverRateLogMassCholesky(         iFormation,iParameter,jParameter), &
                  &                                        self%cholesky                   (massLow ,iFormation,iParameter,jParameter), &
                  &                                        self%cholesky                   (massHigh,iFormation,iParameter,jParameter)  &
                  &                                       )
          end do
       end do
       ! Diagonal elements of the Cholesky matrix are actually log₁₀ of the value.
       do iParameter=1,3
          cholesky_(iParameter,iParameter)=10.0d0**cholesky_(iParameter,iParameter)
       end do
       ! Evaluate the vector of means and the Cholesky matrix.
       means   =vector(means_   )
       cholesky=matrix(cholesky_)
       ! Construct random variables.
       randoms =vector(                                                              &
            &          [                                                             &
            &           node%hostTree%randomNumberGenerator_%standardNormalSample(), &
            &           node%hostTree%randomNumberGenerator_%standardNormalSample(), &
            &           node%hostTree%randomNumberGenerator_%standardNormalSample()  &
            &          ]                                                             &
            &         )
       ! Compute the parameters.
       sample= means    &
            & +cholesky &
            & *randoms
       ! Extract the parameters.
       sample_                  =sample
       self%powerLawIndexEarly  =10.0d0**sample_(1)
       self%powerLawIndexLate   =10.0d0**sample_(2)
       self%timeZeroLogarithmic_=        sample_(3)
       ! Store the parameters for this tree.
       call node%hostTree%properties%set('treeMAHHearin2021PowerLawIndexEarly',self%powerLawIndexEarly  )
       call node%hostTree%properties%set('treeMAHHearin2021PowerLawIndexLate' ,self%powerLawIndexLate   )
       call node%hostTree%properties%set('treeMAHHearin2021Log10TimeZero'     ,self%timeZeroLogarithmic_)
    end if
    return
  end subroutine hearin2021StochasticSample

  double precision function hearin2021StochasticFractionLate(self,node)
    !% Return the fraction of late-forming halos.
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(darkMatterHaloMassAccretionHistoryHearin2021Stochastic), intent(inout) :: self
    type (treeNode                                              ), intent(inout) :: node
    class(nodeComponentBasic                                    ), pointer       :: basic
    
    basic                            => node%basic  (                                                                           )
    hearin2021StochasticFractionLate =  self%sigmoid(log10(basic%mass()),13.0d0,0.5d0,self%fractionLateLow,self%fractionLateHigh)
    return
  end function hearin2021StochasticFractionLate
  
  double precision function hearin2021StochasticPowerLawIndexEarly(self,node)
    !% Return the early power law index for the given node.
    implicit none
    class(darkMatterHaloMassAccretionHistoryHearin2021Stochastic), intent(inout) :: self
    type (treeNode                                              ), intent(inout) :: node

    call self%sample(node)
    hearin2021StochasticPowerLawIndexEarly=self%powerLawIndexEarly
    return
  end function hearin2021StochasticPowerLawIndexEarly
  
  double precision function hearin2021StochasticPowerLawIndexLate(self,node)
    !% Return the late power law index for the given node.
    implicit none
    class(darkMatterHaloMassAccretionHistoryHearin2021Stochastic), intent(inout) :: self
    type (treeNode                                              ), intent(inout) :: node
    
    call self%sample(node)
    hearin2021StochasticPowerLawIndexLate=self%powerLawIndexLate
    return
  end function hearin2021StochasticPowerLawIndexLate
  
  double precision function hearin2021StochasticLog10TimeZero(self,node)
    !% Return the logarithm of the time of zero power law index for the given node.
    implicit none
    class(darkMatterHaloMassAccretionHistoryHearin2021Stochastic), intent(inout) :: self
    type (treeNode                                              ), intent(inout) :: node
    
    call self%sample(node)
    hearin2021StochasticLog10TimeZero=self%timeZeroLogarithmic_
    return
  end function hearin2021StochasticLog10TimeZero
  
  double precision function hearin2021StochasticTimeMaximum(self,node)
    !% Return the time of maximum mass.
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(darkMatterHaloMassAccretionHistoryHearin2021Stochastic), intent(inout) :: self
    type (treeNode                                              ), intent(inout) :: node
    class(nodeComponentBasic                                    ), pointer       :: basic

    basic                           => node %hostTree%baseNode%basic()
    hearin2021StochasticTimeMaximum =  basic%                  time ()
    return
  end function hearin2021StochasticTimeMaximum
