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

  !+ Contributions to this file made by: Sachi Weerasooriya

  !!{
  Implements a cooling function class which implements cooling from molecular hydrogen using the cooling function of
  \cite{galli_chemistry_1998}.
  !!}

  use :: Tables, only : table1DLinearLinear

  !![
  <coolingFunction name="coolingFunctionMolecularHydrogenGalliPalla">
   <description>
    A cooling function class that computes the cooling function due to molecular hydrogen using the results of
    \cite{galli_chemistry_1998}. For the H--H$_2$ cooling function, the fitting functions from \cite{galli_chemistry_1998} are
    used. For the H$_2^+$--e$^-$ and H--H$_2^+$ cooling functions fitting functions to the results plotted in
    \cite{suchkov_cooling_1978} are used:
    \begin{equation}
    \log_{10}\left({\Lambda(T) \over \hbox{erg s}^{-1} \hbox{cm}^3}\right) = C_0 + C_1 \log_{10} \left({T\over\hbox{K}}\right)
    + C_2 \left[\log_{10} \left({T\over\hbox{K}}\right)\right]^2,
    \label{eq:H2CoolingFunction}
    \end{equation}
    where the coefficients $C_{0-2}$ are given in Table~\ref{tb:H2CoolingFunctionCoefficients}.
    
    \begin{table}
     \begin{center}
      \caption{Coefficients of H$_2^+$ cooling functions as appearing in the fitting function,
      eq.~\protect\ref{eq:H2CoolingFunction}.}
      \label{tb:H2CoolingFunctionCoefficients}
      \begin{tabular}{lrrr}
       \hline
       &amp; \multicolumn{3}{c}{{\normalfont \bfseries Coefficient}} \\
       {\normalfont \bfseries Interaction} &amp; \boldmath{$C_0$} &amp; \boldmath{$C_1$} &amp; \boldmath{$C_2$} \\
       \hline
       H$_2^+$--e$^-$ &amp; -33.33 &amp; 5.565 &amp; -0.4675 \\
       H--H$_2^+$ &amp; -35.28 &amp; 5.862 &amp; -0.5124 \\
       \hline
      \end{tabular}
     \end{center}
    \end{table}
   </description>
  </coolingFunction>
  !!]
  type, extends(coolingFunctionClass) :: coolingFunctionMolecularHydrogenGalliPalla
     !!{
     A cooling function class which implements cooling from molecular hydrogen using the cooling function of \cite{galli_chemistry_1998}.
     !!}
     private
     ! Indices of "chemical species" (includes atoms, atomic ions and electrons also) used in this cooling function.
     integer                                           :: atomicHydrogenIndex                                  , electronIndex                                       , &
          &                                               molecularHydrogenCationIndex                         , molecularHydrogenIndex
     double precision                                  :: coolingFunctionH2PlusElectronTemperatureLogGradient  , coolingFunctionH2PlusHTemperatureLogGradient        , &
         &                                                coolingFunctionLowDensityLimitTemperatureLogGradient , temperaturePrevious1                                , &
         &                                                temperaturePrevious2                                 , temperaturePrevious3                                , &
         &                                                temperaturePrevious4
    double precision                                   :: coolingFunctionLowDensityLimit                       , coolingFunctionRotationalTemperaturePart            , &
         &                                                coolingFunctionVibrationalTemperaturePart            , temperatureCommonPrevious                           , &
         &                                                temperatureHH2PlusPrevious                           , coolingFunctionAtomicHydrogenMolecularHydrogenCation, &
         &                                                temperatureH2PlusElectronPrevious                    , coolingFunctionElectronMolecularHydrogenCation      , &
         &                                                coolingFunctionRotationalTemperatureGradient         , coolingFunctionVibrationalTemperatureGradient
    double precision                                   :: temperatureMinimumInterpolators                      , temperatureMaximumInterpolators
    type            (table1DLinearLinear), allocatable :: interpolatorCoolingFunctionCommon
   contains
     !![
     <methods>
       <method description="Compute the cooling function due to H--H$_2$."       method="coolingFunctionH_H2"           />
       <method description="Compute the cooling function due to H$_2^+$--e$^-$." method="coolingFunctionH2Plus_Electron"/>
       <method description="Compute the cooling function due to H--H$_2^+$."     method="coolingFunctionH_H2Plus"       />
       <method description="Compute common factors."                             method="commonFactors"                 />
     </methods>
     !!]
     procedure :: coolingFunction                    => molecularHydrogenGalliPallaCoolingFunction
     procedure :: coolingFunctionFractionInBand      => molecularHydrogenGalliPallaCoolingFunctionFractionInBand
     procedure :: coolingFunctionTemperatureLogSlope => molecularHydrogenGalliPallaCoolingFunctionTemperatureLogSlope
     procedure :: coolingFunctionDensityLogSlope     => molecularHydrogenGalliPallaCoolingFunctionDensityLogSlope
     procedure :: coolingFunctionH_H2                => molecularHydrogenGalliPallaCoolingFunctionH_H2
     procedure :: coolingFunctionH2Plus_Electron     => molecularHydrogenGalliPallaCoolingFunctionH2Plus_Electron
     procedure :: coolingFunctionH_H2Plus            => molecularHydrogenGalliPallaCoolingFunctionH_H2Plus
     procedure :: commonFactors                      => molecularHydrogenGalliPallaCommonFactors
  end type coolingFunctionMolecularHydrogenGalliPalla

  interface coolingFunctionMolecularHydrogenGalliPalla
     !!{
     Constructors for the \refClass{coolingFunctionMolecularHydrogenGalliPalla} cooling function class.
     !!}
     module procedure molecularHydrogenGalliPallaConstructorParameters
     module procedure molecularHydrogenGalliPallaConstructorInternal
  end interface coolingFunctionMolecularHydrogenGalliPalla

  ! Parameters for Hollenbach & McKee cooling function fits.
  double precision                , parameter :: rotationalLambda1           =9.50d-22
  double precision                , parameter :: rotationalLambda2           =3.00d-24
  double precision                , parameter :: rotationalTemperature1      =0.13d+00
  double precision                , parameter :: rotationalTemperature2      =0.51d+00
  double precision                , parameter :: rotationalExponent1         =3.76d+00
  double precision                , parameter :: rotationalExponent2         =2.10d+00
  double precision                , parameter :: rotationalCoefficient1      =0.12d+00
  double precision                , parameter :: vibrationalLambda1          =6.70d-19
  double precision                , parameter :: vibrationalLambda2          =1.60d-18
  double precision                , parameter :: vibrationalTemperature1     =5.86d+00
  double precision                , parameter :: vibrationalTemperature2     =1.17d+01

  ! Parameters for low-density limit cooling function.
  double precision, dimension(0:4), parameter :: lowDensityLimitCoefficient  =[-103.0000d0,+97.59000d0,-48.05000d+0,+10.8000d0,-0.9032d0]

  ! Parameters for H₂⁺ - e⁻ cooling function.
  double precision, dimension(0:2), parameter :: H2PlusElectronCoefficient   =[ -33.3299d0, +5.56465d0, -4.67461d-1                     ]

  ! Parameters for H₂⁺ - H cooling function.
  double precision, dimension(0:2), parameter :: H2PlusHCoefficient          =[ -35.2804d0, +5.86234d0, -5.12276d-1                     ]

  ! Maximum temperature for this cooling function. Above this we assume molecular hydrogen will be dissociated, and the rate
  ! coefficients are, in any case, not valid in this regime.
  double precision                , parameter :: temperatureMaximum          =1.00d+06
  
  ! Minimum hydrogen density for which we will attempt to compute cooling functions.
  double precision                , parameter :: numberDensityHydrogenMinimum=1.00d-20
  
contains

  function molecularHydrogenGalliPallaConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{coolingFunctionMolecularHydrogenGalliPalla} cooling function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(coolingFunctionMolecularHydrogenGalliPalla)                :: self
    type(inputParameters                           ), intent(inout) :: parameters

    self=coolingFunctionMolecularHydrogenGalliPalla()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function molecularHydrogenGalliPallaConstructorParameters

  function molecularHydrogenGalliPallaConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{coolingFunctionMolecularHydrogenGalliPalla} cooling function class.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    implicit none
    type   (coolingFunctionMolecularHydrogenGalliPalla) :: self
    integer                                             :: status
    
    ! Get the indices of chemicals that will be used.
    self%electronIndex               =Chemicals_Index("Electron"               ,status)
    self%atomicHydrogenIndex         =Chemicals_Index("AtomicHydrogen"         ,status)
    self%molecularHydrogenCationIndex=Chemicals_Index("MolecularHydrogenCation",status)
    self%molecularHydrogenIndex      =Chemicals_Index("MolecularHydrogen"      ,status)
    ! Initialized stored calculations to unphysical values.
    self%temperaturePrevious1             =-     1.0d0
    self%temperaturePrevious2             =-     1.0d0
    self%temperaturePrevious3             =-     1.0d0
    self%temperaturePrevious4             =-     1.0d0
    self%temperatureCommonPrevious        =-     1.0d0
    self%temperatureHH2PlusPrevious       =-     1.0d0
    self%temperatureH2PlusElectronPrevious=-     1.0d0
    self%temperatureMinimumInterpolators  =+huge(0.0d0)
    self%temperatureMaximumInterpolators  =-huge(0.0d0)
    return
  end function molecularHydrogenGalliPallaConstructorInternal

  double precision function molecularHydrogenGalliPallaCoolingFunction(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the cooling function due to molecular hydrogen using the cooling function of \cite{galli_chemistry_1998} (which
    refers to the local thermodynamic equilibrium cooling function of \cite{hollenbach_molecule_1979}). Cooling functions
    involving H$_2^+$ are computed using polynomial fits to the results of \cite{suchkov_cooling_1978} found by Andrew
    Benson by measuring curves from the original paper.
    !!}
    use :: Abundances_Structure         , only : abundances
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    double precision                                            , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                                ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances                        ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass                       ), intent(inout) :: radiation
    !$GLC attributes unused :: node, radiation, gasAbundances

    ! Check if the hydrogen density is sufficiently large and temperature is below the maximum.
    if (numberDensityHydrogen > numberDensityHydrogenMinimum .and. temperature < temperatureMaximum) then
       molecularHydrogenGalliPallaCoolingFunction=+self%coolingFunctionH_H2           (numberDensityHydrogen,temperature,chemicalDensities) & ! H   - H₂ cooling function.
            &                                     +self%coolingFunctionH2Plus_Electron(                      temperature,chemicalDensities) & ! H₂⁺ - e⁻ cooling function.
            &                                     +self%coolingFunctionH_H2Plus       (                      temperature,chemicalDensities)   ! H₂⁺ - H  cooling function.
    else
       ! No density, return zero.
       molecularHydrogenGalliPallaCoolingFunction=0.0d0
    end if
    return
  end function molecularHydrogenGalliPallaCoolingFunction

  double precision function molecularHydrogenGalliPallaCoolingFunctionFractionInBand(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation,energyLow,energyHigh)
    !!{
    Return the fraction of the cooling function due to emission in the given band. This is currently unsupported.
    !!}
    use :: Abundances_Structure         , only : abundances
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    use :: Error                        , only : Error_Report
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    double precision                                            , intent(in   ) :: numberDensityHydrogen, temperature, &
         &                                                                         energyLow            , energyHigh
    type            (abundances                                ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances                        ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass                       ), intent(inout) :: radiation
    !$GLC attributes unused :: self, node, numberDensityHydrogen, temperature, gasAbundances, chemicalDensities, radiation, energyLow, energyHigh

    molecularHydrogenGalliPallaCoolingFunctionFractionInBand=0.0d0
    call Error_Report('fraction in band is not supported'//{introspection:location})
    return
  end function molecularHydrogenGalliPallaCoolingFunctionFractionInBand

  double precision function molecularHydrogenGalliPallaCoolingFunctionDensityLogSlope(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the gradient with respect to density of the cooling function due to molecular hydrogen using the cooling function
    of \cite{galli_chemistry_1998}.
    !!}
    use :: Abundances_Structure         , only : abundances
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    double precision                                            , intent(in   ) :: numberDensityHydrogen         , temperature
    type            (abundances                                ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances                        ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass                       ), intent(inout) :: radiation
    double precision                                                            :: coolingFunction               , coolingFunctionLocalThermodynamicEquilibrium  , &
         &                                                                         coolingFunctionLowDensityLimit, numberDensityCriticalOverNumberDensityHydrogen, &
         &                                                                         coolingFunctionCumulative
    !$GLC attributes unused :: node, radiation, gasAbundances

    ! Check if the hydrogen density is sufficiently large and temperature is below the maximum.
    if (numberDensityHydrogen > numberDensityHydrogenMinimum .and. temperature < temperatureMaximum) then
       ! Initialize cumulative cooling function to zero.
       coolingFunctionCumulative=0.0d0
       ! H - H₂ cooling function.
       coolingFunction          =+self%coolingFunctionH_H2 (numberDensityHydrogen,temperature,chemicalDensities)
       coolingFunctionCumulative=+coolingFunctionCumulative                                                      &
            &                    +coolingFunction
       call self%commonFactors(numberDensityHydrogen,temperature &
            &,numberDensityCriticalOverNumberDensityHydrogen,coolingFunctionLocalThermodynamicEquilibrium &
            &,coolingFunctionLowDensityLimit)
       if (coolingFunctionLowDensityLimit > 0.0d0) then
          molecularHydrogenGalliPallaCoolingFunctionDensityLogSlope=+coolingFunction                                    &
               &                                                    /numberDensityHydrogen                              &
               &                                                    *(                                                  &
               &                                                      +  1.0d0                                          &
               &                                                      +  numberDensityCriticalOverNumberDensityHydrogen &
               &                                                      /(                                                &
               &                                                        +1.0d0                                          &
               &                                                        +numberDensityCriticalOverNumberDensityHydrogen &
               &                                                       )                                                &
               &                                                     )
       else
          molecularHydrogenGalliPallaCoolingFunctionDensityLogSlope=+2.0d0                 &
               &                                                    *coolingFunction       &
               &                                                    /numberDensityHydrogen
       end if
       ! H₂⁺ - e⁻ cooling function.
       coolingFunction                                          =+self%coolingFunctionH2Plus_Electron(temperature,chemicalDensities)
       coolingFunctionCumulative                                =+coolingFunctionCumulative                                 &
            &                                                    +coolingFunction
       molecularHydrogenGalliPallaCoolingFunctionDensityLogSlope=+molecularHydrogenGalliPallaCoolingFunctionDensityLogSlope &
            &                                                    +2.0d0                                                     &
            &                                                    *coolingFunction                                           &
            &                                                    /numberDensityHydrogen
       ! H - H₂⁺ cooling function.
       coolingFunction                                          =+self%coolingFunctionH_H2Plus(temperature,chemicalDensities)
       coolingFunctionCumulative                                =+coolingFunctionCumulative                                 &
            &                                                    +coolingFunction
       molecularHydrogenGalliPallaCoolingFunctionDensityLogSlope=+molecularHydrogenGalliPallaCoolingFunctionDensityLogSlope &
            &                                                    +2.0d0                                                     &
            &                                                    *coolingFunction                                           &
            &                                                    /numberDensityHydrogen
       ! Convert to logarithmic slope.
       if (coolingFunctionCumulative /= 0.0d0)                                                                                     &
            & molecularHydrogenGalliPallaCoolingFunctionDensityLogSlope=+molecularHydrogenGalliPallaCoolingFunctionDensityLogSlope &
            &                                                           *numberDensityHydrogen                                     &
            &                                                           /coolingFunctionCumulative
    else
       ! No density, return zero.
       molecularHydrogenGalliPallaCoolingFunctionDensityLogSlope=+0.0d0
    end if
    return
  end function molecularHydrogenGalliPallaCoolingFunctionDensityLogSlope

  double precision function molecularHydrogenGalliPallaCoolingFunctionTemperatureLogSlope(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the gradient with respect to temperature of the cooling function due to molecular hydrogen using the cooling
    function of \cite{galli_chemistry_1998}.
    !!}
    use :: Abundances_Structure         , only : abundances
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    use :: Numerical_Constants_Prefixes , only : milli
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    double precision                                            , intent(in   ) :: numberDensityHydrogen                                          , temperature
    type            (abundances                                ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances                        ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass                       ), intent(inout) :: radiation
    double precision                                                            :: coolingFunction                                                , coolingFunctionLocalThermodynamicEquilibrium   , &
         &                                                                         coolingFunctionLocalThermodynamicEquilibriumTemperatureGradient, coolingFunctionLowDensityLimit                 , &
         &                                                                         coolingFunctionLowDensityLimitTemperatureGradient              , logarithmic10Temperature                       , &
         &                                                                         numberDensityCriticalOverNumberDensityHydrogen                 , temperatureThousand                            , &
         &                                                                         coolingFunctionCumulative
    !$GLC attributes unused :: node, gasAbundances, radiation

    ! Check if the hydrogen density is sufficiently large and temperature is below the maximum.
    if (numberDensityHydrogen > numberDensityHydrogenMinimum .and. temperature < temperatureMaximum) then
       coolingFunctionCumulative=0.0d0
       ! H - H₂ cooling function.
       coolingFunction          =+self%coolingFunctionH_H2(numberDensityHydrogen,temperature,chemicalDensities)
       coolingFunctionCumulative=+coolingFunctionCumulative &
            &                    +coolingFunction
       call self%commonFactors(                                                &
            &                  numberDensityHydrogen                         , &
            &                  temperature                                   , &
            &                  numberDensityCriticalOverNumberDensityHydrogen, &
            &                  coolingFunctionLocalThermodynamicEquilibrium  , &
            &                  coolingFunctionLowDensityLimit                  &
            &                 )
       ! Check if we need to recompute the low density limit cooling function gradient.
       if (temperature /= self%temperaturePrevious1) then
          logarithmic10Temperature=log10(temperature)
          self%coolingFunctionLowDensityLimitTemperatureLogGradient=                             &
               &                    +                          lowDensityLimitCoefficient(1)     &
               &                    +logarithmic10Temperature*(lowDensityLimitCoefficient(2)     &
               &                    +logarithmic10Temperature*(lowDensityLimitCoefficient(3)     &
               &                    +logarithmic10Temperature*(lowDensityLimitCoefficient(4))))
          self%temperaturePrevious1=temperature
       end if
       if (coolingFunctionLowDensityLimit > 0.0d0) then
          ! Get the temperature gradient of the low density limit cooling function.
          coolingFunctionLowDensityLimitTemperatureGradient=+coolingFunctionLowDensityLimit                            &
               &                                            /temperature                                               &
               &                                            *self%coolingFunctionLowDensityLimitTemperatureLogGradient
          ! Check if we need to compute gradients of rotational and vibrational cooling functions.
          if (temperature /= self%temperaturePrevious2) then
             ! Get the temperature gradient of the LTE cooling function.
             temperatureThousand                          =+temperature & ! Convert to units of 1,000K.
                  &                                        *milli
             self%coolingFunctionRotationalTemperatureGradient =          &
                  & +(                                                    &
                  &   +  rotationalLambda1                                &
                  &   *  temperatureThousand**(rotationalExponent1-1.0d0) &
                  &   /(                                                  &
                  &     +1.0d0                                            &
                  &     +rotationalCoefficient1                           &
                  &     *temperatureThousand** rotationalExponent2        &
                  &    )                                                  &
                  &  )                                                    &
                  & *exp(                                                 &
                  &      -(                                               &
                  &        +rotationalTemperature1                        &
                  &        /temperatureThousand                           &
                  &       )**3                                            &
                  &     )                                                 &
                  & *(                                                    &
                  &   +3.0d0                                              &
                  &   *(                                                  &
                  &     +rotationalTemperature1                           &
                  &     /temperatureThousand                              &
                  &    )**3                                               &
                  &   +rotationalExponent1                                &
                  &   -rotationalExponent2                                &
                  &   *rotationalCoefficient1                             &
                  &   *temperatureThousand**rotationalExponent2           &
                  &   /(                                                  &
                  &     +1.0d0                                            &
                  &     +rotationalCoefficient1                           &
                  &     *temperatureThousand**rotationalExponent2         &
                  &    )                                                  &
                  &  )                                                    &
                  & +rotationalLambda2                                    &
                  & /temperatureThousand                                  &
                  & *rotationalTemperature2                               &
                  & /temperatureThousand                                  &
                  & *exp(                                                 &
                  &      -rotationalTemperature2                          &
                  &      /temperatureThousand                             &
                  &     )
             self%coolingFunctionVibrationalTemperatureGradient=          &
                  & +(                                                    &
                  &   +vibrationalLambda1                                 &
                  &   *vibrationalTemperature1                            &
                  &   /temperatureThousand                                &
                  &   *exp(                                               &
                  &        -vibrationalTemperature1                       &
                  &        /temperatureThousand                           &
                  &       )                                               &
                  &   +vibrationalLambda2                                 &
                  &   *vibrationalTemperature2                            &
                  &   /temperatureThousand                                &
                  &   *exp(                                               &
                  &        -vibrationalTemperature2                       &
                  &        /temperatureThousand                           &
                  &       )                                               &
                  &  )                                                    &
                  & /temperatureThousand
             self%temperaturePrevious2=temperature
          end if
          coolingFunctionLocalThermodynamicEquilibriumTemperatureGradient=            &
               & +milli                                                               &
               & *(                                                                   &
               &   +self%coolingFunctionRotationalTemperatureGradient                 &
               &   +self%coolingFunctionVibrationalTemperatureGradient                &
               &  )                                                                   &
               & /numberDensityHydrogen
          ! Get the temperature gradient of the cooling function.
          molecularHydrogenGalliPallaCoolingFunctionTemperatureLogSlope  =            &
               & +coolingFunction                                                     &
               & *(                                                                   &
               &   +coolingFunctionLocalThermodynamicEquilibriumTemperatureGradient   &
               &   /coolingFunctionLocalThermodynamicEquilibrium                      &
               &   -(                                                                 &
               &     +1.0d0                                                           &
               &     /coolingFunctionLowDensityLimit                                  &
               &     /(                                                               &
               &       +1.0d0                                                         &
               &       +numberDensityCriticalOverNumberDensityHydrogen                &
               &      )                                                               &
               &    )                                                                 &
               &   *(                                                                 &
               &     +coolingFunctionLocalThermodynamicEquilibriumTemperatureGradient &
               &     -numberDensityCriticalOverNumberDensityHydrogen                  &
               &     *coolingFunctionLowDensityLimitTemperatureGradient               &
               &    )                                                                 &
               &  )
       else
          molecularHydrogenGalliPallaCoolingFunctionTemperatureLogSlope  =            &
               & +coolingFunction                                                     &
               & *self%coolingFunctionLowDensityLimitTemperatureLogGradient           &
               & /temperature
       end if
       ! H₂⁺ - e⁻ cooling function.
       if (temperature > 1.0d3 .and. temperature < 1.0d4) then
          coolingFunction          =+self%coolingFunctionH2Plus_Electron(temperature,chemicalDensities)
          coolingFunctionCumulative=+coolingFunctionCumulative &
               &                    +coolingFunction
          if (temperature /= self%temperaturePrevious3) then
             logarithmic10Temperature                                = &
                  & +log10(temperature)
             self%coolingFunctionH2PlusElectronTemperatureLogGradient= &
                  & +H2PlusElectronCoefficient(1)                      &
                  & +H2PlusElectronCoefficient(2)                      &
                  & *logarithmic10Temperature
             self%temperaturePrevious3=temperature
          end if
          molecularHydrogenGalliPallaCoolingFunctionTemperatureLogSlope=        &
               & +molecularHydrogenGalliPallaCoolingFunctionTemperatureLogSlope &
               & +coolingFunction                                               &
               & /temperature                                                   &
               & *self%coolingFunctionH2PlusElectronTemperatureLogGradient
       end if
       ! H - H₂⁺ cooling function.
       if (temperature > 1.0d3 .and. temperature < 1.0d4) then
          coolingFunction          =+self%coolingFunctionH_H2Plus(temperature,chemicalDensities)
          coolingFunctionCumulative=+coolingFunctionCumulative &
               &                    +coolingFunction
          if (temperature /= self%temperaturePrevious4) then
             logarithmic10Temperature                         = &
                  & +log10(temperature)
             self%coolingFunctionH2PlusHTemperatureLogGradient= &
                  & +H2PlusHCoefficient(1)                      &
                  & +H2PlusHCoefficient(2)                      &
                  & *logarithmic10Temperature
             self%temperaturePrevious4=temperature
          end if
          molecularHydrogenGalliPallaCoolingFunctionTemperatureLogSlope=        &
               & +molecularHydrogenGalliPallaCoolingFunctionTemperatureLogSlope &
               & +coolingFunction                                               &
               & /temperature                                                   &
               & *self%coolingFunctionH2PlusHTemperatureLogGradient
       end if
       ! Convert to logarithmic slope.
       if (coolingFunctionCumulative /= 0.0d0)                                &
            & molecularHydrogenGalliPallaCoolingFunctionTemperatureLogSlope=  &
            &  +molecularHydrogenGalliPallaCoolingFunctionTemperatureLogSlope &
            &  *temperature                                                   &
            &  /coolingFunctionCumulative
    else
       ! No density, return zero.
       molecularHydrogenGalliPallaCoolingFunctionTemperatureLogSlope=0.0d0
    end if
    return
  end function molecularHydrogenGalliPallaCoolingFunctionTemperatureLogSlope

  subroutine molecularHydrogenGalliPallaCommonFactors(self,numberDensityHydrogen,temperature,numberDensityCriticalOverNumberDensityHydrogen,coolingFunctionLocalThermodynamicEquilibrium,coolingFunctionLowDensityLimit)
    !!{
    Compute the ratio of critical number density to the hydrogen number density for use in molecular hydrogen cooling functions.
    !!}
    use :: Numerical_Constants_Prefixes, only : milli
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout)             :: self
    double precision                                            , intent(in   )             :: numberDensityHydrogen                                  , temperature
    double precision                                            , intent(  out)             :: coolingFunctionLocalThermodynamicEquilibrium           , coolingFunctionLowDensityLimit, &
         &                                                                                     numberDensityCriticalOverNumberDensityHydrogen
    integer                                                     , parameter                 :: temperaturePointsPerDecade                    =100
    double precision                                            , parameter                 :: coolingFunctionMinimum                        =1.0d-300
    double precision                                            , allocatable, dimension(:) :: logTemperatures                                        , temperaturesThousand
    integer                                                                                 :: countTemperatures
    double precision                                                                        :: logarithmic10Temperature

    ! Check if solutions must be updated.
    if (temperature /= self%temperatureCommonPrevious) then
       ! Build interpolation tables if necessary.
       if     (                                                    &
            &   temperature < self%temperatureMinimumInterpolators &
            &  .or.                                                &
            &   temperature > self%temperatureMaximumInterpolators &
            & ) then
          if (allocated(self%interpolatorCoolingFunctionCommon)) deallocate(self%interpolatorCoolingFunctionCommon)
          allocate(self%interpolatorCoolingFunctionCommon)
          self%temperatureMinimumInterpolators=min(self%temperatureMinimumInterpolators,temperature/2.0d0)
          self%temperatureMaximumInterpolators=max(self%temperatureMaximumInterpolators,temperature*2.0d0)
          countTemperatures=int(log10(self%temperatureMaximumInterpolators/self%temperatureMinimumInterpolators)*dble(temperaturePointsPerDecade))+1
          call self%interpolatorCoolingFunctionCommon%create(log10(self%temperatureMinimumInterpolators),log10(self%temperatureMaximumInterpolators),countTemperatures,tableCount=3)
          logTemperatures     =self%interpolatorCoolingFunctionCommon%xs()
          temperaturesThousand=10.0d0**logTemperatures/1000.0d0
          ! The expression from Galli & Palla (1998), assumes an equilibrium (1:3) ratio of para:ortho.
          call self%interpolatorCoolingFunctionCommon%populate(                                                                  &
               &                                               log(                                                              &
               &                                                   10.0d0**(                                                     &
               &                                                                              lowDensityLimitCoefficient(0)      &
               &                                                            +logTemperatures*(lowDensityLimitCoefficient(1)      &
               &                                                            +logTemperatures*(lowDensityLimitCoefficient(2)      &
               &                                                            +logTemperatures*(lowDensityLimitCoefficient(3)      &
               &                                                            +logTemperatures*(lowDensityLimitCoefficient(4)))))  &
               &                                                           )                                                     &
               &                                                  )                                                            , &
               &                                               table=1                                                           &
               &                                              )
          ! Rotational and vibrational cooling functions from Hollenbach & McKee (1979; their equations 6.37 and 6.38).          
          call self%interpolatorCoolingFunctionCommon%populate(                                                                  &
               &                                               log(                                                              &
               &                                                   max(                                                          &
               &                                                       +coolingFunctionMinimum                                 , &
               &                                                       +rotationalLambda1                                        &
               &                                                       *temperaturesThousand**rotationalExponent1                &
               &                                                       /(                                                        &
               &                                                         +1.0d0                                                  &
               &                                                         +rotationalCoefficient1                                 &
               &                                                         *temperaturesThousand**rotationalExponent2              &
               &                                                        )                                                        &
               &                                                       *exp(                                                     &
               &                                                            -(                                                   &
               &                                                              +rotationalTemperature1                            &
               &                                                              /temperaturesThousand                              &
               &                                                             )**3                                                &
               &                                                           )                                                     &
               &                                                       +rotationalLambda2                                        &
               &                                                       *exp(                                                     &
               &                                                            -rotationalTemperature2                              &
               &                                                            /temperaturesThousand                                &
               &                                                           )                                                     &
               &                                                      )                                                          &
               &                                                  )                                                            , &
               &                                               table=2                                                           &
               &                                              )
          call self%interpolatorCoolingFunctionCommon%populate(                                                                  &
               &                                               log(                                                              &
               &                                                   max(                                                          &
               &                                                       +coolingFunctionMinimum                                 , &
               &                                                       +vibrationalLambda1                                       &
               &                                                       *exp(                                                     &
               &                                                            -vibrationalTemperature1                             &
               &                                                            /temperaturesThousand                                &
               &                                                           )                                                     &
               &                                                       +vibrationalLambda2                                       &
               &                                                       *exp(                                                     &
               &                                                            -vibrationalTemperature2                             &
               &                                                            /temperaturesThousand                                &
               &                                                           )                                                     &
               &                                                      )                                                          &
               &                                                  )                                                            , &                                             
               &                                               table=3                                                           &
               &                                              )
       end if
       ! Interpolate in tabulated results.
       logarithmic10Temperature                      =log10(                                                                temperature         )
       self%coolingFunctionLowDensityLimit           =exp  (self%interpolatorCoolingFunctionCommon%interpolate(logarithmic10Temperature,table=1))
       self%coolingFunctionRotationalTemperaturePart =exp  (self%interpolatorCoolingFunctionCommon%interpolate(logarithmic10Temperature,table=2))
       self%coolingFunctionVibrationalTemperaturePart=exp  (self%interpolatorCoolingFunctionCommon%interpolate(logarithmic10Temperature,table=3))
       ! Record the temperature used.
       self%temperatureCommonPrevious=temperature
    end if
    ! Get the low density limit cooling function.
    coolingFunctionLowDensityLimit              =+self%coolingFunctionLowDensityLimit
    ! Compute the LTE cooling function.
    coolingFunctionLocalThermodynamicEquilibrium=+(                                                &
         &                                         +self%coolingFunctionRotationalTemperaturePart  &
         &                                         +self%coolingFunctionVibrationalTemperaturePart &
         &                                        )                                                &
         &                                       /numberDensityHydrogen
    ! Compute the number density ratio.
    if (coolingFunctionLowDensityLimit > 0.0d0) then
       numberDensityCriticalOverNumberDensityHydrogen=+coolingFunctionLocalThermodynamicEquilibrium &
            &                                         /coolingFunctionLowDensityLimit
    else
       ! Unphysical value, should not be used.
       numberDensityCriticalOverNumberDensityHydrogen=-1.0d0
    end if
    return
  end subroutine molecularHydrogenGalliPallaCommonFactors

  double precision function molecularHydrogenGalliPallaCoolingFunctionH2Plus_Electron(self,temperature,chemicalDensities)
    !!{
    Compute the cooling function due to H$_2^+$--e$^-$ interactions.
    !!}
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self
    double precision                                            , intent(in   ) :: temperature
    type            (chemicalAbundances                        ), intent(in   ) :: chemicalDensities
    double precision                                                            :: electronDensity               , logarithmic10Temperature, &
         &                                                                         molecularHydrogenCationDensity

    ! Set to zero if species are not present.
    if     (                                       &
         &   self%molecularHydrogenCationIndex < 0 &
         &  .or.                                   &
         &   self%electronIndex                < 0 &
         & ) then
       molecularHydrogenGalliPallaCoolingFunctionH2Plus_Electron=0.0d0
       return
    end if
    ! Get the relevant densities.
    electronDensity               =chemicalDensities%abundance(self%electronIndex               )
    molecularHydrogenCationDensity=chemicalDensities%abundance(self%molecularHydrogenCationIndex)
    ! H₂⁺ - e⁻ cooling function.
    if (temperature > 1.0d3 .and. temperature < 1.0d4) then
       ! Recompute the temperature dependent part if necessary.
       if (temperature /= self%temperatureH2PlusElectronPrevious) then
          logarithmic10Temperature                           = &
               & +log10(temperature)
          self%coolingFunctionElectronMolecularHydrogenCation=                                                  &
               & +10.0d0**(                                                                                     &
               &           +H2PlusElectronCoefficient(0)*logarithmic10Temperature**0 &
               &           +H2PlusElectronCoefficient(1)*logarithmic10Temperature**1 &
               &           +H2PlusElectronCoefficient(2)*logarithmic10Temperature**2 &
               &          )
          self%temperatureH2PlusElectronPrevious=temperature
       end if
       molecularHydrogenGalliPallaCoolingFunctionH2Plus_Electron=+self%coolingFunctionElectronMolecularHydrogenCation &
            &                                                    *molecularHydrogenCationDensity                      &
            &                                                    *electronDensity
    else
       molecularHydrogenGalliPallaCoolingFunctionH2Plus_Electron=0.0d0
    end if

    return
  end function molecularHydrogenGalliPallaCoolingFunctionH2Plus_Electron

  double precision function molecularHydrogenGalliPallaCoolingFunctionH_H2Plus(self,temperature,chemicalDensities)
    !!{
    Compute the cooling function due to H--H$_2^+$ interactions.
    !!}
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self
    double precision                                            , intent(in   ) :: temperature
    type            (chemicalAbundances                        ), intent(in   ) :: chemicalDensities
    double precision                                                            :: atomicHydrogenDensity         , logarithmic10Temperature, &
         &                                                                         molecularHydrogenCationDensity

    ! Set to zero if species are not present.
    if     (                                       &
         &   self%molecularHydrogenCationIndex < 0 &
         &  .or.                                   &
         &   self%atomicHydrogenIndex          < 0 &
         & ) then
       molecularHydrogenGalliPallaCoolingFunctionH_H2Plus=0.0d0
       return
    end if
    ! Get the relevant densities.
    atomicHydrogenDensity         =chemicalDensities%abundance(self%atomicHydrogenIndex         )
    molecularHydrogenCationDensity=chemicalDensities%abundance(self%molecularHydrogenCationIndex)
    ! H - H₂⁺ cooling function.
    if (temperature > 1.0d3 .and. temperature < 1.0d4) then
       ! Recompute the temperature dependent part if necessary.
       if (temperature /= self%temperatureHH2PlusPrevious) then
          logarithmic10Temperature                            =               &
               & +log10(temperature)
          self%coolingFunctionAtomicHydrogenMolecularHydrogenCation=          &
               & +10.0d0**(                                                   &
               &           +H2PlusHCoefficient(0)*logarithmic10Temperature**0 &
               &           +H2PlusHCoefficient(1)*logarithmic10Temperature**1 &
               &           +H2PlusHCoefficient(2)*logarithmic10Temperature**2 &
               &          )
          self%temperatureHH2PlusPrevious=temperature
       end if
       molecularHydrogenGalliPallaCoolingFunctionH_H2Plus=+self%coolingFunctionAtomicHydrogenMolecularHydrogenCation &
            &                                             *atomicHydrogenDensity                                     &
            &                                             *molecularHydrogenCationDensity
    else
       molecularHydrogenGalliPallaCoolingFunctionH_H2Plus=+0.0d0
    end if
    return
  end function molecularHydrogenGalliPallaCoolingFunctionH_H2Plus

  double precision function molecularHydrogenGalliPallaCoolingFunctionH_H2(self,numberDensityHydrogen,temperature,chemicalDensities)
    !!{
    Compute the cooling function due to H--H$_2$ interactions.
    !!}
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self
    double precision                                            , intent(in   ) :: numberDensityHydrogen                         , temperature
    type            (chemicalAbundances                        ), intent(in   ) :: chemicalDensities
    double precision                                                            :: atomicHydrogenDensity                         , coolingFunctionLocalThermodynamicEquilibrium, &
         &                                                                         coolingFunctionLowDensityLimit                , molecularHydrogenDensity                    , &
         &                                                                         numberDensityCriticalOverNumberDensityHydrogen

    ! Set to zero if species are not present.
    if     (                                 &
         &   self%molecularHydrogenIndex < 0 &
         &  .or.                             &
         &   self%atomicHydrogenIndex    < 0 &
         & ) then
       molecularHydrogenGalliPallaCoolingFunctionH_H2=0.0d0
       return
    end if
    ! Get the relevant densities.
    atomicHydrogenDensity   =chemicalDensities%abundance(self%atomicHydrogenIndex   )
    molecularHydrogenDensity=chemicalDensities%abundance(self%molecularHydrogenIndex)
    ! Compute the cooling function.
    call self%commonFactors(                                                &
         &                  numberDensityHydrogen                         , &
         &                  temperature                                   , &
         &                  numberDensityCriticalOverNumberDensityHydrogen, &
         &                  coolingFunctionLocalThermodynamicEquilibrium  , &
         &                  coolingFunctionLowDensityLimit                  &
         &                 )
    if (coolingFunctionLowDensityLimit > 0.0d0) then
       molecularHydrogenGalliPallaCoolingFunctionH_H2=+coolingFunctionLocalThermodynamicEquilibrium     &
            &                                         /(                                                &
            &                                           +1.0d0                                          &
            &                                           +numberDensityCriticalOverNumberDensityHydrogen &
            &                                          )
    else
       molecularHydrogenGalliPallaCoolingFunctionH_H2=+coolingFunctionLowDensityLimit
    end if
    ! Scale to the density of hydrogen and molecular hydrogen.
    molecularHydrogenGalliPallaCoolingFunctionH_H2=+molecularHydrogenGalliPallaCoolingFunctionH_H2 &
         &                                         *molecularHydrogenDensity                       &
         &                                         *atomicHydrogenDensity
    return
  end function molecularHydrogenGalliPallaCoolingFunctionH_H2
