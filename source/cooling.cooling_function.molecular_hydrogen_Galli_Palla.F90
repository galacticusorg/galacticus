!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Implements a cooling function class which implements cooling from molecular hydrogen using the cooling function of
  !% \cite{galli_chemistry_1998}.
  
  !# <coolingFunction name="coolingFunctionMolecularHydrogenGalliPalla">
  !#  <description>Cooling function class which implements cooling from molecular hydrogen using the cooling function of \cite{galli_chemistry_1998}.</description>
  !# </coolingFunction>
  type, extends(coolingFunctionClass) :: coolingFunctionMolecularHydrogenGalliPalla
     !% A cooling function class which implements cooling from molecular hydrogen using the cooling function of \cite{galli_chemistry_1998}.
     private
     ! Indices of "chemical species" (includes atoms, atomic ions and electrons also) used in this cooling function.
     integer          :: atomicHydrogenIndex                                 , electronIndex                                       , &
          &              molecularHydrogenCationIndex                        , molecularHydrogenIndex
     double precision :: coolingFunctionH2PlusElectronTemperatureLogGradient , coolingFunctionH2PlusHTemperatureLogGradient        , &
         &               coolingFunctionLowDensityLimitTemperatureLogGradient, temperaturePrevious1                                , &
         &               temperaturePrevious2                                , temperaturePrevious3                                , &
         &               temperaturePrevious4                                
    double precision  :: coolingFunctionLowDensityLimit                      , coolingFunctionRotationalTemperaturePart            , &
         &               coolingFunctionVibrationalTemperaturePart           , temperatureCommonPrevious                           , &
         &               temperatureHH2PlusPrevious                          , coolingFunctionAtomicHydrogenMolecularHydrogenCation, &
         &               temperatureH2PlusElectronPrevious                   , coolingFunctionElectronMolecularHydrogenCation
   contains
     !@ <objectMethods>
     !@   <object>coolingFunctionMolecularHydrogenGalliPalla</object>
     !@   <objectMethod>
     !@     <method>coolingFunctionH_H2</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ numberDensityHydrogen\argin, \doublezero\ temperature\argin,\textcolor{red}{\textless type(chemicalAbundances)\textgreater} chemicalDensities\argin</arguments>
     !@     <description>Compute the cooling function due to H--H$_2$.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>coolingFunctionH2Plus_Electron</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ temperature\argin,\textcolor{red}{\textless type(chemicalAbundances)\textgreater} chemicalDensities\argin</arguments>
     !@     <description>Compute the cooling function due to H$_2^+$--e$^-$.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>coolingFunctionH_H2Plus</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ temperature\argin,\textcolor{red}{\textless type(chemicalAbundances)\textgreater} chemicalDensities\argin</arguments>
     !@     <description>Compute the cooling function due to H--H$_2^+$.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>commonFactors</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ numberDensityHydrogen\argin, \doublezero\ temperature\argin, \doublezero\ numberDensityCriticalOverNumberDensityHydrogen\argout, \doublezero\ coolingFunctionLocalThermodynamicEquilibrium\argout, \doublezero\ coolingFunctionLowDensityLimit\argout</arguments>
     !@     <description>Compute common factors.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                       molecularHydrogenGalliPallaDestructor
     procedure :: coolingFunction                    => molecularHydrogenGalliPallaCoolingFunction
     procedure :: coolingFunctionTemperatureLogSlope => molecularHydrogenGalliPallaCoolingFunctionTemperatureLogSlope
     procedure :: coolingFunctionDensityLogSlope     => molecularHydrogenGalliPallaCoolingFunctionDensityLogSlope
     procedure :: coolingFunctionH_H2                => molecularHydrogenGalliPallaCoolingFunctionH_H2
     procedure :: coolingFunctionH2Plus_Electron     => molecularHydrogenGalliPallaCoolingFunctionH2Plus_Electron
     procedure :: coolingFunctionH_H2Plus            => molecularHydrogenGalliPallaCoolingFunctionH_H2Plus
     procedure :: commonFactors                      => molecularHydrogenGalliPallaCommonFactors
  end type coolingFunctionMolecularHydrogenGalliPalla

  interface coolingFunctionMolecularHydrogenGalliPalla
     !% Constructors for the ``molecular hydrogen (Galli \& Palla)'' cooling function class.
     module procedure molecularHydrogenGalliPallaConstructorParameters
     module procedure molecularHydrogenGalliPallaConstructorInternal
  end interface coolingFunctionMolecularHydrogenGalliPalla

  ! Parameters for Hollenbach & McKee cooling function fits.
  double precision                , parameter :: molecularHydrogenGalliPallaRotationalLambda1         =9.50d-22
  double precision                , parameter :: molecularHydrogenGalliPallaRotationalLambda2         =3.00d-24
  double precision                , parameter :: molecularHydrogenGalliPallaRotationalTemperature1    =0.13d+00
  double precision                , parameter :: molecularHydrogenGalliPallaRotationalTemperature2    =0.51d+00
  double precision                , parameter :: molecularHydrogenGalliPallaRotationalExponent1       =3.76d+00
  double precision                , parameter :: molecularHydrogenGalliPallaRotationalExponent2       =2.10d+00
  double precision                , parameter :: molecularHydrogenGalliPallaRotationalCoefficient1    =0.12d+00
  double precision                , parameter :: molecularHydrogenGalliPallaVibrationalLambda1        =6.70d-19
  double precision                , parameter :: molecularHydrogenGalliPallaVibrationalLambda2        =1.60d-18
  double precision                , parameter :: molecularHydrogenGalliPallaVibrationalTemperature1   =5.86d+00
  double precision                , parameter :: molecularHydrogenGalliPallaVibrationalTemperature2   =1.17d+01

  ! Parameters for low-density limit cooling function.
  double precision, dimension(0:4), parameter :: molecularHydrogenGalliPallaLowDensityLimitCoefficient=[-103.0000d0,+97.59000d0,-48.05000d+0,+10.8000d0,-0.9032d0]

  ! Parameters for H_2^+ - e^- cooling function.
  double precision, dimension(0:2), parameter :: molecularHydrogenGalliPallaH2PlusElectronCoefficient =[ -33.3299d0, +5.56465d0, -4.67461d-1                     ]

  ! Parameters for H_2^+ - H cooling function.
  double precision, dimension(0:2), parameter :: molecularHydrogenGalliPallaH2PlusHCoefficient        =[ -35.2804d0, +5.86234d0, -5.12276d-1                     ]

contains

  function molecularHydrogenGalliPallaConstructorParameters(parameters)
    !% Constructor for the ``molecular hydrogen (Galli \& Palla)'' cooling function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(coolingFunctionMolecularHydrogenGalliPalla)                :: molecularHydrogenGalliPallaConstructorParameters
    type(inputParameters                           ), intent(in   ) :: parameters
  
    molecularHydrogenGalliPallaConstructorParameters=molecularHydrogenGalliPallaConstructorInternal()
    return
  end function molecularHydrogenGalliPallaConstructorParameters
  
  function molecularHydrogenGalliPallaConstructorInternal()
    !% Internal constructor for the ``molecular hydrogen (Galli \& Palla)'' cooling function class.
    implicit none
    type(coolingFunctionMolecularHydrogenGalliPalla) :: molecularHydrogenGalliPallaConstructorInternal
    
    ! Get the indices of chemicals that will be used.
    molecularHydrogenGalliPallaConstructorInternal%electronIndex               =Chemicals_Index("Electron"               )
    molecularHydrogenGalliPallaConstructorInternal%atomicHydrogenIndex         =Chemicals_Index("AtomicHydrogen"         )
    molecularHydrogenGalliPallaConstructorInternal%molecularHydrogenCationIndex=Chemicals_Index("MolecularHydrogenCation")
    molecularHydrogenGalliPallaConstructorInternal%molecularHydrogenIndex      =Chemicals_Index("MolecularHydrogen"      )
    ! Initialized stored calculations to unphysical values.
    molecularHydrogenGalliPallaConstructorInternal%temperaturePrevious1             =-1.0d0
    molecularHydrogenGalliPallaConstructorInternal%temperaturePrevious2             =-1.0d0
    molecularHydrogenGalliPallaConstructorInternal%temperaturePrevious3             =-1.0d0
    molecularHydrogenGalliPallaConstructorInternal%temperaturePrevious4             =-1.0d0
    molecularHydrogenGalliPallaConstructorInternal%temperatureCommonPrevious        =-1.0d0
    molecularHydrogenGalliPallaConstructorInternal%temperatureHH2PlusPrevious       =-1.0d0
    molecularHydrogenGalliPallaConstructorInternal%temperatureH2PlusElectronPrevious=-1.0d0
    return
  end function molecularHydrogenGalliPallaConstructorInternal
  
  subroutine molecularHydrogenGalliPallaDestructor(self)
    !% Destructor for the ``molecular hydrogen (Galli \& Palla)'' cooling function class.
    implicit none
    type(coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine molecularHydrogenGalliPallaDestructor

  double precision function molecularHydrogenGalliPallaCoolingFunction(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the cooling function due to molecular hydrogen using the cooling function of \cite{galli_chemistry_1998} (which
    !% refers to the local thermodynamic equilibrium cooling function of \cite{hollenbach_molecule_1979}). Cooling functions
    !% involving H$_2^+$ are computed using polynomial fits to the results of \cite{suchkov_cooling_1978} found by Andrew
    !% Benson by measuring curves from the original paper.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self
    double precision                                            , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                                ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances                        ), intent(in   ) :: chemicalDensities
    type            (radiationStructure                        ), intent(in   ) :: radiation

    ! Check if the hydrogen density is positive.
    if (numberDensityHydrogen > 0.0d0) then
       molecularHydrogenGalliPallaCoolingFunction=+self%coolingFunctionH_H2           (numberDensityHydrogen,temperature,chemicalDensities) & ! H   - H₂ cooling function.
            &                                     +self%coolingFunctionH2Plus_Electron(                      temperature,chemicalDensities) & ! H₂⁺ - e⁻ cooling function.
            &                                     +self%coolingFunctionH_H2Plus       (                      temperature,chemicalDensities)   ! H₂⁺ - H  cooling function.
    else
       ! No density, return zero.
       molecularHydrogenGalliPallaCoolingFunction=0.0d0
    end if
    return
  end function molecularHydrogenGalliPallaCoolingFunction

  double precision function molecularHydrogenGalliPallaCoolingFunctionDensityLogSlope(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the gradient with respect to density of the cooling function due to molecular hydrogen using the cooling function
    !% of \cite{galli_chemistry_1998}.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self
    double precision                                            , intent(in   ) :: numberDensityHydrogen         , temperature
    type            (abundances                                ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances                        ), intent(in   ) :: chemicalDensities
    type            (radiationStructure                        ), intent(in   ) :: radiation
    double precision                                                            :: coolingFunction               , coolingFunctionLocalThermodynamicEquilibrium  , &
         &                                                                         coolingFunctionLowDensityLimit, numberDensityCriticalOverNumberDensityHydrogen, &
         &                                                                         coolingFunctionCumulative

    ! Check if the hydrogen density is positive.
    if (numberDensityHydrogen > 0.0d0) then
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

  double precision function molecularHydrogenGalliPallaCoolingFunctionTemperatureLogSlope(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the gradient with respect to temperature of the cooling function due to molecular hydrogen using the cooling
    !% function of \cite{galli_chemistry_1998}.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    use Numerical_Constants_Prefixes
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self
    double precision                                            , intent(in   ) :: numberDensityHydrogen                                          , temperature
    type            (abundances                                ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances                        ), intent(in   ) :: chemicalDensities
    type            (radiationStructure                        ), intent(in   ) :: radiation
    double precision                                                            :: coolingFunction                                                , coolingFunctionLocalThermodynamicEquilibrium   , &
         &                                                                         coolingFunctionLocalThermodynamicEquilibriumTemperatureGradient, coolingFunctionLowDensityLimit                 , &
         &                                                                         coolingFunctionLowDensityLimitTemperatureGradient              , coolingFunctionRotationalTemperatureGradient   , &
         &                                                                         coolingFunctionVibrationalTemperatureGradient                  , logarithmic10Temperature                       , &
         &                                                                         numberDensityCriticalOverNumberDensityHydrogen                 , temperatureThousand                            , &
         &                                                                         coolingFunctionCumulative

    ! Check if the hydrogen density is positive.
    if (numberDensityHydrogen > 0.0d0) then
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
          self%coolingFunctionLowDensityLimitTemperatureLogGradient=                                                        &
               &                    +                          molecularHydrogenGalliPallaLowDensityLimitCoefficient(1)     &
               &                    +logarithmic10Temperature*(molecularHydrogenGalliPallaLowDensityLimitCoefficient(2)     &
               &                    +logarithmic10Temperature*(molecularHydrogenGalliPallaLowDensityLimitCoefficient(3)     &
               &                    +logarithmic10Temperature*(molecularHydrogenGalliPallaLowDensityLimitCoefficient(4))))
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
             coolingFunctionRotationalTemperatureGradient =                                          &
                  & +(                                                                               &
                  &   +  molecularHydrogenGalliPallaRotationalLambda1                                &
                  &   *  temperatureThousand**(molecularHydrogenGalliPallaRotationalExponent1-1.0d0) &
                  &   /(                                                                             &
                  &     +1.0d0                                                                       &
                  &     +molecularHydrogenGalliPallaRotationalCoefficient1                           &
                  &     *temperatureThousand** molecularHydrogenGalliPallaRotationalExponent2        &
                  &    )                                                                             &
                  &  )                                                                               &
                  & *exp(                                                                            &
                  &      -(                                                                          &
                  &        +molecularHydrogenGalliPallaRotationalTemperature1                        &
                  &        /temperatureThousand                                                      &
                  &       )**3                                                                       &
                  &     )                                                                            &
                  & *(                                                                               &
                  &   +3.0d0                                                                         &
                  &   *(                                                                             &
                  &     +molecularHydrogenGalliPallaRotationalTemperature1                           &
                  &     /temperatureThousand                                                         &
                  &    )**3                                                                          &
                  &   +molecularHydrogenGalliPallaRotationalExponent1                                &
                  &   -molecularHydrogenGalliPallaRotationalExponent2                                &
                  &   *molecularHydrogenGalliPallaRotationalCoefficient1                             &
                  &   *temperatureThousand**molecularHydrogenGalliPallaRotationalExponent2           &
                  &   /(                                                                             &
                  &     +1.0d0                                                                       &
                  &     +molecularHydrogenGalliPallaRotationalCoefficient1                           &
                  &     *temperatureThousand**molecularHydrogenGalliPallaRotationalExponent2         &
                  &    )                                                                             &
                  &  )                                                                               &
                  & +molecularHydrogenGalliPallaRotationalLambda2                                    &
                  & /temperatureThousand                                                             &
                  & *molecularHydrogenGalliPallaRotationalTemperature2                               &
                  & /temperatureThousand                                                             &
                  & *exp(                                                                            &
                  &      -molecularHydrogenGalliPallaRotationalTemperature2                          &
                  &      /temperatureThousand                                                        &
                  &     )
             coolingFunctionVibrationalTemperatureGradient=                                          &
                  & +(                                                                               &
                  &   +molecularHydrogenGalliPallaVibrationalLambda1                                 &
                  &   *molecularHydrogenGalliPallaVibrationalTemperature1                            &
                  &   /temperatureThousand                                                           &
                  &   *exp(                                                                          &
                  &        -molecularHydrogenGalliPallaVibrationalTemperature1                       &
                  &        /temperatureThousand                                                      &
                  &       )                                                                          &
                  &   +molecularHydrogenGalliPallaVibrationalLambda2                                 &
                  &   *molecularHydrogenGalliPallaVibrationalTemperature2                            &
                  &   /temperatureThousand                                                           &
                  &   *exp(                                                                          &
                  &        -molecularHydrogenGalliPallaVibrationalTemperature2                       &
                  &        /temperatureThousand                                                      &
                  &       )                                                                          &
                  &  )                                                                               &
                  & /temperatureThousand
             self%temperaturePrevious2=temperature
          end if
          coolingFunctionLocalThermodynamicEquilibriumTemperatureGradient=            &
               & +milli                                                               &
               & *(                                                                   &
               &   +coolingFunctionRotationalTemperatureGradient                      &
               &   +coolingFunctionVibrationalTemperatureGradient                     &
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
             logarithmic10Temperature                                =       &
                  & +log10(temperature)
             self%coolingFunctionH2PlusElectronTemperatureLogGradient=       &
                  & +molecularHydrogenGalliPallaH2PlusElectronCoefficient(1) &
                  & +molecularHydrogenGalliPallaH2PlusElectronCoefficient(2) &
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
             logarithmic10Temperature                         =       &
                  & +log10(temperature)
             self%coolingFunctionH2PlusHTemperatureLogGradient=       &
                  & +molecularHydrogenGalliPallaH2PlusHCoefficient(1) &
                  & +molecularHydrogenGalliPallaH2PlusHCoefficient(2) &
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
    !% Compute the ratio of critical number density to the hydrogen number density for use in molecular hydrogen cooling functions.
    use Numerical_Constants_Prefixes
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self
    double precision                                            , intent(in   ) :: numberDensityHydrogen                         , temperature
    double precision                                            , intent(  out) :: coolingFunctionLocalThermodynamicEquilibrium  , coolingFunctionLowDensityLimit                 , &
         &                                                                         numberDensityCriticalOverNumberDensityHydrogen
    double precision                                                            :: logarithmic10Temperature                      , temperatureThousand

    if (temperature /= self%temperatureCommonPrevious) then
       ! The expression from Galli & Palla (1998), assumes an equilibrium (1:3) ratio of para:ortho.
       logarithmic10Temperature=log10(temperature)
       self%coolingFunctionLowDensityLimit=                                                                     &
            & +10.0d0**(                                                                                        &
            &                                      molecularHydrogenGalliPallaLowDensityLimitCoefficient(0)     &
            &           +logarithmic10Temperature*(molecularHydrogenGalliPallaLowDensityLimitCoefficient(1)     &
            &           +logarithmic10Temperature*(molecularHydrogenGalliPallaLowDensityLimitCoefficient(2)     &
            &           +logarithmic10Temperature*(molecularHydrogenGalliPallaLowDensityLimitCoefficient(3)     &
            &           +logarithmic10Temperature*(molecularHydrogenGalliPallaLowDensityLimitCoefficient(4))))) &
            &                                 )
       ! Rotational and vibrational cooling functions from Hollenbach & McKee (1979; their equations 6.37 and 6.38).
       temperatureThousand                           =+temperature & ! Convert to units of 1,000K.
            &                                         *milli
       self%coolingFunctionRotationalTemperaturePart =                               &
            & +molecularHydrogenGalliPallaRotationalLambda1                          &
            & *temperatureThousand**molecularHydrogenGalliPallaRotationalExponent1   &
            & /(                                                                     &
            &   +1.0d0                                                               &
            &   +molecularHydrogenGalliPallaRotationalCoefficient1                   &
            &   *temperatureThousand**molecularHydrogenGalliPallaRotationalExponent2 &
            &  )                                                                     &
            & *exp(                                                                  &
            &      -(                                                                &
            &        +molecularHydrogenGalliPallaRotationalTemperature1              &
            &        /temperatureThousand                                            &
            &       )**3                                                             &
            &     )                                                                  &
            & +molecularHydrogenGalliPallaRotationalLambda2                          &
            & *exp(                                                                  &
            &      -molecularHydrogenGalliPallaRotationalTemperature2                &
            &      /temperatureThousand                                              &
            &     )
       self%coolingFunctionVibrationalTemperaturePart=                               &
            & +molecularHydrogenGalliPallaVibrationalLambda1                         &
            & *exp(                                                                  &
            &      -molecularHydrogenGalliPallaVibrationalTemperature1               &
            &      /temperatureThousand                                              &
            &     )                                                                  &
            & +molecularHydrogenGalliPallaVibrationalLambda2                         &
            & *exp(                                                                  &
            &      -molecularHydrogenGalliPallaVibrationalTemperature2               &
            &      /temperatureThousand                                              &
            &     )
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
    !% Compute the cooling function due to H$_2^+$--e$^-$ interactions.
    use Chemical_Abundances_Structure
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self
    double precision                                            , intent(in   ) :: temperature
    type            (chemicalAbundances                        ), intent(in   ) :: chemicalDensities
    double precision                                                            :: electronDensity               , logarithmic10Temperature, &
         &                                                                         molecularHydrogenCationDensity

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
               &           +molecularHydrogenGalliPallaH2PlusElectronCoefficient(0)*logarithmic10Temperature**0 &
               &           +molecularHydrogenGalliPallaH2PlusElectronCoefficient(1)*logarithmic10Temperature**1 &
               &           +molecularHydrogenGalliPallaH2PlusElectronCoefficient(2)*logarithmic10Temperature**2 &
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
    !% Compute the cooling function due to H--H$_2^+$ interactions.
    use Chemical_Abundances_Structure
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self
    double precision                                            , intent(in   ) :: temperature
    type            (chemicalAbundances                        ), intent(in   ) :: chemicalDensities
    double precision                                                            :: atomicHydrogenDensity         , logarithmic10Temperature, &
         &                                                                         molecularHydrogenCationDensity

    ! Get the relevant densities.
    atomicHydrogenDensity         =chemicalDensities%abundance(self%atomicHydrogenIndex         )
    molecularHydrogenCationDensity=chemicalDensities%abundance(self%molecularHydrogenCationIndex)
    ! H - H₂⁺ cooling function.
    if (temperature > 1.0d3 .and. temperature < 1.0d4) then
       ! Recompute the temperature dependent part if necessary.
       if (temperature /= self%temperatureHH2PlusPrevious) then
          logarithmic10Temperature                            =                                          &
               & +log10(temperature)
          self%coolingFunctionAtomicHydrogenMolecularHydrogenCation=                                     &
               & +10.0d0**(                                                                              &
               &           +molecularHydrogenGalliPallaH2PlusHCoefficient(0)*logarithmic10Temperature**0 &
               &           +molecularHydrogenGalliPallaH2PlusHCoefficient(1)*logarithmic10Temperature**1 &
               &           +molecularHydrogenGalliPallaH2PlusHCoefficient(2)*logarithmic10Temperature**2 &
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
    !% Compute the cooling function due to H--H$_2$ interactions.
    use Chemical_Abundances_Structure
    implicit none
    class           (coolingFunctionMolecularHydrogenGalliPalla), intent(inout) :: self
    double precision                                            , intent(in   ) :: numberDensityHydrogen                         , temperature
    type            (chemicalAbundances                        ), intent(in   ) :: chemicalDensities
    double precision                                                            :: atomicHydrogenDensity                         , coolingFunctionLocalThermodynamicEquilibrium, &
         &                                                                         coolingFunctionLowDensityLimit                , molecularHydrogenDensity                    , &
         &                                                                         numberDensityCriticalOverNumberDensityHydrogen

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
