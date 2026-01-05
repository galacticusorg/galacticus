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

  !!{
  An implementation of Gaunt factors using the \cite{sutherland_accurate_1998} fitting function.
  !!}

  use :: Atomic_Ionization_Potentials, only : atomicIonizationPotentialClass

  !![
  <gauntFactor name="gauntFactorSutherland1998">
   <description>Gaunt factors are computed using the fitting function of \cite{sutherland_accurate_1998}.</description>
  </gauntFactor>
  !!]
  type, extends(gauntFactorClass) :: gauntFactorSutherland1998
     !!{
     A gaunt factor class implementing the fitting function of \cite{sutherland_accurate_1998}.
     !!}
     private
     class(atomicIonizationPotentialClass), pointer :: atomicIonizationPotential_ => null()
   contains
     final     ::          sutherland1998Destructor
     procedure :: total => sutherland1998Total
  end type gauntFactorSutherland1998

  interface gauntFactorSutherland1998
     !!{
     Constructors for the \refClass{gauntFactorSutherland1998} gaunt factor class.
     !!}
     module procedure sutherland1998ConstructorParameters
     module procedure sutherland1998ConstructorInternal
  end interface gauntFactorSutherland1998

  ! Arrays to hold coefficients of the fitting function.
  integer         , parameter                   :: countCoefficient       =  41
  double precision, dimension(countCoefficient) :: energyScalesLogarithmic=[                                                                 &
       &                                                                    -4.000000000d0, -3.80000000d0, -3.6000000000d0, -3.4000000000d0, &
       &                                                                    -3.200000000d0, -3.00000000d0, -2.8000000000d0, -2.6000000000d0, &
       &                                                                    -2.400000000d0, -2.20000000d0, -2.0000000000d0, -1.8000000000d0, &
       &                                                                    -1.600000000d0, -1.40000000d0, -1.2000000000d0, -1.0000000000d0, &
       &                                                                    -0.800000000d0, -0.60000000d0, -0.4000000000d0, -0.2000000000d0, &
       &                                                                    +0.000000000d0, +0.20000000d0, +0.4000000000d0, +0.6000000000d0, &
       &                                                                    +0.800000000d0, +1.00000000d0, +1.2000000000d0, +1.4000000000d0, &
       &                                                                    +1.600000000d0, +1.80000000d0, +2.0000000000d0, +2.2000000000d0, &
       &                                                                    +2.400000000d0, +2.60000000d0, +2.8000000000d0, +3.0000000000d0, &
       &                                                                    +3.200000000d0, +3.40000000d0, +3.6000000000d0, +3.8000000000d0, &
       &                                                                    +4.000000000d0                                                   &
       &                                                                   ]                                                               , &
       &                                           gauntFactorZeroPoint   =[                                                                 &
       &                                                                    +1.113883000d0, +1.11698300d0, +1.1208910000d0, +1.1258100000d0, &
       &                                                                    +1.131995000d0, +1.13975200d0, +1.1494450000d0, +1.1614930000d0, &
       &                                                                    +1.176350000d0, +1.19446700d0, +1.2162200000d0, +1.2418220000d0, &
       &                                                                    +1.271043000d0, +1.30327700d0, +1.3371140000d0, +1.3703950000d0, &
       &                                                                    +1.400285000d0, +1.42365000d0, +1.4376770000d0, +1.4405960000d0, &
       &                                                                    +1.432199000d0, +1.41390900d0, +1.3883010000d0, +1.3583240000d0, &
       &                                                                    +1.326582000d0, +1.29496100d0, +1.2646210000d0, +1.2361820000d0, &
       &                                                                    +1.209926000d0, +1.18594400d0, +1.1642110000d0, +1.1446400000d0, &
       &                                                                    +1.127106000d0, +1.11146600d0, +1.0975660000d0, +1.0852550000d0, &
       &                                                                    +1.074382000d0, +1.06480700d0, +1.0563990000d0, +1.0490420000d0, &
       &                                                                    +1.042637000d0                                                   &
       &                                                                   ]                                                               , &
       &                                           fitCoefficient1        =[                                                                 &
       &                                                                    +0.013480000d0, +0.017445800d0, +0.021856800d0, +0.0275320100d0, &
       &                                                                    +0.034575150d0, +0.043297370d0, +0.053985360d0, +0.0668761900d0, &
       &                                                                    +0.082084880d0, +0.099394300d0, +0.118387900d0, +0.1373789000d0, &
       &                                                                    +0.154441300d0, +0.166680800d0, +0.169900300d0, +0.1604878000d0, &
       &                                                                    +0.135713600d0, +0.095482930d0, +0.043234600d0, -0.0142314700d0, &
       &                                                                    -0.068478680d0, -0.112158900d0, -0.141355900d0, -0.1561924000d0, &
       &                                                                    -0.159659600d0, -0.155614300d0, -0.147298200d0, -0.1368782000d0, &
       &                                                                    -0.125614100d0, -0.114235600d0, -0.103168800d0, -0.0926491700d0, &
       &                                                                    -0.082809610d0, -0.073722290d0, -0.065401160d0, -0.0578381300d0, &
       &                                                                    -0.051006360d0, -0.044856500d0, -0.039312590d0, -0.0343681100d0, &
       &                                                                    -0.029644990d0                                                   &
       &                                                                  ]                                                                , &
       &                                           fitCoefficient2       =[                                                                  &
       &                                                                   +0.0104710000d0, +0.009358014d0, +0.012696960d0, +0.0156791100d0, &
       &                                                                   +0.0195366000d0, +0.024074480d0, +0.029365470d0, +0.0350886800d0, &
       &                                                                   +0.0409547600d0, +0.045592330d0, +0.049375890d0, +0.0455791300d0, &
       &                                                                   +0.0397326900d0, +0.021464950d0, -0.005367417d0, -0.0416953600d0, &
       &                                                                   -0.0821755700d0, -0.118977800d0, -0.142263700d0, -0.1450668000d0, &
       &                                                                   -0.1261692000d0, -0.092231910d0, -0.053753110d0, -0.0204291900d0, &
       &                                                                   +0.0030929600d0, +0.017133450d0, +0.024447280d0, +0.0276528200d0, &
       &                                                                   +0.0286676700d0, +0.028224860d0, +0.027108730d0, +0.0254895900d0, &
       &                                                                   +0.0237081900d0, +0.021728380d0, +0.019877320d0, +0.0179378300d0, &
       &                                                                   +0.0162210100d0, +0.014528280d0, +0.013191270d0, +0.0115311200d0, &
       &                                                                   +0.0000000000d0                                                   &
       &                                                                  ]                                                                , &
       &                                           fitCoefficient3       =[                                                                  &
       &                                                                   -0.0018549720d0, +0.005564917d0, +0.004970250d0, +0.0064291400d0, &
       &                                                                   +0.0075631430d0, +0.008818313d0, +0.009538678d0, +0.0097768030d0, &
       &                                                                   +0.0077292820d0, +0.006305929d0, -0.006327923d0, -0.0097440620d0, &
       &                                                                   -0.0304462400d0, -0.044720620d0, -0.060546570d0, -0.0674670500d0, &
       &                                                                   -0.0613371100d0, -0.038809690d0, -0.004671889d0, +0.0314960000d0, &
       &                                                                   +0.0565621500d0, +0.064131340d0, +0.055539840d0, +0.0392035900d0, &
       &                                                                   +0.0234008200d0, +0.012189720d0, +0.005342576d0, +0.0016914110d0, &
       &                                                                   -0.0007380262d0, -0.001860206d0, -0.002698577d0, -0.0029689890d0, &
       &                                                                   -0.0032996870d0, -0.003085106d0, -0.003232478d0, -0.0028613650d0, &
       &                                                                   -0.0028212210d0, -0.002228346d0, -0.002766919d0, +0.0009223057d0, &
       &                                                                   +0.0000000000d0                                                   &
       &                                                                  ]

contains

  function sutherland1998ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{gauntFactorSutherland1998} gaunt factor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (gauntFactorSutherland1998     )                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(atomicIonizationPotentialClass), pointer       :: atomicIonizationPotential_

    !![
    <objectBuilder class="atomicIonizationPotential" name="atomicIonizationPotential_" source="parameters"/>
    !!]
    self=gauntFactorSutherland1998(atomicIonizationPotential_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="atomicIonizationPotential_"/>
    !!]
    return
  end function sutherland1998ConstructorParameters

  function sutherland1998ConstructorInternal(atomicIonizationPotential_) result(self)
    !!{
    Internal constructor for the \refClass{gauntFactorSutherland1998} gaunt factor class.
    !!}
    implicit none
    type (gauntFactorSutherland1998     )                        :: self
    class(atomicIonizationPotentialClass), intent(in   ), target :: atomicIonizationPotential_
    !![
    <constructorAssign variables="*atomicIonizationPotential_"/>
    !!]

    return
  end function sutherland1998ConstructorInternal

  subroutine sutherland1998Destructor(self)
    !!{
    Destructor for the \refClass{gauntFactorSutherland1998} gaunt factor class.
    !!}
    implicit none
    type(gauntFactorSutherland1998), intent(inout) :: self

    !![
    <objectDestructor name="self%atomicIonizationPotential_"/>
    !!]
    return
  end subroutine sutherland1998Destructor

  double precision function sutherland1998Total(self,atomicNumber,electronNumber,temperature)
    !!{
    Compute thermally averaged Gaunt factors for thermal electron distributions using the tabulations and fits of
    \cite{sutherland_accurate_1998}.
    !!}
    use            :: Arrays_Search               , only : searchArray
    use            :: Error                       , only : Error_Report
    use, intrinsic :: ISO_C_Binding               , only : c_size_t
    use            :: Numerical_Constants_Physical, only : boltzmannsConstant
    use            :: Numerical_Constants_Units   , only : electronVolt
    implicit none
    class           (gauntFactorSutherland1998), intent(inout) :: self
    integer                                    , intent(in   ) :: atomicNumber     , electronNumber
    double precision                           , intent(in   ) :: temperature
    integer         (c_size_t                 )                :: iTable
    double precision                                           :: energyScaleOffset, energyScaleLogarithmic

    ! Return zero for unphysical temperatures
    if (temperature <= 0.0d0) then
       sutherland1998Total=0.0d0
       return
    end if
    ! Validate input.
    if (electronNumber > atomicNumber) call Error_Report('number of electrons exceeds atomic number'//{introspection:location})
    ! Return zero if ionization potential is not available for this ion.
    if (self%atomicIonizationPotential_%potential(atomicNumber,electronNumber) == 0.0d0) then
       sutherland1998Total=0.0d0
    else
       energyScaleLogarithmic=log10(self%atomicIonizationPotential_%potential(atomicNumber,electronNumber)*electronVolt/temperature)-log10(boltzmannsConstant)
       iTable=searchArray(energyScalesLogarithmic,energyScaleLogarithmic)
       if (iTable <=                0) iTable=               1
       if (iTable >= countCoefficient) iTable=countCoefficient
       ! Interpolate to get actual value.
       energyScaleOffset          =+energyScaleLogarithmic          &
            &                      -energyScalesLogarithmic(iTable)
       sutherland1998Total=+gauntFactorZeroPoint           (iTable) &
            &                      +    energyScaleOffset           &
            &                      *(                               &
            &                        +  fitCoefficient1    (iTable) &
            &                        +  energyScaleOffset           &
            &                        *(                             &
            &                          +fitCoefficient2    (iTable) &
            &                          +energyScaleOffset           &
            &                          *fitCoefficient3    (iTable) &
            &                         )                             &
            &                       )
    end if
    return
  end function sutherland1998Total
