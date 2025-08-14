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
  Implements a dark matter profile scale radius class that uses simple power-law scalings.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !![
  <darkMatterProfileScaleRadius name="darkMatterProfileScaleRadiusPowerLaw">
    <description>
      A dark matter profile scale radius class that uses simple power-law scalings. Specifically, the scale radius is given by:
      \begin{equation}
      r_\mathrm{s} = r(\nu) \left(\frac{M}{M_0}\right)^{\alpha(\nu)} (1+z)^{-\beta(\nu)}
      \end{equation}      
      where $r(\nu)$, $\alpha(\nu)$, and $\beta(\nu)$ are sigmoid functions of the peak height, $\nu$, of the form:
      \begin{equation}
       y(x) = y_0+(y_1-y_0)/(1+\exp[-(x-x_\nu)/\Delta x]),
      \end{equation}
      where $r_0=${\normalfont \ttfamily [radiusLow]}, $r_1=${\normalfont \ttfamily [radiusHigh]}, $r_\nu=${\normalfont \ttfamily
      [radiusTransition]}, $\Delta r=${\normalfont \ttfamily [radiusWidth]}, $\alpha_0=${\normalfont \ttfamily [massLow]},
      $\alpha_1=${\normalfont \ttfamily [massHigh]}, $\alpha_\nu=${\normalfont \ttfamily [massTransition]}, $\Delta
      \alpha=${\normalfont \ttfamily [massWidth]}, $\beta_0=${\normalfont \ttfamily [expansionFactorLow]}, $\beta_1=${\normalfont
      \ttfamily [expansionFactorHigh]}, $\beta_\nu=${\normalfont \ttfamily [expansionFactorTransition]}, and $\Delta
      \beta=${\normalfont \ttfamily [expansionFactorWidth]} , plus a random log-normal scatter of {\normalfont \ttfamily
      [scatter]}~dex.
    </description>
  </darkMatterProfileScaleRadius>
  !!]
  type, extends(darkMatterProfileScaleRadiusClass) :: darkMatterProfileScaleRadiusPowerLaw
     !!{
     A dark matter profile scale radius class that assigns dark matter profile scale radii using the energy-based model of
     \cite{johnson_random_2021}.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_           => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_          => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_     => null()
     double precision                                         :: radiusLow                              , radiusHigh          , &
          &                                                      radiusTransition                       , radiusWidth         , &
          &                                                      massLow                                , massHigh            , &
          &                                                      massTransition                         , massWidth           , &
          &                                                      expansionFactorLow                     , expansionFactorHigh , &
          &                                                      expansionFactorTransition              , expansionFactorWidth, &
          &                                                      scatter
   contains
     final     ::           darkMatterProfileScalePowerLawDestructor
     procedure :: radius => darkMatterProfileScalePowerLawRadius
  end type darkMatterProfileScaleRadiusPowerLaw
  
  interface darkMatterProfileScaleRadiusPowerLaw
     !!{
     Constructors for the \refClass{darkMatterProfileScaleRadiusPowerLaw} node operator class.
     !!}
     module procedure darkMatterProfileScalePowerLawConstructorParameters
     module procedure darkMatterProfileScalePowerLawConstructorInternal
  end interface darkMatterProfileScaleRadiusPowerLaw

contains
  
  function darkMatterProfileScalePowerLawConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileScaleRadiusPowerLaw} dark matter profile scale radius class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (darkMatterProfileScaleRadiusPowerLaw)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass             ), pointer       :: cosmologyFunctions_
    class           (criticalOverdensityClass            ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass       ), pointer       :: cosmologicalMassVariance_
    double precision                                                      :: radiusLow                , radiusHigh          , &
          &                                                                  radiusTransition         , radiusWidth         , &
          &                                                                  massLow                  , massHigh            , &
          &                                                                  massTransition           , massWidth           , &
          &                                                                  expansionFactorLow       , expansionFactorHigh , &
          &                                                                  expansionFactorTransition, expansionFactorWidth, &
          &                                                                  scatter
    
    !![
    <inputParameter>
      <name>radiusLow</name>
      <defaultValue>+0.0154d0</defaultValue>
      <source>parameters</source>
      <description>The parameter $r_0$ in the power-law scale radius model.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusHigh</name>
      <defaultValue>+0.0962d0</defaultValue>
      <source>parameters</source>
      <description>The parameter $r_1$ in the power-law scale radius model.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusTransition</name>
      <defaultValue>+1.2137d0</defaultValue>
      <source>parameters</source>
      <description>The parameter $r_\nu$ in the power-law scale radius model.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusWidth</name>
      <defaultValue>+0.5482d0</defaultValue>
      <source>parameters</source>
      <description>The parameter $\Delta r$ in the power-law scale radius model.</description>
    </inputParameter>
    <inputParameter>
      <name>massLow</name>
      <defaultValue>+0.3895d0</defaultValue>
      <source>parameters</source>
      <description>The parameter $\alph_0$ in the power-law scale radius model.</description>
    </inputParameter>
    <inputParameter>
      <name>massHigh</name>
      <defaultValue>+0.2984d0</defaultValue>
      <source>parameters</source>
      <description>The parameter $\alph_1$ in the power-law scale radius model.</description>
    </inputParameter>
    <inputParameter>
      <name>massTransition</name>
      <defaultValue>-0.2583d0</defaultValue>
      <source>parameters</source>
      <description>The parameter $\alph_\nu$ in the power-law scale radius model.</description>
    </inputParameter>
    <inputParameter>
      <name>massWidth</name>
      <defaultValue>+16.6050d0</defaultValue>
      <source>parameters</source>
      <description>The parameter $\Delta \alpha$ in the power-law scale radius model.</description>
    </inputParameter>
    <inputParameter>
      <name>expansionFactorLow</name>
      <defaultValue>-0.6977d0</defaultValue>
      <source>parameters</source>
      <description>The parameter $\beta_0$ in the power-law scale radius model.</description>
    </inputParameter>
    <inputParameter>
      <name>expansionFactorHigh</name>
      <defaultValue>+0.7972d0</defaultValue>
      <source>parameters</source>
      <description>The parameter $\beta_1$ in the power-law scale radius model.</description>
    </inputParameter>
    <inputParameter>
      <name>expansionFactorTransition</name>
      <defaultValue>+0.5395d0</defaultValue>
      <source>parameters</source>
      <description>The parameter $\beta_\nu$ in the power-law scale radius model.</description>
    </inputParameter>
    <inputParameter>
      <name>expansionFactorWidth</name>
      <defaultValue>+0.4282d0</defaultValue>
      <source>parameters</source>
      <description>The parameter $\Delta \beta$ in the power-law scale radius model.</description>
    </inputParameter>
    <inputParameter>
      <name>scatter</name>
      <defaultValue>+0.1513d0</defaultValue>
      <source>parameters</source>
      <description>.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    <objectBuilder class="criticalOverdensity"          name="criticalOverdensity_"          source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"     name="cosmologicalMassVariance_"     source="parameters"/>
    !!]
    self=darkMatterProfileScaleRadiusPowerLaw(radiusLow,radiusHigh,radiusTransition,radiusWidth,massLow,massHigh,massTransition,massWidth,expansionFactorLow,expansionFactorHigh,expansionFactorTransition,expansionFactorWidth,scatter,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function darkMatterProfileScalePowerLawConstructorParameters

  function darkMatterProfileScalePowerLawConstructorInternal(radiusLow,radiusHigh,radiusTransition,radiusWidth,massLow,massHigh,massTransition,massWidth,expansionFactorLow,expansionFactorHigh,expansionFactorTransition,expansionFactorWidth,scatter,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileScaleRadiusPowerLaw} dark matter profile scale radius class.
    !!}
    implicit none
    type            (darkMatterProfileScaleRadiusPowerLaw)                        :: self
    class           (cosmologyFunctionsClass             ), intent(in   ), target :: cosmologyFunctions_
    class           (criticalOverdensityClass            ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass       ), intent(in   ), target :: cosmologicalMassVariance_
    double precision                                      , intent(in   )         :: radiusLow                , radiusHigh          , &
          &                                                                          radiusTransition         , radiusWidth         , &
          &                                                                          massLow                  , massHigh            , &
          &                                                                          massTransition           , massWidth           , &
          &                                                                          expansionFactorLow       , expansionFactorHigh , &
          &                                                                          expansionFactorTransition, expansionFactorWidth, &
          &                                                                          scatter
    !![
    <constructorAssign variables="radiusLow, radiusHigh, radiusTransition, radiusWidth, massLow, massHigh, massTransition, massWidth, expansionFactorLow, expansionFactorHigh, expansionFactorTransition, expansionFactorWidth, scatter, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_"/>
    !!]

    return
  end function darkMatterProfileScalePowerLawConstructorInternal

  subroutine darkMatterProfileScalePowerLawDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileScaleRadiusPowerLaw} dark matter halo profile scale radius class.
    !!}
    implicit none
    type(darkMatterProfileScaleRadiusPowerLaw), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine darkMatterProfileScalePowerLawDestructor

  double precision function darkMatterProfileScalePowerLawRadius(self,node) result(radiusScale)
    !!{
    Initialize dark matter profile scale radii.
    !!}
    implicit none
    class           (darkMatterProfileScaleRadiusPowerLaw), intent(inout) , target  :: self
    type            (treeNode                            ), intent(inout) , target  :: node
    class           (nodeComponentBasic                  )                , pointer :: basic
    double precision                                      , parameter               :: massReference   =1.0d12
    double precision                                                                :: expansionFactor        , peakHeight  , &
         &                                                                             radiusNormalization    , exponentMass, &
         &                                                                             exponentExpansionFactor
    
    basic                   =>  node                          %basic          (                                             )
    peakHeight              =  +self%criticalOverdensity_     %value          (time=basic%time(),mass=basic%mass(),node=node) &
         &                     /self%cosmologicalMassVariance_%rootVariance   (time=basic%time(),mass=basic%mass()          )
    expansionFactor         =  +self%cosmologyFunctions_      %expansionFactor(time=basic%time()                            )
    radiusNormalization     =  +sigmoid(self%radiusLow         ,self%radiusHigh         ,self%radiusTransition         ,self%radiusWidth         ,peakHeight)
    exponentMass            =  +sigmoid(self%massLow           ,self%massHigh           ,self%massTransition           ,self%massWidth           ,peakHeight)
    exponentExpansionFactor =  +sigmoid(self%expansionFactorLow,self%expansionFactorHigh,self%expansionFactorTransition,self%expansionFactorWidth,peakHeight)
    radiusScale             =  + radiusNormalization                                           &
         &                     *(basic%mass         ()/massReference)**exponentMass            &
         &                     * expansionFactor                     **exponentExpansionFactor
    return

  contains

    double precision function sigmoid(y0,y1,xNu,DeltaX,x)
      !!{
      A sigmoid function,
      \begin{equation}
       y(x) = y_0+(y_1-y_0)/(1+\exp[-(x-x_\nu)/\Delta x]).
      \end{equation}
      !!}
      implicit none
      double precision, intent(in   ) :: y0 , y1    , &
           &                             xNu, Deltax, &
           &                             x

      sigmoid=y0+(y1-y0)/(1.0d0+exp(-(x-xNu)/Deltax))
      return
    end function sigmoid
    
  end function darkMatterProfileScalePowerLawRadius
