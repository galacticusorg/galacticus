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
  An implementation of dark matter halo profile concentrations using the \cite{prada_halo_2011} algorithm.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMONFW

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationPrada2011">
   <description>
    A dark matter profile concentration class in which the concentration is computed using a fitting function from
    \cite{prada_halo_2011}:
    \begin{equation}
    c(M,t) = B_0(x) \mathcal{C}(\sigma^\prime),
    \end{equation}
    where
    \begin{eqnarray}
    \sigma^\prime(M,t) &amp;=&amp; B_1(x) \sigma(M,t), \\
    B_0(x) &amp;=&amp; c_\mathrm{min}(x)/c_\mathrm{min}(1.393), \\
    B_1(x) &amp;=&amp; \sigma^{-1}_\mathrm{min}(x)/\sigma^{-1}_\mathrm{min}(1.393), \\
    c_\mathrm{min}(x) &amp;=&amp; c_0 + (c_1-c_0) [\tan^{-1}\{\alpha (x-x_0)\}/\Pi+1/2], \\
    \sigma^{-1}_\mathrm{min}(x) &amp;=&amp; \sigma^{-1}_0 + (\sigma^{-1}_1-\sigma^{-1}_0) [\tan^{-1}\{\beta(x-x_1)\}/\Pi+1/2], \\
    \mathcal{C}(\sigma^\prime) &amp;=&amp; A [(\sigma^\prime)/b)^c+1] \exp(d/\sigma^{\prime 2}), \\
    x &amp;=&amp; (\Omega_\Lambda/\Omega_\mathrm{M})^{1/3} a(t),
    \end{eqnarray}
    with the following parameters (default values taken from \cite{prada_halo_2011} given in []): $A=${\normalfont \ttfamily
    [A]}$=2.881$, $b=${\normalfont \ttfamily [B]}$=1.257$, $c=${\normalfont \ttfamily [C]}$=1.022$, $d=${\normalfont \ttfamily
    [D]}$=0.060$, $c_0=${\normalfont \ttfamily [C0]}$=3.681$, $c_1=${\normalfont \ttfamily [C1]}$=5.033$, $x_0=${\normalfont
    \ttfamily [X0]}$=0.424$, $x_1=${\normalfont \ttfamily [X1]}$=0.526$, $\sigma^{-1}_0=${\normalfont \ttfamily
    [sigma0]}$=1.047$, $\sigma^{-1}_1=${\normalfont \ttfamily [sigma1]}$=1.646$, $\alpha=${\normalfont \ttfamily
    [alpha]}$=6.948$, and $\beta=${\normalfont \ttfamily [beta]}$=7.386$.
   </description>
   <deepCopy>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </stateStorable>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationPrada2011
     !!{
     A dark matter halo profile concentration class implementing the algorithm of \cite{prada_halo_2011}.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_              => null()
     class           (cosmologyParametersClass     ), pointer :: cosmologyParameters_             => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_        => null()
     type            (virialDensityContrastFixed   ), pointer :: virialDensityContrastDefinition_ => null()
     type            (darkMatterProfileDMONFW      ), pointer :: darkMatterProfileDMODefinition_  => null()
     double precision                                         :: A                                         , B            , &
          &                                                      C                                         , D            , &
          &                                                      C0                                        , C1           , &
          &                                                      X0                                        , X1           , &
          &                                                      inverseSigma0                             , inverseSigma1, &
          &                                                      alpha                                     , beta
   contains
     final     ::                                   prada2011Destructor
     procedure :: concentration                  => prada2011Concentration
     procedure :: densityContrastDefinition      => prada2011DensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => prada2011DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationPrada2011

  interface darkMatterProfileConcentrationPrada2011
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationPrada2011} dark matter halo profile concentration class.
     !!}
     module procedure prada2011ConstructorParameters
     module procedure prada2011ConstructorInternal
  end interface darkMatterProfileConcentrationPrada2011

contains

  function prada2011ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily prada2011} dark matter halo profile concentration class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileConcentrationPrada2011)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass               ), pointer       :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass          ), pointer       :: cosmologicalMassVariance_
    double precision                                                         :: A                        , B            , &
         &                                                                      C                        , D            , &
         &                                                                      C0                       , C1           , &
         &                                                                      X0                       , inverseSigma0, &
         &                                                                      inverseSigma1            , alpha        , &
         &                                                                      beta                     , X1

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>A</name>
      <source>parameters</source>
      <variable>A</variable>
      <defaultValue>2.881d0</defaultValue>
      <defaultSource>\cite{prada_halo_2011}</defaultSource>
      <description>The parameter $A$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    </inputParameter>
    <inputParameter>
      <name>B</name>
      <source>parameters</source>
      <variable>B</variable>
      <defaultValue>1.257d0</defaultValue>
      <defaultSource>\cite{prada_halo_2011}</defaultSource>
      <description>The parameter $b$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    </inputParameter>
    <inputParameter>
      <name>C</name>
      <source>parameters</source>
      <variable>C</variable>
      <defaultValue>1.022d0</defaultValue>
      <defaultSource>\cite{prada_halo_2011}</defaultSource>
      <description>The parameter $c$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    </inputParameter>
    <inputParameter>
      <name>D</name>
      <source>parameters</source>
      <variable>D</variable>
      <defaultValue>0.060d0</defaultValue>
      <defaultSource>\cite{prada_halo_2011}</defaultSource>
      <description>The parameter $d$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    </inputParameter>
    <inputParameter>
      <name>C0</name>
      <source>parameters</source>
      <variable>C0</variable>
      <defaultValue>3.681d0</defaultValue>
      <defaultSource>\cite{prada_halo_2011}</defaultSource>
      <description>The parameter $c_0$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    </inputParameter>
    <inputParameter>
      <name>C1</name>
      <source>parameters</source>
      <variable>C1</variable>
      <defaultValue>5.033d0</defaultValue>
      <defaultSource>\cite{prada_halo_2011}</defaultSource>
      <description>The parameter $c_1$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    </inputParameter>
    <inputParameter>
      <name>X0</name>
      <source>parameters</source>
      <variable>X0</variable>
      <defaultValue>0.424d0</defaultValue>
      <defaultSource>\cite{prada_halo_2011}</defaultSource>
      <description>The parameter $x_0$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    </inputParameter>
    <inputParameter>
      <name>X1</name>
      <source>parameters</source>
      <variable>X1</variable>
      <defaultValue>0.526d0</defaultValue>
      <defaultSource>\cite{prada_halo_2011}</defaultSource>
      <description>The parameter $x_1$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    </inputParameter>
    <inputParameter>
      <name>inverseSigma0</name>
      <source>parameters</source>
      <variable>inverseSigma0</variable>
      <defaultValue>1.047d0</defaultValue>
      <defaultSource>\cite{prada_halo_2011}</defaultSource>
      <description>The parameter $\sigma^{-1}_0$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    </inputParameter>
    <inputParameter>
      <name>inverseSigma1</name>
      <source>parameters</source>
      <variable>inverseSigma1</variable>
      <defaultValue>1.646d0</defaultValue>
      <defaultSource>\cite{prada_halo_2011}</defaultSource>
      <description>The parameter $\sigma^{-1}_1$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    </inputParameter>
    <inputParameter>
      <name>alpha</name>
      <source>parameters</source>
      <variable>alpha</variable>
      <defaultValue>6.948d0</defaultValue>
      <defaultSource>\cite{prada_halo_2011}</defaultSource>
      <description>The parameter $\alpha$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <variable>beta</variable>
      <defaultValue>7.386d0</defaultValue>
      <defaultSource>\cite{prada_halo_2011}</defaultSource>
      <description>The parameter $\beta$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"      source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=darkMatterProfileConcentrationPrada2011(A,B,C,D,C0,C1,X0,X1,inverseSigma0,inverseSigma1,alpha,beta,cosmologyFunctions_,cosmologyParameters_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function prada2011ConstructorParameters

  function prada2011ConstructorInternal(A,B,C,D,C0,C1,X0,X1,inverseSigma0,inverseSigma1,alpha,beta,cosmologyFunctions_,cosmologyParameters_,cosmologicalMassVariance_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileConcentrationPrada2011} dark matter halo profile concentration class.
    !!}
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleVirialDensityContrastDefinition
    use :: Virial_Density_Contrast, only : fixedDensityTypeCritical
    implicit none
    type            (darkMatterProfileConcentrationPrada2011           )                        :: self
    double precision                                                            , intent(in   ) :: A                        , B            , &
         &                                                                                         C                        , D            , &
         &                                                                                         C0                       , C1           , &
         &                                                                                         X0                       , inverseSigma0, &
         &                                                                                         inverseSigma1            , alpha        , &
         &                                                                                         beta                     , X1
    class           (cosmologyFunctionsClass                           ), target, intent(in   ) :: cosmologyFunctions_
    class           (cosmologyParametersClass                          ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass                     ), target, intent(in   ) :: cosmologicalMassVariance_
    type            (darkMatterHaloScaleVirialDensityContrastDefinition), pointer       :: darkMatterHaloScaleDefinition_
    !![
    <constructorAssign variables="A,B,C,D,C0,C1,X0,X1,inverseSigma0,inverseSigma1,alpha,beta,*cosmologyFunctions_,*cosmologyParameters_,*cosmologicalMassVariance_"/>
    !!]

    allocate(self%virialDensityContrastDefinition_)
    allocate(self%darkMatterProfileDMODefinition_ )
    allocate(     darkMatterHaloScaleDefinition_  )
    !![
    <referenceConstruct owner="self" object="virialDensityContrastDefinition_">
     <constructor>
      virialDensityContrastFixed                        (                                                                            &amp;
       &amp;                                             densityContrastValue                =200.0d0                              , &amp;
       &amp;                                             densityType                         =fixedDensityTypeCritical             , &amp;
       &amp;                                             turnAroundOverVirialRadius          =2.0d0                                , &amp;
       &amp;                                             cosmologyParameters_                =self%cosmologyParameters_            , &amp;
       &amp;                                             cosmologyFunctions_                 =self%cosmologyFunctions_               &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <referenceConstruct              object="darkMatterHaloScaleDefinition_"  >
     <constructor>
      darkMatterHaloScaleVirialDensityContrastDefinition(                                                                            &amp;
       &amp;                                             cosmologyParameters_                =self%cosmologyParameters_            , &amp;
       &amp;                                             cosmologyFunctions_                 =self%cosmologyFunctions_             , &amp;
       &amp;                                             virialDensityContrast_              =self%virialDensityContrastDefinition_  &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <referenceConstruct owner="self" object="darkMatterProfileDMODefinition_" >
     <constructor>
      darkMatterProfileDMONFW                           (                                                                            &amp;
       &amp;                                             velocityDispersionUseSeriesExpansion=.true.                               , &amp;
       &amp;                                             darkMatterHaloScale_                =darkMatterHaloScaleDefinition_         &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <objectDestructor                name  ="darkMatterHaloScaleDefinition_" />
    !!]
    return
  end function prada2011ConstructorInternal

  subroutine prada2011Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationPrada2011} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationPrada2011), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%cosmologicalMassVariance_"       />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    <objectDestructor name="self%darkMatterProfileDMODefinition_" />
    !!]
    return
  end subroutine prada2011Destructor

  double precision function prada2011Concentration(self,node)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the \cite{prada_halo_2011} algorithm.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout), target  :: self
    type            (treeNode                               ), intent(inout), target  :: node
    class           (nodeComponentBasic                     )               , pointer :: basic
    double precision                                                                  :: massNode, sigmaPrime, &
         &                                                                               timeNode, x

    ! Compute concentration.
    basic                 => node %basic()
    massNode              =  basic%mass ()
    timeNode              =  basic%time ()
    x                     =  (self%cosmologyParameters_%OmegaDarkEnergy()/self%cosmologyParameters_%OmegaMatter())**(1.0d0/3.0d0)*self%cosmologyFunctions_%expansionFactor(timeNode)
    sigmaPrime            =  prada2011B1(self,x)*self%cosmologicalMassVariance_%rootVariance(massNode,timeNode)
    prada2011Concentration=  prada2011B0(self,x)*prada2011C(self,sigmaPrime)
    return
  end function prada2011Concentration

  double precision function prada2011B0(self,x)
    !!{
    The function $B_0(x)$ as defined in eqn.~(18) of \cite{prada_halo_2011}.
    !!}
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout) :: self
    double precision                                         , intent(in   ) :: x

    prada2011B0=prada2011cMin(self,x)/prada2011cMin(self,1.393d0)
  end function prada2011B0

  double precision function prada2011B1(self,x)
    !!{
    The function $B_1(x)$ as defined in eqn.~(18) of \cite{prada_halo_2011}.
    !!}
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout) :: self
    double precision                                         , intent(in   ) :: x

    prada2011B1=prada2011inverseSigmaMin(self,x)/prada2011inverseSigmaMin(self,1.393d0)
  end function prada2011B1

  double precision function prada2011cMin(self,x)
    !!{
    The function $c_\mathrm{min}(x)$ as defined in eqn.~(19) of \cite{prada_halo_2011}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout) :: self
    double precision                                         , intent(in   ) :: x

    prada2011cMin=self%C0+(self%C1-self%C0)*(atan(self%alpha*(x-self%X0))/Pi+0.5d0)
  end function prada2011cMin

  double precision function prada2011inverseSigmaMin(self,x)
    !!{
    The function $\sigma^{-1}_\mathrm{min}(x)$ as defined in eqn.~(20) of \cite{prada_halo_2011}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout) :: self
    double precision                                         , intent(in   ) :: x

    prada2011inverseSigmaMin=self%inverseSigma0+(self%inverseSigma1-self%inverseSigma0)*(atan(self%beta*(x-self%X1))/Pi+0.5d0)
  end function prada2011inverseSigmaMin

  double precision function prada2011C(self,sigmaPrime)
    !!{
    The function $\mathcal{C}(\sigma^\prime)$ as defined in eqn.~(17) of \cite{prada_halo_2011}.
    !!}
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout) :: self
    double precision                                         , intent(in   ) :: sigmaPrime

    prada2011C=self%A*((sigmaPrime/self%B)**self%C+1.0d0)*exp(self%D/sigmaPrime**2)
  end function prada2011C

  function prada2011DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of concentration in the \cite{prada_halo_2011} algorithm.
    !!}
    implicit none
    class(virialDensityContrastClass             ), pointer       :: prada2011DensityContrastDefinition
    class(darkMatterProfileConcentrationPrada2011), intent(inout) :: self

    prada2011DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function prada2011DensityContrastDefinition

  function prada2011DarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of concentration in the
    \cite{prada_halo_2011} algorithm.
    !!}
    implicit none
    class(darkMatterProfileDMOClass              ), pointer       :: prada2011DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationPrada2011), intent(inout) :: self

    prada2011DarkMatterProfileDefinition => self%darkMatterProfileDMODefinition_
    return
  end function prada2011DarkMatterProfileDefinition
