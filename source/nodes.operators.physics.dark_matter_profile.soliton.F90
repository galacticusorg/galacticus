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

  !+    Contributions to this file made by: Yu Zhao

  !!{
     Implements a node operator class that evaluates the \gls{fdm} solitonic core–halo relation, following Equation (15) of \cite{chan_diversity_2022}, with modifications to include time differentiation and integration to track the evolution of the core mass.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Dark_Matter_Particles  , only : darkMatterParticleClass
  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Virial_Density_Contrast, only : virialDensityContrastClass

  !![
  <nodeOperator name="nodeOperatorDarkMatterProfileSoliton">
   <description>
     A node operator class that evaluates the \gls{fdm} solitonic core–halo relation,
     following Equation (15) of \cite{chan_diversity_2022}, with modifications to
     include time differentiation and integration to track the evolution of the core mass.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDarkMatterProfileSoliton
     !!{
     A node operator class implementing the time-evolved solitonic core–halo relation following \cite{chan_diversity_2022}.
     !!}
     private
     class           (darkMatterHaloScaleClass  ), pointer :: darkMatterHaloScale_  => null()
     class           (darkMatterParticleClass   ), pointer :: darkMatterParticle_   => null()
     class           (cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_   => null()
     class           (cosmologyParametersClass  ), pointer :: cosmologyParameters_  => null()
     class           (virialDensityContrastClass), pointer :: virialDensityContrast_=> null()
     integer                                               :: densityCoreAccretionID
     double precision                                      :: massParticle
   contains
     final     ::                                darkMatterProfileSolitonDestructor
     procedure :: nodeTreeInitialize          => darkMatterProfileSolitonNodeTreeInitialize
     procedure :: differentialEvolution       => darkMatterProfileSolitonDifferentialEvolution
     procedure :: differentialEvolutionScales => darkMatterProfileSolitonDifferentialEvolutionScales
  end type nodeOperatorDarkMatterProfileSoliton
  
  interface nodeOperatorDarkMatterProfileSoliton
     !!{
     Constructors for the \refClass{nodeOperatorDarkMatterProfileSoliton} node operator class.
     !!}
     module procedure darkMatterProfileSolitonConstructorParameters
     module procedure darkMatterProfileSolitonConstructorInternal
  end interface nodeOperatorDarkMatterProfileSoliton
  
contains

  function darkMatterProfileSolitonConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDarkMatterProfileSoliton} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorDarkMatterProfileSoliton )                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass             ), pointer       :: darkMatterHaloScale_
    class           (darkMatterParticleClass              ), pointer       :: darkMatterParticle_
    class           (cosmologyFunctionsClass              ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass             ), pointer       :: cosmologyParameters_
    class           (virialDensityContrastClass           ), pointer       :: virialDensityContrast_

    !![
    <objectBuilder class="darkMatterHaloScale"   name="darkMatterHaloScale_"   source="parameters"/>
    <objectBuilder class="darkMatterParticle"    name="darkMatterParticle_"    source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
     self=nodeOperatorDarkMatterProfileSoliton(darkMatterHaloScale_,darkMatterParticle_,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_)
    !![
    <inputParametersValidate source="parameters"            />
    <objectDestructor        name  ="darkMatterHaloScale_"  />
    <objectDestructor        name  ="darkMatterParticle_"   />
    <objectDestructor        name  ="cosmologyFunctions_"   />
    <objectDestructor        name  ="cosmologyParameters_"  />
    <objectDestructor        name  ="virialDensityContrast_"/>
    !!]
    return
  end function darkMatterProfileSolitonConstructorParameters

  function darkMatterProfileSolitonConstructorInternal(darkMatterHaloScale_,darkMatterParticle_,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorDarkMatterProfileSoliton} node operator class.
    !!}
    use :: Dark_Matter_Particles           , only : darkMatterParticleFuzzyDarkMatter
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    type            (nodeOperatorDarkMatterProfileSoliton)                     :: self
    class           (darkMatterHaloScaleClass            ), intent(in), target :: darkMatterHaloScale_
    class           (darkMatterParticleClass             ), intent(in), target :: darkMatterParticle_
    class           (cosmologyFunctionsClass             ), intent(in), target :: cosmologyFunctions_
    class           (cosmologyParametersClass            ), intent(in), target :: cosmologyParameters_
    class           (virialDensityContrastClass          ), intent(in), target :: virialDensityContrast_
    !![
    <constructorAssign variables="*darkMatterHaloScale_,*darkMatterParticle_,*cosmologyFunctions_,*cosmologyParameters_,*virialDensityContrast_"/>
    <addMetaProperty component="basic" name="densityCoreAccretion" id="self%densityCoreAccretionID" isEvolvable="yes" isCreator="yes"/>
    !!]

    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       self%massParticle=+darkMatterParticle__%mass()*kilo
    class default
       call Error_Report('expected a `darkMatterParticleFuzzyDarkMatter` dark matter particle object'//{introspection:location})
    end select
    return
  end function darkMatterProfileSolitonConstructorInternal

  subroutine darkMatterProfileSolitonDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorDarkMatterProfileSoliton} node operator class.
    !!}
    implicit none
    type            (nodeOperatorDarkMatterProfileSoliton), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"  />
    <objectDestructor name="self%darkMatterParticle_"   />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]

    return
  end subroutine darkMatterProfileSolitonDestructor

  subroutine darkMatterProfileSolitonDifferentialEvolutionScales(self,node)
    !!{
    Set absolute ODE solver scale for the energy radiated from the hot halo due to cooling following the model of \cite{benson_galaxy_2010-1}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic
    implicit none
    class           (nodeOperatorDarkMatterProfileSoliton), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , parameter     :: scaleRelative     =+1.0d-1
    class           (nodeComponentBasic                  ), pointer       :: basic

    basic=>node%basic()
    call basic%floatRank0MetaPropertyScale(self%densityCoreAccretionID, scaleRelative)

    return
  end subroutine darkMatterProfileSolitonDifferentialEvolutionScales

  subroutine darkMatterProfileSolitonDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Accumulates an estimate of the energy radiated from the hot halo due to cooling following the model of \cite{benson_galaxy_2010-1}.
    !!}
    use :: Galacticus_Nodes                , only : treeNode           , nodeComponentBasic       , nodeComponentDarkMatterProfile
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Units       , only : electronVolt
    use :: Numerical_Constants_Physical    , only : plancksConstant
    use :: Cosmology_Parameters            , only : hubbleUnitsStandard, hubbleUnitsLittleH
    implicit none
    class           (nodeOperatorDarkMatterProfileSoliton), intent   (inout), target  :: self
    type            (treeNode                            ), intent   (inout), target  :: node
    class           (nodeComponentBasic                  ), pointer                   :: basic
    class           (nodeComponentDarkMatterProfile      ), pointer                   :: darkMatterProfile
    logical                                               , intent   (inout)          :: interrupt
    procedure       (interruptTask                       ), intent   (inout), pointer :: functionInterrupt
    integer                                               , intent   (in   )          :: propertyType
    double precision                                      , parameter                 :: plancksConstantBar=+plancksConstant                                & ! ℏ in units of eV s.
         &                                                                                                  /2.0d0                                          &
         &                                                                                                  /Pi                                             &
         &                                                                                                  /electronVolt
    double precision                                                                  :: massHalo                            , expansionFactor            , &
         &                                                                               redshift                            , concentration              , &
         &                                                                               hubbleConstant                      , hubbleConstantLittle       , &
         &                                                                               OmegaMatter                         , densityMatter              , &
         &                                                                               zeta_0                              , zeta_z                     , &
         &                                                                               radiusVirial                        , radiusScale                , &
         &                                                                               massCoreRate                        , densityScale               , &
         &                                                                               expansionRate                       , zetaRate
    double precision                                      , parameter                 :: alpha             =0.515            , beta               =8.0d6  , &
         &                                                                               gamma             =10.0d0**(-5.73d0)                                 ! Best-fitting parameters from Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
    double precision                                                                  :: A                                   , K0                         , & ! first term of Equation 15 from Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
         &                                                                               K1                                                                   ! second term of Equation 15 from Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    ! Get required components.
    basic             => node%basic            ()
    darkMatterProfile => node%darkMatterProfile()
    ! Extract basic properties of the node.
    expansionFactor=+self             %cosmologyFunctions_ %expansionFactor            (basic%time           ())
    expansionRate  =+self             %cosmologyFunctions_ %expansionRate              (expansionFactor        )
    redshift       =+self             %cosmologyFunctions_ %redshiftFromExpansionFactor(expansionFactor        )
    massHalo       =+basic                                 %mass                       (                       )
    radiusScale    =+darkMatterProfile                     %scale                      (                       )
    radiusVirial   =+self             %darkMatterHaloScale_%radiusVirial               (node                   )
    concentration  =+                                       radiusVirial                                         &
         &          /                                       radiusScale
    densityScale   =+1.0d0                                    &
         &          /4.0d0                                    &
         &          /Pi                                       &
         &          *massHalo                                 &
         &          /radiusScale**3                           &
         &          /(                                        &
         &            +              log(1.0d0+concentration) &
         &            -concentration/   (1.0d0+concentration) &
         &          )
    ! Extract cosmological parameters for later use.
    hubbleConstant      =+self%cosmologyParameters_%HubbleConstant (hubbleUnitsStandard)
    hubbleConstantLittle=+self%cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH )
    OmegaMatter         =+self%cosmologyParameters_%OmegaMatter    (                   )
    densityMatter       =+self%cosmologyParameters_%densityCritical(                   ) &
         &                      *                   OmegaMatter
    zeta_0              =+self%virialDensityContrast_%densityContrast            (massHalo,expansionFactor=1.0d0          )
    zeta_z              =+self%virialDensityContrast_%densityContrast            (massHalo,expansionFactor=expansionFactor)
    zetaRate            =+self%virialDensityContrast_%densityContrastRateofChange(massHalo,expansionFactor=expansionFactor)

    A                   =+self%massParticle         &
         &               /8.0d-23
    K0                  =+beta                      &
         &               *A**(-1.5d0)
    K1                  =(sqrt(                     &
         &                     +zeta_z              &
         &                     /zeta_0              &
         &                    )                     &
         &                *massHalo                 &
         &                /gamma                    &
         &               )**alpha                   &
         &               *A**(1.5d0*(alpha-1.0d0))
    massCoreRate        =-0.5d0                     &
         &               *expansionFactor**(-0.5d0) &
         &               *expansionRate             &
         &               *(K0+K1)                   &
         &               +alpha                     &
         &               *expansionFactor**(-0.5d0) &
         &               *K1                        &
         &               *(+0.5d0                   &
         &                 *zetaRate                &
         &                 /zeta_z                  &
         &                 +basic%accretionRate()   &
         &                 /massHalo                &
         &                )
    call basic%floatRank0MetaPropertyRate(                             &
         &                                self%densityCoreAccretionID, &
         &                                massCoreRate                 &
         &                               )
    return
  end subroutine darkMatterProfileSolitonDifferentialEvolution

  subroutine darkMatterProfileSolitonNodeTreeInitialize(self,node)
    use :: Galacticus_Nodes                , only : treeNode           , nodeComponentBasic       , nodeComponentDarkMatterProfile
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Units       , only : electronVolt
    use :: Numerical_Constants_Physical    , only : plancksConstant
    use :: Cosmology_Parameters            , only : hubbleUnitsStandard, hubbleUnitsLittleH
    implicit none
    class           (nodeOperatorDarkMatterProfileSoliton), intent(inout), target :: self
    type            (treeNode                            ), intent(inout), target :: node
    class           (nodeComponentBasic                  ), pointer               :: basic
    class           (nodeComponentDarkMatterProfile      ), pointer               :: darkMatterProfile
    double precision                                      , parameter             :: plancksConstantBar=+plancksConstant                                & ! ℏ in units of eV s.
         &                                                                                              /2.0d0                                          &
         &                                                                                              /Pi                                             &
         &                                                                                              /electronVolt
    double precision                                                              :: massHalo                            , expansionFactor            , &
         &                                                                           redshift                            , concentration              , &
         &                                                                           hubbleConstant                      , hubbleConstantLittle       , &
         &                                                                           OmegaMatter                         , densityMatter              , &
         &                                                                           zeta_0                              , zeta_z                     , &
         &                                                                           radiusVirial                        , radiusScale                , &
         &                                                                           massCoreNormal                      , densityScale
    double precision                                      , parameter             :: alpha             =0.515            , beta               =8.0d6  , &
         &                                                                           gamma             =10.0d0**(-5.73d0)                                 ! Best-fitting parameters from Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).

    ! Get required components.
    basic             => node%basic            ()
    darkMatterProfile => node%darkMatterProfile()
    ! Extract basic properties of the node.
    expansionFactor=+self             %cosmologyFunctions_% expansionFactor            (basic%time           ())
    redshift       =+self             %cosmologyFunctions_ %redshiftFromExpansionFactor(      expansionFactor  )
    massHalo       =+basic                                 %mass                       (                       )
    radiusScale    =+darkMatterProfile                     %scale                      (                       )
    radiusVirial   =+self             %darkMatterHaloScale_%radiusVirial               (node                   )
    concentration  =+                                       radiusVirial                                         &
         &          /                                       radiusScale
    densityScale   =+1.0d0                                    &
         &          /4.0d0                                    &
         &          /Pi                                       &
         &          *massHalo                                 &
         &          /radiusScale**3                           &
         &          /(                                        &
         &            +              log(1.0d0+concentration) &
         &            -concentration/   (1.0d0+concentration) &
         &          )
    ! Extract cosmological parameters for later use.
    hubbleConstant      =+self%cosmologyParameters_%HubbleConstant (hubbleUnitsStandard)
    hubbleConstantLittle=+self%cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH )
    OmegaMatter         =+self%cosmologyParameters_%OmegaMatter    (                   )
    densityMatter       =+self%cosmologyParameters_%densityCritical(                   ) &
         &                      *                   OmegaMatter
    zeta_0             =+self%virialDensityContrast_%densityContrast(massHalo,expansionFactor=1.0d0          )
    zeta_z             =+self%virialDensityContrast_%densityContrast(massHalo,expansionFactor=expansionFactor)

    massCoreNormal     =+(                          & ! Equation (15) of Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
         &                +beta                     &
         &                *(                        &
         &                  +self%massParticle      &
         &                  /8.0d-23                &
         &                 )**(-1.5d0)              &
         &                +(                        &
         &                  +sqrt(                  &
         &                        +zeta_z           &
         &                        /zeta_0           &
         &                       )                  &
         &                  *massHalo               &
         &                  /gamma                  &
         &                 )**alpha                 &
         &                *(                        &
         &                  +self%massParticle      &
         &                  /8.0d-23                &
         &                 )**(1.5d0*(alpha-1.0d0)) &
         &               )&
         &              /sqrt(expansionFactor)
    call basic%floatRank0MetaPropertySet(                             &
         &                               self%densityCoreAccretionID, &
         &                               massCoreNormal               &
         &                              )
    return
  end subroutine darkMatterProfileSolitonNodeTreeInitialize
