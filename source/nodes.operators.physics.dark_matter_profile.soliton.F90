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
     class           (darkMatterHaloScaleClass  ), pointer   :: darkMatterHaloScale_  => null()
     class           (darkMatterParticleClass   ), pointer   :: darkMatterParticle_   => null()
     class           (cosmologyFunctionsClass   ), pointer   :: cosmologyFunctions_   => null()
     class           (cosmologyParametersClass  ), pointer   :: cosmologyParameters_  => null()
     class           (virialDensityContrastClass), pointer   :: virialDensityContrast_=> null()
     integer                                                 :: massCoreID                     , massHaloID        , &
         &                                                      zeta0ID                        , zetazID           , &
         &                                                      massParticleID                 , expansionFactorID
     double precision                                        :: massParticle                   , massCoreMinimum
     double precision                                        :: alpha       =0.515             , beta      =8.0d6  , &
         &                                                      gamma       =10.0d0**(-5.73d0)                       ! Best-fitting parameters from Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
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
    <addMetaProperty component="darkMatterProfile" name="massCore"        id="self%massCoreID"        isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="zeta0"           id="self%zeta0ID"           isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="zetaz"           id="self%zetazID"           isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="expansionFactor" id="self%expansionFactorID" isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="massParticle"    id="self%massParticleID"    isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="massHalo"        id="self%massHaloID"        isEvolvable="no"  isCreator="yes"/>

    !!]

    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       self%massParticle=+darkMatterParticle__%mass()*kilo
    class default
       call Error_Report('expected a `darkMatterParticleFuzzyDarkMatter` dark matter particle object'//{introspection:location})
    end select

    ! Minimum soliton core mass in FDM.
    ! Using Mc,min ~ 0.25 × Mmin,halo (soliton-dominated limit), from Schive et al. (2014; PRL; 113; 1302; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S).
    ! The minimum halo-mass scaling is mentioned in the abstract of Hui et al. (2017; PRD; 95; 3541; https://ui.adsabs.harvard.edu/abs/2017PhRvD..95d3541H).
    self%massCoreMinimum  =+0.25d0*1.0d7*(self%massParticle/1.0d-22)**(-1.5d0)

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
    Set the absolute ODE solver scale for the solitonic core mass evolution,
    using a fraction of the minimum core mass as reference, following \cite{chan_diversity_2022}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorDarkMatterProfileSoliton), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    class           (nodeComponentDarkMatterProfile      ), pointer       :: darkMatterProfile
    double precision                                                      :: scaleRelative

    ! Set the absolute tolerance scale for ODE integration to 10% of the minimum solitonic core mass.
    scaleRelative=+0.1d0*self%massCoreMinimum

    darkMatterProfile=>node%darkMatterProfile()
    call darkMatterProfile%floatRank0MetaPropertyScale(self%massCoreID, scaleRelative)

    return
  end subroutine darkMatterProfileSolitonDifferentialEvolutionScales

  subroutine darkMatterProfileSolitonDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Time derivative of the solitonic core mass following \cite{chan_diversity_2022}.
    !!}
    use :: Galacticus_Nodes                , only : treeNode           , nodeComponentBasic       , nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorDarkMatterProfileSoliton), intent   (inout), target  :: self
    type            (treeNode                            ), intent   (inout), target  :: node
    class           (nodeComponentBasic                  ), pointer                   :: basic
    class           (nodeComponentDarkMatterProfile      ), pointer                   :: darkMatterProfile
    logical                                               , intent   (inout)          :: interrupt
    procedure       (interruptTask                       ), intent   (inout), pointer :: functionInterrupt
    integer                                               , intent   (in   )          :: propertyType
    double precision                                                                  :: massHalo                            , expansionFactor            , &
         &                                                                               zeta_0                              , zeta_z                     , &
         &                                                                               massCoreRate                        , zetaRate                   , &
         &                                                                               expansionRate
    double precision                                                                  :: A                                   , K0                         , & ! first term of Equation 15 from Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
         &                                                                               K1                                                                   ! second term of Equation 15 from Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    ! Get required components.
    basic             => node%basic            ()
    darkMatterProfile => node%darkMatterProfile()
    ! Extract basic properties of the node.
    expansionFactor     =+self             %cosmologyFunctions_ %expansionFactor            (basic%time           ())
    expansionRate       =+self             %cosmologyFunctions_ %expansionRate              (expansionFactor        )
    massHalo            =+basic                                 %mass                       (                       )
    
    zeta_0              =+self%virialDensityContrast_%densityContrast            (massHalo,expansionFactor=1.0d0          )
    zeta_z              =+self%virialDensityContrast_%densityContrast            (massHalo,expansionFactor=expansionFactor)
    zetaRate            =+self%virialDensityContrast_%densityContrastRateofChange(massHalo,expansionFactor=expansionFactor)

    ! Output the parameters for coreMass analytic calculation.
    call darkMatterProfile%floatRank0MetaPropertySet(self%massParticleID   ,+self%massParticle)
    call darkMatterProfile%floatRank0MetaPropertySet(self%massHaloID       ,+massHalo)
    call darkMatterProfile%floatRank0MetaPropertySet(self%zeta0ID          ,+zeta_0)
    call darkMatterProfile%floatRank0MetaPropertySet(self%zetazID          ,+zeta_z)
    call darkMatterProfile%floatRank0MetaPropertySet(self%expansionFactorID,+expansionFactor)

    A                   =+self%massParticle            &
         &               /8.0d-23
    K0                  =+self%beta                    &
         &               *A**(-1.5d0)
    K1                  =(sqrt(                        &
         &                     +zeta_z                 &
         &                     /zeta_0                 &
         &                    )                        &
         &                *massHalo                    &
         &                /self%gamma                  &
         &               )**self%alpha                 &
         &               *A**(1.5d0*(self%alpha-1.0d0))
    massCoreRate        =-0.5d0                        &
         &               *expansionFactor**(-0.5d0)    &
         &               *expansionRate                &
         &               *(K0+K1)                      &
         &               +self%alpha                   &
         &               *expansionFactor**(-0.5d0)    &
         &               *K1                           &
         &               *(+0.5d0                      &
         &                 *zetaRate                   &
         &                 /zeta_z                     &
         &                 +basic%accretionRate()      &
         &                 /massHalo                   &
         &                )
    call darkMatterProfile%floatRank0MetaPropertyRate(                 &
         &                                self%massCoreID, &
         &                                massCoreRate     &
         &                               )
    return
  end subroutine darkMatterProfileSolitonDifferentialEvolution

  subroutine darkMatterProfileSolitonNodeTreeInitialize(self,node)
    !!{
    Initialize solitonic core properties for a tree node.
    Computes the initial core mass using the analytic core–halo relation,
    records it as the minimum core mass for tolerance scaling, and stores the value in the node’s meta-property database.
    !!}
    use :: Galacticus_Nodes                , only : treeNode           , nodeComponentBasic       , nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorDarkMatterProfileSoliton), intent(inout), target :: self
    type            (treeNode                            ), intent(inout), target :: node
    class           (nodeComponentBasic                  ), pointer               :: basic
    class           (nodeComponentDarkMatterProfile      ), pointer               :: darkMatterProfile
    double precision                                                              :: massHalo                            , expansionFactor            , &
         &                                                                           zeta_0                              , zeta_z                     , &
         &                                                                           massCoreNormal

    ! Get required components.
    basic             => node%basic            ()
    darkMatterProfile => node%darkMatterProfile()
    ! Extract basic properties of the node.
    expansionFactor=+self             %cosmologyFunctions_% expansionFactor            (basic%time           ())
    massHalo       =+basic                                 %mass                       (                       )
    zeta_0             =+self%virialDensityContrast_%densityContrast(massHalo,expansionFactor=1.0d0          )
    zeta_z             =+self%virialDensityContrast_%densityContrast(massHalo,expansionFactor=expansionFactor)

    massCoreNormal     =+(                          & ! Equation (15) of Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
         &                +self%beta                &
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
         &                  /self%gamma             &
         &                 )**self%alpha            &
         &                *(                        &
         &                  +self%massParticle      &
         &                  /8.0d-23                &
         &                 )**(1.5d0*(self%alpha-1.0d0)) &
         &               )&
         &              /sqrt(expansionFactor)
    call darkMatterProfile%floatRank0MetaPropertySet(                 &
         &                               self%massCoreID, &
         &                               massCoreNormal   &
         &                              )
    return
  end subroutine darkMatterProfileSolitonNodeTreeInitialize
