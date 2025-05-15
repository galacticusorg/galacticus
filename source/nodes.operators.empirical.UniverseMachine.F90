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

!+    Contributions to this file made by: Charles Gannon, Andrew Benson.

  !!{  
  Implements a node operator that inserts an empirical model of the formation history of a galaxy. Mass evolution is modeled
  using the \textsc{UniverseMachine} \citep{behroozi_universemachine_2019} correlation between galaxy growth and dark matter halo
  assembly.
  !!}

  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Cosmology_Functions    , only : cosmologyFunctions        , cosmologyFunctionsClass
  use :: Virial_Density_Contrast, only : virialDensityContrastClass, virialDensityContrastBryanNorman1998
  
  !![
  <nodeOperator name="nodeOperatorEmpiricalGalaxyUniverseMachine">
   <description>
     A node operator that inserts an empirical model of the formation history of a galaxy. Mass evolution is modeled using the
     \textsc{UniverseMachine} \citep{behroozi_universemachine_2019} correlation between galaxy growth and dark matter halo
     assembly. The \textsc{UniverseMachine} fits are used only for redshifts less than {\normalfont \ttfamily [redshiftMaximum]},
     and for halo masses above {\normalfont \ttfamily [massHaloMinimum]}. Outside of those ranges, no galaxy is inserted.
   </description>
   <deepCopy>
    <functionClass variables="virialDensityContrastDefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="virialDensityContrastDefinition_"/>
   </stateStorable>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorEmpiricalGalaxyUniverseMachine
     !!{     
     A {\normalfont \ttfamily nodeOperator} class that inserts an empirical model of the formation history of a galaxy.  At each
     time step and during mergers, the mass of the central galaxy is computed using the stellar mass--halo mass relation using
     \textsc{UniverseMachine} \citep{behroozi_universemachine_2019} fits.
     !!}
     private
     double precision                                       :: massStellarFinal                          , fractionMassSpheroid, fractionMassDisk           , &
          &                                                    epsilon_0                                 , epsilon_a           , epsilon_lna     , epsilon_z, &
          &                                                    M_0                                       , M_a                 , M_lna           , M_z      , &
          &                                                    alpha_0                                   , alpha_a             , alpha_lna       , alpha_z  , &
          &                                                    beta_0                                    , beta_a              , beta_z                     , &
          &                                                    gamma_0                                   , gamma_a             , gamma_z                    , &
          &                                                    delta_0                                                                                      , &
          &                                                    redshiftMaximum                           , massHaloMinimum
     logical                                                :: setFinalStellarMass                       , hasDisk             , hasSpheroid
     class  (cosmologyParametersClass            ), pointer :: cosmologyParameters_             => null()
     class  (cosmologyFunctionsClass             ), pointer :: cosmologyFunctions_              => null()
     class  (virialDensityContrastClass          ), pointer :: virialDensityContrast_           => null()
     type   (virialDensityContrastBryanNorman1998), pointer :: virialDensityContrastDefinition_ => null()
   contains
     !![
     <methods>
       <method method="scaling"                     description="Compute the scaling of \textsc{UniverseMachine} parameters with redshift."/>
       <method method="stellarMassHaloMassRelation" description="Evaluate the stellar mass--halo mass relation."                           />
       <method method="update"                      description="Update the stellar mass of the galaxy."                                   />
     </methods>
     !!]
     final     ::                                        empiricalGalaxyUniverseMachineDestructor
     procedure :: scaling                             => empiricalGalaxyUniverseMachineScaling
     procedure :: stellarMassHaloMassRelation         => empiricalGalaxyUniverseMachineStellarMassHaloMassRelation
     procedure :: update                              => empiricalGalaxyUniverseMachineUpdate
     procedure :: nodeInitialize                      => empiricalGalaxyUniverseMachineNodeInitialize
     procedure :: differentialEvolutionAnalytics      => empiricalGalaxyUniverseMachineAnalytics
     procedure :: differentialEvolutionSolveAnalytics => empiricalGalaxyUniverseMachineSolveAnalytics
     procedure :: nodesMerge                          => empiricalGalaxyUniverseMachineNodesMerge
  end type nodeOperatorEmpiricalGalaxyUniverseMachine
  
  interface nodeOperatorEmpiricalGalaxyUniverseMachine
     !!{
     Constructors for the {\normalfont \ttfamily empiricalGalaxyUniverseMachine} {\normalfont \ttfamily nodeOperator} class.
     !!}
     module procedure empiricalGalaxyUniverseMachineConstructorParameters
     module procedure empiricalGalaxyUniverseMachineConstructorInternal
  end interface nodeOperatorEmpiricalGalaxyUniverseMachine
  
contains

  function empiricalGalaxyUniverseMachineConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily empiricalGalaxyUniverseMachine} {\normalfont \ttfamily nodeOperator} class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorEmpiricalGalaxyUniverseMachine)                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    double precision                                                            :: massStellarFinal      , fractionMassSpheroid, fractionMassDisk           , &
          &                                                                        epsilon_0             , epsilon_a           , epsilon_lna     , epsilon_z, &
          &                                                                        M_0                   , M_a                 , M_lna           , M_z      , &
          &                                                                        alpha_0               , alpha_a             , alpha_lna       , alpha_z  , &
          &                                                                        beta_0                , beta_a              , beta_z                     , &
          &                                                                        gamma_0               , gamma_a             , gamma_z                    , &
          &                                                                        delta_0                                                                  , &
          &                                                                        redshiftMaximum       , massHaloMinimum
    class           (cosmologyParametersClass                 ), pointer        :: cosmologyParameters_
    class           (cosmologyFunctionsClass                  ), pointer        :: cosmologyFunctions_
    class           (virialDensityContrastClass               ), pointer        :: virialDensityContrast_

    !![ 
    <inputParameter>
      <name>massStellarFinal</name>
      <source>parameters</source>
      <description>
        If positive, rescale the \textsc{UniverseMachine} fitting functions to match this final mass. A negative value indicates
        the final mass of the galaxy will be determined by \textsc{UniverseMachine} fits.
      </description>
      <defaultValue>-1.0d0</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>fractionMassSpheroid</name>
      <source>parameters</source>
      <description>
        Sets the fraction of galaxy mass belonging to the spheroid component.
      </description>
      <defaultValue>1.0d0</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>fractionMassDisk</name>
      <source>parameters</source>
      <description>
        Sets the fraction of galaxy mass belonging to the disk component.
      </description>
      <defaultValue>0.0d0</defaultValue>
    </inputParameter>       
    <inputParameter>
      <name>epsilon_0</name>
      <source>parameters</source>
      <description>Parameter $\epsilon_0$ of the \textsc{UniverseMachine} fits.</description>
      <defaultValue>-1.435d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
    </inputParameter>
    <inputParameter>
      <name>epsilon_a</name>
      <source>parameters</source>
      <defaultValue>+1.831d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $\epsilon_a$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>epsilon_lna</name>
      <source>parameters</source>
      <defaultValue>+1.368d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $\epsilon_{\ln a}$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>epsilon_z</name>
      <source>parameters</source>
      <defaultValue>-0.217d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $\epsilon_0$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>M_0</name>
      <source>parameters</source>
      <description>Parameter $M_0$ of the \textsc{UniverseMachine} fits.</description>
      <defaultValue>+12.04d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
    </inputParameter>
    <inputParameter>
      <name>M_a</name>
      <source>parameters</source>
      <defaultValue>+4.556d0</defaultValue>
       <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
     <description>Parameter $M_a$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>M_lna</name>
      <source>parameters</source>
      <defaultValue>+4.417d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $M_{\ln a}$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>M_z</name>
      <source>parameters</source>
      <defaultValue>-0.731d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $M_z$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>alpha_0</name>
      <source>parameters</source>
      <defaultValue>+1.963d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $\alpha_0$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>alpha_a</name>
      <source>parameters</source>
      <defaultValue>-2.316d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $\alpha_a$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>alpha_lna</name>
      <source>parameters</source>
      <defaultValue>-1.732d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $\alpha_{\ln a}$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>alpha_z</name>
      <source>parameters</source>
      <defaultValue>+0.178d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $\alpha_z$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>beta_0</name>
      <source>parameters</source>
      <defaultValue>+0.482d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $\beta_0$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>beta_a</name>
      <source>parameters</source>
      <defaultValue>-0.841d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $\beta_a$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>beta_z</name>
      <source>parameters</source>
      <defaultValue>-0.471d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $\beta_0$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>delta_0</name>
      <source>parameters</source>
      <defaultValue>+0.411d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $\delta_0$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>gamma_0</name>
      <source>parameters</source>
      <defaultValue>-1.034d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $\gamma_0$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>gamma_a</name>
      <source>parameters</source>
      <defaultValue>-3.100d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $\gamma_a$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter>
    <inputParameter>
      <name>gamma_z</name>
      <source>parameters</source>
      <defaultValue>-1.055d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>Parameter $\gamma_z$ of the \textsc{UniverseMachine} fits.</description>
    </inputParameter> 
    <inputParameter>
      <name>redshiftMaximum</name>
      <source>parameters</source>
      <defaultValue>10.0d0</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>The maximum redshift at which UniverseMachine fits will be applied.</description>
    </inputParameter> 
    <inputParameter>
      <name>massHaloMinimum</name>
      <source>parameters</source>
      <defaultValue>1.0d10</defaultValue>
      <defaultSource>\cite[][Table J1]{behroozi_universemachine_2019}</defaultSource>
      <description>The minimum halo mass at which UniverseMachine fits will be applied.</description>
    </inputParameter> 
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self=nodeOperatorEmpiricalGalaxyUniverseMachine(                                                                            &
         &                                          massStellarFinal    ,fractionMassSpheroid, fractionMassDisk               , &
         &                                          epsilon_0           ,epsilon_a           ,epsilon_lna           ,epsilon_z, &
         &                                          M_0                 ,M_a                 ,M_lna                 ,M_z      , &
         &                                          alpha_0             ,alpha_a             ,alpha_lna             ,alpha_z  , &
         &                                          beta_0              ,beta_a              ,beta_z                          , &
         &                                          gamma_0             ,gamma_a             ,gamma_z                         , &
         &                                          delta_0                                                                   , &
         &                                          redshiftMaximum     ,massHaloMinimum                                      , &
         &                                          cosmologyParameters_,cosmologyFunctions_ ,virialDensityContrast_            &
         &                                         )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>  
    !!]
    return
  end function empiricalGalaxyUniverseMachineConstructorParameters

  function empiricalGalaxyUniverseMachineConstructorInternal(                                                                            &
         &                                                   massStellarFinal    ,fractionMassSpheroid,fractionMassDisk                , &
         &                                                   epsilon_0           ,epsilon_a           ,epsilon_lna           ,epsilon_z, &
         &                                                   M_0                 ,M_a                 ,M_lna                 ,M_z      , &
         &                                                   alpha_0             ,alpha_a             ,alpha_lna             ,alpha_z  , &
         &                                                   beta_0              ,beta_a              ,beta_z                          , &
         &                                                   gamma_0             ,gamma_a             ,gamma_z                         , &
         &                                                   delta_0                                                                   , &
         &                                                   redshiftMaximum     ,massHaloMinimum                                      , &
         &                                                   cosmologyParameters_,cosmologyFunctions_ ,virialDensityContrast_            &
         &                                                  ) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily empiricalGalaxyUniverseMachine} {\normalfont \ttfamily nodeOperator} class.
    !!}
    use :: Error               , only : Error_Report
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    type            (nodeOperatorEmpiricalGalaxyUniverseMachine)                        :: self
    double precision                                            , intent(in)            :: massStellarFinal      , fractionMassSpheroid, fractionMassDisk, &
         &                                                                                 epsilon_0             , epsilon_a           , epsilon_lna     , &
         &                                                                                 epsilon_z             , M_0                 , M_a             , &
         &                                                                                 M_lna                 , M_z                 , alpha_0         , &
         &                                                                                 alpha_a               , alpha_lna           , alpha_z         , &
         &                                                                                 beta_0                , beta_a              , beta_z          , &
         &                                                                                 delta_0               , gamma_0             , gamma_a         , &
         &                                                                                 gamma_z                                                       , &
         &                                                                                 redshiftMaximum       , massHaloMinimum
    class           (cosmologyParametersClass                  ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass                   ), intent(in   ), target :: cosmologyFunctions_
    class           (virialDensityContrastClass                ), intent(in   ), target :: virialDensityContrast_
    !![
    <constructorAssign variables="massStellarFinal, fractionMassSpheroid, fractionMassDisk, epsilon_0, epsilon_a, epsilon_lna, epsilon_z, M_0, M_a, M_lna, M_z, alpha_0, alpha_a, alpha_lna, alpha_z, beta_0, beta_a, beta_z, delta_0, gamma_0, gamma_a, gamma_z, redshiftMaximum, massHaloMinimum, *cosmologyParameters_, *cosmologyFunctions_, *virialDensityContrast_"/>
    !!]
    
    self%setFinalStellarMass=massStellarFinal     >= 0.0d0
    self%hasDisk            =fractionMassDisk     >  0.0d0
    self%hasSpheroid        =fractionMassSpheroid >  0.0d0
    ! Validate fractions.
    if (.not.Values_Agree(fractionMassDisk+fractionMassSpheroid,1.0d0,absTol=1.0d-3)) call Error_Report('disk and spheroid fractions do not sum to 1'//{introspection:location})
    if (fractionMassDisk     < 0.0d0) call Error_Report('negative disk fraction is unphysical'    //{introspection:location})
    if (fractionMassSpheroid < 0.0d0) call Error_Report('negative spheroid fraction is unphysical'//{introspection:location})
    ! Create virial density contrast definition.
    allocate(self%virialDensityContrastDefinition_)
    !![
    <referenceConstruct isResult="yes" owner="self" object="virialDensityContrastDefinition_" constructor="virialDensityContrastBryanNorman1998(allowUnsupportedCosmology=.false.,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_)"/>
    !!]
    return
  end function empiricalGalaxyUniverseMachineConstructorInternal

  subroutine empiricalGalaxyUniverseMachineDestructor(self)    
    !!{
    Destructor for the {\normalfont \ttfamily empiricalGalaxyUniverseMachine} {\normalfont \ttfamily nodeOperator} class.
    !!}
    implicit none
    type(nodeOperatorEmpiricalGalaxyUniverseMachine), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%virialDensityContrast_"          />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    !!]
    return
  end subroutine empiricalGalaxyUniverseMachineDestructor

  double precision function empiricalGalaxyUniverseMachineScaling(self,redshift,y0,ya,ylna,yz) result(parameterScaled)
    !!{
    Implements the scaling relations provided in equations J3--J8 of \cite{behroozi_universemachine_2019}.
    !!}
    implicit none
    class           (nodeOperatorEmpiricalGalaxyUniverseMachine), intent(in) :: self
    double precision                                            , intent(in) :: y0             , ya, &
         &                                                                      ylna           , yz, &
         &                                                                      redshift
    double precision                                                         :: expansionFactor

    expansionFactor=+self%cosmologyFunctions_%expansionFactorFromRedshift(redshift)
    parameterScaled=+y0                               &
         &          +ya  *   (+expansionFactor-1.0d0) &
         &          -ylna*log(+expansionFactor      ) &
         &          +yz  *     redshift      
    return
  end function empiricalGalaxyUniverseMachineScaling

  double precision function empiricalGalaxyUniverseMachineStellarMassHaloMassRelation(self,massHalo,redshift) result(massStellar)
    !!{
    Implements the stellar mass--halo mass relationship provided in equation J1 of \cite{behroozi_universemachine_2019}.
    !!}
    implicit none  
    class           (nodeOperatorEmpiricalGalaxyUniverseMachine), intent(in) :: self 
    double precision                                            , intent(in) :: massHalo          , redshift
    double precision                                            , parameter  :: expMaximum=100.0d0
    double precision                                                         :: MLog10            , gammalog10, &
         &                                                                      M1                , epsilon   , &
         &                                                                      alpha             , beta      , &
         &                                                                      delta             , gamma     , &
         &                                                                      x                 , smfm1Log10, &
         &                                                                      smfm1             , expd      , &
         &                                                                      expab             , ab
    
    MLog10     =self%scaling(redshift,self%M_0      ,self%M_a      ,self%M_lna      ,self%M_z      )
    gammalog10 =self%scaling(redshift,self%gamma_0  ,self%gamma_a  ,     0.0d0      ,self%gamma_z  )
    epsilon    =self%scaling(redshift,self%epsilon_0,self%epsilon_a,self%epsilon_lna,self%epsilon_z)
    alpha      =self%scaling(redshift,self%alpha_0  ,self%alpha_a  ,self%alpha_lna  ,self%alpha_z  )
    beta       =self%scaling(redshift,self%beta_0   ,self%beta_a   ,     0.0d0      ,self%beta_z   )
    delta      =self%scaling(redshift,self%delta_0  ,     0.0d0    ,     0.0d0      ,     0.0d0    )
    M1         =+10.0d0**MLog10
    gamma      =+10.0d0**gammalog10
    x          =+log10(          &
         &             +massHalo &
         &             /M1       &
         &            )
    expd       =+0.50d0  &
         &      *(       &
         &        +x     &
         &        /delta &                     
         &       )**2
    expab      =+(       &
         &        +alpha &
         &        -beta  &
         &       )       &
         &      *x
    if (expab > expMaximum) then
       ab  =+beta *x
    else
       ab  =+alpha*x              &
         &  -log10(               &
         &         +1.0d0         &
         &         +10.0d0**expab &
         &        )   
       end if
    smfm1Log10 =+epsilon     &
         &      +ab          &
         &      +gamma       &
         &      *exp  (      &
         &             -expd &
         &            )         
    smfm1      =+10.0d0**smfm1Log10
    massStellar=+M1*smfm1 
    return
  end function empiricalGalaxyUniverseMachineStellarMassHaloMassRelation

  subroutine empiricalGalaxyUniverseMachineUpdate(self,node)
    !!{
    Updates the stellar mass of the node.
    !!} 
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , nodeComponentSpheroid, nodeComponentDisk
    implicit none
    class           (nodeOperatorEmpiricalGalaxyUniverseMachine), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    type            (treeNode                                  ), pointer       :: nodeRoot     
    class           (nodeComponentBasic                        ), pointer       :: basic      , basicRoot
    class           (nodeComponentSpheroid                     ), pointer       :: spheroid
    class           (nodeComponentDisk                         ), pointer       :: disk
    double precision                                                            :: redshift   , redshiftRoot   , &
         &                                                                         massHalo   , massHaloRoot   , &
         &                                                                         massStellar, massStellarRoot

    ! Insert a galaxy on the main branch of the tree only.
    if (.not.node%isOnMainBranch()) return
    ! Do not insert galaxies below the halo mass limit, or before the earliest time.
    basic       => node   %basic                                            (                 )
    redshift    =  self   %cosmologyFunctions_%redshiftFromExpansionFactor(                     &
        &           self  %cosmologyFunctions_%expansionFactor             (                    &
        &            basic%time                                             (                 ) &
        &                                                                  )                    &
        &                                                                 )
    massHalo    =  Dark_Matter_Profile_Mass_Definition                    (                                                                                                       &
         &                                                                                        node                                                                          , &
         &                                                                                        self%virialDensityContrastDefinition_%densityContrast(                          &
         &                                                                                                                                              basic%mass            (), &
         &                                                                                                                                              basic%timeLastIsolated()  &
         &                                                                                                                                             )                        , &
         &                                                                 cosmologyParameters_  =self%cosmologyParameters_                                                     , &
         &                                                                 cosmologyFunctions_   =self%cosmologyFunctions_                                                      , &
         &                                                                 virialDensityContrast_=self%virialDensityContrast_                                                     &
         &                                                                ) 
    if     (                                 &
         &   redshift > self%redshiftMaximum &
         &  .or.                             &
         &   massHalo < self%massHaloMinimum &
         & ) return
    ! Find the stellar mass from the stellar mass-halo mass relation.
    massStellar =  self   %stellarMassHaloMassRelation                      (massHalo,redshift)     
    ! If necessary, rescale the stellar mass to ensure that we match the requested final stellar mass.
    if (self%setFinalStellarMass) then 
      ! Compute the stellar mass at the root node.
      nodeRoot        => node%hostTree%nodeBase
      basicRoot       => nodeRoot   %basic                                            (                         )     
      massHaloRoot    =  Dark_Matter_Profile_Mass_Definition                          (                                                                                                           &
         &                                                                                                    nodeRoot                                                                          , &
         &                                                                                                    self%virialDensityContrastDefinition_%densityContrast(                              &
         &                                                                                                                                                          basicRoot%mass            (), &
         &                                                                                                                                                          basicRoot%timeLastIsolated()  &
         &                                                                                                                                                         )                            , &
         &                                                                             cosmologyParameters_  =self%cosmologyParameters_                                                         , &
         &                                                                             cosmologyFunctions_   =self%cosmologyFunctions_                                                          , &
         &                                                                             virialDensityContrast_=self%virialDensityContrast_                                                         &
         &                                                                            ) 
      redshiftRoot    =  self       %cosmologyFunctions_%redshiftFromExpansionFactor(                             &
          &               self      %cosmologyFunctions_%expansionFactor             (                            &
          &                basicRoot%time                                             (                         ) &
          &                                                                          )                            &
          &                                                                         )
      massStellarRoot =  self       %stellarMassHaloMassRelation                      (massHaloRoot,redshiftRoot)     
      ! Rescale the stellar mass to ensure that we match the requested final stellar mass.
      massStellar     =+     massStellar      &
        &              /     massStellarRoot  &
        &              *self%massStellarFinal  
    end if
    ! Ensure stellar mass remains non-negative.
    massStellar=max(massStellar,0.0d0)
    ! Set the stellar mass, partitioned between disk and spheroid.
    if (self%hasSpheroid) then
      spheroid => node    %spheroid      (                                     )
      call        spheroid%massStellarSet(self%fractionMassSpheroid*massStellar) 
    end if 
    if (self%hasDisk    ) then
      disk     => node    %disk          (                                     )
      call        disk    %massStellarSet(self%fractionMassDisk    *massStellar) 
    end if
    return
  end subroutine empiricalGalaxyUniverseMachineUpdate

  subroutine empiricalGalaxyUniverseMachineNodeInitialize(self,node)
    !!{
    Initialize the galaxy.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid
    implicit none
    class(nodeOperatorEmpiricalGalaxyUniverseMachine), intent(inout), target  :: self
    type (treeNode                                  ), intent(inout), target  :: node
    class(nodeComponentSpheroid                     ),                pointer :: spheroid
    class(nodeComponentDisk                         ),                pointer :: disk
    
    if (.not.node%isOnMainBranch().or.associated(node%firstChild)) return 
    if (self%hasSpheroid) spheroid => node%spheroid(autoCreate=.true.)
    if (self%hasDisk    ) disk     => node%disk    (autoCreate=.true.)
    call self%update(node)
    return
  end subroutine empiricalGalaxyUniverseMachineNodeInitialize
  
  subroutine empiricalGalaxyUniverseMachineAnalytics(self,node)
    !!{
    Mark disk and spheroid stellar masses as analytically-solvable.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid
    implicit none
    class(nodeOperatorEmpiricalGalaxyUniverseMachine), intent(inout) :: self
    type (treeNode                                  ), intent(inout) :: node
    class(nodeComponentDisk                         ), pointer       :: disk
    class(nodeComponentSpheroid                     ), pointer       :: spheroid

    disk     => node%disk    ()
    spheroid => node%spheroid()
    select type (disk)
    type is (nodeComponentDisk)
       ! Disk does not exist.
    class default
       call disk    %massStellarAnalytic()
    end select
    select type (spheroid)
    type is (nodeComponentSpheroid)
       ! Spheroid does not exist.
    class default
       call spheroid%massStellarAnalytic()     
    end select
    return
  end subroutine empiricalGalaxyUniverseMachineAnalytics
  
  subroutine empiricalGalaxyUniverseMachineSolveAnalytics(self,node,time)
    !!{
    Update galactic properties at each timestep
    !!}
    implicit none
    class           (nodeOperatorEmpiricalGalaxyUniverseMachine), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node        
    double precision                                            , intent(in   ) :: time
    !$GLC attributes unused :: time
    
    call self%update(node)
    return
  end subroutine empiricalGalaxyUniverseMachineSolveAnalytics
  
  subroutine empiricalGalaxyUniverseMachineNodesMerge(self,node)
    !!{
    Update galactic properties at each merger.
    !!}
    implicit none
    class(nodeOperatorEmpiricalGalaxyUniverseMachine), intent(inout) :: self
    type (treeNode                                  ), intent(inout) :: node        

    call self%update(node)
    return
  end subroutine empiricalGalaxyUniverseMachineNodesMerge
