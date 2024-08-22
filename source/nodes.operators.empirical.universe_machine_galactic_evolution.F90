!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  Implements a node operator class that inserts an empirical model of the formation history of a galaxy.
  Mass evolution is modeled using of the {\normalfont \ttfamily [UniverseMachine]} \citep{Behroozi_2019} 
  correlation between Galaxy Growth and Dark Matter Halo Assembly.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass

  !![
  <nodeOperator name="nodeOperatorEmpiricalGalaxyUniverseMachine">
   <description>
    A node operator class that inserts an empirical model of the formation history of a galaxy. 
    Mass evolution is modeled using of the {\normalfont \ttfamily [UniverseMachine]} \citep{Behroozi_2019} correlation between Galaxy Growth and Dark Matter Halo Assembly.    
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorEmpiricalGalaxyUniverseMachine
     !!{     
     A node operator class that inserts an empirical model of the formation history of a galaxy. 
     At each timestep and during mergers, the mass of the central galaxy is computed using the 
     Stellar Mass - Halo Mass using {\normalfont \ttfamily [UniverseMachine]} \citep{Behroozi_2019} fits.
     !!}
     private
     double precision                                                        :: epsilon_0, epsilon_a, epsilon_lna, &
         &                                                                      epsilon_z, M_0      , M_a        , &
         &                                                                      M_lna    , M_z      , alpha_0    , &
         &                                                                      alpha_a  , alpha_lna, alpha_z    , &
         &                                                                      beta_0   , beta_a   , beta_z     , &
         &                                                                      delta_0  , gamma_0  , gamma_a    , &
         &                                                                      gamma_z 
    class(cosmologyFunctionsClass),                                  pointer :: cosmologyFunctions_               => null()
   contains
     final     ::                                        empiricalGalaxyUniverseMachineDestructor
     procedure :: nodeInitialize                      => empiricalGalaxyUniverseMachineNodeInitialize
     procedure :: differentialEvolutionSolveAnalytics => empiricalGalaxyUniverseMachineSolveAnalytics
     procedure :: nodesMerge                          => empiricalGalaxyUniverseMachineNodesMerge
  end type nodeOperatorEmpiricalGalaxyUniverseMachine
  
  interface nodeOperatorEmpiricalGalaxyUniverseMachine
     !!{
     Constructors for the {\normalfont \ttfamily empiricalGalaxyUniverseMachine} node operator class.
     !!}
     module procedure empiricalGalaxyUniverseMachineConstructorParameters
     module procedure empiricalGalaxyUniverseMachineConstructorInternal
  end interface nodeOperatorEmpiricalGalaxyUniverseMachine
  
contains

  function empiricalGalaxyUniverseMachineConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily empiricalGalaxyUniverseMachine} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorEmpiricalGalaxyUniverseMachine)                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    double precision                                                            :: epsilon_0, epsilon_a, epsilon_lna, &
         &                                                                         epsilon_z, M_0      , M_a        , &
         &                                                                         M_lna    , M_z      , alpha_0    , &
         &                                                                         alpha_a  , alpha_lna, alpha_z    , &
         &                                                                         beta_0   , beta_a   , beta_z     , &
         &                                                                         delta_0  , gamma_0  , gamma_a    , &
         &                                                                         gamma_z
    class           (cosmologyFunctionsClass                  ), pointer        :: cosmologyFunctions_
    
    !![
    <inputParameter>
      <name>epsilon_0</name>
      <source>parameters</source>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\Epsilon_0$ (see Table J1 of \cite{Behroozi_2019})</description>
      <defaultValue>-1.435d0</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>epsilon_a</name>
      <source>parameters</source>
      <defaultValue>+1.831d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\Epsilon_a$  (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>epsilon_lna</name>
      <source>parameters</source>
      <defaultValue>+1.368d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\Epsilon_{\ln a} (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>epsilon_z</name>
      <source>parameters</source>
      <defaultValue>-0.217d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\Epsilon_0$ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>M_0</name>
      <source>parameters</source>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $M_0$ (see Table J1 of \cite{Behroozi_2019})</description>
      <defaultValue>+12.04d0</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>M_a</name>
      <source>parameters</source>
      <defaultValue>+4.556d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>M_lna</name>
      <source>parameters</source>
      <defaultValue>+4.417d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $M_{\ln a}$ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>M_z</name>
      <source>parameters</source>
      <defaultValue>-0.731d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $M_z$ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>alpha_0</name>
      <source>parameters</source>
      <defaultValue>+1.963d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\alpha_0$ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>alpha_a</name>
      <source>parameters</source>
      <defaultValue>-2.316d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\alpha_a$ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>alpha_lna</name>
      <source>parameters</source>
      <defaultValue>-1.732d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\alpha_{\ln a}$ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>alpha_z</name>
      <source>parameters</source>
      <defaultValue>+0.178d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\alpha_z$ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>beta_0</name>
      <source>parameters</source>
      <defaultValue>+0.482d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\beta_0$ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>beta_a</name>
      <source>parameters</source>
      <defaultValue>-0.841d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\beta_a$ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>beta_z</name>
      <source>parameters</source>
      <defaultValue>-0.471d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\beta_0$ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>delta_0</name>
      <source>parameters</source>
      <defaultValue>+0.411d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\delta_0$ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>gamma_0</name>
      <source>parameters</source>
      <defaultValue>-1.034d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\gamma_0$ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>gamma_a</name>
      <source>parameters</source>
      <defaultValue>-3.100d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\gamma_a$ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <inputParameter>
      <name>gamma_z</name>
      <source>parameters</source>
      <defaultValue>-1.055d0</defaultValue>
      <description>{\normalfont \ttfamily empiricalGalaxyUniverseMachine} Parameter $\gamma_z$ (see Table J1 of \cite{Behroozi_2019})</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=nodeOperatorEmpiricalGalaxyUniverseMachine(                          &
      &                                  epsilon_0              , &
      &                                  epsilon_a              , &
      &                                  epsilon_lna            , &
      &                                  epsilon_z              , &
      &                                  M_0                    , & 
      &                                  M_a                    , &
      &                                  M_lna                  , &
      &                                  M_z                    , &
      &                                  alpha_0                , &
      &                                  alpha_a                , &
      &                                  alpha_lna              , &
      &                                  alpha_z                , &
      &                                  beta_0                 , &
      &                                  beta_a                 , &
      &                                  beta_z                 , &
      &                                  delta_0                , &
      &                                  gamma_0                , &
      &                                  gamma_a                , &
      &                                  gamma_z                , &
      &                                  cosmologyFunctions_      &
      &                                 )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>  
    !!]
    return
  end function empiricalGalaxyUniverseMachineConstructorParameters

  function empiricalGalaxyUniverseMachineConstructorInternal(epsilon_0, epsilon_a, epsilon_lna, epsilon_z          , &
      &                                                      M_0      , M_a      , M_lna      , M_z                , &
      &                                                      alpha_0  , alpha_a  , alpha_lna  , alpha_z            , &
      &                                                      beta_0   , beta_a   , beta_z     , delta_0            , &
      &                                                      gamma_0  , gamma_a  , gamma_z    , cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily empiricalGalaxyUniverseMachine} node operator class.
    !!}
    implicit none
    type            (nodeOperatorEmpiricalGalaxyUniverseMachine)                     :: self
    double precision                                            , intent(in)         :: epsilon_0, epsilon_a, epsilon_lna, &
         &                                                                              epsilon_z, M_0      , M_a        , &
         &                                                                              M_lna    , M_z      , alpha_0    , &
         &                                                                              alpha_a  , alpha_lna, alpha_z    , &
         &                                                                              beta_0   , beta_a   , beta_z     , &
         &                                                                              delta_0  , gamma_0  , gamma_a    , &
         &                                                                              gamma_z
    class           (cosmologyFunctionsClass                   ), intent(in), target :: cosmologyFunctions_

    !![
    <constructorAssign variables="epsilon_0, epsilon_a, epsilon_lna, epsilon_z, M_0, M_a, M_lna, M_z, alpha_0, alpha_a, alpha_lna, alpha_z, beta_0, beta_a, beta_z, delta_0, gamma_0, gamma_a, gamma_z, *cosmologyFunctions_"/>
    !!]
    
    return
  end function empiricalGalaxyUniverseMachineConstructorInternal

  subroutine empiricalGalaxyUniverseMachineDestructor(self)    
    !!{
    Destructor for the {\normalfont \ttfamily empiricalGalaxyUniverseMachine} node operator class.
    !!}
    implicit none
    type(nodeOperatorEmpiricalGalaxyUniverseMachine), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]

    return
  end subroutine empiricalGalaxyUniverseMachineDestructor

  double precision function universeMachineScaling(z, y0, ya, ylna, yz) 
    !!{
    Implements the scaling relations provided in equations J3-J8 of \citep{Behroozi_2019}
    !!}
    implicit none
    double precision, intent(in) :: z   , y0, ya, &
         &                          ylna, yz   
    double precision             :: a

    a                      =+1.0d0  &
         &                 /(       &
         &                   +1.0d0 &
         &                   +z     &
         &                  )       
                             
    universeMachineScaling=+y0     &
         &                 +ya     &
         &                 *(      & 
         &                   +a    &
         &                   -1    &
         &                  )      &
         &                 -ylna   &
         &                 *log(a) &
         &                 +yz     &
         &                 *z      
    return
  end function universeMachineScaling

  double precision function universeMachineSMHM(self, haloMass, z)
    !!{
    Implements the Stellar Mass - Halo Mass relationship provided in equation J1 of \citep{Behroozi_2019}
    !!}
    implicit none  
    class           (nodeOperatorEmpiricalGalaxyUniverseMachine), intent(in), target :: self
    double precision                                            , intent(in)         :: haloMass  , z
    double precision                                                                 :: MLog10    , gammalog10, M1  , &
        &                                                                               epsilon   , alpha     , beta, &
        &                                                                               delta     , gamma     , x   , &
        &                                                                               smfm1Log10, smfm1     , powa, &
        &                                                                               powb      , expd  
    
    MLog10             =universeMachineScaling(z, self%M_0      , self%M_a      , self%M_lna      , self%M_z      )
    gammalog10         =universeMachineScaling(z, self%gamma_0  , self%gamma_a  , 0.0d0           , self%gamma_z  )

    epsilon            =universeMachineScaling(z, self%epsilon_0, self%epsilon_a, self%epsilon_lna, self%epsilon_z)
    alpha              =universeMachineScaling(z, self%alpha_0  , self%alpha_a  , self%alpha_lna  , self%alpha_z  )
    beta               =universeMachineScaling(z, self%beta_0   , self%beta_a   , 0.0d0           , self%beta_z   )
    delta              =universeMachineScaling(z, self%delta_0  , 0.0d0         , 0.0d0           , 0.0d0         )

    M1                 =+10.0d0**(MLog10    )
    gamma              =+10.0d0**(gammalog10)

    x                  =+log10(          &
      &                        +haloMass &
      &                        /M1       &
      &                       )             

    powa               =+10.0d0**(-alpha*x  )
    powb               =+10.0d0**(-beta*x   )

    expd               =-0.50d0  &
      &                 *(       &
      &                   +x     &
      &                   /delta &                     
      &                  )**2    

    smfm1Log10         =+epsilon     &
     &                  -log10(      &
     &                         +powa &
     &                         +powb &
     &                        )      &
     &                  +gamma       &
     &                  *exp  (      &
     &                         +expd &
     &                        )         

    smfm1              =10.0d0**(smfm1Log10 )

    universeMachineSMHM=+M1    &
      &                 *smfm1 

    return
 
  end function universeMachineSMHM

  subroutine empiricalGalaxyUniverseMachineNodeUpdate(self, node)
    !!{
    Updates the stellar mass of the node.
    !!} 
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSpheroid
    implicit none
    class           (nodeOperatorEmpiricalGalaxyUniverseMachine), intent(inout), target  :: self
    type            (treeNode                                  ), intent(inout), target  :: node
    class           (nodeComponentBasic                        ),                pointer :: basic  
    class           (nodeComponentSpheroid                     )               , pointer :: spheroid
    double precision                                                                     :: z, basicMass, massStellar


    if (.not.node%isOnMainBranch()) return
 
    basic           =>node %basic                                          (                   )

    z               = self %cosmologyFunctions_%redshiftFromExpansionFactor(    &
         &                         self %cosmologyFunctions_%expansionFactor(   &
         &                         basic                    %time            () &
         &                                                                  )   &
         &                                                                  )    
    
    basicMass       = basic              %mass                              (                  )

    massStellar     = universeMachineSMHM                                   (self, basicMass, z)
    spheroid        =>node               %spheroid                          (                  )
    
    call spheroid%    massStellarSet                                        (massStellar       )
  
    return
  end subroutine empiricalGalaxyUniverseMachineNodeUpdate

  subroutine empiricalGalaxyUniverseMachineNodeInitialize(self,node)
    !!{
    Initialize the galaxy
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid
    implicit none
    class        (nodeOperatorEmpiricalGalaxyUniverseMachine), intent(inout), target  :: self
    type         (treeNode                                  ), intent(inout), target  :: node
    class        (nodeComponentSpheroid                     )               , pointer :: spheroid

    spheroid=>node%spheroid(autoCreate=.true.)

    call empiricalGalaxyUniverseMachineNodeUpdate(self, node)
    return
  end subroutine empiricalGalaxyUniverseMachineNodeInitialize

  subroutine empiricalGalaxyUniverseMachineSolveAnalytics(self,node,time)
    !!{
    Update galactic properties at each timestep
    !!}
    implicit none
    class           (nodeOperatorEmpiricalGalaxyUniverseMachine), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node        
    double precision                                            , intent(in   ) :: time

    call empiricalGalaxyUniverseMachineNodeUpdate(self, node)
    return
  end subroutine empiricalGalaxyUniverseMachineSolveAnalytics

  subroutine empiricalGalaxyUniverseMachineNodesMerge(self,node)
    !!{
    Update galactic properties at each merger
    !!}
    implicit none
    class        (nodeOperatorEmpiricalGalaxyUniverseMachine), intent(inout) :: self
    type         (treeNode                                  ), intent(inout) :: node        

    call empiricalGalaxyUniverseMachineNodeUpdate(self, node)
    return
  end subroutine empiricalGalaxyUniverseMachineNodesMerge
