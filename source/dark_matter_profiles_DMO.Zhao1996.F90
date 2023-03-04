!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  An implementation of \cite{zhao_analytical_1996} dark matter halo profiles.
  !!}

  !![
  <enumeration>
   <name>specialCase</name>
   <description>Special cases for {\normalfont \ttfamily zhao1996} dark matter halo profile class.</description>
   <entry label="general"    />
   <entry label="coredNFW"   />
   <entry label="gamma0_5NFW"/>
   <entry label="NFW"        />
   <entry label="gamma1_5NFW"/>
  </enumeration>
  !!]

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOZhao1996">
   <description>
    A dark matter profile DMO class which implements the \cite{zhao_analytical_1996} density profile
    \begin{equation}
      \rho_\mathrm{dark matter}(r) \propto \left({r\over r_\mathrm{s}}\right)^{-\gamma} \left(1+\left[{r\over r_\mathrm{s}}\right]^\alpha \right)^{-(\beta-\gamma)/\alpha},
    \end{equation}
    normalized such that the total mass of the \gls{node} is enclosed with the virial radius and with the scale length    
    $r_\mathrm{s}$. Numerical solutions are implemented for the case of general $(\alpha,\beta,\gamma)$, with some analytic
    solution for certain special cases.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOZhao1996
     !!{
     A dark matter halo profile class implementing \cite{zhao_analytical_1996} dark matter halos.
     !!}
     private
     type            (enumerationSpecialCaseType) :: specialCase
     double precision                             :: alpha      , beta, &
          &                                          gamma
   contains
     !![
     <methods>
       <method description="Reset memoized calculations."                                                                                              method="calculationReset"/>
       <method description="Provides the values of the $(\alpha,\beta,\gamma)$ exponents in \cite{zhao_analytical_1996} dark matter density profiles." method="exponents"       />
       <method description="Compute the scale radius of the profile."                                                                                  method="scaleRadius"     />
       <method description="Compute the normalization of the profile."                                                                                 method="normalization"   />
       <method description="Compute the unnormalized mass within the given scale-free radius"                                                          method="massUnnormalized"/>
     </methods>
     !!]
     final     ::                                      zhao1996Destructor
     procedure :: autoHook                          => zhao1996AutoHook
     procedure :: calculationReset                  => zhao1996CalculationReset
     procedure :: exponents                         => zhao1996Exponents
     procedure :: scaleRadius                       => zhao1996ScaleRadius
     procedure :: normalization                     => zhao1996Normalization
     procedure :: massUnnormalized                  => zhao1996MassUnnormalized
     procedure :: density                           => zhao1996Density
     procedure :: densityLogSlope                   => zhao1996DensityLogSlope
     procedure :: enclosedMass                      => zhao1996EnclosedMass
     procedure :: radiusEnclosingMass               => zhao1996RadiusEnclosingMass
     procedure :: radiusEnclosingDensity            => zhao1996RadiusEnclosingDensity
     procedure :: potential                         => zhao1996Potential
     procedure :: circularVelocity                  => zhao1996CircularVelocity
     procedure :: circularVelocityMaximum           => zhao1996CircularVelocityMaximum
     procedure :: radiusCircularVelocityMaximum     => zhao1996RadiusCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => zhao1996RadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => zhao1996RadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => zhao1996RotationNormalization
     procedure :: energy                            => zhao1996Energy
     procedure :: kSpace                            => zhao1996KSpace
     procedure :: freefallRadius                    => zhao1996FreefallRadius
     procedure :: freefallRadiusIncreaseRate        => zhao1996FreefallRadiusIncreaseRate
     procedure :: radialMoment                      => zhao1996RadialMoment
  end type darkMatterProfileDMOZhao1996

  interface darkMatterProfileDMOZhao1996
     !!{
     Constructors for the {\normalfont \ttfamily zhao1996} dark matter halo profile class.
     !!}
     module procedure zhao1996ConstructorParameters
     module procedure zhao1996ConstructorInternal
  end interface darkMatterProfileDMOZhao1996

  ! Sub-module scope variables used in numerical solutions.
  class           (darkMatterProfileDMOZhao1996), pointer :: self_
  type            (treeNode                    ), pointer :: node_
  double precision                                        :: time_, radiusFreefall_
  !$omp threadprivate(self_,node_,time_,radiusFreefall_)
  
contains

  function zhao1996ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily zhao1996} dark matter halo profile class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOZhao1996)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass    ), pointer       :: darkMatterHaloScale_
    double precision                                              :: alpha               , beta, &
         &                                                           gamma
    
    !![
    <inputParameter>
      <name>alpha</name>
      <source>parameters</source>
      <description>The parameter $\alpha$ of the \cite{zhao_analytical_1996} dark matter density profile.</description>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <description>The parameter $\beta$ of the \cite{zhao_analytical_1996} dark matter density profile.</description>
    </inputParameter>
    <inputParameter>
      <name>gamma</name>
      <source>parameters</source>
      <description>The parameter $\gamma$ of the \cite{zhao_analytical_1996} dark matter density profile.</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOZhao1996(alpha,beta,gamma,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function zhao1996ConstructorParameters

  function zhao1996ConstructorInternal(alpha,beta,gamma,darkMatterHaloScale_) result(self)
    !!{
    Generic constructor for the {\normalfont \ttfamily zhao1996} dark matter halo profile class.
    !!}
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    type            (darkMatterProfileDMOZhao1996)                        :: self
    class           (darkMatterHaloScaleClass    ), intent(in   ), target :: darkMatterHaloScale_
    double precision                              , intent(in   )         :: alpha               , beta, &
         &                                                                   gamma
    !![
    <constructorAssign variables="alpha, beta, gamma, *darkMatterHaloScale_"/>
    !!]

    ! Detect special cases.
    if      (                                         &
         &    Values_Agree(alpha,1.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(beta ,3.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(gamma,1.0d0,absTol=1.0d-6) &
         &  ) then
       ! The "NFW" profile.
       self%specialCase=specialCaseNFW
    else if (                                         &
         &    Values_Agree(alpha,1.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(beta ,3.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(gamma,0.0d0,absTol=1.0d-6) &
         &  ) then
       ! The "cored NFW" profile.
       self%specialCase=specialCaseCoredNFW
    else if (                                         &
         &    Values_Agree(alpha,1.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(beta ,3.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(gamma,0.5d0,absTol=1.0d-6) &
         &  ) then
       ! The "γ=1/2 NFW" profile.
       self%specialCase=specialCaseGamma0_5NFW
    else if (                                         &
         &    Values_Agree(alpha,1.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(beta ,3.0d0,absTol=1.0d-6) &
         &   .and.                                    &
         &    Values_Agree(gamma,1.5d0,absTol=1.0d-6) &
         &  ) then
       ! The "γ=3/2 NFW" profile.
       self%specialCase=specialCaseGamma1_5NFW
    else
       ! Use general solutions.
       self%specialCase=specialCaseGeneral
    end if
    return
  end function zhao1996ConstructorInternal

  subroutine zhao1996AutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileDMOZhao1996), intent(inout) :: self

    call calculationResetEvent%attach(self,zhao1996CalculationReset,openMPThreadBindingAllLevels,label='darkMatterProfileDMOZhao1996')
    return
  end subroutine zhao1996AutoHook

  subroutine zhao1996Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily zhao1996} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOZhao1996), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine zhao1996Destructor

  subroutine zhao1996CalculationReset(self,node)
    !!{
    Reset the dark matter profile calculation.
    !!}
    implicit none
    class(darkMatterProfileDMOZhao1996), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    self%genericLastUniqueID                         =node%uniqueID()
    self%genericEnclosedMassRadiusMinimum            =+huge(0.0d0)
    self%genericEnclosedMassRadiusMaximum            =-huge(0.0d0)
    self%genericVelocityDispersionRadialRadiusMinimum=+huge(0.0d0)
    self%genericVelocityDispersionRadialRadiusMaximum=-huge(0.0d0)
    if (allocated(self%genericVelocityDispersionRadialVelocity)) deallocate(self%genericVelocityDispersionRadialVelocity)
    if (allocated(self%genericVelocityDispersionRadialRadius  )) deallocate(self%genericVelocityDispersionRadialRadius  )
    if (allocated(self%genericEnclosedMassMass                )) deallocate(self%genericEnclosedMassMass                )
    if (allocated(self%genericEnclosedMassRadius              )) deallocate(self%genericEnclosedMassRadius              )
    return
  end subroutine zhao1996CalculationReset

  subroutine zhao1996Exponents(self,node,alpha,beta,gamma)
    !!{
    Compute the exponents of the {\normalfont \ttfamily zhao1996} dark matter halo profile.
    !!}
    implicit none
    class           (darkMatterProfileDMOZhao1996  ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(  out) :: alpha, beta, &
         &                                                             gamma

    alpha=self%alpha
    beta =self%beta
    gamma=self%gamma
    return
  end subroutine zhao1996Exponents
  
  double precision function zhao1996ScaleRadius(self,node)
    !!{
    Compute the scale radius of the {\normalfont \ttfamily zhao1996} dark matter halo profile.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(darkMatterProfileDMOZhao1996  ), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile

    darkMatterProfile   => node             %darkMatterProfile()
    zhao1996ScaleRadius =  darkMatterProfile%scale            ()
    return
  end function zhao1996ScaleRadius
  
  double precision function zhao1996Normalization(self,node)
    !!{
    Returns the normalization of the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (darkMatterProfileDMOZhao1996), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    class           (nodeComponentBasic          ), pointer       :: basic
    double precision                                              :: radiusScale, radiusVirialScaleFree

    basic                 =>  node %basic                                (                          )
    radiusScale           =   self                      %scaleRadius     (node                      )
    radiusVirialScaleFree =  +self %darkMatterHaloScale_%radiusVirial    (node                      )    &
         &                   /                           radiusScale
    zhao1996Normalization =  +basic                     %mass            (                          )    &
         &                   /self                      %massUnnormalized(node,radiusVirialScaleFree)    &
         &                   /                           radiusScale                                 **3
    return
  end function zhao1996Normalization

  double precision function zhao1996MassUnnormalized(self,node,radiusScaleFree)
    !!{
    Returns the unnormalized mass in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radiusScaleFree}.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    implicit none
    class           (darkMatterProfileDMOZhao1996), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radiusScaleFree
    double precision                              , parameter     :: radiusScaleFreeTiny=1.0d-3
    double precision                                              :: alpha                     , beta, &
         &                                                           gamma

    select case (self%specialCase%ID)
    case (specialCaseGeneral%ID)
       call self%exponents(node,alpha,beta,gamma)
       zhao1996MassUnnormalized=+4.0d0                                                                                                            &
            &                   *Pi                                                                                                               &
            &                   *radiusScaleFree**(3.0d0-gamma)                                                                                   &
            &                   *Hypergeometric_2F1([(3.0d0-gamma)/alpha,(beta-gamma)/alpha],[1.0d0+(3.0d0-gamma)/alpha],-radiusScaleFree**alpha) &
            &                   /                    (3.0d0-gamma)
    case (specialCaseNFW%ID)
       if (radiusScaleFree <radiusScaleFreeTiny) then
          ! Use series solution for small radii.
          zhao1996MassUnnormalized=+ 2.0d0      *Pi*radiusScaleFree**2 &
               &                   - 8.0d0/3.0d0*Pi*radiusScaleFree**3 &
               &                   + 3.0d0      *Pi*radiusScaleFree**4 &
               &                   -16.0d0/5.0d0*Pi*radiusScaleFree**5 &
               &                   +10.0d0/3.0d0*Pi*radiusScaleFree**6
       else
          ! Use full solution.
          zhao1996MassUnnormalized=+4.0d0                         &
               &                   *Pi                            &
               &                   *(                             &
               &                     +log(+1.0d0+radiusScaleFree) &
               &                     -           radiusScaleFree  &
               &                     /   (+1.0d0+radiusScaleFree) &
               &                   )
       end if
    case (specialCaseGamma0_5NFW%ID)
       if (radiusScaleFree <radiusScaleFreeTiny) then
          ! Use series solution for small radii.
          zhao1996MassUnnormalized=+  8.0d0/ 5.0d0*Pi*radiusScaleFree**2.5d0 &
               &                   - 20.0d0/ 7.0d0*Pi*radiusScaleFree**3.5d0 &
               &                   + 35.0d0/ 9.0d0*Pi*radiusScaleFree**4.5d0 &
               &                   -105.0d0/22.0d0*Pi*radiusScaleFree**5.5d0
       else
          ! Use full solution.
          zhao1996MassUnnormalized=+(                                      &
               &                     -8.0d0                                &
               &                     *Pi                                   &
               &                     *sqrt(radiusScaleFree)                &
               &                     *(3.0d0+4.0d0*radiusScaleFree)        &
               &                    )                                      &
               &                   /(                                      &
               &                     +3.0d0                                &
               &                     *(1.0d0+      radiusScaleFree)**1.5d0 &
               &                     )                                     &
               &                   +8.0d0                                  &
               &                   *Pi                                     &
               &                   *asinh(sqrt(radiusScaleFree))
       end if
    case (specialCaseGamma1_5NFW%ID)
       if (radiusScaleFree <radiusScaleFreeTiny) then
          ! Use series solution for small radii.
          zhao1996MassUnnormalized=+  8.0d0/  3.0d0*Pi*radiusScaleFree**1.5d0 &
               &                   - 12.0d0/  5.0d0*Pi*radiusScaleFree**2.5d0 &
               &                   + 15.0d0/  7.0d0*Pi*radiusScaleFree**3.5d0 &
               &                   - 35.0d0/ 18.0d0*Pi*radiusScaleFree**4.5d0 &
               &                   +315.0d0/176.0d0*Pi*radiusScaleFree**5.5d0
       else
          ! Use full solution.
          zhao1996MassUnnormalized=+8.0d0                                                  &
               &                   *Pi                                                     &
               &                   *(                                                      &
               &                     -      sqrt(radiusScaleFree/(1.0d0+radiusScaleFree))  &
               &                     +asinh(sqrt(radiusScaleFree                        )) &
               &                    )
       end if
    case (specialCaseCoredNFW%ID)
       if (radiusScaleFree <radiusScaleFreeTiny) then
          ! Use series solution for small radii.
          zhao1996MassUnnormalized=+ 4.0d0/3.0d0*Pi*radiusScaleFree**3 &
               &                   - 3.0d0      *Pi*radiusScaleFree**4 &
               &                   +24.0d0/5.0d0*Pi*radiusScaleFree**5 &
               &                   -20.0d0/3.0d0*Pi*radiusScaleFree**6
       else
          ! Use full solution.
          zhao1996MassUnnormalized=+4.0d0                                  &
               &                   *Pi                                     &
               &                   *(                                      &
               &                     +log(+1.0d0+      radiusScaleFree)    &
               &                     -                 radiusScaleFree     &
               &                     *   (+2.0d0+3.0d0*radiusScaleFree)    &
               &                     /     2.0d0                           &
               &                     /   (+1.0d0+      radiusScaleFree)**2 &
               &                    )
       end if
    case default
       zhao1996MassUnnormalized=+0.0d0
       call Error_Report('unknown special case'//{introspection:location})
    end select
    return
  end function zhao1996MassUnnormalized

  double precision function zhao1996Density(self,node,radius)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (darkMatterProfileDMOZhao1996), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius
    double precision                                              :: alpha          , beta       , &
         &                                                           gamma          , radiusScale, &
         &                                                           radiusScaleFree

    radiusScale    = self%scaleRadius(node)
    radiusScaleFree=+     radius            &
         &          /     radiusScale
    select case (self%specialCase%ID)
    case (specialCaseGeneral%ID)
       call self%exponents(node,alpha,beta,gamma)
       zhao1996Density=+self%normalization(node) &
            &          /  radiusScaleFree**gamma &
            &          /(                        &
            &            +1.0d0                  &
            &            +radiusScaleFree**alpha &
            &           )**((beta-gamma)/alpha)
    case (specialCaseNFW%ID)
       zhao1996Density=+self%normalization(node) &
            &          /  radiusScaleFree        &
            &          /(                        &
            &            +1.0d0                  &
            &            +radiusScaleFree        &
            &           )**2
    case (specialCaseCoredNFW%ID)
       zhao1996Density=+self%normalization(node) &
            &          /(                        &
            &            +1.0d0                  &
            &            +radiusScaleFree        &
            &           )**3
    case (specialCaseGamma0_5NFW%ID)
       zhao1996Density=+self%normalization(node) &
            &          /sqrt(radiusScaleFree)    &
            &          /(                        &
            &            +1.0d0                  &
            &            +radiusScaleFree        &
            &           )**2.5d0
    case (specialCaseGamma1_5NFW%ID)
       zhao1996Density=+self%normalization(node) &
            &          /  radiusScaleFree**1.5d0 &
            &          /(                        &
            &            +1.0d0                  &
            &            +radiusScaleFree        &
            &           )**1.5d0
    case default
       zhao1996Density=+0.0d0
       call Error_Report('unknown special case'//{introspection:location})
    end select
    return
  end function zhao1996Density

  double precision function zhao1996DensityLogSlope(self,node,radius)
    !!{
    Returns the logarithmic slope of the density profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily
    radius} (given in units of Mpc).
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (darkMatterProfileDMOZhao1996), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius
    double precision                                              :: alpha          , beta       , &
         &                                                           gamma          , radiusScale, &
         &                                                           radiusScaleFree

    radiusScale    = self%scaleRadius(node)
    radiusScaleFree=+     radius            &
         &          /     radiusScale
    select case (self%specialCase%ID)
    case (specialCaseGeneral%ID)
       call self%exponents(node,alpha,beta,gamma)
       zhao1996DensityLogSlope=-(+gamma+beta*radiusScaleFree**alpha) &
            &                  /(+1.0d0+     radiusScaleFree**alpha)
    case (specialCaseNFW%ID)
       zhao1996DensityLogSlope=-(+1.0d0+3.0d0*radiusScaleFree) &
            &                  /(+1.0d0+      radiusScaleFree)
    case (specialCaseCoredNFW%ID)
       zhao1996DensityLogSlope=-(      +3.0d0*radiusScaleFree) &
            &                  /(+1.0d0+      radiusScaleFree)
    case (specialCaseGamma0_5NFW%ID)
       zhao1996DensityLogSlope=-(+0.5d0+3.0d0*radiusScaleFree) &
            &                  /(+1.0d0+      radiusScaleFree)
    case (specialCaseGamma1_5NFW%ID)
       zhao1996DensityLogSlope=-(+1.5d0+3.0d0*radiusScaleFree) &
            &                  /(+1.0d0+      radiusScaleFree)
    case default
       zhao1996DensityLogSlope=+0.0d0
       call Error_Report('unknown special case'//{introspection:location})
    end select
    return
  end function zhao1996DensityLogSlope

  double precision function zhao1996EnclosedMass(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOZhao1996), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius
    double precision                                              :: radiusScale, radiusScaleFree

    radiusScale         = self%scaleRadius     (node)
    radiusScaleFree     =+     radius                                    &
         &               /     radiusScale
    zhao1996EnclosedMass=+self%normalization   (node                )    &
         &               *self%scaleRadius     (node                )**3 &
         &               *self%massUnnormalized(node,radiusScaleFree)
    return
  end function zhao1996EnclosedMass

  double precision function zhao1996Potential(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    use :: Galactic_Structure_Options      , only : enumerationStructureErrorCodeType, structureErrorCodeSuccess
    use :: Galacticus_Nodes                , only : nodeComponentBasic
    use :: Gamma_Functions                 , only : Gamma_Function
    use :: Numerical_Constants_Math        , only : Pi
    use :: Hypergeometric_Functions        , only : Hypergeometric_pFq_Regularized
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOZhao1996     ), intent(inout)           :: self
    type            (treeNode                         ), intent(inout), target   :: node
    double precision                                   , intent(in   )           :: radius
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    class           (nodeComponentBasic               )               , pointer  :: basic
    double precision                                                             :: alpha          , beta       , &
         &                                                                          gamma          , radiusScale, &
         &                                                                          radiusScaleFree

    if (present(status)) status=structureErrorCodeSuccess
    call self%exponents(node,alpha,beta,gamma)
    basic                 =>  node%basic       (    )
    radiusScale           =   self%scaleRadius (node)
    radiusScaleFree       =  +     radius             &
         &                   /     radiusScale
    zhao1996Potential     =  +4.0d0                                                                                                                                                                    &
         &                   *Pi                                                                                                                                                                       &
         &                   *radiusScaleFree**(2.0d0      -gamma)                                                                                                                                     &
         &                   /                 (3.0d0      -gamma)                                                                                                                                     &
         &                   *  Gamma_Function((3.0d0+alpha-gamma)/alpha)                                                                                                                              &
         &                   /  Gamma_Function((3.0d0      -gamma)/alpha)                                                                                                                              &
         &                   *(                                                                                                                                                                        &
         &                     +Gamma_Function((2.0d0      -gamma)/alpha)*Hypergeometric_pFq_Regularized([(2.0d0-gamma)/alpha,(beta-gamma)/alpha],[(2.0d0+alpha-gamma)/alpha],-radiusScaleFree**alpha) &
         &                     -Gamma_Function((3.0d0      -gamma)/alpha)*Hypergeometric_pFq_Regularized([(3.0d0-gamma)/alpha,(beta-gamma)/alpha],[(3.0d0+alpha-gamma)/alpha],-radiusScaleFree**alpha) &
         &                    )                                                                                                                                                                        &
         &                   *gravitationalConstantGalacticus                                                                                                                                          &
         &                   *self%normalization(node)                                                                                                                                                 &
         &                   *self %scaleRadius (node)**2
    return
  end function zhao1996Potential

  double precision function zhao1996CircularVelocity(self,node,radius)
    !!{
    Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileDMOZhao1996), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius

    if (radius > 0.0d0) then
       zhao1996CircularVelocity=sqrt(gravitationalConstantGalacticus*self%enclosedMass(node,radius)/radius)
    else
       zhao1996CircularVelocity=0.0d0
    end if
    return
  end function zhao1996CircularVelocity

  double precision function zhao1996CircularVelocityMaximum(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMOZhao1996  ), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node

    zhao1996CircularVelocityMaximum=self%circularVelocityMaximumNumerical(node)
    return
  end function zhao1996CircularVelocityMaximum

  double precision function zhao1996RadiusCircularVelocityMaximum(self,node)
    !!{
    Returns the radius (in Mpc) at which the maximum circular velocity occurs in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMOZhao1996  ), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node

    zhao1996RadiusCircularVelocityMaximum=self%radiusCircularVelocityMaximumNumerical(node)
    return
  end function zhao1996RadiusCircularVelocityMaximum

  double precision function zhao1996RadialVelocityDispersion(self,node,radius)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Polylogarithms                  , only : Polylogarithm_2
    implicit none
    class           (darkMatterProfileDMOZhao1996  ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    double precision                                , parameter     :: radiusScaleFreeTiny               =1.0d-3
    double precision                                                :: radiusScale                              , radiusScaleFree, &
         &                                                             velocityDispersionSquaredScaleFree
    
    select case (self%specialCase%ID)
    case (specialCaseGeneral%ID)
       zhao1996RadialVelocityDispersion=self%radialVelocityDispersionNumerical(node,radius)
    case (specialCaseNFW%ID)       
       radiusScale    = self%scaleRadius(node)
       radiusScaleFree=+     radius            &
            &          /     radiusScale
       if (radiusScaleFree < radiusScaleFreeTiny) then
          ! Use series solution for small radii.
          velocityDispersionSquaredScaleFree=+11.0d0*Pi/15.0d0*radiusScaleFree**4                                                     &
               &                             +       Pi/ 6.0d0*radiusScaleFree**3*(-101.0d0+12.0d0*Pi**2-12.0d0*log(radiusScaleFree)) &
               &                             + 2.0d0*Pi/ 3.0d0*radiusScaleFree**2*(- 59.0d0+ 6.0d0*Pi**2- 6.0d0*log(radiusScaleFree)) &
               &                             +       Pi       *radiusScaleFree   *(- 23.0d0+ 2.0d0*Pi**2- 2.0d0*log(radiusScaleFree))
       else
          ! Use full solution.
          velocityDispersionSquaredScaleFree=+(                                                                                                 &
               &                               +2.0d0                                                                                           &
               &                               *Pi                                                                                              &
               &                               *(                                                                                               &
               &                                 +radiusScaleFree                                                                               &
               &                                 *(                                                                                             &
               &                                   -1.0d0                                                                                       &
               &                                   +radiusScaleFree                                                                             &
               &                                   *(                                                                                           &
               &                                     -        9.0d0                                                                             &
               &                                     -        7.0d0*radiusScaleFree                                                             &
               &                                     +Pi**2*(+1.0d0+radiusScaleFree)**2                                                         &
               &                                    )                                                                                           &
               &                                  )                                                                                             &
               &                                 +radiusScaleFree**4*log(+1.0d0+1.0d0/radiusScaleFree)                                          &
               &                                 +                   log(+1.0d0+      radiusScaleFree)                                          &
               &                                 +radiusScaleFree                                                                               &
               &                                 *(                                                                                             &
               &                                   -(                                                                                           &
               &                                     +            radiusScaleFree*(+1.0d0+2.0d0*radiusScaleFree)   *log(      +radiusScaleFree) &
               &                                    )                                                                                           &
               &                                   +                                                                log(+1.0d0+radiusScaleFree) &
               &                                   *(                                                                                           &
               &                                     -2.0d0-4.0d0*radiusScaleFree*(+2.0d0+      radiusScaleFree)                                &
               &                                     +3.0d0      *radiusScaleFree*(+1.0d0+      radiusScaleFree)**2*log(+1.0d0+radiusScaleFree) &
               &                                    )                                                                                           &
               &                                  )                                                                                             &
               &                                 +                 6.0d0                                                                        &
               &                                 *                       radiusScaleFree **2                                                    &
               &                                 *               (+1.0d0+radiusScaleFree)**2                                                    &
               &                                 *PolyLogarithm_2(      -radiusScaleFree)                                                       &
               &                                )                                                                                               &
               &                              )                                                                                                 &
               &                             /radiusScaleFree
       end if
       zhao1996RadialVelocityDispersion=+sqrt(                                                  &
            &                                 +     gravitationalConstantGalacticus             &
            &                                 *self%normalization                     (node)    &
            &                                 *     radiusScale                             **2 &
            &                                 *     velocityDispersionSquaredScaleFree          &
            &                                )
    case (specialCaseCoredNFW%ID)       
       radiusScale    = self%scaleRadius(node)
       radiusScaleFree=+     radius            &
            &          /     radiusScale
       if (radiusScaleFree < radiusScaleFreeTiny) then
          ! Use series solution for small radii.
          velocityDispersionSquaredScaleFree=+(119.0d0*Pi-12.0d0*Pi**3)/ 6.0d0                    &
               &                             +(119.0d0*Pi-12.0d0*Pi**3)/ 2.0d0*radiusScaleFree    &
               &                             +(353.0d0*Pi-36.0d0*Pi**3)/ 6.0d0*radiusScaleFree**2 &
               &                             +(121.0d0*Pi-12.0d0*Pi**3)/ 6.0d0*radiusScaleFree**3 &
               &                                         - 9.0d0*Pi**4 /20.0d0*radiusScaleFree**4
       else
          ! Use full solution.
          velocityDispersionSquaredScaleFree=+(                                                                                    &
               &                               +(                                                                                  &
               &                                 +Pi                                                                               &
               &                                 *(                                                                                &
               &                                   +95.0d0                                                                         &
               &                                   -12.0d0*Pi**2                                                                   &
               &                                   *(1.0d0+      radiusScaleFree)**4                                               &
               &                                   +       2.0d0*radiusScaleFree                                                   &
               &                                   *(                                                                              &
               &                                     +  130.0d0                                                                    &
               &                                     +    9.0d0*radiusScaleFree                                                    &
               &                                     *(                                                                            &
               &                                       + 13.0d0                                                                    &
               &                                       +  4.0d0*radiusScaleFree                                                    &
               &                                      )                                                                            &
               &                                    )                                                                              &
               &                                  )                                                                                &
               &                                )                                                                                  &
               &                               /6.0d0                                                                              &
               &                               /                         (1.0d0+      radiusScaleFree                         )**4 &
               &                               +2.0d0*Pi*log             (1.0d0+      radiusScaleFree                         )    &
               &                               *                         (2.0d0+9.0d0*radiusScaleFree+6.0d0*radiusScaleFree**2)    &
               &                               /                                      radiusScaleFree                              &
               &                               /                         (1.0d0+      radiusScaleFree                         )**2 &
               &                               - 6.0d0*Pi*log            (1.0d0+      radiusScaleFree                         )**2 &
               &                               -12.0d0*Pi*Polylogarithm_2(     -      radiusScaleFree                         )    &
               &                              )                                                                                    &
               &                             *                           (1.0d0+      radiusScaleFree                         )**3
       end if
       zhao1996RadialVelocityDispersion=+sqrt(                                                  &
            &                                 +     gravitationalConstantGalacticus             &
            &                                 *self%normalization                     (node)    &
            &                                 *     radiusScale                             **2 &
            &                                 *     velocityDispersionSquaredScaleFree          &
            &                                )
    case (specialCaseGamma0_5NFW%ID)       
       radiusScale    = self%scaleRadius(node)
       radiusScaleFree=+     radius            &
            &          /     radiusScale
       if (radiusScaleFree < radiusScaleFreeTiny) then
          ! Use series solution for small radii.
          velocityDispersionSquaredScaleFree=+16.0d0*Pi/  3.0d0*sqrt(radiusScaleFree      )*(- 11.0d0+  16.0d0*log(2.0d0)) &
               &                             +16.0d0*Pi/ 15.0d0*     radiusScaleFree**1.5d0*(-139.0d0+ 200.0d0*log(2.0d0)) &
               &                             + 2.0d0*Pi/  7.0d0*     radiusScaleFree**2.5d0*(-387.0d0+ 560.0d0*log(2.0d0)) &
               &                             + 4.0d0*Pi/189.0d0*     radiusScaleFree**3.5d0*(-887.0d0+1260.0d0*log(2.0d0))
       else
          ! Use full solution.
          velocityDispersionSquaredScaleFree=+(                                                             &
               &                               +8.0d0                                                       &
               &                               *Pi                                                          &
               &                               *(                                                           &
               &                                 - 6.0d0*sqrt(radiusScaleFree   *(+1.0d0+radiusScaleFree))  &
               &                                 -38.0d0*sqrt(radiusScaleFree**3*(+1.0d0+radiusScaleFree))  &
               &                                 -57.0d0*sqrt(radiusScaleFree**5*(+1.0d0+radiusScaleFree))  &
               &                                 -24.0d0*sqrt(radiusScaleFree**7*(+1.0d0+radiusScaleFree))  &
               &                                 -6.0d0                                                     &
               &                                 *(+1.0d0+      radiusScaleFree                        )**2 &
               &                                 *(+1.0d0+2.0d0*radiusScaleFree                        )    &
               &                                 *(-1.0d0+8.0d0*radiusScaleFree*(1.0d0+radiusScaleFree))    &
               &                                 *asinh(sqrt(radiusScaleFree))                              &
               &                                 -24.0d0                                                    &
               &                                 *        radiusScaleFree **1.5d0                           &
               &                                 *(+1.0d0+radiusScaleFree)**3.5d0                           &
               &                                 *log(radiusScaleFree/(16.0d0*(1.0d0+radiusScaleFree)))     &
               &                                 +24.0d0                                                    &
               &                                 *(                                                         &
               &                                   +      sqrt(radiusScaleFree**3*(1.0d0+radiusScaleFree))  &
               &                                   +3.0d0*sqrt(radiusScaleFree**5*(1.0d0+radiusScaleFree))  &
               &                                   +3.0d0*sqrt(radiusScaleFree**7*(1.0d0+radiusScaleFree))  &
               &                                   +      sqrt(radiusScaleFree**9*(1.0d0+radiusScaleFree))  &
               &                                  )                                                         &
               &                                 *(                                                         &
               &                                   +log(       radiusScaleFree)                             &
               &                                   +log(+1.0d0+radiusScaleFree)                             &
               &                                  )                                                         &
               &                                )                                                           &
               &                              )                                                             &
               &                             /(                                                             &
               &                               +9.0d0                                                       &
               &                               *        radiusScaleFree                                     &
               &                               *(+1.0d0+radiusScaleFree)                                    &
               &                              )
       end if
       zhao1996RadialVelocityDispersion=+sqrt(                                                  &
            &                                 +     gravitationalConstantGalacticus             &
            &                                 *self%normalization                     (node)    &
            &                                 *     radiusScale                             **2 &
            &                                 *     velocityDispersionSquaredScaleFree          &
            &                                )
    case (specialCaseGamma1_5NFW%ID)       
       radiusScale    = self%scaleRadius(node)
       radiusScaleFree=+     radius            &
            &          /     radiusScale
       if (radiusScaleFree < radiusScaleFreeTiny) then
          ! Use series solution for small radii.
          velocityDispersionSquaredScaleFree=+8.0d0*Pi/   3.0d0*sqrt(radiusScaleFree       )                                                                 &
               &                             -      Pi/3150.0d0*     radiusScaleFree**3.5d0 *(-19861.0d0+60480.0d0*log(2.0d0)-7560.0d0*log(radiusScaleFree)) &
               &                             -      Pi/ 175.0d0*     radiusScaleFree**2.5d0 *(- 8683.0d0+13440.0d0*log(2.0d0)-1680.0d0*log(radiusScaleFree)) &
               &                             -4.0d0*Pi/  75.0d0*     radiusScaleFree**1.5d0 *(-  817.0d0+  960.0d0*log(2.0d0)- 120.0d0*log(radiusScaleFree))
       else
          ! Use full solution.
          velocityDispersionSquaredScaleFree=+(                                                                     &
               &                               +8.0d0                                                               &
               &                               *Pi                                                                  &
               &                               *        radiusScaleFree **1.5d0                                     &
               &                               *(+1.0d0+radiusScaleFree)**1.5d0                                     &
               &                               *(                                                                   &
               &                                 +(                                                                 &
               &                                   +    2.0d0                                                       &
               &                                   *(                                                               &
               &                                     +  1.0d0                                                       &
               &                                     +        2.0d0*radiusScaleFree                                 &
               &                                     *(-1.0d0+4.0d0*radiusScaleFree+8.0d0*radiusScaleFree**2)       &
               &                                    )                                                               &
               &                                   *asinh(sqrt(radiusScaleFree))                                    &
               &                                  )                                                                 &
               &                                 /sqrt(radiusScaleFree**5*(1.0d0+radiusScaleFree))                  &
               &                                 +(                                                                 &
               &                                   -2.0d0                                                           &
               &                                   +radiusScaleFree                                                 &
               &                                   *(                                                               &
               &                                     + 5.0d0                                                        &
               &                                     +12.0d0*radiusScaleFree                                        &
               &                                     -32.0d0*radiusScaleFree   *(+1.0d0+radiusScaleFree)*log(2.0d0) &
               &                                    )                                                               &
               &                                   +   4.0d0*radiusScaleFree**2*(+1.0d0+radiusScaleFree)            &
               &                                   *(                                                               &
               &                                     +      log(       radiusScaleFree)                             &
               &                                     -5.0d0*log(+1.0d0+radiusScaleFree)                             &
               &                                    )                                                               &
               &                                  )                                                                 &
               &                                 /(                                                                 &
               &                                   +         radiusScaleFree**2*(+1.0d0+radiusScaleFree)            &
               &                                  )                                                                 &
               &                                )                                                                   &
               &                              )                                                                     &
               &                             /5.0d0
       end if
       zhao1996RadialVelocityDispersion=+sqrt(                                                  &
            &                                 +     gravitationalConstantGalacticus             &
            &                                 *self%normalization                     (node)    &
            &                                 *     radiusScale                             **2 &
            &                                 *     velocityDispersionSquaredScaleFree          &
            &                                )
    case default
       zhao1996RadialVelocityDispersion=+0.0d0
       call Error_Report('unknown special case'//{introspection:location})
    end select
    return
  end function zhao1996RadialVelocityDispersion

  double precision function zhao1996RadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !!{
    Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily
    specificAngularMomentum} (given in units of km s$^{-1}$ Mpc)
    !!}
    implicit none
    class           (darkMatterProfileDMOZhao1996), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: specificAngularMomentum

    zhao1996RadiusFromSpecificAngularMomentum=self%radiusFromSpecificAngularMomentumNumerical(node,specificAngularMomentum)
    return
  end function zhao1996RadiusFromSpecificAngularMomentum

  double precision function zhao1996RotationNormalization(self,node)
    !!{
    Return the normalization of the rotation velocity vs. specific angular momentum relation.
    !!}
    implicit none
    class           (darkMatterProfileDMOZhao1996), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                                              :: radiusVirial
    
    radiusVirial                 =+self%darkMatterHaloScale_%radiusVirial(node                                        )
    zhao1996RotationNormalization=+self                     %radialMoment(node,moment=2.0d0,radiusMaximum=radiusVirial) &
         &                        /self                     %radialMoment(node,moment=3.0d0,radiusMaximum=radiusVirial)
    return
  end function zhao1996RotationNormalization

  double precision function zhao1996Energy(self,node)
    !!{
    Return the energy of an Zhao1996 halo density profile.
    !!}
    implicit none
    class(darkMatterProfileDMOZhao1996), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    zhao1996Energy=self%energyNumerical(node)
    return
  end function zhao1996Energy

  double precision function zhao1996KSpace(self,node,waveNumber)
    !!{
    Returns the Fourier transform of the Zhao1996 density profile at the specified {\normalfont \ttfamily waveNumber} (given in Mpc$^{-1}$), using the
    expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    !!}
    implicit none
    class           (darkMatterProfileDMOZhao1996), intent(inout)          :: self
    type            (treeNode                    ), intent(inout), target  :: node
    double precision                              , intent(in   )          :: waveNumber

    zhao1996KSpace=self%kSpaceNumerical(node,wavenumber)
    return
  end function zhao1996KSpace

  double precision function zhao1996FreefallRadius(self,node,time)
    !!{
    Returns the freefall radius in the Zhao1996 density profile at the specified {\normalfont \ttfamily time} (given in Gyr).
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileDMOZhao1996), intent(inout), target :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: time
    double precision                              , parameter             :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-3
    type            (rootFinder                  )                        :: finder

    self_  => self
    node_  => node
    time_  =  time
    finder =  rootFinder(                                                             &
         &               rootFunction                 =rootRadiusFreefall           , &
         &               toleranceAbsolute            =toleranceAbsolute            , &
         &               toleranceRelative            =toleranceRelative            , &
         &               rangeExpandDownward          =0.5d0                        , &
         &               rangeExpandUpward            =2.0d0                        , &
         &               rangeExpandType              =rangeExpandMultiplicative    , &
         &               rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
         &               rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
         &              )    
    zhao1996FreefallRadius=finder%find(rootGuess=self%darkMatterHaloScale_%radiusVirial(node))
    return
  end function zhao1996FreefallRadius

  double precision function rootRadiusFreefall(radiusFreefall)
    !!{
    Root function used in finding the radius corresponding to a given freefall time.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    double precision            , intent(in   ) :: radiusFreefall
    type            (integrator)                :: integrator_

    radiusFreefall_   =+radiusFreefall
    integrator_       = integrator              (integrandTimeFreefall,toleranceRelative=1.0d-3)
    rootRadiusFreefall=+integrator_   %integrate(0.0d0                ,radiusFreefall          ) &
         &             -time_
    return
  end function rootRadiusFreefall

  double precision function integrandTimeFreefall(radius)
    !!{
    Integrand for freefall time in the halo.
    !!}
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: potentialDifference

    potentialDifference=-self_%potential(node_,radiusFreefall_) &
         &              +self_%potential(node_,radius         )
    if (potentialDifference < 0.0d0) then
       integrandTimeFreefall=+Mpc_per_km_per_s_To_Gyr   &
            &                /sqrt(                     &
            &                      -2.0d0               &
            &                      *potentialDifference &
            &                     )
    else
       ! Avoid floating point errors arising from rounding errors.
       integrandTimeFreefall=0.0d0
    end if
    return
  end function integrandTimeFreefall

  double precision function zhao1996FreefallRadiusIncreaseRate(self,node,time)
    !!{
    Returns the rate of increase of the freefall radius in the Zhao1996 density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    use :: Numerical_Differentiation, only : differentiator
    implicit none
    class           (darkMatterProfileDMOZhao1996), intent(inout), target :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: time
    double precision                              , parameter             :: timeLogarithmicStep=0.1d0
    type            (differentiator              )                        :: differentiator_

    self_                              =>  self
    node_                              =>  node
    differentiator_                    =   differentiator            (freefallRadiusEvaluate                    )
    zhao1996FreefallRadiusIncreaseRate =  +differentiator_%derivative(log(time)             ,timeLogarithmicStep) &
         &                                /                               time
    return
  end function zhao1996FreefallRadiusIncreaseRate

  double precision function freefallRadiusEvaluate(timeLogarithmic)
    !!{
    GSL-callable function to evaluate the freefall radius of the dark matter profile.
    !!}
    implicit none
    double precision, intent(in   ), value :: timeLogarithmic

    freefallRadiusEvaluate=self_%freefallRadiusNumerical(node_,exp(timeLogarithmic))
    return
  end function freefallRadiusEvaluate

  double precision function zhao1996RadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given
    in units of Mpc).
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (darkMatterProfileDMOZhao1996), intent(inout)           :: self
    type            (treeNode                    ), intent(inout)           :: node
    double precision                              , intent(in   )           :: moment
    double precision                              , intent(in   ), optional :: radiusMinimum    , radiusMaximum
    double precision                                                        :: radiusScale      , radiusScaleFree  , &
         &                                                                     radialMomentUpper, radialMomentLower, &
         &                                                                     alpha            , beta             , &
         &                                                                     gamma

    call self%exponents(node,alpha,beta,gamma)
    radiusScale           =   self                     %scaleRadius (node)
    if (present(radiusMinimum)) then
       radiusScaleFree  =+radiusMinimum &
            &            /radiusScale
       radialMomentLower= radialMomentIndefinite(radiusScaleFree)
    else
       radialMomentLower=0.0d0
       if (alpha <= 0.0d0 .or. 1.0d0+moment <= gamma) call Error_Report('radial moment is undefined'//{introspection:location})
    end if
    if (present(radiusMaximum)) then
       radiusScaleFree  =+radiusMaximum &
            &            /radiusScale
       radialMomentUpper= radialMomentIndefinite(radiusScaleFree)
    else
       radialMomentUpper=0.0d0
       call Error_Report('radial moment is not implemented'//{introspection:location})
    end if
    zhao1996RadialMoment=+(                           &
         &                 +radialMomentUpper         &
         &                 -radialMomentLower         &
         &                )                           &
         &               *self%normalization(node)    &
         &               *radiusScale**(moment+1.0d0)
    return

  contains

    double precision function radialMomentIndefinite(radiusScaleFree)
      !!{
      Compute the indefinite radial moment.
      !!}
      use :: Hypergeometric_Functions, only : Hypergeometric_2F1
      implicit none
      double precision, intent(in   ) :: radiusScaleFree

      radialMomentIndefinite=+radiusScaleFree**   (1.0d0+moment-gamma)                                                                                       &
           &                 /                    (1.0d0+moment-gamma)                                                                                       &
           &                 *Hypergeometric_2F1([(1.0d0+moment-gamma)/alpha,(beta-gamma)/alpha],[1.0d0+(1.0d0+moment-gamma)/alpha],-radiusScaleFree**alpha)
      return
    end function radialMomentIndefinite
    
  end function zhao1996RadialMoment

  double precision function zhao1996RadiusEnclosingDensity(self,node,density)
    !!{
    Compute the radius enclosing a given density for {\normalfont \ttfamily Zhao1996} dark matter halo profiles.
    !!}
    implicit none
    class           (darkMatterProfileDMOZhao1996), intent(inout), target :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: density

    zhao1996RadiusEnclosingDensity=self%radiusEnclosingDensityNumerical(node,density)
    return
  end function zhao1996RadiusEnclosingDensity

  double precision function zhao1996RadiusEnclosingMass(self,node,mass)
    !!{
    Compute the radius enclosing a given mass for {\normalfont \ttfamily Zhao1996} dark matter halo profiles.
    !!}
    implicit none
    class           (darkMatterProfileDMOZhao1996), intent(inout), target :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: mass

    zhao1996RadiusEnclosingMass=self%radiusEnclosingMassNumerical(node,mass)
    return
  end function zhao1996RadiusEnclosingMass
