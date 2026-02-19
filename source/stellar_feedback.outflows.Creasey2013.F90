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
  Implementation of the \cite{creasey_how_2013} stellar feedback model.
  !!}

  use :: Star_Formation_Rate_Surface_Density_Disks, only : starFormationRateSurfaceDensityDisksClass

  !![
  <stellarFeedbackOutflows name="stellarFeedbackOutflowsCreasey2013">
   <description>
    A stellar feedback outflow class which implements the model of \cite{creasey_how_2013}. Specifically, the outflow rate is:
    \begin{equation}
    \dot{M}_\mathrm{outflow} = {\dot{E}_\mathrm{SN} \over E_\mathrm{SN} \dot{M}_\star} \int_0^\infty \beta_0
    \Sigma_{g,1}^{-\mu}(r) f_\mathrm{g}^\nu(r) \dot{\Sigma}_\star(r) 2 \pi r \mathrm{d}r,
    \end{equation}
    where $\Sigma_{g,1}(r)$ is the surface density of gas in units of $M_\odot$ pc$^{-2}$, $f_\mathrm{g}(r)$ is the gas
    fraction, $\dot{\Sigma}_\star(r)$ is the surface density of star formation rate, $\dot{M}_\star$ is the total star
    formation rate in the disk, $\dot{E}_\mathrm{SN}$ is the current energy input rate from supernovae, $E_\mathrm{SN}$ is the
    total energy input per unit mass from a stellar population after infinite time, $\beta_0=${\normalfont \ttfamily [beta0]},
    $\mu=${\normalfont \ttfamily [mu]}, and $\nu=${\normalfont \ttfamily [nu]}.
   </description>
  </stellarFeedbackOutflows>
  !!]
  type, extends(stellarFeedbackOutflowsClass) :: stellarFeedbackOutflowsCreasey2013
     !!{
     Implementation of the \cite{creasey_how_2013} stellar feedback model.
     !!}
     private
     class           (starFormationRateSurfaceDensityDisksClass), pointer :: starFormationRateSurfaceDensityDisks_ => null()
     double precision                                                     :: nu                                             , mu, &
          &                                                                  beta0
   contains
     final     ::                creasey2013Destructor
     procedure :: outflowRate => creasey2013OutflowRate
  end type stellarFeedbackOutflowsCreasey2013

  interface stellarFeedbackOutflowsCreasey2013
     !!{
     Constructors for the creasey2013 fraction stellar feedback class.
     !!}
     module procedure creasey2013ConstructorParameters
     module procedure creasey2013ConstructorInternal
  end interface stellarFeedbackOutflowsCreasey2013

contains

  function creasey2013ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{creasey_how_2013} stellar feedback class which takes a parameter set as input.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (stellarFeedbackOutflowsCreasey2013       )                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (starFormationRateSurfaceDensityDisksClass), pointer       :: starFormationRateSurfaceDensityDisks_
    double precision                                                           :: mu                                   , nu, &
         &                                                                        beta0

    !![
    <inputParameter>
      <name>mu</name>
      <source>parameters</source>
      <defaultValue>1.15d0</defaultValue>
      <defaultSource>\citep{creasey_how_2013}</defaultSource>
      <description>The parameter $\mu$ appearing in the \cite{creasey_how_2013} model for supernovae feedback.</description>
    </inputParameter>
    <inputParameter>
      <name>nu</name>
      <source>parameters</source>
      <defaultValue>0.16d0</defaultValue>
      <defaultSource>\citep{creasey_how_2013}</defaultSource>
      <description>The parameter $\nu$ appearing in the \cite{creasey_how_2013} model for supernovae feedback.</description>
    </inputParameter>
    <inputParameter>
      <name>beta0</name>
      <source>parameters</source>
      <defaultValue>13.0d0</defaultValue>
      <defaultSource>\citep{creasey_how_2013}</defaultSource>
      <description>The parameter $\beta_0$ appearing in the \cite{creasey_how_2013} model for supernovae feedback.</description>
    </inputParameter>
    <objectBuilder class="starFormationRateSurfaceDensityDisks" name="starFormationRateSurfaceDensityDisks_" source="parameters"/>
    !!]
    self=stellarFeedbackOutflowsCreasey2013(mu,nu,beta0,starFormationRateSurfaceDensityDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateSurfaceDensityDisks_"/>
    !!]
    return
  end function creasey2013ConstructorParameters

  function creasey2013ConstructorInternal(mu,nu,beta0,starFormationRateSurfaceDensityDisks_) result(self)
    !!{
    Internal constructor for the \refClass{stellarFeedbackOutflowsCreasey2013} stellar feedback class.
    !!}
    implicit none
    type            (stellarFeedbackOutflowsCreasey2013       )                        :: self
    class           (starFormationRateSurfaceDensityDisksClass), intent(in   ), target :: starFormationRateSurfaceDensityDisks_
    double precision                                           , intent(in   )         :: mu                                   , nu, &
         &                                                                                beta0
    !![
    <constructorAssign variables="mu, nu, beta0, *starFormationRateSurfaceDensityDisks_"/>
    !!]

    return
  end function creasey2013ConstructorInternal

  subroutine creasey2013Destructor(self)
    !!{
    Destructor for the \refClass{stellarFeedbackOutflowsCreasey2013} feedback in disks class.
    !!}
    implicit none
    type(stellarFeedbackOutflowsCreasey2013), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateSurfaceDensityDisks_"/>
    !!]
    return
  end subroutine creasey2013Destructor

  subroutine creasey2013OutflowRate(self,component,rateStarFormation,rateEnergyInput,rateOutflowEjective,rateOutflowExpulsive)
    !!{
    Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic disk {\normalfont \ttfamily component} using
    the model of \cite{creasey_how_2013}. The outflow rate is given by
    \begin{equation}
    \dot{M}_\mathrm{outflow} = \int_0^\infty \beta_0 \Sigma_{g,1}^{-\mu}(r) f_\mathrm{g}^\nu(r) \dot{\Sigma}_\star(r) 2 \pi r \mathrm{d}r,
    \end{equation}
    where $\Sigma_{g,1}(r)$ is the surface density of gas in units of $M_\odot$ pc$^{-2}$, $f_\mathrm{g}(r)$ is the gas
    fraction, $\dot{\Sigma}_\star(r)$ is the surface density of star formation rate, $\beta_0=${\normalfont \ttfamily [beta0]},
    $\mu=${\normalfont \ttfamily [mu]}, and $\nu=${\normalfont \ttfamily [nu]}.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk                     , coordinateSystemCylindrical, massTypeGaseous, massTypeStellar
    use :: Galacticus_Nodes          , only : nodeComponentDisk                     , nodeComponentSpheroid
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Numerical_Constants_Math  , only : Pi
    use :: Numerical_Integration     , only : integrator
    use :: Stellar_Feedback          , only : feedbackEnergyInputAtInfinityCanonical
    implicit none
    class           (stellarFeedbackOutflowsCreasey2013), intent(inout) :: self
    class           (nodeComponent                     ), intent(inout) :: component
    double precision                                    , intent(in   ) :: rateEnergyInput               , rateStarFormation
    double precision                                    , intent(  out) :: rateOutflowEjective           , rateOutflowExpulsive
    class           (massDistributionClass             ), pointer       :: massDistributionGaseous       , massDistributionStellar
    double precision                                    , parameter     :: radiusInnerDimensionless=0.0d0, radiusOuterDimensionless=10.0d0
    double precision                                                    :: radiusScale                   , massGas                        , &
         &                                                                 radiusInner                   , radiusOuter                    , &
         &                                                                 massStellar
    type            (integrator                        )                :: integrator_

    ! Nothing to do if no energy input.
    if (rateEnergyInput <= 0.0d0) then
       rateOutflowEjective =0.0d0
       rateOutflowExpulsive=0.0d0
       return
    end if
    ! Get the disk properties.
    select type (component)
    class is (nodeComponentDisk    )
       massGas    =component%massGas    ()
       massStellar=component%massStellar()
       radiusScale=component%radius     ()
    class default
       massGas    =0.0d0
       massStellar=0.0d0
       radiusScale=0.0d0
       call Error_Report('unsupported component'//{introspection:location})
    end select
    ! Return immediately for a null component.
    if (massGas <= 0.0d0 .or. massStellar <= 0.0d0 .or. radiusScale <= 0.0d0) then
       rateOutflowEjective =+0.0d0
       rateOutflowExpulsive=+0.0d0
       return
    end if
    ! Compute suitable limits for the integration.
    radiusInner=radiusScale*radiusInnerDimensionless
    radiusOuter=radiusScale*radiusOuterDimensionless
    ! Compute the outflow rate.
    massDistributionGaseous => component%hostNode%massDistribution(componentType=componentTypeDisk,massType=massTypeGaseous)
    massDistributionStellar => component%hostNode%massDistribution(componentType=componentTypeDisk,massType=massTypeStellar)
    integrator_             =  integrator(outflowRateIntegrand,toleranceRelative=1.0d-3)
    rateOutflowEjective     =  +2.0d0                                          &
         &                     *Pi                                             &
         &                     *self%beta0                                     &
         &                     *integrator_%integrate(radiusInner,radiusOuter) &
         &                     /rateStarFormation                              &
         &                     *rateEnergyInput                                &
         &                     /feedbackEnergyInputAtInfinityCanonical
    rateOutflowExpulsive    =  +0.0d0
    !![
    <objectDestructor name="massDistributionGaseous"/>
    <objectDestructor name="massDistributionStellar"/>
    !!]    
    return

  contains

    double precision function outflowRateIntegrand(radius)
      !!{
      Integrand function for the ``Creasey et al. (2012)'' supernovae feedback calculation.
      !!}
      use :: Coordinates                 , only : coordinateCylindrical, assignment(=)
      use :: Numerical_Constants_Prefixes, only : mega
      implicit none
      double precision                       , intent(in   ) :: radius
      double precision                                       :: fractionGas      , densitySurfaceRateStarFormation, &
           &                                                    densitySurfaceGas, densitySurfaceStellar
      type            (coordinateCylindrical)                :: coordinates

      coordinates=[radius,0.0d0,0.0d0]
      ! Get gas surface density.
      densitySurfaceGas    =massDistributionGaseous%surfaceDensity(coordinates)
      ! Get stellar surface density.
      densitySurfaceStellar=massDistributionStellar%surfaceDensity(coordinates)
      ! Compute the gas fraction.
      fractionGas=+  densitySurfaceGas     &
           &      /(                       &
           &        +densitySurfaceGas     &
           &        +densitySurfaceStellar &
           &       )
      ! Convert gas surface density to units of M☉ pc⁻².
      densitySurfaceGas=+densitySurfaceGas &
           &            /mega**2
      ! Get the surface density of star formation rate.
      densitySurfaceRateStarFormation=self%starFormationRateSurfaceDensityDisks_%rate(component%hostNode,radius)
      ! Compute the outflow rate.
      outflowRateIntegrand=+densitySurfaceGas              **(-self%mu) &
           &               *fractionGas                    **  self%nu  &
           &               *densitySurfaceRateStarFormation             &
           &               *radius
      return
    end function outflowRateIntegrand

  end subroutine creasey2013OutflowRate

