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
  Implementation of a the Kennicutt-Schmidt star formation rate surface density for galactic disks.
  !!}

  use :: Kind_Numbers, only : kind_int8

  !![
  <starFormationRateSurfaceDensityDisks name="starFormationRateSurfaceDensityDisksKennicuttSchmidt">
   <description>
    A star formation rate surface density class which assumes that the Kennicutt-Schmidt law holds
    \citep{schmidt_rate_1959,kennicutt_global_1998}:
    \begin{equation}
    \dot{\Sigma}_\star = A \left({\Sigma_\mathrm{H} \over M_\odot \hbox{pc}^{-2}} \right)^N,
    \end{equation}
    where $A=${\normalfont \ttfamily [normalization]} and $N=${\normalfont \ttfamily [exponent]} are parameters. Optionally, if
    the {\normalfont \ttfamily [truncate]} parameter is set to true, then the star formation rate is truncated below a critical
    surface density such that
    \begin{equation}
    \dot{\Sigma}_\star = \left\{ \begin{array}{ll} A \left({\Sigma_\mathrm{H} \over M_\odot \hbox{pc}^{-2}} \right)^N &amp;
    \hbox{ if } \Sigma_\mathrm{gas,disk} &gt; \Sigma_\mathrm{crit} \\ A \left({\Sigma_\mathrm{H} \over M_\odot \hbox{pc}^{-2}}
    \right)^N \left(\Sigma_\mathrm{gas,disk}/\Sigma_\mathrm{crit}\right)^\alpha &amp; \hbox{ otherwise.} \end{array} \right.
    \end{equation}
    Here, $\alpha=${\normalfont \ttfamily [exponentTruncated]} and $\Sigma_\mathrm{crit}$ is a critical surface density for
    star formation which we specify as
    \begin{equation}
    \Sigma_\mathrm{crit} = {q_\mathrm{crit} \kappa \sigma_\mathrm{gas} \over \pi \G},
    \end{equation}
    where $\kappa$ is the epicyclic frequency in the disk, $\sigma_\mathrm{gas}$ is the velocity dispersion of gas in the disk
    and $q_\mathrm{crit}=${\normalfont \ttfamily [toomreParameterCritical]} is a dimensionless constant of order unity which
    controls where the critical density occurs. We assume that $\sigma_\mathrm{gas}$ is a constant equal to {\normalfont
    \ttfamily [velocityDispersionDiskGas]} and that the disk has a flat rotation curve such that $\kappa = \sqrt{2} V/R$.
   </description>
  </starFormationRateSurfaceDensityDisks>
  !!]
  type, extends(starFormationRateSurfaceDensityDisksClass) :: starFormationRateSurfaceDensityDisksKennicuttSchmidt
     !!{
     Implementation of a Kennicutt-Schmidt star formation rate surface density for galactic disks.
     !!}
     private
     double precision            :: normalization               , exponent                 , &
          &                         exponentTruncated           , velocityDispersionDiskGas, &
          &                         toomreParameterCritical
     logical                     :: truncate
     integer         (kind_int8) :: lastUniqueID
     logical                     :: factorsComputed
     double precision            :: surfaceDensityCriticalFactor, hydrogenMassFraction
   contains
     !![
     <methods>
       <method description="Reset memoized calculations." method="calculationReset" />
     </methods>
     !!]
     final     ::                     kennicuttSchmidtDestructor
     procedure :: autoHook         => kennicuttSchmidtAutoHook
     procedure :: calculationReset => kennicuttSchmidtCalculationReset
     procedure :: rate             => kennicuttSchmidtRate
  end type starFormationRateSurfaceDensityDisksKennicuttSchmidt

  interface starFormationRateSurfaceDensityDisksKennicuttSchmidt
     !!{
     Constructors for the \refClass{starFormationRateSurfaceDensityDisksKennicuttSchmidt} star formation surface density rate in disks class.
     !!}
     module procedure kennicuttSchmidtConstructorParameters
     module procedure kennicuttSchmidtConstructorInternal
  end interface starFormationRateSurfaceDensityDisksKennicuttSchmidt

contains

  function kennicuttSchmidtConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationRateSurfaceDensityDisksKennicuttSchmidt} star formation surface density rate in disks class which takes a parameter set as input.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (starFormationRateSurfaceDensityDisksKennicuttSchmidt)                :: self
    type            (inputParameters                                     ), intent(inout) :: parameters
    double precision                                                                      :: normalization          , exponent                 , &
         &                                                                                   exponentTruncated      , velocityDispersionDiskGas, &
         &                                                                                   toomreParameterCritical
    logical                                                                               :: truncate

    !![
    <inputParameter>
      <name>normalization</name>
      <defaultSource>\citep{kennicutt_global_1998}</defaultSource>
      <defaultValue>0.147d0</defaultValue>
      <description>The normalization of the Kennicutt-Schmidt star formation law [$M_\odot$ Gyr$^{-1}$pc$^{-2}$].</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>exponent</name>
      <defaultSource>\citep{kennicutt_global_1998}</defaultSource>
      <defaultValue>1.400d0</defaultValue>
      <description>The exponent in the Kennicutt-Schmidt star formation law.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>truncate</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether or not to truncate star formation below a critical surface density in disks.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>exponentTruncated</name>
      <defaultValue>6.0d0</defaultValue>
      <description>The exponent of the $\Sigma_\mathrm{gas}/\Sigma_\mathrm{crit}$ term used in truncating the Kennicutt-Schmidt star formation law.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>velocityDispersionDiskGas</name>
      <defaultSource>\citep{leroy_star_2008}</defaultSource>
      <defaultValue>10.0d0</defaultValue>
      <description>The velocity dispersion of gas in disks.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>toomreParameterCritical</name>
      <defaultSource>\citep{kennicutt_star_1989}</defaultSource>
      <defaultValue>0.4d0</defaultValue>
      <description>The critical Toomre parameter for star formation in disks.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=starFormationRateSurfaceDensityDisksKennicuttSchmidt(normalization,exponent,truncate,exponentTruncated,velocityDispersionDiskGas,toomreParameterCritical)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function kennicuttSchmidtConstructorParameters

  function kennicuttSchmidtConstructorInternal(normalization,exponent,truncate,exponentTruncated,velocityDispersionDiskGas,toomreParameterCritical) result(self)
    !!{
    Internal constructor for the \refClass{starFormationRateSurfaceDensityDisksKennicuttSchmidt} star formation surface density rate from disks class.
    !!}
    use :: Numerical_Constants_Prefixes, only : mega
    implicit none
    type            (starFormationRateSurfaceDensityDisksKennicuttSchmidt)                :: self
    double precision                                                      , intent(in   ) :: normalization          , exponent                 , &
         &                                                                                   exponentTruncated      , velocityDispersionDiskGas, &
         &                                                                                   toomreParameterCritical
    logical                                                               , intent(in   ) :: truncate
    !![
    <constructorAssign variables="normalization, exponent, truncate, exponentTruncated, velocityDispersionDiskGas, toomreParameterCritical"/>
    !!]

    self%lastUniqueID   =-1_kind_int8
    self%factorsComputed=.false.
    ! Renormalize the Kennicutt-Schmidt relation to our internal units.
    self%normalization=+self%normalization                &
         &             *mega**(2.0d0-2.0d0*self%exponent)
    return
  end function kennicuttSchmidtConstructorInternal

  subroutine kennicuttSchmidtAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(starFormationRateSurfaceDensityDisksKennicuttSchmidt), intent(inout) :: self

    call calculationResetEvent%attach(self,kennicuttSchmidtCalculationReset,openMPThreadBindingAllLevels,label='starFormationRateSurfaceDensityDisksKennicutt\Schmidt')
    return
  end subroutine kennicuttSchmidtAutoHook

  subroutine kennicuttSchmidtDestructor(self)
    !!{
    Destructor for the kennicuttSchmidt cooling radius class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(starFormationRateSurfaceDensityDisksKennicuttSchmidt), intent(inout) :: self

    if (calculationResetEvent%isAttached(self,kennicuttSchmidtCalculationReset)) call calculationResetEvent%detach(self,kennicuttSchmidtCalculationReset)
    return
  end subroutine kennicuttSchmidtDestructor

  subroutine kennicuttSchmidtCalculationReset(self,node,uniqueID)
    !!{
    Reset the Kennicutt-Schmidt relation calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (starFormationRateSurfaceDensityDisksKennicuttSchmidt), intent(inout) :: self
    type   (treeNode                                            ), intent(inout) :: node
    integer(kind_int8                                           ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%factorsComputed=.false.
    self%lastUniqueID   =uniqueID
    return
  end subroutine kennicuttSchmidtCalculationReset

  double precision function kennicuttSchmidtRate(self,node,radius)
    !!{
    Returns the star formation rate surface density  (in $M_\odot$ Gyr$^{-1}$ Mpc$^{-2}$) for star formation in the galactic disk of {\normalfont \ttfamily node}. The disk is assumed to obey the Kennicutt-Schmidt law:
    \begin{equation}
    \Sigma_\star = A \left(x_\mathrm{H} {\Sigma_\mathrm{gas}\over M_\odot \hbox{pc}^{-2}}\right)^N,
    \end{equation}
    where $A=${\normalfont \ttfamily [normalization]} and $N=${\normalfont \ttfamily
    [exponent]}. Optionally, star formation is truncated for gas surface densities below a critical density of:
    \begin{equation}
    \Sigma_\mathrm{crit} = {q_\mathrm{crit} \kappa \sigma_\mathrm{gas} \over \pi \G},
    \end{equation}
    where $\kappa$ is the epicyclic frequency in the disk, $\sigma_\mathrm{gas}$ is the velocity dispersion of gas in the disk and
    $q_\mathrm{crit}=${\normalfont \ttfamily [toomreParameterCritical]} is a dimensionless constant of order unity which controls where the critical
    density occurs. $\sigma_\mathrm{gas}$ is assumed to be a constant equal to {\normalfont \ttfamily [velocityDispersionDiskGas]} and the disk is
    assumed to have a flat rotation curve such that $\kappa = \sqrt{2} V/R$.
    !!}
    use :: Abundances_Structure            , only : abundances
    use :: Coordinates                     , only : coordinateCylindrical         , assignment(=)
    use :: Galactic_Structure_Options      , only : componentTypeDisk             , massTypeGaseous
    use :: Galacticus_Nodes                , only : nodeComponentDisk
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (starFormationRateSurfaceDensityDisksKennicuttSchmidt), intent(inout) :: self
    type            (treeNode                                            ), intent(inout) :: node
    double precision                                                      , intent(in   ) :: radius
    class           (nodeComponentDisk                                   ), pointer       :: disk
    class           (massDistributionClass                               ), pointer       :: massDistribution_
    type            (abundances                                          ), save          :: abundancesFuel
    !$omp threadprivate(abundancesFuel)
    double precision                                                                      :: surfaceDensityCritical, massGas, &
         &                                                                                   surfaceDensityGas
    type            (coordinateCylindrical                               )                :: coordinates

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Check if factors have been precomputed.
    if (.not.self%factorsComputed) then
       ! Get the disk properties.
       disk   => node%disk   ()
       massGas = disk%massGas()
       ! Find the hydrogen fraction in the disk gas of the fuel supply.
       abundancesFuel=disk%abundancesGas()
       call abundancesFuel%massToMassFraction(massGas)
       self%hydrogenMassFraction=abundancesFuel%hydrogenMassFraction()
       ! Compute the constant factor appearing in the critical density.
       self%surfaceDensityCriticalFactor=+self%toomreParameterCritical     &
            &                            *sqrt(2.0d0)                      &
            &                            *self%velocityDispersionDiskGas   &
            &                            *disk%velocity                 () &
            &                            /Pi                               &
            &                            /gravitationalConstant_internal
       ! Record that factors have now been computed.
       self%factorsComputed=.true.
    end if
    ! Get gas surface density.
    coordinates       =  [radius,0.0d0,0.0d0]
    massDistribution_ => node             %massDistribution(componentType=componentTypeDisk,massType=massTypeGaseous)
    surfaceDensityGas =  massDistribution_%surfaceDensity  (              coordinates                               )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Compute the star formation rate surface density.
    kennicuttSchmidtRate=+self%normalization          &
         &               *(                           &
         &                 +self%hydrogenMassFraction &
         &                 *surfaceDensityGas         &
         &                )**self%exponent
    ! Check if we are applying a truncation radius.
    if (self%truncate) then
       ! Always return zero star formation rate at zero radius, as critical density will be infinite.
       if (radius <= 0.0d0) then
          kennicuttSchmidtRate=0.0d0
          return
       end if
       ! Compute the critical density for star formation.
       surfaceDensityCritical=+self%surfaceDensityCriticalFactor &
            &                 /radius
       ! Check if gas is above the critical density. Truncate it if not.
       if (surfaceDensityGas < surfaceDensityCritical)        &
            & kennicuttSchmidtRate=+kennicuttSchmidtRate      &
            &                      *(                         &
            &                        +surfaceDensityGas       &
            &                        /surfaceDensityCritical  &
            &                       )**self%exponentTruncated
    end if
    return
  end function kennicuttSchmidtRate
