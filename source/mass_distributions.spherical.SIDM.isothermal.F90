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
  Provides a mass distribution implementing the ``isothermal'' approximation to the effects of SIDM based on the model of \cite{jiang_semi-analytic_2023}.
  !!}

  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Numerical_Interpolation, only : interpolator
  use            :: Numerical_ODE_Solvers  , only : odeSolver

  !![
  <massDistribution name="massDistributionSphericalSIDMIsothermal">
   <description>
      A mass distribution class for self-interacting dark matter following the ``isothermal'' model of \cite{jiang_semi-analytic_2023}. This
      model assumes that the dark matter within the interaction radius, $r_1$, has thermalized and can therefore be described by a
      constant velocity dispersion, $\sigma_0$. Under this assumption the spherical Jeans equation has a solution of the form:
      \begin{equation}
      \rho(r) = \rho_0 \exp\left[-\frac{\phi(r)}{\sigma_0^2}\right],
      \end{equation}
      where $\rho(r)$ is the density $\rho_0$ is the density at $r=0$, and the gravitational potential satisfies \citep{jiang_semi-analytic_2023}:
      \begin{equation}
      \nabla^2 \phi(r) = 4 \pi \mathrm{G} \rho_0 \exp \left( - \frac{\phi(r)}{\sigma_0^2} \right).
      \end{equation}
      This second-order differential equation is solved using the boundary conditions $\phi(r=0)=0$ and
      $\mathrm{d}\phi/\mathrm{d}r(r=0)=0$. The values of $\rho_0$ and $\sigma_0$ are then found by minimizing a function      
      \begin{equation}
      \delta^2(\rho_0,\sigma_0) = \left[ \frac{\rho(r_1)}{\rho^\prime(r_1)} - 1 \right]^2 + \left[ \frac{M(r_1)}{M^\prime(r_1)} - 1 \right]^2,
      \end{equation}
      where $M(r)$ is the mass contained within radius $r$, and primes indicate the profile prior to SIDM thermalization.

      This can be expressed in a convenient dimensionless form. We define $x=r/r_1$, $y=\rho/\rho_1$, $z=\sigma/\sigma_1$, where
      \begin{equation}
       \sigma_1^2 = \frac{4 \pi}{3} \mathrm{G} \rho_1 r_1^2 \xi,
      \end{equation}
      and we define $\xi$ through the relation:
      \begin{equation}
       M_1 = \xi \frac{4 \pi}{3} \rho_1 r_1^3.
      \end{equation}
      Using these definitions we can define a dimensionless potential, $\Phi(r) = \phi(r) / \sigma_1^2$. The above differential
      equation can then be written as
      \begin{equation}
      \nabla^{\prime 2} \Phi = \frac{3}{\xi} y_0 \exp\left[ - \frac{\Phi}{z_0^2} \right] ,
      \end{equation}
      where $\nabla^{\prime 2}$ indicates the Laplacian with respect to coordinate $x$. Written in this form it is straightforward
      to see that this equation has three parameters, $\xi$, $y_0$, and $z_0$. The value of $\xi$ is determined from the initial
      (pre-thermalization) density profile. We then have two constraints at $x=1$, namely $y=1$ and $m=M/M_1=1$. We can solve for
      the values of $y_0$ and $z_0$ which satisfy these constraints for a given $\xi$. As a result, we can tabulate solutions
      $y_0(\xi)$ and $z_0(\xi)$ which are applicable to any initial density profile and depend only on the effective slope of the
      density profile inside $r_1$, since if $\rho \propto r^\alpha$ then $\xi = 1/(1+\alpha/3)$, such that $\alpha=0$ (the
      largest physically-allowed value of $\alpha$) implies $\xi=1$.
   </description>
  </massDistribution>
  !!]
  type, extends(massDistributionSphericalSIDM) :: massDistributionSphericalSIDMIsothermal
     !!{
     A mass distribution implementing the ``isothermal'' approximation to the effects of SIDM based on the model of \cite{jiang_semi-analytic_2023}.
     !!}
     private
     double precision                                                :: velocityDispersionCentral     , radiusInteraction_                    , &
          &                                                             densityInteraction            , massInteraction                       , &
          &                                                             velocityDispersionInteraction
     type            (interpolator  ), allocatable                   :: densityCentralDimensionless   , velocityDispersionCentralDimensionless, &
          &                                                             interpolatorRadiiDimensionless, interpolatorXi
     double precision                                                :: xiTabulatedMinimum            , xiTabulatedMaximum
     double precision                , allocatable, dimension( : ,:) :: densityProfileDimensionless   , massProfileDimensionless
     double precision                , allocatable, dimension( :   ) :: radiiDimensionless
     integer         (c_size_t      )                                :: indexXi
     double precision                             , dimension(0:1  ) :: factorsXi
   contains
     !![
     <methods>
       <method method="tabulateSolutions" description="Tabulate solutions for the isothermal core of a SIDM halo."/>
       <method method="computeSolution"   description="Compute a solution for the isothermal core of a SIDM halo."/>
     </methods>
     !!]
     final     ::                           sphericalSIDMIsothermalDestructor
     procedure :: computeSolution        => sphericalSIDMIsothermalComputeSolution
     procedure :: tabulateSolutions      => sphericalSIDMIsothermalTabulateSolutions
     procedure :: density                => sphericalSIDMIsothermalDensity
     procedure :: densityGradientRadial  => sphericalSIDMIsothermalDensityGradientRadial
     procedure :: massEnclosedBySphere   => sphericalSIDMIsothermalMassEnclosedBySphere
     procedure :: potential              => sphericalSIDMIsothermalPotential
     procedure :: potentialIsAnalytic    => sphericalSIDMIsothermalPotentialIsAnalytic
  end type massDistributionSphericalSIDMIsothermal

  interface massDistributionSphericalSIDMIsothermal
     !!{
     Constructors for the \refClass{massDistributionSphericalSIDMIsothermal} mass distribution class.
     !!}
     module procedure sphericalSIDMIsothermalConstructorParameters
     module procedure sphericalSIDMIsothermalConstructorInternal
  end interface massDistributionSphericalSIDMIsothermal

  ! Number of properties in ODE.
  integer         (c_size_t ), parameter   :: propertyCount=2

  ! Submodule-scope variables.
  double precision                         :: xi_            , y0_, &
       &                                      z0_
  type            (odeSolver), allocatable :: odeSolver_
  !$omp threadprivate(xi_,y0_,z0_,odeSolver_)

contains

  function sphericalSIDMIsothermalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalSIDMIsothermal} mass distribution class which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSphericalSIDMIsothermal)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (massDistributionClass                  ), pointer       :: massDistribution_
    class           (darkMatterParticleClass                ), pointer       :: darkMatterParticle_
    type            (varying_string                         )                :: componentType      , massType, &
         &                                                                      nonAnalyticSolver
    double precision                                                         :: timeAge

    !![
    <inputParameter>
      <name>timeAge</name>
      <source>parameters</source>
      <description>The age of the halo (in Gyr).</description>
    </inputParameter>
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <inputParameter>
      <name>componentType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The component type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The mass type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="massDistribution"   name="massDistribution_"   source="parameters"/>
    <objectBuilder class="darkMatterParticle" name="darkMatterParticle_" source="parameters"/>
    !!]
    select type (massDistribution_)
    class is (massDistributionSpherical)
       self=massDistributionSphericalSIDMIsothermal(timeAge,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),massDistribution_,darkMatterParticle_,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    class default
       call Error_Report('a spherically-symmetric mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"  />
    <objectDestructor name="darkMatterParticle_"/>
    !!]
    return
  end function sphericalSIDMIsothermalConstructorParameters

  function sphericalSIDMIsothermalConstructorInternal(timeAge,nonAnalyticSolver,massDistribution_,darkMatterParticle_,componentType,massType) result(self)
    !!{
    Internal constructor for the \refClass{massDistributionSphericalSIDMIsothermal} mass distribution class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleSelfInteractingDarkMatter
    implicit none
    type            (massDistributionSphericalSIDMIsothermal)                          :: self
    double precision                                         , intent(in   )           :: timeAge
    class           (massDistributionSpherical              ), intent(in   ), target   :: massDistribution_
    class           (darkMatterParticleClass                ), intent(in   ), target   :: darkMatterParticle_
    type            (enumerationNonAnalyticSolversType      ), intent(in   )           :: nonAnalyticSolver
    type            (enumerationComponentTypeType           ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType                ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="timeAge, nonAnalyticSolver, componentType, massType, *massDistribution_, *darkMatterParticle_"/>
    !!]

    ! Validate the dark matter particle type.
    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       ! This is as expected.
    class default
       call Error_Report('this class expects a self-interacting dark matter particle'//{introspection:location})
    end select
    self%dimensionless     =.false.
    self%xiTabulatedMinimum=+huge(0.0d0)
    self%xiTabulatedMaximum=-huge(0.0d0)
    call self%computeSolution()
    return
  end function sphericalSIDMIsothermalConstructorInternal

  subroutine sphericalSIDMIsothermalDestructor(self)
    !!{
    Destructor for the \refClass{massDistributionSphericalSIDMIsothermal} class.
    !!}
    implicit none
    type(massDistributionSphericalSIDMIsothermal), intent(inout) :: self
    
    !![
    <objectDestructor name="self%massDistribution_"  />
    <objectDestructor name="self%darkMatterParticle_"/>
    !!]
    return
  end subroutine sphericalSIDMIsothermalDestructor

  subroutine sphericalSIDMIsothermalTabulateSolutions(self,xiRequired)
    !!{
    Tabulate solutions for $y_0(\xi)$, $z_0(\xi)$.
    !!}
    use :: Display                   , only : displayIndent  , displayUnindent    , displayMessage, verbosityLevelWorking, &
         &                                    displayCounter , displayCounterClear
    use :: File_Utilities            , only : Directory_Make , File_Exists        , File_Lock     , File_Path            , &
         &                                    File_Unlock    , lockDescriptor
    use :: Input_Paths               , only : inputPath      , pathTypeDataDynamic
    use :: HDF5_Access               , only : hdf5Access
    use :: IO_HDF5                   , only : hdf5Object
    use :: Numerical_Ranges          , only : Make_Range     , rangeTypeLinear
    use :: Multidimensional_Minimizer, only : multiDMinimizer
    implicit none
    class           (massDistributionSphericalSIDMIsothermal), intent(inout)                             :: self
    double precision                                         , intent(in   )                             :: xiRequired
    integer                                                  , parameter                                 :: countXiPerUnit      = 100
    integer                                                  , parameter                                 :: countRadii          =1000
    double precision                                         , parameter                                 :: Y0Minimum           =1.0d+0, Y0Maximum           =1.0d+6
    double precision                                         , parameter                                 :: Z0Minimum           =0.1d+0, Z0Maximum           =3.0d+0
    double precision                                         , parameter                                 :: xiMinimum           =1.1d+0, xiMaximum           =1.0d+1
    double precision                                         , parameter                                 :: x1                  =1.0d+0
    double precision                                         , parameter                                 :: odeToleranceAbsolute=1.0d-9, odeToleranceRelative=1.0d-9
    double precision                                         , dimension(propertyCount+1  )              :: properties                 , propertyScales
    double precision                                         , dimension(propertyCount    )              :: locationMinimum
    double precision                                         , dimension(              :  ), allocatable :: xi                         , y0                         , &
         &                                                                                                  z0
    double precision                                                                                     :: x
    type            (multiDMinimizer                        )                              , allocatable :: minimizer_
    integer                                                                                              :: countXi                    , count
    integer         (c_size_t                               )                                            :: i                          , j                          , &
         &                                                                                                  iteration
    logical                                                                                              :: converged                  , retabulate
    type            (varying_string                         )                                            :: fileName
    type            (hdf5Object                             )                                            :: file
    type            (lockDescriptor                         )                                            :: fileLock
    character       (len=16                                 )                                            :: labelXiMinimum             , labelXiMaximum
 
    ! Return immediately if solutions have been tabulated with sufficient extent already.
    if     (                                       &
         &   xiRequired >= self%xiTabulatedMinimum &
         &  .and.                                  &
         &   xiRequired <= self%xiTabulatedMaximum &
         & ) return
    ! Deallocate existing table if necessary.
    if (allocated(self%radiiDimensionless                    )) deallocate(self%radiiDimensionless                    )
    if (allocated(self%densityProfileDimensionless           )) deallocate(self%densityProfileDimensionless           )
    if (allocated(self%massProfileDimensionless              )) deallocate(self%massProfileDimensionless              )
    if (allocated(self%interpolatorXi                        )) deallocate(self%interpolatorXi                        )
    if (allocated(self%interpolatorRadiiDimensionless        )) deallocate(self%interpolatorRadiiDimensionless        )
    if (allocated(self%densityCentralDimensionless           )) deallocate(self%densityCentralDimensionless           )
    if (allocated(self%velocityDispersionCentralDimensionless)) deallocate(self%velocityDispersionCentralDimensionless)
    ! By default assume that we do need to retabulate.
    retabulate=.true.
    ! Construct a file name for the table.
    fileName=inputPath(pathTypeDataDynamic)// &
         &   'darkMatter/'                 // &
         &   self%objectType()             // &
         &   '.hdf5'
    call Directory_Make(char(File_Path(char(fileName))))
    if (File_Exists(fileName)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
       ! Restore tables from file.
       !$ call hdf5Access%set()
       file=hdf5Object(char(fileName))
       call file%readDataset('xi'                         ,     xi                         )
       call file%readDataset('radii'                      ,self%radiiDimensionless         )
       call file%readDataset('y0'                         ,     y0                         )
       call file%readDataset('z0'                         ,     z0                         )
       call file%readDataset('densityProfileDimensionless',self%densityProfileDimensionless)
       call file%readDataset('massProfileDimensionless'   ,self%massProfileDimensionless   )
       !$ call hdf5Access%unset()
       self%xiTabulatedMinimum=xi(      1 )
       self%xiTabulatedMaximum=xi(size(xi))
       ! Check if the table is sufficient.
       retabulate= xiRequired < self%xiTabulatedMinimum &
            &     .or.                                  &
            &      xiRequired > self%xiTabulatedMaximum
       call File_Unlock(fileLock)
    end if
    ! Retabulate now if necessary.
    if (retabulate) then
       if (allocated(     xi                         )) deallocate(     xi                         )
       if (allocated(     y0                         )) deallocate(     y0                         )
       if (allocated(     z0                         )) deallocate(     z0                         )
       if (allocated(self%radiiDimensionless         )) deallocate(self%radiiDimensionless         )
       if (allocated(self%densityProfileDimensionless)) deallocate(self%densityProfileDimensionless)
       if (allocated(self%massProfileDimensionless   )) deallocate(self%massProfileDimensionless   )
       ! Set extent for tabulation.
       self%xiTabulatedMinimum=min(1.0d0*xiRequired,xiMinimum)
       self%xiTabulatedMaximum=max(1.1d0*xiRequired,xiMaximum)
       write (labelXiMinimum,'(f5.2)') self%xiTabulatedMinimum
       write (labelXiMaximum,'(f5.2)') self%xiTabulatedMaximum
       call displayIndent ('tabulating isothermal SIDM density profile solutions'                                ,verbosityLevelWorking)
       call displayMessage('range: '//trim(adjustl(labelXiMinimum))//' < ξ < '//trim(adjustl(labelXiMaximum))//'',verbosityLevelWorking)
       ! Construct ranges of the parameter ξ to span.
       countXi=int((self%xiTabulatedMaximum-self%xiTabulatedMinimum)*dble(countXiPerUnit))+1
       allocate(     xi                         (           countXi))
       allocate(     y0                         (           countXi))
       allocate(     z0                         (           countXi))
       allocate(self%radiiDimensionless         (countRadii        ))
       allocate(self%densityProfileDimensionless(countRadii,countXi))
       allocate(self%massProfileDimensionless   (countRadii,countXi))
       xi                     =Make_Range(self%xiTabulatedMinimum,self%xiTabulatedMaximum,countXi   ,rangeTypeLinear)
       self%radiiDimensionless=Make_Range(     0.0d0             ,     1.0d0             ,countRadii,rangeTypeLinear)
       ! Set absolute property scales for ODE solving.
       propertyScales=1.0d0
       ! Start parallel region to solve for halo structure at each value of ξ.
       count=0
       call displayCounter(count,isNew=.true.,verbosity=verbosityLevelWorking)
       !$omp parallel private(i,j,x,properties,locationMinimum,iteration,converged,minimizer_)
       !! Allocate and construct objects needed by each thread.
       allocate(odeSolver_)
       allocate(minimizer_)
       odeSolver_=odeSolver      (propertyCount+1,sphericalSIDMIsothermalDimensionlessODEs     ,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative,scale=propertyScales)
       minimizer_=multiDMinimizer(propertyCount  ,sphericalSIDMIsothermalDimensionlessFitMetric                                                                                                   )
       !$omp do schedule(dynamic)
       do i=1,countXi
          xi_=xi(i)
          ! Seek the low-density solution.
          call minimizer_%set(x=[0.0d0,1.0d0],stepSize=[0.01d0,0.01d0])
          iteration=0
          converged=.false.
          do while (.not.converged .and. iteration < 100)
             call minimizer_%iterate()
             iteration=iteration+1
             converged=minimizer_%testSize(toleranceAbsolute=1.0d-12)
          end do
          locationMinimum=minimizer_%x()
          y0(i)=exp(locationMinimum(1))
          z0(i)=    locationMinimum(2)
          ! Tabulate solutions for density and mass.
          do j=1,countRadii
             x         =0.0d0
             properties=0.0d0
             call odeSolver_%solve(x,self%radiiDimensionless(j),properties)
             self%densityProfileDimensionless(j,i)=+y0(i)                 &
                  &                                *exp(                  &
                  &                                     -properties(1)    &
                  &                                     /z0(i)        **2 &
                  &                                    )
             self%massProfileDimensionless   (j,i)=+     properties(3)
          end do
          !$omp atomic
          count=count+1
          call displayCounter(int(100.0d0*dble(count)/dble(countXi)),isNew=.false.,verbosity=verbosityLevelWorking)
       end do
       !$omp end do
       call displayCounterClear(verbosityLevelWorking)
       deallocate(odeSolver_)
       deallocate(minimizer_)
       !$omp end parallel
       ! Write the data to file.
       call File_Lock(char(fileName),fileLock,lockIsShared=.false.)
       !$ call hdf5Access%set()
       file=hdf5Object(char(fileName),overWrite=.true.,readOnly=.false.)
       call file%writeDataset(     xi                          ,'xi'                         )
       call file%writeDataset(self%radiiDimensionless          ,'radii'                      )
       call file%writeDataset(     y0                          ,'y0'                         )
       call file%writeDataset(     z0                          ,'z0'                         )
       call file%writeDataset(self%densityProfileDimensionless ,'densityProfileDimensionless')
       call file%writeDataset(self%massProfileDimensionless    ,'massProfileDimensionless'   )
       !$ call hdf5Access%unset()
       call File_Unlock(fileLock)
       call displayUnindent('done',verbosityLevelWorking)
    end if
    ! Build the interpolators.
    allocate(self%interpolatorXi                        )
    allocate(self%interpolatorRadiiDimensionless        )
    allocate(self%densityCentralDimensionless           )
    allocate(self%velocityDispersionCentralDimensionless)
    self%densityCentralDimensionless           =interpolator(     xi                ,y0)
    self%velocityDispersionCentralDimensionless=interpolator(     xi                ,z0)
    self%interpolatorXi                        =interpolator(     xi                   )
    self%interpolatorRadiiDimensionless        =interpolator(self%radiiDimensionless   )
    return
  end subroutine sphericalSIDMIsothermalTabulateSolutions
  
  double precision function sphericalSIDMIsothermalDimensionlessFitMetric(propertiesCentral)
    !!{
    Evaluate the fit metric.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:)               :: propertiesCentral
    double precision, parameter                                 :: x1               =1.0d0
    double precision               , dimension(propertyCount+1) :: properties
    double precision                                            :: x                      , y1, &
         &                                                         m1

    ! Extract current parameters to submodule-scope.
    y0_=exp(propertiesCentral(1))
    z0_=    propertiesCentral(2)
    ! Solve the ODE to x₁.
    x         =0.0d0
    properties=0.0d0
    call odeSolver_%solve(x,x1,properties)
    ! Extract density and mass at x₁.
    y1     =+y0_                   &
         &  *exp(                  &
         &       -properties(1)    &
         &       /z0_          **2 &
         &      )
    m1     =+     properties(3)
    ! Evaluate the fit metric.
    sphericalSIDMIsothermalDimensionlessFitMetric=+(y1-1.0d0)**2 &
         &                               +(m1-1.0d0)**2
    return
  end function sphericalSIDMIsothermalDimensionlessFitMetric
  
  integer function sphericalSIDMIsothermalDimensionlessODEs(x,properties,propertiesRateOfChange)
    !!{
    Define the dimensionless ODE system to solve for isothermal self-interacting dark matter cores.
    !!}
    use :: Interface_GSL, only : GSL_Success
    implicit none
    double precision, intent(in   )               :: x
    double precision, intent(in   ), dimension(:) :: properties
    double precision, intent(  out), dimension(:) :: propertiesRateOfChange
    double precision                              :: y

    ! Compute the dimensionless density.
    y                               =+y0_                           &
           &                         *exp(                          &
           &                              -max(properties(1),0.0d0) &
           &                              /z0_**2                   &
           &                             )
    ! Evaluate the ODE.
    propertiesRateOfChange       (1)=+properties(2)
    propertiesRateOfChange       (2)=+3.0d0                         &
         &                           /xi_                           &
         &                           *y
    if (x > 0.0d0)                                                  &
         & propertiesRateOfChange(2)=+propertiesRateOfChange(2)     &
         &                           -2.0d0                         &
         &                           *properties            (2)     &
         &                           /x
    propertiesRateOfChange       (3)=+3.0d0                         &
         &                           /xi_                           &
         &                           *x**2                          &
         &                           *y
    sphericalSIDMIsothermalDimensionlessODEs=GSL_Success
    return
  end function sphericalSIDMIsothermalDimensionlessODEs
    
  subroutine sphericalSIDMIsothermalComputeSolution(self)
    !!{
    Compute a solution for the isothermal core of an SIDM halo.
    !!}
    use :: Coordinates                     , only : coordinateSpherical           , assignment(=)
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionSphericalSIDMIsothermal), intent(inout) :: self
    integer                                                  , parameter     :: countTable                   =1000
    double precision                                         , parameter     :: odeToleranceAbsolute         =1.0d-3, odeToleranceRelative     =1.0d-3
    double precision                                                         :: densityCentral                      , velocityDispersionCentral       , &
         &                                                                      densityInteraction                  , massInteraction                 , &
         &                                                                      radiusInteraction                   , xi                              , &
         &                                                                      velocityDispersionInteraction
    type            (coordinateSpherical                    )                :: coordinatesInteraction
    
    ! Find the interaction radius.
    radiusInteraction            =self%radiusInteraction                     (                      )
    coordinatesInteraction=[radiusInteraction,0.0d0,0.0d0]
    ! Properties of the original density profile at the interaction radius.
    densityInteraction           =self%massDistribution_%density             (coordinatesInteraction)
    massInteraction              =self%massDistribution_%massEnclosedBySphere(radiusInteraction     )
    ! Find the velocity dispersion scale to be applied to the dimensionless solutions.
    velocityDispersionInteraction=sqrt(gravitationalConstant_internal*massInteraction/radiusInteraction)
    ! Compute the ξ parameter.
    xi                           =+massInteraction       &
         &                        *3.0d0                 &
         &                        /4.0d0                 &
         &                        /Pi                    &
         &                        /densityInteraction    &
         &                        /radiusInteraction **3
    ! Ensure dimensionless solutions have been tabulated.
    call self%tabulateSolutions(xi)
    ! Find the properties at the halo center.
    densityCentral               =self%densityCentralDimensionless           %interpolate(xi)*densityInteraction
    velocityDispersionCentral    =self%velocityDispersionCentralDimensionless%interpolate(xi)*velocityDispersionInteraction
    ! Store properties of current profile.
    self%radiusInteraction_           =radiusInteraction
    self%densityInteraction           =densityInteraction
    self%massInteraction              =massInteraction
    self%velocityDispersionInteraction=velocityDispersionInteraction
    self%velocityDispersionCentral    =velocityDispersionCentral
    ! Compute interpolating factors in ξ.
    call self%interpolatorXi%linearFactors(xi,self%indexXi,self%factorsXi)
    return
  end subroutine sphericalSIDMIsothermalComputeSolution

  double precision function sphericalSIDMIsothermalDensity(self,coordinates) result(density)
    !!{
    Compute the density at the specified {\normalfont \ttfamily coordinates} for the {\normalfont \ttfamily sphericalSIDMIsothermal}
    mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalSIDMIsothermal), intent(inout)  :: self
    class           (coordinate                             ), intent(in   )  :: coordinates
    integer         (c_size_t                               )                 :: i            , j, &
         &                                                                       indexRadius
    double precision                                         , dimension(0:1) :: factorsRadius

    if (coordinates%rSpherical() > self%radiusInteraction()) then
       density=self%massDistribution_%density(coordinates)
    else
       call self%interpolatorRadiiDimensionless%linearFactors(coordinates%rSpherical()/self%radiusInteraction_,indexRadius,factorsRadius)
       density=0.0d0
       do i=0,1
          do j=0,1
             density=+     density                                                   &
                  &  +self%densityProfileDimensionless(indexRadius+i,self%indexXi+j) &
                  &  *     factorsRadius              (           i                ) &
                  &  *self%factorsXi                  (                           j)
          end do
       end do
       density=+     density            &
            &  *self%densityInteraction
    end if
    return
  end function sphericalSIDMIsothermalDensity

  double precision function sphericalSIDMIsothermalDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a truncated spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalSIDMIsothermal), intent(inout) , target   :: self
    class           (coordinate                             ), intent(in   )            :: coordinates
    logical                                                  , intent(in   ) , optional :: logarithmic
    integer         (c_size_t                               )                           :: indexRadius
    double precision                                         , dimension(0:1)           :: factorsRadius
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    if (coordinates%rSpherical() > self%radiusInteraction()) then
       densityGradient=self%massDistribution_%densityGradientRadial(coordinates,logarithmic=logarithmic)
    else
       call self%interpolatorRadiiDimensionless%linearFactors(coordinates%rSpherical()/self%radiusInteraction_,indexRadius,factorsRadius)
       if (indexRadius > 1) then
          densityGradient=+log(self%densityProfileDimensionless(indexRadius+1,self%indexXi+0)/self%densityProfileDimensionless(indexRadius+0,self%indexXi+0)) &
               &          /log(self%radiiDimensionless         (indexRadius+1               )/self%radiiDimensionless         (indexRadius+0               ))
          if (.not.logarithmic_)                                         &
               densityGradient=+            densityGradient              &
               &               *self       %density        (coordinates) &
               &               /coordinates%rSpherical     (           )
       else
          densityGradient=+0.0d0
       end if
    end if
    return
  end function sphericalSIDMIsothermalDensityGradientRadial
  
  double precision function sphericalSIDMIsothermalMassEnclosedBySphere(self,radius) result(mass)
    !!{   
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for the {\normalfont \ttfamily sphericalSIDMIsothermal}
    mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalSIDMIsothermal), intent(inout) , target :: self
    double precision                                         , intent(in   )          :: radius
    integer         (c_size_t                               )                         :: i            , j, &
         &                                                                               indexRadius
    double precision                                         , dimension(0:1)         :: factorsRadius

    if (radius > self%radiusInteraction()) then
       mass=self%massDistribution_%massEnclosedBySphere(radius)
    else
       call self%interpolatorRadiiDimensionless%linearFactors(radius/self%radiusInteraction_,indexRadius,factorsRadius)
       mass=0.0d0
       do i=0,1
          do j=0,1
             mass   =+     mass                                                   &
                  &  +self%massProfileDimensionless(indexRadius+i,self%indexXi+j) &
                  &  *     factorsRadius           (            i               ) &
                  &  *self%factorsXi               (                           j)
          end do
       end do
       mass   =+     mass            &
            &  *self%massInteraction
    end if
    return
  end function sphericalSIDMIsothermalMassEnclosedBySphere

  logical function sphericalSIDMIsothermalPotentialIsAnalytic(self) result(isAnalytic)
    !!{
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionSphericalSIDMIsothermal), intent(inout) :: self

    isAnalytic=.true.
    return
  end function sphericalSIDMIsothermalPotentialIsAnalytic

  double precision function sphericalSIDMIsothermalPotential(self,coordinates,status) result(potential)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in an burkert mass distribution.
    !!}
    use :: Coordinates               , only : coordinateSpherical      , assignment(=)
    use :: Galactic_Structure_Options, only : structureErrorCodeSuccess
    implicit none
    class(massDistributionSphericalSIDMIsothermal), intent(inout), target   :: self
    class(coordinate                             ), intent(in   )           :: coordinates
    type (enumerationStructureErrorCodeType      ), intent(  out), optional :: status
    type (coordinateSpherical                    )                          :: coordinatesInteraction
    
    if (present(status)) status=structureErrorCodeSuccess
    if (coordinates%rSpherical() > self%radiusInteraction()) then
       potential             =+     self%massDistribution_%potential                (coordinates           )
    else
       coordinatesInteraction=[self%radiusInteraction_,0.0d0,0.0d0]
       potential             =+     self%massDistribution_%potential                (coordinatesInteraction)    &
            &                 -     self                  %velocityDispersionCentral                        **2 &
            &                 *log(                                                                             &
            &                      +self                  %density                  (coordinates           )    &
            &                      /self                  %density                  (coordinatesInteraction)    &
            &                     )
    end if
    return
  end function sphericalSIDMIsothermalPotential
