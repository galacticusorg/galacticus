!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

  !+    Contributions to this file made by: Andrew Benson, Christoph Behrens, Xiaolong Du.

  !% An implementation of cosmological density field mass variance computed using a filtered power spectrum.

  !# <cosmologicalMassVariance name="cosmologicalMassVarianceFilteredPower">
  !#  <description>
  !#   Mass variance of cosmological density fields computed from a filtered power spectrum:
  !#   \begin{equation}
  !#    \sigma^2(M) = {1 \over 2 \pi^2} \int_0^\infty P(k) T^2(k) W^2(k) k^2 \mathrm{d}k
  !#   \end{equation}
  !#   where $P(k)$ is the primordial power spectrum (see \refPhysics{powerSpectrumPrimordial}), $T(k)$ is the transfer function
  !#   (see \refPhysics{transferFunction}), and $W(k)$ is the power spectrum variance window function (see
  !#   \refPhysics{powerSpectrumWindowFunction}).
  !#
  !#   The normalization of the mass variance is specified via the {\normalfont \ttfamily [sigma\_8]} parameter, which defines the
  !#   linear theory root-variance of the density field in spheres of radii $8h^{-1}$Mpc. Note that when computing the
  !#   normalization of the power spectrum to match the specified value of $\sigma_8$ a top-hat real-space window function is
  !#   used (as per the definition of $\sigma_8$), unless a different window function is explicitly defined via the {\normalfont
  !#   \ttfamily [powerSpectrumWindowFunctionTopHat]} parameter.
  !#
  !#   The mass variance, $\sigma(M)$, is found by integration over the linear theory power spectrum, with the specified power
  !#   spectrum window function. The fractional tolerance for this integration can be set via the {\normalfont \ttfamily
  !#   [tolerance]} parameter. (The normalization of $\sigma(M)$ to give the desired $\sigma_8$ always uses a top-hat window
  !#   function. For this integration the tolerance can be set via the {\normalfont \ttfamily [toleranceTopHat]} parameter.) This
  !#   is tabulated across the required range.
  !#
  !#   Cubic spline interpolation is then used to interpolate in this table to give $\sigma(M)$ at any required value of $M$. The
  !#   tabulation is always forced to be monotonically decreasing with $M$. However, the interpolation is not necessarily
  !#   monotonic---for example in cases where $\sigma(M)$ becomes constant or close to constant as a function of $M$ the
  !#   interpolation can become non-monotonic over some ranges of $M$. If strict monotonicity is required set {\normalfont
  !#   \ttfamily [monotonicInterpolation]}={\normalfont \ttfamily true}. This causes a monotonic spline interpolator to be used
  !#   instead which gaurantees monotonicity.
  !#  </description>
  !# </cosmologicalMassVariance>
  use :: Cosmology_Functions                 , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters                , only : cosmologyParametersClass
  use :: Linear_Growth                       , only : linearGrowthClass
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredClass
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionClass
  use :: Tables                              , only : table1DLinearCSpline

  !# <stateStorable class="uniqueTable"/>

  type :: uniqueTable
     !% Type used to store unique values of the mass variance.
     private
     double precision, allocatable, dimension(:) :: rootVariance
     integer         , allocatable, dimension(:) :: index
  end type uniqueTable

  type, extends(cosmologicalMassVarianceClass) :: cosmologicalMassVarianceFilteredPower
     !% A cosmological mass variance class computing variance from a filtered power spectrum.
     private
     class           (cosmologyParametersClass               ), pointer                   :: cosmologyParameters_                => null()
     class           (cosmologyFunctionsClass                ), pointer                   :: cosmologyFunctions_                 => null()
     class           (powerSpectrumPrimordialTransferredClass), pointer                   :: powerSpectrumPrimordialTransferred_ => null()
     class           (linearGrowthClass                      ), pointer                   :: linearGrowth_                       => null()
     class           (powerSpectrumWindowFunctionClass       ), pointer                   :: powerSpectrumWindowFunction_        => null(), powerSpectrumWindowFunctionTopHat_  => null()
     logical                                                                              :: initialized                                  , nonMonotonicIsFatal
     double precision                                                                     :: tolerance                                    , toleranceTopHat                              , &
          &                                                                                  sigma8Value                                  , sigmaNormalization                           , &
          &                                                                                  massMinimum                                  , massMaximum                                  , &
          &                                                                                  timeMinimum                                  , timeMaximum                                  , &
          &                                                                                  timeMinimumLogarithmic                       , timeLogarithmicDeltaInverse
     double precision                                         , allocatable, dimension(:) :: times
     class           (table1DLinearCSpline                   ), allocatable, dimension(:) :: rootVarianceTable
     type            (varying_string                         )                            :: fileName
     ! Unique values in the variance table and their corresponding indices.
     type            (uniqueTable                            ), allocatable, dimension(:) :: rootVarianceUniqueTable
     logical                                                                              :: monotonicInterpolation                       , growthIsMassDependent_
   contains
     !# <methods>
     !#   <method description="Tabulate cosmological mass variance." method="retabulate" />
     !#   <method description="Compute the interpolating factors in time." method="interpolantsTime" />
     !#   <method description="Write the tabulated mass variance to file." method="fileWrite" />
     !#   <method description="Read the tabulated mass variance from file." method="fileRead" />
     !#   <method description="Return true if the table must be remade." method="remakeTable" />
     !# </methods>
     final     ::                                        filteredPowerDestructor
     procedure :: sigma8                              => filteredPowerSigma8
     procedure :: powerNormalization                  => filteredPowerPowerNormalization
     procedure :: rootVariance                        => filteredPowerRootVariance
     procedure :: rootVarianceLogarithmicGradient     => filteredPowerRootVarianceLogarithmicGradient
     procedure :: rootVarianceLogarithmicGradientTime => filteredPowerRootVarianceLogarithmicGradientTime
     procedure :: rootVarianceAndLogarithmicGradient  => filteredPowerRootVarianceAndLogarithmicGradient
     procedure :: mass                                => filteredPowerMass
     procedure :: retabulate                          => filteredPowerRetabulate
     procedure :: interpolantsTime                    => filteredPowerInterpolantsTime
     procedure :: growthIsMassDependent               => filteredPowerGrowthIsMassDependent
     procedure :: fileWrite                           => filteredPowerFileWrite
     procedure :: fileRead                            => filteredPowerFileRead
     procedure :: remakeTable                         => filteredPowerRemakeTable
  end type cosmologicalMassVarianceFilteredPower

  interface cosmologicalMassVarianceFilteredPower
     !% Constructors for the {\normalfont \ttfamily filteredPower} cosmological mass variance class.
     module procedure filteredPowerConstructorParameters
     module procedure filteredPowerConstructorInternal
  end interface cosmologicalMassVarianceFilteredPower

  ! Number of points per decade to use in tabulation of σ(M).
  integer                         , parameter :: filteredPowerTablePointsPerDecade=10, filteredPowerTimePointsPerDecade=100

  ! Module-scope time used in integrals.
  double precision                            :: filteredPowerTime
  !$omp threadprivate(filteredPowerTime)

contains

  function filteredPowerConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily filteredPower} cosmological mass variance class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (cosmologicalMassVarianceFilteredPower  )                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (cosmologyParametersClass               ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class           (powerSpectrumPrimordialTransferredClass), pointer       :: powerSpectrumPrimordialTransferred_
    class           (powerSpectrumWindowFunctionClass       ), pointer       :: powerSpectrumWindowFunction_       , powerSpectrumWindowFunctionTopHat_
    class           (linearGrowthClass                      ), pointer       :: linearGrowth_
    double precision                                                         :: sigma8Value                        , tolerance                         , &
         &                                                                      toleranceTopHat
    logical                                                                  :: monotonicInterpolation             , nonMonotonicIsFatal

    ! Check and read parameters.
    !# <inputParameter>
    !#   <name>sigma_8</name>
    !#   <source>parameters</source>
    !#   <variable>sigma8Value</variable>
    !#   <defaultValue>0.8111d0</defaultValue>
    !#   <defaultSource>(\citealt{planck_collaboration_planck_2018}; TT,TE,EE$+$lowE$+$lensing)</defaultSource>
    !#   <description>The fractional mass fluctuation in the linear density field at the present day in spheres of radius 8~Mpc/h.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>toleranceTopHat</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d-6</defaultValue>
    !#   <description>The relative tolerance to use in integrating over the linear power spectrum using a top-hat (real space) window function to compute the cosmological mass variance.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>tolerance</name>
    !#   <source>parameters</source>
    !#   <defaultValue>4.0d-6</defaultValue>
    !#   <description>The relative tolerance to use in integrating over the linear power spectrum to compute the cosmological mass variance.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>nonMonotonicIsFatal</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>If true any non-monotonicity in the tabulated $\sigma(M)$ is treated as a fatal error. Otherwise a only a warning is issued.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>monotonicInterpolation</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true use a monotonic cubic spline interpolator to interpolate in the $\sigma(M)$ table. Otherwise use a standard cubic spline interpoltor. Use of the monotionic interpolator can be helpful is $\sigma(M)$ must be strictly monotonic but becomes a very weak function of $M$ at low masses.</description>
    !# </inputParameter>
    !# <objectBuilder    class="cosmologyParameters"                name="cosmologyParameters_"                source="parameters"                                                  />
    !# <objectBuilder    class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"                                                  />
    !# <objectBuilder    class="powerSpectrumPrimordialTransferred" name="powerSpectrumPrimordialTransferred_" source="parameters"                                                  />
    !# <objectBuilder    class="powerSpectrumWindowFunction"        name="powerSpectrumWindowFunction_"        source="parameters"                                                  />
    !# <objectBuilder    class="linearGrowth"                       name="linearGrowth_"                       source="parameters"                                                  />
    if (parameters%isPresent('powerSpectrumWindowFunctionTopHat')) then
       !# <objectBuilder class="powerSpectrumWindowFunction"        name="powerSpectrumWindowFunctionTopHat_"  source="parameters" parameterName="powerSpectrumWindowFunctionTopHat"/>
    else
       nullify(powerSpectrumWindowFunctionTopHat_)
    end if
    !# <conditionalCall>
    !#  <call>self=filteredPowerConstructorInternal(sigma8Value,tolerance,toleranceTopHat,nonMonotonicIsFatal,monotonicInterpolation,cosmologyParameters_,cosmologyFunctions_,linearGrowth_,powerSpectrumPrimordialTransferred_,powerSpectrumWindowFunction_{conditions})</call>
    !#  <argument name="powerSpectrumWindowFunctionTopHat_" value="powerSpectrumWindowFunctionTopHat_" parameterPresent="parameters" parameterName="powerSpectrumWindowFunctionTopHat"/>
    !# </conditionalCall>
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor    name="cosmologyParameters_"               />
    !# <objectDestructor    name="cosmologyFunctions_"                />
    !# <objectDestructor    name="linearGrowth_"                      />
    !# <objectDestructor    name="powerSpectrumPrimordialTransferred_"/>
    !# <objectDestructor    name="powerSpectrumWindowFunction_"       />
    if (parameters%isPresent('powerSpectrumWindowFunctionTopHat')) then
       !# <objectDestructor name="powerSpectrumWindowFunctionTopHat_" />
    end if
    return
  end function filteredPowerConstructorParameters

  function filteredPowerConstructorInternal(sigma8,tolerance,toleranceTopHat,nonMonotonicIsFatal,monotonicInterpolation,cosmologyParameters_,cosmologyFunctions_,linearGrowth_,powerSpectrumPrimordialTransferred_,powerSpectrumWindowFunction_,powerSpectrumWindowFunctionTopHat_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily filteredPower} linear growth class.
    use :: File_Utilities                 , only : Directory_Make                   , File_Path
    use :: Galacticus_Paths               , only : galacticusPath                   , pathTypeDataDynamic
    use :: Power_Spectrum_Window_Functions, only : powerSpectrumWindowFunctionTopHat
    implicit none
    type            (cosmologicalMassVarianceFilteredPower  )                                  :: self
    double precision                                         , intent(in   )                   :: tolerance                          , toleranceTopHat                   , &
         &                                                                                        sigma8
    logical                                                  , intent(in   )                   :: nonMonotonicIsFatal                , monotonicInterpolation
    class           (cosmologyParametersClass               ), intent(in   ), target           :: cosmologyParameters_
    class           (cosmologyFunctionsClass                ), intent(in   ), target           :: cosmologyFunctions_
    class           (powerSpectrumPrimordialTransferredClass), intent(in   ), target           :: powerSpectrumPrimordialTransferred_
    class           (powerSpectrumWindowFunctionClass       ), intent(in   ), target           :: powerSpectrumWindowFunction_
    class           (powerSpectrumWindowFunctionClass       ), intent(in   ), target, optional :: powerSpectrumWindowFunctionTopHat_
    class           (linearGrowthClass                      ), intent(in   ), target           :: linearGrowth_
    !# <constructorAssign variables="tolerance, toleranceTopHat, nonMonotonicIsFatal, monotonicInterpolation, *cosmologyParameters_, *cosmologyFunctions_, *linearGrowth_, *powerSpectrumPrimordialTransferred_, *powerSpectrumWindowFunction_, *powerSpectrumWindowFunctionTopHat_"/>

    if (.not.present(powerSpectrumWindowFunctionTopHat_)) then
       allocate(powerSpectrumWindowFunctionTopHat :: self%powerSpectrumWindowFunctionTopHat_)
       select type (powerSpectrumWindowFunctionTopHat__ => self%powerSpectrumWindowFunctionTopHat_)
       type is (powerSpectrumWindowFunctionTopHat)
          !# <referenceConstruct isResult="yes" object="powerSpectrumWindowFunctionTopHat__" constructor="powerSpectrumWindowFunctionTopHat(cosmologyParameters_)"/>  
       end select
    end if
    self%sigma8Value           =sigma8
    self%initialized           =.false.
    self%growthIsMassDependent_=self%powerSpectrumPrimordialTransferred_%growthIsWavenumberDependent()
    self%fileName              =galacticusPath(pathTypeDataDynamic)              // &
         &                      'largeScaleStructure/'                           // &
         &                      self%objectType      (                          )// &
         &                      '_'                                              // &
         &                      self%hashedDescriptor(includeSourceDigest=.true.)// &
         &                      '.hdf5'
    call Directory_Make(File_Path(self%fileName))
    return
  end function filteredPowerConstructorInternal

  subroutine filteredPowerDestructor(self)
    !% Destructor for the {\normalfont \ttfamily filteredPower} linear growth class.
    implicit none
    type   (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    integer                                                       :: i

    !# <objectDestructor name="self%cosmologyParameters_"               />
    !# <objectDestructor name="self%cosmologyFunctions_"                />
    !# <objectDestructor name="self%linearGrowth_"                      />
    !# <objectDestructor name="self%powerSpectrumPrimordialTransferred_"/>
    !# <objectDestructor name="self%powerSpectrumWindowFunction_"       />
    !# <objectDestructor name="self%powerSpectrumWindowFunctionTopHat_" />
    if (self%initialized) then
       do i=1,size(self%rootVarianceTable)
          call self%rootVarianceTable(i)%destroy()
       end do
    end if
    return
  end subroutine filteredPowerDestructor

  double precision function filteredPowerPowerNormalization(self)
    !% Return the normalization of the power spectrum.
    implicit none
    class(cosmologicalMassVarianceFilteredPower), intent(inout) :: self

    call self%retabulate()
    filteredPowerPowerNormalization=self%sigmaNormalization**2
    return
  end function filteredPowerPowerNormalization

  double precision function filteredPowerSigma8(self)
    !% Return the value of $\sigma_8$.
    implicit none
    class(cosmologicalMassVarianceFilteredPower), intent(inout) :: self

    filteredPowerSigma8=self%sigma8Value
    return
  end function filteredPowerSigma8

  double precision function filteredPowerRootVariance(self,mass,time)
    !% Return the root-variance of the cosmological density field in a spherical region containing the given {\normalfont
    !% \ttfamily mass} on average.
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: mass, time
    double precision                                                       :: h
    integer                                                                :: i

    call self%retabulate      (mass,time)
    if (self%growthIsMassDependent_) then
       call self%interpolantsTime(time,i,h)
       filteredPowerRootVariance=+self%rootVarianceTable(i  )%interpolate(mass)*(1.0d0-h) &
            &                    +self%rootVarianceTable(i+1)%interpolate(mass)*       h
    else
       filteredPowerRootVariance=+self%rootVarianceTable(1  )%interpolate(mass)           &
            &                    *self%linearGrowth_         %value      (time)
    end if
    return
  end function filteredPowerRootVariance

  double precision function filteredPowerRootVarianceLogarithmicGradient(self,mass,time)
    !% Return the logairhtmic gradient with respect to mass of the root-variance of the cosmological density field in a spherical
    !% region containing the given {\normalfont \ttfamily mass} on average.
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: mass        , time
    double precision                                                       :: rootVariance

    call self%rootVarianceAndLogarithmicGradient(mass,time,rootVariance,filteredPowerRootVarianceLogarithmicGradient)
    return
  end function filteredPowerRootVarianceLogarithmicGradient

  double precision function filteredPowerRootVarianceLogarithmicGradientTime(self,mass,time)
    !% Return the logarithmic gradient with respect to time of the root-variance of the cosmological density field in a spherical
    !% region containing the given {\normalfont \ttfamily mass} on average.
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: mass, time
    double precision                                                       :: h
    integer                                                                :: i

    call self%retabulate(mass,time)
    if (self%growthIsMassDependent_) then
       call self%interpolantsTime(time,i,h)
       filteredPowerRootVarianceLogarithmicGradientTime=+(                                                         &
            &                                             -self%rootVarianceTable(i  )%interpolate(mass)           &
            &                                             +self%rootVarianceTable(i+1)%interpolate(mass)           &
            &                                            )                                                         &
            &                                           /(                                                         &
            &                                             +self%rootVarianceTable(i  )%interpolate(mass)*(1.0d0-h) &
            &                                             +self%rootVarianceTable(i+1)%interpolate(mass)*       h  &
            &                                            )                                                         &
            &                                           *self%timeLogarithmicDeltaInverse
    else
       filteredPowerRootVarianceLogarithmicGradientTime=+self%linearGrowth_       %logarithmicDerivativeExpansionFactor (time) &
            &                                           *self%cosmologyFunctions_ %expansionRate                       (       &
            &                                            self%cosmologyFunctions_ %expansionFactor                      (time) &
            &                                                                                                          )       &
            &                                           *                                                                time
    end if
    return
  end function filteredPowerRootVarianceLogarithmicGradientTime

  subroutine filteredPowerRootVarianceAndLogarithmicGradient(self,mass,time,rootVariance,rootVarianceLogarithmicGradient)
    !% Return the value and logarithmic gradient with respect to mass of the root-variance of the cosmological density field in a
    !% spherical region containing the given {\normalfont \ttfamily mass} on average.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: mass        , time
    double precision                                       , intent(  out) :: rootVariance, rootVarianceLogarithmicGradient
    double precision                                                       :: wavenumber  , rootVarianceGradient           , &
         &                                                                    interpolant , h
    integer                                                                :: i           , j

    call self%retabulate(mass,time)
    if (self%growthIsMassDependent_) then
       call self%interpolantsTime(time,i,h)
    else
       i=1
       h=0.0d0
    end if
    rootVariance        =0.0d0
    rootVarianceGradient=0.0d0
    do j=i,i+1
       if (j == i) then
          interpolant=1.0d0-h
       else
          interpolant=      h
       end if
       if (interpolant == 0.0d0) cycle
       if (self%powerSpectrumWindowFunction_%amplitudeIsMassIndependent()) then
          ! For the case of a constant window function amplitude the logarithmic gradient can be found analytically.
          wavenumber          =+self%powerSpectrumWindowFunction_       %wavenumberMaximum(                mass    )
          rootVarianceGradient=+rootVarianceGradient                                                                    &
               &               -self%powerSpectrumPrimordialTransferred_%power            (wavenumber,self%times(j))    &
               &               *self%sigmaNormalization                                                             **2 &
               &               *self%powerSpectrumWindowFunction_       %value            (wavenumber,     mass    )**2 &
               &               *                                                           wavenumber               **3 &
               &               /12.0d0                                                                                  &
               &               /Pi                                                                                  **2 &
               &               *interpolant
       else
          ! Compute the gradient by interpolation in the tabulated relation.
          rootVarianceGradient=+rootVarianceGradient                                &
               &               +self%rootVarianceTable(j)%interpolateGradient(mass) &
               &               *                                              mass  &
               &               *interpolant
       end if
       rootVariance=+rootVariance                                &
            &       +self%rootVarianceTable(j)%interpolate(mass) &
            &       *interpolant
    end do
    if (self%powerSpectrumWindowFunction_%amplitudeIsMassIndependent()) then
       rootVarianceLogarithmicGradient=+rootVarianceGradient    &
            &                          /rootVariance        **2
    else
       rootVarianceLogarithmicGradient=+rootVarianceGradient    &
            &                          /rootVariance
    end if
    ! Scale by the linear growth factor if growth is not mass-dependent.
    if (.not.self%growthIsMassDependent_) rootVariance=+rootVariance                   &
         &                                             *self%linearGrowth_%value(time)
    return
  end subroutine filteredPowerRootVarianceAndLogarithmicGradient

  double precision function filteredPowerMass(self,rootVariance,time)
    !% Return the mass corrresponding to the given {\normalfont \ttfamily } root-variance of the cosmological density field.
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: rootVariance   , time
    double precision                                                       :: h              , hTime     , &
         &                                                                    interpolantTime
    integer                                                                :: i              , iBoundLeft, &
         &                                                                    iBoundRight    , k         , &
         &                                                                    j

    ! If the requested root-variance is below the lowest value tabulated, attempt to tabulate to higher mass (lower
    ! root-variance).
    call self%retabulate(time=time)
    do while (rootVariance < self%rootVarianceTable(1)%y(-1))
       call self%retabulate(self%rootVarianceTable(1)%x(-1)*2.0d0,time)
    end do
    ! Get interpolants in time.
    if (self%growthIsMassDependent_) then
       call self%interpolantsTime(time,k,hTime)
    else
       k    =1
       hTime=0.0d0
    end if
    ! If sigma exceeds the highest value tabulated, simply return the lowest tabulated mass.
    if (rootVariance > self%rootVarianceTable(k)%y(1)) then
       filteredPowerMass=self%rootVarianceTable(k)%x(1)
    else
       ! Iterate over times.
       filteredPowerMass=0.0d0
       do j=k,k+1
          ! Compute interpolating factor.
          if (j == k) then
             interpolantTime=1.0d0-hTime
          else
             interpolantTime=      hTime
          end if
          if (interpolantTime == 0.0d0) cycle
          ! Find the largest mass corresponding to this sigma.
          iBoundLeft =1
          iBoundRight=size(self%rootVarianceUniqueTable(j)%rootVariance)
          do while (iBoundLeft+1 < iBoundRight)
             i=int((iBoundLeft+iBoundRight)/2)
             if (self%rootVarianceUniqueTable(j)%rootVariance(i) < rootVariance) then
                iBoundRight=i
             else
                iBoundLeft =i
             end if
          end do
          i                =self%rootVarianceUniqueTable(j)%index(iBoundRight)
          h                =+(     rootVariance               -self%rootVarianceTable(j)%y(i)) &
               &            /(self%rootVarianceTable(j)%y(i-1)-self%rootVarianceTable(j)%y(i))
          filteredPowerMass=+filteredPowerMass                                    &
               &            +exp(                                                 &
               &                 +log(self%rootVarianceTable(j)%x(i  ))*(1.0d0-h) &
               &                 +log(self%rootVarianceTable(j)%x(i-1))*       h  &
               &                )                                                 &
               &            *interpolantTime
       end do
    end if
    return
  end function filteredPowerMass

  subroutine filteredPowerRetabulate(self,mass,time)
    !% Tabulate the cosmological mass variance.
    use :: Cosmology_Parameters    , only : hubbleUnitsLittleH
    use :: File_Utilities          , only : File_Lock                , File_Unlock                      , lockDescriptor
    use :: Galacticus_Display      , only : Galacticus_Display_Indent, Galacticus_Display_Message       , Galacticus_Display_Unindent, verbosityWorking
    use :: Galacticus_Error        , only : Galacticus_Error_Report  , Galacticus_Warn
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Ranges        , only : Make_Range               , rangeTypeLogarithmic
    use :: Tables                  , only : table1DLogarithmicCSpline, table1DLogarithmicMonotoneCSpline
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout)               :: self
    double precision                                       , intent(in   ), optional     :: mass                      , time
    ! Radius for σ(M) normalization in Mpc/h.
    double precision                                       , parameter                   :: radiusNormalization =8.0d0
    integer                                                                              :: i                         , rootVarianceTableCount , &
         &                                                                                  j                         , rootVarianceUniqueCount, &
         &                                                                                  rootVarianceTimeCount     , k                      , &
         &                                                                                  countNewLower             , countNewUpper
    double precision                                                                     :: sigma                     , smoothingMass          , &
         &                                                                                  massMinimum               , massMaximum            , &
         &                                                                                  timeMinimum               , timeMaximum
    logical                                                , allocatable  , dimension(:) :: rootVarianceIsUnique
    type            (varying_string                       ), save                        :: message
    character       (len=12                               )                              :: label                     , labelLow               , &
         &                                                                                  labelHigh                 , labelTarget
    type            (lockDescriptor                       ), save                        :: fileLock
    ! The variables "message", and "fileLock" are saved (and made threadprivate) as their destructors are expensive, and this
    ! functions gets called a lot.
    !$omp threadprivate(message,fileLock)

    if (self%remakeTable(mass,time)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(self%fileName),fileLock,lockIsShared=.true.)
       call self%fileRead()
       call File_Unlock(fileLock,sync=.false.)
    end if
    if (self%remakeTable(mass,time)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(self%fileName),fileLock,lockIsShared=.false.)
       ! Try again to read the file - another process/thread may have already created the file in which case we may not need to do so again.
       call self%fileRead()
       if (self%remakeTable(mass,time)) then
          ! Compute the mass at which the mass variance is normalized.
          smoothingMass=+(                                                               &
               &          +4.0d0                                                         &
               &          /3.0d0                                                         &
               &          *Pi                                                            &
               &         )                                                               &
               &        *  self%cosmologyParameters_%OmegaMatter    (                  ) &
               &        *  self%cosmologyParameters_%densityCritical(                  ) &
               &        *(                                                               &
               &          +radiusNormalization                                           &
               &          /self%cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH) &
               &         )**3
          ! Determine the normalization of the power spectrum.
          self%sigmaNormalization=+self%sigma8Value                                                                &
               &                  /rootVariance(time_=self%cosmologyFunctions_%cosmicTime(1.0d0),useTopHat=.true.)
          ! Find suitable range of masses to tabulate.
          if (present(mass)) then
             countNewLower=0
             countNewUpper=0
             if (self%initialized) then
                massMinimum     =min(mass/10.0d0,self%massMinimum)
                massMaximum     =max(mass*10.0d0,self%massMaximum)
             else
                self%massMinimum=    mass
                self%massMaximum=    mass
                massMinimum     =    mass/10.0d0
                massMaximum     =    mass*10.0d0
             end if
             ! Determine how many points the table must be extended by in each direction to span the new required range.
             if (self%massMinimum > massMinimum) countNewLower=int(+log10(self%massMinimum/massMinimum)*dble(filteredPowerTablePointsPerDecade)+1.0d0)
             if (self%massMaximum < massMaximum) countNewUpper=int(-log10(self%massMaximum/massMaximum)*dble(filteredPowerTablePointsPerDecade)+1.0d0)
             ! Adjust the limits of the table by an integer number of steps.
             self%massMinimum=self%massMinimum/10.0d0**(dble(countNewLower)/dble(filteredPowerTablePointsPerDecade))
             self%massMaximum=self%massMaximum*10.0d0**(dble(countNewUpper)/dble(filteredPowerTablePointsPerDecade))
          else if (.not.self%initialized) then
             ! No mass was given, but the tables are not initialized. Must provide some mass range.
             self%massMinimum=1.0d10
             self%massMaximum=1.0d15
          end if
          rootVarianceTableCount=int(                                         &
               &                     +log10(                                  &
               &                            +self%massMaximum                 &
               &                            /self%massMinimum                 &
               &                           )                                  &
               &                     *dble(filteredPowerTablePointsPerDecade) &
               &                    )
          ! Find suitable range of times to tabulate.
          if (self%growthIsMassDependent_) then
             if (present(time)) then
                countNewLower=0
                countNewUpper=0
                if (self%initialized) then
                   timeMinimum     =min(time/2.0d0,self%timeMinimum)
                   timeMaximum     =max(time*2.0d0,self%timeMaximum)
                else
                   self%timeMinimum=    time
                   self%timeMaximum=    time
                   timeMinimum     =    time/2.0d0
                   timeMaximum     =    time*2.0d0
                end if
                ! Determine how many points the table must be extended by in each direction to span the new required range.
                if (self%timeMinimum > timeMinimum) countNewLower=int(+log10(self%timeMinimum/timeMinimum)*dble(filteredPowerTimePointsPerDecade)+1.0d0)
                if (self%timeMaximum < timeMaximum) countNewUpper=int(-log10(self%timeMaximum/timeMaximum)*dble(filteredPowerTimePointsPerDecade)+1.0d0)
                ! Adjust the limits of the table by an integer number of steps.
                self%timeMinimum=self%timeMinimum/10.0d0**(dble(countNewLower)/dble(filteredPowerTimePointsPerDecade))
                self%timeMaximum=self%timeMaximum*10.0d0**(dble(countNewUpper)/dble(filteredPowerTimePointsPerDecade))
             else if (.not.self%initialized) then
                ! No time was given, but the tables are not initialized. Must provide some time range.
                self%timeMinimum=self%cosmologyFunctions_%cosmicTime(0.5d0)
                self%timeMaximum=self%cosmologyFunctions_%cosmicTime(1.0d0)
             end if
             rootVarianceTimeCount =int(                                         &
                  &                     +log10(                                  &
                  &                            +self%timeMaximum                 &
                  &                            /self%timeMinimum                 &
                  &                           )                                  &
                  &                     *dble(filteredPowerTimePointsPerDecade ) &
                  &                    )
             self%timeMinimumLogarithmic     =                              log(                 self%timeMinimum)
             self%timeLogarithmicDeltaInverse=dble(rootVarianceTimeCount-1)/log(self%timeMaximum/self%timeMinimum)
          else
             ! Growth of the transferred power spectrum is independent of mass - we can therefore tabulate σ(M) at a single epoch
             ! and use the linear growth factor to transform it to other epochs.
             self%timeMinimum                =self%cosmologyFunctions_%cosmicTime(1.0d0)
             self%timeMaximum                =self%cosmologyFunctions_%cosmicTime(1.0d0)
             rootVarianceTimeCount           =1
             self%timeMinimumLogarithmic     =0.0d0
             self%timeLogarithmicDeltaInverse=0.0d0
          end if
          if (allocated(self%times                  )) deallocate(self%times                  )
          if (allocated(self%rootVarianceUniqueTable)) deallocate(self%rootVarianceUniqueTable)
          allocate(self%rootVarianceUniqueTable(rootVarianceTimeCount))
          allocate(self%times                  (rootVarianceTimeCount))
          if (self%growthIsMassDependent_) then
             self%times=Make_Range(self%timeMinimum,self%timeMaximum,rootVarianceTimeCount,rangeTypeLogarithmic)
          else
             self%times=self%timeMinimum
          end if
          ! Allocate table grid.
          if (allocated(self%rootVarianceTable)) then
             do i=1,size(self%rootVarianceTable)
                call self%rootVarianceTable(i)%destroy()
             end do
             deallocate(self%rootVarianceTable)
          end if
          if (self%monotonicInterpolation) then
             allocate(table1DLogarithmicMonotoneCSpline :: self%rootVarianceTable(rootVarianceTimeCount))
          else
             allocate(table1DLogarithmicCSpline         :: self%rootVarianceTable(rootVarianceTimeCount))
          end if
          call Galacticus_Display_Indent("retabulating σ(M)",verbosityWorking)
          write    (labelLow   ,'(e9.2)') self%massMinimum
          write    (labelHigh  ,'(e9.2)') self%massMaximum
          if (present(mass)) then
             write (labelTarget,'(e9.2)')      mass
          else
             labelTarget="unspecified"
          end if
          call Galacticus_Display_Message("mass range: "//labelLow//" < "//labelTarget//" < "//labelHigh//" M☉" ,verbosityWorking)
          write    (labelLow   ,'(f9.4)') self%timeMinimum
          write    (labelHigh  ,'(f9.4)') self%timeMaximum
          if (present(time)) then
             write (labelTarget,'(f9.4)')      time
          else
             labelTarget="unspecified"
          end if
          call Galacticus_Display_Message("time range: "//labelLow//" < "//labelTarget//" < "//labelHigh//" Gyr",verbosityWorking)
          do k=1,rootVarianceTimeCount
             call self%rootVarianceTable(k)%create(self%massMinimum,self%massMaximum,rootVarianceTableCount)
             allocate(rootVarianceIsUnique(rootVarianceTableCount))
             rootVarianceIsUnique=.true.
             ! Compute σ(M) at each tabulated point.
             massMinimum=-1.0d0
             do i=1,rootVarianceTableCount
                smoothingMass=+self        %rootVarianceTable(k)%x(                                    i)
                sigma        =+rootVariance                       (time_=self%times(k),useTopHat=.false.) &
                     &        *self%sigmaNormalization
                ! Enforce monotonicity.
                if (i > 1) then
                   if (sigma >= self%rootVarianceTable(k)%y(i-1)) then
                      massMinimum            =smoothingMass
                      rootVarianceIsUnique(i)=.false.
                   end if
                   sigma=min(sigma,self%rootVarianceTable(k)%y(i-1))
                end if
                ! Store the value.
                call self%rootVarianceTable(k)%populate(sigma,i,computeSpline=(i == rootVarianceTableCount))
             end do
             ! Find unique values in the variance table.
             rootVarianceUniqueCount=count(rootVarianceIsUnique)
             allocate(self%rootVarianceUniqueTable(k)%rootVariance(rootVarianceUniqueCount))
             allocate(self%rootVarianceUniqueTable(k)%index       (rootVarianceUniqueCount))
             j=1
             do i=1,rootVarianceTableCount
                if (rootVarianceIsUnique(i)) then
                   self%rootVarianceUniqueTable(k)%rootVariance(j)=self%rootVarianceTable(k)%y(i)
                   self%rootVarianceUniqueTable(k)%index       (j)=i
                   j                                              =j+1
                end if
             end do
             deallocate(rootVarianceIsUnique)
             ! Abort or warn if σ(M) has no increase below some mass scale.
             if (massMinimum > 0.0d0) then
                if (self%nonMonotonicIsFatal) then
                   message=         ""
                else
                   message=         "WARNING: "
                end if
                write (label,'(e12.6)') massMinimum
                message=         "σ(M) is non-increasing below mass M="//label//"M☉"
                write (label,'(e12.6)') self%times(k)
                message=message//" at time t="//label//"Gyr."//char(10)
                if (self%nonMonotonicIsFatal) then
                   call Galacticus_Error_Report(message//{introspection:location})
                else
                   message=message//"         If problems occur consider not attempting to model structure below this mass scale."
                   call Galacticus_Warn(message)
                end if
             end if
          end do
          call Galacticus_Display_Unindent("done",verbosityWorking)
          ! Table is now initialized.
          self%initialized=.true.
          ! Store file.
          call self%fileWrite()
       end if
       call File_Unlock(fileLock)
    end if
    return

  contains

    double precision function rootVariance(time_,useTopHat)
      !% Compute the root-variance of mass in spheres enclosing the given {\normalfont \ttfamily mass} from the power spectrum.
      use :: Numerical_Constants_Math, only : Pi
      use :: Numerical_Integration   , only : integrator, GSL_Integ_Gauss15
      implicit none
      double precision            , intent(in   ) :: time_
      logical                     , intent(in   ) :: useTopHat
      double precision            , parameter     :: wavenumberBAO    =5.0d0 ! The wavenumber above which baryon acoustic oscillations are small - used to split the integral allowing the oscillating part to be handled robustly.
      double precision                            :: topHatRadius           , wavenumberMaximum, &
           &                                         wavenumberMinimum      , integrandLow     , &
           &                                         integrandHigh
      type            (integrator)                :: integrator_

      filteredPowerTime=time_
      topHatRadius     =(                                             &
           &             +(                                           &
           &               +3.0d0                                     &
           &               /4.0d0                                     &
           &               /Pi                                        &
           &              )                                           &
           &             *smoothingMass                               &
           &             /self%cosmologyParameters_%OmegaMatter    () &
           &             /self%cosmologyParameters_%densityCritical() &
           &            )**(1.0d0/3.0d0)
      wavenumberMinimum=0.0d0
      ! The integral over the power spectrum is split at a wavenumber corresponding to the smallest scale at which BAO features
      ! are significant (unless the upper limit of the integral is already below that wavenumber). This allows the oscillatory
      ! part of the integral to be computed more accurately, without affecting the non-oscillatory part at larger wavenumbers, and
      ! leads to an overall more accurate and robust determination of σ(M).
      if (useTopHat) then
         integrator_=integrator(varianceIntegrandTopHat,toleranceRelative=+self%tolerance,integrationRule=GSL_Integ_Gauss15)
         wavenumberMaximum=min(1.0d3/topHatRadius,self%powerSpectrumWindowFunctionTopHat_%wavenumberMaximum(smoothingMass))
         if (wavenumberMaximum > wavenumberBAO) then
            integrandLow =   integrator_%integrate(wavenumberMinimum,wavenumberBAO    )
            integrandHigh=   integrator_%integrate(wavenumberBAO    ,wavenumberMaximum)
            rootVariance =+(                                                            &
                 &          +integrandLow                                               &
                 &          +integrandHigh                                              &
                 &         )                                                            &
                 &        /2.0d0                                                        &
                 &        /Pi**2
         else
            rootVariance =+  integrator_%integrate(wavenumberMinimum,wavenumberMaximum) &
                 &       /2.0d0                                                         &
                 &       /Pi**2
         end if
      else
         integrator_=integrator(varianceIntegrand,toleranceRelative=+self%tolerance,integrationRule=GSL_Integ_Gauss15)
         wavenumberMaximum=min(1.0d3/topHatRadius,self%powerSpectrumWindowFunction_      %wavenumberMaximum(smoothingMass))
         if (wavenumberMaximum > wavenumberBAO) then
            integrandLow =   integrator_%integrate(wavenumberMinimum,wavenumberBAO    )
            integrandHigh=   integrator_%integrate(wavenumberBAO    ,wavenumberMaximum)
            rootVariance =+(                                                            &
                 &          +integrandLow                                               &
                 &          +integrandHigh                                              &
                 &         )                                                            &
                 &        /2.0d0                                                        &
                 &        /Pi**2
         else
            rootVariance =+  integrator_%integrate(wavenumberMinimum,wavenumberMaximum) &
                 &        /2.0d0                                                        &
                 &        /Pi**2
         end if
      end if
      rootVariance=sqrt(rootVariance)
      return
    end function rootVariance

    double precision function varianceIntegrand(wavenumber)
      !% Integrand function used in computing the variance in (real space) top-hat spheres from the power spectrum.
      implicit none
      double precision, intent(in   ) :: wavenumber

      ! Return power spectrum multiplied by window function and volume element in k-space. Factors of 2 and π are included
      ! elsewhere.
      varianceIntegrand=+  self%powerSpectrumPrimordialTransferred_%power(wavenumber,filteredPowerTime) &
           &            *(                                                                              &
           &              +self%powerSpectrumWindowFunction_       %value(wavenumber,smoothingMass    ) &
           &              *                                               wavenumber                    &
           &             )**2
      return
    end function varianceIntegrand

    double precision function varianceIntegrandTopHat(wavenumber)
      !% Integrand function used in computing the variance in (real space) top-hat spheres from the power spectrum.
      implicit none
      double precision, intent(in   ) :: wavenumber

      ! Return power spectrum multiplied by window function and volume element in k-space. Factors of 2 and π are included
      ! elsewhere.
      varianceIntegrandTopHat=+  self%powerSpectrumPrimordialTransferred_%power(wavenumber,filteredPowerTime) &
           &                  *(                                                                              &
           &                    +self%powerSpectrumWindowFunctionTopHat_ %value(wavenumber,smoothingMass    ) &
           &                    *                                               wavenumber                    &
           &                   )**2
      return
    end function varianceIntegrandTopHat

  end subroutine filteredPowerRetabulate

  subroutine filteredPowerInterpolantsTime(self,time,i,h)
    !% Compute interoplants in time.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    integer                                                , intent(  out) :: i
    double precision                                       , intent(  out) :: h

    h=(log(time)-self%timeMinimumLogarithmic)*self%timeLogarithmicDeltaInverse+1.0d0
    i=  int (h)
    h=h-dble(i)
    if (i == size(self%times)) then
       ! Requested time must exactly equal the maximum tabulated time.
       i=size(self%times)-1
       h=1.0d0
    else if (i < 1) then
       call Galacticus_Error_Report('interpolant out of range'//{introspection:location})
    end if
    return
  end subroutine filteredPowerInterpolantsTime

  logical function filteredPowerGrowthIsMassDependent(self)
    !% Return true if the growth rate of the variance is mass-dependent.
    implicit none
    class(cosmologicalMassVarianceFilteredPower), intent(inout) :: self

    filteredPowerGrowthIsMassDependent=self%growthIsMassDependent_
    return
  end function filteredPowerGrowthIsMassDependent

  subroutine filteredPowerFileRead(self)
    !% Read tabulated data on mass variance from file.
    use :: File_Utilities    , only : File_Exists
    use :: Galacticus_Display, only : Galacticus_Display_Message, verbosityWorking
    use :: IO_HDF5           , only : hdf5Access                , hdf5Object
    use :: Tables            , only : table1DLogarithmicCSpline , table1DLogarithmicMonotoneCSpline
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout)               :: self
    double precision                                       , dimension(:  ), allocatable :: massTmp        , timesTmp
    double precision                                       , dimension(:,:), allocatable :: rootVarianceTmp, rootVarianceUniqueTmp
    integer                                                , dimension(:  ), allocatable :: uniqueSizeTmp
    integer                                                , dimension(:,:), allocatable :: indexTmp
    type            (hdf5Object                           )                              :: dataFile
    integer                                                                              :: i

    ! Return immediately if the file does not exist.
    if (.not.File_Exists(char(self%fileName))) return
    call Galacticus_Display_Message('reading σ(M) data from: '//self%fileName,verbosityWorking)
    !$ call hdf5Access%set()
    call dataFile%openFile     (char(self%fileName)          ,overWrite                       =.false.)
    call dataFile%readDataset  ('times'                      ,     timesTmp                           )
    call dataFile%readDataset  ('mass'                       ,     massTmp                            )
    call dataFile%readDataset  ('rootVariance'               ,     rootVarianceTmp                    )
    call dataFile%readDataset  ('rootVarianceUnique'         ,     rootVarianceUniqueTmp              )
    call dataFile%readDataset  ('indexUnique'                ,     indexTmp                           )
    call dataFile%readDataset  ('uniqueSize'                 ,     uniqueSizeTmp                      )
    call dataFile%readAttribute('sigmaNormalization'         ,self%sigmaNormalization                 )
    call dataFile%readAttribute('massMinimum'                ,self%massMinimum                        )
    call dataFile%readAttribute('massMaximum'                ,self%massMaximum                        )
    call dataFile%readAttribute('timeMinimum'                ,self%timeMinimum                        )
    call dataFile%readAttribute('timeMaximum'                ,self%timeMaximum                        )
    call dataFile%readAttribute('timeMinimumLogarithmic'     ,self%timeMinimumLogarithmic             )
    call dataFile%readAttribute('timeLogarithmicDeltaInverse',self%timeLogarithmicDeltaInverse        )
    call dataFile%close        (                                                                      )
    !$ call hdf5Access%unset()
    if (allocated(self%times                  )) deallocate(self%times                  )
    if (allocated(self%rootVarianceTable      )) deallocate(self%rootVarianceTable      )
    if (allocated(self%rootVarianceUniqueTable)) deallocate(self%rootVarianceUniqueTable)
    allocate(self%times                  (size(timesTmp)))
    allocate(self%rootVarianceUniqueTable(size(timesTmp)))
    if (self%monotonicInterpolation) then
       allocate(table1DLogarithmicMonotoneCSpline :: self%rootVarianceTable(size(timesTmp)))
    else
       allocate(table1DLogarithmicCSpline         :: self%rootVarianceTable(size(timesTmp)))
    end if
    self%times=timesTmp
    do i=1,size(self%times)
       allocate(self%rootVarianceUniqueTable(i)%rootVariance(uniqueSizeTmp(i)))
       allocate(self%rootVarianceUniqueTable(i)%index       (uniqueSizeTmp(i)))
       call self%rootVarianceTable(i)%create  (self%massMinimum,self%massMaximum,size(massTmp))
       call self%rootVarianceTable(i)%populate(rootVarianceTmp(:,i))
       self%rootVarianceUniqueTable(i)%rootVariance=rootVarianceUniqueTmp(1:uniqueSizeTmp(i),i)
       self%rootVarianceUniqueTable(i)%index       =indexTmp             (1:uniqueSizeTmp(i),i)
    end do
    deallocate(rootVarianceTmp      )
    deallocate(rootVarianceUniqueTmp)
    deallocate(indexTmp             )
    deallocate(massTmp              )
    deallocate(timesTmp             )
    deallocate(uniqueSizeTmp        )
    self%initialized=.true.
    return
  end subroutine filteredPowerFileRead

  subroutine filteredPowerFileWrite(self)
    !% Write tabulated data on mass variance to file.
    use :: Galacticus_Display, only : Galacticus_Display_Message, verbosityWorking
    use :: HDF5              , only : hsize_t
    use :: IO_HDF5           , only : hdf5Access                , hdf5Object
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout)               :: self
    double precision                                       , dimension(:  ), allocatable :: massTmp
    double precision                                       , dimension(:,:), allocatable :: rootVarianceTmp, rootVarianceUniqueTmp
    integer                                                , dimension(:  ), allocatable :: uniqueSizeTmp
    integer                                                , dimension(:,:), allocatable :: indexTmp
    type            (hdf5Object                           )                              :: dataFile
    integer                                                                              :: i

    call Galacticus_Display_Message('writing σ(M) data to: '//self%fileName,verbosityWorking)
    ! Prepare data.
    allocate(massTmp              (self%rootVarianceTable(1)%size()                 ))
    allocate(rootVarianceTmp      (self%rootVarianceTable(1)%size(),size(self%times)))
    allocate(rootVarianceUniqueTmp(self%rootVarianceTable(1)%size(),size(self%times)))
    allocate(indexTmp             (self%rootVarianceTable(1)%size(),size(self%times)))
    allocate(uniqueSizeTmp        (                                 size(self%times)))
    rootVarianceUniqueTmp   ( :                , : )=     0.0d0
    indexTmp                ( :                , : )=     0
    massTmp                 ( :                    )=     self%rootVarianceTable      (1)%xs          ()
    do i=1,size(self%times)
       rootVarianceTmp      ( :                ,i:i)=     self%rootVarianceTable      (i)%ys          ()
       uniqueSizeTmp        (                   i  )=size(self%rootVarianceUniqueTable(i)%index         )
       rootVarianceUniqueTmp(1:uniqueSizeTmp(i),i  )=     self%rootVarianceUniqueTable(i)%rootVariance
       indexTmp             (1:uniqueSizeTmp(i),i  )=     self%rootVarianceUniqueTable(i)%index
    end do
    ! Open the data file.
    !$ call hdf5Access%set()
    call dataFile%openFile      (char(self%fileName)             ,overWrite                    =.true.,objectsOverwritable=.true.,chunkSize=100_hsize_t,compressionLevel=9)
    call dataFile%writeDataset  (self%times                      ,'times'                                                                                                 )
    call dataFile%writeDataset  (     massTmp                    ,'mass'                                                                                                  )
    call dataFile%writeDataset  (     rootVarianceTmp            ,'rootVariance'                                                                                          )
    call dataFile%writeDataset  (     rootVarianceUniqueTmp      ,'rootVarianceUnique'                                                                                    )
    call dataFile%writeDataset  (     indexTmp                   ,'indexUnique'                                                                                           )
    call dataFile%writeDataset  (     uniqueSizeTmp              ,'uniqueSize'                                                                                            )
    call dataFile%writeAttribute(self%sigmaNormalization         ,'sigmaNormalization'                                                                                    )
    call dataFile%writeAttribute(self%massMinimum                ,'massMinimum'                                                                                           )
    call dataFile%writeAttribute(self%massMaximum                ,'massMaximum'                                                                                           )
    call dataFile%writeAttribute(self%timeMinimum                ,'timeMinimum'                                                                                           )
    call dataFile%writeAttribute(self%timeMaximum                ,'timeMaximum'                                                                                           )
    call dataFile%writeAttribute(self%timeMinimumLogarithmic     ,'timeMinimumLogarithmic'                                                                                )
    call dataFile%writeAttribute(self%timeLogarithmicDeltaInverse,'timeLogarithmicDeltaInverse'                                                                           )
    call dataFile%close         (                                                                                                                                         )
    !$ call hdf5Access%unset()
    deallocate(rootVarianceTmp      )
    deallocate(rootVarianceUniqueTmp)
    deallocate(indexTmp             )
    deallocate(massTmp              )
    deallocate(uniqueSizeTmp        )
    return
  end subroutine filteredPowerFileWrite

  logical function filteredPowerRemakeTable(self,mass,time)
    !% Determine if the table should be remade.
    implicit none
    class           (cosmologicalMassVarianceFilteredPower  ), intent(inout)           :: self
    double precision                                         , intent(in   ), optional :: mass, time

    if (self%initialized) then
       filteredPowerRemakeTable=.false.
       if (present(mass)) then
          filteredPowerRemakeTable=(                          &
               &                     mass < self%massMinimum  &
               &                    .or.                      &
               &                     mass > self%massMaximum  &
               &                   )
       end if
       if (present(time).and.self%growthIsMassDependent_) then
          filteredPowerRemakeTable=(                          &
               &                     filteredPowerRemakeTable &
               &                    .or.                      &
               &                     time < self%timeMinimum  &
               &                    .or.                      &
               &                     time > self%timeMaximum  &
               &                   )
       end if
    else
       filteredPowerRemakeTable=.true.
    end if
    return
  end function filteredPowerRemakeTable
