!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of cosmological density field mass variance computed using a filtered power spectrum.

  !# <cosmologicalMassVariance name="cosmologicalMassVarianceFilteredPower" defaultThreadPrivate="yes">
  !#  <description>Mass variance of cosmological density fields computed from a filtered power spectrum.</description>
  !# </cosmologicalMassVariance>
  use Tables

  type, extends(cosmologicalMassVarianceClass) :: cosmologicalMassVarianceFilteredPower
     !% A cosmological mass variance class computing variance from a filtered power spectrum.
     private
     logical                                     :: initialized
     double precision                            :: tolerance        , toleranceTopHat   , &
          &                                         sigma8Value      , sigmaNormalization, &
          &                                         massMinimum      , massMaximum
     type            (table1DLogarithmicCSpline) :: rootVarianceTable
   contains
     !@ <objectMethods>
     !@   <object>cosmologicalMassVarianceFilteredPower</object>
     !@   <objectMethod>
     !@     <method>retabulate</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ time\argin</arguments>
     !@     <description>Tabulate cosmological mass variance.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                       filteredPowerDestructor
     procedure :: stateStore                         => filteredPowerStateStore
     procedure :: stateRestore                       => filteredPowerStateRestore
     procedure :: sigma8                             => filteredPowerSigma8
     procedure :: powerNormalization                 => filteredPowerPowerNormalization
     procedure :: rootVariance                       => filteredPowerRootVariance
     procedure :: rootVarianceLogarithmicGradient    => filteredPowerRootVarianceLogarithmicGradient
     procedure :: rootVarianceAndLogarithmicGradient => filteredPowerRootVarianceAndLogarithmicGradient
     procedure :: mass                               => filteredPowerMass
     procedure :: retabulate                         => filteredPowerRetabulate
  end type cosmologicalMassVarianceFilteredPower

  interface cosmologicalMassVarianceFilteredPower
     !% Constructors for the {\normalfont \ttfamily filteredPower} cosmological mass variance class.
     module procedure filteredPowerConstructorParameters
     module procedure filteredPowerConstructorInternal
  end interface cosmologicalMassVarianceFilteredPower

  ! Number of points per decade to use in tabulation of σ(M).
  integer, parameter :: filteredPowerTablePointsPerDecade=10

contains

  function filteredPowerConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily filteredPower} cosmological mass variance class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type            (cosmologicalMassVarianceFilteredPower)                :: filteredPowerConstructorParameters
    type            (inputParameters                      ), intent(in   ) :: parameters
    double precision                                                       :: sigma_8                           , tolerance, &
         &                                                                    toleranceTopHat
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>sigma_8</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.817d0</defaultValue>
    !#   <defaultSource>(\citealt{hinshaw_nine-year_2012}; CMB$+H_0+$BAO)</defaultSource>
    !#   <description>The fractional mass fluctuation in the linear density field at the present day in spheres of radius 8~Mpc/h.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>toleranceTopHat</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d-6</defaultValue>
    !#   <description>The relative tolerance to use in integrating over the linear power spectrum using a top-hat (real space) window function to compute the cosmological mass variance.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>tolerance</name>
    !#   <source>parameters</source>
    !#   <defaultValue>4.0d-6</defaultValue>
    !#   <description>The relative tolerance to use in integrating over the linear power spectrum to compute the cosmological mass variance.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    ! Construct the instance.
    filteredPowerConstructorParameters=filteredPowerConstructorInternal(sigma_8,tolerance,toleranceTopHat)
    return
  end function filteredPowerConstructorParameters

  function filteredPowerConstructorInternal(sigma8,tolerance,toleranceTopHat)
    !% Internal constructor for the {\normalfont \ttfamily filteredPower} linear growth class.
    implicit none
    type            (cosmologicalMassVarianceFilteredPower)                :: filteredPowerConstructorInternal
    double precision                                       , intent(in   ) :: tolerance                       , toleranceTopHat, &
         &                                                                    sigma8

    filteredPowerConstructorInternal%sigma8Value    =sigma8
    filteredPowerConstructorInternal%tolerance      =tolerance
    filteredPowerConstructorInternal%toleranceTopHat=toleranceTopHat
    filteredPowerConstructorInternal%initialized    =.false.
    filteredPowerConstructorInternal%massMinimum    =1.0d06
    filteredPowerConstructorInternal%massMaximum    =1.0d15
    return
  end function filteredPowerConstructorInternal

  subroutine filteredPowerDestructor(self)
    !% Destructor for the {\normalfont \ttfamily filteredPower} linear growth class.
    implicit none
    type (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    
    if (self%initialized) call self%rootVarianceTable%destroy()
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
  
  double precision function filteredPowerRootVariance(self,mass)
    !% Return the root-variance of the cosmological density field in a spherical region containing the given {\normalfont
    !% \ttfamily mass} on average.
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: mass

    call self%retabulate(mass)
    filteredPowerRootVariance=self%rootVarianceTable%interpolate(mass)
    return
  end function filteredPowerRootVariance
  
  double precision function filteredPowerRootVarianceLogarithmicGradient(self,mass)
    !% Return the logairhtmic gradient with respect to mass of the root-variance of the cosmological density field in a spherical
    !% region containing the given {\normalfont \ttfamily mass} on average.
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: mass

    call self%retabulate(mass)
    filteredPowerRootVarianceLogarithmicGradient=+self%rootVarianceTable%interpolateGradient(mass) &
         &                                       /self%rootVarianceTable%interpolate        (mass) &
         &                                       *                                           mass
    return
  end function filteredPowerRootVarianceLogarithmicGradient
  
  subroutine filteredPowerRootVarianceAndLogarithmicGradient(self,mass,rootVariance,rootVarianceLogarithmicGradient)
    !% Return the value and logairhtmic gradient with respect to mass of the root-variance of the cosmological density field in a
    !% spherical region containing the given {\normalfont \ttfamily mass} on average.
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: mass
    double precision                                       , intent(  out) :: rootVariance, rootVarianceLogarithmicGradient

    call self%retabulate(mass)
    rootVariance                   =+self%rootVarianceTable%interpolate        (mass)
    rootVarianceLogarithmicGradient=+self%rootVarianceTable%interpolateGradient(mass) &
         &                          /     rootVariance                                &
         &                          *                                           mass
    return
  end subroutine filteredPowerRootVarianceAndLogarithmicGradient

  double precision function filteredPowerMass(self,rootVariance)
    !% Return the mass corrresponding to the given {\normalfont \ttfamily } root-variance of the cosmological density field.
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: rootVariance
    double precision                                                       :: h
    integer                                                                :: i
    
    ! If the requested root-variance is below the lowest value tabulated, attempt to tabulate to higher mass (lower
    ! root-variance).
    if (.not.self%initialized) call self%retabulate()
    do while (rootVariance < self%rootVarianceTable%y(-1))
       call self%retabulate(self%rootVarianceTable%x(-1)*2.0d0)
    end do
    ! If sigma exceeds the highest value tabulated, simply return the lowest tabulated mass.
    if (rootVariance > self%rootVarianceTable%y(1)) then
       filteredPowerMass=self%rootVarianceTable%x(1)
    else
       ! Find the largest mass corresponding to this sigma.
       i=self%rootVarianceTable%size()
       do while (i > 1 .and. self%rootVarianceTable%y(i-1) < rootVariance)
          i=i-1
       end do
       h                =+(     rootVariance            -self%rootVarianceTable%y(i)) &
            &            /(self%rootVarianceTable%y(i-1)-self%rootVarianceTable%y(i))
       filteredPowerMass=exp(                                              &
            &                +log(self%rootVarianceTable%x(i  ))*(1.0d0-h) &
            &                +log(self%rootVarianceTable%x(i-1))*       h  &
            &               )
    end if
    return
  end function filteredPowerMass
 
  subroutine filteredPowerRetabulate(self,mass)
    !% Tabulate the cosmological mass variance.
    use Numerical_Constants_Math
    use Cosmology_Parameters
    use Power_Spectra_Primordial_Transferred
    implicit none
    class           (cosmologicalMassVarianceFilteredPower  ), intent(inout)           :: self
    double precision                                         , intent(in   ), optional :: mass
    class           (cosmologyParametersClass               ), pointer                 :: cosmologyParameters_
    class           (powerSpectrumPrimordialTransferredClass), pointer                 :: powerSpectrumPrimordialTransferred_
    ! Radius for σ(M) normalization in Mpc/h.
    double precision                                         , parameter               :: radiusNormalization                =8.0d0
    integer                                                                            :: i                                        , rootVarianceTableCount
    double precision                                                                   :: sigma                                    , smoothingMass
    logical                                                                            :: remakeTable
    
    ! Check if we need to recompute our table.
    if (self%initialized) then
       if (present(mass)) then
          remakeTable=(                         &
               &        mass < self%massMinimum &
               &       .or.                     &
               &        mass > self%massMaximum &
               &      )
       else
          remakeTable=.false.
       end if
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       ! Get required objects.
       cosmologyParameters_                => cosmologyParameters               ()
       powerSpectrumPrimordialTransferred_ => powerSpectrumPrimordialTransferred()
       ! Compute the mass at which the mass variance is normalized.
       smoothingMass=+(                                                          &
            &          +4.0d0                                                    &
            &          /3.0d0                                                    &
            &          *Pi                                                       &
            &         )                                                          &
            &        *  cosmologyParameters_%OmegaMatter    (                  ) &
            &        *  cosmologyParameters_%densityCritical(                  ) &
            &        *(                                                          &
            &          +radiusNormalization                                      &
            &          /cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH) &
            &         )**3
       ! Determine the normalization of the power spectrum.
       self%sigmaNormalization=+self%sigma8Value               &
            &                  /rootVariance(useTopHat=.true.)
       ! Find suitable range of masses to tabulate.
       if (present(mass)) then
          self%massMinimum=min(self%massMinimum,mass/10.0d0)
          self%massMaximum=max(self%massMaximum,mass*10.0d0)
       end if
       rootVarianceTableCount=int(                                         &
            &                     +log10(                                  &
            &                            +self%massMaximum                 &
            &                            /self%massMinimum                 &
            &                           )                                  &
            &                     *dble(filteredPowerTablePointsPerDecade) &
            &                    )
       ! Allocate table grid.
       call self%rootVarianceTable%destroy(                                                        )
       call self%rootVarianceTable%create (self%massMinimum,self%massMaximum,rootVarianceTableCount)
       ! Compute sigma(M) at each tabulated point.
       do i=1,rootVarianceTableCount
          smoothingMass=+self        %rootVarianceTable%x(                i)
          sigma        =+rootVariance                    (useTopHat=.false.) &
               &        *self%sigmaNormalization
          ! Enforce monotonicity.
          if (i > 1) sigma=min(sigma,self%rootVarianceTable%y(i-1))
          ! Store the value.
          call self%rootVarianceTable%populate(sigma,i,computeSpline=(i == rootVarianceTableCount))
       end do
       ! Table is now initialized.
       self%initialized=.true.
    end if
    return
    
  contains
  
    double precision function rootVariance(useTopHat)
      !% Compute the root-variance of mass in spheres enclosing the given {\normalfont \ttfamily mass} from the power spectrum.
      use, intrinsic :: ISO_C_Binding
      use               Numerical_Constants_Math
      use               Numerical_Integration
      use               Power_Spectrum_Window_Functions
      implicit none
      logical                                     , intent(in   ) :: useTopHat
      double precision                                            :: topHatRadius        , wavenumberMaximum, &
           &                                                         wavenumberMinimum
      type            (c_ptr                     )                :: parameterPointer
      type            (fgsl_function             )                :: integrandFunction
      type            (fgsl_integration_workspace)                :: integrationWorkspace
      
      topHatRadius     =(                                        &
           &             +(                                      &
           &               +3.0d0                                &
           &               /4.0d0                                &
           &               /Pi                                   &
           &              )                                      &
           &             *smoothingMass                          &
           &             /cosmologyParameters_%OmegaMatter    () &
           &             /cosmologyParameters_%densityCritical() &
           &            )**(1.0d0/3.0d0)
      wavenumberMinimum=    0.0d0/topHatRadius
      wavenumberMaximum=min(1.0d3/topHatRadius,Power_Spectrum_Window_Function_Wavenumber_Maximum(smoothingMass))
      if (useTopHat) then
         rootVariance=+Integrate(                                              &
              &                  wavenumberMinimum                           , &
              &                  wavenumberMaximum                           , &
              &                  varianceIntegrandTopHat                     , &
              &                  parameterPointer                            , &
              &                  integrandFunction                           , &
              &                  integrationWorkspace                        , &
              &                  toleranceAbsolute      =0.0d0               , &
              &                  toleranceRelative      =self%toleranceTopHat, &
              &                  integrationRule        =FGSL_Integ_Gauss15    &
              &                 )                                              &
              &       /2.0d0                                                   &
              &       /Pi**2
      else
         rootVariance=+Integrate(                                              &
              &                  wavenumberMinimum                           , &
              &                  wavenumberMaximum                           , &
              &                  varianceIntegrand                           , &
              &                  parameterPointer                            , &
              &                  integrandFunction                           , &
              &                  integrationWorkspace                        , &
              &                  toleranceAbsolute      =0.0d0               , &
              &                  toleranceRelative      =self%tolerance      , &
              &                  integrationRule        =FGSL_Integ_Gauss15    &
              &                  )                                             &
              &       /2.0d0                                                   &
              &       /Pi**2
      end if
      call Integrate_Done(integrandFunction,integrationWorkspace)
      rootVariance=sqrt(rootVariance)
      return
    end function rootVariance
    
    function varianceIntegrand(wavenumber,parameterPointer) bind(c)
      !% Integrand function used in compute the variance in (real space) top-hat spheres from the power spectrum.
      use, intrinsic :: ISO_C_Binding
      use               Power_Spectrum_Window_Functions
      implicit none
      real(kind=c_double)        :: varianceIntegrand
      real(kind=c_double), value :: wavenumber
      type(     c_ptr   ), value :: parameterPointer

      ! Return power spectrum multiplied by window function and volume element in k-space. Factors of 2 and Pi are included
      ! elsewhere.
      varianceIntegrand=+  powerSpectrumPrimordialTransferred_%power(wavenumber              ) &
           &            *(                                                                     &
           &              +Power_Spectrum_Window_Function           (wavenumber,smoothingMass) &
           &              *                                          wavenumber                &
           &             )**2
      return
    end function varianceIntegrand
    
    function varianceIntegrandTopHat(wavenumber,parameterPointer) bind(c)
      !% Integrand function used in compute the variance in (real space) top-hat spheres from the power spectrum.
      use, intrinsic :: ISO_C_Binding
      use               Power_Spectrum_Window_Functions_Top_Hat
      implicit none
      real(kind=c_double)        :: varianceIntegrandTopHat
      real(kind=c_double), value :: wavenumber
      type(     c_ptr   ), value :: parameterPointer
      
      ! Return power spectrum multiplied by window function and volume element in k-space. Factors of 2 and Pi are included
      ! elsewhere.
      varianceIntegrandTopHat=+  powerSpectrumPrimordialTransferred_   %power(wavenumber              ) &
           &                  *(                                                                        &
           &                    +Power_Spectrum_Window_Function_Top_Hat      (wavenumber,smoothingMass) &
           &                    *                                             wavenumber                &
           &                   )**2
      return
    end function varianceIntegrandTopHat

  end subroutine filteredPowerRetabulate
      
  subroutine filteredPowerStateStore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use FGSL
    implicit none
    class  (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    integer                                       , intent(in   ) :: stateFile
    type   (fgsl_file                            ), intent(in   ) :: fgslStateFile

    write (stateFile) self%massMinimum,self%massMaximum
    return
  end subroutine filteredPowerStateStore

  subroutine filteredPowerStateRestore(self,stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use FGSL
    implicit none
    class  (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    integer                                       , intent(in   ) :: stateFile
    type   (fgsl_file                            ), intent(in   ) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) self%massMinimum,self%massMaximum
    self%initialized=.false.
    call self%retabulate()
    return
  end subroutine filteredPowerStateRestore
