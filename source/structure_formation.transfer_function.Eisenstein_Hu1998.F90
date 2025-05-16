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
  Implements a transfer function class using the fitting function of \cite{eisenstein_baryonic_1998}.
  !!}

  use :: Cosmology_Functions  , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters , only : cosmologyParametersClass
  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !![
  <transferFunction name="transferFunctionEisensteinHu1998">
   <description>
    Provides the \cite{eisenstein_baryonic_1998} fitting function to the transfer function. The effective number of neutrino
    species and the summed mass (in electron volts) of all neutrino species are specified via the {\normalfont \ttfamily
   neutrinoNumberEffective} and {\normalfont \ttfamily neutrinoMassSummed} parameters respectively.
   </description>
  </transferFunction>
  !!]
  type, extends(transferFunctionClass) :: transferFunctionEisensteinHu1998
     !!{
     The {\normalfont \ttfamily eisensteinHu1998} transfer function class.
     !!}
     private
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_      => null()
     class           (darkMatterParticleClass ), pointer :: darkMatterParticle_      => null()
     double precision                                    :: temperatureCMB27                  , fractionBaryons            , &
          &                                                 transferFunctionBaryons           , transferFunctionDarkMatter , &
          &                                                 DtransferFunctionBaryons          , DtransferFunctionDarkMatter, &
          &                                                 omegaMatter                       , omegaDarkMatter            , &
          &                                                 alphaDarkMatter                   , betaDarkMatter             , &
          &                                                 alphaBaryons                      , betaBaryons                , &
          &                                                 wavenumberEquality                , wavenumberSilk             , &
          &                                                 densityMatter                     , densityBaryons             , &
          &                                                 bNode                             , time                       , &
          &                                                 distanceSoundHorizon              , wavenumberPrevious
   contains
     !![
     <methods>
       <method description="Compute common factors needed by \cite{eisenstein_baryonic_1998} transfer function calculations." method="computeFactors" />
     </methods>
     !!]
     final     ::                          eisensteinHu1998Destructor
     procedure :: value                 => eisensteinHu1998Value
     procedure :: logarithmicDerivative => eisensteinHu1998LogarithmicDerivative
     procedure :: computeFactors        => eisensteinHu1998ComputeFactors
     procedure :: halfModeMass          => eisensteinHu1998HalfModeMass
     procedure :: quarterModeMass       => eisensteinHu1998QuarterModeMass
     procedure :: epochTime             => eisensteinHu1998EpochTime
  end type transferFunctionEisensteinHu1998

  interface transferFunctionEisensteinHu1998
     !!{
     Constructors for the \refClass{transferFunctionEisensteinHu1998} transfer function class.
     !!}
     module procedure eisensteinHu1998ConstructorParameters
     module procedure eisensteinHu1998ConstructorInternal
  end interface transferFunctionEisensteinHu1998

contains

  function eisensteinHu1998ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{transferFunctionEisensteinHu1998} transfer function class
    which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (transferFunctionEisensteinHu1998)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyParametersClass        ), pointer       :: cosmologyParameters_
    class           (darkMatterParticleClass         ), pointer       :: darkMatterParticle_
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    ! Call the internal constructor.
    self=transferFunctionEisensteinHu1998(darkMatterParticle_,cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="darkMatterParticle_" />
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function eisensteinHu1998ConstructorParameters

  function eisensteinHu1998ConstructorInternal(darkMatterParticle_,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{transferFunctionEisensteinHu1998} transfer function class.
    !!}
    use :: Cosmology_Parameters , only : hubbleUnitsLittleH
    use :: Dark_Matter_Particles, only : darkMatterParticleCDM
    use :: Error                , only : Error_Report
    implicit none
    type            (transferFunctionEisensteinHu1998)                        :: self
    class           (darkMatterParticleClass         ), intent(in   ), target :: darkMatterParticle_
    class           (cosmologyParametersClass        ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    double precision                                                          :: redshiftEquality                      , redshiftComptonDrag                       , &
         &                                                                       b1                                    , b2                                        , &
         &                                                                       a1                                    , a2                                        , &
         &                                                                       b1Drag                                , b2Drag                                    , &
         &                                                                       y                                     , Gy                                        , &
         &                                                                       ratioBaryonToPhotonMomentumDensityDrag, ratioBaryonToPhotonMomentumDensityEquality
    !![
    <constructorAssign variables="*darkMatterParticle_, *cosmologyParameters_, *cosmologyFunctions_"/>
    !!]

    ! Require that the dark matter be cold dark matter.
    select type (darkMatterParticle_)
       class is (darkMatterParticleCDM)
          ! Cold dark matter particle - this is as expected.
       class default
       call Error_Report('transfer function expects a cold dark matter particle'//{introspection:location})
    end select
    ! Compute the epoch - the transfer function is assumed to be for z=0.
    self%time=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(0.0d0))
    ! Present day CMB temperature [in units of 2.7K].
    self%temperatureCMB27                     =+self%cosmologyParameters_%temperatureCMB() &
         &                                     /2.7d0
    ! Compute various density factors.
    self%omegaMatter                          =+self%cosmologyParameters_%OmegaMatter()
    self%omegaDarkMatter                      =+self%cosmologyParameters_%OmegaMatter()-self%cosmologyParameters_%OmegaBaryon   (                  )
    self%fractionBaryons                      =+self%cosmologyParameters_%OmegaBaryon()/self%cosmologyParameters_%OmegaMatter   (                  )
    self%densityMatter                        =+self%cosmologyParameters_%OmegaMatter()*self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2    
    self%densityBaryons                       =+self%cosmologyParameters_%OmegaBaryon()*self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2    
    ! Redshift of matter-radiation equality.
    redshiftEquality                          =+2.50d4                   &
         &                                     *self%densityMatter       &
         &                                     /self%temperatureCMB27**4
    ! Wavenumber scale imprinted by matter-radiation equality.
    self%wavenumberEquality                   =+7.46d-2                  &
         &                                     *self%densityMatter       &
         &                                     /self%temperatureCMB27**2
    ! Compute redshift at which baryons are released from Compton drag of photons (eqn. 2)
    b1Drag                                    =+  0.313d0                     &
         &                                     /  self%densityMatter**0.419d0 &
         &                                     *(                             &
         &                                       +1.000d0                     &
         &                                       +0.607d0                     &
         &                                       *self%densityMatter**0.674d0 &
         &                                      )
    b2Drag                                    =+0.238d0                       &
         &                                     *  self%densityMatter**0.223d0
    redshiftComptonDrag                       =+1291.0d0                      &
         &                                     *  self%densityMatter**0.251d0 &
         &                                     *(                             &
         &                                       +1.0d0                       &
         &                                       +b1Drag                      &
         &                                       *self%densityBaryons**b2Drag &
         &                                      )                             &
         &                                     /(                             &
         &                                       +1.000d0                     &
         &                                       +0.659d0                     &
         &                                       *self%densityMatter**0.828d0 &
         &                                      )
    ! Compute ratio of baryon to photon momentum density at redshifts of interest.
    ratioBaryonToPhotonMomentumDensityEquality=+31.5d0*self%densityBaryons/self%temperatureCMB27**4/(redshiftEquality   /1.0d3)
    ratioBaryonToPhotonMomentumDensityDrag    =+31.5d0*self%densityBaryons/self%temperatureCMB27**4/(redshiftComptonDrag/1.0d3)
    ! Compute the comoving distance that a sound wave can propagate prior to redshiftComptonDrag (i.e. sound horizon; eq. 6)
    self%distanceSoundHorizon                 =+2.0d0                                                                                                  &
         &                                     /3.0d0                                                                                                  &
         &                                     /self%wavenumberEquality                                                                                &
         &                                     *sqrt(6.0d0/ratioBaryonToPhotonMomentumDensityEquality)                                                 &
         &                                     *log(                                                                                                   &
         &                                          +(                                                                                                 &
         &                                                  +sqrt(+1.0d0                                     +ratioBaryonToPhotonMomentumDensityDrag)  &
         &                                                  +sqrt(+ratioBaryonToPhotonMomentumDensityEquality+ratioBaryonToPhotonMomentumDensityDrag)  &
         &                                           )                                                                                                 &
         &                                          /(+1.0d0+sqrt(+ratioBaryonToPhotonMomentumDensityEquality                                       )) &
         &                                         )
    ! Find the wavenumber for Silk damping.
    self%wavenumberSilk                       =+1.6d0                                               &
         &                                     *                      self%densityBaryons **0.52d0  &
         &                                     *                      self%densityMatter  **0.73d0  &
         &                                     *(+1.0d0+1.0d0/(10.4d0*self%densityMatter )**0.95d0)
    ! Evaluate coefficients of the CDM transfer function (eqn. 11 & 12).
    a1                                        =+(46.9d0*self%densityMatter)**0.670d0*(1.0d0+1.0d0/(32.1d0*self%densityMatter)**0.532d0)
    a2                                        =+(12.0d0*self%densityMatter)**0.424d0*(1.0d0+1.0d0/(45.0d0*self%densityMatter)**0.582d0)
    self%alphaDarkMatter                      =+a1**(-self%fractionBaryons   ) &
         &                                     *a2**(-self%fractionBaryons**3)    
    b1                                        =+0.944d0/(+1.0d0+1.0d0/(+458.0d0*self%densityMatter)**0.708d0)
    b2                                        =+1.0d0/(+0.395d0*self%densityMatter)**0.0266d0
    self%betaDarkMatter                       =+1.0d0/(+1.0d0+b1*((self%omegaDarkMatter/self%omegaMatter)**b2-1.0d0))
    ! Compute the factor G(y) (eq. 15).
    y                                         =+(1.0d0+redshiftEquality   ) &
         &                                     /(1.0d0+redshiftComptonDrag)
    Gy                                        =+y                                                                    &
         &                                     *(                                                                    &
         &                                       -  6.0d0         *      sqrt(1.0d0+y)                               &
         &                                       +(+2.0d0+3.0d0*y)*log((+sqrt(1.0d0+y)+1.0d0)/(sqrt(1.0d0+y)-1.0d0)) &
         &                                     )
    ! Compute the normalization of the baryonic transfer function (eq. 14).
    self%alphaBaryons                         =+2.07d0*self%wavenumberEquality*self%distanceSoundHorizon*(+1.0d0+ratioBaryonToPhotonMomentumDensityDrag)**(-3.0d0/4.0d0)*Gy
    ! Compute node shift factors.
    ! Evaluate velocity term characteristic scale (eq. 24).
    self%betaBaryons                          =+0.50d0+self%fractionBaryons+(+3.0d0-2.0d0*self%fractionBaryons)*sqrt((17.2d0*self%densityMatter)**2+1.0d0)
    ! Evaluate node shift (eq. 23).
    self%bNode                                =+8.41d0*self%densityMatter**0.435d0
    ! Initialize wavenumber for which results were computed to a non-physical value.
    self%wavenumberPrevious                   =-1.0d0
    return
  end function eisensteinHu1998ConstructorInternal

  subroutine eisensteinHu1998Destructor(self)
    !!{
    Destructor for the eisensteinHu1998 transfer function class.
    !!}
    implicit none
    type(transferFunctionEisensteinHu1998), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterParticle_" />
    !!]
    return
  end subroutine eisensteinHu1998Destructor

  subroutine eisensteinHu1998ComputeFactors(self,wavenumber)
    !!{
    Compute common factors required by {\normalfont \ttfamily eisensteinHu1998} transfer function class.
    !!}
    use :: Numerical_Constants_Math, only : e
    implicit none
    class           (transferFunctionEisensteinHu1998), intent(inout) :: self
    double precision                                  , intent(in   ) :: wavenumber
    double precision                                                  :: wavenumberScaleFree, DwavenumberScaleFree, &
         &                                                               f                  , Df                  , &
         &                                                               T0t                , DT0t                , &
         &                                                               T0t1bc             , DT0t1bc             , &
         &                                                               T0t11              , DT0t11              , &
         &                                                               st                 , Dst                 , &
         &                                                               termBaryons1       , DtermBaryons1       , &
         &                                                               termBaryons2       , DtermBaryons2       , &
         &                                                               termBaryons3       , DtermBaryons3       , &
         &                                                               C                  , DC                  , &
         &                                                               C11                , DC11                , &
         &                                                               C1bc               , DC1bc
    
    ! If called again with the same wavenumber, return without recomputing result.
    if (wavenumber == self%wavenumberPrevious) return
    wavenumberScaleFree=+     wavenumber         &
         &              /self%wavenumberEquality &
         &              /13.41d0
    ! Evaluate the CDM transfer function.   
    ! Compute interpolating factor (eq. 18).
    f                               =+1.0d0/(1.0d0+(wavenumber*self%distanceSoundHorizon/5.4d0)**4)
    ! Evaluate fitting coefficient (eq. 20).
    C                               =+14.2d0/self%alphaDarkMatter+386.0d0/(+1.0d0+69.9d0*wavenumberScaleFree**1.08d0)
    ! Evaluate the transfer function (eq. 19).
    T0t                             =+log(e+1.8d0*self%betaDarkMatter*wavenumberScaleFree)/(log(e+1.8d0*self%betaDarkMatter*wavenumberScaleFree)+C   *wavenumberScaleFree**2)
    ! Interpolate to get the transfer function (eq. 17).
    C1bc                            =+14.2d0                     +386.0d0/(+1.0d0+69.9d0*wavenumberScaleFree**1.08d0)
    T0t1bc                          =+log(e+1.8d0*self%betaDarkMatter*wavenumberScaleFree)/(log(e+1.8d0*self%betaDarkMatter*wavenumberScaleFree)+C1bc*wavenumberScaleFree**2)
    self%transferFunctionDarkMatter =+        f *T0t1bc &
         &                           +(+1.0d0-f)*T0t
    ! Compute the baryonic transfer function.
    ! Compute node shift (eq. 22).
    st                              =self%distanceSoundHorizon/(+1.0d0+(self%bNode/wavenumber/self%distanceSoundHorizon)**3)**(1.0d0/3.0d0)
    ! Baryonic transfer function (eq. 21).
    C11                             =+14.2d0+386.0d0/(+1.0d0+69.9d0*wavenumberScaleFree**1.08d0)
    T0t11                           =+log(e+1.8d0*wavenumberScaleFree)/(log(e+1.8d0*wavenumberScaleFree)+C11*wavenumberScaleFree**2)
    termBaryons1                    =+T0t11            /(1.0d0+(                 wavenumber*self%distanceSoundHorizon/5.2d0)**2)
    termBaryons2                    =+self%alphaBaryons/(1.0d0+(self%betaBaryons/wavenumber/self%distanceSoundHorizon      )**3)*exp(-(wavenumber/self%wavenumberSilk)**1.4d0)
    termBaryons3                    =+sin(wavenumber*st) &
         &                           /   (wavenumber*st)
    self%transferFunctionBaryons    =+(termBaryons1+termBaryons2)*termBaryons3
    ! Evaluate derivatives of functons.
    DwavenumberScaleFree            =+1.0d0                   &
         &                           /self%wavenumberEquality &
         &                           /13.41d0
    Df                              =- 4.0d0*  wavenumber**3*(self%distanceSoundHorizon/5.4d0)**4     &
         &                           /(1.0d0+(+wavenumber   * self%distanceSoundHorizon/5.4d0)**4)**2
    DC                              =-386.0d0*69.9d0*1.08d0*wavenumberScaleFree**0.08d0/(+1.0d0+69.9d0*wavenumberScaleFree**1.08d0)**2
    DC1bc                           =-386.0d0*69.9d0*1.08d0*wavenumberScaleFree**0.08d0/(+1.0d0+69.9d0*wavenumberScaleFree**1.08d0)**2
    DC11                            =-386.0d0*69.9d0*1.08d0*wavenumberScaleFree**0.08d0/(+1.0d0+69.9d0*wavenumberScaleFree**1.08d0)**2
    DT0t11                          =+DT0(C11 ,DC11 ,1.0d0              )
    DT0t                            =+DT0(C   ,DC   ,self%betaDarkMatter)
    DT0t1bc                         =+DT0(C1bc,DC1bc,self%betaDarkMatter)
    Dst                             =+self%bNode**3/(+1.0d0+(self%bNode/wavenumber/self%distanceSoundHorizon)**3)**(4.0d0/3.0d0)/wavenumber**4/self%distanceSoundHorizon**2
    DtermBaryons1                   =-2.0d0*wavenumber*self%distanceSoundHorizon**2*T0t11/5.2d0**2/(1.0d0+(wavenumber*self%distanceSoundHorizon/5.2d0)**2)**2+DT0t11/(1.0d0+(wavenumber*self%distanceSoundHorizon/5.2d0)**2)
    DtermBaryons2                   =termBaryons2*(3.0d0/(1.0d0+(wavenumber*self%distanceSoundHorizon/self%betaBaryons)**3)-1.4d0*(wavenumber/self%wavenumberSilk)**1.4d0)/wavenumber
    DtermBaryons3                   =(cos(wavenumber*st)/(wavenumber*st)-sin(wavenumber*st)/(wavenumber*st)**2)*(st+wavenumber*Dst)
    self%DtransferFunctionDarkMatter=+Df*T0t1bc+       f *DT0t1bc &
         &                           -Df*T0t   +(1.0d0-f)*DT0t
    self%DtransferFunctionBaryons   =+(DtermBaryons1+DtermBaryons2)* termBaryons3 &
         &                           +( termBaryons1+ termBaryons2)*DtermBaryons3
    ! Record the wavenumber for which factors were computed.
    self%wavenumberPrevious        =+wavenumber
    return
    
  contains
    
    double precision function DT0(C,DC,beta)
      !!{
      Evaluate derivatices of the $T_0$ factors.
      !!}
      implicit none
      double precision, intent(in   ) :: C, DC, beta

      DT0    =+(                                                                       &
           &    -(                                                                     &
           &      +(                                                                   &
           &        +(                                                                 &
           &          +2.0d0*C*wavenumberScaleFree                                     &
           &          +(  1.8d0                    *beta)                              &
           &          /(e+1.8d0*wavenumberScaleFree*beta)                              &
           &         )                                                                 &
           &        *                        log(e+1.8d0*wavenumberScaleFree*beta)     &
           &       )                                                                   &
           &      /(C*wavenumberScaleFree**2+log(e+1.8d0*wavenumberScaleFree*beta))**2 &
           &     )                                                                     &
           &    +     (  1.8d0                    *beta)                               &
           &    /(                                                                     &
           &          (e+1.8d0*wavenumberScaleFree*beta)                               &
           &      *(C*wavenumberScaleFree**2+log(e+1.8d0*wavenumberScaleFree*beta))    &
           &     )                                                                     &
           &   )                                                                       &
           &  *DwavenumberScaleFree                                                    &
           &  -wavenumberScaleFree**2                                                  &
           &  *                          log(e+1.8d0*wavenumberScaleFree*beta)         &
           &  /(C*wavenumberScaleFree**2+log(e+1.8d0*wavenumberScaleFree*beta))**2     &
           &  *DC                                                                      &
           &  *DwavenumberScaleFree
      return
    end function DT0

  end subroutine eisensteinHu1998ComputeFactors

  double precision function eisensteinHu1998Value(self,wavenumber)
    !!{
    Return the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionEisensteinHu1998), intent(inout) :: self
    double precision                                  , intent(in   ) :: wavenumber

    ! Compute common factors.
    call self%computeFactors(wavenumber)
    ! Compute the total transfer function.
    eisensteinHu1998Value=+self%fractionBaryons                 *self%transferFunctionBaryons    &
         &                +self%omegaDarkMatter/self%omegaMatter*self%transferFunctionDarkMatter
    return
  end function eisensteinHu1998Value

  double precision function eisensteinHu1998LogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionEisensteinHu1998), intent(inout) :: self
    double precision                                  , intent(in   ) :: wavenumber
    double precision                                                  :: transferFunction

    ! Compute common factors.
    call self%computeFactors(wavenumber)
    ! Get the transfer function itself.
    transferFunction=self%value(wavenumber)
    ! Compute the logarithmic derivative.
    eisensteinHu1998LogarithmicDerivative=+(                                                                        &
         &                                  +self%fractionBaryons                 *self%DtransferFunctionBaryons    &
         &                                  +self%omegaDarkMatter/self%omegaMatter*self%DtransferFunctionDarkMatter &
         &                                 )                                                                        &
         &                                *wavenumber                                                               &
         &                                /transferFunction
    return
  end function eisensteinHu1998LogarithmicDerivative

  double precision function eisensteinHu1998HalfModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    to a \gls{cdm} transfer function. Not supported in this implementation.
    !!}
    use :: Error, only : Error_Report, errorStatusFail
    implicit none
    class  (transferFunctionEisensteinHu1998), intent(inout), target   :: self
    integer                                  , intent(  out), optional :: status
    !$GLC attributes unused :: self

    eisensteinHu1998HalfModeMass=0.0d0
    if (present(status)) then
       status=errorStatusFail
    else
       call Error_Report('not supported by this implementation'//{introspection:location})
    end if
    return
  end function eisensteinHu1998HalfModeMass

  double precision function eisensteinHu1998QuarterModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of four relative
    to a \gls{cdm} transfer function. Not supported in this implementation.
    !!}
    use :: Error, only : Error_Report, errorStatusFail
    implicit none
    class  (transferFunctionEisensteinHu1998), intent(inout), target   :: self
    integer                                  , intent(  out), optional :: status
    !$GLC attributes unused :: self

    eisensteinHu1998QuarterModeMass=0.0d0
    if (present(status)) then
       status=errorStatusFail
    else
       call Error_Report('not supported by this implementation'//{introspection:location})
    end if
    return
  end function eisensteinHu1998QuarterModeMass

  double precision function eisensteinHu1998EpochTime(self)
    !!{
    Return the cosmic time at the epoch at which this transfer function is defined.
    !!}
    implicit none
    class(transferFunctionEisensteinHu1998), intent(inout) :: self

    eisensteinHu1998EpochTime=self%time
    return
  end function eisensteinHu1998EpochTime
