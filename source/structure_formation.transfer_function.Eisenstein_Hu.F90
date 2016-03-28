!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.
  
  !% Contains a module which implements a transfer function class using the fitting function of
  !% \cite{eisenstein_power_1999}.

  use Cosmology_Parameters

  !# <transferFunction name="transferFunctionEisensteinHu1999">
  !#  <description>Provides the \cite{eisenstein_power_1999} fitting function to the transfer function. The effective number of neutrino species and the summed mass (in electron volts) of all neutrino species are specified via the {\normalfont \ttfamily neutrinoNumberEffective} and {\normalfont \ttfamily neutrinoMassSummed} parameters respectively.</description>
  !# </transferFunction>
  type, extends(transferFunctionClass) :: transferFunctionEisensteinHu1999
     !% The ``{\normalfont \ttfamily eisensteinHu1999}'' transfer function class.
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_
     double precision                                    :: temperatureCMB27    , distanceSoundWave      , &
          &                                                 neutrinoMassFraction, neutrinoNumberEffective, &
          &                                                 neutrinoFactor      , betaDarkMatter         , &
          &                                                 neutrinoMassSummed
   contains
     !@ <objectMethods>
     !@   <object>transferFunctionEisensteinHu1999</object>
     !@   <objectMethod>
     !@     <method>computeFactors</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ wavenumber\argin, \doublezero\ wavenumberEffective\argout, \doublezero\ wavenumberNeutrino\argout, \doublezero\ L\argout, \doublezero\ C\argout</arguments>
     !@     <description>Compute common factors needed by \cite{eisenstein_power_1999} transfer function calculations.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                          eisensteinHu1999Destructor
     procedure :: value                 => eisensteinHu1999Value
     procedure :: logarithmicDerivative => eisensteinHu1999LogarithmicDerivative
     procedure :: computeFactors        => eisensteinHu1999ComputeFactors
     procedure :: halfModeMass          => eisensteinHu1999HalfModeMass
     procedure :: descriptor            => eisensteinHu1999Descriptor
  end type transferFunctionEisensteinHu1999

  interface transferFunctionEisensteinHu1999
     !% Constructors for the ``{\normalfont \ttfamily eisensteinHu1999}'' transfer function class.
     module procedure eisensteinHu1999ConstructorParameters
     module procedure eisensteinHu1999ConstructorInternal
  end interface transferFunctionEisensteinHu1999

contains

  function eisensteinHu1999ConstructorParameters(parameters)
    !% Constructor for the ``{\normalfont \ttfamily eisensteinHu1999}'' transfer function class
    !% which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type            (transferFunctionEisensteinHu1999)                :: eisensteinHu1999ConstructorParameters
    type            (inputParameters                 ), intent(in   ) :: parameters
    class           (cosmologyParametersClass        ), pointer       :: cosmologyParameters_
    double precision                                                  :: neutrinoNumberEffective             , neutrinoMassSummed
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>neutrinoNumberEffective</name>
    !#   <source>parameters</source>
    !#   <defaultValue>3.046d0</defaultValue>
    !#   <defaultSource>\citep{mangano_relic_2005}</defaultSource>
    !#   <description>The effective number of neutrino species.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>neutrinoMassSummed</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The summed mass (in electron volts) of all neutrino species.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    ! Call the internal constructor.
    eisensteinHu1999ConstructorParameters=eisensteinHu1999ConstructorInternal(neutrinoNumberEffective,neutrinoMassSummed,cosmologyParameters_)
    return
  end function eisensteinHu1999ConstructorParameters

  function eisensteinHu1999ConstructorInternal(neutrinoNumberEffective,neutrinoMassSummed,cosmologyParameters_)
    !% Internal constructor for the ``{\normalfont \ttfamily eisensteinHu1999}'' transfer function class.
    implicit none
    type            (transferFunctionEisensteinHu1999)                                  :: eisensteinHu1999ConstructorInternal
    double precision                                  , intent(in   )                   :: neutrinoNumberEffective           , neutrinoMassSummed
    class           (cosmologyParametersClass        ), intent(in   ), target, optional :: cosmologyParameters_
    double precision                                                                    :: redshiftEquality                   , redshiftComptonDrag   , &
         &                                                                                 b1                                 , b2                    , &
         &                                                                                 massFractionBaryonic               , massFractionDarkMatter, &
         &                                                                                 expansionFactorRatio               , massFractionMatter    , &
         &                                                                                 massFractionBaryonsNeutrinos       , suppressionDarkMatter , &
         &                                                                                 suppressionMatter

    ! Determine the cosmological parameters to use.
    if (present(cosmologyParameters_)) then
       eisensteinHu1999ConstructorInternal%cosmologyParameters_ => cosmologyParameters_
    else
       eisensteinHu1999ConstructorInternal%cosmologyParameters_ => cosmologyParameters()
    end if
    ! Present day CMB temperature [in units of 2.7K].
    eisensteinHu1999ConstructorInternal%temperatureCMB27       =+eisensteinHu1999ConstructorInternal%cosmologyParameters_%temperatureCMB(                  )    &
         &                                                       /2.7d0
    ! Redshift of matter-radiation equality.
    redshiftEquality                                           =+2.50d4                                                                                         &
         &                                                      *eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaMatter   (                  )    &
         &                                                      *eisensteinHu1999ConstructorInternal%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                                                      /eisensteinHu1999ConstructorInternal%temperatureCMB27                                       **4
    ! Compute redshift at which baryons are released from Compton drag of photons (eqn. 2)
    b1                                                         =+0.313d0                                                                                            &
         &                                                      /(                                                                                                  &
         &                                                        +  eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaMatter   (                  )    &
         &                                                        *  eisensteinHu1999ConstructorInternal%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                                                       )**0.419d0                                                                                         &
         &                                                      *(                                                                                                  &
         &                                                        +1.000d0                                                                                          &
         &                                                        +0.607d0                                                                                          &
         &                                                        *(                                                                                                &
         &                                                          +eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaMatter   (                  )    &
         &                                                          *eisensteinHu1999ConstructorInternal%cosmologyParameters_%hubbleConstant(hubbleUnitsLittleH)**2 &
         &                                                         )**0.674d0                                                                                       &
         &                                                       )
    b2                                                         =+0.238d0                                                                                            &
         &                                                      *(                                                                                                  &
         &                                                        +  eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaMatter   (                  )    &
         &                                                        *  eisensteinHu1999ConstructorInternal%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                                                       )**0.223d0
    redshiftComptonDrag                                        =+1291.0d0                                                                                           &
         &                                                      *(                                                                                                  &
         &                                                        +  eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaMatter   (                  )    &
         &                                                        *  eisensteinHu1999ConstructorInternal%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                                                       )**0.251d0                                                                                         &
         &                                                      *(                                                                                                  &
         &                                                        +1.0d0                                                                                            &
         &                                                        +b1                                                                                               &
         &                                                        *(                                                                                                &
         &                                                          +eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaBaryon   (                  )    &
         &                                                          *eisensteinHu1999ConstructorInternal%cosmologyParameters_%hubbleConstant(hubbleUnitsLittleH)**2 &
         &                                                         )**b2                                                                                            &
         &                                                       )                                                                                                  &
         &                                                      /(                                                                                                  &
         &                                                        +1.000d0                                                                                          &
         &                                                        +0.659d0                                                                                          &
         &                                                        *(                                                                                                &
         &                                                          +eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaMatter   (                  )    &
         &                                                          *eisensteinHu1999ConstructorInternal%cosmologyParameters_%hubbleConstant(hubbleUnitsLittleH)**2 &
         &                                                         )**0.828d0                                                                                       &
         &                                                       )
    ! Relative expansion factor between previous two computed redshifts.
    expansionFactorRatio                                       =+(1.0d0+redshiftEquality   ) &
         &                                                      /(1.0d0+redshiftComptonDrag)
    ! Compute the comoving distance that a sound wave can propagate prior to redshiftComptonDrag (i.e. sound horizon; eq. 4)
    eisensteinHu1999ConstructorInternal%distanceSoundWave      =+44.5d0                                                                                                 &
         &                                                      *log (                                                                                                  &
         &                                                            +9.83d0                                                                                           &
         &                                                            /  eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaMatter   (                  )    &
         &                                                            /  eisensteinHu1999ConstructorInternal%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                                                           )                                                                                                  &
         &                                                      /sqrt(                                                                                                  &
         &                                                            + 1.0d0                                                                                           &
         &                                                            +10.0d0                                                                                           &
         &                                                            *(                                                                                                &
         &                                                              +eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaBaryon   (                  )    &
         &                                                              *eisensteinHu1999ConstructorInternal%cosmologyParameters_%hubbleConstant(hubbleUnitsLittleH)**2 &
         &                                                             )**0.75d0                                                                                        &
         &                                                           )
    ! Specify properties of neutrinos. Mass fraction formula is from Komatsu et al. (2007; http://adsabs.harvard.edu/abs/2010arXiv1001.4538K).
    eisensteinHu1999ConstructorInternal%neutrinoMassSummed     =+neutrinoMassSummed
    eisensteinHu1999ConstructorInternal%neutrinoMassFraction   =+neutrinoMassSummed                                                                           &
         &                                                      /94.0d0                                                                                         &
         &                                                      /eisensteinHu1999ConstructorInternal%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                                                      /eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaMatter   (                  )
    eisensteinHu1999ConstructorInternal%neutrinoNumberEffective=+neutrinoNumberEffective
    ! Compute baryonic and cold dark matter fractions.
    massFractionBaryonic                                       =+  eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaBaryon() &
         &                                                      /  eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaMatter()
    massFractionDarkMatter                                     =+(                                                                        &
         &                                                        +eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaMatter() &
         &                                                        -eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaBaryon() &
         &                                                       )                                                                        &
         &                                                      /  eisensteinHu1999ConstructorInternal%cosmologyParameters_%OmegaMatter()
    ! Total matter fraction.
    massFractionMatter                                         =+massFractionBaryonic   &
         &                                                      +massFractionDarkMatter
    ! Baryonic + neutrino fraction.
    massFractionBaryonsNeutrinos                               =+eisensteinHu1999ConstructorInternal%neutrinoMassFraction &
         &                                                      +massFractionBaryonic
    ! Compute small scale suppression factor (eqn. 15).
    suppressionDarkMatter                                      =+0.25d0*(5.0d0-sqrt(1.0d0+24.0d0*massFractionDarkMatter))
    suppressionMatter                                          =+0.25d0*(5.0d0-sqrt(1.0d0+24.0d0*massFractionMatter    ))
    eisensteinHu1999ConstructorInternal%neutrinoFactor         =+(                                                                      &
         &                                                        +massFractionDarkMatter                                               &
         &                                                        /massFractionMatter                                                   &
         &                                                       )                                                                      &
         &                                                      *(                                                                      &
         &                                                        +(5.0d0-2.0d0*(suppressionDarkMatter+suppressionMatter))              &
         &                                                        /(5.0d0-4.0d0*                       suppressionMatter )              &
         &                                                       )                                                                      &
         &                                                      *(                                                                      &
         &                                                        +(                                                                    &
         &                                                          +1.000d0                                                            &
         &                                                          -0.533d0                                                            &
         &                                                          *massFractionBaryonsNeutrinos                                       &
         &                                                          +0.126d0                                                            &
         &                                                          *massFractionBaryonsNeutrinos**3                                    &
         &                                                         )                                                                    &
         &                                                        *(                                                                    &
         &                                                          +1.0d0                                                              &
         &                                                          +expansionFactorRatio                                               &
         &                                                         )**(                                                                 &
         &                                                             +suppressionMatter                                               &
         &                                                             -suppressionDarkMatter                                           &
         &                                                            )                                                                 &
         &                                                        /(                                                                    &
         &                                                          +1.000d0                                                            &
         &                                                          -0.193d0                                                            &
         &                                                          *sqrt(                                                              &
         &                                                                +eisensteinHu1999ConstructorInternal%neutrinoMassFraction     &
         &                                                                *eisensteinHu1999ConstructorInternal%neutrinoNumberEffective  &
         &                                                               )                                                              &
         &                                                          +0.169d0                                                            &
         &                                                          *eisensteinHu1999ConstructorInternal%neutrinoMassFraction           &
         &                                                          *eisensteinHu1999ConstructorInternal%neutrinoNumberEffective**0.2d0 &
         &                                                         )                                                                    &
         &                                                       )                                                                      &
         &                                                      *(                                                                      &
         &                                                        +1.0d0                                                                &
         &                                                        +0.5d0                                                                &
         &                                                        *(                                                                    &
         &                                                          +suppressionDarkMatter                                              &
         &                                                          -suppressionMatter                                                  &
         &                                                         )                                                                    &
         &                                                        *(                                                                    &
         &                                                          +1.0d0                                                              &
         &                                                          +1.0d0                                                              &
         &                                                          /(3.0d0-4.0d0*suppressionDarkMatter)                                &
         &                                                          /(7.0d0-4.0d0*suppressionMatter    )                                &
         &                                                         )                                                                    &
         &                                                        /(                                                                    &
         &                                                          +1.0d0                                                              &
         &                                                          +expansionFactorRatio                                               &
         &                                                         )                                                                    &
         &                                                       )
    eisensteinHu1999ConstructorInternal%betaDarkMatter      =+1.0d0                                                                     & ! Eqn. 21.
         &                                                   /(                                                                         &
         &                                                     +1.0d0                                                                   &
         &                                                     -0.949d0                                                                 &
         &                                                     *massFractionBaryonsNeutrinos                                            &
         &                                                    )
    return
  end function eisensteinHu1999ConstructorInternal

  subroutine eisensteinHu1999Destructor(self)
    !% Destructor for the eisensteinHu1999 transfer function class.
    implicit none
    type(transferFunctionEisensteinHu1999), intent(inout) :: self

    if     (                                                       &
         &   associated(self%cosmologyParameters_                ) &
         &  .and.                                                  &
         &              self%cosmologyParameters_%isFinalizable()  &
         & ) deallocate(self%cosmologyParameters_                )
    return
  end subroutine eisensteinHu1999Destructor

  subroutine eisensteinHu1999ComputeFactors(self,wavenumber,wavenumberEffective,wavenumberNeutrino,L,C)
    !% Compute common factors required by ``{\normalfont \ttfamily eisensteinHu1999}'' transfer function class.
    implicit none
    class           (transferFunctionEisensteinHu1999), intent(in   ) :: self
    double precision                                  , intent(in   ) :: wavenumber
    double precision                                  , intent(  out) :: wavenumberEffective    , wavenumberNeutrino, &
         &                                                               L                      , C
    double precision                                                  :: wavenumberScaleFree
    double precision                                                  :: shapeParameterEffective
    
    ! Compute rescaled shape parameter (eqn. 16)
    shapeParameterEffective=+self%cosmologyParameters_%OmegaMatter   (                  )    &
         &                  *self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                  *(                                                               &
         &                    +       sqrt(self%neutrinoFactor)                              &
         &                    +(1.0d0-sqrt(self%neutrinoFactor))                             &
         &                    /(                                                             &
         &                      +1.0d0                                                       &
         &                      +(                                                           &
         &                        +0.43d0                                                    &
         &                        *wavenumber                                                &
         &                        *self%distanceSoundWave                                    &
         &                       )**4                                                        &
         &                     )                                                             &
         &                   )
    wavenumberEffective    =+wavenumber                                                      &
         &                  *self%temperatureCMB27**2                                        &
         &                  /shapeParameterEffective
    L                      =+log(                                                            & ! Eqn. 19.
         &                       +exp(1.0d0)                                                 &
         &                       +1.84d0                                                     &
         &                       *self%betaDarkMatter                                        &
         &                       *sqrt(self%neutrinoFactor)                                  &
         &                       *wavenumberEffective                                        &
         &                      )
    C                      =+ 14.4d0                                                         & ! Eqn. 20.
         &                  +325.0d0                                                         &
         &                  /(                                                               &
         &                    + 1.0d0                                                        &
         &                    +60.5d0                                                        &
         &                    *wavenumberEffective**1.11d0                                   &
         &                   )
    ! Compute wavenumbers needed for horizon scale modes.
    if     (                                    &
         &   self%neutrinoMassFraction >  0.0d0 &
         &  .and.                               &
         &   self%neutrinoMassFraction <= 0.3d0 &
         & ) then
       ! Compute effective q.
       wavenumberScaleFree =+wavenumber                                                      &
            &               *self%temperatureCMB27                                       **2 &
            &               /self%cosmologyParameters_%OmegaMatter   (                  )    &
            &               /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2
       wavenumberNeutrino =+3.92d0                             &
            &              *wavenumberScaleFree                &
            &              *sqrt(self%neutrinoNumberEffective) &
            &              /     self%neutrinoMassFraction
    else
       wavenumberNeutrino =+0.0d0
    end if
    return
  end subroutine eisensteinHu1999ComputeFactors
  
  double precision function eisensteinHu1999Value(self,wavenumber)
    !% Return the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionEisensteinHu1999), intent(inout) :: self
    double precision                                  , intent(in   ) :: wavenumber
    double precision                                                  :: wavenumberEffective, L                      , &
         &                                                               C                  , wavenumberNeutrino     , &
         &                                                               suppressionNeutrino

    ! Compute common factors.
    call self%computeFactors(wavenumber,wavenumberEffective,wavenumberNeutrino,L,C)
    ! Evaluate the transfer function.
    eisensteinHu1999Value  =+  L                                                             & ! Zero baryon form of the transfer function (eqn. 18).
         &                  /(                                                               &
         &                    +L                                                             &
         &                    +C                                                             &
         &                    *wavenumberEffective**2                                        &
         &                   )                              
    ! Apply correction for scales close to horizon.
    if     (                                    &
         &   self%neutrinoMassFraction >  0.0d0 &
         &  .and.                               &
         &   self%neutrinoMassFraction <= 0.3d0 &
         & ) then
       suppressionNeutrino=+1.0d0                                                                    &
            &              +(                                                                        &
            &                +1.2d0                                                                  &
            &                *self%neutrinoMassFraction   ** 0.64d0                                  &
            &                *self%neutrinoNumberEffective**(0.30d0+0.6d0*self%neutrinoMassFraction) &
            &               )                                                                        &
            &              /(                                                                        &
            &                +wavenumberNeutrino**(-1.6d0)                                           &
            &                +wavenumberNeutrino**(+0.8d0)                                           &
            &               )
       eisensteinHu1999Value=eisensteinHu1999Value*suppressionNeutrino
    end if
    return
  end function eisensteinHu1999Value

  double precision function eisensteinHu1999LogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionEisensteinHu1999), intent(inout) :: self
    double precision                                  , intent(in   ) :: wavenumber
    double precision                                                  :: wavenumberEffective   , L                     , &
         &                                                               C                     , wavenumberNeutrino    , &
         &                                                               transferFunction      , dLdwavenumberEffective, &
         &                                                               dCdwavenumberEffective

    ! Compute common factors.
    call self%computeFactors(wavenumber,wavenumberEffective,wavenumberNeutrino,L,C)
    ! Get the transfer function itself.
    transferFunction=self%value(wavenumber)
    ! Compute derivatives of transfer function terms.
    dCdwavenumberEffective                                 =-325.00d0&
         &                                * 60.50d0&
         &                                *  1.11d0&
         &                                /(&
         &                                  + 1.0d0&
         &                                  +60.5d0&
         &                                  *wavenumberEffective**1.11d0&
         &                                 )**2&
         &                                *  wavenumberEffective**0.11d0    
    dLdwavenumberEffective                                 =+1.0d0&
         &                                /(&
         &                                  +exp(1.0d0)&
         &                                  /1.84d0                                                     &
         &                                  /self%betaDarkMatter                                        &
         &                                  /sqrt(self%neutrinoFactor)&
         &                                  +wavenumberEffective&
         &                                 )
    ! Compute logarithmic derivative of transfer function.
    eisensteinHu1999LogarithmicDerivative=+(                                                                                                        &
         &                                  +dLdwavenumberEffective                                                                                 &
         &                                  /(                            + L                    + C                    *wavenumberEffective**2)    &
         &                                  - L                                                                                                     &
         &                                  *(+2.0d0*C*wavenumberEffective+dLdwavenumberEffective+dCdwavenumberEffective*wavenumberEffective**2)    &
         &                                  /(                            + L                    + C                    *wavenumberEffective**2)**2 &
         &                                 )                                                                                                        &
         &                                 *wavenumberEffective                                                                                     &
         &                                 /transferFunction
    ! Apply correction for scales close to horizon.
    if     (                                    &
         &   self%neutrinoMassFraction >  0.0d0 &
         &  .and.                               &
         &   self%neutrinoMassFraction <= 0.3d0 &
         & ) eisensteinHu1999LogarithmicDerivative=eisensteinHu1999LogarithmicDerivative+0.8d0*(3.0d0/(1.0d0+wavenumberNeutrino**2.4d0)-1.0d0)
    return
  end function eisensteinHu1999LogarithmicDerivative

  double precision function eisensteinHu1999HalfModeMass(self)
    !% Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    !% to a \gls{cdm} transfer function. Not supported in this implementation.
    use Galacticus_Error
    implicit none
    class(transferFunctionEisensteinHu1999), intent(inout) :: self

    call Galacticus_Error_Report('eisensteinHu1999HalfModeMass','not supported by this implementation')
    return
  end function eisensteinHu1999HalfModeMass

  subroutine eisensteinHu1999Descriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class    (transferFunctionEisensteinHu1999), intent(inout) :: self
    type     (inputParameters                 ), intent(inout) :: descriptor
    type     (node                            ), pointer       :: parameterNode
    type     (inputParameters                 )                :: subParameters
    character(len=10                          )                :: parameterLabel

    call descriptor%addParameter("transferFunctionMethod","eisensteinHu1999")
    parameterNode => descriptor%node("transferFunctionMethod")
    subParameters=inputParameters(parameterNode)
    write (parameterLabel,'(f10.6)') self%neutrinoMassSummed
    call subParameters%addParameter("neutrinoMassSummed"     ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(f10.6)') self%neutrinoNumberEffective
    call subParameters%addParameter("neutrinoNumberEffective",trim(adjustl(parameterLabel)))
    call self%cosmologyParameters_%descriptor(subParameters)
    return
  end subroutine eisensteinHu1999Descriptor
