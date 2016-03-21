!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of dark matter halo profile concentrations using the
  !% \cite{diemer_universal_2014} algorithm.

  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationDiemerKravtsov2014">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{diemer_universal_2014}.</description>
  !# </darkMatterProfileConcentration>
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationDiemerKravtsov2014
     !% A dark matter halo profile concentration class implementing the algorithm of
     !% \cite{diemer_universal_2014}.
     private
     double precision           :: kappa, phi0, phi1, eta0, eta1, alpha, beta, scatter
     type            (fgsl_rng) :: clonedPseudoSequenceObject, pseudoSequenceObject
     logical                    :: resetSequence             , resetSequenceSnapshot
   contains
     final     ::                                diemerKravtsov2014Destructor
     procedure :: concentration               => diemerKravtsov2014Concentration
     procedure :: concentrationMean           => diemerKravtsov2014ConcentrationMean
     procedure :: densityContrastDefinition   => diemerKravtsov2014DensityContrastDefinition
     procedure :: darkMatterProfileDefinition => diemerKravtsov2014DarkMatterProfileDefinition
     procedure :: stateStore                  => diemerKravtsov2014StateStore
     procedure :: stateRestore                => diemerKravtsov2014StateRestore
     procedure :: stateSnapshot               => diemerKravtsov2014StateSnapshot
  end type darkMatterProfileConcentrationDiemerKravtsov2014

  interface darkMatterProfileConcentrationDiemerKravtsov2014
     !% Constructors for the {\normalfont \ttfamily diemerKravtsov2014} dark matter halo profile concentration class.
     module procedure diemerKravtsov2014ConstructorParameters
     module procedure diemerKravtsov2014ConstructorInternal
  end interface darkMatterProfileConcentrationDiemerKravtsov2014

contains

  function diemerKravtsov2014ConstructorParameters(parameters)
    !% Default constructor for the {\normalfont \ttfamily diemerKravtsov2014} dark matter halo
    !% profile concentration class.
    implicit none
    type(darkMatterProfileConcentrationDiemerKravtsov2014)                :: diemerKravtsov2014ConstructorParameters
    type(inputParameters                                 ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>kappa</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%kappa</variable>
    !#   <defaultValue>0.69d0</defaultValue>
    !#   <description>The parameter $\kappa$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>phi0</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%phi0</variable>
    !#   <defaultValue>6.58d0</defaultValue>
    !#   <description>The parameter $\phi_0$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>phi1</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%phi1</variable>
    !#   <defaultValue>1.37d0</defaultValue>
    !#   <description>The parameter $\phi_1$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>eta0</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%eta0</variable>
    !#   <defaultValue>6.82d0</defaultValue>
    !#   <description>The parameter $\eta_0$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>eta1</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%eta1</variable>
    !#   <defaultValue>1.42d0</defaultValue>
    !#   <description>The parameter $\eta_1$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>alpha</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%alpha</variable>
    !#   <defaultValue>1.12d0</defaultValue>
    !#   <attachedTo>module</attachedTo>
    !#   <description>The parameter $\alpha$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>beta</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%beta</variable>
    !#   <defaultValue>1.69d0</defaultValue>
    !#   <description>The parameter $\beta$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>scatter</name>
    !#   <source>parameters</source>
    !#   <variable>diemerKravtsov2014ConstructorParameters%scatter</variable>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The scatter (in dex) to assume in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    diemerKravtsov2014ConstructorParameters%resetSequence=.true.
    return
  end function diemerKravtsov2014ConstructorParameters

  function diemerKravtsov2014ConstructorInternal(kappa,phi0,phi1,eta0,eta1,alpha,beta,scatter)
    !% Constructor for the {\normalfont \ttfamily diemerKravtsov2014} dark matter halo profile
    !% concentration class.
    use Galacticus_Error
    implicit none
    type            (darkMatterProfileConcentrationDiemerKravtsov2014)                :: diemerKravtsov2014ConstructorInternal
    double precision                                                  , intent(in   ) :: kappa, phi0, phi1, eta0, eta1, alpha, beta, scatter
    
    diemerKravtsov2014ConstructorInternal%kappa        =kappa
    diemerKravtsov2014ConstructorInternal%phi0         =phi0
    diemerKravtsov2014ConstructorInternal%phi1         =phi1
    diemerKravtsov2014ConstructorInternal%eta0         =eta0
    diemerKravtsov2014ConstructorInternal%eta1         =eta1
    diemerKravtsov2014ConstructorInternal%alpha        =alpha
    diemerKravtsov2014ConstructorInternal%beta         =beta
    diemerKravtsov2014ConstructorInternal%scatter      =scatter
    diemerKravtsov2014ConstructorInternal%resetSequence=.true.
    return
  end function diemerKravtsov2014ConstructorInternal
  
  subroutine diemerKravtsov2014Destructor(self)
    !% Destructor for the {\normalfont \ttfamily diemerKravtsov2014} dark matter halo profile
    !% concentration class.
    implicit none
    type(darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine diemerKravtsov2014Destructor

  double precision function diemerKravtsov2014Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node}
    !% using the \cite{diemer_universal_2014} algorithm.
    use Gaussian_Random
    implicit none
    class(darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout)          :: self
    type (treeNode                                        ), intent(inout), pointer :: node

    ! Get the mean concentration.
    diemerKravtsov2014Concentration=self%concentrationMean(node)
    ! Add scatter if necessary.
    if (self%scatter > 0.0d0)                                      &
         &  diemerKravtsov2014Concentration                        &
         & =diemerKravtsov2014Concentration                        &
         & *10.0d0**Gaussian_Random_Get(                           &
         &                              self%pseudoSequenceObject, &
         &                              self%scatter             , &
         &                              self%resetSequence         &
         &                             )
    return
  end function diemerKravtsov2014Concentration

  double precision function diemerKravtsov2014ConcentrationMean(self,node)
    !% Return the mean concentration of the dark matter halo profile of {\normalfont \ttfamily node}
    !% using the \cite{diemer_universal_2014} algorithm.
    use Numerical_Constants_Math
    use Cosmological_Mass_Variance
    use Critical_Overdensities
    use Cosmology_Parameters
    use Power_Spectra
    implicit none
    class           (darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout)          :: self
    type            (treeNode                                        ), intent(inout), pointer :: node
    class           (nodeComponentBasic                              )               , pointer :: basic
    class           (cosmologyParametersClass                        )               , pointer :: cosmologyParameters_
    class           (criticalOverdensityClass                        )               , pointer :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                   )               , pointer :: cosmologicalMassVariance_
    class           (powerSpectrumClass                              )               , pointer :: powerSpectrum_
    double precision                                                                           :: radiusHaloLagrangian     , peakHeight        , &
         &                                                                                        wavenumber               , powerSpectrumSlope, &
         &                                                                                        concentrationMinimum     , peakHeightMinimum
    
    cosmologyParameters_                => cosmologyParameters      ()
    criticalOverdensity_                => criticalOverdensity      ()
    cosmologicalMassVariance_           => cosmologicalMassVariance ()
    powerSpectrum_                      => powerSpectrum            ()
    basic                               => node               %basic()
    radiusHaloLagrangian                =  +(                                                   &
         &                                   +3.0d0                                             &
         &                                   *basic%mass()                                      &
         &                                   /4.0d0                                             &
         &                                   /Pi                                                &
         &                                   /cosmologyParameters_%densityCritical()            &
         &                                   /cosmologyParameters_%OmegaMatter    ()            &
         &                                  )**(1.0d0/3.0d0)
    peakHeight                          = +criticalOverdensity_     %value       (basic%time()) &
         &                                /cosmologicalMassVariance_%rootVariance(basic%mass())
    wavenumber                          = +self%kappa                                           &
         &                                *2.0d0                                                &
         &                                *Pi                                                   &
         &                                /radiusHaloLagrangian
    powerSpectrumSlope                  = +powerSpectrum_%powerLogarithmicDerivative(wavenumber)
    concentrationMinimum                = +self%phi0                                            &
         &                                +self%phi1                                            &
         &                                *powerSpectrumSlope
    peakHeightMinimum                   = +self%eta0                                            &
         &                                +self%eta1                                            &
         &                                *powerSpectrumSlope
    diemerKravtsov2014ConcentrationMean = +0.5d0                                                &
         &                                *concentrationMinimum                                 &
         &                                *(                                                    &
         &                                  +(peakHeight/peakHeightMinimum)**(-self%alpha)      &
         &                                  +(peakHeight/peakHeightMinimum)**(+self%beta )      &
         &                                 )
    return
  end function diemerKravtsov2014ConcentrationMean

  function diemerKravtsov2014DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of concentration in the
    !% \cite{diemer_universal_2014} algorithm.
    implicit none
    class  (virialDensityContrastClass                      ), pointer                     :: diemerKravtsov2014DensityContrastDefinition
    class  (darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout)               :: self
    class  (virialDensityContrastClass                      ), allocatable  , target, save :: densityContrastDefinition
    logical                                                                         , save :: densityContrastDefinitionInitialized=.false.
    !$omp threadprivate(densityContrastDefinition,densityContrastDefinitionInitialized)
    
    if (.not.densityContrastDefinitionInitialized) then
       allocate(virialDensityContrastFixed :: densityContrastDefinition)
       select type (densityContrastDefinition)
       type is (virialDensityContrastFixed)
          densityContrastDefinition=virialDensityContrastFixed(200.0d0,virialDensityContrastFixedDensityTypeCritical)
          call densityContrastDefinition%makeIndestructible()
       end select
       densityContrastDefinitionInitialized=.true.
    end if
    diemerKravtsov2014DensityContrastDefinition => densityContrastDefinition
    return
  end function diemerKravtsov2014DensityContrastDefinition

  function diemerKravtsov2014DarkMatterProfileDefinition(self)
    !% Return a dark matter density profile object defining that used in the definition of concentration in the
    !% \cite{diemer_universal_2014} algorithm.
    use Dark_Matter_Halo_Scales
    implicit none
    class  (darkMatterProfileClass                            ), pointer                     :: diemerKravtsov2014DarkMatterProfileDefinition
    class  (darkMatterProfileConcentrationDiemerKravtsov2014  ), intent(inout)               :: self
    class  (darkMatterHaloScaleVirialDensityContrastDefinition), pointer                     :: darkMatterHaloScaleDefinition
    class  (darkMatterProfileClass                            ), allocatable  , target, save :: densityProfileDefinition
    logical                                                                           , save :: densityProfileDefinitionInitialized=.false.
    !$omp threadprivate(densityProfileDefinition,densityProfileDefinitionInitialized)
    
    if (.not.densityProfileDefinitionInitialized) then
       allocate(darkMatterProfileNFW                               :: densityProfileDefinition     )
       allocate(darkMatterHaloScaleVirialDensityContrastDefinition :: darkMatterHaloScaleDefinition)
       select type (densityProfileDefinition)
       type is (darkMatterProfileNFW)
          select type (darkMatterHaloScaleDefinition)
          type is (darkMatterHaloScaleVirialDensityContrastDefinition)
             darkMatterHaloScaleDefinition=darkMatterHaloScaleVirialDensityContrastDefinition(self%densityContrastDefinition())
             densityProfileDefinition     =darkMatterProfileNFW                              (darkMatterHaloScaleDefinition   )
             call densityProfileDefinition%makeIndestructible()
          end select
       end select
       densityProfileDefinitionInitialized=.true.
    end if
    diemerKravtsov2014DarkMatterProfileDefinition => densityProfileDefinition
    return
  end function diemerKravtsov2014DarkMatterProfileDefinition

  subroutine diemerKravtsov2014StateSnapshot(self)
    !% Write the tablulation state to file.
    implicit none
    class(darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout) :: self

    if (.not.self%resetSequence) self%clonedPseudoSequenceObject=FGSL_Rng_Clone(self%pseudoSequenceObject)
    self%resetSequenceSnapshot=self%resetSequence
    return
  end subroutine diemerKravtsov2014StateSnapshot

  subroutine diemerKravtsov2014StateStore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use Pseudo_Random
    implicit none
    class  (darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout) :: self
    integer                                                  , intent(in   ) :: stateFile
    type   (fgsl_file                                       ), intent(in   ) :: fgslStateFile

    write (stateFile) self%resetSequenceSnapshot
    if (.not.self%resetSequenceSnapshot) call Pseudo_Random_Store(self%clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine diemerKravtsov2014StateStore

  subroutine diemerKravtsov2014StateRestore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use Pseudo_Random
    implicit none
    class  (darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout) :: self
    integer                                                  , intent(in   ) :: stateFile
    type   (fgsl_file                                       ), intent(in   ) :: fgslStateFile

    read (stateFile) self%resetSequence
    if (.not.self%resetSequence) call Pseudo_Random_Retrieve(self%pseudoSequenceObject,fgslStateFile)
   return
  end subroutine diemerKravtsov2014StateRestore
