!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of dark matter halo profile concentrations using the \cite{diemer_universal_2014} algorithm.

  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationDiemerKravtsov2014">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{diemer_universal_2014}.</description>
  !# </darkMatterProfileConcentration>

  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationDiemerKravtsov2014
     !% A dark matter halo profile concentration class implementing the algorithm of \cite{diemer_universal_2014}.
     private
     double precision :: kappa, phi0, phi1, eta0, eta1, alpha, beta
   contains
     procedure :: concentration               => diemerKravtsov2014Concentration
     procedure :: densityContrastDefinition   => diemerKravtsov2014DensityContrastDefinition
     procedure :: darkMatterProfileDefinition => diemerKravtsov2014DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationDiemerKravtsov2014

  interface darkMatterProfileConcentrationDiemerKravtsov2014
     !% Constructors for the {\tt diemerKravtsov2014} dark matter halo profile concentration class.
     module procedure diemerKravtsov2014DefaultConstructor
     module procedure diemerKravtsov2014Constructor
  end interface darkMatterProfileConcentrationDiemerKravtsov2014

  ! Initialization status.
  logical          :: diemerKravtsov2014Initialized=.false.

  ! Default parameters.
  double precision :: darkMatterProfileConcentrationDiemerKravtsov2014Kappa, darkMatterProfileConcentrationDiemerKravtsov2014Phi0 , &
       &              darkMatterProfileConcentrationDiemerKravtsov2014Phi1 , darkMatterProfileConcentrationDiemerKravtsov2014Eta0 , &
       &              darkMatterProfileConcentrationDiemerKravtsov2014Eta1 , darkMatterProfileConcentrationDiemerKravtsov2014Alpha, &
       &              darkMatterProfileConcentrationDiemerKravtsov2014Beta

contains

  function diemerKravtsov2014DefaultConstructor()
    !% Default constructor for the {\tt diemerKravtsov2014} dark matter halo profile concentration class.
    use Input_Parameters
    implicit none
    type(darkMatterProfileConcentrationDiemerKravtsov2014), target :: diemerKravtsov2014DefaultConstructor

    if (.not.diemerKravtsov2014Initialized) then
       !$omp critical(diemerKravtsov2014DefaultInitialize)
       if (.not.diemerKravtsov2014Initialized) then
          ! Get parameters of the model.
          !@ <inputParameter>
          !@   <name>darkMatterProfileConcentrationDiemerKravtsov2014Kappa</name>
          !@   <defaultValue>0.69</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $\kappa$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("darkMatterProfileConcentrationDiemerKravtsov2014Kappa",darkMatterProfileConcentrationDiemerKravtsov2014Kappa,defaultValue=0.69d0)
          !@ <inputParameter>
          !@   <name>darkMatterProfileConcentrationDiemerKravtsov2014Phi0</name>
          !@   <defaultValue>6.58</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $\phi_0$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("darkMatterProfileConcentrationDiemerKravtsov2014Phi0",darkMatterProfileConcentrationDiemerKravtsov2014Phi0,defaultValue=6.58d0)
          !@ <inputParameter>
          !@   <name>darkMatterProfileConcentrationDiemerKravtsov2014Phi1</name>
          !@   <defaultValue>1.37</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $\phi_1$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("darkMatterProfileConcentrationDiemerKravtsov2014Phi1",darkMatterProfileConcentrationDiemerKravtsov2014Phi1,defaultValue=1.37d0)
          !@ <inputParameter>
          !@   <name>darkMatterProfileConcentrationDiemerKravtsov2014Eta0</name>
          !@   <defaultValue>6.82</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $\eta_0$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("darkMatterProfileConcentrationDiemerKravtsov2014Eta0",darkMatterProfileConcentrationDiemerKravtsov2014Eta0,defaultValue=6.82d0)
          !@ <inputParameter>
          !@   <name>darkMatterProfileConcentrationDiemerKravtsov2014Eta1</name>
          !@   <defaultValue>1.42</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $\eta_1$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("darkMatterProfileConcentrationDiemerKravtsov2014Eta1",darkMatterProfileConcentrationDiemerKravtsov2014Eta1,defaultValue=1.42d0)
          !@ <inputParameter>
          !@   <name>darkMatterProfileConcentrationDiemerKravtsov2014Alpha</name>
          !@   <defaultValue>1.12</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $\alpha$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("darkMatterProfileConcentrationDiemerKravtsov2014Alpha",darkMatterProfileConcentrationDiemerKravtsov2014Alpha,defaultValue=1.12d0)
          !@ <inputParameter>
          !@   <name>darkMatterProfileConcentrationDiemerKravtsov2014Beta</name>
          !@   <defaultValue>1.69</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $\beta$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("darkMatterProfileConcentrationDiemerKravtsov2014Beta",darkMatterProfileConcentrationDiemerKravtsov2014Beta,defaultValue=1.69d0)
          ! Record that method is now initialized.
          diemerKravtsov2014Initialized=.true.
       end if
       !$omp end critical(diemerKravtsov2014DefaultInitialize)
    end if
    ! Construct the object.
    diemerKravtsov2014DefaultConstructor                                                         &
         & =diemerKravtsov2014Constructor(                                                       &
         &                                darkMatterProfileConcentrationDiemerKravtsov2014Kappa, &
         &                                darkMatterProfileConcentrationDiemerKravtsov2014Phi0 , &
         &                                darkMatterProfileConcentrationDiemerKravtsov2014Phi1 , &
         &                                darkMatterProfileConcentrationDiemerKravtsov2014Eta0 , &
         &                                darkMatterProfileConcentrationDiemerKravtsov2014Eta1 , &
         &                                darkMatterProfileConcentrationDiemerKravtsov2014Alpha, &
         &                                darkMatterProfileConcentrationDiemerKravtsov2014Beta   &
         &                               )
    return
  end function diemerKravtsov2014DefaultConstructor

  function diemerKravtsov2014Constructor(kappa,phi0,phi1,eta0,eta1,alpha,beta)
    !% Constructor for the {\tt diemerKravtsov2014} dark matter halo profile concentration class.
    use Galacticus_Error
    implicit none
    type            (darkMatterProfileConcentrationDiemerKravtsov2014)                :: diemerKravtsov2014Constructor
    double precision                                                  , intent(in   ) :: kappa, phi0, phi1, eta0, eta1, alpha, beta
    
    diemerKravtsov2014Constructor%kappa=kappa
    diemerKravtsov2014Constructor%phi0 =phi0
    diemerKravtsov2014Constructor%phi1 =phi1
    diemerKravtsov2014Constructor%eta0 =eta0
    diemerKravtsov2014Constructor%eta1 =eta1
    diemerKravtsov2014Constructor%alpha=alpha
    diemerKravtsov2014Constructor%beta =beta
    return
  end function diemerKravtsov2014Constructor
  
  double precision function diemerKravtsov2014Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\tt node} using the \cite{diemer_universal_2014} algorithm.
    use Numerical_Constants_Math
    use Power_Spectra
    use Critical_Overdensity
    use Cosmology_Parameters
    implicit none
    class           (darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout)          :: self
    type            (treeNode                                        ), intent(inout), pointer :: node
    class           (nodeComponentBasic                              )               , pointer :: basic
    class           (cosmologyParametersClass                        )               , pointer :: cosmologyParameters_
    double precision                                                                           :: radiusHaloLagrangian, peakHeight        , &
         &                                                                                        wavenumber          , powerSpectrumSlope, &
         &                                                                                        concentrationMinimum, peakHeightMinimum
    
    cosmologyParameters_           => cosmologyParameters      ()
    basic                          => node               %basic()
    radiusHaloLagrangian           =  +(                                                &
         &                              +3.0d0                                          &
         &                              *basic%mass()                                   &
         &                              /4.0d0                                          &
         &                              /Pi                                             &
         &                              /cosmologyParameters_%densityCritical()         &
         &                              /cosmologyParameters_%OmegaMatter    ()         &
         &                             )**(1.0d0/3.0d0)
    peakHeight                     =                                                    &
         &                            +Critical_Overdensity_for_Collapse(basic%time())  &
         &                            /Cosmological_Mass_Root_Variance  (basic%mass())
    wavenumber                     =                                                    &
         &                            +self%kappa                                       &
         &                            *2.0d0                                            &
         &                            *Pi                                               &
         &                            /radiusHaloLagrangian
    powerSpectrumSlope             = +Power_Spectrum_Logarithmic_Derivative(wavenumber)
    concentrationMinimum           =                                                    &
         &                           +self%phi0                                         &
         &                           +self%phi1                                         &
         &                           *powerSpectrumSlope
    peakHeightMinimum              =                                                    &
         &                           +self%eta0                                         &
         &                           +self%eta1                                         &
         &                           *powerSpectrumSlope
    diemerKravtsov2014Concentration=                                                    &
         &                           +0.5d0                                             &
         &                           *concentrationMinimum                              &
         &                           *(                                                 &
         &                             +(peakHeight/peakHeightMinimum)**(-self%alpha)   &
         &                             +(peakHeight/peakHeightMinimum)**(+self%beta )   &
         &                            )
    return
  end function diemerKravtsov2014Concentration

  function diemerKravtsov2014DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of concentration in the
    !% \cite{diemer_universal_2014} algorithm.
    implicit none
    class(virialDensityContrastClass                      ), pointer       :: diemerKravtsov2014DensityContrastDefinition
    class(darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout) :: self
    
    allocate(virialDensityContrastFixed :: diemerKravtsov2014DensityContrastDefinition)
    select type (diemerKravtsov2014DensityContrastDefinition)
    type is (virialDensityContrastFixed)
       diemerKravtsov2014DensityContrastDefinition=virialDensityContrastFixed(200.0d0,virialDensityContrastFixedDensityTypeCritical)
    end select    
    return
  end function diemerKravtsov2014DensityContrastDefinition

  function diemerKravtsov2014DarkMatterProfileDefinition(self)
    !% Return a dark matter density profile object defining that used in the definition of concentration in the
    !% \cite{diemer_universal_2014} algorithm.
    use Dark_Matter_Halo_Scales
    implicit none
    class(darkMatterProfileClass                            ), pointer       :: diemerKravtsov2014DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationDiemerKravtsov2014  ), intent(inout) :: self
    class(darkMatterHaloScaleVirialDensityContrastDefinition), pointer       :: darkMatterHaloScaleDefinition

    allocate(darkMatterProfileNFW                               :: diemerKravtsov2014DarkMatterProfileDefinition)
    allocate(darkMatterHaloScaleVirialDensityContrastDefinition :: darkMatterHaloScaleDefinition                )
    select type (diemerKravtsov2014DarkMatterProfileDefinition)
    type is (darkMatterProfileNFW)
       select type (darkMatterHaloScaleDefinition)
       type is (darkMatterHaloScaleVirialDensityContrastDefinition)
          darkMatterHaloScaleDefinition                =darkMatterHaloScaleVirialDensityContrastDefinition(self%densityContrastDefinition())
          diemerKravtsov2014DarkMatterProfileDefinition=darkMatterProfileNFW                              (darkMatterHaloScaleDefinition   )
       end select
    end select
    return
  end function diemerKravtsov2014DarkMatterProfileDefinition
