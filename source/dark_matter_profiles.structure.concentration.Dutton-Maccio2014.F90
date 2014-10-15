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

  !% An implementation of dark matter halo profile concentrations using the \cite{dutton_cold_2014} algorithm.

  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationDuttonMaccio2014">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{dutton_cold_2014}.</description>
  !# </darkMatterProfileConcentration>

  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationDuttonMaccio2014
     !% A dark matter halo profile concentration class implementing the algorithm of \cite{dutton_cold_2014}.
     private
     double precision :: a1,a2,a3,a4,b1,b2
     integer          :: densityContrastMethod, densityProfileMethod
   contains
     procedure :: concentration               => duttonMaccio2014Concentration
     procedure :: densityContrastDefinition   => duttonMaccio2014DensityContrastDefinition
     procedure :: darkMatterProfileDefinition => duttonMaccio2014DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationDuttonMaccio2014

  interface darkMatterProfileConcentrationDuttonMaccio2014
     !% Constructors for the {\tt duttonMaccio2014} dark matter halo profile concentration class.
     module procedure duttonMaccio2014DefaultConstructor
     module procedure duttonMaccio2014Constructor
     module procedure duttonMaccio2014UserDefinedConstructor
  end interface darkMatterProfileConcentrationDuttonMaccio2014

  ! Initialization status.
  logical                            :: duttonMaccio2014Initialized=.false.

  ! Type of fit used for the concentration-mass relation.
  type   (varying_string)            :: duttonMaccio2014FitType

  ! Parameters for user-defined fit.
  double precision                 :: duttonMaccio2014A1, duttonMaccio2014A2, duttonMaccio2014A3, duttonMaccio2014A4, &
       &                              duttonMaccio2014B1, duttonMaccio2014B2
  ! Density contrast methods.
  integer                , parameter :: duttonMaccio2014DensityContrastMethod200   =0
  integer                , parameter :: duttonMaccio2014DensityContrastMethodVirial=1

  ! Density profile methods.
  integer                , parameter :: duttonMaccio2014DensityProfileMethodNFW    =0
  integer                , parameter :: duttonMaccio2014DensityProfileMethodEinasto=1

contains

  function duttonMaccio2014DefaultConstructor()
    !% Default constructor for the {\tt duttonMaccio2014} dark matter halo profile concentration class.
    use Input_Parameters
    implicit none
    type(darkMatterProfileConcentrationDuttonMaccio2014), target :: duttonMaccio2014DefaultConstructor

    if (.not.duttonMaccio2014Initialized) then
       !$omp critical(duttonMaccio2014DefaultInitialize)
       if (.not.duttonMaccio2014Initialized) then
          ! Get parameters of the model.
          !@ <inputParameter>
          !@   <name>duttonMaccio2014FitType</name>
          !@   <defaultValue>nfwVirial</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Option controlling which of \cite{dutton_cold_2014}'s fits to the halo concentration--mass relation to use.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("duttonMaccio2014FitType",duttonMaccio2014FitType,defaultValue="nfwVirial")
          ! Check for user-defined fit.
          if (duttonMaccio2014FitType == "userDefined") then
             ! Get user-defined parameters.
             !@ <inputParameter>
             !@   <name>duttonMaccio2014A1</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Parameter $a_1$ in the \cite{dutton_cold_2014} halo concentration--mass relation.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter("duttonMaccio2014A1",duttonMaccio2014A1)
             !@ <inputParameter>
             !@   <name>duttonMaccio2014A2</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Parameter $a_2$ in the \cite{dutton_cold_2014} halo concentration--mass relation.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter("duttonMaccio2014A2",duttonMaccio2014A2)
             !@ <inputParameter>
             !@   <name>duttonMaccio2014A3</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Parameter $a_3$ in the \cite{dutton_cold_2014} halo concentration--mass relation.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter("duttonMaccio2014A3",duttonMaccio2014A3)
             !@ <inputParameter>
             !@   <name>duttonMaccio2014A4</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Parameter $a_4$ in the \cite{dutton_cold_2014} halo concentration--mass relation.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter("duttonMaccio2014A4",duttonMaccio2014A4)
             !@ <inputParameter>
             !@   <name>duttonMaccio2014B1</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Parameter $b_1$ in the \cite{dutton_cold_2014} halo concentration--mass relation.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter("duttonMaccio2014B1",duttonMaccio2014B1)
             !@ <inputParameter>
             !@   <name>duttonMaccio2014B2</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Parameter $b_2$ in the \cite{dutton_cold_2014} halo concentration--mass relation.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter("duttonMaccio2014B2",duttonMaccio2014B2)
          end if
          ! Record that method is now initialized.
          duttonMaccio2014Initialized=.true.
       end if
       !$omp end critical(duttonMaccio2014DefaultInitialize)
    end if
    ! Construct the object.
    if (duttonMaccio2014FitType == "userDefined") then
       duttonMaccio2014DefaultConstructor=duttonMaccio2014UserDefinedConstructor(duttonMaccio2014A1,duttonMaccio2014A2,duttonMaccio2014A3,duttonMaccio2014A4,duttonMaccio2014B1,duttonMaccio2014B2)
   else
       duttonMaccio2014DefaultConstructor=duttonMaccio2014Constructor(char(duttonMaccio2014FitType))
    end if
    return
  end function duttonMaccio2014DefaultConstructor

  function duttonMaccio2014Constructor(fitType)
    !% Constructor for the {\tt duttonMaccio2014} dark matter halo profile concentration class.
    use Galacticus_Error
    implicit none
    type     (darkMatterProfileConcentrationDuttonMaccio2014)                :: duttonMaccio2014Constructor
    character(len=*                                         ), intent(in   ) :: fitType

    select case (fitType)
    case ('nfwVirial' )
       duttonMaccio2014Constructor%a1                   =+0.537d0
       duttonMaccio2014Constructor%a2                   =+1.025d0
       duttonMaccio2014Constructor%a3                   =-0.718d0
       duttonMaccio2014Constructor%a4                   =+1.080d0
       duttonMaccio2014Constructor%b1                   =-0.097d0
       duttonMaccio2014Constructor%b2                   =+0.024d0
       duttonMaccio2014Constructor%densityContrastMethod=duttonMaccio2014DensityContrastMethodVirial
       duttonMaccio2014Constructor%densityProfileMethod =duttonMaccio2014DensityProfileMethodNFW
    case ('nfw200'    )
       duttonMaccio2014Constructor%a1                   =+0.520d0
       duttonMaccio2014Constructor%a2                   =+0.905d0
       duttonMaccio2014Constructor%a3                   =-0.617d0
       duttonMaccio2014Constructor%a4                   =+1.210d0
       duttonMaccio2014Constructor%b1                   =-0.101d0
       duttonMaccio2014Constructor%b2                   =+0.026d0
       duttonMaccio2014Constructor%densityContrastMethod=duttonMaccio2014DensityContrastMethod200
       duttonMaccio2014Constructor%densityProfileMethod =duttonMaccio2014DensityProfileMethodNFW
    case ('einasto200')
       duttonMaccio2014Constructor%a1                   =+0.459d0
       duttonMaccio2014Constructor%a2                   =+0.977d0
       duttonMaccio2014Constructor%a3                   =-0.490d0
       duttonMaccio2014Constructor%a4                   =+1.303d0
       duttonMaccio2014Constructor%b1                   =-0.130d0
       duttonMaccio2014Constructor%b2                   =+0.029d0
       duttonMaccio2014Constructor%densityContrastMethod=duttonMaccio2014DensityContrastMethod200
       duttonMaccio2014Constructor%densityProfileMethod =duttonMaccio2014DensityProfileMethodEinasto
    case default
       call Galacticus_Error_Report('duttonMaccio2014Constructor','unrecognized fit type [available types are: nfwVirial, nfw200, einasto200]')
    end select
    return
  end function duttonMaccio2014Constructor
  
  function duttonMaccio2014UserDefinedConstructor(a1,a2,a3,a4,b1,b2)
    !% Constructor for the {\tt duttonMaccio2014} dark matter halo profile concentration class with user defined parameters.
    use Galacticus_Error
    implicit none
    type           (darkMatterProfileConcentrationDuttonMaccio2014)                :: duttonMaccio2014UserDefinedConstructor
    double precision                                               , intent(in   ) :: a1, a2, a3, a4, b1, b2

    duttonMaccio2014UserDefinedConstructor%a1=a1
    duttonMaccio2014UserDefinedConstructor%a2=a2
    duttonMaccio2014UserDefinedConstructor%a3=a3
    duttonMaccio2014UserDefinedConstructor%a4=a4
    duttonMaccio2014UserDefinedConstructor%b1=b1
    duttonMaccio2014UserDefinedConstructor%b2=b2
    return
  end function duttonMaccio2014UserDefinedConstructor
  
  double precision function duttonMaccio2014Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\tt node} using the \cite{dutton_cold_2014} algorithm.
    use Cosmology_Functions
    implicit none
    class           (darkMatterProfileConcentrationDuttonMaccio2014), intent(inout)          :: self
    type            (treeNode                                      ), intent(inout), pointer :: node
    class           (nodeComponentBasic                            )               , pointer :: basic
    class           (cosmologyFunctionsClass                       )               , pointer :: cosmologyFunctions_
    double precision                                                , parameter              :: littleHubbleConstantDuttonMaccio2014= 0.671d0
    double precision                                                , parameter              :: massNormalization                   =12.000d0
    double precision                                                                         :: redshift                                     , logarithmHaloMass, &
         &                                                                                      parameterA                                   , parameterB

    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions()
    ! Get the basic component.
    basic               => node%basic()
    ! Compute the concentration.
    logarithmHaloMass            =+log10(                                                       &
         &                                littleHubbleConstantDuttonMaccio2014                  &
         &                               *basic%mass()                                          &
         &                              )                                                       &
         &                        -massNormalization
    redshift                     =cosmologyFunctions_%redshiftFromExpansionFactor(              &
         &                        cosmologyFunctions_%expansionFactor             (             &
         &                                                                         basic%time() &
         &                                                                        )             &
         &                                                                       )
    parameterA                   =self%a1+(self%a2-self%a1)*exp(self%a3*redshift**self%a4)
    parameterB                   =self%b1+self%b2*redshift
    duttonMaccio2014Concentration=10.0d0**(parameterA+parameterB*logarithmHaloMass)
    return
  end function duttonMaccio2014Concentration

  function duttonMaccio2014DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of concentration in the \cite{dutton_cold_2014} algorithm.
    implicit none
    class(virialDensityContrastClass                    ), pointer       :: duttonMaccio2014DensityContrastDefinition
    class(darkMatterProfileConcentrationDuttonMaccio2014), intent(inout) :: self
    
    select case (self%densityContrastMethod)
    case (duttonMaccio2014DensityContrastMethod200   )
       allocate(virialDensityContrastFixed                         :: duttonMaccio2014DensityContrastDefinition)
       select type (duttonMaccio2014DensityContrastDefinition)
       type is (virialDensityContrastFixed)
          duttonMaccio2014DensityContrastDefinition=virialDensityContrastFixed                        (200.0d0,virialDensityContrastFixedDensityTypeMean)
       end select
    case (duttonMaccio2014DensityContrastMethodVirial)
       allocate(virialDensityContrastSphericalCollapseMatterLambda :: duttonMaccio2014DensityContrastDefinition)
       select type (duttonMaccio2014DensityContrastDefinition)
       type is (virialDensityContrastSphericalCollapseMatterLambda)
          duttonMaccio2014DensityContrastDefinition=virialDensityContrastSphericalCollapseMatterLambda(                                                 )
       end select
    end select
    return
  end function duttonMaccio2014DensityContrastDefinition

  function duttonMaccio2014DarkMatterProfileDefinition(self)
    !% Return a dark matter density profile object defining that used in the definition of concentration in the
    !% \cite{diemer_universal_2014} algorithm.
    use Dark_Matter_Halo_Scales
    implicit none
    class(darkMatterProfileClass                            ), pointer       :: duttonMaccio2014DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationDuttonMaccio2014    ), intent(inout) :: self
    class(darkMatterHaloScaleVirialDensityContrastDefinition), pointer       :: darkMatterHaloScaleDefinition

    allocate(darkMatterHaloScaleVirialDensityContrastDefinition :: darkMatterHaloScaleDefinition)
    select type (darkMatterHaloScaleDefinition)
    type is (darkMatterHaloScaleVirialDensityContrastDefinition)
       select case (self%densityProfileMethod)
       case (duttonMaccio2014DensityProfileMethodNFW    )
          allocate(darkMatterProfileNFW     :: duttonMaccio2014DarkMatterProfileDefinition)
          select type (duttonMaccio2014DarkMatterProfileDefinition)
          type is (darkMatterProfileNFW    )
             darkMatterHaloScaleDefinition              =darkMatterHaloScaleVirialDensityContrastDefinition(self%densityContrastDefinition())
             duttonMaccio2014DarkMatterProfileDefinition=darkMatterProfileNFW                              (darkMatterHaloScaleDefinition   )
          end select
       case (duttonMaccio2014DensityProfileMethodEinasto)
          allocate(darkMatterProfileEinasto :: duttonMaccio2014DarkMatterProfileDefinition)
          select type (duttonMaccio2014DarkMatterProfileDefinition)
          type is (darkMatterProfileEinasto)
             darkMatterHaloScaleDefinition              =darkMatterHaloScaleVirialDensityContrastDefinition(self%densityContrastDefinition())
             duttonMaccio2014DarkMatterProfileDefinition=darkMatterProfileEinasto                          (darkMatterHaloScaleDefinition   )
          end select
       end select
    end select
    return
  end function duttonMaccio2014DarkMatterProfileDefinition
