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

  !% An implementation of dark matter halo profile concentrations using the \cite{navarro_structure_1996} algorithm.

  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationNFW1996">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{navarro_structure_1996}.</description>
  !# </darkMatterProfileConcentration>

  ! Module global variable used in root finding.
  double precision :: nfw1996TargetValue
  !$omp threadprivate(nfw1996TargetValue)
  ! Parameters of the fit.
  logical          :: nfw1996Initialized=.false.
  double precision :: nfw1996ConcentrationC     , nfw1996ConcentrationF

  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationNFW1996
     !% A dark matter halo profile concentration class implementing the algorithm of \cite{navarro_structure_1996}.
     private
     double precision :: F, C
   contains
     procedure :: concentration               => nfw1996Concentration
     procedure :: densityContrastDefinition   => nfw1996DensityContrastDefinition
     procedure :: darkMatterProfileDefinition => nfw1996DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationNFW1996

  interface darkMatterProfileConcentrationNFW1996
     !% Constructors for the {\tt nfw1996} dark matter halo profile concentration class.
     module procedure nfw1996DefaultConstructor
     module procedure nfw1996Constructor
  end interface darkMatterProfileConcentrationNFW1996

contains

  function nfw1996DefaultConstructor()
    !% Default constructor for the {\tt nfw1996} dark matter halo profile concentration class.
    use Input_Parameters
    implicit none
    type(darkMatterProfileConcentrationNFW1996), target :: nfw1996DefaultConstructor

    if (.not.nfw1996Initialized) then
       !$omp critical(nfw1996DefaultInitialize)
       if (.not.nfw1996Initialized) then
          ! Get parameters of the model.
          !@ <inputParameter>
          !@   <name>nfw1996ConcentrationF</name>
          !@   <defaultValue>0.01 \cite{navarro_structure_1996}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $C$ appearing in the halo concentration algorithm of \cite{navarro_structure_1996}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("nfw1996ConcentrationF",nfw1996ConcentrationF,defaultValue=   0.01d0)
          !@ <inputParameter>
          !@   <name>nfw1996ConcentrationC</name>
          !@   <defaultValue>2000 \cite{navarro_structure_1996}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $f$ appearing in the halo concentration algorithm of \cite{navarro_structure_1996}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("nfw1996ConcentrationC",nfw1996ConcentrationC,defaultValue=2000.0d0)
          ! Record that method is now initialized.
          nfw1996Initialized=.true.
       end if
       !$omp end critical(nfw1996DefaultInitialize)
    end if
    ! Construct the object.
    nfw1996DefaultConstructor=nfw1996Constructor(nfw1996ConcentrationF,nfw1996ConcentrationC)
    return
  end function nfw1996DefaultConstructor

  function nfw1996Constructor(F,C)
    !% Constructor for the {\tt nfw1996} dark matter halo profile concentration class.
    implicit none
    type            (darkMatterProfileConcentrationNFW1996)                :: nfw1996Constructor
    double precision                                       , intent(in   ) :: F                 , C

    nfw1996Constructor%F=F
    nfw1996Constructor%C=C
    return
  end function nfw1996Constructor

  double precision function nfw1996Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\tt node} using the \cite{navarro_structure_1996} algorithm.
    use Cosmology_Functions
    use Power_Spectra
    use Critical_Overdensity
    use Root_Finder
    use Virial_Density_Contrast
    implicit none
    class           (darkMatterProfileConcentrationNFW1996), intent(inout)          :: self
    type            (treeNode                             ), intent(inout), pointer :: node
    double precision                                       , parameter              :: fitParameterNuHalf         =0.47693628d0
    double precision                                       , parameter              :: toleranceAbsolute          =0.0d0       , toleranceRelative      =1.0d-6
    class           (nodeComponentBasic                   )               , pointer :: basic
    class           (cosmologyFunctionsClass              )               , pointer :: cosmologyFunctions_
    class           (virialDensityContrastClass           )               , pointer :: virialDensityContrast_
    type            (rootFinder                           ), save                   :: finder
    !$omp threadprivate(finder)
    double precision                                                                :: collapseCriticalOverdensity             , collapseExpansionFactor       , &
         &                                                                             collapseMass                            , collapseOverdensity           , &
         &                                                                             collapseTime                            , expansionFactor               , &
         &                                                                             nodeMass                                , nodeTime

    ! Get default objects.
    cosmologyFunctions_    => cosmologyFunctions   ()
    virialDensityContrast_ => virialDensityContrast()
    ! Get the basic component.
    basic               => node%basic()
    ! Get the properties of the node.
    nodeMass                   =basic%mass()
    nodeTime                   =basic%time()
    expansionFactor            =cosmologyFunctions_%expansionFactor(nodeTime)
    ! Compute the mass of a progenitor as defined by NFW.
    collapseMass               =self%F*nodeMass
    ! Find the time of collapse for this progenitor.
    collapseCriticalOverdensity=sqrt(2.0d0*fitParameterNuHalf**2*(Cosmological_Mass_Root_Variance(collapseMass)**2-Cosmological_Mass_Root_Variance(nodeMass)**2))+Critical_Overdensity_for_Collapse(nodeTime)
    collapseTime               =Time_of_Collapse(collapseCriticalOverdensity)
    collapseExpansionFactor    =cosmologyFunctions_%expansionFactor(collapseTime               )
    ! Compute the overdensity of the progenitor at collapse using the scaling given by NFW.
    collapseOverdensity        =self%C*(expansionFactor/collapseExpansionFactor)**3
    ! Find the ratio of this overdensity to that at for the present node.
    nfw1996TargetValue         =collapseOverdensity/virialDensityContrast_%densityContrast(nodeTime)
    ! Initialize our root finder.
    if (.not.finder%isInitialized()) then
       call finder%rangeExpand (                                               &
            &                   rangeExpandDownward=0.5d0                    , &
            &                   rangeExpandUpward  =2.0d0                    , &
            &                   rangeExpandType    =rangeExpandMultiplicative  &
            &                  )
       call finder%rootFunction(nfw1996RootFunction                )
       call finder%tolerance   (toleranceAbsolute,toleranceRelative)
    end if
    ! Find the concentration.
    nfw1996Concentration=finder%find(rootRange=[1.0d0,20.0d0])
    return
  end function nfw1996Concentration

  double precision function nfw1996RootFunction(concentration)
    !% Root function used in finding concentrations in the \cite{navarro_structure_1996} method.
    implicit none
    double precision, intent(in   ) :: concentration

    nfw1996RootFunction=concentration**3/(log(1.0d0+concentration)-concentration/(1.0d0+concentration))/3.0d0-nfw1996TargetValue
    return
  end function nfw1996RootFunction

  function nfw1996DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of concentration in the \cite{navarro_structure_1996} algorithm.
    implicit none
    class(virialDensityContrastClass           ), pointer       :: nfw1996DensityContrastDefinition
    class(darkMatterProfileConcentrationNfw1996), intent(inout) :: self
    
    allocate(virialDensityContrastFixed :: nfw1996DensityContrastDefinition)
    select type (nfw1996DensityContrastDefinition)
    type is (virialDensityContrastFixed)
      nfw1996DensityContrastDefinition=virialDensityContrastFixed(200.0d0,virialDensityContrastFixedDensityTypeCritical)
    end select
    return
  end function nfw1996DensityContrastDefinition
  
  function nfw1996DarkMatterProfileDefinition(self)
    !% Return a dark matter density profile object defining that used in the definition of concentration in the
    !% \cite{navarro_structure_1996} algorithm.
    use Dark_Matter_Halo_Scales
    implicit none
    class(darkMatterProfileClass                            ), pointer       :: nfw1996DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationNFW1996             ), intent(inout) :: self
    class(darkMatterHaloScaleVirialDensityContrastDefinition), pointer       :: darkMatterHaloScaleDefinition

    allocate(darkMatterProfileNFW                               :: nfw1996DarkMatterProfileDefinition)
    allocate(darkMatterHaloScaleVirialDensityContrastDefinition :: darkMatterHaloScaleDefinition     )
    select type (nfw1996DarkMatterProfileDefinition)
    type is (darkMatterProfileNFW)
       select type (darkMatterHaloScaleDefinition)
       type is (darkMatterHaloScaleVirialDensityContrastDefinition)
          darkMatterHaloScaleDefinition     =darkMatterHaloScaleVirialDensityContrastDefinition(self%densityContrastDefinition())
          nfw1996DarkMatterProfileDefinition=darkMatterProfileNFW                              (darkMatterHaloScaleDefinition   )
       end select
    end select
    return
  end function nfw1996DarkMatterProfileDefinition

