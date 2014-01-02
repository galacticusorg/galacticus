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

!% Contains a module which implements the \cite{navarro_structure_1996} NFW halo concentration algorithm.

module Dark_Matter_Profiles_Concentrations_NFW1996
  !% Implements the \cite{navarro_structure_1996} NFW halo concentration algorithm.
  implicit none
  private
  public :: Dark_Matter_Concentrations_NFW1996_Initialize

  ! Module global variable used in root finding.
  double precision :: targetValue
  !$omp threadprivate(targetValue)
  ! Parameters of the fit.
  double precision :: nfw96ConcentrationC, nfw96ConcentrationF

contains

  !# <darkMatterConcentrationMethod>
  !#  <unitName>Dark_Matter_Concentrations_NFW1996_Initialize</unitName>
  !# </darkMatterConcentrationMethod>
  subroutine Dark_Matter_Concentrations_NFW1996_Initialize(darkMatterConcentrationMethod,Dark_Matter_Profile_Concentration_Get)
    !% Initializes the ``NFW1996'' halo concentration module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                           ), intent(in   )          :: darkMatterConcentrationMethod
    procedure(Dark_Matter_Profile_Concentration_NFW1996), intent(inout), pointer :: Dark_Matter_Profile_Concentration_Get

    if (darkMatterConcentrationMethod == 'NFW1996') then
       ! Return a pointer to our implementation.
       Dark_Matter_Profile_Concentration_Get => Dark_Matter_Profile_Concentration_NFW1996
       ! Get parameters of the model.
       !@ <inputParameter>
       !@   <name>nfw96ConcentrationF</name>
       !@   <defaultValue>0.01 \cite{navarro_structure_1996}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $C$ appearing in the halo concentration algorithm of \cite{navarro_structure_1996}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("nfw96ConcentrationF",nfw96ConcentrationF,defaultValue=   0.01d0)
       !@ <inputParameter>
       !@   <name>nfw96ConcentrationC</name>
       !@   <defaultValue>2000 \cite{navarro_structure_1996}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $f$ appearing in the halo concentration algorithm of \cite{navarro_structure_1996}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("nfw96ConcentrationC",nfw96ConcentrationC,defaultValue=2000.0d0)
    end if
    return
  end subroutine Dark_Matter_Concentrations_NFW1996_Initialize

  double precision function Dark_Matter_Profile_Concentration_NFW1996(thisNode)
    !% Returns the concentration of the dark matter profile of {\tt thisNode} using the method of \cite{navarro_structure_1996}.
    use Galacticus_Nodes
    use Power_Spectra
    use Cosmology_Functions
    use Critical_Overdensity
    use Root_Finder
    use Virial_Density_Contrast
    implicit none
    type            (treeNode               ), intent(inout), pointer :: thisNode
    double precision                         , parameter              :: fitParameterNuHalf         =0.47693628d0
    double precision                         , parameter              :: toleranceAbsolute          =0.0d0       , toleranceRelative      =1.0d-6
    class           (nodeComponentBasic     )               , pointer :: thisBasicComponent
    class           (cosmologyFunctionsClass)               , pointer :: cosmologyFunctionsDefault
    type            (rootFinder             ), save                   :: finder
    !$omp threadprivate(finder)
    double precision                                                  :: collapseCriticalOverdensity             , collapseExpansionFactor       , &
         &                                                               collapseMass                            , collapseOverdensity           , &
         &                                                               collapseTime                            , expansionFactor               , &
         &                                                               nodeMass                                , nodeTime

    ! Get the default cosmology functions object.
    cosmologyFunctionsDefault => cosmologyFunctions()
    ! Get the basic component.
    thisBasicComponent         => thisNode%basic()
    ! Get the properties of the node.
    nodeMass                   =thisBasicComponent%mass()
    nodeTime                   =thisBasicComponent%time()
    expansionFactor            =cosmologyFunctionsDefault%expansionFactor(nodeTime)
    ! Compute the mass of a progenitor as defined by NFW.
    collapseMass               =nfw96ConcentrationF*nodeMass
    ! Find the time of collapse for this progenitor.
    collapseCriticalOverdensity=sqrt(2.0d0*fitParameterNuHalf**2*(Cosmological_Mass_Root_Variance(collapseMass)**2-Cosmological_Mass_Root_Variance(nodeMass)**2))+Critical_Overdensity_for_Collapse(nodeTime)
    collapseTime               =Time_of_Collapse(collapseCriticalOverdensity)
    collapseExpansionFactor    =cosmologyFunctionsDefault%expansionFactor(collapseTime               )
    ! Compute the overdensity of the progenitor at collapse using the scaling given by NFW.
    collapseOverdensity        =nfw96ConcentrationC*(expansionFactor/collapseExpansionFactor)**3
    ! Find the ratio of this overdensity to that at for the present node.
    targetValue                =collapseOverdensity/Halo_Virial_Density_Contrast(nodeTime)
    ! Initialize our root finder.
    if (.not.finder%isInitialized()) then
       call finder%rangeExpand (                                               &
            &                   rangeExpandDownward=0.5d0                    , &
            &                   rangeExpandUpward  =2.0d0                    , &
            &                   rangeExpandType    =rangeExpandMultiplicative  &
            &                  )
       call finder%rootFunction(NFW_Concentration_Function_Root    )
       call finder%tolerance   (toleranceAbsolute,toleranceRelative)
    end if
    ! Find the concentration.
    Dark_Matter_Profile_Concentration_NFW1996=finder%find(rootRange=[1.0d0,20.0d0])
    return
  end function Dark_Matter_Profile_Concentration_NFW1996

  double precision function NFW_Concentration_Function_Root(concentration)
    !% Root function used in finding concentrations in the \cite{navarro_structure_1996} method.
    implicit none
    double precision, intent(in   ) :: concentration

    NFW_Concentration_Function_Root=concentration**3/(log(1.0d0+concentration)-concentration/(1.0d0+concentration))/3.0d0-targetValue
    return
  end function NFW_Concentration_Function_Root

end module Dark_Matter_Profiles_Concentrations_NFW1996
