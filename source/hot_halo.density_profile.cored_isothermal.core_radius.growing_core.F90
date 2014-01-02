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

!% Contains a module which implements a calculation of the core radius in the hot halo density profile that is a fraction of the
!% dark matter profile scale radius, but which grows as gas is depleted to keep the density at the virial radius constant (similar
!% to the algorithm of \citep{cole_hierarchical_2000}).

module Hot_Halo_Density_Cored_Isothermal_Core_Radii_Growing_Core
  !% Implements a calculation of the core radius in the hot halo density profile that is a fraction of the
  !% dark matter profile scale radius, but which grows as gas is depleted to keep the density at the virial radius constant (similar
  !% to the algorithm of \citep{cole_hierarchical_2000}).
  use Galacticus_Nodes
  use Tables
  implicit none
  private
  public :: Hot_Halo_Density_Cored_Isothermal_Core_Radii_GC_Initialize,&
       & Hot_Halo_Density_Cored_Isothermal_Core_Radius_GC_State_Store,&
       & Hot_Halo_Density_Cored_Isothermal_Core_Radius_GC_State_Retrieve

  ! Parameters of the model.
  double precision                                        :: isothermalCoreRadiusOverScaleRadius        , isothermalCoreRadiusOverVirialRadiusMaximum

  ! Minimum and maximum radii (in units of the virial radius) allowed for cores.
  double precision                                        :: coreRadiusMaximum                          , coreRadiusMinimum

  ! Tabulation of core radius vs. virial density relation.
  integer                                   , parameter   :: coreRadiusTablePointsPerDecade     =100
  integer                                                 :: coreRadiusTableCount
  logical                                                 :: coreRadiusTableInitialized         =.false.
  type            (table1DLogarithmicLinear)              :: coreRadiusTable
  class           (table1D                 ), allocatable :: coreRadiusTableInverse

contains

  !# <hotHaloCoredIsothermalCoreRadiiMethod>
  !#  <unitName>Hot_Halo_Density_Cored_Isothermal_Core_Radii_GC_Initialize</unitName>
  !# </hotHaloCoredIsothermalCoreRadiiMethod>
  subroutine Hot_Halo_Density_Cored_Isothermal_Core_Radii_GC_Initialize(hotHaloCoredIsothermalCoreRadiiMethod&
       &,Hot_Halo_Density_Cored_Isothermal_Core_Radius_Get)
    !% Initializes the ``growing core'' cored isothermal hot halo profile core radius module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type     (varying_string                                            ), intent(in   )          :: hotHaloCoredIsothermalCoreRadiiMethod
    procedure(Hot_Halo_Density_Cored_Isothermal_Core_Radius_Growing_Core), intent(inout), pointer :: Hot_Halo_Density_Cored_Isothermal_Core_Radius_Get

    if (hotHaloCoredIsothermalCoreRadiiMethod == 'growingCore') then
       Hot_Halo_Density_Cored_Isothermal_Core_Radius_Get => Hot_Halo_Density_Cored_Isothermal_Core_Radius_Growing_Core
       !@ <inputParameter>
       !@   <name>isothermalCoreRadiusOverScaleRadius</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The core radius in the ``cored isothermal'' hot halo density profile in units of the dark matter profile scale radius.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('isothermalCoreRadiusOverScaleRadius',isothermalCoreRadiusOverScaleRadius,defaultValue=0.1d0)
       !@ <inputParameter>
       !@   <name>isothermalCoreRadiusOverVirialRadiusMaximum</name>
       !@   <defaultValue>10</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum core radius in the ``cored isothermal'' hot halo density profile in units of the virial radius.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('isothermalCoreRadiusOverVirialRadiusMaximum',isothermalCoreRadiusOverVirialRadiusMaximum,defaultValue=10.0d0)
       ! Ensure that the dark matter profile supports the scale property.
       if (.not.defaultDarkMatterProfileComponent%scaleIsGettable())                                                          &
            & call Galacticus_Error_Report                                                                                    &
            &      (                                                                                                          &
            &       'Hot_Halo_Density_Cored_Isothermal_Core_Radii_GC_Initialize'                                            , &
            &       'method requires a dark matter profile component that provides a gettable "scale" property.'//            &
            &       Galacticus_Component_List(                                                                                &
            &                                 'darkMatterProfile'                                                           , &
            &                                  defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)  &
            &                                )                                                                                &
            &      )
    end if
    return
  end subroutine Hot_Halo_Density_Cored_Isothermal_Core_Radii_GC_Initialize

  double precision function Hot_Halo_Density_Cored_Isothermal_Core_Radius_Growing_Core(thisNode)
    !% Returns the radius (in Mpc) of the core radius in a cored isothermal hot halo density profile. Assumes that the radius is
    !% a fraction of the dark matter profile scale radius, but which grows as gas is depleted to keep the density at the virial
    !% radius constant (similar to the algorithm of \citep{cole_hierarchical_2000}).
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: thisNode
    class           (nodeComponentBasic            )               , pointer :: thisBasicComponent
    class           (nodeComponentHotHalo          )               , pointer :: thisHotHaloComponent
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    class           (cosmologyParametersClass      )               , pointer :: thisCosmologyParameters
    double precision                                , save                   :: hotGasFractionSaved                        , isothermalCoreRadiusOverVirialRadiusInitialSaved, &
         &                                                                      isothermalCoreRadiusOverVirialRadiusSaved
    !$omp threadprivate(isothermalCoreRadiusOverVirialRadiusInitialSaved,hotGasFractionSaved,isothermalCoreRadiusOverVirialRadiusSaved)
    double precision                                                         :: hotGasFraction                             , isothermalCoreRadiusOverVirialRadius            , &
         &                                                                      isothermalCoreRadiusOverVirialRadiusInitial, targetValue
    logical                                                                  :: makeTable

    ! Get components.
    thisBasicComponent             => thisNode%basic            ()
    thisHotHaloComponent           => thisNode%hotHalo          ()
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile()

    ! Get the default cosmology.
    thisCosmologyParameters => cosmologyParameters()
    ! Find the fraction of gas in the hot halo relative to that expected from the universal baryon fraction.
    hotGasFraction=(thisHotHaloComponent%mass()/thisBasicComponent%mass())*(thisCosmologyParameters%OmegaMatter()/thisCosmologyParameters%OmegaBaryon())

    ! Return an arbitrary value for empty halos.
    if (hotGasFraction <= 0.0d0) then
       Hot_Halo_Density_Cored_Isothermal_Core_Radius_Growing_Core=isothermalCoreRadiusOverScaleRadius
       return
    end if

    ! Comptue the desired core radius (in units of the virial radius) for a fully populated halo.
    isothermalCoreRadiusOverVirialRadiusInitial=isothermalCoreRadiusOverScaleRadius*thisDarkMatterProfileComponent%scale()&
         &/Dark_Matter_Halo_Virial_Radius(thisNode)

    ! Check if the initial core radius and hot gas fraction equal the previously stored values.
    if (.not.(isothermalCoreRadiusOverVirialRadiusInitial == isothermalCoreRadiusOverVirialRadiusInitialSaved .and.&
         & hotGasFraction == hotGasFractionSaved)) then

       ! Create a tabulation of core radius vs. virial density factor if necessary.
       !$omp critical (Hot_Halo_Density_Growing_Core_Interpolation)
       makeTable=(.not.coreRadiusTableInitialized)
       if (.not.makeTable) makeTable=(isothermalCoreRadiusOverVirialRadiusInitial < coreRadiusTable%x(-1))
       if (makeTable) then
          coreRadiusMinimum   =min(isothermalCoreRadiusOverScaleRadius,isothermalCoreRadiusOverVirialRadiusInitial)
          coreRadiusMaximum   =isothermalCoreRadiusOverVirialRadiusMaximum
          coreRadiusTableCount=int(log10(coreRadiusMaximum/coreRadiusMinimum)*dble(coreRadiusTablePointsPerDecade))+1
          call coreRadiusTable%destroy (                                                          )
          call coreRadiusTable%create  (coreRadiusMaximum,coreRadiusMinimum,coreRadiusTableCount  )
          call coreRadiusTable%populate(Growing_Core_Virial_Density_Function(coreRadiusTable%xs()))
          call coreRadiusTable%reverse (coreRadiusTableInverse                                    )
          coreRadiusTableInitialized=.true.
       end if
       !$omp end critical (Hot_Halo_Density_Growing_Core_Interpolation)

       ! Compute the target value of the function giving the density at the virial radius per unit gas mass.
       targetValue=Growing_Core_Virial_Density_Function(isothermalCoreRadiusOverVirialRadiusInitial)*hotGasFraction

       ! Interpolate to get the required core radius.
       !$omp critical (Hot_Halo_Density_Growing_Core_Interpolation)
       if      (hotGasFraction >= 1.0d0                                             ) then
          isothermalCoreRadiusOverVirialRadius=isothermalCoreRadiusOverVirialRadiusInitial
       else if (targetValue    <= coreRadiusTable%y(1)) then
          isothermalCoreRadiusOverVirialRadius=coreRadiusTable%x(1)
       else
          isothermalCoreRadiusOverVirialRadius=coreRadiusTableInverse%interpolate(targetValue)
       end if
       !$omp end critical (Hot_Halo_Density_Growing_Core_Interpolation)
       isothermalCoreRadiusOverVirialRadiusInitialSaved=isothermalCoreRadiusOverVirialRadiusInitial
       hotGasFractionSaved                             =hotGasFraction
       isothermalCoreRadiusOverVirialRadiusSaved       =isothermalCoreRadiusOverVirialRadius
    end if
    ! Compute the resulting core radius.
    Hot_Halo_Density_Cored_Isothermal_Core_Radius_Growing_Core=isothermalCoreRadiusOverVirialRadiusSaved*Dark_Matter_Halo_Virial_Radius(thisNode)
    return
  end function Hot_Halo_Density_Cored_Isothermal_Core_Radius_Growing_Core

  elemental double precision function Growing_Core_Virial_Density_Function(radiusOverVirialRadius)
    !% Returns the function $(1+r_{\rm c}^2)[1-r_{\rm c} \tan^{-1}(1/r_{\rm c}]$ which is proportional to the density at the
    !% virial radius of a cored isothermal profile with core radius $r_{\rm c}$ (in units of the virial radius) per unit mass.
    implicit none
    double precision, intent(in   ) :: radiusOverVirialRadius

    Growing_Core_Virial_Density_Function=(1.0d0+radiusOverVirialRadius**2)*(1.0d0-radiusOverVirialRadius*atan(1.0d0/radiusOverVirialRadius))
    return
  end function Growing_Core_Virial_Density_Function

  !# <galacticusStateStoreTask>
  !#  <unitName>Hot_Halo_Density_Cored_Isothermal_Core_Radius_GC_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Hot_Halo_Density_Cored_Isothermal_Core_Radius_GC_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    write (stateFile) coreRadiusMinimum,coreRadiusMaximum
    return
  end subroutine Hot_Halo_Density_Cored_Isothermal_Core_Radius_GC_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Hot_Halo_Density_Cored_Isothermal_Core_Radius_GC_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Hot_Halo_Density_Cored_Isothermal_Core_Radius_GC_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) coreRadiusMinimum,coreRadiusMaximum
    ! Force retabulation on next evaluation.
    coreRadiusTableInitialized=.false.
    return
  end subroutine Hot_Halo_Density_Cored_Isothermal_Core_Radius_GC_State_Retrieve

end module Hot_Halo_Density_Cored_Isothermal_Core_Radii_Growing_Core
