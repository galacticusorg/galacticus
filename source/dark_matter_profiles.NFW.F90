!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements NFW halo profiles.

module Dark_Matter_Profiles_NFW
  !% Implements NFW halo profiles.
  use Galacticus_Nodes
  use FGSL
  use Kind_Numbers
  use Tables
  implicit none
  private
  public :: Dark_Matter_Profile_NFW_Initialize, Dark_Matter_Profiles_NFW_State_Store, Dark_Matter_Profiles_NFW_State_Retrieve,&
       & Dark_Matter_Profile_NFW_Reset

  ! Minimum and maximum concentrations to tabulate.
  double precision                            :: concentrationMinimum= 1.0d0
  double precision                            :: concentrationMaximum=20.0d0
  ! Minimum and maximum radii to tabulate.
  double precision                            :: radiusMinimum       = 1.0d-3,freefallRadiusMinimum=1.0d-3
  double precision                            :: radiusMaximum       = 1.0d+2,freefallRadiusMaximum=1.0d+2
  double precision                            :: specificAngularMomentumMinimum,freefallTimeMinimum
  double precision                            :: specificAngularMomentumMaximum,freefallTimeMaximum
  ! Number of points per decade of concentration in NFW tabulations.
  integer,          parameter                 :: nfwTablePointsPerDecade        =100
  integer,          parameter                 :: nfwInverseTablePointsPerDecade =100
  integer,          parameter                 :: nfwFreefallTablePointsPerDecade=100
  ! Tables of NFW properties.
  logical                                     :: nfwTableInitialized=.false.,nfwInverseTableInitialized=.false.,nfwFreefallTableInitialized=.false.
  integer                                     :: nfwTableNumberPoints,nfwInverseTableNumberPoints,nfwFreefallTableNumberPoints
  double precision, allocatable, dimension(:) :: nfwRadius ,nfwSpecificAngularMomentum,nfwFreefallTime,nfwFreefallRadius
  integer                 , parameter         :: nfwConcentrationEnergyIndex=1,nfwConcetrationRotationNormalizationIndex=2
  type(table1DLogarithmicLinear)                    :: nfwConcentrationTable

  ! Interpolator variables.
  type(fgsl_interp)                           :: interpolationInverseObject      ,interpolationFreefallObject
  type(fgsl_interp_accel)                     :: interpolationInverseAccelerator ,interpolationFreefallAccelerator
  logical                                     :: interpolationInverseReset=.true.,interpolationFreefallReset=.true.

  ! Module variables used in integrations.
  double precision                            :: concentrationParameter,radiusStart

  ! Record of unique ID of node which we last computed results for.
  integer(kind=kind_int8)                     :: lastUniqueID=-1
  !$omp threadprivate(lastUniqueID)

  ! Record of whether or not quantities have been computed.
  logical :: specificAngularMomentumScalingsComputed=.false.
  !$omp threadprivate(specificAngularMomentumScalingsComputed)

  ! Stored values of computed quantities.
  double precision :: specificAngularMomentumLengthScale,specificAngularMomentumScale
  !$omp threadprivate(specificAngularMomentumLengthScale,specificAngularMomentumScale)

contains

  !# <darkMatterProfileMethod>
  !#  <unitName>Dark_Matter_Profile_NFW_Initialize</unitName>
  !# </darkMatterProfileMethod>
  subroutine Dark_Matter_Profile_NFW_Initialize(darkMatterProfileMethod,Dark_Matter_Profile_Density_Get&
       &,Dark_Matter_Profile_Energy_Get ,Dark_Matter_Profile_Energy_Growth_Rate_Get&
       &,Dark_Matter_Profile_Rotation_Normalization_Get ,Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get&
       &,Dark_Matter_Profile_Circular_Velocity_Get ,Dark_Matter_Profile_Potential_Get,Dark_Matter_Profile_Enclosed_Mass_Get&
       &,Dark_Matter_Profile_kSpace_Get,Dark_Matter_Profile_Freefall_Radius_Get &
       &,Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Get)
    !% Initializes the ``NFW'' halo profile module.
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    type(varying_string),                 intent(in)    :: darkMatterProfileMethod
    procedure(Dark_Matter_Profile_Density_NFW), pointer, intent(inout) :: Dark_Matter_Profile_Density_Get
    procedure(Dark_Matter_Profile_Energy_NFW), pointer, intent(inout) :: Dark_Matter_Profile_Energy_Get
    procedure(Dark_Matter_Profile_Energy_Growth_Rate_NFW), pointer, intent(inout) :: Dark_Matter_Profile_Energy_Growth_Rate_Get 
    procedure(Dark_Matter_Profile_Rotation_Normalization_NFW), pointer, intent(inout) :: Dark_Matter_Profile_Rotation_Normalization_Get
    procedure(Radius_from_Specific_Angular_Momentum_NFW), pointer, intent(inout) :: Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get 
    procedure(Dark_Matter_Profile_Circular_Velocity_NFW), pointer, intent(inout) :: Dark_Matter_Profile_Circular_Velocity_Get
    procedure(Dark_Matter_Profile_Potential_NFW), pointer, intent(inout) :: Dark_Matter_Profile_Potential_Get
    procedure(Dark_Matter_Profile_Enclosed_Mass_NFW), pointer, intent(inout) :: Dark_Matter_Profile_Enclosed_Mass_Get 
    procedure(Dark_Matter_Profile_kSpace_NFW), pointer, intent(inout) :: Dark_Matter_Profile_kSpace_Get
    procedure(Dark_Matter_Profile_Freefall_Radius_NFW), pointer, intent(inout) :: Dark_Matter_Profile_Freefall_Radius_Get 
    procedure(Dark_Matter_Profile_Freefall_Radius_Increase_Rate_NFW), pointer, intent(inout) :: Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Get
    
    if (darkMatterProfileMethod == 'NFW') then
       Dark_Matter_Profile_Density_Get                               => Dark_Matter_Profile_Density_NFW
       Dark_Matter_Profile_Energy_Get                                => Dark_Matter_Profile_Energy_NFW
       Dark_Matter_Profile_Energy_Growth_Rate_Get                    => Dark_Matter_Profile_Energy_Growth_Rate_NFW
       Dark_Matter_Profile_Rotation_Normalization_Get                => Dark_Matter_Profile_Rotation_Normalization_NFW
       Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get => Radius_from_Specific_Angular_Momentum_NFW
       Dark_Matter_Profile_Circular_Velocity_Get                     => Dark_Matter_Profile_Circular_Velocity_NFW
       Dark_Matter_Profile_Potential_Get                             => Dark_Matter_Profile_Potential_NFW
       Dark_Matter_Profile_Enclosed_Mass_Get                         => Dark_Matter_Profile_Enclosed_Mass_NFW       
       Dark_Matter_Profile_kSpace_Get                                => Dark_Matter_Profile_kSpace_NFW
       Dark_Matter_Profile_Freefall_Radius_Get                       => Dark_Matter_Profile_Freefall_Radius_NFW
       Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Get         => Dark_Matter_Profile_Freefall_Radius_Increase_Rate_NFW
        ! Ensure that the dark matter profile component supports a "scale" property. Since we've been called with a treeNode to
       ! process, it should have been initialized by now.
       if (.not.defaultDarkMatterProfileComponent%scaleIsGettable()) call&
            & Galacticus_Error_Report('Dark_Matter_Profile_NFW_Initialize','NFW dark matter profile requires a dark matter&
            & profile component with a gettable "scale" property')
       ! Initialize the tabulations.
       call Dark_Matter_Profile_NFW_Tabulate
       call Dark_Matter_Profile_NFW_Inverse_Angular_Momentum
    end if
    return
  end subroutine Dark_Matter_Profile_NFW_Initialize

  !# <calculationResetTask>
  !# <unitName>Dark_Matter_Profile_NFW_Reset</unitName>
  !# </calculationResetTask>
  subroutine Dark_Matter_Profile_NFW_Reset(thisNode)
    !% Reset the cooling radius calculation.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    specificAngularMomentumScalingsComputed=.false.
    lastUniqueID                           =thisNode%uniqueID()
    return
  end subroutine Dark_Matter_Profile_NFW_Reset

  subroutine Dark_Matter_Profile_NFW_Tabulate(concentration)
    !% Tabulate properties of the NFW halo profile which must be computed numerically.
    use Memory_Management
    use Numerical_Interpolation
    implicit none
    double precision, intent(in), optional :: concentration
    integer                                :: iConcentration
    logical                                :: retabulate
    double precision                       :: tableConcentration

    !$omp critical (NFW_Interpolation)
    retabulate=.not.nfwTableInitialized
    if (present(concentration)) then
       if (concentration < concentrationMinimum) then
          concentrationMinimum=0.5d0*concentration
          retabulate=.true.
       end if
       if (concentration > concentrationMaximum) then
          concentrationMaximum=2.0d0*concentration
          retabulate=.true.
       end if
    end if
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       nfwTableNumberPoints=int(dlog10(concentrationMaximum/concentrationMinimum)*dble(nfwTablePointsPerDecade))+1
       call nfwConcentrationTable%destroy()
       call nfwConcentrationTable%create(concentrationMinimum,concentrationMaximum,nfwTableNumberPoints,2)
       ! Loop over concentrations and populate tables.
       do iConcentration=1,nfwTableNumberPoints
          tableConcentration=nfwConcentrationTable%x(iConcentration)
          call nfwConcentrationTable%populate(NFW_Profile_Energy(tableConcentration),iConcentration,table&
               &=nfwConcentrationEnergyIndex)
          call nfwConcentrationTable%populate(tableConcentration /Angular_Momentum_NFW_Scale_Free(tableConcentration)&
               &,iConcentration,table=nfwConcetrationRotationNormalizationIndex)
       end do
       ! Specify that tabulation has been made.
       nfwTableInitialized=.true.
    end if
    !$omp end critical (NFW_Interpolation)
    return
  end subroutine Dark_Matter_Profile_NFW_Tabulate

  subroutine Dark_Matter_Profile_NFW_Inverse_Angular_Momentum(specificAngularMomentum)
    !% Tabulates the specific angular momentum vs. radius in an NFW profile for rapid inversion.
    use Numerical_Ranges
    use Memory_Management
    use Numerical_Interpolation
    implicit none
    double precision, intent(in), optional :: specificAngularMomentum
    integer                                :: iRadius
    logical                                :: retabulate

    !$omp critical (NFW_Inverse_Interpolation)
    retabulate=.not.nfwInverseTableInitialized
    ! If the table has not yet been made, compute and store the specific angular momenta corresponding to the minimum and maximum
    ! radii that will be tabulated by default.
    if (retabulate) then
       specificAngularMomentumMinimum=Specific_Angular_Momentum_NFW_Scale_Free(radiusMinimum)
       specificAngularMomentumMaximum=Specific_Angular_Momentum_NFW_Scale_Free(radiusMaximum)
    end if
    if (present(specificAngularMomentum)) then
       do while (specificAngularMomentum < specificAngularMomentumMinimum)
          radiusMinimum=0.5d0*radiusMinimum
          specificAngularMomentumMinimum=Specific_Angular_Momentum_NFW_Scale_Free(radiusMinimum)
          retabulate=.true.
       end do
       do while (specificAngularMomentum > specificAngularMomentumMaximum)
          radiusMaximum=2.0d0*radiusMaximum
          specificAngularMomentumMaximum=Specific_Angular_Momentum_NFW_Scale_Free(radiusMaximum)
          retabulate=.true.
       end do
    end if
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       nfwInverseTableNumberPoints=int(dlog10(radiusMaximum/radiusMinimum)*dble(nfwInverseTablePointsPerDecade))+1
       if (allocated(nfwRadius)) then
          call Dealloc_Array(nfwRadius                 )
          call Dealloc_Array(nfwSpecificAngularMomentum)
       end if
       call Alloc_Array(nfwRadius                 ,[nfwInverseTableNumberPoints])
       call Alloc_Array(nfwSpecificAngularMomentum,[nfwInverseTableNumberPoints])
       ! Create a range of radii.
       nfwRadius=Make_Range(radiusMinimum,radiusMaximum,nfwInverseTableNumberPoints,rangeType=rangeTypeLogarithmic)
       ! Loop over radii and populate tables.
       do iRadius=1,nfwInverseTableNumberPoints
          nfwSpecificAngularMomentum(iRadius)=Specific_Angular_Momentum_NFW_Scale_Free(nfwRadius(iRadius))
       end do
       ! Ensure interpolations get reset.
       call Interpolate_Done(interpolationInverseObject,interpolationInverseAccelerator,interpolationInverseReset)
       interpolationInverseReset=.true.
       ! Specify that tabulation has been made.
       nfwInverseTableInitialized=.true.
    end if
    !$omp end critical (NFW_Inverse_Interpolation)
    return
  end subroutine Dark_Matter_Profile_NFW_Inverse_Angular_Momentum

  double precision function Dark_Matter_Profile_Density_NFW(thisNode,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given
    !% in units of Mpc).
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    double precision                     , intent(in   )          :: radius
    class(nodeComponentBasic            ),                pointer :: thisBasicComponent
    class(nodeComponentDarkMatterProfile),                pointer :: thisDarkMatterProfileComponent
    double precision                                              :: scaleRadius,radiusOverScaleRadius,virialRadiusOverScaleRadius

    thisBasicComponent             => thisNode%basic            (                 )
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)
    scaleRadius                    =thisDarkMatterProfileComponent%scale()
    radiusOverScaleRadius          =radius                                  /scaleRadius
    virialRadiusOverScaleRadius    =Dark_Matter_Halo_Virial_Radius(thisNode)/scaleRadius
    Dark_Matter_Profile_Density_NFW=Density_NFW_Scale_Free(radiusOverScaleRadius,virialRadiusOverScaleRadius)&
         &*thisBasicComponent%mass()/scaleRadius**3
    return
  end function Dark_Matter_Profile_Density_NFW
  
  double precision function Dark_Matter_Profile_Enclosed_Mass_NFW(thisNode,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc).
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type (treeNode)                      , intent(inout), pointer :: thisNode
    double precision                     , intent(in    )         :: radius
    class(nodeComponentBasic            ),                pointer :: thisBasicComponent
    class(nodeComponentDarkMatterProfile),                pointer :: thisDarkMatterProfileComponent
    double precision                                              :: scaleRadius,radiusOverScaleRadius,virialRadiusOverScaleRadius

    thisBasicComponent             => thisNode%basic            (                 )
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)
    scaleRadius                    =thisDarkMatterProfileComponent%scale()
    radiusOverScaleRadius          =radius                                  /scaleRadius
    virialRadiusOverScaleRadius    =Dark_Matter_Halo_Virial_Radius(thisNode)/scaleRadius
    Dark_Matter_Profile_Enclosed_Mass_NFW=Enclosed_Mass_NFW_Scale_Free(radiusOverScaleRadius,virialRadiusOverScaleRadius) &
         &*thisBasicComponent%mass()
    return
  end function Dark_Matter_Profile_Enclosed_Mass_NFW

  double precision function Dark_Matter_Profile_Potential_NFW(thisNode,radius)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc).
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    double precision                     , intent(in   )          :: radius
    class(nodeComponentDarkMatterProfile),                pointer :: thisDarkMatterProfileComponent
    double precision                     , parameter              :: radiusSmall=1.0d-10
    double precision                                              :: radiusOverScaleRadius,virialRadiusOverScaleRadius,radiusTerm

    thisDarkMatterProfileComponent   => thisNode%darkMatterProfile(autoCreate=.true.)
    radiusOverScaleRadius            =radius                                  /thisDarkMatterProfileComponent%scale()
    virialRadiusOverScaleRadius      =Dark_Matter_Halo_Virial_Radius(thisNode)/thisDarkMatterProfileComponent%scale()
    if (radiusOverScaleRadius < radiusSmall) then
       ! Use a series solution for very small radii.
       radiusTerm=1.0d0-0.5d0*radiusOverScaleRadius
    else
       ! Use the full expression for larger radii.
       radiusTerm=dlog(1.0d0+radiusOverScaleRadius)/radiusOverScaleRadius
    end if
    Dark_Matter_Profile_Potential_NFW=(-1.0d0-virialRadiusOverScaleRadius*(radiusTerm-dlog(1.0d0+virialRadiusOverScaleRadius)&
         &/virialRadiusOverScaleRadius)/(dlog(1.0d0 +virialRadiusOverScaleRadius)-virialRadiusOverScaleRadius/(1.0d0&
         &+virialRadiusOverScaleRadius))) *Dark_Matter_Halo_Virial_Velocity(thisNode)**2
    return
  end function Dark_Matter_Profile_Potential_NFW
  
  double precision function Dark_Matter_Profile_Circular_Velocity_NFW(thisNode,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc). For an NFW halo this is independent of radius and therefore equal to the virial velocity.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Physical
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius

    if (radius > 0.0d0) then
       Dark_Matter_Profile_Circular_Velocity_NFW=dsqrt(gravitationalConstantGalacticus&
            &*Dark_Matter_Profile_Enclosed_Mass_NFW(thisNode,radius)/radius)
    else
       Dark_Matter_Profile_Circular_Velocity_NFW=0.0d0
    end if
    return
  end function Dark_Matter_Profile_Circular_Velocity_NFW
  
  double precision function Radius_from_Specific_Angular_Momentum_NFW(thisNode,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\tt thisNode} at which a circular orbit has the given {\tt specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc). For an NFW halo, the circular velocity is constant (and therefore equal to the virial
    !% velocity). Therefore, $r = j/V_{\rm virial}$ where $j$(={\tt specificAngularMomentum}) is the specific angular momentum and
    !% $r$ the required radius.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    double precision                     , intent(in)             :: specificAngularMomentum
    class(nodeComponentDarkMatterProfile),                pointer :: thisDarkMatterProfileComponent
    double precision                                             :: specificAngularMomentumScaleFree

    ! Return immediately with zero radius for non-positive specific angular momenta.
    if (specificAngularMomentum <= 0.0d0) then
       Radius_from_Specific_Angular_Momentum_NFW=0.0d0
       return
    end if

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Dark_Matter_Profile_NFW_Reset(thisNode)

    ! Check if scalings are already computed. Compute and store if not.
    if (.not.specificAngularMomentumScalingsComputed) then
       ! Flag that scale quantities are now computed.
       specificAngularMomentumScalingsComputed=.true.
      
       ! Get the dark matter profile.
       thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

       ! Get the scale radius.
       specificAngularMomentumLengthScale=thisDarkMatterProfileComponent%scale()

       ! Get the specific angular momentum scale.
       specificAngularMomentumScale=specificAngularMomentumLengthScale*Dark_Matter_Profile_Circular_Velocity_NFW(thisNode&
            &,specificAngularMomentumLengthScale)
    end if

    ! Compute the specific angular momentum in scale free units (using the scale length for distances the sqrt(G M(r_scale) /
    ! r_scale) for velocities).
    specificAngularMomentumScaleFree=specificAngularMomentum/specificAngularMomentumScale

    ! Ensure that the interpolations exist and extend sufficiently far.
    call Dark_Matter_Profile_NFW_Inverse_Angular_Momentum(specificAngularMomentumScaleFree)

    ! Interpolate to get the dimensionless radius at which this specific angular momentum is found.
    !$omp critical(NFW_Inverse_Interpolation)
    Radius_from_Specific_Angular_Momentum_NFW=Interpolate(nfwInverseTableNumberPoints,nfwSpecificAngularMomentum,nfwRadius &
         &,interpolationInverseObject,interpolationInverseAccelerator,specificAngularMomentumScaleFree,reset &
         &=interpolationInverseReset)
    !$omp end critical(NFW_Inverse_Interpolation)

    ! Convert to a physical radius.
    Radius_from_Specific_Angular_Momentum_NFW=Radius_from_Specific_Angular_Momentum_NFW*specificAngularMomentumLengthScale

    return
  end function Radius_from_Specific_Angular_Momentum_NFW
  
  double precision function Dark_Matter_Profile_Rotation_Normalization_NFW(thisNode)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    class(nodeComponentDarkMatterProfile),                pointer :: thisDarkMatterProfileComponent
    double precision                                              :: concentration

    ! Get components.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    ! Find the concentration parameter of this halo.
    concentration=Dark_Matter_Halo_Virial_Radius(thisNode)/thisDarkMatterProfileComponent%scale()

    ! Ensure that the interpolations exist and extend sufficiently far.
    call Dark_Matter_Profile_NFW_Tabulate(concentration)

    ! Find the energy by interpolation.
    !$omp critical(NFW_Interpolation)
    Dark_Matter_Profile_Rotation_Normalization_NFW=nfwConcentrationTable%interpolate(concentration,table&
         &=nfwConcentrationEnergyIndex)/Dark_Matter_Halo_Virial_Radius(thisNode)
    !$omp end critical(NFW_Interpolation)
    return
  end function Dark_Matter_Profile_Rotation_Normalization_NFW
  
  double precision function Dark_Matter_Profile_Energy_NFW(thisNode)
    !% Return the energy of an NFW halo density profile.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    class(nodeComponentDarkMatterProfile),                pointer :: thisDarkMatterProfileComponent
    class(nodeComponentBasic            ),                pointer :: thisBasicComponent
    double precision                                              :: concentration

    ! Get components.
    thisBasicComponent             => thisNode%basic            (                 )
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    ! Find the concentration parameter of this halo.
    concentration=Dark_Matter_Halo_Virial_Radius(thisNode)/thisDarkMatterProfileComponent%scale()

    ! Ensure that the interpolations exist and extend sufficiently far.
    call Dark_Matter_Profile_NFW_Tabulate(concentration)

    ! Find the energy by interpolation.
    !$omp critical(NFW_Interpolation)
    Dark_Matter_Profile_Energy_NFW=nfwConcentrationTable%interpolate(concentration,table=nfwConcentrationEnergyIndex)&
         &*thisBasicComponent%mass()*Dark_Matter_Halo_Virial_Velocity(thisNode)**2
    !$omp end critical(NFW_Interpolation)
    return
  end function Dark_Matter_Profile_Energy_NFW
  
  double precision function Dark_Matter_Profile_Energy_Growth_Rate_NFW(thisNode)
    !% Return the rate of change of the energy of an NFW halo density profile.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    class(nodeComponentDarkMatterProfile),                pointer :: thisDarkMatterProfileComponent
    class(nodeComponentBasic            ),                pointer :: thisBasicComponent
    double precision                                              :: concentration,energy,energyGradient

    ! Get components.
    thisBasicComponent             => thisNode%basic            (                 )
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    ! Find the concentration parameter of this halo.
    concentration=Dark_Matter_Halo_Virial_Radius(thisNode)/thisDarkMatterProfileComponent%scale()

    ! Ensure that the interpolations exist and extend sufficiently far.
    call Dark_Matter_Profile_NFW_Tabulate(concentration)

    ! Find the energy gradient by interpolation.
    !$omp critical(NFW_Interpolation)
    energy        =nfwConcentrationTable%interpolate        (concentration,table=nfwConcentrationEnergyIndex)
    energyGradient=nfwConcentrationTable%interpolateGradient(concentration,table=nfwConcentrationEnergyIndex)
    !$omp end critical(NFW_Interpolation)

    Dark_Matter_Profile_Energy_Growth_Rate_NFW=Dark_Matter_Profile_Energy_NFW(thisNode)&
         &*(thisBasicComponent%accretionRate()/thisBasicComponent%mass()+2.0d0 &
         &*Dark_Matter_Halo_Virial_Velocity_Growth_Rate(thisNode)/Dark_Matter_Halo_Virial_Velocity(thisNode)+(energyGradient&
         &*concentration/energy)*(Dark_Matter_Halo_Virial_Radius_Growth_Rate(thisNode)/Dark_Matter_Halo_Virial_Radius(thisNode)&
         &-thisDarkMatterProfileComponent%scaleGrowthRate()/thisDarkMatterProfileComponent%scale()))

    return
  end function Dark_Matter_Profile_Energy_Growth_Rate_NFW

  double precision function Angular_Momentum_NFW_Scale_Free(concentration)
    !% Returns the total angular momentum (in units of the virial mass times scale radius times [assumed constant] rotation speed)
    !% in an NFW dark matter profile with given {\tt concentration}. This is given by:
    !% \begin{equation}
    !% J = \left. \int_0^c 4 \pi x^3 \rho(x) \d x \right/ \int_0^c 4 \pi x^2 \rho(x) \d x,
    !% \end{equation}
    !% where $x$ is radius in units of the scale radius and $c$ is concentration. This can be evaluated to give
    !% \begin{equation}
    !% J = \left. \left[ 1 + c - 2 \ln (1+c) - {1 \over 1+c} \right] \right/ \left[ \ln(1+c)-{c\over 1+c} \right].
    !% \end{equation}
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in) :: concentration

    Angular_Momentum_NFW_Scale_Free=(1.0d0+concentration-2.0d0*dlog(1.0d0+concentration)-1.0d0/(1.0d0+concentration)) &
         &/(dlog(1.0d0+concentration)-concentration/(1.0d0+concentration))
    return
  end function Angular_Momentum_NFW_Scale_Free

  double precision function Specific_Angular_Momentum_NFW_Scale_Free(radius)
    !% Returns the specific angular momentum, normalized to unit scale length and unit velocity at the scale radius, at position
    !% {\tt radius} (in units of the scale radius) in an NFW profile.
    implicit none
    double precision, intent(in) :: radius

    Specific_Angular_Momentum_NFW_Scale_Free=dsqrt(radius*Enclosed_Mass_NFW_Scale_Free(radius,1.0d0))
    return
  end function Specific_Angular_Momentum_NFW_Scale_Free

  double precision function Enclosed_Mass_NFW_Scale_Free(radius,concentration)
    !% Returns the enclosed mass (in units of the virial mass) in an NFW dark matter profile with given {\tt concentration} at the
    !% given {\tt radius} (given in units of the scale radius).
    implicit none
    double precision, intent(in) :: radius,concentration
    double precision, parameter  :: minimumRadiusForExactSolution=1.0d-7
    ! Precomputed NFW normalization factor for unit concentration.
    double precision, parameter  :: nfwNormalizationFactorUnitConcentration=1.0d0/(dlog(2.0d0)-0.5d0)
    ! Precomputed NFW normalization factor for unit radius.
    double precision, parameter  :: nfwNormalizationFactorUnitRadius       =dlog(2.0d0)-0.5d0
    double precision, save       :: concentrationPrevious=-1.0d0,nfwNormalizationFactorPrevious
    !$omp threadprivate(concentrationPrevious,nfwNormalizationFactorPrevious)

    if (radius == 1.0d0) then
       Enclosed_Mass_NFW_Scale_Free=nfwNormalizationFactorUnitRadius
    else if (radius >= minimumRadiusForExactSolution) then
       Enclosed_Mass_NFW_Scale_Free=(dlog(1.0d0+radius)-radius/(1.0d0+radius))
    else
       Enclosed_Mass_NFW_Scale_Free=(radius**2)*(0.5d0+radius*(-2.0d0/3.0d0+radius*(0.75d0+radius*(-0.8d0))))
    end if
    ! Check if we were called with a different concentration compared to the previous call.
    if (concentration /= concentrationPrevious) then
       ! We were, so recompute the normalization factor.
       if (concentration == 1.0d0) then
          nfwNormalizationFactorPrevious=nfwNormalizationFactorUnitConcentration
       else
          nfwNormalizationFactorPrevious=1.0d0/(dlog(1.0d0+concentration)-concentration/(1.0d0 &
               &+concentration))
       end if
       concentrationPrevious=concentration
    end if
    Enclosed_Mass_NFW_Scale_Free=Enclosed_Mass_NFW_Scale_Free*nfwNormalizationFactorPrevious
    return
  end function Enclosed_Mass_NFW_Scale_Free
  
  double precision function Density_NFW_Scale_Free(radius,concentration)
    !% Returns the density (in units such that the virial mass and scale length are unity) in an NFW dark matter profile with
    !% given {\tt concentration} at the given {\tt radius} (given in units of the scale radius).
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in) :: radius,concentration

    Density_NFW_Scale_Free=1.0d0/(dlog(1.0d0+concentration)-concentration/(1.0d0+concentration))/radius/(1.0d0+radius)**2/4.0d0/Pi
    return
  end function Density_NFW_Scale_Free
  
  double precision function NFW_Profile_Energy(concentration)
    !% Computes the total energy of an NFW profile halo of given {\tt concentration} using the methods of
    !% \citeauthor{cole_hierarchical_2000}~(\citeyear{cole_hierarchical_2000}; their Appendix~A).
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Math
    use Numerical_Integration
    implicit none
    double precision,                intent(in) :: concentration
    type(c_ptr)                                 :: parameterPointer
    type(fgsl_function)                         :: integrandFunction
    type(fgsl_integration_workspace)            :: integrationWorkspace
    double precision                            :: radiusMinimum,radiusMaximum,potentialEnergyIntegral,potentialEnergy&
         &,kineticEnergyIntegral,kineticEnergy,jeansEquationIntegral

    ! Compute the potential energy.
    radiusMinimum=0.0d0
    radiusMaximum=concentration
    concentrationParameter=concentration
    potentialEnergyIntegral=Integrate(radiusMinimum,radiusMaximum,Potential_Energy_Integrand&
         &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    potentialEnergy=-0.5d0*(1.0d0/concentration+potentialEnergyIntegral)

    ! Compute the velocity dispersion at the virial radius.
    radiusMinimum=concentration
    radiusMaximum=100.0d0*concentration
    concentrationParameter=concentration
    jeansEquationIntegral=Integrate(radiusMinimum,radiusMaximum,Jeans_Equation_Integrand&
         &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)

    ! Compute the kinetic energy.
    radiusMinimum=0.0d0
    radiusMaximum=concentration
    concentrationParameter=concentration
    kineticEnergyIntegral=Integrate(radiusMinimum,radiusMaximum,Kinetic_Energy_Integrand&
         &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    kineticEnergy=2.0d0*Pi*(jeansEquationIntegral*concentration**3+kineticEnergyIntegral)
    
    ! Compute the total energy.
    NFW_Profile_Energy=(potentialEnergy+kineticEnergy)*concentration

    return
  end function NFW_Profile_Energy
  
  function Potential_Energy_Integrand(radius,parameterPointer) bind(c)
    !% Integrand for NFW profile potential energy.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double)          :: Potential_Energy_Integrand
    real(c_double), value   :: radius
    type(c_ptr),    value   :: parameterPointer
    
    Potential_Energy_Integrand=(Enclosed_Mass_NFW_Scale_Free(radius,concentrationParameter)/radius)**2
    return
  end function Potential_Energy_Integrand
  
  function Kinetic_Energy_Integrand(radius,parameterPointer) bind(c)
    !% Integrand for NFW profile kinetic energy.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double)          :: Kinetic_Energy_Integrand
    real(c_double), value   :: radius
    type(c_ptr),    value   :: parameterPointer

    Kinetic_Energy_Integrand=Enclosed_Mass_NFW_Scale_Free(radius,concentrationParameter)*Density_NFW_Scale_Free(radius&
         &,concentrationParameter)*radius
    return
  end function Kinetic_Energy_Integrand

  function Jeans_Equation_Integrand(radius,parameterPointer) bind(c)
    !% Integrand for NFW profile Jeans equation.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double)          :: Jeans_Equation_Integrand
    real(c_double), value   :: radius
    type(c_ptr),    value   :: parameterPointer
    
    Jeans_Equation_Integrand=Enclosed_Mass_NFW_Scale_Free(radius,concentrationParameter)*Density_NFW_Scale_Free(radius &
         &,concentrationParameter)/radius**2
    return
  end function Jeans_Equation_Integrand

  double precision function Dark_Matter_Profile_kSpace_NFW(thisNode,waveNumber)
    !% Returns the Fourier transform of the NFW density profile at the specified {\tt waveNumber} (given in Mpc$^{-1}$), using the
    !% expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Exponential_Integrals
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    double precision                     , intent(in   )          :: waveNumber
    class(nodeComponentDarkMatterProfile),                pointer :: thisDarkMatterProfileComponent
    double precision                                              :: radiusScale,waveNumberScaleFree,concentration
    
    ! Get components.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    ! Get the scale radius.
    radiusScale=thisDarkMatterProfileComponent%scale()

    ! Compute the concentration parameter.
    concentration=Dark_Matter_Halo_Virial_Radius(thisNode)/radiusScale

    ! Get the dimensionless wavenumber.
    waveNumberScaleFree=waveNumber*radiusScale

    ! Compute the Fourier transformed profile.
    Dark_Matter_Profile_kSpace_NFW=(                                                                                                                  &
         & +dsin(              waveNumberScaleFree)*(Sine_Integral  ((1.0d0+concentration)*waveNumberScaleFree)-Sine_Integral  (waveNumberScaleFree)) &
         & -dsin(concentration*waveNumberScaleFree)/(1.0d0+concentration)/waveNumberScaleFree                                                         &
         & +dcos(              waveNumberScaleFree)*(Cosine_Integral((1.0d0+concentration)*waveNumberScaleFree)-Cosine_Integral(waveNumberScaleFree)) &
         &                         )                                                                                                                  &
         & /(dlog(1.0d0+concentration)-concentration/(1.0d0+concentration))
    return
  end function Dark_Matter_Profile_kSpace_NFW

  double precision function Dark_Matter_Profile_Freefall_Radius_NFW(thisNode,time)
    !% Returns the freefall radius in the NFW density profile at the specified {\tt time} (given in Gyr).
    use Galacticus_Nodes
    use Numerical_Interpolation
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    double precision                     , intent(in   )          :: time
    class(nodeComponentDarkMatterProfile),                pointer :: thisDarkMatterProfileComponent
    double precision                                              :: freefallTimeScaleFree,radiusScale,concentration&
         &,velocityScale,timeScale

    ! For non-positive freefall times, return a zero freefall radius immediately.
    if (time <= 0.0d0) then
       Dark_Matter_Profile_Freefall_Radius_NFW=0.0d0
       return
    end if

    ! Get components.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    ! Get the scale radius.
    radiusScale=thisDarkMatterProfileComponent%scale()

    ! Get the concentration.
    concentration=Dark_Matter_Halo_Virial_Radius(thisNode)/radiusScale

    ! Get the virial velocity.
    velocityScale=Dark_Matter_Halo_Virial_Velocity(thisNode)

    ! Compute time scale.
    timeScale=Mpc_per_km_per_s_To_Gyr*radiusScale/velocityScale/dsqrt(concentration/(dlog(1.0d0+concentration)-concentration&
         &/(1.0d0+concentration)))

    ! Compute dimensionless time.
    freefallTimeScaleFree=time/timeScale

    ! Ensure table is sufficiently extensive.
    call Dark_Matter_Profile_NFW_Freefall_Tabulate(freefallTimeScaleFree)

    ! Interpolate to get the freefall radius.
    !$omp critical(NFW_Freefall_Interpolation)
    Dark_Matter_Profile_Freefall_Radius_NFW=Interpolate(nfwFreefallTableNumberPoints,nfwFreefallTime,nfwFreefallRadius &
         &,interpolationFreefallObject,interpolationFreefallAccelerator,freefallTimeScaleFree,reset=interpolationFreefallReset)&
         &*radiusScale
    !$omp end critical(NFW_Freefall_Interpolation)

    return
  end function Dark_Matter_Profile_Freefall_Radius_NFW
  
  double precision function Dark_Matter_Profile_Freefall_Radius_Increase_Rate_NFW(thisNode,time)
    !% Returns the rate of increase of the freefall radius in the NFW density profile at the specified {\tt time} (given in
    !% Gyr).
    use Galacticus_Nodes
    use Numerical_Interpolation
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    double precision                     , intent(in   )          :: time
    class(nodeComponentDarkMatterProfile),                pointer :: thisDarkMatterProfileComponent
    double precision                                              :: freefallTimeScaleFree,radiusScale,concentration&
         &,velocityScale,timeScale

    ! For non-positive freefall times, return the limiting value for small radii.
    if (time <= 0.0d0) then
       Dark_Matter_Profile_Freefall_Radius_Increase_Rate_NFW=0.0d0
       return
    end if

    ! Get components.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    ! Get the scale radius.
    radiusScale=thisDarkMatterProfileComponent%scale()

    ! Get the concentration.
    concentration=Dark_Matter_Halo_Virial_Radius(thisNode)/radiusScale

    ! Get the virial velocity.
    velocityScale=Dark_Matter_Halo_Virial_Velocity(thisNode)

    ! Compute time scale.
    timeScale=Mpc_per_km_per_s_To_Gyr*radiusScale/velocityScale/dsqrt(concentration/(dlog(1.0d0+concentration)-concentration&
         &/(1.0d0+concentration)))

    ! Compute dimensionless time.
    freefallTimeScaleFree=time/timeScale

    ! Ensure table is sufficiently extensive.
    call Dark_Matter_Profile_NFW_Freefall_Tabulate(freefallTimeScaleFree)

    ! Interpolate to get the freefall radius growth rate.
    !$omp critical(NFW_Freefall_Interpolation)
    Dark_Matter_Profile_Freefall_Radius_Increase_Rate_NFW=Interpolate_Derivative(nfwFreefallTableNumberPoints,nfwFreefallTime,nfwFreefallRadius &
         &,interpolationFreefallObject,interpolationFreefallAccelerator,freefallTimeScaleFree,reset=interpolationFreefallReset)&
         &*radiusScale/timeScale
    !$omp end critical(NFW_Freefall_Interpolation)
    return
  end function Dark_Matter_Profile_Freefall_Radius_Increase_Rate_NFW

  subroutine Dark_Matter_Profile_NFW_Freefall_Tabulate(freefallTimeScaleFree)
    !% Tabulates the freefall time vs. freefall radius for NFW halos.
    use Numerical_Ranges
    use Memory_Management
    use Numerical_Interpolation
    implicit none
    double precision, intent(in) :: freefallTimeScaleFree
    logical                      :: retabulate
    integer                      :: iRadius

    !$omp critical (NFW_Freefall_Interpolation)
    retabulate=.not.nfwFreefallTableInitialized
    ! If the table has not yet been made, compute and store the freefall corresponding to the minimum and maximum
    ! radii that will be tabulated by default.
    if (retabulate) then
       freefallTimeMinimum=Freefall_Time_Scale_Free(freefallRadiusMinimum)
       freefallTimeMaximum=Freefall_Time_Scale_Free(freefallRadiusMaximum)
    end if
    do while (freefallTimeScaleFree < freefallTimeMinimum)
       freefallRadiusMinimum=0.5d0*freefallRadiusMinimum
       freefallTimeMinimum=Freefall_Time_Scale_Free(freefallRadiusMinimum)
       retabulate=.true.
    end do
    do while (freefallTimeScaleFree > freefallTimeMaximum)
       freefallRadiusMaximum=2.0d0*freefallRadiusMaximum
       freefallTimeMaximum=Freefall_Time_Scale_Free(freefallRadiusMaximum)
       retabulate=.true.
    end do
    if (retabulate) then

       ! Decide how many points to tabulate and allocate table arrays.
       nfwFreefallTableNumberPoints=int(dlog10(freefallRadiusMaximum/freefallRadiusMinimum)*dble(nfwFreefallTablePointsPerDecade))+1
       if (allocated(nfwFreefallRadius)) then
          call Dealloc_Array(nfwFreefallRadius)
          call Dealloc_Array(nfwFreefallTime  )
       end if
       call Alloc_Array(nfwFreefallRadius,[nfwFreefallTableNumberPoints])
       call Alloc_Array(nfwFreefallTime  ,[nfwFreefallTableNumberPoints])
       ! Create a range of radii.
       nfwFreefallRadius=Make_Range(freefallRadiusMinimum,freefallRadiusMaximum,nfwFreefallTableNumberPoints,rangeType&
            &=rangeTypeLogarithmic)
       ! Loop over radii and populate tables.
       do iRadius=1,nfwFreefallTableNumberPoints
          nfwFreefallTime(iRadius)=Freefall_Time_Scale_Free(nfwFreefallRadius(iRadius))
       end do
       ! Ensure interpolations get reset.
       call Interpolate_Done(interpolationFreefallObject,interpolationFreefallAccelerator,interpolationFreefallReset)
       interpolationFreefallReset=.true.
       ! Specify that tabulation has been made.
       nfwFreefallTableInitialized=.true.
    end if
    !$omp end critical (NFW_Freefall_Interpolation)
    return
  end subroutine Dark_Matter_Profile_NFW_Freefall_Tabulate
  
  double precision function Freefall_Time_Scale_Free(radius)
    !% Compute the freefall time in a scale-free NFW halo.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    implicit none
    double precision,                intent(in) :: radius
    double precision,                parameter  :: radiusSmall=4.0d-6
    type(c_ptr)                                 :: parameterPointer
    type(fgsl_function)                         :: integrandFunction
    type(fgsl_integration_workspace)            :: integrationWorkspace
    double precision                            :: radiusEnd

    if (radius > radiusSmall) then
       ! Use the full solution.
       radiusStart=radius
       radiusEnd  =0.0d0
       Freefall_Time_Scale_Free=Integrate(radiusEnd,radiusStart,Freefall_Time_Scale_Free_Integrand,parameterPointer&
            &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
       call Integrate_Done(integrandFunction,integrationWorkspace)
    else
       ! Use an approximation here, found by taking series expansions of the logarithms in the integrand and keeping only the
       ! first order terms.
       Freefall_Time_Scale_Free=2.0d0*dsqrt(radius)
    end if
    return
  end function Freefall_Time_Scale_Free
  
  function Freefall_Time_Scale_Free_Integrand(radius,parameterPointer) bind(c)
    !% Integrand function used for finding the free-fall time in NFW halos.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double)            :: Freefall_Time_Scale_Free_Integrand
    real(c_double), value     :: radius
    type(c_ptr),    value     :: parameterPointer
    real(c_double), parameter :: radiusSmall        =1.0d-6
    real(c_double), parameter :: radiusSmallFraction=1.0d-3
    real(c_double)            :: x
    
     if (radius < radiusSmall) then
       ! Use a series approximation for small radii.
       Freefall_Time_Scale_Free_Integrand=dlog(1.0d0+radiusStart)/radiusStart-1.0d0+radius*(0.5d0-radius/3.0d0)
    else if (radius > radiusStart*(1.0d0-radiusSmallFraction)) then
       ! Use a series approximation for radii close to the initial radius.
       x=1.0d0-radius/radiusStart
       Freefall_Time_Scale_Free_Integrand=(1.0d0/(1.0d0+radiusStart)-dlog(1.0d0+radiusStart)/radiusStart)*x+(0.5d0*radiusStart&
            &/(1.0d0+radiusStart)**2+(radiusStart-(1.0d0+radiusStart)*dlog(1.0d0+radiusStart))/radiusStart/(1.0d0+radiusStart))*x&
            &**2
    else
       ! Use full expression for larger radii.
       Freefall_Time_Scale_Free_Integrand=dlog(1.0d0+radiusStart)/radiusStart-dlog(1.0d0+radius)/radius
    end if
    Freefall_Time_Scale_Free_Integrand=1.0d0/dsqrt(-2.0d0*Freefall_Time_Scale_Free_Integrand)
    return
  end function Freefall_Time_Scale_Free_Integrand
  
  !# <galacticusStateStoreTask>
  !#  <unitName>Dark_Matter_Profiles_NFW_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Dark_Matter_Profiles_NFW_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) concentrationMinimum,concentrationMaximum,radiusMinimum,radiusMaximum,freefallRadiusMinimum,freefallRadiusMaximum
    return
  end subroutine Dark_Matter_Profiles_NFW_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Dark_Matter_Profiles_NFW_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Dark_Matter_Profiles_NFW_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) concentrationMinimum,concentrationMaximum,radiusMinimum,radiusMaximum,freefallRadiusMinimum,freefallRadiusMaximum
    ! Retabulate.
    nfwTableInitialized        =.false.
    nfwInverseTableInitialized =.false.
    nfwFreefallTableInitialized=.false.
    call Dark_Matter_Profile_NFW_Tabulate
    call Dark_Matter_Profile_NFW_Inverse_Angular_Momentum
    return
  end subroutine Dark_Matter_Profiles_NFW_State_Retrieve
  
end module Dark_Matter_Profiles_NFW
