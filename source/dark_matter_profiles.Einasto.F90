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

!% Contains a module which implements Einasto halo profiles.

module Dark_Matter_Profiles_Einasto
  !% Implements Einasto halo profiles.
  use Galacticus_Nodes
  use FGSL
  implicit none
  private
  public :: Dark_Matter_Profile_Einasto_Initialize, Dark_Matter_Profiles_Einasto_State_Store, Dark_Matter_Profiles_Einasto_State_Retrieve

  ! Module scope variables used in integrations.
  double precision                                                              :: alphaParameter                                               , concentrationParameter                                         , &
       &                                                                           radiusStart                                                  , wavenumberParameter

  ! Tables for specific angular momentum vs. radius table
  double precision                                                              :: angularMomentumTableRadiusMinimum                    =1.0d-3
  double precision                                                              :: angularMomentumTableRadiusMaximum                    =20.0d+0
  integer                                                           , parameter :: angularMomentumTableRadiusPointsPerDecade            =100
  double precision                                                              :: angularMomentumTableAlphaMinimum                     =0.1d+0
  double precision                                                              :: angularMomentumTableAlphaMaximum                     =0.3d+0
  integer                                                           , parameter :: angularMomentumTableAlphaPointsPerUnit               =100
  logical                                                                       :: angularMomentumTableInitialized                      =.false.
  integer                                                                       :: angularMomentumTableAlphaCount                               , angularMomentumTableRadiusCount
  double precision                   , allocatable, dimension(:  )              :: angularMomentumTableAlpha                                    , angularMomentumTableRadius
  double precision                   , allocatable, dimension(:,:)              :: angularMomentumTable
  type            (fgsl_interp      )                                           :: angularMomentumTableRadiusInterpolationObject
  type            (fgsl_interp_accel)                                           :: angularMomentumTableAlphaInterpolationAccelerator            , angularMomentumTableRadiusInterpolationAccelerator
  logical                                                                       :: angularMomentumTableAlphaInterpolationReset          =.true. , angularMomentumTableRadiusInterpolationReset            =.true.

  ! Tables for freefall time vs. radius table
  double precision                                                              :: freefallRadiusTableRadiusMinimum                     =1.0d-3
  double precision                                                              :: freefallRadiusTableRadiusMaximum                     =20.0d+0
  integer                                                           , parameter :: freefallRadiusTableRadiusPointsPerDecade             =10
  double precision                                                              :: freefallRadiusTableAlphaMinimum                      =0.1d+0
  double precision                                                              :: freefallRadiusTableAlphaMaximum                      =0.3d+0
  integer                                                           , parameter :: freefallRadiusTableAlphaPointsPerUnit                =30
  logical                                                                       :: freefallRadiusTableInitialized                       =.false.
  integer                                                                       :: freefallRadiusTableAlphaCount                                , freefallRadiusTableRadiusCount
  double precision                   , allocatable, dimension(:  )              :: freefallRadiusTableAlpha                                     , freefallRadiusTableRadius
  double precision                   , allocatable, dimension(:,:)              :: freefallRadiusTable
  double precision                                                              :: freefallTimeMaximum                                          , freefallTimeMinimum
  type            (fgsl_interp      )                                           :: freefallRadiusTableRadiusInterpolationObject
  type            (fgsl_interp_accel)                                           :: freefallRadiusTableAlphaInterpolationAccelerator             , freefallRadiusTableRadiusInterpolationAccelerator
  logical                                                                       :: freefallRadiusTableAlphaInterpolationReset           =.true. , freefallRadiusTableRadiusInterpolationReset             =.true.

  ! Tables for energy as a function of concentration and alpha.
  double precision                                                              :: energyTableConcentrationMinimum                      =2.0d0
  double precision                                                              :: energyTableConcentrationMaximum                      =20.0d0
  integer                                                           , parameter :: energyTableConcentrationPointsPerDecade              =100
  double precision                                                              :: energyTableAlphaMinimum                              =0.1d0
  double precision                                                              :: energyTableAlphaMaximum                              =0.3d0
  integer                                                           , parameter :: energyTableAlphaPointsPerUnit                        =100
  logical                                                                       :: energyTableInitialized                               =.false.
  integer                                                                       :: energyTableAlphaCount                                        , energyTableConcentrationCount
  double precision                   , allocatable, dimension(:  )              :: energyTableAlpha                                             , energyTableConcentration
  double precision                   , allocatable, dimension(:,:)              :: energyTable
  type            (fgsl_interp      )                                           :: energyTableAlphaInterpolationObject                          , energyTableConcentrationInterpolationObject
  type            (fgsl_interp_accel)                                           :: energyTableAlphaInterpolationAccelerator                     , energyTableConcentrationInterpolationAccelerator
  logical                                                                       :: energyTableAlphaInterpolationReset                   =.true. , energyTableConcentrationInterpolationReset              =.true.

  ! Tables for specific Fourier transform of density profile as a function of alpha and radius.
  double precision                                                              :: fourierProfileTableConcentrationMinimum              =2.0d0
  double precision                                                              :: fourierProfileTableConcentrationMaximum              =20.0d0
  integer                                                           , parameter :: fourierProfileTableConcentrationPointsPerDecade      =10
  double precision                                                              :: fourierProfileTableWavenumberMinimum                 =1.0d-3
  double precision                                                              :: fourierProfileTableWavenumberMaximum                 =1.0d+3
  integer                                                           , parameter :: fourierProfileTableWavenumberPointsPerDecade         =10
  double precision                                                              :: fourierProfileTableAlphaMinimum                      =0.1d+0
  double precision                                                              :: fourierProfileTableAlphaMaximum                      =0.3d+0
  integer                                                           , parameter :: fourierProfileTableAlphaPointsPerUnit                =100
  logical                                                                       :: fourierProfileTableInitialized                       =.false.
  integer                                                                       :: fourierProfileTableAlphaCount                                , fourierProfileTableConcentrationCount                          , &
       &                                                                           fourierProfileTableWavenumberCount
  double precision                   , allocatable, dimension(:    )            :: fourierProfileTableAlpha                                     , fourierProfileTableConcentration                               , &
       &                                                                           fourierProfileTableWavenumber
  double precision                   , allocatable, dimension(:,:,:)            :: fourierProfileTable
  type            (fgsl_interp      )                                           :: fourierProfileTableWavenumberInterpolationObject
  type            (fgsl_interp_accel)                                           :: fourierProfileTableAlphaInterpolationAccelerator             , fourierProfileTableConcentrationInterpolationAccelerator       , &
       &                                                                           fourierProfileTableWavenumberInterpolationAccelerator
  logical                                                                       :: fourierProfileTableAlphaInterpolationReset           =.true. , fourierProfileTableConcentrationInterpolationReset      =.true., &
       &                                                                           fourierProfileTableWavenumberInterpolationReset      =.true.

contains

  !# <darkMatterProfileMethod>
  !#  <unitName>Dark_Matter_Profile_Einasto_Initialize</unitName>
  !# </darkMatterProfileMethod>
  subroutine Dark_Matter_Profile_Einasto_Initialize(darkMatterProfileMethod,Dark_Matter_Profile_Density_Get&
       &,Dark_Matter_Profile_Energy_Get ,Dark_Matter_Profile_Energy_Growth_Rate_Get &
       &,Dark_Matter_Profile_Rotation_Normalization_Get ,Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get &
       &,Dark_Matter_Profile_Circular_Velocity_Get ,Dark_Matter_Profile_Potential_Get,Dark_Matter_Profile_Enclosed_Mass_Get &
       &,Dark_Matter_Profile_kSpace_Get,Dark_Matter_Profile_Freefall_Radius_Get &
       &,Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Get)
    !% Initializes the ``Einasto'' halo profile module.
    use ISO_Varying_String
    use Galacticus_Error
    use Array_Utilities
    implicit none
    type     (varying_string                                           ), intent(in   )          :: darkMatterProfileMethod
    procedure(Dark_Matter_Profile_Density_Einasto                      ), intent(inout), pointer :: Dark_Matter_Profile_Density_Get
    procedure(Dark_Matter_Profile_Energy_Einasto                       ), intent(inout), pointer :: Dark_Matter_Profile_Energy_Get
    procedure(Dark_Matter_Profile_Energy_Growth_Rate_Einasto           ), intent(inout), pointer :: Dark_Matter_Profile_Energy_Growth_Rate_Get
    procedure(Dark_Matter_Profile_Rotation_Normalization_Einasto       ), intent(inout), pointer :: Dark_Matter_Profile_Rotation_Normalization_Get
    procedure(Radius_from_Specific_Angular_Momentum_Einasto            ), intent(inout), pointer :: Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get
    procedure(Dark_Matter_Profile_Circular_Velocity_Einasto            ), intent(inout), pointer :: Dark_Matter_Profile_Circular_Velocity_Get
    procedure(Dark_Matter_Profile_Potential_Einasto                    ), intent(inout), pointer :: Dark_Matter_Profile_Potential_Get
    procedure(Dark_Matter_Profile_Enclosed_Mass_Einasto                ), intent(inout), pointer :: Dark_Matter_Profile_Enclosed_Mass_Get
    procedure(Dark_Matter_Profile_kSpace_Einasto                       ), intent(inout), pointer :: Dark_Matter_Profile_kSpace_Get
    procedure(Dark_Matter_Profile_Freefall_Radius_Einasto              ), intent(inout), pointer :: Dark_Matter_Profile_Freefall_Radius_Get
    procedure(Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Einasto), intent(inout), pointer :: Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Get

    if (darkMatterProfileMethod == 'Einasto') then
       Dark_Matter_Profile_Density_Get                               => Dark_Matter_Profile_Density_Einasto
       Dark_Matter_Profile_Energy_Get                                => Dark_Matter_Profile_Energy_Einasto
       Dark_Matter_Profile_Energy_Growth_Rate_Get                    => Dark_Matter_Profile_Energy_Growth_Rate_Einasto
       Dark_Matter_Profile_Rotation_Normalization_Get                => Dark_Matter_Profile_Rotation_Normalization_Einasto
       Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get => Radius_from_Specific_Angular_Momentum_Einasto
       Dark_Matter_Profile_Circular_Velocity_Get                     => Dark_Matter_Profile_Circular_Velocity_Einasto
       Dark_Matter_Profile_Potential_Get                             => Dark_Matter_Profile_Potential_Einasto
       Dark_Matter_Profile_Enclosed_Mass_Get                         => Dark_Matter_Profile_Enclosed_Mass_Einasto
       Dark_Matter_Profile_kSpace_Get                                => Dark_Matter_Profile_kSpace_Einasto
       Dark_Matter_Profile_Freefall_Radius_Get                       => Dark_Matter_Profile_Freefall_Radius_Einasto
       Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Get         => Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Einasto
       ! Ensure that the dark matter profile component supports both "scale" and "shape" properties. Since we've been called with
       ! a treeNode to process, it should have been initialized by now.
       if     (                                                                                                                 &
            &  .not.(                                                                                                           &
            &         defaultDarkMatterProfileComponent%scaleIsGettable()                                                       &
            &        .and.                                                                                                      &
            &         defaultDarkMatterProfileComponent%shapeIsGettable()                                                       &
            &       )                                                                                                           &
            & ) call Galacticus_Error_Report                                                                                    &
            &        (                                                                                                          &
            &         'Dark_Matter_Profile_Einasto_Initialize'                                                                , &
            &         'Einasto dark matter profile requires a dark matter profile component that supports gettable '//          &
            &         '"scale" and "shape" properties.'                                                             //          &
            &         Galacticus_Component_List(                                                                                &
            &                                   'darkMatterProfile'                                                           , &
            &                                    defaultDarkMatterProfileComponent%shapeAttributeMatch(requireGettable=.true.)  &
            &                                   .intersection.                                                                  &
            &                                    defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)  &
            &                                  )                                                                                &
            &        )
    end if
    return
  end subroutine Dark_Matter_Profile_Einasto_Initialize

  double precision function Dark_Matter_Profile_Density_Einasto(thisNode,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given
    !% in units of Mpc).
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: thisNode
    double precision                                , intent(in   )          :: radius
    class           (nodeComponentBasic            )               , pointer :: thisBasicComponent
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: alpha                         , radiusOverScaleRadius      , &
         &                                                                      scaleRadius                   , virialRadiusOverScaleRadius

    ! Get components.
    thisBasicComponent             => thisNode%basic            (                 )
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    scaleRadius                =thisDarkMatterProfileComponent%scale()
    alpha                      =thisDarkMatterProfileComponent%shape()
    radiusOverScaleRadius      =radius                                  /scaleRadius
    virialRadiusOverScaleRadius=Dark_Matter_Halo_Virial_Radius(thisNode)/scaleRadius
    Dark_Matter_Profile_Density_Einasto=Density_Einasto_Scale_Free(radiusOverScaleRadius,virialRadiusOverScaleRadius,alpha)&
         &*thisBasicComponent%mass()/scaleRadius**3
    return
  end function Dark_Matter_Profile_Density_Einasto

  double precision function Dark_Matter_Profile_Enclosed_Mass_Einasto(thisNode,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc).
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: thisNode
    double precision                                , intent(in   )          :: radius
    class           (nodeComponentBasic            )               , pointer :: thisBasicComponent
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: alpha                         , radiusOverScaleRadius      , &
         &                                                                      scaleRadius                   , virialRadiusOverScaleRadius

    ! Get components.
    thisBasicComponent             => thisNode%basic            (                 )
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    scaleRadius                =thisDarkMatterProfileComponent%scale()
    alpha                      =thisDarkMatterProfileComponent%shape()
    radiusOverScaleRadius      =radius                                  /scaleRadius
    virialRadiusOverScaleRadius=Dark_Matter_Halo_Virial_Radius(thisNode)/scaleRadius
    Dark_Matter_Profile_Enclosed_Mass_Einasto=Enclosed_Mass_Einasto_Scale_Free(radiusOverScaleRadius,virialRadiusOverScaleRadius,alpha)&
         &*thisBasicComponent%mass()
    return
  end function Dark_Matter_Profile_Enclosed_Mass_Einasto

  double precision function Dark_Matter_Profile_Circular_Velocity_Einasto(thisNode,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc).
    use Numerical_Constants_Physical
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: radius

    if (radius > 0.0d0) then
       Dark_Matter_Profile_Circular_Velocity_Einasto=sqrt(gravitationalConstantGalacticus&
            &*Dark_Matter_Profile_Enclosed_Mass_Einasto(thisNode,radius)/radius)
    else
       Dark_Matter_Profile_Circular_Velocity_Einasto=0.0d0
    end if
    return
  end function Dark_Matter_Profile_Circular_Velocity_Einasto

  double precision function Dark_Matter_Profile_Potential_Einasto(thisNode,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc).
    use Dark_Matter_Halo_Scales
    use Dark_Matter_Profiles_Error_Codes
    use Numerical_Constants_Physical
    implicit none
    type            (treeNode                      ), intent(inout), pointer  :: thisNode
    double precision                                , intent(in   )           :: radius
    integer                                         , intent(  out), optional :: status
    class           (nodeComponentBasic            )               , pointer  :: thisBasicComponent
    class           (nodeComponentDarkMatterProfile)               , pointer  :: thisDarkMatterProfileComponent
    double precision                                                          :: alpha                         , radiusOverScaleRadius      , &
         &                                                                       scaleRadius                   , virialRadiusOverScaleRadius

    ! Assume success.
    if (present(status)) status=darkMatterProfileSuccess
    ! Get components.
    thisBasicComponent             => thisNode%basic            (                 )
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    scaleRadius                =thisDarkMatterProfileComponent%scale()
    alpha                      =thisDarkMatterProfileComponent%shape()
    radiusOverScaleRadius      =radius                                  /scaleRadius
    virialRadiusOverScaleRadius=Dark_Matter_Halo_Virial_Radius(thisNode)/scaleRadius
    Dark_Matter_Profile_Potential_Einasto=                                                                   &
         & ( Potential_Einasto_Scale_Free(radiusOverScaleRadius      ,virialRadiusOverScaleRadius,alpha)     &
         &  -Potential_Einasto_Scale_Free(virialRadiusOverScaleRadius,virialRadiusOverScaleRadius,alpha)     &
         &  -1.0d0/                                                   virialRadiusOverScaleRadius            &
         & )                                                                                                 &
         & *gravitationalConstantGalacticus*thisBasicComponent%mass()/scaleRadius
    return
  end function Dark_Matter_Profile_Potential_Einasto

  double precision function Radius_from_Specific_Angular_Momentum_Einasto(thisNode,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\tt thisNode} at which a circular orbit has the given {\tt specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc).
    use Numerical_Constants_Physical
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: thisNode
    double precision                                , intent(in   )          :: specificAngularMomentum
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: alpha                           , scaleRadius, &
         &                                                                      specificAngularMomentumScaleFree

    ! Get components.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    ! Get the scale radius of the halo.
    scaleRadius=thisDarkMatterProfileComponent%scale()
    ! Get the shape parameter of the halo.
    alpha      =thisDarkMatterProfileComponent%shape()
    ! Compute the specific angular momentum in scale free units.
    specificAngularMomentumScaleFree=specificAngularMomentum/sqrt(gravitationalConstantGalacticus*scaleRadius&
         &*Dark_Matter_Profile_Enclosed_Mass_Einasto(thisNode,scaleRadius))
    ! Compute the corresponding radius.
    Radius_from_Specific_Angular_Momentum_Einasto=scaleRadius*Radius_from_Specific_Angular_Momentum_Scale_Free(alpha&
         &,specificAngularMomentumScaleFree)
    return
  end function Radius_from_Specific_Angular_Momentum_Einasto

  double precision function Radius_from_Specific_Angular_Momentum_Scale_Free(alpha,specificAngularMomentumScaleFree)
    !% Comptue the radius at which a circular orbit has the given {\tt specificAngularMomentumScaleFree} in a scale free Einasto
    !% profile.
    use Numerical_Interpolation
    implicit none
    double precision, intent(in   )  :: alpha , specificAngularMomentumScaleFree
    integer         , dimension(0:1) :: jAlpha
    double precision, dimension(0:1) :: hAlpha
    integer                          :: iAlpha

    ! Return immediately for zero angular momentum.
    if (specificAngularMomentumScaleFree <= 0.0d0) then
       Radius_from_Specific_Angular_Momentum_Scale_Free=0.0d0
       return
    end if

    ! Ensure the table exists and is sufficiently tabulated.
    call Radius_from_Specific_Angular_Momentum_Table_Make(alpha,specificAngularMomentumScaleFree)

    ! Get interpolating factors in alpha.
    !$omp critical (Einasto_Interpolate_Specific_Angular_Momentum)
    jAlpha(0)=Interpolate_Locate(angularMomentumTableAlphaCount,angularMomentumTableAlpha&
         &,angularMomentumTableAlphaInterpolationAccelerator,alpha,reset=angularMomentumTableAlphaInterpolationReset)
    !$omp end critical (Einasto_Interpolate_Specific_Angular_Momentum)
    jAlpha(1)=jAlpha(0)+1
    hAlpha=Interpolate_Linear_Generate_Factors(angularMomentumTableAlphaCount,angularMomentumTableAlpha,jAlpha(0),alpha)

    ! Interpolate in specific angular momentum to get radius.
    Radius_from_Specific_Angular_Momentum_Scale_Free=0.0d0
    do iAlpha=0,1
       Radius_from_Specific_Angular_Momentum_Scale_Free=                       &
            &  Radius_from_Specific_Angular_Momentum_Scale_Free                &
            & +Interpolate( angularMomentumTableRadiusCount                    &
            &              ,angularMomentumTable(:,jAlpha(iAlpha))             &
            &              ,angularMomentumTableRadius                         &
            &              ,angularMomentumTableRadiusInterpolationObject      &
            &              ,angularMomentumTableRadiusInterpolationAccelerator &
            &              ,specificAngularMomentumScaleFree                   &
            &              ,reset=angularMomentumTableRadiusInterpolationReset &
            &             )                                                    &
            & *hAlpha(iAlpha)
    end do

    return
  end function Radius_from_Specific_Angular_Momentum_Scale_Free

  subroutine Radius_from_Specific_Angular_Momentum_Table_Make(alphaRequired,specificAngularMomentumRequired)
    !% Create a tabulation of the relation between specific angular momentum and radius in an Einasto profile.
    use Numerical_Interpolation
    use Numerical_Ranges
    use Gamma_Functions
    use Memory_Management
    implicit none
    double precision, intent(in   ) :: alphaRequired, specificAngularMomentumRequired
    integer                         :: iAlpha       , iRadius
    logical                         :: makeTable
    double precision                :: alpha        , enclosedMass                   , &
         &                             radius

    !$omp critical (Einasto_Interpolate_Specific_Angular_Momentum)
    ! Always check if we need to make the table.
    makeTable=.true.
    do while (makeTable)
       ! Assume table does not need remaking.
       makeTable=.false.
       ! Check for uninitialized table.
       if (.not.angularMomentumTableInitialized) then
          makeTable=.true.
          ! Check for alpha out of range.
       else if (alphaRequired < angularMomentumTableAlpha(1) .or. alphaRequired >&
            & angularMomentumTableAlpha(angularMomentumTableAlphaCount)) then
          makeTable=.true.
          ! Compute the range of tabulation and number of points to use.
          angularMomentumTableAlphaMinimum =min(angularMomentumTableAlphaMinimum,0.9d0*alphaRequired)
          angularMomentumTableAlphaMaximum =max(angularMomentumTableAlphaMaximum,1.1d0*alphaRequired)
          ! Check for angular momentum below minimum tabulated value.
       else if (any(specificAngularMomentumRequired < angularMomentumTable(1,:))) then
          makeTable=.true.
          angularMomentumTableRadiusMinimum=0.5d0*angularMomentumTableRadiusMinimum
          ! Check for angular momentum above maximum tabulated value.
       else if (any(specificAngularMomentumRequired > angularMomentumTable(angularMomentumTableRadiusCount,:))) then
          makeTable=.true.
          angularMomentumTableRadiusMaximum=2.0d0*angularMomentumTableRadiusMaximum
       end if
       ! Remake the table if necessary.
       if (makeTable) then
          ! Allocate arrays to the appropriate sizes.
          angularMomentumTableAlphaCount =int(      (angularMomentumTableAlphaMaximum -angularMomentumTableAlphaMinimum )&
               &*dble(angularMomentumTableAlphaPointsPerUnit   ))+1
          angularMomentumTableRadiusCount=int(log10(angularMomentumTableRadiusMaximum/angularMomentumTableRadiusMinimum)&
               &*dble(angularMomentumTableRadiusPointsPerDecade))+1
          if (allocated(angularMomentumTableAlpha )) call Dealloc_Array(angularMomentumTableAlpha )
          if (allocated(angularMomentumTableRadius)) call Dealloc_Array(angularMomentumTableRadius)
          if (allocated(angularMomentumTable      )) call Dealloc_Array(angularMomentumTable      )
          call Alloc_Array(angularMomentumTableAlpha ,[                                angularMomentumTableAlphaCount])
          call Alloc_Array(angularMomentumTableRadius,[angularMomentumTableRadiusCount                               ])
          call Alloc_Array(angularMomentumTable      ,[angularMomentumTableRadiusCount,angularMomentumTableAlphaCount])
          ! Create ranges of alpha and radius.
          angularMomentumTableAlpha =Make_Range(angularMomentumTableAlphaMinimum ,angularMomentumTableAlphaMaximum &
               &,angularMomentumTableAlphaCount ,rangeType=rangeTypeLinear     )
          angularMomentumTableRadius=Make_Range(angularMomentumTableRadiusMinimum,angularMomentumTableRadiusMaximum&
               &,angularMomentumTableRadiusCount,rangeType=rangeTypeLogarithmic)
          ! Tabulate the radius vs. specific angular momentum relation.
          do iAlpha=1,angularMomentumTableAlphaCount
             alpha=angularMomentumTableAlpha(iAlpha)
             do iRadius=1,angularMomentumTableRadiusCount
                radius=angularMomentumTableRadius(iRadius)
                enclosedMass= Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*radius**alpha/alpha) &
                     &       /Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0              /alpha)
                angularMomentumTable(iRadius,iAlpha)=sqrt(enclosedMass*radius)
             end do
          end do
          ! Reset interpolators.
          call Interpolate_Done(angularMomentumTableRadiusInterpolationObject,angularMomentumTableRadiusInterpolationAccelerator &
               &,angularMomentumTableRadiusInterpolationReset)
          call Interpolate_Done(interpolationAccelerator=angularMomentumTableAlphaInterpolationAccelerator &
               &,reset=angularMomentumTableAlphaInterpolationReset)
          angularMomentumTableRadiusInterpolationReset=.true.
          angularMomentumTableAlphaInterpolationReset =.true.
          ! Flag that the table is now initialized.
          angularMomentumTableInitialized=.true.
       end if
    end do
    !$omp end critical (Einasto_Interpolate_Specific_Angular_Momentum)
    return
  end subroutine Radius_from_Specific_Angular_Momentum_Table_Make

  double precision function Dark_Matter_Profile_Rotation_Normalization_Einasto(thisNode)
    !% Return the rotation normalization of an Einasto halo density profile.
    use Numerical_Constants_Math
    use Gamma_Functions
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: thisNode
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: alpha                         , scaleRadius, &
         &                                                                      virialRadiusOverScaleRadius

    ! Get components.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    ! Get scale radius, shape and concentration.
    scaleRadius                =thisDarkMatterProfileComponent%scale()
    alpha                      =thisDarkMatterProfileComponent%shape()
    virialRadiusOverScaleRadius=Dark_Matter_Halo_Virial_Radius(thisNode)/scaleRadius

    ! Compute the rotation normalization.
    Dark_Matter_Profile_Rotation_Normalization_Einasto=                                                         &
         &  (4.0d0/Pi/scaleRadius)                                                                              &
         & *(2.0d0/alpha)**(1.0d0/alpha)                                                                        &
         & *Gamma_Function                         (3.0d0/alpha                                               ) &
         & /Gamma_Function                         (4.0d0/alpha                                               ) &
         & *Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*virialRadiusOverScaleRadius**alpha/alpha) &
         & /Gamma_Function_Incomplete_Complementary(4.0d0/alpha,2.0d0*virialRadiusOverScaleRadius**alpha/alpha)

    return
  end function Dark_Matter_Profile_Rotation_Normalization_Einasto

  double precision function Dark_Matter_Profile_Energy_Einasto(thisNode)
    !% Return the energy of an Einasto halo density profile.
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    implicit none
    type            (treeNode                      ), intent(inout) , pointer :: thisNode
    class           (nodeComponentBasic            )                , pointer :: thisBasicComponent
    class           (nodeComponentDarkMatterProfile)                , pointer :: thisDarkMatterProfileComponent
    integer                                         , dimension(0:1)          :: jAlpha
    double precision                                , dimension(0:1)          :: hAlpha
    integer                                                                   :: iAlpha
    double precision                                                          :: alpha                         , scaleRadius, &
         &                                                                       virialRadiusOverScaleRadius

    ! Get components.
    thisBasicComponent             => thisNode%basic            (                 )
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    ! Get scale radius, shape parameter and concentration.
    scaleRadius                =thisDarkMatterProfileComponent%scale()
    alpha                      =thisDarkMatterProfileComponent%shape()
    virialRadiusOverScaleRadius=Dark_Matter_Halo_Virial_Radius(thisNode)/scaleRadius

    ! Ensure the table exists and is sufficiently tabulated.
    call Energy_Table_Make(virialRadiusOverScaleRadius,alpha)

    !$omp critical(Einasto_Interpolation)
    ! Get interpolating factors in alpha.
    jAlpha(0)=Interpolate_Locate(energyTableAlphaCount,energyTableAlpha,energyTableAlphaInterpolationAccelerator,alpha,reset&
         &=energyTableAlphaInterpolationReset)
    jAlpha(1)=jAlpha(0)+1
    hAlpha=Interpolate_Linear_Generate_Factors(energyTableAlphaCount,energyTableAlpha,jAlpha(0),alpha)

    ! Find the energy by interpolation.
    Dark_Matter_Profile_Energy_Einasto=0.0d0
    do iAlpha=0,1
       Dark_Matter_Profile_Energy_Einasto=Dark_Matter_Profile_Energy_Einasto+Interpolate(energyTableConcentrationCount&
            &,energyTableConcentration,energyTable(:,jAlpha(iAlpha)),energyTableConcentrationInterpolationObject &
            &,energyTableConcentrationInterpolationAccelerator,virialRadiusOverScaleRadius,reset&
            &=energyTableConcentrationInterpolationReset)*hAlpha(iAlpha)
    end do

    ! Scale to dimensionful units.
    Dark_Matter_Profile_Energy_Einasto=Dark_Matter_Profile_Energy_Einasto*thisBasicComponent%mass() &
         &*Dark_Matter_Halo_Virial_Velocity(thisNode)**2
    !$omp end critical(Einasto_Interpolation)
    return
  end function Dark_Matter_Profile_Energy_Einasto

  double precision function Dark_Matter_Profile_Energy_Growth_Rate_Einasto(thisNode)
    !% Return the energy of an Einasto halo density profile.
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    implicit none
    type            (treeNode                      ), intent(inout) , pointer :: thisNode
    class           (nodeComponentBasic            )                , pointer :: thisBasicComponent
    class           (nodeComponentDarkMatterProfile)                , pointer :: thisDarkMatterProfileComponent
    integer                                         , dimension(0:1)          :: jAlpha
    double precision                                , dimension(0:1)          :: hAlpha
    integer                                                                   :: iAlpha
    double precision                                                          :: alpha                         , energy     , &
         &                                                                       energyGradient                , scaleRadius, &
         &                                                                       virialRadiusOverScaleRadius

    ! Get components.
    thisBasicComponent             => thisNode%basic            (                 )
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    ! Get scale radius, shape parameter and concentration.
    scaleRadius                =thisDarkMatterProfileComponent%scale()
    alpha                      =thisDarkMatterProfileComponent%shape()
    virialRadiusOverScaleRadius=Dark_Matter_Halo_Virial_Radius(thisNode)/scaleRadius

    ! Ensure the table exists and is sufficiently tabulated.
    call Energy_Table_Make(virialRadiusOverScaleRadius,alpha)

    !$omp critical(Einasto_Interpolation)
    ! Get interpolating factors in alpha.
    jAlpha(0)=Interpolate_Locate(energyTableAlphaCount,energyTableAlpha,energyTableAlphaInterpolationAccelerator,alpha,reset&
         &=energyTableAlphaInterpolationReset)
    jAlpha(1)=jAlpha(0)+1
    hAlpha=Interpolate_Linear_Generate_Factors(energyTableAlphaCount,energyTableAlpha,jAlpha(0),alpha)

    ! Find the energy gradient by interpolation.
    energy        =0.0d0
    energyGradient=0.0d0
    do iAlpha=0,1
       energy        =energy        +Interpolate           (energyTableConcentrationCount,energyTableConcentration,energyTable(:&
            &,jAlpha(iAlpha)),energyTableConcentrationInterpolationObject ,energyTableConcentrationInterpolationAccelerator&
            &,virialRadiusOverScaleRadius,reset=energyTableConcentrationInterpolationReset)*hAlpha(iAlpha)
       energyGradient=energyGradient+Interpolate_Derivative(energyTableConcentrationCount,energyTableConcentration,energyTable(:&
            &,jAlpha(iAlpha)),energyTableConcentrationInterpolationObject ,energyTableConcentrationInterpolationAccelerator&
            &,virialRadiusOverScaleRadius,reset=energyTableConcentrationInterpolationReset)*hAlpha(iAlpha)
    end do
    !$omp end critical(Einasto_Interpolation)

    ! Compute the energy growth rate.
    Dark_Matter_Profile_Energy_Growth_Rate_Einasto=Dark_Matter_Profile_Energy_Einasto(thisNode)&
         &*(thisBasicComponent%accretionRate()/thisBasicComponent%mass()+2.0d0 &
         &*Dark_Matter_Halo_Virial_Velocity_Growth_Rate(thisNode)/Dark_Matter_Halo_Virial_Velocity(thisNode)+(energyGradient&
         &*virialRadiusOverScaleRadius/energy)*(Dark_Matter_Halo_Virial_Radius_Growth_Rate(thisNode)&
         &/Dark_Matter_Halo_Virial_Radius(thisNode)-thisDarkMatterProfileComponent%scaleGrowthRate()&
         &/thisDarkMatterProfileComponent%scale()))

    return
  end function Dark_Matter_Profile_Energy_Growth_Rate_Einasto

  subroutine Energy_Table_Make(concentrationRequired,alphaRequired)
    !% Create a tabulation of the energy of Einasto profiles as a function of their concentration of $\alpha$ parameter.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Interpolation
    use Numerical_Integration
    use Numerical_Ranges
    use Memory_Management
    use Numerical_Constants_Math
    implicit none
    double precision                            , intent(in   ) :: alphaRequired          , concentrationRequired
    integer                                                     :: iAlpha                 , iConcentration
    logical                                                     :: makeTable
    double precision                                            :: alpha                  , concentration        , &
         &                                                         jeansEquationIntegral  , kineticEnergy        , &
         &                                                         kineticEnergyIntegral  , potentialEnergy      , &
         &                                                         potentialEnergyIntegral, radiusMaximum        , &
         &                                                         radiusMinimum
    type            (c_ptr                     )                :: parameterPointer
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace

    !$omp critical (Einasto_Interpolation)
    ! Assume table does not need remaking.
    makeTable=.false.
    ! Check for uninitialized table.
    if (.not.energyTableInitialized) makeTable=.true.
    ! Check for alpha out of range.
    if (alphaRequired < energyTableAlphaMinimum .or. alphaRequired > energyTableAlphaMaximum) then
       makeTable=.true.
       ! Compute the range of tabulation and number of points to use.
       energyTableAlphaMinimum =min(energyTableAlphaMinimum,0.9d0*alphaRequired)
       energyTableAlphaMaximum =max(energyTableAlphaMaximum,1.1d0*alphaRequired)
    end if
    ! Check for concentration below minimum tabulated value.
    if (concentrationRequired < energyTableConcentrationMinimum .or. concentrationRequired > energyTableConcentrationMaximum) then
       makeTable=.true.
       energyTableConcentrationMinimum=min(energyTableConcentrationMinimum,0.5d0*energyTableConcentrationMinimum)
       energyTableConcentrationMaximum=max(energyTableConcentrationMaximum,2.0d0*energyTableConcentrationMaximum)
    end if
    ! Remake the table if necessary.
    if (makeTable) then
       ! Allocate arrays to the appropriate sizes.
       energyTableAlphaCount        =int(      (energyTableAlphaMaximum        -energyTableAlphaMinimum        ) &
            &*dble(energyTableAlphaPointsPerUnit          ))+1
       energyTableConcentrationCount=int(log10(energyTableConcentrationMaximum/energyTableConcentrationMinimum) &
            &*dble(energyTableConcentrationPointsPerDecade))+1
       if (allocated(energyTableAlpha        )) call Dealloc_Array(energyTableAlpha        )
       if (allocated(energyTableConcentration)) call Dealloc_Array(energyTableConcentration)
       if (allocated(energyTable             )) call Dealloc_Array(energyTable             )
       call Alloc_Array(energyTableAlpha        ,[                              energyTableAlphaCount])
       call Alloc_Array(energyTableConcentration,[energyTableConcentrationCount                      ])
       call Alloc_Array(energyTable             ,[energyTableConcentrationCount,energyTableAlphaCount])
       ! Create ranges of alpha and concentration.
       energyTableAlpha        =Make_Range(energyTableAlphaMinimum        ,energyTableAlphaMaximum         &
            &,energyTableAlphaCount        ,rangeType=rangeTypeLinear     )
       energyTableConcentration=Make_Range(energyTableConcentrationMinimum,energyTableConcentrationMaximum &
            &,energyTableConcentrationCount,rangeType=rangeTypeLogarithmic)
       ! Tabulate the radius vs. specific angular momentum relation.
       do iAlpha=1,energyTableAlphaCount
          alpha=energyTableAlpha(iAlpha)
          do iConcentration=1,energyTableConcentrationCount
             concentration=energyTableConcentration(iConcentration)

             ! Compute the potential energy.
             radiusMinimum         =0.0d0
             radiusMaximum         =concentration
             concentrationParameter=concentration
             alphaParameter        =alpha
             potentialEnergyIntegral=Integrate(radiusMinimum,radiusMaximum,Potential_Energy_Integrand_Einasto&
                  &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
             call Integrate_Done(integrandFunction,integrationWorkspace)
             potentialEnergy=-0.5d0*(1.0d0/concentration+potentialEnergyIntegral)

             ! Compute the velocity dispersion at the virial radius.
             radiusMinimum         =        concentration
             radiusMaximum         =100.0d0*concentration
             concentrationParameter=        concentration
             alphaParameter        =alpha
             jeansEquationIntegral=Integrate(radiusMinimum,radiusMaximum,Jeans_Equation_Integrand_Einasto&
                  &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
             call Integrate_Done(integrandFunction,integrationWorkspace)

             ! Compute the kinetic energy.
             radiusMinimum         =0.0d0
             radiusMaximum         =concentration
             concentrationParameter=concentration
             alphaParameter        =alpha
             kineticEnergyIntegral=Integrate(radiusMinimum,radiusMaximum,Kinetic_Energy_Integrand_Einasto&
                  &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
             call Integrate_Done(integrandFunction,integrationWorkspace)
             kineticEnergy=2.0d0*Pi*(jeansEquationIntegral*concentration**3+kineticEnergyIntegral)

             ! Compute the total energy.
             energyTable(iConcentration,iAlpha)=(potentialEnergy+kineticEnergy)*concentration

          end do
       end do
       ! Reset interpolators.
       call Interpolate_Done(energyTableConcentrationInterpolationObject,energyTableConcentrationInterpolationAccelerator &
            &,energyTableConcentrationInterpolationReset)
       call Interpolate_Done(energyTableAlphaInterpolationObject        ,energyTableAlphaInterpolationAccelerator         &
            &,energyTableAlphaInterpolationReset        )
       energyTableConcentrationInterpolationReset=.true.
       energyTableAlphaInterpolationReset        =.true.
       ! Flag that the table is now initialized.
       energyTableInitialized=.true.
    end if
    !$omp end critical (Einasto_Interpolation)
    return
  end subroutine Energy_Table_Make

  function Potential_Energy_Integrand_Einasto(radius,parameterPointer) bind(c)
    !% Integrand for Einasto profile potential energy.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(kind=c_double)        :: Potential_Energy_Integrand_Einasto
    real(kind=c_double), value :: radius
    type(c_ptr        ), value :: parameterPointer

    Potential_Energy_Integrand_Einasto=(Enclosed_Mass_Einasto_Scale_Free(radius,concentrationParameter,alphaParameter)/radius)**2
    return
  end function Potential_Energy_Integrand_Einasto

  function Kinetic_Energy_Integrand_Einasto(radius,parameterPointer) bind(c)
    !% Integrand for Einasto profile kinetic energy.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(kind=c_double)        :: Kinetic_Energy_Integrand_Einasto
    real(kind=c_double), value :: radius
    type(c_ptr        ), value :: parameterPointer

    Kinetic_Energy_Integrand_Einasto=Enclosed_Mass_Einasto_Scale_Free(radius,concentrationParameter,alphaParameter)&
         &*Density_Einasto_Scale_Free(radius,concentrationParameter,alphaParameter)*radius
    return
  end function Kinetic_Energy_Integrand_Einasto

  function Jeans_Equation_Integrand_Einasto(radius,parameterPointer) bind(c)
    !% Integrand for Einasto profile Jeans equation.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(kind=c_double)        :: Jeans_Equation_Integrand_Einasto
    real(kind=c_double), value :: radius
    type(c_ptr        ), value :: parameterPointer

    Jeans_Equation_Integrand_Einasto=Enclosed_Mass_Einasto_Scale_Free(radius,concentrationParameter,alphaParameter)&
         &*Density_Einasto_Scale_Free(radius ,concentrationParameter,alphaParameter)/radius**2
    return
  end function Jeans_Equation_Integrand_Einasto

  double precision function Enclosed_Mass_Einasto_Scale_Free(radius,concentration,alpha)
    !% Returns the enclosed mass (in units of the virial mass) in an Einasto dark matter profile with given {\tt concentration} at the
    !% given {\tt radius} (given in units of the scale radius).
    use Gamma_Functions
    implicit none
    double precision, intent(in   ) :: alpha, concentration, radius

    if (radius >= concentration) then
       Enclosed_Mass_Einasto_Scale_Free=1.0d0
    else
       Enclosed_Mass_Einasto_Scale_Free=                                                             &
            &  Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*radius       **alpha/alpha) &
            & /Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*concentration**alpha/alpha)
    end if
    return
  end function Enclosed_Mass_Einasto_Scale_Free

  double precision function Density_Einasto_Scale_Free(radius,concentration,alpha)
    !% Returns the density (in units such that the virial mass and scale length are unity) in an Einasto dark matter profile with
    !% given {\tt concentration} and {\tt alpha} at the given {\tt radius} (given in units of the scale radius).
    use Numerical_Constants_Math
    use Gamma_Functions
    implicit none
    double precision, intent(in   ) :: alpha               , concentration, radius
    double precision                :: densityNormalization

    densityNormalization= (alpha/4.0d0/Pi)                                                                      &
         &               *    ((2.0d0/alpha)                   **(3.0d0/alpha)                                ) &
         &               *exp(-2.0d0/alpha                                                                   ) &
         &               /Gamma_Function                         (3.0d0/alpha                                 ) &
         &               /Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*concentration**alpha/alpha)
    Density_Einasto_Scale_Free=densityNormalization*exp(-(2.0d0/alpha)*(radius**alpha-1.0d0))
    return
  end function Density_Einasto_Scale_Free

  double precision function Potential_Einasto_Scale_Free(radius,concentration,alpha)
    !% Returns the gravitational potential (in units where the virial mass and scale radius are unity) in an Einasto dark matter
    !% profile with given {\tt concentration} and {\tt alpha} at the given {\tt radius} (given in units of the scale radius).
    use Gamma_Functions
    implicit none
    double precision, intent(in   ) :: alpha, concentration, radius
       if (radius <= 0.0d0) then
         Potential_Einasto_Scale_Free=                                                                  &
               & -((2.0d0/alpha)**(1.0d0/alpha))                                                        &
               & *2.0d0                                                                                 &
               & *Gamma_Function                         (2.0d0/alpha                                 ) &
               & /Gamma_Function                         (3.0d0/alpha                                 ) &
               & /Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*concentration**alpha/alpha)

       else
          Potential_Einasto_Scale_Free=                                                                                 &
               & -(                                                                                                     &
               &           Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*radius       **alpha/alpha)/radius &
               &   +((2.0d0/alpha)**(1.0d0/alpha))                                                                      &
               &   *(1.0d0+Gamma_Function_Incomplete              (2.0d0/alpha,2.0d0*radius       **alpha/alpha))       &
               &   *       Gamma_Function                         (2.0d0/alpha                                 )        &
               &   /       Gamma_Function                         (3.0d0/alpha                                 )        &
               &  )                                                                                                     &
               & /         Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*concentration**alpha/alpha)
       end if
    return
   end function Potential_Einasto_Scale_Free

  double precision function Dark_Matter_Profile_kSpace_Einasto(thisNode,wavenumber)
    !% Returns the Fourier transform of the Einasto density profile at the specified {\tt waveNumber} (given in Mpc$^{-1}$)).
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    implicit none
    type            (treeNode                      )                , intent(inout), pointer :: thisNode
    double precision                                                , intent(in   )          :: wavenumber
    class           (nodeComponentDarkMatterProfile)                               , pointer :: thisDarkMatterProfileComponent
    integer                                         , dimension(0:1)                         :: jAlpha                        , jConcentration
    double precision                                , dimension(0:1)                         :: hAlpha                        , hConcentration
    integer                                                                                  :: iAlpha                        , iConcentration
    double precision                                                                         :: alpha                         , scaleRadius        , &
         &                                                                                      virialRadiusOverScaleRadius   , wavenumberScaleFree

    ! Get components.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    ! Get scale radius, shape parameter and concentration.
    scaleRadius                =thisDarkMatterProfileComponent%scale()
    alpha                      =thisDarkMatterProfileComponent%shape()
    virialRadiusOverScaleRadius=Dark_Matter_Halo_Virial_Radius(thisNode)/scaleRadius
    wavenumberScaleFree        =wavenumber*scaleRadius

    ! Ensure the table exists and is sufficiently tabulated.
    call Fourier_Profile_Table_Make(wavenumberScaleFree,virialRadiusOverScaleRadius,alpha)

    !$omp critical(Einasto_Fourier_Interpolation)
    ! Get interpolating factors in alpha.
    jAlpha(0)=Interpolate_Locate(fourierProfileTableAlphaCount,fourierProfileTableAlpha&
         &,fourierProfileTableAlphaInterpolationAccelerator,alpha,reset=fourierProfileTableAlphaInterpolationReset)
    jAlpha(1)=jAlpha(0)+1
    hAlpha=Interpolate_Linear_Generate_Factors(fourierProfileTableAlphaCount,fourierProfileTableAlpha,jAlpha(0),alpha)

    ! Get interpolating factors in concentration.
    jConcentration(0)=Interpolate_Locate(fourierProfileTableConcentrationCount,fourierProfileTableConcentration &
         &,fourierProfileTableConcentrationInterpolationAccelerator,virialRadiusOverScaleRadius,reset &
         &=fourierProfileTableConcentrationInterpolationReset)
    jConcentration(1)=jConcentration(0)+1
    hConcentration=Interpolate_Linear_Generate_Factors(fourierProfileTableConcentrationCount,fourierProfileTableConcentration&
         &,jConcentration(0),virialRadiusOverScaleRadius)

    ! Find the Fourier profile by interpolation.
    Dark_Matter_Profile_kSpace_Einasto=0.0d0
    do iAlpha=0,1
       do iConcentration=0,1
          Dark_Matter_Profile_kSpace_Einasto=Dark_Matter_Profile_kSpace_Einasto+Interpolate(fourierProfileTableWavenumberCount &
               &,fourierProfileTableWavenumber,fourierProfileTable(:,jConcentration(iConcentration),jAlpha(iAlpha))&
               &,fourierProfileTableWavenumberInterpolationObject ,fourierProfileTableWavenumberInterpolationAccelerator&
               &,wavenumberScaleFree,reset =fourierProfileTableWavenumberInterpolationReset)*hAlpha(iAlpha)*hConcentration(iConcentration)
       end do
    end do
    !$omp end critical(Einasto_Fourier_Interpolation)
    return
  end function Dark_Matter_Profile_kSpace_Einasto

  subroutine Fourier_Profile_Table_Make(wavenumberRequired,concentrationRequired,alphaRequired)
    !% Create a tabulation of the Fourier transform of Einasto profiles as a function of their $\alpha$ parameter and
    !% dimensionless wavenumber.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Interpolation
    use Numerical_Integration
    use Numerical_Ranges
    use Memory_Management
    use Galacticus_Display
    implicit none
    double precision                            , intent(in   ) :: alphaRequired       , concentrationRequired, wavenumberRequired
    integer                                                     :: iAlpha              , iConcentration       , iWavenumber       , &
         &                                                         percentage
    logical                                                     :: makeTable
    double precision                                            :: alpha               , concentration        , radiusMaximum     , &
         &                                                         radiusMinimum       , wavenumber
    type            (c_ptr                     )                :: parameterPointer
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace

    !$omp critical (Einasto_Fourier_Interpolation)
    ! Assume table does not need remaking.
    makeTable=.false.
    ! Check for uninitialized table.
    if (.not.fourierProfileTableInitialized) makeTable=.true.
    ! Check for alpha out of range.
    if (alphaRequired         < fourierProfileTableAlphaMinimum         .or. alphaRequired         > fourierProfileTableAlphaMaximum        ) then
       makeTable=.true.
       ! Compute the range of tabulation and number of points to use.
       fourierProfileTableAlphaMinimum        =min(fourierProfileTableAlphaMinimum     ,0.9d0*alphaRequired           )
       fourierProfileTableAlphaMaximum        =max(fourierProfileTableAlphaMaximum     ,1.1d0*alphaRequired           )
    end if
    ! Check for concentration out of range.
    if (concentrationRequired < fourierProfileTableConcentrationMinimum .or. concentrationRequired > fourierProfileTableConcentrationMaximum ) then
       makeTable=.true.
       ! Compute the range of tabulation and number of points to use.
       fourierProfileTableConcentrationMinimum=min(fourierProfileTableConcentrationMinimum,0.5d0*concentrationRequired)
       fourierProfileTableConcentrationMaximum=max(fourierProfileTableConcentrationMaximum,2.0d0*concentrationRequired)
    end if
    ! Check for wavenumber below minimum tabulated value.
    if (wavenumberRequired    < fourierProfileTableWavenumberMinimum    .or. wavenumberRequired    > fourierProfileTableWavenumberMaximum    ) then
       makeTable=.true.
       fourierProfileTableWavenumberMinimum   =min(fourierProfileTableWavenumberMinimum   ,0.5d0*wavenumberRequired   )
       fourierProfileTableWavenumberMaximum   =max(fourierProfileTableWavenumberMaximum   ,2.0d0*wavenumberRequired   )
    end if
    ! Remake the table if necessary.
    if (makeTable) then
       ! Display a message.
       call Galacticus_Display_Indent('Constructing Einasto profile Fourier transform lookup table...',verbosityInfo)
       ! Allocate arrays to the appropriate sizes.
       fourierProfileTableAlphaCount        =int(      (fourierProfileTableAlphaMaximum        -fourierProfileTableAlphaMinimum        ) &
            &*dble(fourierProfileTableAlphaPointsPerUnit          ))+1
       fourierProfileTableConcentrationCount=int(log10(fourierProfileTableConcentrationMaximum/fourierProfileTableConcentrationMinimum) &
            &*dble(fourierProfileTableConcentrationPointsPerDecade))+1
       fourierProfileTableWavenumberCount   =int(log10(fourierProfileTableWavenumberMaximum   /fourierProfileTableWavenumberMinimum   ) &
            &*dble(fourierProfileTableWavenumberPointsPerDecade   ))+1
       if (allocated(fourierProfileTableAlpha        )) call Dealloc_Array(fourierProfileTableAlpha        )
       if (allocated(fourierProfileTableConcentration)) call Dealloc_Array(fourierProfileTableConcentration)
       if (allocated(fourierProfileTableWavenumber   )) call Dealloc_Array(fourierProfileTableWavenumber   )
       if (allocated(fourierProfileTable             )) call Dealloc_Array(fourierProfileTable             )
       call Alloc_Array(fourierProfileTableAlpha        ,[                                                                         fourierProfileTableAlphaCount])
       call Alloc_Array(fourierProfileTableConcentration,[                                   fourierProfileTableConcentrationCount                              ])
       call Alloc_Array(fourierProfileTableWavenumber   ,[fourierProfileTableWavenumberCount                                                                    ])
       call Alloc_Array(fourierProfileTable             ,[fourierProfileTableWavenumberCount,fourierProfileTableConcentrationCount,fourierProfileTableAlphaCount])
       ! Create ranges of alpha and wavenumber.
       fourierProfileTableAlpha        =Make_Range(fourierProfileTableAlphaMinimum        ,fourierProfileTableAlphaMaximum         &
            &,fourierProfileTableAlphaCount        ,rangeType=rangeTypeLinear     )
       fourierProfileTableConcentration=Make_Range(fourierProfileTableConcentrationMinimum,fourierProfileTableConcentrationMaximum &
            &,fourierProfileTableConcentrationCount,rangeType=rangeTypeLogarithmic)
       fourierProfileTableWavenumber   =Make_Range(fourierProfileTableWavenumberMinimum   ,fourierProfileTableWavenumberMaximum    &
            &,fourierProfileTableWavenumberCount   ,rangeType=rangeTypeLogarithmic)
       ! Tabulate the Fourier profile.
       do iAlpha=1,fourierProfileTableAlphaCount
          alpha=fourierProfileTableAlpha(iAlpha)
          do iConcentration=1,fourierProfileTableConcentrationCount
             concentration=fourierProfileTableConcentration(iConcentration)

             ! Show progress.
             percentage=int(100.0d0*dble((iAlpha-1)*fourierProfileTableConcentrationCount+iConcentration-1)&
                  &/dble(fourierProfileTableAlphaCount*fourierProfileTableConcentrationCount))
             call Galacticus_Display_Counter(percentage,iAlpha == 1 .and. iConcentration == 1,verbosityWorking)

             do iWavenumber=1,fourierProfileTableWavenumberCount
                wavenumber=fourierProfileTableWavenumber(iWavenumber)
                ! Compute the potential fourierProfile.
                radiusMinimum         =0.0d0
                radiusMaximum         =concentration
                wavenumberParameter   =wavenumber
                alphaParameter        =alpha
                concentrationParameter=concentration
                fourierProfileTable(iWavenumber,iConcentration,iAlpha)=Integrate(radiusMinimum,radiusMaximum&
                     &,Fourier_Profile_Integrand_Einasto ,parameterPointer,integrandFunction,integrationWorkspace&
                     &,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-2,maxIntervals=10000)
                call Integrate_Done(integrandFunction,integrationWorkspace)

             end do
          end do
       end do
       call Galacticus_Display_Counter_Clear(verbosityWorking)
       ! Reset interpolators.
       call Interpolate_Done(fourierProfileTableWavenumberInterpolationObject,fourierProfileTableWavenumberInterpolationAccelerator &
            &,fourierProfileTableWavenumberInterpolationReset)
       call Interpolate_Done(interpolationAccelerator=fourierProfileTableAlphaInterpolationAccelerator         &
            &,reset=fourierProfileTableAlphaInterpolationReset        )
       call Interpolate_Done(interpolationAccelerator=fourierProfileTableConcentrationInterpolationAccelerator &
            &,reset=fourierProfileTableConcentrationInterpolationReset)
       fourierProfileTableWavenumberInterpolationReset   =.true.
       fourierProfileTableAlphaInterpolationReset        =.true.
       fourierProfileTableConcentrationInterpolationReset=.true.
       ! Flag that the table is now initialized.
       fourierProfileTableInitialized=.true.
       ! Display a message.
       call Galacticus_Display_Unindent('done',verbosityInfo)
    end if
    !$omp end critical (Einasto_Fourier_Interpolation)
    return
  end subroutine Fourier_Profile_Table_Make

  function Fourier_Profile_Integrand_Einasto(radius,parameterPointer) bind(c)
    !% Integrand for Einasto Fourier profile.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Math
    implicit none
    real(kind=c_double)        :: Fourier_Profile_Integrand_Einasto
    real(kind=c_double), value :: radius
    type(c_ptr        ), value :: parameterPointer

    Fourier_Profile_Integrand_Einasto=4.0d0*Pi*radius*sin(wavenumberParameter*radius)*Density_Einasto_Scale_Free(radius&
         &,concentrationParameter,alphaParameter)/wavenumberParameter
    return
  end function Fourier_Profile_Integrand_Einasto

  double precision function Dark_Matter_Profile_Freefall_Radius_Einasto(thisNode,time)
    !% Returns the freefall radius in the Einasto density profile at the specified {\tt time} (given in Gyr).
    use Numerical_Interpolation
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    implicit none
    type            (treeNode                      )                , intent(inout), pointer :: thisNode
    double precision                                                , intent(in   )          :: time
    class           (nodeComponentBasic            )                               , pointer :: thisBasicComponent
    class           (nodeComponentDarkMatterProfile)                               , pointer :: thisDarkMatterProfileComponent
    integer                                         , dimension(0:1)                         :: jAlpha
    double precision                                , dimension(0:1)                         :: hAlpha
    integer                                                                                  :: iAlpha
    double precision                                                                         :: alpha                         , freefallTimeScaleFree, &
         &                                                                                      radiusScale                   , timeScale            , &
         &                                                                                      velocityScale

    ! For non-positive freefall times, return a zero freefall radius immediately.
    if (time <= 0.0d0) then
       Dark_Matter_Profile_Freefall_Radius_Einasto=0.0d0
       return
    end if

    ! Get components.
    thisBasicComponent             => thisNode%basic            (                 )
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    ! Get the shape parameter.
    alpha        =thisDarkMatterProfileComponent%shape()

    ! Get the scale radius.
    radiusScale  =thisDarkMatterProfileComponent%scale()

    ! Get the velocity scale.
    velocityScale=sqrt(gravitationalConstantGalacticus*thisBasicComponent%mass()/radiusScale)

    ! Compute time scale.
    timeScale=Mpc_per_km_per_s_To_Gyr*radiusScale/velocityScale

    ! Compute dimensionless time.
    freefallTimeScaleFree=time/timeScale

    ! Ensure table is sufficiently extensive.
    call Dark_Matter_Profile_Einasto_Freefall_Tabulate(freefallTimeScaleFree,alpha)

    ! Interpolate to get the freefall radius.
    !$omp critical(Einasto_Freefall_Interpolation)
    ! Get interpolating factors in alpha.
    jAlpha(0)=Interpolate_Locate(freefallRadiusTableAlphaCount,freefallRadiusTableAlpha&
         &,freefallRadiusTableAlphaInterpolationAccelerator,alpha,reset=freefallRadiusTableAlphaInterpolationReset)
    jAlpha(1)=jAlpha(0)+1
    hAlpha=Interpolate_Linear_Generate_Factors(freefallRadiusTableAlphaCount,freefallRadiusTableAlpha,jAlpha(0),alpha)

    Dark_Matter_Profile_Freefall_Radius_Einasto=0.0d0
    do iAlpha=0,1
       Dark_Matter_Profile_Freefall_Radius_Einasto=Interpolate(freefallRadiusTableRadiusCount,freefallRadiusTable(:&
            &,jAlpha(iAlpha)),freefallRadiusTableRadius ,freefallRadiusTableRadiusInterpolationObject&
            &,freefallRadiusTableRadiusInterpolationAccelerator ,freefallTimeScaleFree,reset&
            &=freefallRadiusTableRadiusInterpolationReset)*hAlpha(iAlpha)
    end do
    Dark_Matter_Profile_Freefall_Radius_Einasto=Dark_Matter_Profile_Freefall_Radius_Einasto*radiusScale
    !$omp end critical(Einasto_Freefall_Interpolation)
    return
  end function Dark_Matter_Profile_Freefall_Radius_Einasto

  double precision function Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Einasto(thisNode,time)
    !% Returns the rate of increase of the freefall radius in the Einasto density profile at the specified {\tt time} (given in
    !% Gyr).
    use Numerical_Interpolation
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    implicit none
    type            (treeNode                      )                , intent(inout), pointer :: thisNode
    double precision                                                , intent(in   )          :: time
    class           (nodeComponentBasic            )                               , pointer :: thisBasicComponent
    class           (nodeComponentDarkMatterProfile)                               , pointer :: thisDarkMatterProfileComponent
    integer                                         , dimension(0:1)                         :: jAlpha
    double precision                                , dimension(0:1)                         :: hAlpha
    integer                                                                                  :: iAlpha
    double precision                                                                         :: alpha                         , freefallTimeScaleFree, &
         &                                                                                      radiusScale                   , timeScale            , &
         &                                                                                      velocityScale

    ! For non-positive freefall times, return a zero freefall radius immediately.
    if (time <= 0.0d0) then
       Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Einasto=0.0d0
       return
    end if

    ! Get components.
    thisBasicComponent             => thisNode%basic            (                 )
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)

    ! Get the shape parameter.
    alpha        =thisDarkMatterProfileComponent%shape()

    ! Get the scale radius.
    radiusScale  =thisDarkMatterProfileComponent%scale()

    ! Get the velocity scale.
    velocityScale=sqrt(gravitationalConstantGalacticus*thisBasicComponent%mass()/radiusScale)

    ! Compute time scale.
    timeScale=Mpc_per_km_per_s_To_Gyr*radiusScale/velocityScale

    ! Compute dimensionless time.
    freefallTimeScaleFree=time/timeScale

    ! Ensure table is sufficiently extensive.
    call Dark_Matter_Profile_Einasto_Freefall_Tabulate(freefallTimeScaleFree,alpha)

    ! Interpolate to get the freefall radius.
    !$omp critical(Einasto_Freefall_Interpolation)
    ! Get interpolating factors in alpha.
    jAlpha(0)=Interpolate_Locate(freefallRadiusTableAlphaCount,freefallRadiusTableAlpha&
         &,freefallRadiusTableAlphaInterpolationAccelerator,alpha,reset=freefallRadiusTableAlphaInterpolationReset)
    jAlpha(1)=jAlpha(0)+1
    hAlpha=Interpolate_Linear_Generate_Factors(freefallRadiusTableAlphaCount,freefallRadiusTableAlpha,jAlpha(0),alpha)

    Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Einasto=0.0d0
    do iAlpha=0,1
       Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Einasto=Interpolate_Derivative(freefallRadiusTableRadiusCount,freefallRadiusTable(:&
            &,jAlpha(iAlpha)),freefallRadiusTableRadius ,freefallRadiusTableRadiusInterpolationObject&
            &,freefallRadiusTableRadiusInterpolationAccelerator ,freefallTimeScaleFree,reset&
            &=freefallRadiusTableRadiusInterpolationReset)*hAlpha(iAlpha)
    end do
    Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Einasto=Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Einasto*radiusScale/timeScale
    !$omp end critical(Einasto_Freefall_Interpolation)
    return
  end function Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Einasto

  subroutine Dark_Matter_Profile_Einasto_Freefall_Tabulate(freefallTimeScaleFree,alphaRequired)
    !% Tabulates the freefall time vs. freefall radius for Einasto halos.
    use Galacticus_Display
    use Numerical_Ranges
    use Memory_Management
    use Numerical_Interpolation
    implicit none
    double precision, intent(in   ) :: alphaRequired, freefallTimeScaleFree
    logical                         :: retabulate
    integer                         :: iAlpha       , iRadius              , percentage
    double precision                :: alpha

    !$omp critical (Einasto_Freefall_Interpolation)
    retabulate=.not.freefallRadiusTableInitialized
    ! If the table has not yet been made, compute and store the freefall times corresponding to the minimum and maximum
    ! radii that will be tabulated by default.
    if (retabulate) then
       freefallTimeMinimum=Freefall_Time_Scale_Free(freefallRadiusTableRadiusMinimum,alphaRequired)
       freefallTimeMaximum=Freefall_Time_Scale_Free(freefallRadiusTableRadiusMaximum,alphaRequired)
    end if
    do while (freefallTimeScaleFree < freefallTimeMinimum)
       freefallRadiusTableRadiusMinimum=0.5d0*freefallRadiusTableRadiusMinimum
       freefallTimeMinimum=Freefall_Time_Scale_Free(freefallRadiusTableRadiusMinimum,alphaRequired)
       retabulate=.true.
    end do
    do while (freefallTimeScaleFree > freefallTimeMaximum)
       freefallRadiusTableRadiusMaximum=2.0d0*freefallRadiusTableRadiusMaximum
       freefallTimeMaximum=Freefall_Time_Scale_Free(freefallRadiusTableRadiusMaximum,alphaRequired)
       retabulate=.true.
    end do
    ! Check for alpha out of range.
    if (alphaRequired < freefallRadiusTableAlphaMinimum .or. alphaRequired > freefallRadiusTableAlphaMaximum) then
       retabulate=.true.
       ! Compute the range of tabulation.
       freefallRadiusTableAlphaMinimum=min(freefallRadiusTableAlphaMinimum,0.9d0*alphaRequired)
       freefallRadiusTableAlphaMaximum=max(freefallRadiusTableAlphaMaximum,1.1d0*alphaRequired)
    end if

    if (retabulate) then
       ! Display a message.
       call Galacticus_Display_Indent('Constructing Einasto profile freefall radius lookup table...',verbosityWorking)
       ! Decide how many points to tabulate and allocate table arrays.
       freefallRadiusTableRadiusCount=int(log10(freefallRadiusTableRadiusMaximum/freefallRadiusTableRadiusMinimum)*dble(freefallRadiusTableRadiusPointsPerDecade))+1
       freefallRadiusTableAlphaCount =int(      (freefallRadiusTableAlphaMaximum -freefallRadiusTableAlphaMinimum )*dble(freefallRadiusTableAlphaPointsPerUnit   ))+1
       if (allocated(freefallRadiusTableRadius)) then
          call Dealloc_Array(freefallRadiusTableAlpha )
          call Dealloc_Array(freefallRadiusTableRadius)
          call Dealloc_Array(freefallRadiusTable      )
       end if
       call Alloc_Array(freefallRadiusTableAlpha ,[                               freefallRadiusTableAlphaCount])
       call Alloc_Array(freefallRadiusTableRadius,[freefallRadiusTableRadiusCount                              ])
       call Alloc_Array(freefallRadiusTable      ,[freefallRadiusTableRadiusCount,freefallRadiusTableAlphaCount])
       ! Create a range of radii and alpha.
       freefallRadiusTableAlpha =Make_Range(freefallRadiusTableAlphaMinimum ,freefallRadiusTableAlphaMaximum ,freefallRadiusTableAlphaCount ,rangeType=rangeTypeLinear     )
       freefallRadiusTableRadius=Make_Range(freefallRadiusTableRadiusMinimum,freefallRadiusTableRadiusMaximum,freefallRadiusTableRadiusCount,rangeType=rangeTypeLogarithmic)
       ! Loop over radii and alpha and populate tables.
       do iAlpha=1,freefallRadiusTableAlphaCount
          alpha=freefallRadiusTableAlpha(iAlpha)
          do iRadius=1,freefallRadiusTableRadiusCount
             ! Show progress.
             percentage=int(100.0d0*dble((iAlpha-1)*freefallRadiusTableRadiusCount+iRadius-1)&
                  &/dble(freefallRadiusTableAlphaCount*freefallRadiusTableRadiusCount))
             call Galacticus_Display_Counter(percentage,iAlpha == 1 .and. iRadius == 1,verbosityWorking)
             ! Compute the freefall radius.
             freefallRadiusTable(iRadius,iAlpha)=Freefall_Time_Scale_Free(freefallRadiusTableRadius(iRadius),alpha)
          end do
       end do
       ! Ensure interpolations get reset.
       call Interpolate_Done(                                                                            &
            &                interpolationAccelerator=freefallRadiusTableAlphaInterpolationAccelerator,  &
            &                reset                   =freefallRadiusTableAlphaInterpolationReset         &
            &               )
       call Interpolate_Done(                                                                            &
            &                interpolationObject     =freefallRadiusTableRadiusInterpolationObject,      &
            &                interpolationAccelerator=freefallRadiusTableRadiusInterpolationAccelerator, &
            &                reset                   =freefallRadiusTableRadiusInterpolationReset        &
            &               )
       freefallRadiusTableAlphaInterpolationReset =.true.
       freefallRadiusTableRadiusInterpolationReset=.true.
       ! Store the minimum and maximum tabulated freefall times across all alpha values.
       freefallTimeMinimum=maxval(freefallRadiusTable(                             1,:))
       freefallTimeMaximum=minval(freefallRadiusTable(freefallRadiusTableRadiusCount,:))
       ! Display a message.
       call Galacticus_Display_Unindent('...done',verbosityWorking)
       ! Specify that tabulation has been made.
       freefallRadiusTableInitialized=.true.
    end if
    !$omp end critical (Einasto_Freefall_Interpolation)
    return
  end subroutine Dark_Matter_Profile_Einasto_Freefall_Tabulate

  double precision function Freefall_Time_Scale_Free(radius,alpha)
    !% Compute the freefall time in a scale-free Einasto halo.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    implicit none
    double precision                            , intent(in   ) :: alpha               , radius
    type            (c_ptr                     )                :: parameterPointer
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    double precision                                            :: radiusEnd

    radiusStart   =radius
    radiusEnd     =0.0d0
    alphaParameter=alpha
    Freefall_Time_Scale_Free=Integrate(radiusEnd,radiusStart,Freefall_Time_Scale_Free_Integrand_Einasto,parameterPointer&
         &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    return
  end function Freefall_Time_Scale_Free

  function Freefall_Time_Scale_Free_Integrand_Einasto(radius,parameterPointer) bind(c)
    !% Integrand function used for finding the free-fall time in Einasto halos.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(kind=c_double)        :: Freefall_Time_Scale_Free_Integrand_Einasto
    real(kind=c_double), value :: radius
    type(c_ptr        ), value :: parameterPointer

    Freefall_Time_Scale_Free_Integrand_Einasto= 1.0d0                                                                   &
         &                                     /sqrt(                                                                  &
         &                                             2.0d0                                                            &
         &                                            *(                                                                &
         &                                               Potential_Einasto_Scale_Free(radiusStart,1.0d0,alphaParameter) &
         &                                              -Potential_Einasto_Scale_Free(radius     ,1.0d0,alphaParameter) &
         &                                             )                                                                &
         &                                           )
    return
  end function Freefall_Time_Scale_Free_Integrand_Einasto

  !# <galacticusStateStoreTask>
  !#  <unitName>Dark_Matter_Profiles_Einasto_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Dark_Matter_Profiles_Einasto_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    write (stateFile) angularMomentumTableRadiusMinimum,angularMomentumTableRadiusMaximum,angularMomentumTableAlphaMinimum &
         &,angularMomentumTableAlphaMaximum,energyTableConcentrationMinimum,energyTableConcentrationMaximum &
         &,energyTableAlphaMinimum,energyTableAlphaMaximum,fourierProfileTableWavenumberMinimum &
         &,fourierProfileTableWavenumberMaximum,fourierProfileTableAlphaMinimum,fourierProfileTableAlphaMaximum &
         &,fourierProfileTableConcentrationMinimum,fourierProfileTableConcentrationMaximum,freefallRadiusTableRadiusMinimum&
         &,freefallRadiusTableRadiusMaximum,freefallRadiusTableAlphaMinimum ,freefallRadiusTableAlphaMaximum

    return
  end subroutine Dark_Matter_Profiles_Einasto_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Dark_Matter_Profiles_Einasto_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Dark_Matter_Profiles_Einasto_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) angularMomentumTableRadiusMinimum,angularMomentumTableRadiusMaximum,angularMomentumTableAlphaMinimum &
         &,angularMomentumTableAlphaMaximum,energyTableConcentrationMinimum,energyTableConcentrationMaximum &
         &,energyTableAlphaMinimum,energyTableAlphaMaximum,fourierProfileTableWavenumberMinimum &
         &,fourierProfileTableWavenumberMaximum,fourierProfileTableAlphaMinimum,fourierProfileTableAlphaMaximum &
         &,fourierProfileTableConcentrationMinimum,fourierProfileTableConcentrationMaximum,freefallRadiusTableRadiusMinimum &
         &,freefallRadiusTableRadiusMaximum,freefallRadiusTableAlphaMinimum,freefallRadiusTableAlphaMaximum
    ! Retabulate.
    angularMomentumTableInitialized=.false.
    energyTableInitialized         =.false.
    fourierProfileTableInitialized =.false.
    freefallRadiusTableInitialized =.false.
    return
  end subroutine Dark_Matter_Profiles_Einasto_State_Retrieve

end module Dark_Matter_Profiles_Einasto
