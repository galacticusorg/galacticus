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

  !% An implementation of ``Einasto'' dark matter halo profiles.

  !# <darkMatterProfile name="darkMatterProfileEinasto">
  !#  <description>``Einasto'' dark matter halo profiles</description>
  !# </darkMatterProfile>

  use Dark_Matter_Halo_Scales
  use Tables
  use Kind_Numbers

  type, extends(darkMatterProfileClass) :: darkMatterProfileEinasto
     !% A dark matter halo profile class implementing ``Einasto'' dark matter halos.
     private
     ! Tables for specific angular momentum vs. radius table
     double precision                                                   :: angularMomentumTableRadiusMinimum                    
     double precision                                                   :: angularMomentumTableRadiusMaximum                    
     double precision                                                   :: angularMomentumTableAlphaMinimum                     
     double precision                                                   :: angularMomentumTableAlphaMaximum                     
     logical                                                            :: angularMomentumTableInitialized                      
     integer                                                            :: angularMomentumTableAlphaCount                       , angularMomentumTableRadiusCount
     double precision                   , allocatable, dimension(:  )   :: angularMomentumTableAlpha                            , angularMomentumTableRadius
     double precision                   , allocatable, dimension(:,:)   :: angularMomentumTable
     type            (fgsl_interp      )                                :: angularMomentumTableRadiusInterpolationObject
     type            (fgsl_interp_accel)                                :: angularMomentumTableAlphaInterpolationAccelerator    , angularMomentumTableRadiusInterpolationAccelerator
     logical                                                            :: angularMomentumTableAlphaInterpolationReset          , angularMomentumTableRadiusInterpolationReset            
     ! Tables for freefall time vs. radius table
     double precision                                                   :: freefallRadiusTableRadiusMinimum                     
     double precision                                                   :: freefallRadiusTableRadiusMaximum                     
     double precision                                                   :: freefallRadiusTableAlphaMinimum                      
     double precision                                                   :: freefallRadiusTableAlphaMaximum                      
     logical                                                            :: freefallRadiusTableInitialized                       
     integer                                                            :: freefallRadiusTableAlphaCount                        , freefallRadiusTableRadiusCount
     double precision                   , allocatable, dimension(:  )   :: freefallRadiusTableAlpha                             , freefallRadiusTableRadius
     double precision                   , allocatable, dimension(:,:)   :: freefallRadiusTable
     double precision                                                   :: freefallTimeMaximum                                  , freefallTimeMinimum
     type            (fgsl_interp      )                                :: freefallRadiusTableRadiusInterpolationObject
     type            (fgsl_interp_accel)                                :: freefallRadiusTableAlphaInterpolationAccelerator     , freefallRadiusTableRadiusInterpolationAccelerator
     logical                                                            :: freefallRadiusTableAlphaInterpolationReset           , freefallRadiusTableRadiusInterpolationReset             
     ! Tables for energy as a function of concentration and alpha.
     double precision                                                   :: energyTableConcentrationMinimum                      
     double precision                                                   :: energyTableConcentrationMaximum                      
     double precision                                                   :: energyTableAlphaMinimum                              
     double precision                                                   :: energyTableAlphaMaximum                              
     logical                                                            :: energyTableInitialized                               
     integer                                                            :: energyTableAlphaCount                                , energyTableConcentrationCount
     double precision                   , allocatable, dimension(:  )   :: energyTableAlpha                                     , energyTableConcentration
     double precision                   , allocatable, dimension(:,:)   :: energyTable
     type            (fgsl_interp      )                                :: energyTableAlphaInterpolationObject                  , energyTableConcentrationInterpolationObject
     type            (fgsl_interp_accel)                                :: energyTableAlphaInterpolationAccelerator             , energyTableConcentrationInterpolationAccelerator
     logical                                                            :: energyTableAlphaInterpolationReset                   , energyTableConcentrationInterpolationReset              
     ! Tables for specific Fourier transform of density profile as a function of alpha and radius.
     double precision                                                   :: fourierProfileTableConcentrationMinimum              

     double precision                                                   :: fourierProfileTableConcentrationMaximum              
     double precision                                                   :: fourierProfileTableWavenumberMinimum                 
     double precision                                                   :: fourierProfileTableWavenumberMaximum                 
     double precision                                                   :: fourierProfileTableAlphaMinimum                      
     double precision                                                   :: fourierProfileTableAlphaMaximum                      
     logical                                                            :: fourierProfileTableInitialized                       
     integer                                                            :: fourierProfileTableAlphaCount                        , fourierProfileTableConcentrationCount                   , &
          &                                                                fourierProfileTableWavenumberCount
     double precision                   , allocatable, dimension(:    ) :: fourierProfileTableAlpha                             , fourierProfileTableConcentration                        , &
       &                                                                   fourierProfileTableWavenumber
     double precision                   , allocatable, dimension(:,:,:) :: fourierProfileTable
     type            (fgsl_interp      )                                :: fourierProfileTableWavenumberInterpolationObject
     type            (fgsl_interp_accel)                                :: fourierProfileTableAlphaInterpolationAccelerator     , fourierProfileTableConcentrationInterpolationAccelerator, &
          &                                                                fourierProfileTableWavenumberInterpolationAccelerator
     logical                                                            :: fourierProfileTableAlphaInterpolationReset           , fourierProfileTableConcentrationInterpolationReset      , &
          &                                                                fourierProfileTableWavenumberInterpolationReset      
     ! Pointer to object setting halo scales.
     class(darkMatterHaloScaleClass    ), pointer                       :: scale
   contains
     !@ <objectMethods>
     !@   <object>darkMatterProfileEinasto</object>
     !@   <objectMethod>
     !@     <method>densityScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin, \doublezero\ concentration\argin, \doublezero\ alpha\argin</arguments>
     !@     <description>Returns the density (in units such that the virial mass and scale length are unity) in an Einasto dark matter profile with given {\tt concentration} and {\tt alpha} at the given {\tt radius} (given in units of the scale radius).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>enclosedMassScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin, \doublezero\ concentration\argin, \doublezero\ alpha\argin</arguments>
     !@     <description>Returns the enclosed mass (in units of the virial mass) in an Einasto dark matter profile with given {\tt concentration} at the given {\tt radius} (given in units of the scale radius).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>energyTableMake</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ concentrationRequired\argin, \doublezero\ alphaRequired\argin</arguments>
     !@     <description>Create a tabulation of the energy of Einasto profiles as a function of their concentration of $\alpha$ parameter.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>fourierProfileTableMake</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ wavenumberRequired\argin, \doublezero\ concentrationRequired\argin, \doublezero\ alphaRequired\argin</arguments>
     !@     <description>Create a tabulation of the Fourier transform of Einasto profiles as a function of their $\alpha$ parameter and dimensionless wavenumber.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>freefallTabulate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ freefallTimeScaleFree\argin, \doublezero\ alphaRequired\argin</arguments>
     !@     <description>Tabulates the freefall time vs. freefall radius for Einasto halos.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>freefallTimeScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin, \doublezero\ alpha\argin</arguments>
     !@     <description>Compute the freefall time in a scale-free Einasto halo.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>potentialScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin, \doublezero\ concentration\argin, \doublezero\ alpha\argin</arguments>
     !@     <description>Returns the gravitational potential (in units where the virial mass and scale radius are unity) in an Einasto dark matter profile with given {\tt concentration} and {\tt alpha} at the given {\tt radius} (given in units of the scale radius).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radiusFromSpecificAngularMomentumScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ alpha\argin, \doublezero\ specificAngularMomentumScaleFree\argin</arguments>
     !@     <description> Comptue the radius at which a circular orbit has the given {\tt specificAngularMomentumScaleFree} in a scale free Einasto profile.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radiusFromSpecificAngularMomentumTableMake</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ alphaRequired\argin, \doublezero\ specificAngularMomentumRequired\argin</arguments>
     !@     <description>Create a tabulation of the relation between specific angular momentum and radius in an Einasto profile.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final                                                      einastoDestructor
     procedure :: calculationReset                           => einastoCalculationReset
     procedure :: stateStore                                 => einastoStateStore
     procedure :: stateRestore                               => einastoStateRestore
     procedure :: density                                    => einastoDensity
     procedure :: enclosedMass                               => einastoEnclosedMass
     procedure :: potential                                  => einastoPotential
     procedure :: circularVelocity                           => einastoCircularVelocity
     procedure :: circularVelocityMaximum                    => einastoCircularVelocityMaximum
     procedure :: radiusFromSpecificAngularMomentum          => einastoRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization                      => einastoRotationNormalization
     procedure :: energy                                     => einastoEnergy
     procedure :: energyGrowthRate                           => einastoEnergyGrowthRate
     procedure :: kSpace                                     => einastoKSpace
     procedure :: freefallRadius                             => einastoFreefallRadius
     procedure :: freefallRadiusIncreaseRate                 => einastoFreefallRadiusIncreaseRate
     procedure :: densityScaleFree                           => einastoDensityScaleFree
     procedure :: enclosedMassScaleFree                      => einastoEnclosedMassScaleFree
     procedure :: freefallTimeScaleFree                      => einastoFreefallTimeScaleFree
     procedure :: potentialScaleFree                         => einastoPotentialScaleFree
     procedure :: radiusFromSpecificAngularMomentumScaleFree => einastoRadiusFromSpecificAngularMomentumScaleFree
     procedure :: radiusFromSpecificAngularMomentumTableMake => einastoRadiusFromSpecificAngularMomentumTableMake
     procedure :: freefallTabulate                           => einastoFreefallTabulate
     procedure :: energyTableMake                            => einastoEnergyTableMake
     procedure :: fourierProfileTableMake                    => einastoFourierProfileTableMake
  end type darkMatterProfileEinasto

  interface darkMatterProfileEinasto
     !% Constructors for the {\tt einasto} dark matter halo profile class.
     module procedure einastoDefaultConstructor
     module procedure einastoConstructor
  end interface darkMatterProfileEinasto

  ! Granularity parameters for tabulations.
  integer, parameter :: einastoAngularMomentumTableRadiusPointsPerDecade      =100
  integer, parameter :: einastoAngularMomentumTableAlphaPointsPerUnit         =100
  integer, parameter :: einastoFreefallRadiusTableRadiusPointsPerDecade       = 10
  integer, parameter :: einastoFreefallRadiusTableAlphaPointsPerUnit          = 30
  integer, parameter :: einastoEnergyTableConcentrationPointsPerDecade        =100
  integer, parameter :: einastoEnergyTableAlphaPointsPerUnit                  =100
  integer, parameter :: einastoFourierProfileTableConcentrationPointsPerDecade= 10
  integer, parameter :: einastoFourierProfileTableWavenumberPointsPerDecade   = 10
  integer, parameter :: einastoFourierProfileTableAlphaPointsPerUnit          =100

contains

  function einastoDefaultConstructor()
    !% Default constructor for the {\tt einasto} dark matter halo profile class.
    use Input_Parameters
    implicit none
    type(darkMatterProfileEinasto), target :: einastoDefaultConstructor

    einastoDefaultConstructor=einastoConstructor(darkMatterHaloScale())
    return
  end function einastoDefaultConstructor

  function einastoConstructor(scale)
    !% Generic constructor for the {\tt einasto} dark matter halo profile class.
    use Galacticus_Error
    use Array_Utilities
    implicit none
    type (darkMatterProfileEinasto), target :: einastoConstructor
    class(darkMatterHaloScaleClass), target :: scale

    ! Initialize table states.
    einastoConstructor%angularMomentumTableRadiusMinimum                 = 1.0d-3
    einastoConstructor%angularMomentumTableRadiusMaximum                 =20.0d+0
    einastoConstructor%angularMomentumTableAlphaMinimum                  = 0.1d+0
    einastoConstructor%angularMomentumTableAlphaMaximum                  = 0.3d+0
    einastoConstructor%angularMomentumTableInitialized                   =.false.
    einastoConstructor%angularMomentumTableAlphaInterpolationReset       =.true.
    einastoConstructor%angularMomentumTableRadiusInterpolationReset      =.true.
    einastoConstructor%freefallRadiusTableRadiusMinimum                  = 1.0d-3
    einastoConstructor%freefallRadiusTableRadiusMaximum                  =20.0d+0
    einastoConstructor%freefallRadiusTableAlphaMinimum                   = 0.1d+0
    einastoConstructor%freefallRadiusTableAlphaMaximum                   = 0.3d+0
    einastoConstructor%freefallRadiusTableInitialized                    =.false.
    einastoConstructor%freefallRadiusTableAlphaInterpolationReset        =.true.
    einastoConstructor%freefallRadiusTableRadiusInterpolationReset       =.true.
    einastoConstructor%energyTableConcentrationMinimum                   = 2.0d0
    einastoConstructor%energyTableConcentrationMaximum                   =20.0d0
    einastoConstructor%energyTableAlphaMinimum                           = 0.1d0
    einastoConstructor%energyTableAlphaMaximum                           = 0.3d0
    einastoConstructor%energyTableInitialized                            =.false.
    einastoConstructor%energyTableAlphaInterpolationReset                =.true.
    einastoConstructor%energyTableConcentrationInterpolationReset        =.true.
    einastoConstructor%fourierProfileTableConcentrationMinimum           = 2.0d0
    einastoConstructor%fourierProfileTableConcentrationMaximum           =20.0d0
    einastoConstructor%fourierProfileTableWavenumberMinimum              = 1.0d-3
    einastoConstructor%fourierProfileTableWavenumberMaximum              = 1.0d+3
    einastoConstructor%fourierProfileTableAlphaMinimum                   = 0.1d+0
    einastoConstructor%fourierProfileTableAlphaMaximum                   = 0.3d+0
    einastoConstructor%fourierProfileTableInitialized                    =.false.
    einastoConstructor%fourierProfileTableAlphaInterpolationReset        =.true.
    einastoConstructor%fourierProfileTableConcentrationInterpolationReset=.true.
    einastoConstructor%fourierProfileTableWavenumberInterpolationReset   =.true.
    ! Store a pointer to the provided scale object.
    einastoConstructor%scale => scale
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
         &         'einastoConstructor'                                                                                    , &
         &         'Einasto dark matter profile requires a dark matter profile component that supports gettable '//          &
         &         '"scale" and "shape" properties.'                                                             //          &
         &         Galacticus_Component_List(                                                                                &
         &                                   'darkMatterProfile'                                                           , &
         &                                    defaultDarkMatterProfileComponent%shapeAttributeMatch(requireGettable=.true.)  &
         &                                   .intersection.                                                                  &
         &                                    defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)  &
         &                                  )                                                                                &
         &        )
    return
  end function einastoConstructor
  
  subroutine einastoDestructor(self)
    !% Destructor for the {\tt einasto} dark matter halo profile class.
    use Numerical_Interpolation
    implicit none
    type(darkMatterProfileEinasto), intent(inout) :: self

    call Interpolate_Done(                                                                                        &
         &                interpolationObject     =self%angularMomentumTableRadiusInterpolationObject           , &
         &                interpolationAccelerator=self%angularMomentumTableRadiusInterpolationAccelerator      , &
         &                reset                   =self%angularMomentumTableRadiusInterpolationReset              &
         &               )
    call Interpolate_Done(                                                                                        &
         &                interpolationAccelerator=self%angularMomentumTableAlphaInterpolationAccelerator       , &
         &                reset                   =self%angularMomentumTableAlphaInterpolationReset               &
         &                )
    call Interpolate_Done(                                                                                        &
         &                interpolationObject     =self%energyTableConcentrationInterpolationObject             , &
         &                interpolationAccelerator=self%energyTableConcentrationInterpolationAccelerator        , &
         &                reset                   =self%energyTableConcentrationInterpolationReset                &
         &               )
    call Interpolate_Done(                                                                                        &
         &                interpolationObject     =self%energyTableAlphaInterpolationObject                     , &
         &                interpolationAccelerator=self%energyTableAlphaInterpolationAccelerator                , &
         &                reset                   =self%energyTableAlphaInterpolationReset                        &
         &                  )
    call Interpolate_Done(                                                                                        &
         &                interpolationObject     =self%fourierProfileTableWavenumberInterpolationObject        , &
         &                interpolationAccelerator=self%fourierProfileTableWavenumberInterpolationAccelerator   , &
         &                reset                   =self%fourierProfileTableWavenumberInterpolationReset           &
         &               )
    call Interpolate_Done(                                                                                        &
         &                interpolationAccelerator=self%fourierProfileTableAlphaInterpolationAccelerator        , &
         &                reset                   =self%fourierProfileTableAlphaInterpolationReset                &
         &                  )
    call Interpolate_Done(                                                                                        &
         &                interpolationAccelerator=self%fourierProfileTableConcentrationInterpolationAccelerator, &
         &                reset                   =self%fourierProfileTableConcentrationInterpolationReset        &
         &               )
    call Interpolate_Done(                                                                                        &
         &                interpolationAccelerator=self%freefallRadiusTableAlphaInterpolationAccelerator        , &
         &                reset                   =self%freefallRadiusTableAlphaInterpolationReset                &
         &               )
    call Interpolate_Done(                                                                                        &
         &                interpolationObject     =self%freefallRadiusTableRadiusInterpolationObject            , &
         &                interpolationAccelerator=self%freefallRadiusTableRadiusInterpolationAccelerator       , &
         &                reset                   =self%freefallRadiusTableRadiusInterpolationReset               &
         &               )
    if (self%scale%isFinalizable()) deallocate(self%scale)
    return
  end subroutine einastoDestructor
  
  subroutine einastoCalculationReset(self,thisNode)
    !% Reset the dark matter profile calculation.
    implicit none
    class(darkMatterProfileEinasto), intent(inout)          :: self
    type (treeNode                ), intent(inout), pointer :: thisNode

    call self%scale%calculationReset(thisNode)
    return
  end subroutine einastoCalculationReset

  double precision function einastoDensity(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\tt node} at the given {\tt radius} (given
    !% in units of Mpc).
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileEinasto      ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: radius
    class           (nodeComponentBasic            )               , pointer :: thisBasicComponent
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: alpha                         , radiusOverScaleRadius      , &
         &                                                                      scaleRadius                   , virialRadiusOverScaleRadius

    ! Get components.
    thisBasicComponent             => node%basic            (                 )
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                =thisDarkMatterProfileComponent%scale()
    alpha                      =thisDarkMatterProfileComponent%shape()
    radiusOverScaleRadius      =radius                       /scaleRadius
    virialRadiusOverScaleRadius=self%scale%virialRadius(node)/scaleRadius
    einastoDensity=self%densityScaleFree(radiusOverScaleRadius,virialRadiusOverScaleRadius,alpha)&
         &*thisBasicComponent%mass()/scaleRadius**3
    return
  end function einastoDensity

  double precision function einastoEnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\tt node} at the given {\tt radius} (given in
    !% units of Mpc).
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileEinasto      ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: radius
    class           (nodeComponentBasic            )               , pointer :: thisBasicComponent
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: alpha                         , radiusOverScaleRadius      , &
         &                                                                      scaleRadius                   , virialRadiusOverScaleRadius

    ! Get components.
    thisBasicComponent             => node%basic            (                 )
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                =thisDarkMatterProfileComponent%scale()
    alpha                      =thisDarkMatterProfileComponent%shape()
    radiusOverScaleRadius      =radius                       /scaleRadius
    virialRadiusOverScaleRadius=self%scale%virialRadius(node)/scaleRadius
    einastoEnclosedMass=self%enclosedMassScaleFree(radiusOverScaleRadius,virialRadiusOverScaleRadius,alpha)&
         &*thisBasicComponent%mass()
    return
  end function einastoEnclosedMass

  double precision function einastoCircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\tt node} at the given {\tt radius} (given in
    !% units of Mpc).
    use Numerical_Constants_Physical
    implicit none
    class           (darkMatterProfileEinasto), intent(inout)          :: self
    type            (treeNode                ), intent(inout), pointer :: node
    double precision                          , intent(in   )          :: radius

    if (radius > 0.0d0) then
       einastoCircularVelocity=sqrt(gravitationalConstantGalacticus&
            &*self%enclosedMass(node,radius)/radius)
    else
       einastoCircularVelocity=0.0d0
    end if
    return
  end function einastoCircularVelocity

  double precision function einastoCircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\tt node}.
    use Numerical_Constants_Physical
    use Root_Finder
    implicit none
    class           (darkMatterProfileEinasto      ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfile
    double precision                                , parameter              :: toleranceRelative    =1.0d-3
    double precision                                                         :: alpha                       , radiusScale, &
         &                                                                      radiusPeak
    type            (rootFinder                    )                         :: finder

    ! Get the shape parameter for this halo.
    thisDarkMatterProfile => node                 %darkMatterProfile(autoCreate=.true.)
    alpha                 =  thisDarkMatterProfile%shape            (                 )
    radiusScale           =  thisDarkMatterProfile%scale            (                 )
    ! Solve for the radius (in units of the scale radius) at which the rotation curve peaks.
    call finder%tolerance   (                                                       &
         &                   toleranceRelative  =toleranceRelative                  &
         &                  )
    call finder%rangeExpand (                                                       &
         &                   rangeExpandUpward  =2.0d0                            , &
         &                   rangeExpandDownward=0.5d0                            , &
         &                   rangeExpandType    =rangeExpandMultiplicative          &
         &                  )
    call finder%rootFunction(                                                       &
         &                                       einastoCircularVelocityPeakRadius  &
         &                  )
    radiusPeak=finder%find(rootGuess=radiusScale)
    ! Find the peak velocity.
    einastoCircularVelocityMaximum=self%circularVelocity(node,radiusPeak*radiusScale)
    return
    
  contains
    
    double precision function einastoCircularVelocityPeakRadius(radius)
      !% Computes the derivative of the square of circular velocity for an Einasto density profile.
      use Gamma_Functions
      implicit none
      double precision, intent(in   ) :: radius

      einastoCircularVelocityPeakRadius=                                         &
           & +        2.0d0               **(      +3.0d0/alpha)                 &
           & *        radius              **(-2.0d0+      alpha)                 &
           & *(       radius**alpha/alpha)**(-1.0d0+3.0d0/alpha)                 &
           & * exp(                                                              &
           &       -(                                                            &
           &         +2.0d0                                                      &
           &         *radius**alpha                                              &
           &        )                                                            &
           &       /alpha                                                        &
           &      )                                                              &
           & -Gamma_Function_Incomplete_Complementary(                           &
           &                                          3.0d0              /alpha, &
           &                                          2.0d0*radius**alpha/alpha  &
           &                                         )                           &
           & /       radius              **  2
     return
    end function einastoCircularVelocityPeakRadius

  end function einastoCircularVelocityMaximum

  double precision function einastoPotential(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\tt node} at the given {\tt radius} (given in
    !% units of Mpc).
    use Dark_Matter_Halo_Scales
    use Dark_Matter_Profiles_Error_Codes
    use Numerical_Constants_Physical
    implicit none
    class           (darkMatterProfileEinasto      ), intent(inout)           :: self
    type            (treeNode                      ), intent(inout), pointer  :: node
    double precision                                , intent(in   )           :: radius
    integer                                         , intent(  out), optional :: status
    class           (nodeComponentBasic            )               , pointer  :: thisBasicComponent
    class           (nodeComponentDarkMatterProfile)               , pointer  :: thisDarkMatterProfileComponent
    double precision                                                          :: alpha                         , radiusOverScaleRadius      , &
         &                                                                       scaleRadius                   , virialRadiusOverScaleRadius

    ! Assume success.
    if (present(status)) status=darkMatterProfileSuccess
    ! Get components.
    thisBasicComponent             => node%basic            (                 )
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                =thisDarkMatterProfileComponent%scale()
    alpha                      =thisDarkMatterProfileComponent%shape()
    radiusOverScaleRadius      =radius                       /scaleRadius
    virialRadiusOverScaleRadius=self%scale%virialRadius(node)/scaleRadius
    einastoPotential=                                                                                   &
         & ( self%potentialScaleFree(radiusOverScaleRadius      ,virialRadiusOverScaleRadius,alpha)     &
         &  -self%potentialScaleFree(virialRadiusOverScaleRadius,virialRadiusOverScaleRadius,alpha)     &
         &  -1.0d0/                                              virialRadiusOverScaleRadius            &
         & )                                                                                            &
         & *gravitationalConstantGalacticus*thisBasicComponent%mass()/scaleRadius
    return
  end function einastoPotential

  double precision function einastoRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\tt node} at which a circular orbit has the given {\tt specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc).
    use Numerical_Constants_Physical
    implicit none
    class           (darkMatterProfileEinasto      ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: specificAngularMomentum
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: alpha                           , scaleRadius, &
         &                                                                      specificAngularMomentumScaleFree

    ! Get components.
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)
    ! Get the scale radius of the halo.
    scaleRadius=thisDarkMatterProfileComponent%scale()
    ! Get the shape parameter of the halo.
    alpha      =thisDarkMatterProfileComponent%shape()
    ! Compute the specific angular momentum in scale free units.
    specificAngularMomentumScaleFree=specificAngularMomentum/sqrt(gravitationalConstantGalacticus*scaleRadius&
         &*self%enclosedMass(node,scaleRadius))
    ! Compute the corresponding radius.
    einastoRadiusFromSpecificAngularMomentum=scaleRadius*self%radiusFromSpecificAngularMomentumScaleFree(alpha&
         &,specificAngularMomentumScaleFree)
    return
  end function einastoRadiusFromSpecificAngularMomentum

  double precision function einastoRadiusFromSpecificAngularMomentumScaleFree(self,alpha,specificAngularMomentumScaleFree)
    !% Comptue the radius at which a circular orbit has the given {\tt specificAngularMomentumScaleFree} in a scale free Einasto
    !% profile.
    use Numerical_Interpolation
    implicit none
    class           (darkMatterProfileEinasto), intent(inout)  :: self
    double precision                          , intent(in   )  :: alpha , specificAngularMomentumScaleFree
    integer                                   , dimension(0:1) :: jAlpha
    double precision                          , dimension(0:1) :: hAlpha
    integer                                                    :: iAlpha

    ! Return immediately for zero angular momentum.
    if (specificAngularMomentumScaleFree <= 0.0d0) then
       einastoRadiusFromSpecificAngularMomentumScaleFree=0.0d0
       return
    end if

    ! Ensure the table exists and is sufficiently tabulated.
    call self%radiusFromSpecificAngularMomentumTableMake(alpha,specificAngularMomentumScaleFree)

    ! Get interpolating factors in alpha.
    jAlpha(0)=Interpolate_Locate(self%angularMomentumTableAlphaCount,self%angularMomentumTableAlpha&
         &,self%angularMomentumTableAlphaInterpolationAccelerator,alpha,reset=self%angularMomentumTableAlphaInterpolationReset)
    jAlpha(1)=jAlpha(0)+1
    hAlpha=Interpolate_Linear_Generate_Factors(self%angularMomentumTableAlphaCount,self%angularMomentumTableAlpha,jAlpha(0),alpha)

    ! Interpolate in specific angular momentum to get radius.
    einastoRadiusFromSpecificAngularMomentumScaleFree=0.0d0
    do iAlpha=0,1
       einastoRadiusFromSpecificAngularMomentumScaleFree=                      &
            &  einastoRadiusFromSpecificAngularMomentumScaleFree               &
            & +Interpolate( self%angularMomentumTableRadiusCount                    &
            &              ,self%angularMomentumTable(:,jAlpha(iAlpha))             &
            &              ,self%angularMomentumTableRadius                         &
            &              ,self%angularMomentumTableRadiusInterpolationObject      &
            &              ,self%angularMomentumTableRadiusInterpolationAccelerator &
            &              ,specificAngularMomentumScaleFree                   &
            &              ,reset=self%angularMomentumTableRadiusInterpolationReset &
            &             )                                                    &
            & *hAlpha(iAlpha)
    end do
    return
  end function einastoRadiusFromSpecificAngularMomentumScaleFree

  subroutine einastoRadiusFromSpecificAngularMomentumTableMake(self,alphaRequired,specificAngularMomentumRequired)
    !% Create a tabulation of the relation between specific angular momentum and radius in an Einasto profile.
    use Numerical_Interpolation
    use Numerical_Ranges
    use Gamma_Functions
    use Memory_Management
    implicit none
    class           (darkMatterProfileEinasto), intent(inout) :: self
    double precision                          , intent(in   ) :: alphaRequired, specificAngularMomentumRequired
    integer                                                   :: iAlpha       , iRadius
    logical                                                   :: makeTable
    double precision                                          :: alpha        , enclosedMass                   , &
         &                                                       radius

    ! Always check if we need to make the table.
    makeTable=.true.
    do while (makeTable)
       ! Assume table does not need remaking.
       makeTable=.false.
       ! Check for uninitialized table.
       if (.not.self%angularMomentumTableInitialized) then
          makeTable=.true.
          ! Check for alpha out of range.
       else if (alphaRequired < self%angularMomentumTableAlpha(1) .or. alphaRequired >&
            & self%angularMomentumTableAlpha(self%angularMomentumTableAlphaCount)) then
          makeTable=.true.
          ! Compute the range of tabulation and number of points to use.
          self%angularMomentumTableAlphaMinimum=min(self%angularMomentumTableAlphaMinimum,0.9d0*alphaRequired)
          self%angularMomentumTableAlphaMaximum=max(self%angularMomentumTableAlphaMaximum,1.1d0*alphaRequired)
          ! Check for angular momentum below minimum tabulated value.
       else if (any(specificAngularMomentumRequired < self%angularMomentumTable(1,:))) then
          makeTable=.true.
          self%angularMomentumTableRadiusMinimum=0.5d0*self%angularMomentumTableRadiusMinimum
          ! Check for angular momentum above maximum tabulated value.
       else if (any(specificAngularMomentumRequired > self%angularMomentumTable(self%angularMomentumTableRadiusCount,:))) then
          makeTable=.true.
          self%angularMomentumTableRadiusMaximum=2.0d0*self%angularMomentumTableRadiusMaximum
       end if
       ! Remake the table if necessary.
       if (makeTable) then
          ! Allocate arrays to the appropriate sizes.
          self%angularMomentumTableAlphaCount =int(      (self%angularMomentumTableAlphaMaximum -self%angularMomentumTableAlphaMinimum )&
               &*dble(einastoAngularMomentumTableAlphaPointsPerUnit   ))+1
          self%angularMomentumTableRadiusCount=int(log10(self%angularMomentumTableRadiusMaximum/self%angularMomentumTableRadiusMinimum)&
               &*dble(einastoAngularMomentumTableRadiusPointsPerDecade))+1
          if (allocated(self%angularMomentumTableAlpha )) call Dealloc_Array(self%angularMomentumTableAlpha )
          if (allocated(self%angularMomentumTableRadius)) call Dealloc_Array(self%angularMomentumTableRadius)
          if (allocated(self%angularMomentumTable      )) call Dealloc_Array(self%angularMomentumTable      )
          call Alloc_Array(self%angularMomentumTableAlpha ,[                                self%angularMomentumTableAlphaCount])
          call Alloc_Array(self%angularMomentumTableRadius,[self%angularMomentumTableRadiusCount                               ])
          call Alloc_Array(self%angularMomentumTable      ,[self%angularMomentumTableRadiusCount,self%angularMomentumTableAlphaCount])
          ! Create ranges of alpha and radius.
          self%angularMomentumTableAlpha =Make_Range(self%angularMomentumTableAlphaMinimum ,self%angularMomentumTableAlphaMaximum &
               &,self%angularMomentumTableAlphaCount ,rangeType=rangeTypeLinear     )
          self%angularMomentumTableRadius=Make_Range(self%angularMomentumTableRadiusMinimum,self%angularMomentumTableRadiusMaximum&
               &,self%angularMomentumTableRadiusCount,rangeType=rangeTypeLogarithmic)
          ! Tabulate the radius vs. specific angular momentum relation.
          do iAlpha=1,self%angularMomentumTableAlphaCount
             alpha=self%angularMomentumTableAlpha(iAlpha)
             do iRadius=1,self%angularMomentumTableRadiusCount
                radius=self%angularMomentumTableRadius(iRadius)
                enclosedMass= Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*radius**alpha/alpha) &
                     &       /Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0              /alpha)
                self%angularMomentumTable(iRadius,iAlpha)=sqrt(enclosedMass*radius)
             end do
          end do
          ! Reset interpolators.
          call Interpolate_Done(self%angularMomentumTableRadiusInterpolationObject,self%angularMomentumTableRadiusInterpolationAccelerator &
               &,self%angularMomentumTableRadiusInterpolationReset)
          call Interpolate_Done(interpolationAccelerator=self%angularMomentumTableAlphaInterpolationAccelerator &
               &,reset=self%angularMomentumTableAlphaInterpolationReset)
          self%angularMomentumTableRadiusInterpolationReset=.true.
          self%angularMomentumTableAlphaInterpolationReset =.true.
          ! Flag that the table is now initialized.
          self%angularMomentumTableInitialized=.true.
       end if
    end do
    return
  end subroutine einastoRadiusFromSpecificAngularMomentumTableMake

  double precision function einastoRotationNormalization(self,node)
    !% Return the rotation normalization of an Einasto halo density profile.
    use Numerical_Constants_Math
    use Gamma_Functions
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileEinasto      ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: alpha                         , scaleRadius, &
         &                                                                      virialRadiusOverScaleRadius

    ! Get components.
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)

    ! Get scale radius, shape and concentration.
    scaleRadius                =thisDarkMatterProfileComponent%scale()
    alpha                      =thisDarkMatterProfileComponent%shape()
    virialRadiusOverScaleRadius=self%scale%virialRadius(node)/scaleRadius

    ! Compute the rotation normalization.
    einastoRotationNormalization=                                                                               &
         &  (4.0d0/Pi/scaleRadius)                                                                              &
         & *(2.0d0/alpha)**(1.0d0/alpha)                                                                        &
         & *Gamma_Function                         (3.0d0/alpha                                               ) &
         & /Gamma_Function                         (4.0d0/alpha                                               ) &
         & *Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*virialRadiusOverScaleRadius**alpha/alpha) &
         & /Gamma_Function_Incomplete_Complementary(4.0d0/alpha,2.0d0*virialRadiusOverScaleRadius**alpha/alpha)

    return
  end function einastoRotationNormalization

  double precision function einastoEnergy(self,node)
    !% Return the energy of an Einasto halo density profile.
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    implicit none
    class           (darkMatterProfileEinasto      ), intent(inout)           :: self
    type            (treeNode                      ), intent(inout) , pointer :: node
    class           (nodeComponentBasic            )                , pointer :: thisBasicComponent
    class           (nodeComponentDarkMatterProfile)                , pointer :: thisDarkMatterProfileComponent
    integer                                         , dimension(0:1)          :: jAlpha
    double precision                                , dimension(0:1)          :: hAlpha
    integer                                                                   :: iAlpha
    double precision                                                          :: alpha                         , scaleRadius, &
         &                                                                       virialRadiusOverScaleRadius

    ! Get components.
    thisBasicComponent             => node%basic            (                 )
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)

    ! Get scale radius, shape parameter and concentration.
    scaleRadius                =thisDarkMatterProfileComponent%scale()
    alpha                      =thisDarkMatterProfileComponent%shape()
    virialRadiusOverScaleRadius=self%scale%virialRadius(node)/scaleRadius

    ! Ensure the table exists and is sufficiently tabulated.
    call self%energyTableMake(virialRadiusOverScaleRadius,alpha)

    ! Get interpolating factors in alpha.
    jAlpha(0)=Interpolate_Locate(self%energyTableAlphaCount,self%energyTableAlpha,self%energyTableAlphaInterpolationAccelerator,alpha,reset&
         &=self%energyTableAlphaInterpolationReset)
    jAlpha(1)=jAlpha(0)+1
    hAlpha=Interpolate_Linear_Generate_Factors(self%energyTableAlphaCount,self%energyTableAlpha,jAlpha(0),alpha)

    ! Find the energy by interpolation.
    einastoEnergy=0.0d0
    do iAlpha=0,1
       einastoEnergy=einastoEnergy+Interpolate(self%energyTableConcentrationCount&
            &,self%energyTableConcentration,self%energyTable(:,jAlpha(iAlpha)),self%energyTableConcentrationInterpolationObject &
            &,self%energyTableConcentrationInterpolationAccelerator,virialRadiusOverScaleRadius,reset&
            &=self%energyTableConcentrationInterpolationReset)*hAlpha(iAlpha)
    end do

    ! Scale to dimensionful units.
    einastoEnergy=einastoEnergy*thisBasicComponent%mass() &
         &*self%scale%virialVelocity(node)**2
    return
  end function einastoEnergy

  double precision function einastoEnergyGrowthRate(self,node)
    !% Return the energy of an Einasto halo density profile.
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    implicit none
    class           (darkMatterProfileEinasto      ), intent(inout)           :: self
    type            (treeNode                      ), intent(inout) , pointer :: node
    class           (nodeComponentBasic            )                , pointer :: thisBasicComponent
    class           (nodeComponentDarkMatterProfile)                , pointer :: thisDarkMatterProfileComponent
    integer                                         , dimension(0:1)          :: jAlpha
    double precision                                , dimension(0:1)          :: hAlpha
    integer                                                                   :: iAlpha
    double precision                                                          :: alpha                         , energy     , &
         &                                                                       energyGradient                , scaleRadius, &
         &                                                                       virialRadiusOverScaleRadius

    ! Get components.
    thisBasicComponent             => node%basic            (                 )
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)

    ! Get scale radius, shape parameter and concentration.
    scaleRadius                =thisDarkMatterProfileComponent%scale()
    alpha                      =thisDarkMatterProfileComponent%shape()
    virialRadiusOverScaleRadius=self%scale%virialRadius(node)/scaleRadius

    ! Ensure the table exists and is sufficiently tabulated.
    call self%energyTableMake(virialRadiusOverScaleRadius,alpha)

    ! Get interpolating factors in alpha.
    jAlpha(0)=Interpolate_Locate(self%energyTableAlphaCount,self%energyTableAlpha,self%energyTableAlphaInterpolationAccelerator,alpha,reset&
         &=self%energyTableAlphaInterpolationReset)
    jAlpha(1)=jAlpha(0)+1
    hAlpha=Interpolate_Linear_Generate_Factors(self%energyTableAlphaCount,self%energyTableAlpha,jAlpha(0),alpha)

    ! Find the energy gradient by interpolation.
    energy        =0.0d0
    energyGradient=0.0d0
    do iAlpha=0,1
       energy        =energy        +Interpolate           (self%energyTableConcentrationCount,self%energyTableConcentration,self%energyTable(:&
            &,jAlpha(iAlpha)),self%energyTableConcentrationInterpolationObject ,self%energyTableConcentrationInterpolationAccelerator&
            &,virialRadiusOverScaleRadius,reset=self%energyTableConcentrationInterpolationReset)*hAlpha(iAlpha)
       energyGradient=energyGradient+Interpolate_Derivative(self%energyTableConcentrationCount,self%energyTableConcentration,self%energyTable(:&
            &,jAlpha(iAlpha)),self%energyTableConcentrationInterpolationObject ,self%energyTableConcentrationInterpolationAccelerator&
            &,virialRadiusOverScaleRadius,reset=self%energyTableConcentrationInterpolationReset)*hAlpha(iAlpha)
    end do

    ! Compute the energy growth rate.
    einastoEnergyGrowthRate=self%energy(node)&
         &*(thisBasicComponent%accretionRate()/thisBasicComponent%mass()+2.0d0 &
         &*self%scale%virialVelocityGrowthRate(node)/self%scale%virialVelocity(node)+(energyGradient&
         &*virialRadiusOverScaleRadius/energy)*(self%scale%virialRadiusGrowthRate(node)&
         &/self%scale%virialRadius(node)-thisDarkMatterProfileComponent%scaleGrowthRate()&
         &/thisDarkMatterProfileComponent%scale()))

    return
  end function einastoEnergyGrowthRate

  subroutine einastoEnergyTableMake(self,concentrationRequired,alphaRequired)
    !% Create a tabulation of the energy of Einasto profiles as a function of their concentration of $\alpha$ parameter.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Interpolation
    use Numerical_Integration
    use Numerical_Ranges
    use Memory_Management
    use Numerical_Constants_Math
    implicit none
    class           (darkMatterProfileEinasto  ), intent(inout) :: self
    double precision                            , intent(in   ) :: alphaRequired          , concentrationRequired
    integer                                                     :: iAlpha                 , iConcentration
    logical                                                     :: makeTable
    double precision                                            :: alpha                  , concentration        , &
         &                                                         jeansEquationIntegral  , kineticEnergy        , &
         &                                                         kineticEnergyIntegral  , potentialEnergy      , &
         &                                                         potentialEnergyIntegral, radiusMaximum        , &
         &                                                         radiusMinimum          , concentrationParameter, &
         &                                                         alphaParameter
    type            (c_ptr                     )                :: parameterPointer
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace

    ! Assume table does not need remaking.
    makeTable=.false.
    ! Check for uninitialized table.
    if (.not.self%energyTableInitialized) makeTable=.true.
    ! Check for alpha out of range.
    if (alphaRequired < self%energyTableAlphaMinimum .or. alphaRequired > self%energyTableAlphaMaximum) then
       makeTable=.true.
       ! Compute the range of tabulation and number of points to use.
       self%energyTableAlphaMinimum =min(self%energyTableAlphaMinimum,0.9d0*alphaRequired)
       self%energyTableAlphaMaximum =max(self%energyTableAlphaMaximum,1.1d0*alphaRequired)
    end if
    ! Check for concentration below minimum tabulated value.
    if (concentrationRequired < self%energyTableConcentrationMinimum .or. concentrationRequired > self%energyTableConcentrationMaximum) then
       makeTable=.true.
       self%energyTableConcentrationMinimum=min(self%energyTableConcentrationMinimum,0.5d0*self%energyTableConcentrationMinimum)
       self%energyTableConcentrationMaximum=max(self%energyTableConcentrationMaximum,2.0d0*self%energyTableConcentrationMaximum)
    end if
    ! Remake the table if necessary.
    if (makeTable) then
       ! Allocate arrays to the appropriate sizes.
       self%energyTableAlphaCount        =int(      (self%energyTableAlphaMaximum        -self%energyTableAlphaMinimum        ) &
            &*dble(einastoEnergyTableAlphaPointsPerUnit          ))+1
       self%energyTableConcentrationCount=int(log10(self%energyTableConcentrationMaximum/self%energyTableConcentrationMinimum) &
            &*dble(einastoEnergyTableConcentrationPointsPerDecade))+1
       if (allocated(self%energyTableAlpha        )) call Dealloc_Array(self%energyTableAlpha        )
       if (allocated(self%energyTableConcentration)) call Dealloc_Array(self%energyTableConcentration)
       if (allocated(self%energyTable             )) call Dealloc_Array(self%energyTable             )
       call Alloc_Array(self%energyTableAlpha        ,[                              self%energyTableAlphaCount])
       call Alloc_Array(self%energyTableConcentration,[self%energyTableConcentrationCount                      ])
       call Alloc_Array(self%energyTable             ,[self%energyTableConcentrationCount,self%energyTableAlphaCount])
       ! Create ranges of alpha and concentration.
       self%energyTableAlpha        =Make_Range(self%energyTableAlphaMinimum        ,self%energyTableAlphaMaximum         &
            &,self%energyTableAlphaCount        ,rangeType=rangeTypeLinear     )
       self%energyTableConcentration=Make_Range(self%energyTableConcentrationMinimum,self%energyTableConcentrationMaximum &
            &,self%energyTableConcentrationCount,rangeType=rangeTypeLogarithmic)
       ! Tabulate the radius vs. specific angular momentum relation.
       do iAlpha=1,self%energyTableAlphaCount
          alpha=self%energyTableAlpha(iAlpha)
          do iConcentration=1,self%energyTableConcentrationCount
             concentration=self%energyTableConcentration(iConcentration)

             ! Compute the potential energy.
             radiusMinimum         =0.0d0
             radiusMaximum         =concentration
             concentrationParameter=concentration
             alphaParameter        =alpha
             potentialEnergyIntegral=Integrate(radiusMinimum,radiusMaximum,einastoPotentialEnergyIntegrand&
                  &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
             call Integrate_Done(integrandFunction,integrationWorkspace)
             potentialEnergy=-0.5d0*(1.0d0/concentration+potentialEnergyIntegral)

             ! Compute the velocity dispersion at the virial radius.
             radiusMinimum         =        concentration
             radiusMaximum         =100.0d0*concentration
             concentrationParameter=        concentration
             alphaParameter        =alpha
             jeansEquationIntegral=Integrate(radiusMinimum,radiusMaximum,einastoJeansEquationIntegrand&
                  &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
             call Integrate_Done(integrandFunction,integrationWorkspace)

             ! Compute the kinetic energy.
             radiusMinimum         =0.0d0
             radiusMaximum         =concentration
             concentrationParameter=concentration
             alphaParameter        =alpha
             kineticEnergyIntegral=Integrate(radiusMinimum,radiusMaximum,einastoKineticEnergyIntegrand&
                  &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
             call Integrate_Done(integrandFunction,integrationWorkspace)
             kineticEnergy=2.0d0*Pi*(jeansEquationIntegral*concentration**3+kineticEnergyIntegral)

             ! Compute the total energy.
             self%energyTable(iConcentration,iAlpha)=(potentialEnergy+kineticEnergy)*concentration

          end do
       end do
       ! Reset interpolators.
       call Interpolate_Done(self%energyTableConcentrationInterpolationObject,self%energyTableConcentrationInterpolationAccelerator &
            &,self%energyTableConcentrationInterpolationReset)
       call Interpolate_Done(self%energyTableAlphaInterpolationObject        ,self%energyTableAlphaInterpolationAccelerator         &
            &,self%energyTableAlphaInterpolationReset        )
       self%energyTableConcentrationInterpolationReset=.true.
       self%energyTableAlphaInterpolationReset        =.true.
       ! Flag that the table is now initialized.
       self%energyTableInitialized=.true.
    end if
    return
    
  contains
    
    function einastoPotentialEnergyIntegrand(radius,parameterPointer) bind(c)
      !% Integrand for Einasto profile potential energy.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real(kind=c_double)        :: einastoPotentialEnergyIntegrand
      real(kind=c_double), value :: radius
      type(c_ptr        ), value :: parameterPointer
      
      einastoPotentialEnergyIntegrand=(self%enclosedMassScaleFree(radius,concentrationParameter,alphaParameter)/radius)**2
      return
    end function einastoPotentialEnergyIntegrand

    function einastoKineticEnergyIntegrand(radius,parameterPointer) bind(c)
      !% Integrand for Einasto profile kinetic energy.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real(kind=c_double)        :: einastoKineticEnergyIntegrand
      real(kind=c_double), value :: radius
      type(c_ptr        ), value :: parameterPointer
      
      einastoKineticEnergyIntegrand=self%enclosedMassScaleFree(radius,concentrationParameter,alphaParameter)&
           &*self%densityScaleFree(radius,concentrationParameter,alphaParameter)*radius
      return
    end function einastoKineticEnergyIntegrand
    
    function einastoJeansEquationIntegrand(radius,parameterPointer) bind(c)
      !% Integrand for Einasto profile Jeans equation.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real(kind=c_double)        :: einastoJeansEquationIntegrand
      real(kind=c_double), value :: radius
      type(c_ptr        ), value :: parameterPointer
      
      einastoJeansEquationIntegrand=self%enclosedMassScaleFree(radius,concentrationParameter,alphaParameter)&
           &*self%densityScaleFree(radius ,concentrationParameter,alphaParameter)/radius**2
      return
    end function einastoJeansEquationIntegrand
    
  end subroutine einastoEnergyTableMake

  double precision function einastoEnclosedMassScaleFree(self,radius,concentration,alpha)
    !% Returns the enclosed mass (in units of the virial mass) in an Einasto dark matter profile with given {\tt concentration} at the
    !% given {\tt radius} (given in units of the scale radius).
    use Gamma_Functions
    implicit none
    class           (darkMatterProfileEinasto), intent(inout) :: self
    double precision                          , intent(in   ) :: alpha, concentration, radius

    if (radius >= concentration) then
       einastoEnclosedMassScaleFree=1.0d0
    else
       einastoEnclosedMassScaleFree=                                                                 &
            &  Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*radius       **alpha/alpha) &
            & /Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*concentration**alpha/alpha)
    end if
    return
  end function einastoEnclosedMassScaleFree

  double precision function einastoDensityScaleFree(self,radius,concentration,alpha)
    !% Returns the density (in units such that the virial mass and scale length are unity) in an Einasto dark matter profile with
    !% given {\tt concentration} and {\tt alpha} at the given {\tt radius} (given in units of the scale radius).
    use Numerical_Constants_Math
    use Gamma_Functions
    implicit none
    class           (darkMatterProfileEinasto), intent(inout) :: self
    double precision                          , intent(in   ) :: alpha               , concentration, radius
    double precision                                          :: densityNormalization

    densityNormalization= (alpha/4.0d0/Pi)                                                                      &
         &               *    ((2.0d0/alpha)                   **(3.0d0/alpha)                                ) &
         &               *exp(-2.0d0/alpha                                                                   ) &
         &               /Gamma_Function                         (3.0d0/alpha                                 ) &
         &               /Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*concentration**alpha/alpha)
    einastoDensityScaleFree=densityNormalization*exp(-(2.0d0/alpha)*(radius**alpha-1.0d0))
    return
  end function einastoDensityScaleFree

  double precision function einastoPotentialScaleFree(self,radius,concentration,alpha)
    !% Returns the gravitational potential (in units where the virial mass and scale radius are unity) in an Einasto dark matter
    !% profile with given {\tt concentration} and {\tt alpha} at the given {\tt radius} (given in units of the scale radius).
    use Gamma_Functions
    implicit none
    class           (darkMatterProfileEinasto), intent(inout) :: self
    double precision                          , intent(in   ) :: alpha, concentration, radius

    if (radius <= 0.0d0) then
       einastoPotentialScaleFree=                                                                    &
            & -((2.0d0/alpha)**(1.0d0/alpha))                                                        &
            & *2.0d0                                                                                 &
            & *Gamma_Function                         (2.0d0/alpha                                 ) &
            & /Gamma_Function                         (3.0d0/alpha                                 ) &
            & /Gamma_Function_Incomplete_Complementary(3.0d0/alpha,2.0d0*concentration**alpha/alpha)
       
    else
       einastoPotentialScaleFree=                                                                                    &
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
  end function einastoPotentialScaleFree

  double precision function einastoKSpace(self,node,wavenumber)
    !% Returns the Fourier transform of the Einasto density profile at the specified {\tt waveNumber} (given in Mpc$^{-1}$)).
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    implicit none
    class           (darkMatterProfileEinasto      )                , intent(inout)          :: self
    type            (treeNode                      )                , intent(inout), pointer :: node
    double precision                                                , intent(in   )          :: wavenumber
    class           (nodeComponentDarkMatterProfile)                               , pointer :: thisDarkMatterProfileComponent
    integer                                         , dimension(0:1)                         :: jAlpha                        , jConcentration
    double precision                                , dimension(0:1)                         :: hAlpha                        , hConcentration
    integer                                                                                  :: iAlpha                        , iConcentration
    double precision                                                                         :: alpha                         , scaleRadius        , &
         &                                                                                      virialRadiusOverScaleRadius   , wavenumberScaleFree

    ! Get components.
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)

    ! Get scale radius, shape parameter and concentration.
    scaleRadius                =thisDarkMatterProfileComponent%scale()
    alpha                      =thisDarkMatterProfileComponent%shape()
    virialRadiusOverScaleRadius=self%scale%virialRadius(node)/scaleRadius
    wavenumberScaleFree        =wavenumber*scaleRadius

    ! Ensure the table exists and is sufficiently tabulated.
    call self%fourierProfileTableMake(wavenumberScaleFree,virialRadiusOverScaleRadius,alpha)

    ! Get interpolating factors in alpha.
    jAlpha(0)=Interpolate_Locate(self%fourierProfileTableAlphaCount,self%fourierProfileTableAlpha&
         &,self%fourierProfileTableAlphaInterpolationAccelerator,alpha,reset=self%fourierProfileTableAlphaInterpolationReset)
    jAlpha(1)=jAlpha(0)+1
    hAlpha=Interpolate_Linear_Generate_Factors(self%fourierProfileTableAlphaCount,self%fourierProfileTableAlpha,jAlpha(0),alpha)

    ! Get interpolating factors in concentration.
    jConcentration(0)=Interpolate_Locate(self%fourierProfileTableConcentrationCount,self%fourierProfileTableConcentration &
         &,self%fourierProfileTableConcentrationInterpolationAccelerator,virialRadiusOverScaleRadius,reset &
         &=self%fourierProfileTableConcentrationInterpolationReset)
    jConcentration(1)=jConcentration(0)+1
    hConcentration=Interpolate_Linear_Generate_Factors(self%fourierProfileTableConcentrationCount,self%fourierProfileTableConcentration&
         &,jConcentration(0),virialRadiusOverScaleRadius)

    ! Find the Fourier profile by interpolation.
    einastoKSpace=0.0d0
    do iAlpha=0,1
       do iConcentration=0,1
          einastoKSpace=einastoKSpace+Interpolate(self%fourierProfileTableWavenumberCount &
               &,self%fourierProfileTableWavenumber,self%fourierProfileTable(:,jConcentration(iConcentration),jAlpha(iAlpha))&
               &,self%fourierProfileTableWavenumberInterpolationObject ,self%fourierProfileTableWavenumberInterpolationAccelerator&
               &,wavenumberScaleFree,reset =self%fourierProfileTableWavenumberInterpolationReset)*hAlpha(iAlpha)*hConcentration(iConcentration)
       end do
    end do
    return
  end function einastoKSpace

  subroutine einastoFourierProfileTableMake(self,wavenumberRequired,concentrationRequired,alphaRequired)
    !% Create a tabulation of the Fourier transform of Einasto profiles as a function of their $\alpha$ parameter and
    !% dimensionless wavenumber.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Interpolation
    use Numerical_Integration
    use Numerical_Ranges
    use Memory_Management
    use Galacticus_Display
    implicit none
    class           (darkMatterProfileEinasto  ), intent(inout) :: self
    double precision                            , intent(in   ) :: alphaRequired       , concentrationRequired, wavenumberRequired
    integer                                                     :: iAlpha              , iConcentration       , iWavenumber       , &
         &                                                         percentage
    logical                                                     :: makeTable
    double precision                                            :: alpha               , concentration        , radiusMaximum     , &
         &                                                         radiusMinimum       , wavenumber           , wavenumberParameter, &
         &                                                         concentrationParameter, alphaParameter
    type            (c_ptr                     )                :: parameterPointer
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace

    ! Assume table does not need remaking.
    makeTable=.false.
    ! Check for uninitialized table.
    if (.not.self%fourierProfileTableInitialized) makeTable=.true.
    ! Check for alpha out of range.
    if (alphaRequired         < self%fourierProfileTableAlphaMinimum         .or. alphaRequired         > self%fourierProfileTableAlphaMaximum        ) then
       makeTable=.true.
       ! Compute the range of tabulation and number of points to use.
       self%fourierProfileTableAlphaMinimum        =min(self%fourierProfileTableAlphaMinimum     ,0.9d0*alphaRequired           )
       self%fourierProfileTableAlphaMaximum        =max(self%fourierProfileTableAlphaMaximum     ,1.1d0*alphaRequired           )
    end if
    ! Check for concentration out of range.
    if (concentrationRequired < self%fourierProfileTableConcentrationMinimum .or. concentrationRequired > self%fourierProfileTableConcentrationMaximum ) then
       makeTable=.true.
       ! Compute the range of tabulation and number of points to use.
       self%fourierProfileTableConcentrationMinimum=min(self%fourierProfileTableConcentrationMinimum,0.5d0*concentrationRequired)
       self%fourierProfileTableConcentrationMaximum=max(self%fourierProfileTableConcentrationMaximum,2.0d0*concentrationRequired)
    end if
    ! Check for wavenumber below minimum tabulated value.
    if (wavenumberRequired    < self%fourierProfileTableWavenumberMinimum    .or. wavenumberRequired    > self%fourierProfileTableWavenumberMaximum    ) then
       makeTable=.true.
       self%fourierProfileTableWavenumberMinimum   =min(self%fourierProfileTableWavenumberMinimum   ,0.5d0*wavenumberRequired   )
       self%fourierProfileTableWavenumberMaximum   =max(self%fourierProfileTableWavenumberMaximum   ,2.0d0*wavenumberRequired   )
    end if
    ! Remake the table if necessary.
    if (makeTable) then
       ! Display a message.
       call Galacticus_Display_Indent('Constructing Einasto profile Fourier transform lookup table...',verbosityInfo)
       ! Allocate arrays to the appropriate sizes.
       self%fourierProfileTableAlphaCount        =int(      (self%fourierProfileTableAlphaMaximum        -self%fourierProfileTableAlphaMinimum        ) &
            &*dble(einastoFourierProfileTableAlphaPointsPerUnit          ))+1
       self%fourierProfileTableConcentrationCount=int(log10(self%fourierProfileTableConcentrationMaximum/self%fourierProfileTableConcentrationMinimum) &
            &*dble(einastoFourierProfileTableConcentrationPointsPerDecade))+1
       self%fourierProfileTableWavenumberCount   =int(log10(self%fourierProfileTableWavenumberMaximum   /self%fourierProfileTableWavenumberMinimum   ) &
            &*dble(einastoFourierProfileTableWavenumberPointsPerDecade   ))+1
       if (allocated(self%fourierProfileTableAlpha        )) call Dealloc_Array(self%fourierProfileTableAlpha        )
       if (allocated(self%fourierProfileTableConcentration)) call Dealloc_Array(self%fourierProfileTableConcentration)
       if (allocated(self%fourierProfileTableWavenumber   )) call Dealloc_Array(self%fourierProfileTableWavenumber   )
       if (allocated(self%fourierProfileTable             )) call Dealloc_Array(self%fourierProfileTable             )
       call Alloc_Array(self%fourierProfileTableAlpha        ,[                                                                         self%fourierProfileTableAlphaCount])
       call Alloc_Array(self%fourierProfileTableConcentration,[                                   self%fourierProfileTableConcentrationCount                              ])
       call Alloc_Array(self%fourierProfileTableWavenumber   ,[self%fourierProfileTableWavenumberCount                                                                    ])
       call Alloc_Array(self%fourierProfileTable             ,[self%fourierProfileTableWavenumberCount,self%fourierProfileTableConcentrationCount,self%fourierProfileTableAlphaCount])
       ! Create ranges of alpha and wavenumber.
       self%fourierProfileTableAlpha        =Make_Range(self%fourierProfileTableAlphaMinimum        ,self%fourierProfileTableAlphaMaximum         &
            &,self%fourierProfileTableAlphaCount        ,rangeType=rangeTypeLinear     )
       self%fourierProfileTableConcentration=Make_Range(self%fourierProfileTableConcentrationMinimum,self%fourierProfileTableConcentrationMaximum &
            &,self%fourierProfileTableConcentrationCount,rangeType=rangeTypeLogarithmic)
       self%fourierProfileTableWavenumber   =Make_Range(self%fourierProfileTableWavenumberMinimum   ,self%fourierProfileTableWavenumberMaximum    &
            &,self%fourierProfileTableWavenumberCount   ,rangeType=rangeTypeLogarithmic)
       ! Tabulate the Fourier profile.
       do iAlpha=1,self%fourierProfileTableAlphaCount
          alpha=self%fourierProfileTableAlpha(iAlpha)
          do iConcentration=1,self%fourierProfileTableConcentrationCount
             concentration=self%fourierProfileTableConcentration(iConcentration)

             ! Show progress.
             percentage=int(100.0d0*dble((iAlpha-1)*self%fourierProfileTableConcentrationCount+iConcentration-1)&
                  &/dble(self%fourierProfileTableAlphaCount*self%fourierProfileTableConcentrationCount))
             call Galacticus_Display_Counter(percentage,iAlpha == 1 .and. iConcentration == 1,verbosityWorking)

             do iWavenumber=1,self%fourierProfileTableWavenumberCount
                wavenumber=self%fourierProfileTableWavenumber(iWavenumber)
                ! Compute the potential fourierProfile.
                radiusMinimum         =0.0d0
                radiusMaximum         =concentration
                wavenumberParameter   =wavenumber
                alphaParameter        =alpha
                concentrationParameter=concentration
                self%fourierProfileTable(iWavenumber,iConcentration,iAlpha)=Integrate(radiusMinimum,radiusMaximum&
                     &,einastoFourierProfileIntegrand ,parameterPointer,integrandFunction,integrationWorkspace&
                     &,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-2,maxIntervals=10000)
                call Integrate_Done(integrandFunction,integrationWorkspace)

             end do
          end do
       end do
       call Galacticus_Display_Counter_Clear(verbosityWorking)
       ! Reset interpolators.
       call Interpolate_Done(self%fourierProfileTableWavenumberInterpolationObject,self%fourierProfileTableWavenumberInterpolationAccelerator &
            &,self%fourierProfileTableWavenumberInterpolationReset)
       call Interpolate_Done(interpolationAccelerator=self%fourierProfileTableAlphaInterpolationAccelerator         &
            &,reset=self%fourierProfileTableAlphaInterpolationReset        )
       call Interpolate_Done(interpolationAccelerator=self%fourierProfileTableConcentrationInterpolationAccelerator &
            &,reset=self%fourierProfileTableConcentrationInterpolationReset)
       self%fourierProfileTableWavenumberInterpolationReset   =.true.
       self%fourierProfileTableAlphaInterpolationReset        =.true.
       self%fourierProfileTableConcentrationInterpolationReset=.true.
       ! Flag that the table is now initialized.
       self%fourierProfileTableInitialized=.true.
       ! Display a message.
       call Galacticus_Display_Unindent('done',verbosityInfo)
    end if
    return
    
  contains

    function einastoFourierProfileIntegrand(radius,parameterPointer) bind(c)
      !% Integrand for Einasto Fourier profile.
      use, intrinsic :: ISO_C_Binding
      use Numerical_Constants_Math
      implicit none
      real(kind=c_double)        :: einastoFourierProfileIntegrand
      real(kind=c_double), value :: radius
      type(c_ptr        ), value :: parameterPointer
      
      einastoFourierProfileIntegrand=4.0d0*Pi*radius*sin(wavenumberParameter*radius)*self%densityScaleFree(radius&
           &,concentrationParameter,alphaParameter)/wavenumberParameter
      return
    end function einastoFourierProfileIntegrand
    
  end subroutine einastoFourierProfileTableMake

  double precision function einastoFreefallRadius(self,node,time)
    !% Returns the freefall radius in the Einasto density profile at the specified {\tt time} (given in Gyr).
    use Numerical_Interpolation
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    implicit none
    class           (darkMatterProfileEinasto      )                , intent(inout)          :: self
    type            (treeNode                      )                , intent(inout), pointer :: node
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
       einastoFreefallRadius=0.0d0
       return
    end if

    ! Get components.
    thisBasicComponent             => node%basic            (                 )
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)

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
    call self%freefallTabulate(freefallTimeScaleFree,alpha)

    ! Interpolate to get the freefall radius.
    ! Get interpolating factors in alpha.
    jAlpha(0)=Interpolate_Locate(self%freefallRadiusTableAlphaCount,self%freefallRadiusTableAlpha&
         &,self%freefallRadiusTableAlphaInterpolationAccelerator,alpha,reset=self%freefallRadiusTableAlphaInterpolationReset)
    jAlpha(1)=jAlpha(0)+1
    hAlpha=Interpolate_Linear_Generate_Factors(self%freefallRadiusTableAlphaCount,self%freefallRadiusTableAlpha,jAlpha(0),alpha)

    einastoFreefallRadius=0.0d0
    do iAlpha=0,1
       einastoFreefallRadius=einastoFreefallRadius+Interpolate(self%freefallRadiusTableRadiusCount,self%freefallRadiusTable(:&
            &,jAlpha(iAlpha)),self%freefallRadiusTableRadius ,self%freefallRadiusTableRadiusInterpolationObject&
            &,self%freefallRadiusTableRadiusInterpolationAccelerator ,freefallTimeScaleFree,reset&
            &=self%freefallRadiusTableRadiusInterpolationReset)*hAlpha(iAlpha)
    end do
    einastoFreefallRadius=einastoFreefallRadius*radiusScale
    return
  end function einastoFreefallRadius

  double precision function einastoFreefallRadiusIncreaseRate(self,node,time)
    !% Returns the rate of increase of the freefall radius in the Einasto density profile at the specified {\tt time} (given in
    !% Gyr).
    use Numerical_Interpolation
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    implicit none
    class           (darkMatterProfileEinasto      )                , intent(inout)          :: self
    type            (treeNode                      )                , intent(inout), pointer :: node
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
       einastoFreefallRadiusIncreaseRate=0.0d0
       return
    end if

    ! Get components.
    thisBasicComponent             => node%basic            (                 )
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)

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
    call self%freefallTabulate(freefallTimeScaleFree,alpha)

    ! Interpolate to get the freefall radius.
    ! Get interpolating factors in alpha.
    jAlpha(0)=Interpolate_Locate(self%freefallRadiusTableAlphaCount,self%freefallRadiusTableAlpha&
         &,self%freefallRadiusTableAlphaInterpolationAccelerator,alpha,reset=self%freefallRadiusTableAlphaInterpolationReset)
    jAlpha(1)=jAlpha(0)+1
    hAlpha=Interpolate_Linear_Generate_Factors(self%freefallRadiusTableAlphaCount,self%freefallRadiusTableAlpha,jAlpha(0),alpha)

    einastoFreefallRadiusIncreaseRate=0.0d0
    do iAlpha=0,1
       einastoFreefallRadiusIncreaseRate=einastoFreefallRadiusIncreaseRate+Interpolate_Derivative(self%freefallRadiusTableRadiusCount,self%freefallRadiusTable(:&
            &,jAlpha(iAlpha)),self%freefallRadiusTableRadius ,self%freefallRadiusTableRadiusInterpolationObject&
            &,self%freefallRadiusTableRadiusInterpolationAccelerator ,freefallTimeScaleFree,reset&
            &=self%freefallRadiusTableRadiusInterpolationReset)*hAlpha(iAlpha)
    end do
    einastoFreefallRadiusIncreaseRate=einastoFreefallRadiusIncreaseRate*radiusScale/timeScale
    return
  end function einastoFreefallRadiusIncreaseRate

  subroutine einastoFreefallTabulate(self,freefallTimeScaleFree,alphaRequired)
    !% Tabulates the freefall time vs. freefall radius for Einasto halos.
    use Galacticus_Display
    use Numerical_Ranges
    use Memory_Management
    use Numerical_Interpolation
    implicit none
    class           (darkMatterProfileEinasto), intent(inout) :: self
    double precision                          , intent(in   ) :: alphaRequired, freefallTimeScaleFree
    logical                                                   :: retabulate
    integer                                                   :: iAlpha       , iRadius              , percentage
    double precision                                          :: alpha

    retabulate=.not.self%freefallRadiusTableInitialized
    ! If the table has not yet been made, compute and store the freefall times corresponding to the minimum and maximum
    ! radii that will be tabulated by default.
    if (retabulate) then
       self%freefallTimeMinimum=self%freefallTimeScaleFree(self%freefallRadiusTableRadiusMinimum,alphaRequired)
       self%freefallTimeMaximum=self%freefallTimeScaleFree(self%freefallRadiusTableRadiusMaximum,alphaRequired)
    end if
    do while (freefallTimeScaleFree < self%freefallTimeMinimum)
       self%freefallRadiusTableRadiusMinimum=0.5d0*self%freefallRadiusTableRadiusMinimum
       self%freefallTimeMinimum=self%freefallTimeScaleFree(self%freefallRadiusTableRadiusMinimum,alphaRequired)
       retabulate=.true.
    end do
    do while (freefallTimeScaleFree > self%freefallTimeMaximum)
       self%freefallRadiusTableRadiusMaximum=2.0d0*self%freefallRadiusTableRadiusMaximum
       self%freefallTimeMaximum=self%freefallTimeScaleFree(self%freefallRadiusTableRadiusMaximum,alphaRequired)
       retabulate=.true.
    end do
    ! Check for alpha out of range.
    if (alphaRequired < self%freefallRadiusTableAlphaMinimum .or. alphaRequired > self%freefallRadiusTableAlphaMaximum) then
       retabulate=.true.
       ! Compute the range of tabulation.
       self%freefallRadiusTableAlphaMinimum=min(self%freefallRadiusTableAlphaMinimum,0.9d0*alphaRequired)
       self%freefallRadiusTableAlphaMaximum=max(self%freefallRadiusTableAlphaMaximum,1.1d0*alphaRequired)
    end if

    if (retabulate) then
       ! Display a message.
       call Galacticus_Display_Indent('Constructing Einasto profile freefall radius lookup table...',verbosityWorking)
       ! Decide how many points to tabulate and allocate table arrays.
       self%freefallRadiusTableRadiusCount=int(log10(self%freefallRadiusTableRadiusMaximum/self%freefallRadiusTableRadiusMinimum)*dble(einastoFreefallRadiusTableRadiusPointsPerDecade))+1
       self%freefallRadiusTableAlphaCount =int(      (self%freefallRadiusTableAlphaMaximum -self%freefallRadiusTableAlphaMinimum )*dble(einastoFreefallRadiusTableAlphaPointsPerUnit   ))+1
       if (allocated(self%freefallRadiusTableRadius)) then
          call Dealloc_Array(self%freefallRadiusTableAlpha )
          call Dealloc_Array(self%freefallRadiusTableRadius)
          call Dealloc_Array(self%freefallRadiusTable      )
       end if
       call Alloc_Array(self%freefallRadiusTableAlpha ,[                               self%freefallRadiusTableAlphaCount])
       call Alloc_Array(self%freefallRadiusTableRadius,[self%freefallRadiusTableRadiusCount                              ])
       call Alloc_Array(self%freefallRadiusTable      ,[self%freefallRadiusTableRadiusCount,self%freefallRadiusTableAlphaCount])
       ! Create a range of radii and alpha.
       self%freefallRadiusTableAlpha =Make_Range(self%freefallRadiusTableAlphaMinimum ,self%freefallRadiusTableAlphaMaximum ,self%freefallRadiusTableAlphaCount ,rangeType=rangeTypeLinear     )
       self%freefallRadiusTableRadius=Make_Range(self%freefallRadiusTableRadiusMinimum,self%freefallRadiusTableRadiusMaximum,self%freefallRadiusTableRadiusCount,rangeType=rangeTypeLogarithmic)
       ! Loop over radii and alpha and populate tables.
       do iAlpha=1,self%freefallRadiusTableAlphaCount
          alpha=self%freefallRadiusTableAlpha(iAlpha)
          do iRadius=1,self%freefallRadiusTableRadiusCount
             ! Show progress.
             percentage=int(100.0d0*dble((iAlpha-1)*self%freefallRadiusTableRadiusCount+iRadius-1)&
                  &/dble(self%freefallRadiusTableAlphaCount*self%freefallRadiusTableRadiusCount))
             call Galacticus_Display_Counter(percentage,iAlpha == 1 .and. iRadius == 1,verbosityWorking)
             ! Compute the freefall radius.
             self%freefallRadiusTable(iRadius,iAlpha)=self%freefallTimeScaleFree(self%freefallRadiusTableRadius(iRadius),alpha)
          end do
       end do
       ! Ensure interpolations get reset.
       call Interpolate_Done(                                                                            &
            &                interpolationAccelerator=self%freefallRadiusTableAlphaInterpolationAccelerator,  &
            &                reset                   =self%freefallRadiusTableAlphaInterpolationReset         &
            &               )
       call Interpolate_Done(                                                                            &
            &                interpolationObject     =self%freefallRadiusTableRadiusInterpolationObject,      &
            &                interpolationAccelerator=self%freefallRadiusTableRadiusInterpolationAccelerator, &
            &                reset                   =self%freefallRadiusTableRadiusInterpolationReset        &
            &               )
       self%freefallRadiusTableAlphaInterpolationReset =.true.
       self%freefallRadiusTableRadiusInterpolationReset=.true.
       ! Store the minimum and maximum tabulated freefall times across all alpha values.
       self%freefallTimeMinimum=maxval(self%freefallRadiusTable(                             1,:))
       self%freefallTimeMaximum=minval(self%freefallRadiusTable(self%freefallRadiusTableRadiusCount,:))
       ! Display a message.
       call Galacticus_Display_Unindent('...done',verbosityWorking)
       ! Specify that tabulation has been made.
       self%freefallRadiusTableInitialized=.true.
    end if
    return
  end subroutine einastoFreefallTabulate

  double precision function einastoFreefallTimeScaleFree(self,radius,alpha)
    !% Compute the freefall time in a scale-free Einasto halo.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    implicit none
    class           (darkMatterProfileEinasto  ), intent(inout) :: self
    double precision                            , intent(in   ) :: alpha               , radius
    type            (c_ptr                     )                :: parameterPointer
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    double precision                                            :: radiusStart         , radiusEnd, &
         &                                                         alphaParameter

    radiusStart   =radius
    radiusEnd     =0.0d0
    alphaParameter=alpha
    einastoFreefallTimeScaleFree=Integrate(radiusEnd,radiusStart,einastoFreefallTimeScaleFreeIntegrand,parameterPointer&
         &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    return

  contains
    
    function einastoFreefallTimeScaleFreeIntegrand(radius,parameterPointer) bind(c)
      !% Integrand function used for finding the free-fall time in Einasto halos.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real(kind=c_double)        :: einastoFreefallTimeScaleFreeIntegrand
      real(kind=c_double), value :: radius
      type(c_ptr        ), value :: parameterPointer
      
      einastoFreefallTimeScaleFreeIntegrand= 1.0d0                                                                   &
           &                                     /sqrt(                                                                  &
           &                                             2.0d0                                                            &
           &                                            *(                                                                &
           &                                               self%potentialScaleFree(radiusStart,1.0d0,alphaParameter) &
           &                                              -self%potentialScaleFree(radius     ,1.0d0,alphaParameter) &
           &                                             )                                                                &
           &                                           )
      return
    end function einastoFreefallTimeScaleFreeIntegrand

  end function einastoFreefallTimeScaleFree

  subroutine einastoStateStore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    implicit none
    class  (darkMatterProfileEinasto), intent(inout) :: self
    integer                          , intent(in   ) :: stateFile
    type   (fgsl_file               ), intent(in   ) :: fgslStateFile

    write (stateFile) self%angularMomentumTableRadiusMinimum,self%angularMomentumTableRadiusMaximum,self%angularMomentumTableAlphaMinimum &
         &,self%angularMomentumTableAlphaMaximum,self%energyTableConcentrationMinimum,self%energyTableConcentrationMaximum &
         &,self%energyTableAlphaMinimum,self%energyTableAlphaMaximum,self%fourierProfileTableWavenumberMinimum &
         &,self%fourierProfileTableWavenumberMaximum,self%fourierProfileTableAlphaMinimum,self%fourierProfileTableAlphaMaximum &
         &,self%fourierProfileTableConcentrationMinimum,self%fourierProfileTableConcentrationMaximum,self%freefallRadiusTableRadiusMinimum&
         &,self%freefallRadiusTableRadiusMaximum,self%freefallRadiusTableAlphaMinimum ,self%freefallRadiusTableAlphaMaximum

    return
  end subroutine einastoStateStore

  subroutine einastoStateRestore(self,stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    implicit none
    class  (darkMatterProfileEinasto), intent(inout) :: self
    integer                          , intent(in   ) :: stateFile
    type   (fgsl_file               ), intent(in   ) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) self%angularMomentumTableRadiusMinimum,self%angularMomentumTableRadiusMaximum,self%angularMomentumTableAlphaMinimum &
         &,self%angularMomentumTableAlphaMaximum,self%energyTableConcentrationMinimum,self%energyTableConcentrationMaximum &
         &,self%energyTableAlphaMinimum,self%energyTableAlphaMaximum,self%fourierProfileTableWavenumberMinimum &
         &,self%fourierProfileTableWavenumberMaximum,self%fourierProfileTableAlphaMinimum,self%fourierProfileTableAlphaMaximum &
         &,self%fourierProfileTableConcentrationMinimum,self%fourierProfileTableConcentrationMaximum,self%freefallRadiusTableRadiusMinimum &
         &,self%freefallRadiusTableRadiusMaximum,self%freefallRadiusTableAlphaMinimum,self%freefallRadiusTableAlphaMaximum
    ! Retabulate.
    self%angularMomentumTableInitialized=.false.
    self%energyTableInitialized         =.false.
    self%fourierProfileTableInitialized =.false.
    self%freefallRadiusTableInitialized =.false.
    return
  end subroutine einastoStateRestore
