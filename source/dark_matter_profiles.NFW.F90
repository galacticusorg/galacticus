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

  !% An implementation of \cite{navarro_universal_1997} dark matter halo profiles.

  !# <darkMatterProfile name="darkMatterProfileNFW">
  !#  <description>\cite{navarro_universal_1997} dark matter halo profiles</description>
  !# </darkMatterProfile>

  use Dark_Matter_Halo_Scales
  use Tables
  use Kind_Numbers

  type, extends(darkMatterProfileClass) :: darkMatterProfileNFW
     !% A dark matter halo profile class implementing \cite{navarro_universal_1997} dark matter halos.
     private
     ! Minimum and maximum concentrations to tabulate.
     double precision                                        :: concentrationMinimum                   , concentrationMaximum
     ! Minimum and maximum radii to tabulate.
     double precision                                        :: freefallRadiusMinimum                  , radiusMinimum
     double precision                                        :: freefallRadiusMaximum                  , radiusMaximum
     double precision                                        :: freefallTimeMinimum                    , specificAngularMomentumMinimum
     double precision                                        :: freefallTimeMaximum                    , specificAngularMomentumMaximum
     ! Tables of NFW properties.
     logical                                                 :: nfwFreefallTableInitialized            , nfwInverseTableInitialized       , &
          &                                                     nfwTableInitialized                    
     integer                                                 :: nfwFreefallTableNumberPoints           , nfwInverseTableNumberPoints      , &
          &                                                     nfwTableNumberPoints
     type            (table1DLogarithmicLinear)              :: nfwConcentrationTable
     ! Tables.
     type            (table1DLogarithmicLinear)              :: nfwFreeFall                            , nfwSpecificAngularMomentum
     class           (table1D                 ), allocatable :: nfwFreefallInverse                     , nfwSpecificAngularMomentumInverse
     ! Module variables used in integrations.
     double precision                                        :: concentrationParameter                 , radiusStart
     ! Record of unique ID of node which we last computed results for.
     integer         (kind=kind_int8          )              :: lastUniqueID
     ! Record of whether or not quantities have been computed.
     logical                                                 :: specificAngularMomentumScalingsComputed
     ! Stored values of computed quantities.
     double precision                                        :: specificAngularMomentumLengthScale     , specificAngularMomentumScale     , &
          &                                                     concentrationPrevious                  , nfwNormalizationFactorPrevious
     ! Pointer to object setting halo scales.
     class(darkMatterHaloScaleClass           ), pointer     :: scale
   contains
     !@ <objectMethods>
     !@   <object>darkMatterProfileNFW</object>
     !@   <objectMethod>
     !@     <method>densityScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin, \doublezero\ concentration\argin, \doublezero\ alpha\argin</arguments>
     !@     <description>Returns the density (in units such that the virial mass and scale length are unity) in an NFW dark matter profile with given {\tt concentration} and {\tt alpha} at the given {\tt radius} (given in units of the scale radius).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>enclosedMassScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin, \doublezero\ concentration\argin, \doublezero\ alpha\argin</arguments>
     !@     <description>Returns the enclosed mass (in units of the virial mass) in an NFW dark matter profile with given {\tt concentration} at the given {\tt radius} (given in units of the scale radius).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>freefallTabulate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ freefallTimeScaleFree\argin, \doublezero\ alphaRequired\argin</arguments>
     !@     <description>Tabulates the freefall time vs. freefall radius for NFW halos.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>freefallTimeScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin, \doublezero\ alpha\argin</arguments>
     !@     <description>Compute the freefall time in a scale-free NFW halo.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>angularMomentumScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ concentration\argin</arguments>
     !@     <description>Returns the total angular momentum in an NFW dark matter profile with given {\tt concentration}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>inverseAngularMomentum</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ specificAngularMomentum\argin</arguments>
     !@     <description>Tabulates the specific angular momentum vs. radius in an NFW profile for rapid inversion.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>profileEnergy</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ concentration\argin</arguments>
     !@     <description>Computes the total energy of an NFW profile halo of given {\tt concentration}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>specificAngularMomentumScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin</arguments>
     !@     <description>Returns the specific angular momentum, normalized to unit scale length and unit velocity at the scale radius, at position {\tt radius} (in units of the scale radius) in an NFW profile.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>tabulate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ concentration\argin</arguments>
     !@     <description>Tabulate properties of the NFW halo profile which must be computed numerically.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final                                             nfwDestructor
     procedure :: calculationReset                  => nfwCalculationReset
     procedure :: stateStore                        => nfwStateStore
     procedure :: stateRestore                      => nfwStateRestore
     procedure :: density                           => nfwDensity
     procedure :: enclosedMass                      => nfwEnclosedMass
     procedure :: potential                         => nfwPotential
     procedure :: circularVelocity                  => nfwCircularVelocity
     procedure :: circularVelocityMaximum           => nfwCircularVelocityMaximum
     procedure :: radiusFromSpecificAngularMomentum => nfwRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => nfwRotationNormalization
     procedure :: energy                            => nfwEnergy
     procedure :: energyGrowthRate                  => nfwEnergyGrowthRate
     procedure :: kSpace                            => nfwKSpace
     procedure :: freefallRadius                    => nfwFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => nfwFreefallRadiusIncreaseRate
     procedure :: profileEnergy                     => nfwProfileEnergy
     procedure :: specificAngularMomentumScaleFree  => nfwSpecificAngularMomentumScaleFree
     procedure :: angularMomentumScaleFree          => nfwAngularMomentumScaleFree
     procedure :: enclosedMassScaleFree             => nfwEnclosedMassScaleFree
     procedure :: densityScaleFree                  => nfwDensityScaleFree
     procedure :: tabulate                          => nfwTabulate
     procedure :: inverseAngularMomentum            => nfwInverseAngularMomentum
     procedure :: freefallTabulate                  => nfwFreefallTabulate
     procedure :: freefallTimeScaleFree             => nfwFreefallTimeScaleFree
  end type darkMatterProfileNFW

  interface darkMatterProfileNFW
     !% Constructors for the {\tt nfw} dark matter halo profile class.
     module procedure nfwDefaultConstructor
     module procedure nfwConstructor
  end interface darkMatterProfileNFW

  ! Number of points per decade of concentration in NFW tabulations.
  integer, parameter   :: nfwTablePointsPerDecade        =100
  integer, parameter   :: nfwInverseTablePointsPerDecade =100
  integer, parameter   :: nfwFreefallTablePointsPerDecade=100
  ! Indices for tabulated quantities.
  integer, parameter   :: nfwConcentrationEnergyIndex    =  1, nfwConcetrationRotationNormalizationIndex=2

contains

  function nfwDefaultConstructor()
    !% Default constructor for the {\tt nfw} dark matter halo profile class.
    use Input_Parameters
    implicit none
    type(darkMatterProfileNFW), target :: nfwDefaultConstructor

    nfwDefaultConstructor=nfwConstructor(darkMatterHaloScale())
    return
  end function nfwDefaultConstructor

  function nfwConstructor(scale)
    !% Generic constructor for the {\tt nfw} dark matter halo profile class.
    use Galacticus_Error
    implicit none
    type (darkMatterProfileNFW    ), target :: nfwConstructor
    class(darkMatterHaloScaleClass), target :: scale
 
    nfwConstructor%concentrationPrevious       =  -1.0d+0
    nfwConstructor%concentrationMinimum        =   1.0d+0
    nfwConstructor%concentrationMaximum        =  20.0d+0
    nfwConstructor%freefallRadiusMinimum       =   1.0d-3 
    nfwConstructor%radiusMinimum               =   1.0d-3
    nfwConstructor%freefallRadiusMaximum       =   1.0d+2  
    nfwConstructor%radiusMaximum               =   1.0d+2
    nfwConstructor%nfwFreefallTableInitialized =  .false.
    nfwConstructor%nfwInverseTableInitialized  =  .false.
    nfwConstructor%nfwTableInitialized         =  .false.
    nfwConstructor%lastUniqueID                =  -1
    nfwConstructor%scale                       => scale
    ! Ensure that the dark matter profile component supports a "scale" property.
    if (.not.defaultDarkMatterProfileComponent%scaleIsGettable())                                                         &
         & call Galacticus_Error_Report                                                                                   &
         &      (                                                                                                         &
         &       'nfwConstructor'                                                                                       , &
         &       'NFW dark matter profile requires a dark matter profile component with a gettable "scale" property.'//   &
         &       Galacticus_Component_List(                                                                               &
         &                                 'darkMatterProfile'                                                          , &
         &                                 defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)  &
         &                                )                                                                               &
         &      )
    ! Initialize the tabulations.
    call nfwConstructor%tabulate              ()
    call nfwConstructor%inverseAngularMomentum()
    return
  end function nfwConstructor
  
  subroutine nfwDestructor(self)
    !% Destructor for the {\tt nfw} dark matter halo profile class.
    implicit none
    type(darkMatterProfileNFW), intent(inout) :: self

    if (self%nfwFreefallTableInitialized) then
       call self%nfwFreeFall                      %destroy()
       call self%nfwFreeFallInverse               %destroy()
       deallocate(self%nfwFreefallInverse)
    end if
    if (self%nfwInverseTableInitialized ) then
       call self%nfwSpecificAngularMomentum       %destroy()
       call self%nfwSpecificAngularMomentumInverse%destroy()
       deallocate(self%nfwSpecificAngularMomentumInverse)
    end if
    if (self%nfwTableInitialized        ) then
       call self%nfwConcentrationTable            %destroy()
    end if
    if (self%scale%isFinalizable()) deallocate(self%scale)
    return
  end subroutine nfwDestructor
  
  subroutine nfwCalculationReset(self,thisNode)
    !% Reset the dark matter profile calculation.
    implicit none
    class(darkMatterProfileNFW), intent(inout)          :: self
    type (treeNode            ), intent(inout), pointer :: thisNode

    self%specificAngularMomentumScalingsComputed=.false.
    self%lastUniqueID                           =thisNode%uniqueID()
    call self%scale%calculationReset(thisNode)
    return
  end subroutine nfwCalculationReset

  subroutine nfwTabulate(self,concentration)
    !% Tabulate properties of the NFW halo profile which must be computed numerically.
    implicit none
    class           (darkMatterProfileNFW), intent(inout)           :: self
    double precision                      , intent(in   ), optional :: concentration
    integer                                                         :: iConcentration
    logical                                                         :: retabulate
    double precision                                                :: tableConcentration

    retabulate=.not.self%nfwTableInitialized
    if (present(concentration)) then
       if (concentration < self%concentrationMinimum) then
          self%concentrationMinimum=0.5d0*concentration
          retabulate=.true.
       end if
       if (concentration > self%concentrationMaximum) then
          self%concentrationMaximum=2.0d0*concentration
          retabulate=.true.
       end if
    end if
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%nfwTableNumberPoints=int(log10(self%concentrationMaximum/self%concentrationMinimum)*dble(nfwTablePointsPerDecade))+1
       call self%nfwConcentrationTable%destroy()
       call self%nfwConcentrationTable%create(self%concentrationMinimum,self%concentrationMaximum,self%nfwTableNumberPoints,2)
       ! Loop over concentrations and populate tables.
       do iConcentration=1,self%nfwTableNumberPoints
          tableConcentration=self%nfwConcentrationTable%x(iConcentration)
          call self%nfwConcentrationTable%populate(                   self%profileEnergy           (tableConcentration),iConcentration,table=nfwConcentrationEnergyIndex              )
          call self%nfwConcentrationTable%populate(tableConcentration/self%angularMomentumScaleFree(tableConcentration),iConcentration,table=nfwConcetrationRotationNormalizationIndex)
       end do
       ! Specify that tabulation has been made.
       self%nfwTableInitialized=.true.
    end if
    return
  end subroutine nfwTabulate

  subroutine nfwInverseAngularMomentum(self,specificAngularMomentum)
    !% Tabulates the specific angular momentum vs. radius in an NFW profile for rapid inversion.
    implicit none
    class           (darkMatterProfileNFW), intent(inout)           :: self
    double precision                      , intent(in   ), optional :: specificAngularMomentum
    integer                                                         :: iRadius
    logical                                                         :: retabulate

    retabulate=.not.self%nfwInverseTableInitialized
    ! If the table has not yet been made, compute and store the specific angular momenta corresponding to the minimum and maximum
    ! radii that will be tabulated by default.
    if (retabulate) then
       self%specificAngularMomentumMinimum=self%specificAngularMomentumScaleFree(self%radiusMinimum)
       self%specificAngularMomentumMaximum=self%specificAngularMomentumScaleFree(self%radiusMaximum)
    end if
    if (present(specificAngularMomentum)) then
       do while (specificAngularMomentum < self%specificAngularMomentumMinimum)
          self%radiusMinimum                 =0.5d0*self%radiusMinimum
          self%specificAngularMomentumMinimum=self%specificAngularMomentumScaleFree(self%radiusMinimum)
          retabulate=.true.
       end do
       do while (specificAngularMomentum > self%specificAngularMomentumMaximum)
          self%radiusMaximum                 =2.0d0*self%radiusMaximum
          self%specificAngularMomentumMaximum=self%specificAngularMomentumScaleFree(self%radiusMaximum)
          retabulate=.true.
       end do
    end if
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%nfwInverseTableNumberPoints=int(log10(self%radiusMaximum/self%radiusMinimum)*dble(nfwInverseTablePointsPerDecade))+1
       ! Create a range of radii.
       call self%nfwSpecificAngularMomentum%destroy(                                                       )
       call self%nfwSpecificAngularMomentum%create (self%radiusMinimum,self%radiusMaximum,self%nfwInverseTableNumberPoints)
       ! Loop over radii and populate tables.
       do iRadius=1,self%nfwInverseTableNumberPoints
          call self%nfwSpecificAngularMomentum%populate(                                                                                   &
               &                                        self%specificAngularMomentumScaleFree(self%nfwSpecificAngularMomentum%x(iRadius)), &
               &                                        iRadius                                                                            &
               &                                       )
       end do
       call self%nfwSpecificAngularMomentum%reverse(self%nfwSpecificAngularMomentumInverse)
       ! Specify that tabulation has been made.
       self%nfwInverseTableInitialized=.true.
    end if
    return
  end subroutine nfwInverseAngularMomentum

  double precision function nfwDensity(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\tt node} at the given {\tt radius} (given
    !% in units of Mpc).
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileNFW          ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: radius
    class           (nodeComponentBasic            )               , pointer :: thisBasicComponent
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: radiusOverScaleRadius         , scaleRadius, &
         &                                                                      virialRadiusOverScaleRadius

    thisBasicComponent             => node%basic            (                 )
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                    =thisDarkMatterProfileComponent%scale()
    radiusOverScaleRadius          =radius                           /scaleRadius
    virialRadiusOverScaleRadius    =self%scale%virialRadius(node)/scaleRadius
    nfwDensity=self%densityScaleFree(radiusOverScaleRadius,virialRadiusOverScaleRadius)&
         &*thisBasicComponent%mass()/scaleRadius**3
    return
  end function nfwDensity

  double precision function nfwEnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\tt node} at the given {\tt radius} (given in
    !% units of Mpc).
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileNFW          ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: radius
    class           (nodeComponentBasic            )               , pointer :: thisBasicComponent
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: radiusOverScaleRadius         , scaleRadius, &
         &                                                                      virialRadiusOverScaleRadius

    thisBasicComponent             => node%basic            (                 )
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)
    scaleRadius                    =thisDarkMatterProfileComponent%scale()
    radiusOverScaleRadius          =radius                       /scaleRadius
    virialRadiusOverScaleRadius    =self%scale%virialRadius(node)/scaleRadius
    nfwEnclosedMass                =self%enclosedMassScaleFree(radiusOverScaleRadius,virialRadiusOverScaleRadius) &
         &*thisBasicComponent%mass()
    return
  end function nfwEnclosedMass

  double precision function nfwPotential(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\tt node} at the given {\tt radius} (given in
    !% units of Mpc).
    use Dark_Matter_Halo_Scales
    use Dark_Matter_Profiles_Error_Codes
    implicit none
    class           (darkMatterProfileNFW          ), intent(inout)           :: self
    type            (treeNode                      ), intent(inout), pointer  :: node
    double precision                                , intent(in   )           :: radius
    integer                                         , intent(  out), optional :: status
    class           (nodeComponentDarkMatterProfile)               , pointer  :: thisDarkMatterProfileComponent
    double precision                                , parameter               :: radiusSmall                   =1.0d-10
    double precision                                                          :: radiusOverScaleRadius                 , radiusTerm, &
         &                                                                       virialRadiusOverScaleRadius

    if (present(status)) status=darkMatterProfileSuccess
    thisDarkMatterProfileComponent   => node%darkMatterProfile(autoCreate=.true.)
    radiusOverScaleRadius            =radius                           /thisDarkMatterProfileComponent%scale()
    virialRadiusOverScaleRadius      =self%scale%virialRadius(node)/thisDarkMatterProfileComponent%scale()
    if (radiusOverScaleRadius < radiusSmall) then
       ! Use a series solution for very small radii.
       radiusTerm=1.0d0-0.5d0*radiusOverScaleRadius
    else
       ! Use the full expression for larger radii.
       radiusTerm=log(1.0d0+radiusOverScaleRadius)/radiusOverScaleRadius
    end if
    nfwPotential=(-1.0d0-virialRadiusOverScaleRadius*(radiusTerm-log(1.0d0+virialRadiusOverScaleRadius)&
         &/virialRadiusOverScaleRadius)/(log(1.0d0 +virialRadiusOverScaleRadius)-virialRadiusOverScaleRadius/(1.0d0&
         &+virialRadiusOverScaleRadius)))*self%scale%virialVelocity(node)**2
    return
  end function nfwPotential

  double precision function nfwCircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\tt node} at the given {\tt radius} (given in
    !% units of Mpc).
    use Numerical_Constants_Physical
    implicit none
    class           (darkMatterProfileNFW), intent(inout)          :: self
    type            (treeNode            ), intent(inout), pointer :: node
    double precision                      , intent(in   )          :: radius

    if (radius > 0.0d0) then
       nfwCircularVelocity=sqrt(gravitationalConstantGalacticus&
            &*self%enclosedMass(node,radius)/radius)
    else
       nfwCircularVelocity=0.0d0
    end if
    return
  end function nfwCircularVelocity

  double precision function nfwCircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\tt node}.
    use Numerical_Constants_Physical
    implicit none
    class           (darkMatterProfileNFW          ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    ! The radius (in units of the scale radius) at which the rotation speed peaks in an NFW halo.
    double precision                                , parameter              :: radiusMaximum=2.1625815870646097d0
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfile
    double precision                                                         :: scaleRadius

    thisDarkMatterProfile      => node                 %darkMatterProfile(autoCreate=.true.                          )
    scaleRadius                =  thisDarkMatterProfile%scale            (                                           )
    nfwCircularVelocityMaximum =  self                 %circularVelocity (node             ,radiusMaximum*scaleRadius)
    return
  end function nfwCircularVelocityMaximum

  double precision function nfwRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\tt node} at which a circular orbit has the given {\tt specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc). For an NFW halo, the circular velocity is constant (and therefore equal to the virial
    !% velocity). Therefore, $r = j/V_{\rm virial}$ where $j$(={\tt specificAngularMomentum}) is the specific angular momentum and
    !% $r$ the required radius.
    implicit none
    class           (darkMatterProfileNFW          ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: specificAngularMomentum
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: specificAngularMomentumScaleFree

    ! Return immediately with zero radius for non-positive specific angular momenta.
    if (specificAngularMomentum <= 0.0d0) then
       nfwRadiusFromSpecificAngularMomentum=0.0d0
       return
    end if
    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Check if scalings are already computed. Compute and store if not.
    if (.not.self%specificAngularMomentumScalingsComputed) then
       ! Flag that scale quantities are now computed.
       self%specificAngularMomentumScalingsComputed=.true.

       ! Get the dark matter profile.
       thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)

       ! Get the scale radius.
       self%specificAngularMomentumLengthScale=thisDarkMatterProfileComponent%scale()

       ! Get the specific angular momentum scale.
       self%specificAngularMomentumScale=self%specificAngularMomentumLengthScale*self%circularVelocity(node&
            &,self%specificAngularMomentumLengthScale)
    end if

    ! Compute the specific angular momentum in scale free units (using the scale length for distances the sqrt(G M(r_scale) /
    ! r_scale) for velocities).
    specificAngularMomentumScaleFree=specificAngularMomentum/self%specificAngularMomentumScale

    ! Ensure that the interpolations exist and extend sufficiently far.
    call self%inverseAngularMomentum(specificAngularMomentumScaleFree)

    ! Interpolate to get the dimensionless radius at which this specific angular momentum is found.
    nfwRadiusFromSpecificAngularMomentum=self%nfwSpecificAngularMomentumInverse%interpolate(specificAngularMomentumScaleFree)

    ! Convert to a physical radius.
    nfwRadiusFromSpecificAngularMomentum=nfwRadiusFromSpecificAngularMomentum*self%specificAngularMomentumLengthScale
    return
  end function nfwRadiusFromSpecificAngularMomentum

  double precision function nfwRotationNormalization(self,node)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileNFW          ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: concentration

    ! Get components.
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)

    ! Find the concentration parameter of this halo.
    concentration=self%scale%virialRadius(node)/thisDarkMatterProfileComponent%scale()

    ! Ensure that the interpolations exist and extend sufficiently far.
    call self%tabulate(concentration)

    ! Find the energy by interpolation.
    nfwRotationNormalization=self%nfwConcentrationTable%interpolate(concentration,table&
         &=nfwConcentrationEnergyIndex)/self%scale%virialRadius(node)
    return
  end function nfwRotationNormalization

  double precision function nfwEnergy(self,node)
    !% Return the energy of an NFW halo density profile.
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileNFW          ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    class           (nodeComponentBasic            )               , pointer :: thisBasicComponent
    double precision                                                         :: concentration

    ! Get components.
    thisBasicComponent             => node%basic            (                 )
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)

    ! Find the concentration parameter of this halo.
    concentration=self%scale%virialRadius(node)/thisDarkMatterProfileComponent%scale()

    ! Ensure that the interpolations exist and extend sufficiently far.
    call self%tabulate(concentration)

    ! Find the energy by interpolation.
    nfwEnergy=self%nfwConcentrationTable%interpolate(concentration,table=nfwConcentrationEnergyIndex)&
         &*thisBasicComponent%mass()*self%scale%virialVelocity(node)**2
    return
  end function nfwEnergy

  double precision function nfwEnergyGrowthRate(self,node)
    !% Return the rate of change of the energy of an NFW halo density profile.
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileNFW          ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    class           (nodeComponentBasic            )               , pointer :: thisBasicComponent
    double precision                                                         :: concentration                 , energy, &
         &                                                                      energyGradient

    ! Get components.
    thisBasicComponent             => node%basic            (                 )
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)

    ! Find the concentration parameter of this halo.
    concentration=self%scale%virialRadius(node)/thisDarkMatterProfileComponent%scale()

    ! Ensure that the interpolations exist and extend sufficiently far.
    call self%tabulate(concentration)

    ! Find the energy gradient by interpolation.
    energy        =self%nfwConcentrationTable%interpolate        (concentration,table=nfwConcentrationEnergyIndex)
    energyGradient=self%nfwConcentrationTable%interpolateGradient(concentration,table=nfwConcentrationEnergyIndex)

    nfwEnergyGrowthRate=self%energy(node)&
         &*(thisBasicComponent%accretionRate()/thisBasicComponent%mass()+2.0d0 &
         &*self%scale%virialVelocityGrowthRate(node)/self%scale%virialVelocity(node)+(energyGradient&
         &*concentration/energy)*(self%scale%virialRadiusGrowthRate(node)/self%scale%virialRadius(node)&
         &-thisDarkMatterProfileComponent%scaleGrowthRate()/thisDarkMatterProfileComponent%scale()))

    return
  end function nfwEnergyGrowthRate

  double precision function nfwAngularMomentumScaleFree(self,concentration)
    !% Returns the total angular momentum (in units of the virial mass times scale radius times [assumed constant] rotation speed)
    !% in an NFW dark matter profile with given {\tt concentration}. This is given by:
    !% \begin{equation}
    !% J = \left. \int_0^c 4 \pi x^3 \rho(x) \d x \right/ \int_0^c 4 \pi x^2 \rho(x) \d x,
    !% \end{equation}
    !% where $x$ is radius in units of the scale radius and $c$ is concentration. This can be evaluated to give
    !% \begin{equation}
    !% J = \left. \left[ 1 + c - 2 \ln (1+c) - {1 \over 1+c} \right] \right/ \left[ \ln(1+c)-{c\over 1+c} \right].
    !% \end{equation}
    implicit none
    class           (darkMatterProfileNFW), intent(inout) :: self
    double precision                      , intent(in   ) :: concentration

    nfwAngularMomentumScaleFree=(1.0d0+concentration-2.0d0*log(1.0d0+concentration)-1.0d0/(1.0d0+concentration)) &
         &/(log(1.0d0+concentration)-concentration/(1.0d0+concentration))
    return
  end function nfwAngularMomentumScaleFree

  double precision function nfwSpecificAngularMomentumScaleFree(self,radius)
    !% Returns the specific angular momentum, normalized to unit scale length and unit velocity at the scale radius, at position
    !% {\tt radius} (in units of the scale radius) in an NFW profile.
    implicit none
    class           (darkMatterProfileNFW), intent(inout) :: self
    double precision                      , intent(in   ) :: radius

    nfwSpecificAngularMomentumScaleFree=sqrt(radius*self%enclosedMassScaleFree(radius,1.0d0))
    return
  end function nfwSpecificAngularMomentumScaleFree

  double precision function nfwEnclosedMassScaleFree(self,radius,concentration)
    !% Returns the enclosed mass (in units of the virial mass) in an NFW dark matter profile with given {\tt concentration} at the
    !% given {\tt radius} (given in units of the scale radius).
    implicit none
    class           (darkMatterProfileNFW), intent(inout) :: self
    double precision                      , intent(in   ) :: concentration                                                   , radius
    double precision                      , parameter     :: minimumRadiusForExactSolution          =1.0d-7
    ! Precomputed NFW normalization factor for unit concentration.
    double precision                      , parameter     :: nfwNormalizationFactorUnitConcentration=1.0d0/(log(2.0d0)-0.5d0)
    ! Precomputed NFW normalization factor for unit radius.
    double precision                      , parameter     :: nfwNormalizationFactorUnitRadius       =log(2.0d0)-0.5d0

    if (radius == 1.0d0) then
       nfwEnclosedMassScaleFree=nfwNormalizationFactorUnitRadius
    else if (radius >= minimumRadiusForExactSolution) then
       nfwEnclosedMassScaleFree=(log(1.0d0+radius)-radius/(1.0d0+radius))
    else
       nfwEnclosedMassScaleFree=(radius**2)*(0.5d0+radius*(-2.0d0/3.0d0+radius*(0.75d0+radius*(-0.8d0))))
    end if
    ! Check if we were called with a different concentration compared to the previous call.
    if (concentration /= self%concentrationPrevious) then
       ! We were, so recompute the normalization factor.
       if (concentration == 1.0d0) then
          self%nfwNormalizationFactorPrevious=nfwNormalizationFactorUnitConcentration
       else
          self%nfwNormalizationFactorPrevious=1.0d0/(log(1.0d0+concentration)-concentration/(1.0d0 &
               &+concentration))
       end if
       self%concentrationPrevious=concentration
    end if
    nfwEnclosedMassScaleFree=nfwEnclosedMassScaleFree*self%nfwNormalizationFactorPrevious
    return
  end function nfwEnclosedMassScaleFree

  double precision function nfwDensityScaleFree(self,radius,concentration)
    !% Returns the density (in units such that the virial mass and scale length are unity) in an NFW dark matter profile with
    !% given {\tt concentration} at the given {\tt radius} (given in units of the scale radius).
    use Numerical_Constants_Math
    implicit none
    class           (darkMatterProfileNFW), intent(inout) :: self
    double precision                      , intent(in   ) :: concentration, radius

    nfwDensityScaleFree=1.0d0/(log(1.0d0+concentration)-concentration/(1.0d0+concentration))/radius/(1.0d0+radius)**2/4.0d0/Pi
    return
  end function nfwDensityScaleFree

  double precision function nfwProfileEnergy(self,concentration)
    !% Computes the total energy of an NFW profile halo of given {\tt concentration} using the methods of
    !% \citeauthor{cole_hierarchical_2000}~(\citeyear{cole_hierarchical_2000}; their Appendix~A).
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Math
    use Numerical_Integration
    implicit none
    class           (darkMatterProfileNFW      ), intent(inout) :: self
    double precision                            , intent(in   ) :: concentration
    type            (c_ptr                     )                :: parameterPointer
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    double precision                                            :: jeansEquationIntegral  , kineticEnergy         , &
         &                                                         kineticEnergyIntegral  , potentialEnergy       , &
         &                                                         potentialEnergyIntegral, radiusMaximum         , &
         &                                                         radiusMinimum          , concentrationParameter

    ! Compute the potential energy.
    radiusMinimum    =0.0d0
    radiusMaximum    =concentration
    concentrationParameter=concentration
    potentialEnergyIntegral=Integrate(radiusMinimum,radiusMaximum,nfwPotentialEnergyIntegrand&
         &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    potentialEnergy=-0.5d0*(1.0d0/concentration+potentialEnergyIntegral)
    ! Compute the velocity dispersion at the virial radius.
    radiusMinimum=concentration
    radiusMaximum=100.0d0*concentration
    concentrationParameter=concentration
    jeansEquationIntegral=Integrate(radiusMinimum,radiusMaximum,nfwJeansEquationIntegrand&
         &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    ! Compute the kinetic energy.
    radiusMinimum=0.0d0
    radiusMaximum=concentration
    concentrationParameter=concentration
    kineticEnergyIntegral=Integrate(radiusMinimum,radiusMaximum,nfwKineticEnergyIntegrand&
         &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    kineticEnergy=2.0d0*Pi*(jeansEquationIntegral*concentration**3+kineticEnergyIntegral)
    ! Compute the total energy.
    nfwProfileEnergy=(potentialEnergy+kineticEnergy)*concentration
    return

  contains
    
    function nfwPotentialEnergyIntegrand(radius,parameterPointer) bind(c)
      !% Integrand for NFW profile potential energy.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real(kind=c_double)        :: nfwPotentialEnergyIntegrand
      real(kind=c_double), value :: radius
      type(c_ptr        ), value :: parameterPointer
      
      nfwPotentialEnergyIntegrand=(self%enclosedMassScaleFree(radius,concentrationParameter)/radius)**2
      return
    end function nfwPotentialEnergyIntegrand
    
    function nfwKineticEnergyIntegrand(radius,parameterPointer) bind(c)
      !% Integrand for NFW profile kinetic energy.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real(kind=c_double)        :: nfwKineticEnergyIntegrand
      real(kind=c_double), value :: radius
      type(c_ptr        ), value :: parameterPointer
      
      nfwKineticEnergyIntegrand=self%EnclosedMassScaleFree(radius,concentrationParameter)*self%densityScaleFree(radius&
           &,concentrationParameter)*radius
      return
    end function nfwKineticEnergyIntegrand
    
    function nfwJeansEquationIntegrand(radius,parameterPointer) bind(c)
      !% Integrand for NFW profile Jeans equation.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real(kind=c_double)        :: nfwJeansEquationIntegrand
      real(kind=c_double), value :: radius
      type(c_ptr        ), value :: parameterPointer
      
      nfwJeansEquationIntegrand=self%enclosedMassScaleFree(radius,concentrationParameter)*self%densityScaleFree(radius &
           &,concentrationParameter)/radius**2
      return
    end function nfwJeansEquationIntegrand

  end function nfwProfileEnergy
  
  double precision function nfwKSpace(self,node,waveNumber)
    !% Returns the Fourier transform of the NFW density profile at the specified {\tt waveNumber} (given in Mpc$^{-1}$), using the
    !% expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    use Dark_Matter_Halo_Scales
    use Exponential_Integrals
    implicit none
    class           (darkMatterProfileNFW          ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: waveNumber
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: concentration                 , radiusScale, &
         &                                                                      waveNumberScaleFree

    ! Get components.
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)

    ! Get the scale radius.
    radiusScale=thisDarkMatterProfileComponent%scale()

    ! Compute the concentration parameter.
    concentration=self%scale%virialRadius(node)/radiusScale

    ! Get the dimensionless wavenumber.
    waveNumberScaleFree=waveNumber*radiusScale

    ! Compute the Fourier transformed profile.
    nfwKSpace= (                                                                                                                                          &
         &      +sin(              waveNumberScaleFree)*(Sine_Integral  ((1.0d0+concentration)*waveNumberScaleFree)-Sine_Integral  (waveNumberScaleFree)) &
         &      -sin(concentration*waveNumberScaleFree)/(1.0d0+concentration)/waveNumberScaleFree                                                         &
         &      +cos(              waveNumberScaleFree)*(Cosine_Integral((1.0d0+concentration)*waveNumberScaleFree)-Cosine_Integral(waveNumberScaleFree)) &
         &     )                                                                                                                                          &
         &    /(log(1.0d0+concentration)-concentration/(1.0d0+concentration))
    return
  end function nfwKSpace

  double precision function nfwFreefallRadius(self,node,time)
    !% Returns the freefall radius in the NFW density profile at the specified {\tt time} (given in Gyr).
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    implicit none
    class           (darkMatterProfileNFW          ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: time
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: concentration                 , freefallTimeScaleFree, &
         &                                                                      radiusScale                   , timeScale            , &
         &                                                                      velocityScale

    ! For non-positive freefall times, return a zero freefall radius immediately.
    if (time <= 0.0d0) then
       nfwFreefallRadius=0.0d0
       return
    end if

    ! Get components.
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)

    ! Get the scale radius.
    radiusScale=thisDarkMatterProfileComponent%scale()

    ! Get the concentration.
    concentration=self%scale%virialRadius(node)/radiusScale

    ! Get the virial velocity.
    velocityScale=self%scale%virialVelocity(node)

    ! Compute time scale.
    timeScale=Mpc_per_km_per_s_To_Gyr*radiusScale/velocityScale/sqrt(concentration/(log(1.0d0+concentration)-concentration&
         &/(1.0d0+concentration)))

    ! Compute dimensionless time.
    freefallTimeScaleFree=time/timeScale

    ! Ensure table is sufficiently extensive.
    call self%freefallTabulate(freefallTimeScaleFree)

    ! Interpolate to get the freefall radius.
    nfwFreefallRadius=self%nfwFreefallInverse%interpolate(freefallTimeScaleFree)*radiusScale
    return
  end function nfwFreefallRadius

  double precision function nfwFreefallRadiusIncreaseRate(self,node,time)
    !% Returns the rate of increase of the freefall radius in the NFW density profile at the specified {\tt time} (given in
    !% Gyr).
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    implicit none
    class           (darkMatterProfileNFW          ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: time
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: concentration                 , freefallTimeScaleFree, &
         &                                                                      radiusScale                   , timeScale            , &
         &                                                                      velocityScale

    ! For non-positive freefall times, return the limiting value for small radii.
    if (time <= 0.0d0) then
       nfwFreefallRadiusIncreaseRate=0.0d0
       return
    end if

    ! Get components.
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)

    ! Get the scale radius.
    radiusScale=thisDarkMatterProfileComponent%scale()

    ! Get the concentration.
    concentration=self%scale%virialRadius(node)/radiusScale

    ! Get the virial velocity.
    velocityScale=self%scale%virialVelocity(node)

    ! Compute time scale.
    timeScale=Mpc_per_km_per_s_To_Gyr*radiusScale/velocityScale/sqrt(concentration/(log(1.0d0+concentration)-concentration&
         &/(1.0d0+concentration)))

    ! Compute dimensionless time.
    freefallTimeScaleFree=time/timeScale

    ! Ensure table is sufficiently extensive.
    call self%freefallTabulate(freefallTimeScaleFree)

    ! Interpolate to get the freefall radius growth rate.
    nfwFreefallRadiusIncreaseRate=self%nfwFreefallInverse%interpolateGradient(freefallTimeScaleFree)*radiusScale/timeScale
    return
  end function nfwFreefallRadiusIncreaseRate

  subroutine nfwFreefallTabulate(self,freefallTimeScaleFree)
    !% Tabulates the freefall time vs. freefall radius for NFW halos.
    implicit none
    class           (darkMatterProfileNFW), intent(inout) :: self
    double precision                      , intent(in   ) :: freefallTimeScaleFree
    logical                                               :: retabulate
    integer                                               :: iRadius

    retabulate=.not.self%nfwFreefallTableInitialized
    ! If the table has not yet been made, compute and store the freefall corresponding to the minimum and maximum
    ! radii that will be tabulated by default.
    if (retabulate) then
       self%freefallTimeMinimum=self%freefallTimeScaleFree(self%freefallRadiusMinimum)
       self%freefallTimeMaximum=self%freefallTimeScaleFree(self%freefallRadiusMaximum)
    end if
    do while (freefallTimeScaleFree < self%freefallTimeMinimum)
       self%freefallRadiusMinimum=0.5d0*self%freefallRadiusMinimum
       self%freefallTimeMinimum=self%freefallTimeScaleFree(self%freefallRadiusMinimum)
       retabulate=.true.
    end do
    do while (freefallTimeScaleFree > self%freefallTimeMaximum)
       self%freefallRadiusMaximum=2.0d0*self%freefallRadiusMaximum
       self%freefallTimeMaximum=self%freefallTimeScaleFree(self%freefallRadiusMaximum)
       retabulate=.true.
    end do
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%nfwFreefallTableNumberPoints=int(log10(self%freefallRadiusMaximum/self%freefallRadiusMinimum)*dble(nfwFreefallTablePointsPerDecade))+1
       ! Create the table.
       call self%nfwFreefall%destroy(                                                                        )
       call self%nfwFreefall%create (self%freefallRadiusMinimum,self%freefallRadiusMaximum,self%nfwFreefallTableNumberPoints)
       ! Loop over radii and populate tables.
       do iRadius=1,self%nfwFreefallTableNumberPoints
          call self%nfwFreefall%populate(                                                         &
               &                         self%freefallTimeScaleFree(self%nfwFreefall%x(iRadius)), &
               &                                                                       iRadius    &
               &                        )
       end do
       call self%nfwFreefall%reverse(self%nfwFreefallInverse)
       ! Specify that tabulation has been made.
       self%nfwFreefallTableInitialized=.true.
    end if
    return
  end subroutine nfwFreefallTabulate

  double precision function nfwFreefallTimeScaleFree(self,radius)
    !% Compute the freefall time in a scale-free NFW halo.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    implicit none
    class           (darkMatterProfileNFW      ), intent(inout) :: self
    double precision                            , intent(in   ) :: radius
    double precision                            , parameter     :: radiusSmall         =4.0d-6
    type            (c_ptr                     )                :: parameterPointer
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    double precision                                            :: radiusEnd                  , radiusStart

    if (radius > radiusSmall) then
       ! Use the full solution.
       radiusStart=radius
       radiusEnd  =0.0d0
       nfwFreefallTimeScaleFree=Integrate(radiusEnd,radiusStart,nfwFreefallTimeScaleFreeIntegrand,parameterPointer&
            &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
       call Integrate_Done(integrandFunction,integrationWorkspace)
    else
       ! Use an approximation here, found by taking series expansions of the logarithms in the integrand and keeping only the
       ! first order terms.
       nfwFreefallTimeScaleFree=2.0d0*sqrt(radius)
    end if
    return

  contains
    
    function nfwFreefallTimeScaleFreeIntegrand(radius,parameterPointer) bind(c)
      !% Integrand function used for finding the free-fall time in NFW halos.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real(kind=c_double)                   :: nfwFreefallTimeScaleFreeIntegrand
      real(kind=c_double)           , value :: radius
      type(c_ptr        )           , value :: parameterPointer
      real(kind=c_double), parameter        :: radiusSmall                       =1.0d-6
      real(kind=c_double), parameter        :: radiusSmallFraction               =1.0d-3
      real(kind=c_double)                   :: x
      
      if (radius < radiusSmall) then
         ! Use a series approximation for small radii.
         nfwFreefallTimeScaleFreeIntegrand=log(1.0d0+radiusStart)/radiusStart-1.0d0+radius*(0.5d0-radius/3.0d0)
      else if (radius > radiusStart*(1.0d0-radiusSmallFraction)) then
         ! Use a series approximation for radii close to the initial radius.
         x=1.0d0-radius/radiusStart
         nfwFreefallTimeScaleFreeIntegrand=(1.0d0/(1.0d0+radiusStart)-log(1.0d0+radiusStart)/radiusStart)*x+(0.5d0*radiusStart&
              &/(1.0d0+radiusStart)**2+(radiusStart-(1.0d0+radiusStart)*log(1.0d0+radiusStart))/radiusStart/(1.0d0+radiusStart))*x&
              &**2
      else
         ! Use full expression for larger radii.
         nfwFreefallTimeScaleFreeIntegrand=log(1.0d0+radiusStart)/radiusStart-log(1.0d0+radius)/radius
      end if
      nfwFreefallTimeScaleFreeIntegrand=1.0d0/sqrt(-2.0d0*nfwFreefallTimeScaleFreeIntegrand)
      return
    end function nfwFreefallTimeScaleFreeIntegrand

  end function nfwFreefallTimeScaleFree

  subroutine nfwStateStore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    implicit none
    class  (darkMatterProfileNFW), intent(inout) :: self
    integer                      , intent(in   ) :: stateFile
    type   (fgsl_file           ), intent(in   ) :: fgslStateFile

    write (stateFile) self%concentrationMinimum,self%concentrationMaximum,self%radiusMinimum,self%radiusMaximum,self%freefallRadiusMinimum,self%freefallRadiusMaximum
    return
  end subroutine nfwStateStore

  subroutine nfwStateRestore(self,stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    implicit none
    class  (darkMatterProfileNFW), intent(inout) :: self
    integer                      , intent(in   ) :: stateFile
    type   (fgsl_file           ), intent(in   ) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) self%concentrationMinimum,self%concentrationMaximum,self%radiusMinimum,self%radiusMaximum,self%freefallRadiusMinimum,self%freefallRadiusMaximum
    ! Retabulate.
    self%nfwTableInitialized        =.false.
    self%nfwInverseTableInitialized =.false.
    self%nfwFreefallTableInitialized=.false.
    call self%tabulate              ()
    call self%inverseAngularMomentum()
    return
  end subroutine nfwStateRestore
