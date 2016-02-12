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

  !% An implementation of \cite{burkert_structure_1995} dark matter halo profiles.

  !# <darkMatterProfile name="darkMatterProfileBurkert">
  !#  <description>\cite{burkert_structure_1995} dark matter halo profiles</description>
  !# </darkMatterProfile>

  use Dark_Matter_Halo_Scales
  use Tables
  use Kind_Numbers

  type, extends(darkMatterProfileClass) :: darkMatterProfileBurkert
     !% A dark matter halo profile class implementing \cite{burkert_structure_1995} dark matter halos.
     private
     ! Minimum and maximum concentrations to tabulate.
     double precision                                        :: concentrationMinimum                   , concentrationMaximum
     ! Minimum and maximum radii to tabulate.
     double precision                                        :: freefallRadiusMinimum                  , radiusMinimum
     double precision                                        :: freefallRadiusMaximum                  , radiusMaximum
     double precision                                        :: freefallTimeMinimum                    , specificAngularMomentumMinimum
     double precision                                        :: freefallTimeMaximum                    , specificAngularMomentumMaximum
     ! Tables of Burkert properties.
     logical                                                 :: burkertFreefallTableInitialized        , burkertInverseTableInitialized       , &
          &                                                     burkertTableInitialized                
     integer                                                 :: burkertFreefallTableNumberPoints       , burkertInverseTableNumberPoints      , &
          &                                                     burkertTableNumberPoints
     type            (table1DLogarithmicLinear)              :: burkertConcentrationTable
     ! Tables.
     type            (table1DLogarithmicLinear)              :: burkertFreeFall                        , burkertSpecificAngularMomentum
     class           (table1D                 ), allocatable :: burkertFreefallInverse                 , burkertSpecificAngularMomentumInverse
     ! Module variables used in integrations.
     double precision                                        :: concentrationParameter                 , radiusStart
     ! Record of unique ID of node which we last computed results for.
     integer         (kind=kind_int8          )              :: lastUniqueID
     ! Record of whether or not quantities have been computed.
     logical                                                 :: specificAngularMomentumScalingsComputed, maximumVelocityComputed
     ! Stored values of computed quantities.
     double precision                                        :: specificAngularMomentumLengthScale     , specificAngularMomentumScale     , &
          &                                                     concentrationPrevious                  , burkertNormalizationFactorPrevious   , &
          &                                                     maximumVelocityPrevious
     ! Pointer to object setting halo scales.
     class(darkMatterHaloScaleClass           ), pointer     :: scale
   contains
     !@ <objectMethods>
     !@   <object>darkMatterProfileBurkert</object>
     !@   <objectMethod>
     !@     <method>densityScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin, \doublezero\ concentration\argin, \doublezero\ alpha\argin</arguments>
     !@     <description>Returns the density (in units such that the virial mass and scale length are unity) in an Burkert dark matter profile with given {\normalfont \ttfamily concentration} and {\normalfont \ttfamily alpha} at the given {\normalfont \ttfamily radius} (given in units of the scale radius).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>enclosedMassScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin, \doublezero\ concentration\argin, \doublezero\ alpha\argin</arguments>
     !@     <description>Returns the enclosed mass (in units of the virial mass) in an Burkert dark matter profile with given {\normalfont \ttfamily concentration} at the given {\normalfont \ttfamily radius} (given in units of the scale radius).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>freefallTabulate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ freefallTimeScaleFree\argin, \doublezero\ alphaRequired\argin</arguments>
     !@     <description>Tabulates the freefall time vs. freefall radius for Burkert halos.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>freefallTimeScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin, \doublezero\ alpha\argin</arguments>
     !@     <description>Compute the freefall time in a scale-free Burkert halo.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>angularMomentumScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ concentration\argin</arguments>
     !@     <description>Returns the total angular momentum in an Burkert dark matter profile with given {\normalfont \ttfamily concentration}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>inverseAngularMomentum</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ specificAngularMomentum\argin</arguments>
     !@     <description>Tabulates the specific angular momentum vs. radius in an Burkert profile for rapid inversion.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>profileEnergy</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ concentration\argin</arguments>
     !@     <description>Computes the total energy of an Burkert profile halo of given {\normalfont \ttfamily concentration}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>specificAngularMomentumScaleFree</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin</arguments>
     !@     <description>Returns the specific angular momentum, normalized to unit scale length and unit velocity at the scale radius, at position {\normalfont \ttfamily radius} (in units of the scale radius) in an Burkert profile.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>tabulate</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ concentration\argin</arguments>
     !@     <description>Tabulate properties of the Burkert halo profile which must be computed numerically.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final                                             burkertDestructor
     procedure :: calculationReset                  => burkertCalculationReset
     procedure :: stateStore                        => burkertStateStore
     procedure :: stateRestore                      => burkertStateRestore
     procedure :: density                           => burkertDensity
     procedure :: enclosedMass                      => burkertEnclosedMass
     procedure :: potential                         => burkertPotential
     procedure :: circularVelocity                  => burkertCircularVelocity
     procedure :: circularVelocityMaximum           => burkertCircularVelocityMaximum
     procedure :: radiusFromSpecificAngularMomentum => burkertRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => burkertRotationNormalization
     procedure :: energy                            => burkertEnergy
     procedure :: energyGrowthRate                  => burkertEnergyGrowthRate
     procedure :: kSpace                            => burkertKSpace
     procedure :: freefallRadius                    => burkertFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => burkertFreefallRadiusIncreaseRate
     procedure :: profileEnergy                     => burkertProfileEnergy
     procedure :: specificAngularMomentumScaleFree  => burkertSpecificAngularMomentumScaleFree
     procedure :: angularMomentumScaleFree          => burkertAngularMomentumScaleFree
     procedure :: enclosedMassScaleFree             => burkertEnclosedMassScaleFree
     procedure :: densityScaleFree                  => burkertDensityScaleFree
     procedure :: tabulate                          => burkertTabulate
     procedure :: inverseAngularMomentum            => burkertInverseAngularMomentum
     procedure :: freefallTabulate                  => burkertFreefallTabulate
     procedure :: freefallTimeScaleFree             => burkertFreefallTimeScaleFree
  end type darkMatterProfileBurkert

  interface darkMatterProfileBurkert
     !% Constructors for the {\normalfont \ttfamily burkert} dark matter halo profile class.
     module procedure burkertDefaultConstructor
     module procedure burkertConstructor
  end interface darkMatterProfileBurkert

  ! Number of points per decade of concentration in Burkert tabulations.
  integer, parameter   :: burkertTablePointsPerDecade        =100
  integer, parameter   :: burkertInverseTablePointsPerDecade =100
  integer, parameter   :: burkertFreefallTablePointsPerDecade=100
  ! Indices for tabulated quantities.
  integer, parameter   :: burkertConcentrationEnergyIndex    =  1, burkertConcetrationRotationNormalizationIndex=2

contains

  function burkertDefaultConstructor()
    !% Default constructor for the {\normalfont \ttfamily burkert} dark matter halo profile class.
    use Input_Parameters
    implicit none
    type(darkMatterProfileBurkert), target :: burkertDefaultConstructor

    burkertDefaultConstructor=burkertConstructor(darkMatterHaloScale())
    return
  end function burkertDefaultConstructor

  function burkertConstructor(scale)
    !% Generic constructor for the {\normalfont \ttfamily burkert} dark matter halo profile class.
    use Galacticus_Error
    implicit none
    type (darkMatterProfileBurkert    ), target :: burkertConstructor
    class(darkMatterHaloScaleClass), target :: scale
 
    burkertConstructor%concentrationPrevious       =  -1.0d+0
    burkertConstructor%concentrationMinimum        =   1.0d+0
    burkertConstructor%concentrationMaximum        =  20.0d+0
    burkertConstructor%freefallRadiusMinimum       =   1.0d-3 
    burkertConstructor%radiusMinimum               =   1.0d-3
    burkertConstructor%freefallRadiusMaximum       =   1.0d+2  
    burkertConstructor%radiusMaximum               =   1.0d+2
    burkertConstructor%burkertFreefallTableInitialized =  .false.
    burkertConstructor%burkertInverseTableInitialized  =  .false.
    burkertConstructor%burkertTableInitialized         =  .false.
    burkertConstructor%lastUniqueID                =  -1
    burkertConstructor%scale                       => scale
    ! Ensure that the dark matter profile component supports a "scale" property.
    if (.not.defaultDarkMatterProfileComponent%scaleIsGettable())                                                             &
         & call Galacticus_Error_Report                                                                                       &
         &      (                                                                                                             &
         &       'burkertConstructor'                                                                                       , &
         &       'Burkert dark matter profile requires a dark matter profile component with a gettable "scale" property.'//   &
         &       Galacticus_Component_List(                                                                                   &
         &                                 'darkMatterProfile'                                                              , &
         &                                 defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)      &
         &                                )                                                                                   &
         &      )
    ! Initialize the tabulations.
    call burkertConstructor%tabulate              ()
    call burkertConstructor%inverseAngularMomentum()
    return
  end function burkertConstructor
  
  subroutine burkertDestructor(self)
    !% Destructor for the {\normalfont \ttfamily burkert} dark matter halo profile class.
    implicit none
    type(darkMatterProfileBurkert), intent(inout) :: self

    if (self%burkertFreefallTableInitialized) then
       call self%burkertFreeFall                      %destroy()
       call self%burkertFreeFallInverse               %destroy()
       deallocate(self%burkertFreefallInverse)
    end if
    if (self%burkertInverseTableInitialized ) then
       call self%burkertSpecificAngularMomentum       %destroy()
       call self%burkertSpecificAngularMomentumInverse%destroy()
       deallocate(self%burkertSpecificAngularMomentumInverse)
    end if
    if (self%burkertTableInitialized        ) then
       call self%burkertConcentrationTable            %destroy()
    end if
    if (self%scale%isFinalizable()) deallocate(self%scale)
    return
  end subroutine burkertDestructor
  
  subroutine burkertCalculationReset(self,thisNode)
    !% Reset the dark matter profile calculation.
    implicit none
    class(darkMatterProfileBurkert), intent(inout)          :: self
    type (treeNode            ), intent(inout), pointer :: thisNode

    self%specificAngularMomentumScalingsComputed=.false.
    self%maximumVelocityComputed                =.false.
    self%lastUniqueID                           =thisNode%uniqueID()
    call self%scale%calculationReset(thisNode)
    return
  end subroutine burkertCalculationReset

  subroutine burkertTabulate(self,concentration)
    !% Tabulate properties of the Burkert halo profile which must be computed numerically.
    implicit none
    class           (darkMatterProfileBurkert), intent(inout)           :: self
    double precision                      , intent(in   ), optional :: concentration
    integer                                                         :: iConcentration
    logical                                                         :: retabulate
    double precision                                                :: tableConcentration

    retabulate=.not.self%burkertTableInitialized
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
       self%burkertTableNumberPoints=int(log10(self%concentrationMaximum/self%concentrationMinimum)*dble(burkertTablePointsPerDecade))+1
       call self%burkertConcentrationTable%destroy()
       call self%burkertConcentrationTable%create(self%concentrationMinimum,self%concentrationMaximum,self%burkertTableNumberPoints,2)
       ! Loop over concentrations and populate tables.
       do iConcentration=1,self%burkertTableNumberPoints
          tableConcentration=self%burkertConcentrationTable%x(iConcentration)
          call self%burkertConcentrationTable%populate(                   self%profileEnergy           (tableConcentration),iConcentration,table=burkertConcentrationEnergyIndex              )
          call self%burkertConcentrationTable%populate(tableConcentration/self%angularMomentumScaleFree(tableConcentration),iConcentration,table=burkertConcetrationRotationNormalizationIndex)
       end do
       ! Specify that tabulation has been made.
       self%burkertTableInitialized=.true.
    end if
    return
  end subroutine burkertTabulate

  subroutine burkertInverseAngularMomentum(self,specificAngularMomentum)
    !% Tabulates the specific angular momentum vs. radius in an Burkert profile for rapid inversion.
    implicit none
    class           (darkMatterProfileBurkert), intent(inout)           :: self
    double precision                      , intent(in   ), optional :: specificAngularMomentum
    integer                                                         :: iRadius
    logical                                                         :: retabulate

    retabulate=.not.self%burkertInverseTableInitialized
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
       self%burkertInverseTableNumberPoints=int(log10(self%radiusMaximum/self%radiusMinimum)*dble(burkertInverseTablePointsPerDecade))+1
       ! Create a range of radii.
       call self%burkertSpecificAngularMomentum%destroy(                                                       )
       call self%burkertSpecificAngularMomentum%create (self%radiusMinimum,self%radiusMaximum,self%burkertInverseTableNumberPoints)
       ! Loop over radii and populate tables.
       do iRadius=1,self%burkertInverseTableNumberPoints
          call self%burkertSpecificAngularMomentum%populate(                                                                                   &
               &                                        self%specificAngularMomentumScaleFree(self%burkertSpecificAngularMomentum%x(iRadius)), &
               &                                        iRadius                                                                            &
               &                                       )
       end do
       call self%burkertSpecificAngularMomentum%reverse(self%burkertSpecificAngularMomentumInverse)
       ! Specify that tabulation has been made.
       self%burkertInverseTableInitialized=.true.
    end if
    return
  end subroutine burkertInverseAngularMomentum

  double precision function burkertDensity(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given
    !% in units of Mpc).
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileBurkert      ), intent(inout)          :: self
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
    burkertDensity=self%densityScaleFree(radiusOverScaleRadius,virialRadiusOverScaleRadius)*thisBasicComponent%mass()/scaleRadius**3
    return
  end function burkertDensity

  double precision function burkertEnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileBurkert      ), intent(inout)          :: self
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
    burkertEnclosedMass            =self%enclosedMassScaleFree(radiusOverScaleRadius,virialRadiusOverScaleRadius)*thisBasicComponent%mass()
    return
  end function burkertEnclosedMass

  double precision function burkertPotential(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    use Galactic_Structure_Options
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Math
    implicit none
    class           (darkMatterProfileBurkert      ), intent(inout)           :: self
    type            (treeNode                      ), intent(inout), pointer  :: node
    double precision                                , intent(in   )           :: radius
    integer                                         , intent(  out), optional :: status
    class           (nodeComponentDarkMatterProfile)               , pointer  :: thisDarkMatterProfileComponent
    double precision                                , parameter               :: radiusSmall                   =1.0d-10
    double precision                                                          :: radiusOverScaleRadius                 , radiusTerm, &
         &                                                                       virialRadiusOverScaleRadius

    if (present(status)) status=structureErrorCodeSuccess
    thisDarkMatterProfileComponent   => node%darkMatterProfile(autoCreate=.true.)
    radiusOverScaleRadius            =radius                       /thisDarkMatterProfileComponent%scale()
    virialRadiusOverScaleRadius      =self%scale%virialRadius(node)/thisDarkMatterProfileComponent%scale()
    if (radiusOverScaleRadius < radiusSmall) then
       burkertPotential=                                                                          &
            & +(                                                                                  &
            &   -Pi                                                                               &
            &   +2.0d0                                                                            &
            &   /3.0d0                                                                            &
            &   *radiusOverScaleRadius**2                                                         &
            &  )                                                                                  &
            & /(                                                                                  &
            &     -2.0d0                              *atan(      virialRadiusOverScaleRadius   ) &
            &     +2.0d0                              *log (1.0d0+virialRadiusOverScaleRadius   ) &
            &     +                                    log (1.0d0+virialRadiusOverScaleRadius**2) &
            &  )                                                                                  &
            & *virialRadiusOverScaleRadius                                                        &
            & *self%scale%virialVelocity(node)**2
    else
       burkertPotential=                                                                          &
            & +(                                                                                  &
            &   -Pi                                                                               &
            &   *radiusOverScaleRadius                                                            &
            &   +(                                                                                &
            &     +2.0d0*(1.0d0+radiusOverScaleRadius)*atan(      radiusOverScaleRadius         ) &
            &     -2.0d0*(1.0d0+radiusOverScaleRadius)*log (1.0d0+radiusOverScaleRadius         ) &
            &     -      (1.0d0-radiusOverScaleRadius)*log (1.0d0+radiusOverScaleRadius      **2) &
            &    )                                                                                &
            &  )                                                                                  &
            & /radiusOverScaleRadius                                                              &
            & /(                                                                                  &
            &     -2.0d0                              *atan(      virialRadiusOverScaleRadius   ) &
            &     +2.0d0                              *log (1.0d0+virialRadiusOverScaleRadius   ) &
            &     +                                    log (1.0d0+virialRadiusOverScaleRadius**2) &
            &  )                                                                                  &
            & *virialRadiusOverScaleRadius                                                        &
            & *self%scale%virialVelocity(node)**2
    end if
    return
  end function burkertPotential

  double precision function burkertCircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    use Numerical_Constants_Physical
    implicit none
    class           (darkMatterProfileBurkert), intent(inout)          :: self
    type            (treeNode                ), intent(inout), pointer :: node
    double precision                          , intent(in   )          :: radius

    if (radius > 0.0d0) then
       burkertCircularVelocity=sqrt(gravitationalConstantGalacticus*self%enclosedMass(node,radius)/radius)
    else
       burkertCircularVelocity=0.0d0
    end if
    return
  end function burkertCircularVelocity

  double precision function burkertCircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    use Numerical_Constants_Physical
    implicit none
    class           (darkMatterProfileBurkert      ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    ! The radius (in units of the scale radius) at which the rotation speed peaks in a Burkert halo.
    double precision                                , parameter              :: radiusMaximum    =3.2446257246042642d0
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile
    double precision                                                         :: scaleRadius

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Check if maximum velocity is already computed. Compute and store if not.
    if (.not.self%maximumVelocityComputed) then
       darkMatterProfile            => node             %darkMatterProfile(autoCreate=.true.                          )
       scaleRadius                  =  darkMatterProfile%scale            (                                           )
       self%maximumVelocityPrevious =  self             %circularVelocity (node             ,radiusMaximum*scaleRadius)
       self%maximumVelocityComputed =  .true.
    end if
    burkertCircularVelocityMaximum=self%maximumVelocityPrevious
    return
  end function burkertCircularVelocityMaximum

  double precision function burkertRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily
    !% specificAngularMomentum} (given in units of km s$^{-1}$ Mpc)
    implicit none
    class           (darkMatterProfileBurkert      ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: specificAngularMomentum
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: specificAngularMomentumScaleFree

    ! Return immediately with zero radius for non-positive specific angular momenta.
    if (specificAngularMomentum <= 0.0d0) then
       burkertRadiusFromSpecificAngularMomentum=0.0d0
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

    ! Compute the specific angular momentum in scale free units (using the scale length for distances and sqrt(G M(r_scale) /
    ! r_scale) for velocities).
    specificAngularMomentumScaleFree=specificAngularMomentum/self%specificAngularMomentumScale

    ! Ensure that the interpolations exist and extend sufficiently far.
    call self%inverseAngularMomentum(specificAngularMomentumScaleFree)

    ! Interpolate to get the dimensionless radius at which this specific angular momentum is found.
    burkertRadiusFromSpecificAngularMomentum=self%burkertSpecificAngularMomentumInverse%interpolate(specificAngularMomentumScaleFree)

    ! Convert to a physical radius.
    burkertRadiusFromSpecificAngularMomentum=burkertRadiusFromSpecificAngularMomentum*self%specificAngularMomentumLengthScale
    return
  end function burkertRadiusFromSpecificAngularMomentum

  double precision function burkertRotationNormalization(self,node)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileBurkert          ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: concentration

    ! Get components.
    thisDarkMatterProfileComponent => node%darkMatterProfile(autoCreate=.true.)

    ! Find the concentration parameter of this halo.
    concentration=self%scale%virialRadius(node)/thisDarkMatterProfileComponent%scale()

    ! Ensure that the interpolations exist and extend sufficiently far.
    call self%tabulate(concentration)

    ! Find the rotation normalization by interpolation.
    burkertRotationNormalization=self%burkertConcentrationTable%interpolate(concentration,table&
         &=burkertConcetrationRotationNormalizationIndex)/self%scale%virialRadius(node)
    return
  end function burkertRotationNormalization

  double precision function burkertEnergy(self,node)
    !% Return the energy of an Burkert halo density profile.
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileBurkert      ), intent(inout)          :: self
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
    burkertEnergy=self%burkertConcentrationTable%interpolate(concentration,table=burkertConcentrationEnergyIndex)&
         &*thisBasicComponent%mass()*self%scale%virialVelocity(node)**2
    return
  end function burkertEnergy

  double precision function burkertEnergyGrowthRate(self,node)
    !% Return the rate of change of the energy of an Burkert halo density profile.
    use Dark_Matter_Halo_Scales
    implicit none
    class           (darkMatterProfileBurkert      ), intent(inout)          :: self
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
    energy        =self%burkertConcentrationTable%interpolate        (concentration,table=burkertConcentrationEnergyIndex)
    energyGradient=self%burkertConcentrationTable%interpolateGradient(concentration,table=burkertConcentrationEnergyIndex)

    burkertEnergyGrowthRate=self%energy(node)&
         &*(thisBasicComponent%accretionRate()/thisBasicComponent%mass()+2.0d0 &
         &*self%scale%virialVelocityGrowthRate(node)/self%scale%virialVelocity(node)+(energyGradient&
         &*concentration/energy)*(self%scale%virialRadiusGrowthRate(node)/self%scale%virialRadius(node)&
         &-thisDarkMatterProfileComponent%scaleGrowthRate()/thisDarkMatterProfileComponent%scale()))

    return
  end function burkertEnergyGrowthRate

  double precision function burkertAngularMomentumScaleFree(self,concentration)
    !% Returns the total angular momentum (in units of the virial mass times scale radius times [assumed constant] rotation speed)
    !% in an Burkert dark matter profile with given {\normalfont \ttfamily concentration}. This is given by:
    !% \begin{equation}
    !% J = \left. \int_0^c 4 \pi x^3 \rho(x) \d x \right/ \int_0^c 4 \pi x^2 \rho(x) \d x,
    !% \end{equation}
    !% where $x$ is radius in units of the scale radius and $c$ is concentration. This can be evaluated to give
    !% \begin{equation}
    !% J = \left. \left[ 2 \tan^{-1} c + 2 \log(1+c) + \log(1+c^2) - 4c \right] \right/ \left[ 2 \tan^{-1} c - 2 \log(1+c) - \log(1+c^2) \right].
    !% \end{equation}
    implicit none
    class           (darkMatterProfileBurkert), intent(inout) :: self
    double precision                          , intent(in   ) :: concentration

    burkertAngularMomentumScaleFree=+(                                    &
         &                            +2.0d0*atan(      concentration   ) &
         &                            +2.0d0*log (1.0d0+concentration   ) &
         &                            +      log (1.0d0+concentration**2) &
         &                            -4.0d0*           concentration     &
         &                           )                                    &
         &                          /(                                    &
         &                            +2.0d0*atan(      concentration   ) &
         &                            -2.0d0*log (1.0d0+concentration   ) &
         &                            -      log (1.0d0+concentration**2) &
         &                           )
    return
  end function burkertAngularMomentumScaleFree

  double precision function burkertSpecificAngularMomentumScaleFree(self,radius)
    !% Returns the specific angular momentum, normalized to unit scale length and unit velocity at the scale radius, at position
    !% {\normalfont \ttfamily radius} (in units of the scale radius) in an Burkert profile.
    implicit none
    class           (darkMatterProfileBurkert), intent(inout) :: self
    double precision                      , intent(in   ) :: radius

    burkertSpecificAngularMomentumScaleFree=sqrt(radius*self%enclosedMassScaleFree(radius,1.0d0))
    return
  end function burkertSpecificAngularMomentumScaleFree

  double precision function burkertEnclosedMassScaleFree(self,radius,concentration)
    !% Returns the enclosed mass (in units of the virial mass) in an Burkert dark matter profile with given {\normalfont \ttfamily
    !% concentration} at the given {\normalfont \ttfamily radius} (given in units of the scale radius).
    implicit none
    class           (darkMatterProfileBurkert), intent(inout) :: self
    double precision                      , intent(in   ) :: concentration                                                   , radius
    double precision                      , parameter     :: minimumRadiusForExactSolution          =1.0d-4
    ! Precomputed Burkert normalization factor for unit concentration.
    double precision                      , parameter     :: burkertNormalizationFactorUnitConcentration=1.0d0/(3.0d0*log(2.0d0)-2.0d0*atan(1.0d0))

    if (radius < minimumRadiusForExactSolution) then
       ! Use a series solution for small radii.
       burkertEnclosedMassScaleFree=(4.0d0/3.0d0)*radius**3*(1.0d0-radius)
    else
       ! Use the exact solution.
       burkertEnclosedMassScaleFree=(                             &
            &                        +2.0d0*log (1.0d0+radius   ) &
            &                        +      log (1.0d0+radius**2) &
            &                        -2.0d0*atan(      radius   ) &
            &                       )
    end if
    ! Check if we were called with a different concentration compared to the previous call.
    if (concentration /= self%concentrationPrevious) then
       ! We were, so recompute the normalization factor.
       if (concentration == 1.0d0) then
          self%burkertNormalizationFactorPrevious=burkertNormalizationFactorUnitConcentration
       else
          self%burkertNormalizationFactorPrevious=+1.0d0                                &
               &                                  /(                                    &
               &                                    +2.0d0*log (1.0d0+concentration   ) &
               &                                    +      log (1.0d0+concentration**2) &
               &                                    -2.0d0*atan(      concentration   ) &
               &                                   )
       end if
       self%concentrationPrevious=concentration
    end if
    burkertEnclosedMassScaleFree=burkertEnclosedMassScaleFree*self%burkertNormalizationFactorPrevious
    return
  end function burkertEnclosedMassScaleFree

  double precision function burkertDensityScaleFree(self,radius,concentration)
    !% Returns the density (in units such that the virial mass and scale length are unity) in an Burkert dark matter profile with
    !% given {\normalfont \ttfamily concentration} at the given {\normalfont \ttfamily radius} (given in units of the scale radius).
    use Numerical_Constants_Math
    implicit none
    class           (darkMatterProfileBurkert), intent(inout) :: self
    double precision                          , intent(in   ) :: concentration, radius

    burkertDensityScaleFree=+1.0d0                                &
         &                  /(1.0d0+radius   )                    &
         &                  /(1.0d0+radius**2)                    &
         &                  /Pi                                   &
         &                  /(                                    &
         &                    +2.0d0*log (1.0d0+concentration   ) &
         &                    +      log (1.0d0+concentration**2) &
         &                    -2.0d0*atan(      concentration   ) &
         &                   )
    return
  end function burkertDensityScaleFree

  double precision function burkertProfileEnergy(self,concentration)
    !% Computes the total energy of an Burkert profile halo of given {\normalfont \ttfamily concentration} using the methods of
    !% \citeauthor{cole_hierarchical_2000}~(\citeyear{cole_hierarchical_2000}; their Appendix~A).
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Math
    use Numerical_Integration
    implicit none
    class           (darkMatterProfileBurkert      ), intent(inout) :: self
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
    potentialEnergyIntegral=Integrate(radiusMinimum,radiusMaximum,burkertPotentialEnergyIntegrand&
         &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    potentialEnergy=-0.5d0*(1.0d0/concentration+potentialEnergyIntegral)
    ! Compute the velocity dispersion at the virial radius.
    radiusMinimum=concentration
    radiusMaximum=100.0d0*concentration
    concentrationParameter=concentration
    jeansEquationIntegral=Integrate(radiusMinimum,radiusMaximum,burkertJeansEquationIntegrand&
         &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    ! Compute the kinetic energy.
    radiusMinimum=0.0d0
    radiusMaximum=concentration
    concentrationParameter=concentration
    kineticEnergyIntegral=Integrate(radiusMinimum,radiusMaximum,burkertKineticEnergyIntegrand&
         &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    kineticEnergy=2.0d0*Pi*(jeansEquationIntegral*concentration**3+kineticEnergyIntegral)
    ! Compute the total energy.
    burkertProfileEnergy=(potentialEnergy+kineticEnergy)*concentration
    return

  contains
    
    function burkertPotentialEnergyIntegrand(radius,parameterPointer) bind(c)
      !% Integrand for Burkert profile potential energy.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real(kind=c_double)        :: burkertPotentialEnergyIntegrand
      real(kind=c_double), value :: radius
      type(c_ptr        ), value :: parameterPointer
      
      burkertPotentialEnergyIntegrand=(self%enclosedMassScaleFree(radius,concentrationParameter)/radius)**2
      return
    end function burkertPotentialEnergyIntegrand
    
    function burkertKineticEnergyIntegrand(radius,parameterPointer) bind(c)
      !% Integrand for Burkert profile kinetic energy.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real(kind=c_double)        :: burkertKineticEnergyIntegrand
      real(kind=c_double), value :: radius
      type(c_ptr        ), value :: parameterPointer
      
      burkertKineticEnergyIntegrand=self%EnclosedMassScaleFree(radius,concentrationParameter)*self%densityScaleFree(radius&
           &,concentrationParameter)*radius
      return
    end function burkertKineticEnergyIntegrand
    
    function burkertJeansEquationIntegrand(radius,parameterPointer) bind(c)
      !% Integrand for Burkert profile Jeans equation.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real(kind=c_double)        :: burkertJeansEquationIntegrand
      real(kind=c_double), value :: radius
      type(c_ptr        ), value :: parameterPointer
      
      burkertJeansEquationIntegrand=self%enclosedMassScaleFree(radius,concentrationParameter)*self%densityScaleFree(radius &
           &,concentrationParameter)/radius**2
      return
    end function burkertJeansEquationIntegrand

  end function burkertProfileEnergy
  
  double precision function burkertKSpace(self,node,waveNumber)
    !% Returns the Fourier transform of the Burkert density profile at the specified {\normalfont \ttfamily waveNumber} (given in Mpc$^{-1}$), using the
    !% expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    use Dark_Matter_Halo_Scales
    use Exponential_Integrals
    use Numerical_Constants_Math
    implicit none
    class           (darkMatterProfileBurkert      ), intent(inout)          :: self
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
    burkertKSpace=dimag(                                                                                                                                              &
         &              +(                                                                                                                                            &
         &                +       2.0d0         *exp(-dcmplx(0.0d0,1.0d0)*waveNumberScaleFree)*Cosine_Integral     (                             waveNumberScaleFree) &
         &                -       2.0d0         *exp(-dcmplx(0.0d0,1.0d0)*waveNumberScaleFree)*Cosine_Integral     (      (+1.0d0+concentration)*waveNumberScaleFree) &
         &                +dcmplx(1.0d0,1.0d0)                                                                                                                        &
         &                *(                                                                                                                                          &
         &                  -dcmplx(0.0d0,1.0d0)*exp(+                    waveNumberScaleFree)*Pi                                                                     &
         &                  -                    exp(+                    waveNumberScaleFree)*Exponential_Integral(       -1.0d0               *waveNumberScaleFree) & 
         &                  +dcmplx(0.0d0,1.0d0)*exp(-                    waveNumberScaleFree)*Exponential_Integral(       +1.0d0               *waveNumberScaleFree) &
         &                  +                    exp(+                    waveNumberScaleFree)*Exponential_Integral(dcmplx(-1.0d0,concentration)*waveNumberScaleFree) &
         &                  -dcmplx(0.0d0,1.0d0)*exp(-                    waveNumberScaleFree)*Exponential_Integral(dcmplx(+1.0d0,concentration)*waveNumberScaleFree) &
         &                  +dcmplx(1.0d0,1.0d0)*exp(-dcmplx(0.0d0,1.0d0)*waveNumberScaleFree)*Sine_Integral       (                             waveNumberScaleFree) &
         &                  -dcmplx(1.0d0,1.0d0)*exp(-dcmplx(0.0d0,1.0d0)*waveNumberScaleFree)*Sine_Integral       (      (+1.0d0+concentration)*waveNumberScaleFree) &
         &                 )                                                                                                                                          &
         &               )                                                                                                                                            &
         &              /(                                                                                                                                            &
         &                +waveNumberScaleFree                                                                                                                        &
         &                *(                                                                                                                                          &
         &                  -2.0d0*atan(      concentration   )                                                                                                       &
         &                  +2.0d0*log (1.0d0+concentration   )                                                                                                       &
         &                  +      log (1.0d0+concentration**2)                                                                                                       &
         &                 )                                                                                                                                          &
         &               )                                                                                                                                            &
         &             )
    return
  end function burkertKSpace

  double precision function burkertFreefallRadius(self,node,time)
    !% Returns the freefall radius in the Burkert density profile at the specified {\normalfont \ttfamily time} (given in Gyr).
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    implicit none
    class           (darkMatterProfileBurkert      ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: time
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: concentration                 , freefallTimeScaleFree, &
         &                                                                      radiusScale                   , timeScale            , &
         &                                                                      velocityScale

    ! For non-positive freefall times, return a zero freefall radius immediately.
    if (time <= 0.0d0) then
       burkertFreefallRadius=0.0d0
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
    burkertFreefallRadius=self%burkertFreefallInverse%interpolate(freefallTimeScaleFree)*radiusScale
    return
  end function burkertFreefallRadius

  double precision function burkertFreefallRadiusIncreaseRate(self,node,time)
    !% Returns the rate of increase of the freefall radius in the Burkert density profile at the specified {\normalfont \ttfamily time} (given in
    !% Gyr).
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    implicit none
    class           (darkMatterProfileBurkert      ), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: time
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfileComponent
    double precision                                                         :: concentration                 , freefallTimeScaleFree, &
         &                                                                      radiusScale                   , timeScale            , &
         &                                                                      velocityScale

    ! For non-positive freefall times, return the limiting value for small radii.
    if (time <= 0.0d0) then
       burkertFreefallRadiusIncreaseRate=0.0d0
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
    burkertFreefallRadiusIncreaseRate=self%burkertFreefallInverse%interpolateGradient(freefallTimeScaleFree)*radiusScale/timeScale
    return
  end function burkertFreefallRadiusIncreaseRate

  subroutine burkertFreefallTabulate(self,freefallTimeScaleFree)
    !% Tabulates the freefall time vs. freefall radius for Burkert halos.
    implicit none
    class           (darkMatterProfileBurkert), intent(inout) :: self
    double precision                      , intent(in   ) :: freefallTimeScaleFree
    logical                                               :: retabulate
    integer                                               :: iRadius

    retabulate=.not.self%burkertFreefallTableInitialized
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
       self%burkertFreefallTableNumberPoints=int(log10(self%freefallRadiusMaximum/self%freefallRadiusMinimum)*dble(burkertFreefallTablePointsPerDecade))+1
       ! Create the table.
       call self%burkertFreefall%destroy(                                                                        )
       call self%burkertFreefall%create (self%freefallRadiusMinimum,self%freefallRadiusMaximum,self%burkertFreefallTableNumberPoints)
       ! Loop over radii and populate tables.
       do iRadius=1,self%burkertFreefallTableNumberPoints
          call self%burkertFreefall%populate(                                                         &
               &                         self%freefallTimeScaleFree(self%burkertFreefall%x(iRadius)), &
               &                                                                       iRadius    &
               &                        )
       end do
       call self%burkertFreefall%reverse(self%burkertFreefallInverse)
       ! Specify that tabulation has been made.
       self%burkertFreefallTableInitialized=.true.
    end if
    return
  end subroutine burkertFreefallTabulate

  double precision function burkertFreefallTimeScaleFree(self,radius)
    !% Compute the freefall time in a scale-free Burkert halo.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    implicit none
    class           (darkMatterProfileBurkert      ), intent(inout) :: self
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
       burkertFreefallTimeScaleFree=Integrate(radiusEnd,radiusStart,burkertFreefallTimeScaleFreeIntegrand,parameterPointer&
            &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
       call Integrate_Done(integrandFunction,integrationWorkspace)
    else
       ! Use an approximation here, found by taking series expansions of the logarithms in the integrand and keeping only the
       ! first order terms.
       burkertFreefallTimeScaleFree=2.0d0*sqrt(radius)
    end if
    return

  contains
    
    function burkertFreefallTimeScaleFreeIntegrand(radius,parameterPointer) bind(c)
      !% Integrand function used for finding the free-fall time in Burkert halos.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real(kind=c_double)                   :: burkertFreefallTimeScaleFreeIntegrand
      real(kind=c_double)           , value :: radius
      type(c_ptr        )           , value :: parameterPointer
      real(kind=c_double), parameter        :: radiusSmall                       =1.0d-6
      real(kind=c_double), parameter        :: radiusSmallFraction               =1.0d-3
      real(kind=c_double)                   :: x
      
      if (radius < radiusSmall) then
         ! Use a series approximation for small radii.
         burkertFreefallTimeScaleFreeIntegrand=log(1.0d0+radiusStart)/radiusStart-1.0d0+radius*(0.5d0-radius/3.0d0)
      else if (radius > radiusStart*(1.0d0-radiusSmallFraction)) then
         ! Use a series approximation for radii close to the initial radius.
         x=1.0d0-radius/radiusStart
         burkertFreefallTimeScaleFreeIntegrand=(1.0d0/(1.0d0+radiusStart)-log(1.0d0+radiusStart)/radiusStart)*x+(0.5d0*radiusStart&
              &/(1.0d0+radiusStart)**2+(radiusStart-(1.0d0+radiusStart)*log(1.0d0+radiusStart))/radiusStart/(1.0d0+radiusStart))*x&
              &**2
      else
         ! Use full expression for larger radii.
         burkertFreefallTimeScaleFreeIntegrand=log(1.0d0+radiusStart)/radiusStart-log(1.0d0+radius)/radius
      end if
      burkertFreefallTimeScaleFreeIntegrand=1.0d0/sqrt(-2.0d0*burkertFreefallTimeScaleFreeIntegrand)
      return
    end function burkertFreefallTimeScaleFreeIntegrand

  end function burkertFreefallTimeScaleFree

  subroutine burkertStateStore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    implicit none
    class  (darkMatterProfileBurkert), intent(inout) :: self
    integer                      , intent(in   ) :: stateFile
    type   (fgsl_file           ), intent(in   ) :: fgslStateFile

    write (stateFile) self%concentrationMinimum,self%concentrationMaximum,self%radiusMinimum,self%radiusMaximum,self%freefallRadiusMinimum,self%freefallRadiusMaximum
    return
  end subroutine burkertStateStore

  subroutine burkertStateRestore(self,stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    implicit none
    class  (darkMatterProfileBurkert), intent(inout) :: self
    integer                      , intent(in   ) :: stateFile
    type   (fgsl_file           ), intent(in   ) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) self%concentrationMinimum,self%concentrationMaximum,self%radiusMinimum,self%radiusMaximum,self%freefallRadiusMinimum,self%freefallRadiusMaximum
    ! Retabulate.
    self%burkertTableInitialized        =.false.
    self%burkertInverseTableInitialized =.false.
    self%burkertFreefallTableInitialized=.false.
    call self%tabulate              ()
    call self%inverseAngularMomentum()
    return
  end subroutine burkertStateRestore
