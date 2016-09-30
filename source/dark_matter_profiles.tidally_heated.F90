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

  !% An implementation of tidally heated dark matter halo profiles.

  !# <darkMatterProfile name="darkMatterProfileTidallyHeated">
  !#  <description>Tidally heated dark matter halo profiles.</description>
  !# </darkMatterProfile>

  type, extends(darkMatterProfileClass) :: darkMatterProfileTidallyHeated
     !% A dark matter halo profile class implementing tidally heated dark matter halos.
     private
     class(darkMatterProfileClass), pointer :: unheatedProfile
     logical                                :: unimplementedFatal
   contains
     !@ <objectMethods>
     !@   <object>darkMatterProfileTidallyHeated</object>
     !@   <objectMethod>
     !@     <method>radiusInitial</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless *type(treeNode)\textgreater} node\arginout, \doublezero\ radiusFinal\argin</arguments>
     !@     <description>Return the initial radius corresponding to the given final radius in a tidally-heated dark matter halo density profile.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final                                             tidallyHeatedDestructor
     procedure :: calculationReset                  => tidallyHeatedCalculationReset
     procedure :: stateStore                        => tidallyHeatedStateStore
     procedure :: stateRestore                      => tidallyHeatedStateRestore
     procedure :: radiusInitial                     => tidallyHeatedRadiusInitial
     procedure :: density                           => tidallyHeatedDensity
     procedure :: densityLogSlope                   => tidallyHeatedDensityLogSlope
     procedure :: radialMoment                      => tidallyHeatedRadialMoment
     procedure :: enclosedMass                      => tidallyHeatedEnclosedMass
     procedure :: potential                         => tidallyHeatedPotential
     procedure :: circularVelocity                  => tidallyHeatedCircularVelocity
     procedure :: circularVelocityMaximum           => tidallyHeatedCircularVelocityMaximum
     procedure :: radiusFromSpecificAngularMomentum => tidallyHeatedRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => tidallyHeatedRotationNormalization
     procedure :: energy                            => tidallyHeatedEnergy
     procedure :: energyGrowthRate                  => tidallyHeatedEnergyGrowthRate
     procedure :: kSpace                            => tidallyHeatedKSpace
     procedure :: freefallRadius                    => tidallyHeatedFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => tidallyHeatedFreefallRadiusIncreaseRate
  end type darkMatterProfileTidallyHeated

  interface darkMatterProfileTidallyHeated
     !% Constructors for the {\normalfont \ttfamily tidallyHeated} dark matter halo profile class.
     module procedure tidallyHeatedDefaultConstructor
  end interface darkMatterProfileTidallyHeated

  ! Record of whether method is initialized.
  logical                                                   :: tidallyHeatedInitialized                            =.false.

  ! Default settings for this method.
  type            (varying_string                )          :: darkMatterProfileTidallyHeatedUnheatedProfile  
  logical                                                   :: darkMatterProfileTidallyHeatedUnimplementedIsFatal
  
  ! Warnings issued.
  logical                                                   :: tidallyHeatedDensityLogSlopeWarned                  =.false.
  logical                                                   :: tidallyHeatedRadialMomentWarned                     =.false.
  logical                                                   :: tidallyHeatedCircularVelocityMaximumWarned          =.false.
  logical                                                   :: tidallyHeatedRotationNormalizationWarned            =.false.
  logical                                                   :: tidallyHeatedEnergyWarned                           =.false.
  logical                                                   :: tidallyHeatedEnergyGrowthRateWarned                 =.false.
  logical                                                   :: tidallyHeatedKSpaceWarned                           =.false.
  logical                                                   :: tidallyHeatedFreefallRadiusWarned                   =.false.
  logical                                                   :: tidallyHeatedFreefallRadiusIncreaseRateWarned       =.false.

  ! Global variables used in root solving.
  double precision                                          :: tidallyHeatedHeatingNormalized      , tidallyHeatedRadiusFinal, &
       &                                                       tidallyHeatedSpecificAngularMomentum
  type            (treeNode                      ), pointer :: tidallyHeatedNode
  type            (darkMatterProfileTidallyHeated), pointer :: tidallyHeatedSelf
  !$omp threadprivate(tidallyHeatedHeatingNormalized,tidallyHeatedRadiusFinal,tidallyHeatedNode,tidallyHeatedSelf,tidallyHeatedSpecificAngularMomentum)

contains

  function tidallyHeatedDefaultConstructor()
    !% Default constructor for the {\normalfont \ttfamily tidallyHeated} dark matter halo profile class.
    use Galacticus_Error
    use Input_Parameters
    implicit none
    type(darkMatterProfileTidallyHeated), target :: tidallyHeatedDefaultConstructor    

    ! Initialize if necessary.
    if (.not.tidallyHeatedInitialized) then
       !$omp critical(tidallyHeatedInitialization)
       if (.not.tidallyHeatedInitialized) then
          ! Read parameter controlling the unheated profile
          !@ <inputParameter>
          !@   <name>darkMatterProfileTidallyHeatedUnheatedProfile</name>
          !@   <defaultValue>NFW</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The unheated {\normalfont \ttfamily darkMatterProfile} method to which tidal heating should be applied.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>starFormation</group>
          !@ </inputParameter>
          call Get_Input_Parameter('darkMatterProfileTidallyHeatedUnheatedProfile',darkMatterProfileTidallyHeatedUnheatedProfile,defaultValue="NFW")
          !@ <inputParameter>
          !@   <name>darkMatterProfileTidallyHeatedUnimplementedIsFatal</name>
          !@   <defaultValue>true</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     If {\normalfont \ttfamily true}, unimplemented features of the tidally-heated dark matter profile cause fatal errors. Otherwise, a warning is issued and the solution for an unheated profile is returned.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>starFormation</group>
          !@ </inputParameter>
          call Get_Input_Parameter('darkMatterProfileTidallyHeatedUnimplementedIsFatal',darkMatterProfileTidallyHeatedUnimplementedIsFatal,defaultValue=.true.)
          ! Check that tidal heating is available.
          if (.not.defaultSatelliteComponent%tidalHeatingNormalizedIsGettable())                                                   &
            & call Galacticus_Error_Report                                                                                         &
            & (                                                                                                                    &
            &  'tidallyHeatedDefaultConstructor'                                                                                 , &
            &  'This method requires that the "tidalHeatingNormalized" property of the satellite are gettable.'                 // &
            &  Galacticus_Component_List(                                                                                          &
            &                            'satellite'                                                                            ,  &
            &                             defaultSatelliteComponent%tidalHeatingNormalizedAttributeMatch(requireGettable=.true.)   &
            &                           )                                                                                          &
            & )
          ! Record that this method is now initialized
          tidallyHeatedInitialized=.true.
       end if
       !$omp end critical(tidallyHeatedInitialization)
    end if
    ! Construct the default object.
    tidallyHeatedDefaultConstructor%unheatedProfile   => darkMatterProfile(char(darkMatterProfileTidallyHeatedUnheatedProfile))
    tidallyHeatedDefaultConstructor%unimplementedFatal=  darkMatterProfileTidallyHeatedUnimplementedIsFatal
    return
  end function tidallyHeatedDefaultConstructor

  function tidallyHeatedGenericConstructor(unheatedProfile)
    !% Generic constructor for the {\normalfont \ttfamily tidallyHeated} dark matter profile class.
    implicit none
    type (darkMatterProfileTidallyHeated)                        :: tidallyHeatedGenericConstructor
    class(darkMatterProfileClass        ), intent(in   ), target :: unheatedProfile

    ! Construct the object.
    tidallyHeatedGenericConstructor%unheatedProfile => unheatedProfile
    return
  end function tidallyHeatedGenericConstructor

  subroutine tidallyHeatedDestructor(self)
    !% Destructor for the {\normalfont \ttfamily tidallyHeated} dark matter halo profile class.
    implicit none
    type(darkMatterProfileTidallyHeated), intent(inout) :: self

    if (self%unheatedProfile%isFinalizable()) deallocate(self%unheatedProfile)
    return
  end subroutine tidallyHeatedDestructor
  
  subroutine tidallyHeatedCalculationReset(self,thisNode)
    !% Reset the dark matter profile calculation.
    implicit none
    class(darkMatterProfileTidallyHeated), intent(inout)          :: self
    type (treeNode                      ), intent(inout) :: thisNode

    call self%unheatedProfile%calculationReset(thisNode)
    return
  end subroutine tidallyHeatedCalculationReset

  subroutine tidallyHeatedSetGlobalSelf(self)
    !% Set a module-scope pointer to ``{\normalfont \ttfamily self}''.
    implicit none
    class(darkMatterProfileTidallyHeated), intent(inout), target :: self

    tidallyHeatedSelf => self
    return
  end subroutine tidallyHeatedSetGlobalSelf
  
  double precision function tidallyHeatedDensity(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Display
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    implicit none
    class           (darkMatterProfileTidallyHeated), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    class           (nodeComponentSatellite        ), pointer       :: satellite
    double precision                                                :: radiusInitial , massEnclosed, &
         &                                                             densityInitial, jacobian
    
    satellite => node%satellite()
    if (satellite%tidalHeatingNormalized() <= 0.0d0) then
       ! No tidal heating has occurred - fall through to the unheated profile.
       tidallyHeatedDensity=self%unheatedProfile%density(node,radius)
    else
       radiusInitial      =self                %radiusInitial(node,radius       )
       massEnclosed       =self%unheatedProfile%enclosedMass (node,radiusInitial)
       densityInitial     =self%unheatedProfile%density      (node,radiusInitial)
       jacobian           =+1.0d0                                   &
            &              /(                                       &
            &                +(                                     &
            &                  +radius                              &
            &                  /radiusInitial                       &
            &                 )                                 **2 &
            &                +4.0d0                                 &
            &                *satellite%tidalHeatingNormalized()    &
            &                *radiusInitial                         &
            &                *radius                            **2 &
            &                /gravitationalConstantGalacticus       &
            &                /massEnclosed                          &
            &                *(                                     &
            &                  +1.0d0                               &
            &                  -2.0d0                               &
            &                  *Pi                                  &
            &                  *radiusInitial                   **3 &
            &                  *densityInitial                      &
            &                  /massEnclosed                        &
            &                 )                                     &
            &               )
       tidallyHeatedDensity=+densityInitial                         &
            &               *(                                      &
            &                 +radiusInitial                        &
            &                 /radius                               &
            &                )                                  **2 &
            &               *jacobian       
    end if
    return
  end function tidallyHeatedDensity

  double precision function tidallyHeatedDensityLogSlope(self,node,radius)
    !% Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTidallyHeated), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: radius
    class           (nodeComponentSatellite        )               , pointer :: satellite

    satellite => node%satellite()
    if (satellite%tidalHeatingNormalized() <= 0.0d0) then
       ! No tidal heating has occurred - fall through to the unheated profile.
       tidallyHeatedDensityLogSlope=self%unheatedProfile%densityLogSlope(node,radius)
    else
       ! Check if unimplemented features are fatal.
       if (self%unimplementedFatal) then
          call Galacticus_Error_Report('tidallyHeatedDensityLogSlope','density logarithmic slope in tidally heated dark matter profiles is not supported')
       else
          ! Issue a warning and then fall through to the unheated profile.
          if (.not.tidallyHeatedDensityLogSlopeWarned) then
             !$omp critical(tidallyHeatedDarkMatterProfileDensityLogSlopeWarn)
             if (.not.tidallyHeatedDensityLogSlopeWarned) then
                call Galacticus_Display_Message('WARNING: density logarithmic slope in tidally heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
                tidallyHeatedDensityLogSlopeWarned=.true.
             end if
             !$omp end critical(tidallyHeatedDarkMatterProfileDensityLogSlopeWarn)
          end if
       end if
       tidallyHeatedDensityLogSlope=self%unheatedProfile%densityLogSlope(node,radius)
    end if
    return
  end function tidallyHeatedDensityLogSlope

  double precision function tidallyHeatedRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Error
    use Galacticus_Display
    implicit none
    class           (darkMatterProfileTidallyHeated), intent(inout)           :: self
    type            (treeNode                      ), intent(inout), pointer  :: node
    double precision                                , intent(in   )           :: moment
    double precision                                , intent(in   ), optional :: radiusMinimum                 , radiusMaximum

    class           (nodeComponentSatellite        )               , pointer :: satellite

    satellite => node%satellite()
    if (satellite%tidalHeatingNormalized() <= 0.0d0) then
       ! No tidal heating has occurred - fall through to the unheated profile.
       tidallyHeatedRadialMoment=self%unheatedProfile%radialMoment(node,moment,radiusMinimum,radiusMaximum)
    else
       ! Check if unimplemented features are fatal.
       if (self%unimplementedFatal) then
          tidallyHeatedRadialMoment=0.0d0
          call Galacticus_Error_Report('tidallyHeatedRadialMoment','radial moment in tidally heated dark matter profiles is not supported')
       else
          ! Issue a warning and then fall through to the unheated profile.
          if (.not.tidallyHeatedRadialMomentWarned) then
             !$omp critical(tidallyHeatedDarkMatterProfileRadialMomentWarn)
             if (.not.tidallyHeatedRadialMomentWarned) then
                call Galacticus_Display_Message('WARNING: radial moment in tidally heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
                tidallyHeatedRadialMomentWarned=.true.
             end if
             !$omp end critical(tidallyHeatedDarkMatterProfileRadialMomentWarn)
          end if
          tidallyHeatedRadialMoment=self%unheatedProfile%radialMoment(node,moment,radiusMinimum,radiusMaximum)
       end if
    end if
    return 
  end function tidallyHeatedRadialMoment

  double precision function tidallyHeatedEnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    use Galacticus_Display
    implicit none
    class           (darkMatterProfileTidallyHeated), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: radius
    class           (nodeComponentSatellite        ), pointer       :: satellite

    satellite => node%satellite()
    if (satellite%tidalHeatingNormalized() <= 0.0d0) then
       ! No tidal heating has occurred - fall through to the unheated profile.
       tidallyHeatedEnclosedMass=self%unheatedProfile%enclosedMass(node,                        radius )
    else
       tidallyHeatedEnclosedMass=self%unheatedProfile%enclosedMass(node,self%radiusInitial(node,radius))
    end if
    return
  end function tidallyHeatedEnclosedMass

  double precision function tidallyHeatedRadiusInitial(self,node,radiusFinal)
    !% Find the initial radius corresponding to the given {\normalfont \ttfamily radiusFinal} in
    !% the tidally-heated dark matter profile.
    use Root_Finder
    implicit none
    class           (darkMatterProfileTidallyHeated), intent(inout), target  :: self
    type            (treeNode                      ), intent(inout), target  :: node
    double precision                                , intent(in   )          :: radiusFinal
    class           (nodeComponentSatellite        )               , pointer :: satellite
    double precision                                , parameter              :: toleranceAbsolute     =0.0d0, toleranceRelative=1.0d-6
    type            (rootFinder                    ), save                   :: finder
    !$omp threadprivate(finder)
    
    ! Get the degree of tidal heating.
    tidallyHeatedSelf             =>           self
    tidallyHeatedNode             =>           node
    satellite                     =>           node     %satellite             ()
    tidallyHeatedHeatingNormalized=  max(0.0d0,satellite%tidalHeatingNormalized())
    tidallyHeatedRadiusFinal      =  radiusFinal
    ! Initialize the root finder.
    if (.not.finder%isInitialized()) then
       call finder%rootFunction(tidallyHeatedRadiusInitialRoot     )
       call finder%tolerance   (toleranceAbsolute,toleranceRelative)
       call finder%rangeExpand (                                                             &
            &                   rangeExpandUpward            =1.0d0                        , &
            &                   rangeExpandDownward          =0.5d0                        , &
            &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                   rangeExpandType              =rangeExpandMultiplicative      &
            &                  )
    end if
    tidallyHeatedRadiusInitial=finder%find(rootGuess=radiusFinal)
    return
  end function tidallyHeatedRadiusInitial
  
  double precision function tidallyHeatedRadiusInitialRoot(radiusInitial)
    !% Root function used in finding initial radii in tidally-heated dark matter halo profiles.
    use Numerical_Constants_Physical
    implicit none
    double precision, intent(in   ) :: radiusInitial
    double precision                :: massEnclosed
    
    massEnclosed                  =tidallyHeatedSelf%unheatedProfile%enclosedMass(tidallyHeatedNode,radiusInitial)
    tidallyHeatedRadiusInitialRoot=+tidallyHeatedHeatingNormalized     &
         &                         *radiusInitial                  **2 &
         &                         +0.5d0                              &
         &                         *gravitationalConstantGalacticus    &
         &                         *massEnclosed                       &
         &                         *(                                  &
         &                           +1.0d0/tidallyHeatedRadiusFinal   &
         &                           -1.0d0/radiusInitial              &
         &                          )
    return
  end function tidallyHeatedRadiusInitialRoot
    
  double precision function tidallyHeatedPotential(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    !% \ttfamily radius} (given in units of Mpc).
    use FGSL
    use Numerical_Integration
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Physical
    use Galactic_Structure_Options
    implicit none
    class           (darkMatterProfileTidallyHeated), intent(inout)           :: self
    type            (treeNode                      ), intent(inout), pointer  :: node
    double precision                                , intent(in   )           :: radius
    integer                                         , intent(  out), optional :: status    
    class           (nodeComponentSatellite        )               , pointer  :: satellite
    class           (darkMatterHaloScaleClass      )               , pointer  :: darkMatterHaloScale_
    double precision                                , parameter               :: radiusMaximumFactor =1.0d2
    type            (fgsl_function                 )                          :: integrandFunction
    type            (fgsl_integration_workspace    )                          :: integrationWorkspace
    double precision                                                          :: radiusMaximum

    if (present(status)) status=structureErrorCodeSuccess   
    satellite => node%satellite()
    if (satellite%tidalHeatingNormalized() <= 0.0d0) then
       ! No tidal heating has occurred - fall through to the unheated profile.
       tidallyHeatedPotential=self%unheatedProfile%potential(node,radius)
    else
       ! The only option here is to do a numerical integration.
       call tidallyHeatedSetGlobalSelf(self)
       tidallyHeatedNode    => node
       darkMatterHaloScale_ => darkMatterHaloScale()
       radiusMaximum        =  +radiusMaximumFactor                                              &
            &                  *self%radiusInitial(node,darkMatterHaloScale_%virialRadius(node))
       tidallyHeatedPotential=Integrate(                                        &
            &                           radius                                , &
            &                           radiusMaximum                         , &
            &                           tidallyHeatedPotentialIntegrand       , &
            &                           integrandFunction                     , &
            &                           integrationWorkspace                  , &
            &                           toleranceAbsolute              =0.0d+0, &
            &                           toleranceRelative              =1.0d-3  &
            &                          )
       call Integrate_Done(integrandFunction,integrationWorkspace)
    end if
    return
  end function tidallyHeatedPotential

  double precision function tidallyHeatedPotentialIntegrand(radius)
    !% Integrand for gravitational potential in a tidally-heated dark matter profile.
    use Numerical_Constants_Physical
    implicit none
    double precision, intent(in   ) :: radius
    
    tidallyHeatedPotentialIntegrand=-gravitationalConstantGalacticus                             &
         &                          *tidallyHeatedSelf%enclosedMass(tidallyHeatedNode,radius)    &
         &                          /                                                 radius **2
    return
  end function tidallyHeatedPotentialIntegrand
    
  double precision function tidallyHeatedCircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Numerical_Constants_Physical
    implicit none
    class           (darkMatterProfileTidallyHeated), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: radius

    if (radius > 0.0d0) then
       tidallyHeatedCircularVelocity=sqrt(                                 &
            &                             +gravitationalConstantGalacticus &
            &                             *self%enclosedMass(node,radius)  &
            &                             /                       radius   &
            &                            )
    else
       tidallyHeatedCircularVelocity=0.0d0
    end if
    return
  end function tidallyHeatedCircularVelocity

  double precision function tidallyHeatedCircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class(darkMatterProfileTidallyHeated), intent(inout)          :: self
    type (treeNode                      ), intent(inout), pointer :: node
    class(nodeComponentSatellite        )               , pointer :: satellite

    satellite => node%satellite()
    if (satellite%tidalHeatingNormalized() <= 0.0d0) then
       ! No tidal heating has occurred - fall through to the unheated profile.
       tidallyHeatedCircularVelocityMaximum=self%unheatedProfile%circularVelocityMaximum(node)
    else
       ! Check if unimplemented features are fatal.
       if (self%unimplementedFatal) then
          tidallyHeatedCircularVelocityMaximum=0.0d0
          call Galacticus_Error_Report('tidallyHeatedCircularVelocityMaximum','circular velocity maximum in tidally heated dark matter profiles is not supported')
       else
          ! Issue a warning and then fall through to the unheated profile.
          if (.not.tidallyHeatedCircularVelocityMaximumWarned) then
             !$omp critical(tidallyHeatedDarkMatterProfileCircularVelocityMaximumWarn)
             if (.not.tidallyHeatedCircularVelocityMaximumWarned) then
                call Galacticus_Display_Message('WARNING: circular velocity maximum in tidally heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
                tidallyHeatedCircularVelocityMaximumWarned=.true.
             end if
             !$omp end critical(tidallyHeatedDarkMatterProfileCircularVelocityMaximumWarn)
          end if
          tidallyHeatedCircularVelocityMaximum=self%unheatedProfile%circularVelocityMaximum(node)
       end if
    end if
    return
  end function tidallyHeatedCircularVelocityMaximum

  double precision function tidallyHeatedRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc).
    use Root_Finder
    implicit none
    class           (darkMatterProfileTidallyHeated), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: specificAngularMomentum
    class           (nodeComponentSatellite        )               , pointer :: satellite
    double precision                                , parameter              :: toleranceAbsolute     =0.0d0, toleranceRelative=1.0d-6
    type            (rootFinder                    ), save                   :: finder
    !$omp threadprivate(finder)

    ! Handle the trivial case.
    if (specificAngularMomentum <= 0.0d0) then
       tidallyHeatedRadiusFromSpecificAngularMomentum=0.0d0
    else
       ! Compute radius in unheated profile.
       tidallyHeatedRadiusFromSpecificAngularMomentum=self%unheatedProfile%radiusFromSpecificAngularMomentum(node,specificAngularMomentum)
       ! If profile is heated, perform the full calculation.
       satellite => node%satellite()
       if (satellite%tidalHeatingNormalized() > 0.0d0) then
          ! Set global pointers.
          call tidallyHeatedSetGlobalSelf(self)
          tidallyHeatedNode                   => node
          tidallyHeatedSpecificAngularMomentum=  specificAngularMomentum
          ! Initialize the root finder.
          if (.not.finder%isInitialized()) then
             call finder%rootFunction(tidallyHeatedRadiusFromSpecificAngularMomentumRoot)
             call finder%tolerance   (toleranceAbsolute,toleranceRelative)
             call finder%rangeExpand (                                                             &
                  &                   rangeExpandUpward            =2.0d0                        , &
                  &                   rangeExpandDownward          =0.5d0                        , &
                  &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
                  &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
                  &                   rangeExpandType              =rangeExpandMultiplicative      &
                  &                  )
          end if
          tidallyHeatedRadiusFromSpecificAngularMomentum=finder%find(rootGuess=tidallyHeatedRadiusFromSpecificAngularMomentum)
       end if
    end if
    return
  end function tidallyHeatedRadiusFromSpecificAngularMomentum

  double precision function tidallyHeatedRadiusFromSpecificAngularMomentumRoot(radius)
    !% Root function used in finding radii corresponding to given specific angular momenta in
    !% tidally heated dark matter profiles.
    implicit none
    double precision, intent(in   ) :: radius

    tidallyHeatedRadiusFromSpecificAngularMomentumRoot=tidallyHeatedSelf%circularVelocity(tidallyHeatedNode,radius)*radius-tidallyHeatedSpecificAngularMomentum
    return
  end function tidallyHeatedRadiusFromSpecificAngularMomentumRoot
  
  double precision function tidallyHeatedRotationNormalization(self,node)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTidallyHeated), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    class           (nodeComponentSatellite        )               , pointer :: satellite

    satellite => node%satellite()
    if (satellite%tidalHeatingNormalized() <= 0.0d0) then
       ! No tidal heating has occurred - fall through to the unheated profile.
       tidallyHeatedRotationNormalization=self%unheatedProfile%rotationNormalization(node)
    else
       ! Check if unimplemented features are fatal.
       if (self%unimplementedFatal) then
          tidallyHeatedRotationNormalization=0.0d0
          call Galacticus_Error_Report('tidallyHeatedRotationNormalization','rotation normalization in tidally heated dark matter profiles is not supported')
       else
          ! Issue a warning and then fall through to the unheated profile.
          if (.not.tidallyHeatedRotationNormalizationWarned) then
             !$omp critical(tidallyHeatedDarkMatterProfileRotationNormalizationWarn)
             if (.not.tidallyHeatedRotationNormalizationWarned) then
                call Galacticus_Display_Message('WARNING: rotation normalization in tidally heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
                tidallyHeatedRotationNormalizationWarned=.true.
             end if
             !$omp end critical(tidallyHeatedDarkMatterProfileRotationNormalizationWarn)
          end if
          tidallyHeatedRotationNormalization=self%unheatedProfile%rotationNormalization(node)
       end if
    end if
    return
  end function tidallyHeatedRotationNormalization

  double precision function tidallyHeatedEnergy(self,node)
    !% Return the energy of a tidally heated halo density profile.
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTidallyHeated), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    class           (nodeComponentSatellite        )               , pointer :: satellite

    satellite => node%satellite()
    if (satellite%tidalHeatingNormalized() <= 0.0d0) then
       ! No tidal heating has occurred - fall through to the unheated profile.
       tidallyHeatedEnergy=self%unheatedProfile%energy(node)
    else
       ! Check if unimplemented features are fatal.
       if (self%unimplementedFatal) then
          tidallyHeatedEnergy=0.0d0
          call Galacticus_Error_Report('tidallyHeatedEnergy','energy in tidally heated dark matter profiles is not supported')
       else
          ! Issue a warning and then fall through to the unheated profile.
          if (.not.tidallyHeatedEnergyWarned) then
             !$omp critical(tidallyHeatedDarkMatterProfileEnergyWarn)
             if (.not.tidallyHeatedEnergyWarned) then
                call Galacticus_Display_Message('WARNING: energy in tidally heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
                tidallyHeatedEnergyWarned=.true.
             end if
             !$omp end critical(tidallyHeatedDarkMatterProfileEnergyWarn)
          end if
          tidallyHeatedEnergy=self%unheatedProfile%energy(node)
       end if
    end if
    return
  end function tidallyHeatedEnergy

  double precision function tidallyHeatedEnergyGrowthRate(self,node)
    !% Return the rate of change of the energy of a tidally heated halo density profile.
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTidallyHeated), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    class           (nodeComponentSatellite        )               , pointer :: satellite

    satellite => node%satellite()
    if (satellite%tidalHeatingNormalized() <= 0.0d0) then
       ! No tidal heating has occurred - fall through to the unheated profile.
       tidallyHeatedEnergyGrowthRate=self%unheatedProfile%energyGrowthRate(node)
    else
       ! Check if unimplemented features are fatal.
       if (self%unimplementedFatal) then
          tidallyHeatedEnergyGrowthRate=0.0d0
          call Galacticus_Error_Report('tidallyHeatedEnergyGrowthRate','energy growth rate in tidally heated dark matter profiles is not supported')
       else
          ! Issue a warning and then fall through to the unheated profile.
          if (.not.tidallyHeatedEnergyGrowthRateWarned) then
             !$omp critical(tidallyHeatedDarkMatterProfileEnergyGrowthRateWarn)
             if (.not.tidallyHeatedEnergyGrowthRateWarned) then
                call Galacticus_Display_Message('WARNING: energyGrowthRate in tidally heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
                tidallyHeatedEnergyGrowthRateWarned=.true.
             end if
             !$omp end critical(tidallyHeatedDarkMatterProfileEnergyGrowthRateWarn)
          end if
          tidallyHeatedEnergyGrowthRate=self%unheatedProfile%energyGrowthRate(node)
       end if
    end if
    return
  end function tidallyHeatedEnergyGrowthRate
  
  double precision function tidallyHeatedKSpace(self,node,waveNumber)
    !% Returns the Fourier transform of the tidally heated density profile at the specified {\normalfont \ttfamily waveNumber}
    !% (given in Mpc$^{-1}$), using the expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTidallyHeated), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: waveNumber
    class           (nodeComponentSatellite        )               , pointer :: satellite

    satellite => node%satellite()
    if (satellite%tidalHeatingNormalized() <= 0.0d0) then
       ! No tidal heating has occurred - fall through to the unheated profile.
       tidallyHeatedKSpace=self%unheatedProfile%kSpace(node,waveNumber)
    else
       ! Check if unimplemented features are fatal.
       if (self%unimplementedFatal) then
          tidallyHeatedKSpace=0.0d0
          call Galacticus_Error_Report('tidallyHeatedKSpace','Fourier profile in tidally heated dark matter profiles is not supported')
       else
          ! Issue a warning and then fall through to the unheated profile.
          if (.not.tidallyHeatedKSpaceWarned) then
             !$omp critical(tidallyHeatedDarkMatterProfileKSpaceWarn)
             if (.not.tidallyHeatedKSpaceWarned) then
                call Galacticus_Display_Message('WARNING: Fourier profile in tidally heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
                tidallyHeatedKSpaceWarned=.true.
             end if
             !$omp end critical(tidallyHeatedDarkMatterProfileKSpaceWarn)
          end if
          tidallyHeatedKSpace=self%unheatedProfile%kSpace(node,waveNumber)
       end if
    end if
    return
  end function tidallyHeatedKSpace

  double precision function tidallyHeatedFreefallRadius(self,node,time)
    !% Returns the freefall radius in the tidally heated density profile at the specified {\normalfont \ttfamily time} (given in
    !% Gyr).
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTidallyHeated), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: time
    class           (nodeComponentSatellite        )               , pointer :: satellite

    satellite => node%satellite()
    if (satellite%tidalHeatingNormalized() <= 0.0d0) then
       ! No tidal heating has occurred - fall through to the unheated profile.
       tidallyHeatedFreefallRadius=self%unheatedProfile%freefallRadius(node,time)
    else
       ! Check if unimplemented features are fatal.
       if (self%unimplementedFatal) then
          tidallyHeatedFreefallRadius=0.0d0
          call Galacticus_Error_Report('tidallyHeatedFreefallRadius','freefall radius in tidally heated dark matter profiles is not supported')
       else
          ! Issue a warning and then fall through to the unheated profile.
          if (.not.tidallyHeatedFreefallRadiusWarned) then
             !$omp critical(tidallyHeatedDarkMatterProfileFreefallRadiusWarn)
             if (.not.tidallyHeatedFreefallRadiusWarned) then
                call Galacticus_Display_Message('WARNING: freefallRadius in tidally heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
                tidallyHeatedFreefallRadiusWarned=.true.
             end if
             !$omp end critical(tidallyHeatedDarkMatterProfileFreefallRadiusWarn)
          end if
          tidallyHeatedFreefallRadius=self%unheatedProfile%freefallRadius(node,time)
       end if
    end if
    return
  end function tidallyHeatedFreefallRadius

  double precision function tidallyHeatedFreefallRadiusIncreaseRate(self,node,time)
    !% Returns the rate of increase of the freefall radius in the tidally heated density profile at the specified {\normalfont
    !% \ttfamily time} (given in Gyr).
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTidallyHeated), intent(inout)          :: self
    type            (treeNode                      ), intent(inout), pointer :: node
    double precision                                , intent(in   )          :: time
    class           (nodeComponentSatellite        )               , pointer :: satellite

    satellite => node%satellite()
    if (satellite%tidalHeatingNormalized() <= 0.0d0) then
       ! No tidal heating has occurred - fall through to the unheated profile.
       tidallyHeatedFreefallRadiusIncreaseRate=self%unheatedProfile%freefallRadiusIncreaseRate(node,time)
    else
       ! Check if unimplemented features are fatal.
       if (self%unimplementedFatal) then
          tidallyHeatedFreefallRadiusIncreaseRate=0.0d0
          call Galacticus_Error_Report('tidallyHeatedFreefallRadiusIncreaseRate','freefall radius increase rate in tidally heated dark matter profiles is not supported')
       else
          ! Issue a warning and then fall through to the unheated profile.
          if (.not.tidallyHeatedFreefallRadiusIncreaseRateWarned) then
             !$omp critical(tidallyHeatedDarkMatterProfileFreefallRadiusIncreaseRateWarn)
             if (.not.tidallyHeatedFreefallRadiusIncreaseRateWarned) then
                call Galacticus_Display_Message('WARNING: freefall radius increase rate in tidally heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
                tidallyHeatedFreefallRadiusIncreaseRateWarned=.true.
             end if
             !$omp end critical(tidallyHeatedDarkMatterProfileFreefallRadiusIncreaseRateWarn)
          end if
          tidallyHeatedFreefallRadiusIncreaseRate=self%unheatedProfile%freefallRadiusIncreaseRate(node,time)
       end if
    end if
    return
  end function tidallyHeatedFreefallRadiusIncreaseRate

  subroutine tidallyHeatedStateStore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    implicit none
    class  (darkMatterProfileTidallyHeated), intent(inout) :: self
    integer                                , intent(in   ) :: stateFile
    type   (fgsl_file                     ), intent(in   ) :: fgslStateFile

    call self%unheatedProfile%stateStore(stateFile,fgslStateFile)
    return
  end subroutine tidallyHeatedStateStore

  subroutine tidallyHeatedStateRestore(self,stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    implicit none
    class  (darkMatterProfileTidallyHeated), intent(inout) :: self
    integer                                , intent(in   ) :: stateFile
    type   (fgsl_file                     ), intent(in   ) :: fgslStateFile

    call self%unheatedProfile%stateRestore(stateFile,fgslStateFile)
    return
  end subroutine tidallyHeatedStateRestore
