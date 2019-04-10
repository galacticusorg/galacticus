!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
!!    Andrew Benson <abenson@carnegiescience.edu>
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

  !% An implementation of heated dark matter halo profiles.

  !# <darkMatterProfileDMO name="darkMatterProfileDMOHeated">
  !#  <description>Heated dark matter halo profiles.</description>
  !# </darkMatterProfileDMO>

  use Kind_Numbers
  use Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass, darkMatterHaloScale

  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOHeated
     !% A dark matter halo profile class implementing heated dark matter halos.
     private
     class           (darkMatterProfileDMOClass       ), pointer :: darkMatterProfileDMO_        => null()
     class           (darkMatterHaloScaleClass     ), pointer :: darkMatterHaloScale_      => null()
     class           (darkMatterProfileHeatingClass), pointer :: darkMatterProfileHeating_ => null()
     logical                                                  :: unimplementedIsFatal
     integer         (kind=kind_int8               )          :: lastUniqueID
     double precision                                         :: radiusFinalPrevious      , radiusInitialPrevious
   contains
     !@ <objectMethods>
     !@   <object>darkMatterProfileDMOHeated</object>
     !@   <objectMethod>
     !@     <method>radiusInitial</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless *type(treeNode)\textgreater} node\arginout, \doublezero\ radiusFinal\argin</arguments>
     !@     <description>Return the initial radius corresponding to the given final radius in a heated dark matter halo density profile.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final                                             heatedDestructor
     procedure :: calculationReset                  => heatedCalculationReset
     procedure :: radiusInitial                     => heatedRadiusInitial
     procedure :: density                           => heatedDensity
     procedure :: densityLogSlope                   => heatedDensityLogSlope
     procedure :: radiusEnclosingDensity            => heatedRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => heatedRadiusEnclosingMass
     procedure :: radialMoment                      => heatedRadialMoment
     procedure :: enclosedMass                      => heatedEnclosedMass
     procedure :: potential                         => heatedPotential
     procedure :: circularVelocity                  => heatedCircularVelocity
     procedure :: circularVelocityMaximum           => heatedCircularVelocityMaximum
     procedure :: radiusFromSpecificAngularMomentum => heatedRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => heatedRotationNormalization
     procedure :: energy                            => heatedEnergy
     procedure :: energyGrowthRate                  => heatedEnergyGrowthRate
     procedure :: kSpace                            => heatedKSpace
     procedure :: freefallRadius                    => heatedFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => heatedFreefallRadiusIncreaseRate     
  end type darkMatterProfileDMOHeated

  interface darkMatterProfileDMOHeated
     !% Constructors for the {\normalfont \ttfamily heated} dark matter halo profile class.
     module procedure heatedConstructorParameters
     module procedure heatedConstructorInternal
  end interface darkMatterProfileDMOHeated

  ! Warnings issued.
  logical                                            :: heatedRadiusEnclosingDensityWarned           =.false.
  logical                                            :: heatedDensityLogSlopeWarned                  =.false.
  logical                                            :: heatedRadialMomentWarned                     =.false.
  logical                                            :: heatedCircularVelocityMaximumWarned          =.false.
  logical                                            :: heatedRotationNormalizationWarned            =.false.
  logical                                            :: heatedEnergyWarned                           =.false.
  logical                                            :: heatedEnergyGrowthRateWarned                 =.false.
  logical                                            :: heatedKSpaceWarned                           =.false.
  logical                                            :: heatedFreefallRadiusWarned                   =.false.
  logical                                            :: heatedFreefallRadiusIncreaseRateWarned       =.false.

  ! Global variables used in root solving.
  double precision                                   :: heatedSpecificAngularMomentum, heatedRadiusFinal
  type            (treeNode               ), pointer :: heatedNode
  type            (darkMatterProfileDMOHeated), pointer :: heatedSelf
  !$omp threadprivate(heatedRadiusFinal,heatedNode,heatedSelf,heatedSpecificAngularMomentum)

contains

  function heatedConstructorParameters(parameters) result(self)
    !% Default constructor for the {\normalfont \ttfamily heated} dark matter halo profile class.
    use Input_Parameters
    implicit none
    type   (darkMatterProfileDMOHeated       )                :: self    
    type   (inputParameters               ), intent(inout) :: parameters
    class  (darkMatterProfileDMOClass        ), pointer       :: darkMatterProfileDMO_
    class  (darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_
    class  (darkMatterProfileHeatingClass ), pointer       :: darkMatterProfileHeating_
    logical                                                :: unimplementedIsFatal

    !# <inputParameter>
    !#   <name>unimplementedIsFatal</name>
    !#   <defaultValue>.true.</defaultValue>
    !#   <source>parameters</source>
    !#   <description>If {\normalfont \ttfamily true}, unimplemented features of the heated dark matter profile cause fatal errors. Otherwise, a warning is issued and the solution for an unheated profile is returned.</description>
    !#   <type>boolean</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterProfileDMO"        name="darkMatterProfileDMO_"   source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"      name="darkMatterHaloScale_" source="parameters"/>
    !# <objectBuilder class="darkMatterProfileHeating" name="darkMatterProfileHeating_" source="parameters"/>
    self=darkMatterProfileDMOHeated(unimplementedIsFatal,darkMatterProfileDMO_,darkMatterHaloScale_,darkMatterProfileHeating_)
    !# <inputParametersValidate source="parameters"/>

    !# <objectDestructor name="darkMatterProfileDMO_"       />
    !# <objectDestructor name="darkMatterHaloScale_"     />
    !# <objectDestructor name="darkMatterProfileHeating_"/>
    return
  end function heatedConstructorParameters

  function heatedConstructorInternal(unimplementedIsFatal,darkMatterProfileDMO_,darkMatterHaloScale_,darkMatterProfileHeating_) result(self)
    !% Generic constructor for the {\normalfont \ttfamily heated} dark matter profile class.
    implicit none
    type   (darkMatterProfileDMOHeated       )                        :: self
    class  (darkMatterProfileDMOClass        ), intent(in   ), target :: darkMatterProfileDMO_
    class  (darkMatterHaloScaleClass      ), intent(in   ), target :: darkMatterHaloScale_
    class  (darkMatterProfileHeatingClass ), intent(in   ), target :: darkMatterProfileHeating_
    logical                                , intent(in   )         :: unimplementedIsFatal
    !# <constructorAssign variables="unimplementedIsFatal,*darkMatterProfileDMO_,*darkMatterHaloScale_,*darkMatterProfileHeating_"/>

    ! Construct the object.
    self%lastUniqueID       =-1_kind_int8
    self%radiusFinalPrevious=-huge(0.0d0)
    return
  end function heatedConstructorInternal

  subroutine heatedDestructor(self)
    !% Destructor for the {\normalfont \ttfamily heated} dark matter halo profile class.
    implicit none
    type(darkMatterProfileDMOHeated), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterProfileDMO_"           />
    !# <objectDestructor name="self%darkMatterHaloScale_"      />
    !# <objectDestructor name="self%darkMatterProfileHeating_" />
    return
  end subroutine heatedDestructor
  
  subroutine heatedCalculationReset(self,node)
    !% Reset the dark matter profile calculation.
    implicit none
    class(darkMatterProfileDMOHeated), intent(inout) :: self
    type (treeNode               ), intent(inout) :: node

    ! Reset the unheated profile.
    call self%darkMatterProfileDMO_%calculationReset(node)
    ! Reset calculations for this profile.
    self%lastUniqueID       =node%uniqueID()
    self%radiusFinalPrevious=-huge(0.0d0)
    return
  end subroutine heatedCalculationReset

  subroutine heatedSetGlobalSelf(self)
    !% Set a module-scope pointer to ``{\normalfont \ttfamily self}''.
    implicit none
    class(darkMatterProfileDMOHeated), intent(inout), target :: self

    heatedSelf => self
    return
  end subroutine heatedSetGlobalSelf
  
  double precision function heatedDensity(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Display
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout) :: self
    type            (treeNode               ), intent(inout) :: node
    double precision                         , intent(in   ) :: radius
    double precision                                         :: radiusInitial , massEnclosed, &
         &                                                      densityInitial, jacobian

    radiusInitial      =self                %radiusInitial(node,radius       )
    massEnclosed       =self%darkMatterProfileDMO_%enclosedMass (node,radiusInitial)
    densityInitial     =self%darkMatterProfileDMO_%density      (node,radiusInitial)
    jacobian           =+1.0d0                                                                                                 &
         &              /(                                                                                                     &
         &                +(                                                                                                   &
         &                  +radius                                                                                            &
         &                  /radiusInitial                                                                                     &
         &                 )                                                                                               **2 &
         &                +2.0d0                                                                                               &
         &                *radius                                                                                          **2 &
         &                /gravitationalConstantGalacticus                                                                     &
         &                /massEnclosed                                                                                        &
         &                *(                                                                                                   &
         &                  +self%darkMatterProfileHeating_%specificEnergyGradient(node,self%darkMatterProfileDMO_,radiusInitial)    &
         &                  -4.0d0                                                                                             &
         &                  *Pi                                                                                                &
         &                  *radiusInitial                                                                                 **2 &
         &                  *densityInitial                                                                                    &
         &                  *self%darkMatterProfileHeating_%specificEnergy        (node,self%darkMatterProfileDMO_,radiusInitial)    &
         &                  /massEnclosed                                                                                      &
         &                 )                                                                                                   &
         &               )
    heatedDensity=+densityInitial                                                                                       &
         &               *(                                                                                                    &
         &                 +radiusInitial                                                                                      &
         &                 /radius                                                                                             &
         &                )                                                                                                **2 &
         &               *jacobian       
    return
  end function heatedDensity

  double precision function heatedDensityLogSlope(self,node,radius)
    !% Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout) :: self
    type            (treeNode               ), intent(inout) :: node
    double precision                         , intent(in   ) :: radius

    ! Check if unimplemented features are fatal.
    if (self%unimplementedIsFatal) then
       call Galacticus_Error_Report('density logarithmic slope in heated dark matter profiles is not supported'//{introspection:location})
    else
       ! Issue a warning and then fall through to the unheated profile.
       if (.not.heatedDensityLogSlopeWarned) then
          !$omp critical(heatedDarkMatterProfileDensityLogSlopeWarn)
          if (.not.heatedDensityLogSlopeWarned) then
             call Galacticus_Display_Message('WARNING: density logarithmic slope in heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
             heatedDensityLogSlopeWarned=.true.
          end if
          !$omp end critical(heatedDarkMatterProfileDensityLogSlopeWarn)
       end if
    end if
    heatedDensityLogSlope=self%darkMatterProfileDMO_%densityLogSlope(node,radius)
    return
  end function heatedDensityLogSlope

  double precision function heatedRadiusEnclosingDensity(self,node,density)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    use Galacticus_Error
    use Galacticus_Display
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout), target :: self
    type            (treeNode               ), intent(inout), target :: node
    double precision                         , intent(in   )         :: density

    ! Check if unimplemented features are fatal.
    if (self%unimplementedIsFatal) then
       heatedRadiusEnclosingDensity=0.0d0
       call Galacticus_Error_Report('radius enclosing density in heated dark matter profiles is not supported'//{introspection:location})
    else
       ! Issue a warning and then fall through to the unheated profile.
       if (.not.heatedRadiusEnclosingDensityWarned) then
          !$omp critical(heatedDarkMatterProfileRadiusEnclosingDensityWarn)
          if (.not.heatedRadiusEnclosingDensityWarned) then
             call Galacticus_Display_Message('WARNING: radius enclosing density in heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
             heatedRadiusEnclosingDensityWarned=.true.
          end if
          !$omp end critical(heatedDarkMatterProfileRadiusEnclosingDensityWarn)
       end if
       heatedRadiusEnclosingDensity=self%darkMatterProfileDMO_%radiusEnclosingDensity(node,density)
    end if
    return
  end function heatedRadiusEnclosingDensity

  double precision function heatedRadiusEnclosingMass(self,node,mass)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    use Numerical_Constants_Physical
    use Galactic_Structure_Options
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout), target :: self
    type            (treeNode               ), intent(inout), target :: node
    double precision                         , intent(in   )         :: mass
    double precision                                                 :: radiusInitial
    double precision                                                 :: energySpecific

    radiusInitial            =self%darkMatterProfileDMO_%radiusEnclosingMass(node,mass)
    energySpecific           =self%darkMatterProfileHeating_%specificEnergy(node,self%darkMatterProfileDMO_,radiusInitial)
    heatedRadiusEnclosingMass=+1.0d0                                                      &
         &                    /                                                           &
         &                    (                                                           &
         &                     +1.0d0/radiusInitial                                       &
         &                     -2.0d0/gravitationalConstantGalacticus/mass*energySpecific &
         &                    )
    ! If the radius found is negative, which means the intial shell has expanded to infinity, return the largest radius.
    if (heatedRadiusEnclosingMass < 0.0d0) then
       heatedRadiusEnclosingMass=radiusLarge
    end if
    return
  end function heatedRadiusEnclosingMass

  double precision function heatedRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Error
    use Galacticus_Display
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout)           :: self
    type            (treeNode               ), intent(inout)           :: node
    double precision                         , intent(in   )           :: moment
    double precision                         , intent(in   ), optional :: radiusMinimum, radiusMaximum

    ! Check if unimplemented features are fatal.
    if (self%unimplementedIsFatal) then
       heatedRadialMoment=0.0d0
       call Galacticus_Error_Report('radial moment in heated dark matter profiles is not supported'//{introspection:location})
    else
       ! Issue a warning and then fall through to the unheated profile.
       if (.not.heatedRadialMomentWarned) then
          !$omp critical(heatedDarkMatterProfileRadialMomentWarn)
          if (.not.heatedRadialMomentWarned) then
             call Galacticus_Display_Message('WARNING: radial moment in heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
             heatedRadialMomentWarned=.true.
          end if
          !$omp end critical(heatedDarkMatterProfileRadialMomentWarn)
       end if
       heatedRadialMoment=self%darkMatterProfileDMO_%radialMoment(node,moment,radiusMinimum,radiusMaximum)
    end if
    return 
  end function heatedRadialMoment

  double precision function heatedEnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    use Galacticus_Display
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout) :: self
    type            (treeNode               ), intent(inout) :: node
    double precision                         , intent(in   ) :: radius

    heatedEnclosedMass=self%darkMatterProfileDMO_%enclosedMass(node,self%radiusInitial(node,radius))
    return
  end function heatedEnclosedMass

  double precision function heatedRadiusInitial(self,node,radiusFinal)
    !% Find the initial radius corresponding to the given {\normalfont \ttfamily radiusFinal} in
    !% the heated dark matter profile.
    use Root_Finder
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout), target  :: self
    type            (treeNode               ), intent(inout), target  :: node
    double precision                         , intent(in   )          :: radiusFinal
    double precision                         , parameter              :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-6
    type            (rootFinder             ), save                   :: finder
    double precision                                                  :: factorExpand
    !$omp threadprivate(finder)

    ! If profile is unheated, the initial radius equals the final radius.
    if (self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,self%darkMatterProfileDMO_)) then
       heatedRadiusInitial=radiusFinal
       return
    end if
    ! Reset calculations if necessary.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Find the initial radius in the unheated profile.
    if (radiusFinal /= self%radiusFinalPrevious) then
       heatedSelf        => self
       heatedNode        => node
       heatedRadiusFinal =  radiusFinal
       if (.not.finder%isInitialized()) then
          call finder%rootFunction(heatedRadiusInitialRoot                  )
          call finder%tolerance   (toleranceAbsolute      ,toleranceRelative)
       end if
       if (self%radiusFinalPrevious <= -huge(0.0d0) .or. radiusFinal < self%radiusInitialPrevious) then    
          ! No previous solution is available, or the requested final radius is smaller than the previous initial radius. In this
          ! case, our guess for the initial radius is the final radius, and we expand the range downward to find a solution.
          call finder%rangeExpand (                                                             &
               &                   rangeExpandUpward            =1.01d0                       , &
               &                   rangeExpandDownward          =0.50d0                       , &
               &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                   rangeExpandType              =rangeExpandMultiplicative      &
               &                  )
          self%radiusInitialPrevious=finder%find(rootGuess=radiusFinal)
       else
          ! Previous solution exists, and the requested final radius is larger than the previous initial radius. Use the previous
          ! initial radius as a guess for the solution, with range expansion in steps determined by the relative values of the
          ! current and previous final radii. If the current final radius is close to the previous final radius this should give a
          ! guess for the initial radius close to the actual solution.
          if (radiusFinal > self%radiusFinalPrevious) then
             factorExpand=     radiusFinal        /self%radiusFinalPrevious
          else
             factorExpand=self%radiusFinalPrevious/     radiusFinal
          end if
          call finder%rangeExpand (                                                             &
               &                   rangeExpandUpward            =1.0d0*factorExpand           , &
               &                   rangeExpandDownward          =1.0d0/factorExpand           , &
               &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                   rangeExpandType              =rangeExpandMultiplicative      &
               &                  )
          self%radiusInitialPrevious=finder%find(rootGuess=self%radiusInitialPrevious)
       end if
       self%radiusFinalPrevious=radiusFinal
    end if
    heatedRadiusInitial=self%radiusInitialPrevious
    return
  end function heatedRadiusInitial
  
  double precision function heatedRadiusInitialRoot(radiusInitial)
    !% Root function used in finding initial radii in heated dark matter halo profiles.
    use Numerical_Constants_Physical
    implicit none
    double precision, intent(in   ) :: radiusInitial
    double precision                :: massEnclosed
    
    massEnclosed           =+heatedSelf%darkMatterProfileDMO_       %enclosedMass  (heatedNode,                              radiusInitial)
    heatedRadiusInitialRoot=+heatedSelf%darkMatterProfileHeating_%specificEnergy(heatedNode,heatedSelf%darkMatterProfileDMO_,radiusInitial) &
         &                  +0.5d0                                                                                                    &
         &                  *gravitationalConstantGalacticus                                                                          &
         &                  *massEnclosed                                                                                             &
         &                  *(                                                                                                        &
         &                    +1.0d0/heatedRadiusFinal                                                                                &
         &                    -1.0d0/radiusInitial                                                                                    &
         &                   )
    return
  end function heatedRadiusInitialRoot

  double precision function heatedPotential(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    !% \ttfamily radius} (given in units of Mpc).
    use FGSL                        , only : fgsl_function, fgsl_integration_workspace
    use Numerical_Integration
    use Numerical_Constants_Physical
    use Galactic_Structure_Options
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout)           :: self
    type            (treeNode                  ), intent(inout), pointer  :: node
    double precision                            , intent(in   )           :: radius
    integer                                     , intent(  out), optional :: status    
    double precision                            , parameter               :: radiusMaximumFactor =1.0d2
    type            (fgsl_function             )                          :: integrandFunction
    type            (fgsl_integration_workspace)                          :: integrationWorkspace
    double precision                                                      :: radiusMaximum

    if (present(status)) status=structureErrorCodeSuccess   
    ! The only option here is to do a numerical integration.
    call heatedSetGlobalSelf(self)
    heatedNode    => node
    radiusMaximum        =  +radiusMaximumFactor                                                   &
         &                  *self%radiusInitial(node,self%darkMatterHaloScale_%virialRadius(node))
    heatedPotential=Integrate(                                        &
         &                           radius                                , &
         &                           radiusMaximum                         , &
         &                           heatedPotentialIntegrand       , &
         &                           integrandFunction                     , &
         &                           integrationWorkspace                  , &
         &                           toleranceAbsolute              =0.0d+0, &
         &                           toleranceRelative              =1.0d-3  &
         &                          )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
  end function heatedPotential

  double precision function heatedPotentialIntegrand(radius)
    !% Integrand for gravitational potential in a heated dark matter profile.
    use Numerical_Constants_Physical
    implicit none
    double precision, intent(in   ) :: radius
    
    heatedPotentialIntegrand=-gravitationalConstantGalacticus                             &
         &                          *heatedSelf%enclosedMass(heatedNode,radius)    &
         &                          /                                                 radius **2
    return
  end function heatedPotentialIntegrand
    
  double precision function heatedCircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Numerical_Constants_Physical
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout) :: self
    type            (treeNode               ), intent(inout) :: node
    double precision                         , intent(in   ) :: radius

    if (radius > 0.0d0) then
       heatedCircularVelocity=sqrt(                                 &
            &                             +gravitationalConstantGalacticus &
            &                             *self%enclosedMass(node,radius)  &
            &                             /                       radius   &
            &                            )
    else
       heatedCircularVelocity=0.0d0
    end if
    return
  end function heatedCircularVelocity

  double precision function heatedCircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class(darkMatterProfileDMOHeated), intent(inout) :: self
    type (treeNode               ), intent(inout) :: node

    ! Check if unimplemented features are fatal.
    if (self%unimplementedIsFatal) then
       heatedCircularVelocityMaximum=0.0d0
       call Galacticus_Error_Report('circular velocity maximum in heated dark matter profiles is not supported'//{introspection:location})
    else
       ! Issue a warning and then fall through to the unheated profile.
       if (.not.heatedCircularVelocityMaximumWarned) then
          !$omp critical(heatedDarkMatterProfileCircularVelocityMaximumWarn)
          if (.not.heatedCircularVelocityMaximumWarned) then
             call Galacticus_Display_Message('WARNING: circular velocity maximum in heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
             heatedCircularVelocityMaximumWarned=.true.
          end if
          !$omp end critical(heatedDarkMatterProfileCircularVelocityMaximumWarn)
       end if
       heatedCircularVelocityMaximum=self%darkMatterProfileDMO_%circularVelocityMaximum(node)
    end if
    return
  end function heatedCircularVelocityMaximum

  double precision function heatedRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc).
    use Root_Finder
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout)          :: self
    type            (treeNode               ), intent(inout), pointer :: node
    double precision                         , intent(in   )          :: specificAngularMomentum
    double precision                         , parameter              :: toleranceAbsolute     =0.0d0, toleranceRelative=1.0d-6
    type            (rootFinder             ), save                   :: finder
    !$omp threadprivate(finder)

    ! Handle the trivial case.
    if (specificAngularMomentum <= 0.0d0) then
       heatedRadiusFromSpecificAngularMomentum=0.0d0
    else
       ! Compute radius in unheated profile.
       heatedRadiusFromSpecificAngularMomentum=self%darkMatterProfileDMO_%radiusFromSpecificAngularMomentum(node,specificAngularMomentum)
       ! Find radius in the heated profile.
       ! Set global pointers.
       call heatedSetGlobalSelf(self)
       heatedNode                   => node
       heatedSpecificAngularMomentum=  specificAngularMomentum
       ! Initialize the root finder.
       if (.not.finder%isInitialized()) then
          call finder%rootFunction(heatedRadiusFromSpecificAngularMomentumRoot)
          call finder%tolerance   (toleranceAbsolute,toleranceRelative)
          call finder%rangeExpand (                                                             &
               &                   rangeExpandUpward            =2.0d0                        , &
               &                   rangeExpandDownward          =0.5d0                        , &
               &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                   rangeExpandType              =rangeExpandMultiplicative      &
               &                  )
       end if
       heatedRadiusFromSpecificAngularMomentum=finder%find(rootGuess=heatedRadiusFromSpecificAngularMomentum)      
    end if
    return
  end function heatedRadiusFromSpecificAngularMomentum

  double precision function heatedRadiusFromSpecificAngularMomentumRoot(radius)
    !% Root function used in finding radii corresponding to given specific angular momenta in
    !% heated dark matter profiles.
    implicit none
    double precision, intent(in   ) :: radius

    heatedRadiusFromSpecificAngularMomentumRoot=heatedSelf%circularVelocity(heatedNode,radius)*radius-heatedSpecificAngularMomentum
    return
  end function heatedRadiusFromSpecificAngularMomentumRoot
  
  double precision function heatedRotationNormalization(self,node)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class(darkMatterProfileDMOHeated), intent(inout) :: self
    type (treeNode               ), intent(inout) :: node

    ! Check if unimplemented features are fatal.
    if (self%unimplementedIsFatal) then
       heatedRotationNormalization=0.0d0
       call Galacticus_Error_Report('rotation normalization in heated dark matter profiles is not supported'//{introspection:location})
    else       ! Issue a warning and then fall through to the unheated profile.
       if (.not.heatedRotationNormalizationWarned) then
          !$omp critical(heatedDarkMatterProfileRotationNormalizationWarn)
          if (.not.heatedRotationNormalizationWarned) then
             call Galacticus_Display_Message('WARNING: rotation normalization in heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
             heatedRotationNormalizationWarned=.true.
          end if
          !$omp end critical(heatedDarkMatterProfileRotationNormalizationWarn)
       end if
       heatedRotationNormalization=self%darkMatterProfileDMO_%rotationNormalization(node)
    end if
    return
  end function heatedRotationNormalization

  double precision function heatedEnergy(self,node)
    !% Return the energy of a heated halo density profile.
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class(darkMatterProfileDMOHeated), intent(inout) :: self
    type (treeNode               ), intent(inout) :: node

    ! Check if unimplemented features are fatal.
    if (self%unimplementedIsFatal) then
       heatedEnergy=0.0d0
       call Galacticus_Error_Report('energy in heated dark matter profiles is not supported'//{introspection:location})
    else
       ! Issue a warning and then fall through to the unheated profile.
       if (.not.heatedEnergyWarned) then
          !$omp critical(heatedDarkMatterProfileEnergyWarn)
          if (.not.heatedEnergyWarned) then
             call Galacticus_Display_Message('WARNING: energy in heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
             heatedEnergyWarned=.true.
          end if
          !$omp end critical(heatedDarkMatterProfileEnergyWarn)
       end if
       heatedEnergy=self%darkMatterProfileDMO_%energy(node)
    end if
    return
  end function heatedEnergy

  double precision function heatedEnergyGrowthRate(self,node)
    !% Return the rate of change of the energy of a heated halo density profile.
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class(darkMatterProfileDMOHeated), intent(inout)         :: self
    type (treeNode               ), intent(inout), target :: node

    ! Check if unimplemented features are fatal.
    if (self%unimplementedIsFatal) then
       heatedEnergyGrowthRate=0.0d0
       call Galacticus_Error_Report('energy growth rate in heated dark matter profiles is not supported'//{introspection:location})
    else
       ! Issue a warning and then fall through to the unheated profile.
       if (.not.heatedEnergyGrowthRateWarned) then
          !$omp critical(heatedDarkMatterProfileEnergyGrowthRateWarn)
          if (.not.heatedEnergyGrowthRateWarned) then
             call Galacticus_Display_Message('WARNING: energyGrowthRate in heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
             heatedEnergyGrowthRateWarned=.true.
          end if
          !$omp end critical(heatedDarkMatterProfileEnergyGrowthRateWarn)
       end if
       heatedEnergyGrowthRate=self%darkMatterProfileDMO_%energyGrowthRate(node)
    end if
    return
  end function heatedEnergyGrowthRate

  double precision function heatedKSpace(self,node,waveNumber)
    !% Returns the Fourier transform of the heated density profile at the specified {\normalfont \ttfamily waveNumber}
    !% (given in Mpc$^{-1}$), using the expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; eqn.~81).
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout)          :: self
    type            (treeNode               ), intent(inout), pointer :: node
    double precision                         , intent(in   )          :: waveNumber
    
    ! Check if unimplemented features are fatal.
    if (self%unimplementedIsFatal) then
       heatedKSpace=0.0d0
       call Galacticus_Error_Report('Fourier profile in heated dark matter profiles is not supported'//{introspection:location})
    else
       ! Issue a warning and then fall through to the unheated profile.
       if (.not.heatedKSpaceWarned) then
          !$omp critical(heatedDarkMatterProfileKSpaceWarn)
          if (.not.heatedKSpaceWarned) then
             call Galacticus_Display_Message('WARNING: Fourier profile in heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
             heatedKSpaceWarned=.true.
          end if
          !$omp end critical(heatedDarkMatterProfileKSpaceWarn)
       end if
       heatedKSpace=self%darkMatterProfileDMO_%kSpace(node,waveNumber)
    end if
    return
  end function heatedKSpace

  double precision function heatedFreefallRadius(self,node,time)
    !% Returns the freefall radius in the heated density profile at the specified {\normalfont \ttfamily time} (given in
    !% Gyr).
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout) :: self
    type            (treeNode               ), intent(inout) :: node
    double precision                         , intent(in   ) :: time

    ! Check if unimplemented features are fatal.
    if (self%unimplementedIsFatal) then
       heatedFreefallRadius=0.0d0
       call Galacticus_Error_Report('freefall radius in heated dark matter profiles is not supported'//{introspection:location})
    else
       ! Issue a warning and then fall through to the unheated profile.
       if (.not.heatedFreefallRadiusWarned) then
          !$omp critical(heatedDarkMatterProfileFreefallRadiusWarn)
          if (.not.heatedFreefallRadiusWarned) then
             call Galacticus_Display_Message('WARNING: freefallRadius in heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
             heatedFreefallRadiusWarned=.true.
          end if
          !$omp end critical(heatedDarkMatterProfileFreefallRadiusWarn)
       end if
       heatedFreefallRadius=self%darkMatterProfileDMO_%freefallRadius(node,time)
    end if
    return
  end function heatedFreefallRadius

  double precision function heatedFreefallRadiusIncreaseRate(self,node,time)
    !% Returns the rate of increase of the freefall radius in the heated density profile at the specified {\normalfont
    !% \ttfamily time} (given in Gyr).
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileDMOHeated), intent(inout) :: self
    type            (treeNode               ), intent(inout) :: node
    double precision                         , intent(in   ) :: time

    ! Check if unimplemented features are fatal.
    if (self%unimplementedIsFatal) then
       heatedFreefallRadiusIncreaseRate=0.0d0
       call Galacticus_Error_Report('freefall radius increase rate in heated dark matter profiles is not supported'//{introspection:location})
    else
       ! Issue a warning and then fall through to the unheated profile.
       if (.not.heatedFreefallRadiusIncreaseRateWarned) then
          !$omp critical(heatedDarkMatterProfileFreefallRadiusIncreaseRateWarn)
          if (.not.heatedFreefallRadiusIncreaseRateWarned) then
             call Galacticus_Display_Message('WARNING: freefall radius increase rate in heated dark matter profiles is not supported - using unheated profile',verbosity=verbosityWarn)
             heatedFreefallRadiusIncreaseRateWarned=.true.
          end if
          !$omp end critical(heatedDarkMatterProfileFreefallRadiusIncreaseRateWarn)
       end if
       heatedFreefallRadiusIncreaseRate=self%darkMatterProfileDMO_%freefallRadiusIncreaseRate(node,time)
    end if
    return
  end function heatedFreefallRadiusIncreaseRate
