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

!+    Contributions to this file made by: Xiaolong Du, Andrew Benson.

  !% An implementation of exponentially truncated dark matter halo profiles \cite{kazantzidis_2006}.

  use Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass, darkMatterHaloScale

  !# <darkMatterProfile name="darkMatterProfileTruncatedExponential">
  !#  <description>exponentially truncated dark matter halo profiles \cite{kazantzidis_2006}.</description>
  !# </darkMatterProfile>
  type, extends(darkMatterProfileClass) :: darkMatterProfileTruncatedExponential
     !% A dark matter halo profile class implementing exponentially truncated dark matter halos.
     private
     class           (darkMatterProfileClass  ), pointer :: darkMatterProfile_   => null()
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: radiusFractionalDecay      , alpha, &
          &                                                 beta                       , gamma
     logical                                             :: unimplementedIsFatal
     ! Record of unique ID of node which we last computed results for.
     integer         (kind=kind_int8          )          :: lastUniqueID
     ! Stored values of computed quantities.
     double precision                                    :: enclosingMassRadiusPrevious, kappaPrevious
   contains
     final                                             truncatedExponentialDestructor
     procedure :: calculationReset                  => truncatedExponentialCalculationReset
     procedure :: density                           => truncatedExponentialDensity
     procedure :: densityLogSlope                   => truncatedExponentialDensityLogSlope
     procedure :: radiusEnclosingDensity            => truncatedExponentialRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => truncatedExponentialRadiusEnclosingMass
     procedure :: radialMoment                      => truncatedExponentialRadialMoment
     procedure :: enclosedMass                      => truncatedExponentialEnclosedMass
     procedure :: potential                         => truncatedExponentialPotential
     procedure :: circularVelocity                  => truncatedExponentialCircularVelocity
     procedure :: circularVelocityMaximum           => truncatedExponentialCircularVelocityMaximum
     procedure :: radiusFromSpecificAngularMomentum => truncatedExponentialRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => truncatedExponentialRotationNormalization
     procedure :: energy                            => truncatedExponentialEnergy
     procedure :: energyGrowthRate                  => truncatedExponentialEnergyGrowthRate
     procedure :: kSpace                            => truncatedExponentialKSpace
     procedure :: freefallRadius                    => truncatedExponentialFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => truncatedExponentialFreefallRadiusIncreaseRate     
  end type darkMatterProfileTruncatedExponential

  interface darkMatterProfileTruncatedExponential
     !% Constructors for the {\normalfont \ttfamily exponentially truncated} dark matter halo profile class.
     module procedure truncatedExponentialConstructorParameters
     module procedure truncatedExponentialConstructorInternal
  end interface darkMatterProfileTruncatedExponential

contains

  function truncatedExponentialConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily exponentially truncated} dark matter halo profile class which takes a parameter set as input.
    use Galacticus_Error
    use Input_Parameters
    implicit none
    type            (darkMatterProfileTruncatedExponential)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (darkMatterProfileClass               ), pointer       :: darkMatterProfile_
    class           (darkMatterHaloScaleClass             ), pointer       :: darkMatterHaloScale_
    logical                                                                :: unimplementedIsFatal
    double precision                                                       :: radiusFractionalDecay, alpha, &
         &                                                                    beta                 , gamma

    !# <inputParameter>
    !#   <name>unimplementedIsFatal</name>
    !#   <defaultValue>.true.</defaultValue>
    !#   <source>parameters</source>
    !#   <description>If {\normalfont \ttfamily true}, unimplemented features of the exponentially truncated dark matter profile cause fatal errors. Otherwise, the solution for an untruncated profile is returned.</description>
    !#   <type>boolean</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>radiusFractionalDecay</name>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>The truncation scale (in units of the virial radius).</description>
    !#   <type>float</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>alpha</name>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Parameter $\alpha$ in the \cite{kazantzidis_2006} truncated profile.</description>
    !#   <type>float</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>beta</name>
    !#   <defaultValue>3.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Parameter $\beta$ in the \cite{kazantzidis_2006} truncated profile.</description>
    !#   <type>float</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>gamma</name>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Parameter $\gamma$ in the \cite{kazantzidis_2006} truncated profile.</description>
    !#   <type>float</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterProfile"   name="darkMatterProfile_"   source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=darkMatterProfileTruncatedExponential(radiusFractionalDecay,alpha,beta,gamma,unimplementedIsFatal,darkMatterProfile_,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterProfile_"  />
    !# <objectDestructor name="darkMatterHaloScale_"/>
    return
  end function truncatedExponentialConstructorParameters

  function truncatedExponentialConstructorInternal(radiusFractionalDecay,alpha,beta,gamma,unimplementedIsFatal,darkMatterProfile_,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily exponentially truncated} dark matter profile class.
    implicit none
    type            (darkMatterProfileTruncatedExponential)                        :: self
    class           (darkMatterProfileClass               ), intent(in   ), target :: darkMatterProfile_
    class           (darkMatterHaloScaleClass             ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                       , intent(in   )         :: radiusFractionalDecay, alpha, &
         &                                                                            beta                 , gamma
    logical                                                , intent(in   )         :: unimplementedIsFatal
    !# <constructorAssign variables="radiusFractionalDecay,alpha,beta,gamma,unimplementedIsFatal,*darkMatterProfile_,*darkMatterHaloScale_"/>

    self%lastUniqueID=-1_kind_int8
    return
  end function truncatedExponentialConstructorInternal

  subroutine truncatedExponentialDestructor(self)
    !% Destructor for the {\normalfont \ttfamily exponentially truncated} dark matter halo profile class.
    implicit none
    type(darkMatterProfileTruncatedExponential), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterProfile_"   />
    !# <objectDestructor name="self%darkMatterHaloScale_" />
    return
  end subroutine truncatedExponentialDestructor

  subroutine truncatedExponentialCalculationReset(self,node)
    !% Reset the dark matter profile calculation.
    implicit none
    class(darkMatterProfileTruncatedExponential), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node

    self%lastUniqueID               = node%uniqueID()
    self%kappaPrevious              =-huge(0.0d0)
    self%enclosingMassRadiusPrevious=-1.0d0
    call self%darkMatterHaloScale_%calculationReset(node)
    return
  end subroutine truncatedExponentialCalculationReset 

  double precision function truncatedExponentialDensity(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class           (darkMatterProfileTruncatedExponential), intent(inout) :: self
    type            (treeNode                             ), intent(inout) :: node
    double precision                                       , intent(in   ) :: radius
    class           (nodeComponentDarkMatterProfile       ), pointer       :: darkMatterProfile
    double precision                                                       :: radiusVirial     , scaleRadius, &
         &                                                                    concentration    , radiusDecay

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Get the virial radius.
    radiusVirial=self%darkMatterHaloScale_%virialRadius(node)
    if (radius <= radiusVirial) then
       truncatedExponentialDensity=+self%darkMatterProfile_%density(node,radius)
    else
       ! Compute kappa if required.
       if (self%kappaPrevious == -huge(0.0d0)) then
          darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
          scaleRadius       =  darkMatterProfile%scale()
          concentration     =  radiusVirial/scaleRadius
          self%kappaPrevious= -(self%gamma+self%beta*concentration**self%alpha)/(1.0d0+concentration**self%alpha) &
               &              +1.0d0/self%radiusFractionalDecay
       end if
       ! Compute decay scale.
       radiusDecay                =+self%radiusFractionalDecay*radiusVirial
       truncatedExponentialDensity=+self%darkMatterProfile_%density(node,radiusVirial) &
            &                      *(radius/radiusVirial)**self%kappaPrevious          &
            &                      *exp(-(radius-radiusVirial)/radiusDecay)
    end if
    return
  end function truncatedExponentialDensity

  double precision function truncatedExponentialDensityLogSlope(self,node,radius)
    !% Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncatedExponential), intent(inout) :: self
    type            (treeNode                             ), intent(inout) :: node
    double precision                                       , intent(in   ) :: radius

    if (self%unimplementedIsFatal) then
       truncatedExponentialDensityLogSlope=0.0d0
       call Galacticus_Error_Report('density logarithmic slope in exponentially truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedExponentialDensityLogSlope=self%darkMatterProfile_%densityLogSlope(node,radius)
    end if
    return
  end function truncatedExponentialDensityLogSlope

  double precision function truncatedExponentialRadiusEnclosingDensity(self,node,density)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncatedExponential), intent(inout), target :: self
    type            (treeNode                             ), intent(inout), target :: node
    double precision                                       , intent(in   )         :: density
    
    if (self%unimplementedIsFatal) then
       truncatedExponentialRadiusEnclosingDensity=0.0d0
       call Galacticus_Error_Report('radius enclosing density in exponentially truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedExponentialRadiusEnclosingDensity=self%darkMatterProfile_%radiusEnclosingDensity(node,density)
    end if
    return
  end function truncatedExponentialRadiusEnclosingDensity

  double precision function truncatedExponentialRadiusEnclosingMass(self,node,mass)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    use Galacticus_Nodes, only : nodeComponentBasic
    use Root_Finder
    implicit none
    class           (darkMatterProfileTruncatedExponential), intent(inout), target :: self
    type            (treeNode                             ), intent(inout), target :: node
    double precision                                       , intent(in   )         :: mass
    class           (nodeComponentBasic                   ), pointer               :: basic
    type            (rootFinder                           ), save                  :: finder
    !$omp threadprivate(finder)

    if (mass <= 0.0d0) then
       truncatedExponentialRadiusEnclosingMass=0.0d0
       return
    end if
    ! Get basic component.
    basic => node%basic()
    ! If the given mass is smaller than the virial mass, compute the radius from untruncated profile.
    if (mass <= basic%mass()) then
       truncatedExponentialRadiusEnclosingMass=self%darkMatterProfile_%radiusEnclosingMass(node,mass)
    else
       ! Initialize the root finder.
       if (.not.finder%isInitialized()) then
          call finder%rangeExpand (                                                                 &
               &                       rangeExpandDownward          =0.5d0                        , &
               &                       rangeExpandUpward            =2.0d0                        , &
               &                       rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                       rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                       rangeExpandType              =rangeExpandMultiplicative      &
               &                  )
          call finder%rootFunction(truncatedExponentialEnclosedMassRoot                         )
          call finder%tolerance   (toleranceAbsolute=0.0d0             ,toleranceRelative=1.0d-6)
       end if
       ! Check if node differs from previous one for which we performed calculations.
       if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
       if (self%enclosingMassRadiusPrevious < 0.0d0) then
          self%enclosingMassRadiusPrevious=self%darkMatterHaloScale_%virialRadius(node)
       end if
       self%enclosingMassRadiusPrevious=finder%find(rootGuess=self%enclosingMassRadiusPrevious)
       truncatedExponentialRadiusEnclosingMass=self%enclosingMassRadiusPrevious
    end if
    return
    contains
      double precision function truncatedExponentialEnclosedMassRoot(radius)
        !% Root function used in solving for the radius that encloses a given mass.
        implicit none
        double precision, intent(in) :: radius

        truncatedExponentialEnclosedMassRoot=self%enclosedMass(node,radius)-mass
        return
      end function truncatedExponentialEnclosedMassRoot
  end function truncatedExponentialRadiusEnclosingMass

  double precision function truncatedExponentialRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncatedExponential), intent(inout)           :: self
    type            (treeNode                             ), intent(inout)           :: node
    double precision                                       , intent(in   )           :: moment
    double precision                                       , intent(in   ), optional :: radiusMinimum, radiusMaximum

    if (self%unimplementedIsFatal) then
       truncatedExponentialRadialMoment=0.0d0
       call Galacticus_Error_Report('radial moment in exponentially truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedExponentialRadialMoment=self%darkMatterProfile_%radialMoment(node,moment,radiusMinimum,radiusMaximum)
    end if
    return 
  end function truncatedExponentialRadialMoment

  double precision function truncatedExponentialEnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    use FGSL                 , only : fgsl_function, fgsl_integration_workspace
    use Numerical_Integration
    implicit none
    class           (darkMatterProfileTruncatedExponential), intent(inout) :: self
    type            (treeNode                             ), intent(inout) :: node
    double precision                                       , intent(in   ) :: radius
    double precision                                                       :: radiusVirial
    type            (fgsl_function                        )                :: integrandFunction
    type            (fgsl_integration_workspace           )                :: integrationWorkspace

    radiusVirial=self%darkMatterHaloScale_%virialRadius(node)
    if (radius <= radiusVirial) then
       truncatedExponentialEnclosedMass=+self%darkMatterProfile_%enclosedMass(node,radius)
    else
       truncatedExponentialEnclosedMass=+Integrate(                                           &
            &                                             +radiusVirial                     , &
            &                                             +radius                           , &
            &                                              truncatedExponentialMassIntegrand, &
            &                                              integrandFunction                , &
            &                                              integrationWorkspace             , &
            &                                      toleranceAbsolute=+0.0d0                 , &
            &                                      toleranceRelative=+1.0d-9                  &
            &                                     )                                           &
            &                           +self%darkMatterProfile_%enclosedMass(node,radiusVirial)
       call Integrate_Done(integrandFunction,integrationWorkspace)
    end if
    return

  contains

    double precision function truncatedExponentialMassIntegrand(radius)
      !% Integrand for mass in exponentially truncated dark matter profiles.
      use Numerical_Constants_Math
      implicit none
      double precision, intent(in   ) :: radius

      truncatedExponentialMassIntegrand=4.0d0*Pi*radius**2*self%density(node,radius)
      return
    end function truncatedExponentialMassIntegrand

  end function truncatedExponentialEnclosedMass
  
  double precision function truncatedExponentialPotential(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    !% \ttfamily radius} (given in units of Mpc).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncatedExponential), intent(inout)           :: self
    type            (treeNode                             ), intent(inout), pointer  :: node
    double precision                                       , intent(in   )           :: radius
    integer                                                , intent(  out), optional :: status    
 
     if (self%unimplementedIsFatal) then
       truncatedExponentialPotential=0.0d0
       call Galacticus_Error_Report('potential in exponentially truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedExponentialPotential=self%darkMatterProfile_%potential(node,radius,status)
    end if
   return
  end function truncatedExponentialPotential

  double precision function truncatedExponentialCircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncatedExponential), intent(inout) :: self
    type            (treeNode                             ), intent(inout) :: node
    double precision                                       , intent(in   ) :: radius

    if (self%unimplementedIsFatal) then
       truncatedExponentialCircularVelocity=0.0d0
       call Galacticus_Error_Report('circular velocity in exponentially truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedExponentialCircularVelocity=self%darkMatterProfile_%circularVelocity(node,radius)
    end if
    return
  end function truncatedExponentialCircularVelocity

  double precision function truncatedExponentialCircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    use Galacticus_Error
    implicit none
    class(darkMatterProfileTruncatedExponential), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node

    if (self%unimplementedIsFatal) then
       truncatedExponentialCircularVelocityMaximum=0.0d0
       call Galacticus_Error_Report('circular velocity maximum in exponentially truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedExponentialCircularVelocityMaximum=self%darkMatterProfile_%circularVelocityMaximum(node)
    end if
    return
  end function truncatedExponentialCircularVelocityMaximum

  double precision function truncatedExponentialRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncatedExponential), intent(inout)          :: self
    type            (treeNode                             ), intent(inout), pointer :: node
    double precision                                       , intent(in   )          :: specificAngularMomentum

    if (self%unimplementedIsFatal) then
       truncatedExponentialRadiusFromSpecificAngularMomentum=0.0d0
       call Galacticus_Error_Report('radius from specific angular momentum in exponentially truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedExponentialRadiusFromSpecificAngularMomentum=self%darkMatterProfile_%radiusFromSpecificAngularMomentum(node,specificAngularMomentum)
    end if
    return
  end function truncatedExponentialRadiusFromSpecificAngularMomentum

  double precision function truncatedExponentialRotationNormalization(self,node)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Galacticus_Error
    implicit none
    class(darkMatterProfileTruncatedExponential), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node

    if (self%unimplementedIsFatal) then
       truncatedExponentialRotationNormalization=0.0d0
       call Galacticus_Error_Report('rotation normalization in exponentially truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedExponentialRotationNormalization=self%darkMatterProfile_%rotationNormalization(node)
    end if
    return
  end function truncatedExponentialRotationNormalization

  double precision function truncatedExponentialEnergy(self,node)
    !% Return the energy of a truncatedExponential halo density profile.
    use Galacticus_Error
    implicit none
    class(darkMatterProfileTruncatedExponential), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node

    if (self%unimplementedIsFatal) then
       truncatedExponentialEnergy=0.0d0
       call Galacticus_Error_Report('energy in exponentially truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedExponentialEnergy=self%darkMatterProfile_%energy(node)
    end if
    return
  end function truncatedExponentialEnergy

  double precision function truncatedExponentialEnergyGrowthRate(self,node)
    !% Return the rate of change of the energy of a truncatedExponential halo density profile.
    use Galacticus_Error
    implicit none
    class(darkMatterProfileTruncatedExponential), intent(inout)         :: self
    type (treeNode                             ), intent(inout), target :: node

    if (self%unimplementedIsFatal) then
       truncatedExponentialEnergyGrowthRate=0.0d0
       call Galacticus_Error_Report('energy growth rate in exponentially truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedExponentialEnergyGrowthRate=self%darkMatterProfile_%energyGrowthRate(node)
    end if
    return
  end function truncatedExponentialEnergyGrowthRate
  
  double precision function truncatedExponentialKSpace(self,node,waveNumber)
    !% Returns the Fourier transform of the truncatedExponential density profile at the specified {\normalfont \ttfamily waveNumber}
    !% (given in Mpc$^{-1}$).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncatedExponential), intent(inout)          :: self
    type            (treeNode                             ), intent(inout), pointer :: node
    double precision                                       , intent(in   )          :: waveNumber

    if (self%unimplementedIsFatal) then
       truncatedExponentialKSpace=0.0d0
       call Galacticus_Error_Report('k-space in exponentially truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedExponentialKSpace=self%darkMatterProfile_%kSpace(node,waveNumber)
    end if
    return
  end function truncatedExponentialKSpace

  double precision function truncatedExponentialFreefallRadius(self,node,time)
    !% Returns the freefall radius in the truncatedExponential density profile at the specified {\normalfont \ttfamily time} (given in
    !% Gyr).
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncatedExponential), intent(inout) :: self
    type            (treeNode                             ), intent(inout) :: node
    double precision                                       , intent(in   ) :: time

    if (self%unimplementedIsFatal) then
       truncatedExponentialFreefallRadius=0.0d0
       call Galacticus_Error_Report('freefall radius in exponentially truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedExponentialFreefallRadius=self%darkMatterProfile_%freefallRadius(node,time)
    end if
    return
  end function truncatedExponentialFreefallRadius

  double precision function truncatedExponentialFreefallRadiusIncreaseRate(self,node,time)
    !% Returns the rate of increase of the freefall radius in the truncatedExponential density profile at the specified {\normalfont
    !% \ttfamily time} (given in Gyr).
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncatedExponential), intent(inout) :: self
    type            (treeNode                             ), intent(inout) :: node
    double precision                                       , intent(in   ) :: time

    if (self%unimplementedIsFatal) then
       truncatedExponentialFreefallRadiusIncreaseRate=0.0d0
       call Galacticus_Error_Report('freefall radius increase rate in exponentially truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedExponentialFreefallRadiusIncreaseRate=self%darkMatterProfile_%freefallRadiusIncreaseRate(node,time)
    end if
    return
  end function truncatedExponentialFreefallRadiusIncreaseRate
