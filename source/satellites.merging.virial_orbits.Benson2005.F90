!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

  !% An implementation of virial orbits using the \cite{benson_orbital_2005} orbital parameter distribution.

  use Dark_Matter_Halo_Scales
  use Cosmology_Functions

  !# <virialOrbit name="virialOrbitBenson2005">
  !#  <description>Virial orbits using the \cite{benson_orbital_2005} orbital parameter distribution.</description>
  !# </virialOrbit>
  type, extends(virialOrbitClass) :: virialOrbitBenson2005
     !% A virial orbit class using the \cite{benson_orbital_2005} orbital parameter distribution.
     private
     class(darkMatterHaloScaleClass                          ), pointer :: darkMatterHaloScale_   => null()
     class(cosmologyFunctionsClass                           ), pointer :: cosmologyFunctions_    => null()
     class(virialDensityContrastSphericalCollapseMatterLambda), pointer :: virialDensityContrast_
   contains
     final     ::                              benson2005Destructor
     procedure :: orbit                     => benson2005Orbit
     procedure :: densityContrastDefinition => benson2005DensityContrastDefinition
  end type virialOrbitBenson2005

  interface virialOrbitBenson2005
     !% Constructors for the {\normalfont \ttfamily benson2005} virial orbit class.
     module procedure benson2005ConstructorParameters
     module procedure benson2005ConstructorInternal
  end interface virialOrbitBenson2005

contains

  function benson2005ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily benson2005} virial orbits class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type (virialOrbitBenson2005   )                :: self
    type (inputParameters         ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass), pointer       :: darkMatterHaloScale_
    class(cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_

    !# <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"  source="parameters"/>
    self=virialOrbitBenson2005(darkMatterHaloScale_,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function benson2005ConstructorParameters

  function benson2005ConstructorInternal(darkMatterHaloScale_,cosmologyFunctions_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily benson2005} virial orbits class.
    implicit none
    type (virialOrbitBenson2005   )                        :: self
    class(darkMatterHaloScaleClass), intent(in   ), target :: darkMatterHaloScale_
    class(cosmologyFunctionsClass ), intent(in   ), target :: cosmologyFunctions_
    !# <constructorAssign variables="*darkMatterHaloScale_, *cosmologyFunctions_"/>

    allocate(self%virialDensityContrast_)
    select type (virialDensityContrast_ => self%virialDensityContrast_)
    type is (virialDensityContrastSphericalCollapseMatterLambda)
       virialDensityContrast_=virialDensityContrastSphericalCollapseMatterLambda(cosmologyFunctions_)
    end select
    return
  end function benson2005ConstructorInternal

  subroutine benson2005Destructor(self)
    !% Destructor for the {\normalfont \ttfamily benson2005} virial orbits class.
    implicit none
    type(virialOrbitBenson2005), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_" />
    !# <objectDestructor name="self%cosmologyFunctions_"  />
    if     (associated(self%virialDensityContrast_)) &
         &  deallocate(self%virialDensityContrast_)
    return
  end subroutine benson2005Destructor

  function benson2005Orbit(self,node,host,acceptUnboundOrbits)
    !% Return benson2005 orbital parameters for a satellite.
    use Dark_Matter_Profile_Mass_Definitions
    use Galacticus_Error
    implicit none
    type            (keplerOrbit          )                        :: benson2005Orbit
    class           (virialOrbitBenson2005), intent(inout), target :: self
    type            (treeNode             ), intent(inout)         :: host                   , node
    logical                                , intent(in   )         :: acceptUnboundOrbits
    class           (nodeComponentBasic   ), pointer               :: hostBasic              , basic
    double precision                       , parameter             :: pMax                   =1.96797d0                              , &
         &                                                            velocityMax            =3.00000d0
    double precision                       , parameter             :: a                   (9)=[                                        &
         &                                                                                     0.390052d+01,0.247973d+01,0.102373d+02, &
         &                                                                                     0.683922d+00,0.353953d+00,0.107716d+01, &
         &                                                                                     0.509837d+00,0.206204d+00,0.314641d+00  &
         &                                                                                    ]
    double precision                       , parameter             :: boundTolerance         =1.0d-4 !  Tolerence to ensure that orbits are sufficiently bound.
    double precision                                               :: b1                             , b2                        , &
         &                                                            distributionFunction           , energyInternal            , &
         &                                                            uniformRandom                  , velocityRadialInternal    , &
         &                                                            velocityHost                   , velocityTangentialInternal, &
         &                                                            massHost                       , radiusHost                , &
         &                                                            massSatellite                  , radiusHostSelf
    logical                                                        :: foundOrbit

    ! Get basic components.
    basic     => node%basic()
    hostBasic => host%basic()
    ! Find mass, radius, and velocity in the host corresponding to the Benson (2005) virial density contrast definition.
    massHost     =Dark_Matter_Profile_Mass_Definition(host,self%virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%time()),radiusHostSelf,velocityHost)
    massSatellite=Dark_Matter_Profile_Mass_Definition(node,self%virialDensityContrast_%densityContrast(    basic%mass(),    basic%time())                            )
    ! Select an orbit.
    foundOrbit=.false.
    do while(.not.foundOrbit)
       ! Reset the orbit.
       call benson2005Orbit%reset()
       ! Set basic properties of the orbit.
       call benson2005Orbit%massesSet(massSatellite,massHost      )
       call benson2005Orbit%radiusSet(              radiusHostSelf)
       ! Select potential radial and tangential velocities.
       velocityRadialInternal    =node%hostTree%randomNumberGenerator%sample()*velocityMax
       velocityTangentialInternal=node%hostTree%randomNumberGenerator%sample()*velocityMax
       ! Evaluate distribution function for these parameters.
       b1                  =+a(3)                                                &
            &               *exp(-a (4)*( velocityTangentialInternal-a (5))**2)
       b2                  =+a(6)                                                &
            &               *exp(-a (7)*( velocityTangentialInternal-a (8))**2)
       distributionFunction=+a(1)                                                &
            &               *velocityTangentialInternal                          &
            &               *exp(-a (2)*((velocityTangentialInternal-a (9))**2)) &
            &               *exp(-b1   *( velocityRadialInternal    -b2   )**2)
       if (distributionFunction > pMax) call Galacticus_Error_Report('benson2005Orbit','distribution function exceeds expected peak value')
       uniformRandom=pMax*node%hostTree%randomNumberGenerator%sample()
       if (uniformRandom <= distributionFunction) then
          foundOrbit=.true.
          ! If requested, check that the orbit is bound. We require it to have E<-boundTolerance to ensure that it is sufficiently
          ! bound that later rounding errors will not make it appear unbound.
          if (.not.acceptUnboundOrbits) then
             energyInternal=-1.0d0+0.5d0*(velocityRadialInternal**2+velocityTangentialInternal**2)*benson2005Orbit%specificReducedMass()
             foundOrbit=(energyInternal < -boundTolerance)
          end if
       end if
       if (.not.foundOrbit) cycle
       call benson2005Orbit%velocityRadialSet    (velocityRadialInternal    *velocityHost)
       call benson2005Orbit%velocityTangentialSet(velocityTangentialInternal*velocityHost)
       ! Propagate the orbit to the virial radius under the default density contrast definition.
       radiusHost=self%darkMatterHaloScale_%virialRadius(host)
       if (benson2005Orbit%radiusApocenter() >= radiusHost .and. benson2005Orbit%radiusPericenter() <= radiusHost) then
          foundOrbit=.true.
          call benson2005Orbit%propagate(radiusHost  ,infalling=.true.)
          call benson2005Orbit%massesSet(basic%mass(),hostBasic%mass())
       end if
    end do
    return
  end function benson2005Orbit

  function benson2005DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of \cite{benson_orbital_2005} virial orbits.
    implicit none
    class(virialDensityContrastClass), pointer       :: benson2005DensityContrastDefinition
    class(virialOrbitBenson2005     ), intent(inout) :: self
    
    benson2005DensityContrastDefinition => self%virialDensityContrast_
    return
  end function benson2005DensityContrastDefinition
