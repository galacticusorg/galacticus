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

  !% An implementation of virial orbits using the \cite{benson_orbital_2005} orbital parameter distribution.

  use FGSL

  !# <virialOrbit name="virialOrbitBenson2005">
  !#  <description>Virial orbits using the \cite{benson_orbital_2005} orbital parameter distribution.</description>
  !# </virialOrbit>

  type, extends(virialOrbitClass) :: virialOrbitBenson2005
     !% A virial orbit class using the \cite{benson_orbital_2005} orbital parameter distribution.
     private
     type   (fgsl_rng) :: clonedPseudoSequenceObject, pseudoSequenceObject
     logical           :: resetSequence             , resetSequenceSnapshot
   contains
     procedure :: stateSnapshot             => benson2005StateSnapshot
     procedure :: stateStore                => benson2005StateStore
     procedure :: stateRestore              => benson2005StateRestore
     procedure :: orbit                     => benson2005Orbit
     procedure :: densityContrastDefinition => benson2005DensityContrastDefinition
  end type virialOrbitBenson2005

  interface virialOrbitBenson2005
     !% Constructors for the {\tt benson2005} virial orbit class.
     module procedure benson2005Constructor
  end interface virialOrbitBenson2005

contains

  function benson2005Constructor()
    !% Generic constructor for the {\tt benson2005} virial orbits class.
    use Input_Parameters
    implicit none
    type(virialOrbitBenson2005), target :: benson2005Constructor

    benson2005Constructor%resetSequence=.true.
    return
  end function benson2005Constructor

  function benson2005Orbit(self,node,host,acceptUnboundOrbits)
    !% Return benson2005 orbital parameters for a satellite.
    use Dark_Matter_Profile_Mass_Definitions
    use Dark_Matter_Halo_Scales
    use Galacticus_Error
    use Pseudo_Random
    implicit none
    type            (keplerOrbit               )                         :: benson2005Orbit
    class           (virialOrbitBenson2005     ), intent(inout)          :: self
    type            (treeNode                  ), intent(inout), pointer :: host                   , node
    logical                                     , intent(in   )          :: acceptUnboundOrbits
    class           (darkMatterHaloScaleClass  )               , pointer :: darkMatterHaloScale_
    class           (nodeComponentBasic        )               , pointer :: hostBasic              , basic
    class           (virialDensityContrastClass), pointer                :: virialDensityContrast_
    double precision                            , parameter              :: pMax                   =1.96797d0                              , &
         &                                                                  velocityMax            =3.00000d0
    double precision                            , parameter              :: a                   (9)=[                                        &
         &                                                                                           0.390052d+01,0.247973d+01,0.102373d+02, &
         &                                                                                           0.683922d+00,0.353953d+00,0.107716d+01, &
         &                                                                                           0.509837d+00,0.206204d+00,0.314641d+00  &
         &                                                                                          ]
    double precision                            , parameter              :: boundTolerance         =1.0d-4 !  Tolerence to ensure that orbits are sufficiently bound.
    double precision                                                     :: b1                             , b2                        , &
         &                                                                  distributionFunction           , energyInternal            , &
         &                                                                  uniformRandom                  , velocityRadialInternal    , &
         &                                                                  velocityHost                   , velocityTangentialInternal, &
         &                                                                  massHost                       , radiusHost
    logical                                                              :: foundOrbit

    ! Reset the orbit.
    call benson2005Orbit%reset()
    ! Get required objects.
    darkMatterHaloScale_ => darkMatterHaloScale()
    ! Get basic components.
    basic                => node%basic         ()
    hostBasic            => host%basic         ()
    ! Find virial density contrast under Benson (2005) definition.
    virialDensityContrast_ => self %densityContrastDefinition()
    ! Find mass, radius, and velocity in the host corresponding to the Benson (2005) virial density contrast definition.
    massHost=Dark_Matter_Profile_Mass_Definition(host,virialDensityContrast_%densityContrast(hostBasic%time()),radiusHost,velocityHost)
    deallocate(virialDensityContrast_)
    ! Set basic properties of the orbit.
    call benson2005Orbit%massesSet(basic%mass(),hostBasic%mass())
    call benson2005Orbit%radiusSet(radiusHost                   )
    ! Select an orbit.
    foundOrbit=.false.
    do while(.not.foundOrbit)
       ! Select potential radial and tangential velocities.
       velocityRadialInternal    =Pseudo_Random_Get(self%pseudoSequenceObject,self%resetSequence)*velocityMax
       velocityTangentialInternal=Pseudo_Random_Get(self%pseudoSequenceObject,self%resetSequence)*velocityMax
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
       uniformRandom=pMax*Pseudo_Random_Get(self%pseudoSequenceObject,self%resetSequence)
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
       radiusHost=darkMatterHaloScale_%virialRadius(host)
       if (benson2005Orbit%radiusApocenter() >= radiusHost) then
          foundOrbit=.true.
          call benson2005Orbit%propagate(radiusHost,infalling=.true.)
       end if
    end do
    return
  end function benson2005Orbit

  function benson2005DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of \cite{benson_orbital_2005} virial orbits.
    implicit none
    class(virialDensityContrastClass), pointer       :: benson2005DensityContrastDefinition
    class(virialOrbitBenson2005     ), intent(inout) :: self
    
    allocate(virialDensityContrastSphericalCollapseMatterLambda :: benson2005DensityContrastDefinition)
    select type (benson2005DensityContrastDefinition)
    type is (virialDensityContrastSphericalCollapseMatterLambda)
      benson2005DensityContrastDefinition=virialDensityContrastSphericalCollapseMatterLambda()
    end select
    return
  end function benson2005DensityContrastDefinition
  
  subroutine benson2005StateSnapshot(self)
    !% Write the tablulation state to file.
    implicit none
    class(virialOrbitBenson2005), intent(inout) :: self

    if (.not.self%resetSequence) self%clonedPseudoSequenceObject=FGSL_Rng_Clone(self%pseudoSequenceObject)
    self%resetSequenceSnapshot=self%resetSequence
    return
  end subroutine benson2005StateSnapshot

  subroutine benson2005StateStore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use Pseudo_Random
    implicit none
    class  (virialOrbitBenson2005), intent(inout) :: self
    integer                       , intent(in   ) :: stateFile
    type   (fgsl_file            ), intent(in   ) :: fgslStateFile

    write (stateFile) self%resetSequenceSnapshot
    if (.not.self%resetSequenceSnapshot) call Pseudo_Random_Store(self%clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine benson2005StateStore

  subroutine benson2005StateRestore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use Pseudo_Random
    implicit none
    class  (virialOrbitBenson2005), intent(inout) :: self
    integer                       , intent(in   ) :: stateFile
    type   (fgsl_file            ), intent(in   ) :: fgslStateFile

    read (stateFile) self%resetSequence
    if (.not.self%resetSequence) call Pseudo_Random_Retrieve(self%pseudoSequenceObject,fgslStateFile)
   return
  end subroutine benson2005StateRestore
