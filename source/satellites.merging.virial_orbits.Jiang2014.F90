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

  !% An implementation of virial orbits using the \cite{jiang_orbital_2014} orbital parameter distribution.

  use FGSL

  !# <virialOrbit name="virialOrbitJiang2014">
  !#  <description>Virial orbits using the \cite{jiang_orbital_2014} orbital parameter distribution.</description>
  !# </virialOrbit>

  type, extends(virialOrbitClass) :: virialOrbitJiang2014
     !% A virial orbit class using the \cite{jiang_orbital_2014} orbital parameter distribution.
     private
     type   (fgsl_rng) :: clonedPseudoSequenceObject, pseudoSequenceObject
     logical           :: resetSequence             , resetSequenceSnapshot
   contains
     procedure :: stateSnapshot             => jiang2014StateSnapshot
     procedure :: stateStore                => jiang2014StateStore
     procedure :: stateRestore              => jiang2014StateRestore
     procedure :: orbit                     => jiang2014Orbit
     procedure :: densityContrastDefinition => jiang2014DensityContrastDefinition
  end type virialOrbitJiang2014

  interface virialOrbitJiang2014
     !% Constructors for the {\tt jiang2014} virial orbit class.
     module procedure jiang2014Constructor
  end interface virialOrbitJiang2014

contains

  function jiang2014Constructor()
    !% Generic constructor for the {\tt jiang2014} virial orbits class.
    use Input_Parameters
    implicit none
    type(virialOrbitJiang2014), target :: jiang2014Constructor

    jiang2014Constructor%resetSequence=.true.
    return
  end function jiang2014Constructor

  function jiang2014Orbit(self,node,host,acceptUnboundOrbits)
    !% Return jiang2014 orbital parameters for a satellite.
    use Dark_Matter_Profile_Mass_Definitions
    use Dark_Matter_Halo_Scales
    use Galacticus_Error
    use Pseudo_Random
    use Root_Finder
    use Statistics_Distributions
    implicit none
    type            (keplerOrbit               )                                :: jiang2014Orbit
    class           (virialOrbitJiang2014      ), intent(inout)                 :: self
    type            (treeNode                  ), intent(inout), pointer        :: host                   , node
    logical                                     , intent(in   )                 :: acceptUnboundOrbits
    class           (darkMatterHaloScaleClass  )               , pointer        :: darkMatterHaloScale_
    class           (nodeComponentBasic        )               , pointer        :: hostBasic              , basic
    class           (virialDensityContrastClass), pointer                       :: virialDensityContrast_
    double precision                            , parameter                     :: boundTolerance         =1.0d-4 !  Tolerence to ensure that orbits are sufficiently bound.
    double precision                            , parameter    , dimension(3,3) :: B    =reshape(                                  &
         &                                                                                       [                                 &
         &                                                                                        +0.049d0, +0.548d0, +1.229d0,    &
         &                                                                                        +1.044d0, +1.535d0, +3.396d0,    &
         &                                                                                        +2.878d0, +3.946d0, +2.982d0     &
         &                                                                                       ]                               , &
         &                                                                                       [3,3]                             &
         &                                                                                      )
    double precision                            , parameter    , dimension(3,3) :: gamma=reshape(                                  &
         &                                                                                       [                                 &
         &                                                                                        +0.109d0, +0.114d0, +0.110d0,    &
         &                                                                                        +0.098d0, +0.087d0, +0.050d0,    &
         &                                                                                        +0.071d0, +0.030d0, -0.012d0     &
         &                                                                                       ]                               , &
         &                                                                                       [3,3]                             &
         &                                                                                      )
    double precision                            , parameter    , dimension(3,3) :: sigma=reshape(                                  &
         &                                                                                       [                                 &
         &                                                                                        +0.077d0, +0.094d0, +0.072d0,    &
         &                                                                                        +0.073d0, +0.083d0, +0.118d0,    &
         &                                                                                        +0.091d0, +0.139d0, +0.187d0     &
         &                                                                                       ]                               , &
         &                                                                                       [3,3]                             &
         &                                                                                      )
    double precision                            , parameter    , dimension(3,3) :: mu   =reshape(                                  &
         &                                                                                       [                                 &
         &                                                                                        +1.220d0, +1.231d0, +1.254d0,    &
         &                                                                                        +1.181d0, +1.201d0, +1.236d0,    &
         &                                                                                        +1.100d0, +1.100d0, +1.084d0     &
         &                                                                                       ]                               , &
         &                                                                                       [3,3]                             &
         &                                                                                      )
    double precision                            , parameter                     :: toleranceAbsolute=0.0d+0
    double precision                            , parameter                     :: toleranceRelative=1.0d-3
    type            (rootFinder                ), save                          :: totalFinder                    , radialFinder
    !$omp threadprivate(totalFinder,radialFinder)
    double precision                                                            :: velocityHost                   , radiusHost   , &
         &                                                                         massHost                       , massSatellite, &
         &                                                                         energyInternal                 , probabilityRadialNormalization, &
         &                                                                         velocityRadialInternal         , velocityTangentialInternal, &
         &                                                                         velocityTotalInternal          , xRadial                  , &
         &                                                                         xTotal
    logical                                                                     :: foundOrbit
    integer                                                                     :: i                              , j
    type            (distributionVoight        )                                :: voightDistribution

    ! Reset the orbit.
    call jiang2014Orbit%reset()
    ! Get required objects.
    darkMatterHaloScale_ => darkMatterHaloScale()
    ! Get basic components.
    basic                => node%basic         ()
    hostBasic            => host%basic         ()
    ! Find virial density contrast under Jiang et al. (2014) definition.
    virialDensityContrast_ => self %densityContrastDefinition()
    ! Find mass, radius, and velocity in the host and satellite corresponding to the Jiang et al. (2014) virial density contrast definition.
    massHost     =Dark_Matter_Profile_Mass_Definition(host,virialDensityContrast_%densityContrast(hostBasic%time()),radiusHost,velocityHost)
    massSatellite=Dark_Matter_Profile_Mass_Definition(node,virialDensityContrast_%densityContrast(hostBasic%time())                        )
    deallocate(virialDensityContrast_)
    ! Set basic properties of the orbit.
    call jiang2014Orbit%massesSet(massSatellite,massHost)
    call jiang2014Orbit%radiusSet(radiusHost            )
    ! Select parameters appropriate for this host-satellite pair.
    if      (massHost < 10.0d0**12.5d0) then
       i=1
    else if (massHost < 10.0d0**13.5d0) then
       i=2
    else
       i=3
    end if
    if      (massSatellite/massHost < 0.005d0) then
       j=1
    else if (massSatellite/massHost < 0.050d0) then
       j=2
    else
       j=3
    end if
    ! Configure finder objects.
    if (.not.totalFinder %isInitialized()) then
       call totalFinder %rootFunction(totalVelocityCDF )
       call totalFinder %tolerance   (toleranceAbsolute,toleranceRelative)
       call totalFinder %rangeExpand (                                                             &
            &                         rangeExpandUpward            =2.0d0                        , &
            &                         rangeExpandDownward          =0.5d0                        , &
            &                         rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                         rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                         rangeExpandType              =rangeExpandMultiplicative      &
            &                        )
    end if
    if (.not.radialFinder%isInitialized()) then
       call radialFinder%rootFunction(radialVelocityCDF)
       call radialFinder%tolerance   (toleranceAbsolute,toleranceRelative)
       call radialFinder%rangeExpand (                                                             &
            &                         rangeExpandUpward            =2.0d0                        , &
            &                         rangeExpandDownward          =0.5d0                        , &
            &                         rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                         rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                         rangeExpandType              =rangeExpandMultiplicative      &
            &                        )
    end if
    ! Select an orbit.
    foundOrbit=.false.
    do while(.not.foundOrbit)
       ! Solve for the total velocity.
       xTotal                        =Pseudo_Random_Get(self%pseudoSequenceObject,self%resetSequence)
       voightDistribution            =distributionVoight(gamma(i,j),mu(i,j),sigma(i,j),limitLower=0.0d0,limitUpper=2.0d0)
       velocityTotalInternal         =totalFinder %find(rootGuess=1.0d0)
       ! Solve for the radial velocity.
       xRadial                       =0.0d0
       probabilityRadialNormalization=1.0d0
       probabilityRadialNormalization=1.0d0/(radialVelocityCDF(velocityTotalInternal)-radialVelocityCDF(0.0d0))
       xRadial                       =Pseudo_Random_Get(self%pseudoSequenceObject,self%resetSequence)
       velocityRadialInternal        =radialFinder%find(rootGuess=sqrt(2.0d0)*velocityTotalInternal)
       ! Compute tangential velocity.
       velocityTangentialInternal=sqrt(max(0.0d0,velocityTotalInternal**2-velocityRadialInternal**2))
       ! If requested, check that the orbit is bound. We require it to have E<-boundTolerance to ensure that it is sufficiently
       ! bound that later rounding errors will not make it appear unbound.
       foundOrbit=.true.
       if (.not.acceptUnboundOrbits) then
          energyInternal=-1.0d0+0.5d0*velocityTotalInternal**2*jiang2014Orbit%specificReducedMass()
          foundOrbit=(energyInternal < -boundTolerance)
       end if
       if (.not.foundOrbit) cycle
       call jiang2014Orbit%velocityRadialSet    (velocityRadialInternal    *velocityHost)
       call jiang2014Orbit%velocityTangentialSet(velocityTangentialInternal*velocityHost)
       ! Propagate the orbit to the virial radius under the default density contrast definition.
       radiusHost=darkMatterHaloScale_%virialRadius(host)
       foundOrbit=.false.
       if (jiang2014Orbit%radiusApocenter() >= radiusHost .and. jiang2014Orbit%radiusPericenter() <= radiusHost) then
          foundOrbit=.true.
          call jiang2014Orbit%propagate(radiusHost  ,infalling=.true.)
          call jiang2014Orbit%massesSet(basic%mass(),hostBasic%mass())
       end if
    end do
    return

  contains

    double precision function totalVelocityCDF(velocityTotal)
      !% Cumulative distribution function for the total velocity.
      implicit none
      double precision, intent(in   ) :: velocityTotal

      totalVelocityCDF=voightDistribution%cumulative(velocityTotal)-xTotal
      return
    end function totalVelocityCDF

    double precision function radialVelocityCDF(velocityRadial)
      !% Cumulative distribution function for the radial velocity.
      implicit none
      double precision, intent(in   ) :: velocityRadial

      radialVelocityCDF=+probabilityRadialNormalization &
           &            *(                              &
           &              +velocityTotalInternal        &
           &              /B(i,j)                       &
           &              *(                            &
           &                +exp(                       &
           &                     +B(i,j)                &
           &                     *velocityRadial        &
           &                     /velocityTotalInternal &
           &                    )                       &
           &                -1.0d0                      &
           &               )                            &
           &              -velocityRadial               &
           &             )                              &
           &            -xRadial
      return
    end function radialVelocityCDF

  end function jiang2014Orbit

  function jiang2014DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of \cite{jiang_orbital_2014} virial orbits.
    implicit none
    class(virialDensityContrastClass), pointer       :: jiang2014DensityContrastDefinition
    class(virialOrbitJiang2014      ), intent(inout) :: self
    
    allocate(virialDensityContrastFixed :: jiang2014DensityContrastDefinition)
    select type (jiang2014DensityContrastDefinition)
    type is (virialDensityContrastFixed)
      jiang2014DensityContrastDefinition=virialDensityContrastFixed(200.0d0,virialDensityContrastFixedDensityTypeCritical)
    end select
    return
  end function jiang2014DensityContrastDefinition
  
  subroutine jiang2014StateSnapshot(self)
    !% Write the tablulation state to file.
    implicit none
    class(virialOrbitJiang2014), intent(inout) :: self

    if (.not.self%resetSequence) self%clonedPseudoSequenceObject=FGSL_Rng_Clone(self%pseudoSequenceObject)
    self%resetSequenceSnapshot=self%resetSequence
    return
  end subroutine jiang2014StateSnapshot

  subroutine jiang2014StateStore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use Pseudo_Random
    implicit none
    class  (virialOrbitJiang2014), intent(inout) :: self
    integer                       , intent(in   ) :: stateFile
    type   (fgsl_file            ), intent(in   ) :: fgslStateFile

    write (stateFile) self%resetSequenceSnapshot
    if (.not.self%resetSequenceSnapshot) call Pseudo_Random_Store(self%clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine jiang2014StateStore

  subroutine jiang2014StateRestore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use Pseudo_Random
    implicit none
    class  (virialOrbitJiang2014), intent(inout) :: self
    integer                       , intent(in   ) :: stateFile
    type   (fgsl_file            ), intent(in   ) :: fgslStateFile

    read (stateFile) self%resetSequence
    if (.not.self%resetSequence) call Pseudo_Random_Retrieve(self%pseudoSequenceObject,fgslStateFile)
   return
  end subroutine jiang2014StateRestore
