!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program which tests AGN tori misaligned from their galaxies.
!!}

program Test_Dust_Attenuation_AGN_Torus_Misaligned
  !!{RST
  Tests the :galacticus-class:`galacticInclinationMisaligned` galactic inclination class, and its use by the
  :galacticus-class:`dustAttenuationAGNTorus` dust attenuation class.

  The checks are made against results found here independently of the class:

  * the inclinations returned for averaging reproduce moments of the distribution of the orientation of the
    component. For axes isotropic within a cone of half-angle :math:`\psi_\mathrm{max}` about that of the galaxy,
    :math:`\langle\cos^2\theta\rangle = \cos^2 i \langle\cos^2\psi\rangle + \sin^2 i (1-\langle\cos^2\psi\rangle)/2`,
    with :math:`\langle\cos^2\psi\rangle=(1+a+a^2)/3` for :math:`a=\cos\psi_\mathrm{max}`. The quadratures are exact
    for this moment, so it is reproduced to rounding error; for :math:`\psi_\mathrm{max}=90^\circ` it is
    :math:`1/3` whatever the inclination of the galaxy;
  * with no misalignment the component is seen at the inclination of the galaxy;
  * for a galaxy with its own misalignment the component is seen at :math:`|i-\psi|` when the axis of the component
    is tilted towards the observer (:math:`\phi=0`), and at :math:`i+\psi` (folded into :math:`[0,\pi/2]`) when tilted
    away (:math:`\phi=\pi`);
  * the transmission through a misaligned torus is the weighted sum of that through an aligned torus at the
    inclinations returned, and, for a torus whose axis is independent of its galaxy, is the average over
    :math:`\cos\theta` of the transmission through an aligned torus, whatever the inclination of the galaxy.
  !!}
  use :: Black_Hole_Torus_Covering_Fractions, only : blackHoleTorusCoveringFractionFixed
  use :: Display                            , only : displayVerbositySet                , verbosityLevelStandard
  use :: Dust_Attenuation_Descriptors       , only : emissionDescriptor                 , emissionSourceAccretionDisk
  use :: Dust_Attenuations                  , only : dustAttenuationAGNTorus
  use :: Dust_Extinction_Curves             , only : dustExtinctionCurvePowerLaw
  use :: Error                              , only : Error_Handler_Register
  use :: Events_Hooks                       , only : eventsHooksInitialize
  use :: Functions_Global_Utilities         , only : Functions_Global_Set
  use :: Galactic_Inclinations              , only : galacticInclinationFixed           , galacticInclinationMisaligned    , galacticInclinationNull
  use :: Galactic_Structure_Options         , only : componentTypeBlackHole
  use :: Galacticus_Nodes                   , only : nodeClassHierarchyInitialize       , nodeComponentBasic               , treeNode
  use :: Input_Parameters                   , only : inputParameters
  use :: Node_Components                    , only : Node_Components_Initialize         , Node_Components_Thread_Initialize
  use :: Numerical_Constants_Math           , only : Pi
  use :: Unit_Tests                         , only : Assert                             , Unit_Tests_Begin_Group           , Unit_Tests_End_Group   , Unit_Tests_Finish
  implicit none
  ! The V-band wavelength, which is the reference wavelength of the extinction curve.
  double precision                                     , parameter                     :: wavelengthV                    =5.5d3
  ! Half-angles of the cones of misalignment tested, in degrees.
  double precision                                     , parameter    , dimension(3)   :: anglesMisalignmentMaximum      =[0.0d0,40.0d0,90.0d0]
  ! Number of points in the fine grid used for an independent average over cos(θ).
  integer                                              , parameter                     :: countFine                      =200001
  type            (inputParameters                    )                                :: parameters
  type            (treeNode                           ), pointer                       :: node
  class           (nodeComponentBasic                 ), pointer                       :: basic
  type            (dustExtinctionCurvePowerLaw        ), pointer                       :: dustExtinctionCurve_
  type            (galacticInclinationNull            ), pointer                       :: galacticInclinationNull_
  type            (galacticInclinationFixed           ), pointer                       :: galacticInclinationFixed_
  type            (galacticInclinationMisaligned      ), pointer                       :: galacticInclination_
  type            (blackHoleTorusCoveringFractionFixed), pointer                       :: blackHoleTorusCoveringFraction_
  type            (dustAttenuationAGNTorus            ), pointer                       :: torusAligned_                                        , torusMisaligned_
  type            (emissionDescriptor                 ), dimension(2)                  :: descriptors
  double precision                                     , dimension(2)                  :: transmission                                         , transmissionExpected
  double precision                                     , allocatable  , dimension(:)   :: inclinations                                         , weights
  double precision                                     , dimension(countFine)          :: cosineInclinationFine                                , integrand
  double precision                                                                     :: inclinationHost                                      , cosineMisalignmentMaximum, &
       &                                                                                  momentCosineMisalignment
  integer                                                                              :: i                                                    , j                        , &
       &                                                                                  k                                                    , misalignmentDeviateID    , &
       &                                                                                  misalignmentAzimuthID
  character       (len=96                             )                                :: label

  call displayVerbositySet(verbosityLevelStandard)
  call Error_Handler_Register()
  parameters=inputParameters('testSuite/parameters/dustAttenuationAGNTorus.xml')
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! The meta-properties holding the misalignment of each galaxy, created here as `nodeOperatorMisalignmentRandom`
  ! would create them, so that their values can be set directly.
  !![
  <addMetaProperty component="basic" name="misalignmentDeviate" id="misalignmentDeviateID" isCreator="yes" isEvolvable="no"/>
  <addMetaProperty component="basic" name="misalignmentAzimuth" id="misalignmentAzimuthID" isCreator="yes" isEvolvable="no"/>
  !!]
  node  => treeNode      (                 )
  basic => node    %basic(autoCreate=.true.)
  ! A fine grid in cos(θ) for the independent average.
  do i=1,countFine
     cosineInclinationFine(i)=dble(i-1)/dble(countFine-1)
  end do
  ! Parcels: accretion disk continuum at V and at 2μm.
  descriptors%wavelength   =[wavelengthV                ,2.0d4                      ]
  descriptors%componentType=[componentTypeBlackHole     ,componentTypeBlackHole     ]
  descriptors%sourceType   =[emissionSourceAccretionDisk,emissionSourceAccretionDisk]
  allocate(dustExtinctionCurve_           )
  allocate(galacticInclinationNull_       )
  allocate(galacticInclinationFixed_      )
  allocate(blackHoleTorusCoveringFraction_)
  allocate(torusAligned_                  )
  !![
  <referenceConstruct object="dustExtinctionCurve_"            constructor="dustExtinctionCurvePowerLaw        (0.7d0,wavelengthV)"/>
  <referenceConstruct object="galacticInclinationNull_"        constructor="galacticInclinationNull            (                 )"/>
  <referenceConstruct object="galacticInclinationFixed_"       constructor="galacticInclinationFixed           (60.0d0           )"/>
  <referenceConstruct object="blackHoleTorusCoveringFraction_" constructor="blackHoleTorusCoveringFractionFixed(0.5d0            )"/>
  <referenceConstruct object="torusAligned_"                   constructor="dustAttenuationAGNTorus(5.0d0,30.0d0,dustExtinctionCurve_,galacticInclinationNull_,blackHoleTorusCoveringFraction_)"/>
  !!]

  ! Averages over orientation, with no misalignment for each galaxy.
  do j=1,size(anglesMisalignmentMaximum)
     write (label,'(a,f4.1,a)') "Misaligned inclination: averaged, ψ_max=",anglesMisalignmentMaximum(j),"°"
     call Unit_Tests_Begin_Group(trim(label))
     allocate(galacticInclination_)
     allocate(torusMisaligned_    )
     !![
     <referenceConstruct object="galacticInclination_" constructor="galacticInclinationMisaligned(anglesMisalignmentMaximum(j),8,16,galacticInclinationNull_)"                      />
     <referenceConstruct object="torusMisaligned_"     constructor="dustAttenuationAGNTorus(5.0d0,30.0d0,dustExtinctionCurve_,galacticInclination_,blackHoleTorusCoveringFraction_)"/>
     !!]
     call Assert("no inclination is available for each galaxy",galacticInclination_%isAvailable(),.false.)
     cosineMisalignmentMaximum=cos(anglesMisalignmentMaximum(j)*Pi/180.0d0)
     momentCosineMisalignment =(1.0d0+cosineMisalignmentMaximum+cosineMisalignmentMaximum**2)/3.0d0
     do i=0,3
        inclinationHost=dble(i)*Pi/6.0d0
        call galacticInclination_%inclinationsMisaligned(node,inclinationHost,inclinations,weights)
        write (label,'(a,i2,a)') "i=",30*i,"°: weights sum to one"
        call Assert(trim(label),sum(weights),1.0d0,relTol=1.0d-12)
        write (label,'(a,i2,a)') "i=",30*i,"°: <cos²θ> of the cone of misalignment"
        call Assert(trim(label),sum(weights*cos(inclinations)**2),cos(inclinationHost)**2*momentCosineMisalignment+sin(inclinationHost)**2*(1.0d0-momentCosineMisalignment)/2.0d0,absTol=1.0d-12)
        if (anglesMisalignmentMaximum(j) == 0.0d0) then
           write (label,'(a,i2,a)') "i=",30*i,"°: aligned torus is seen at the inclination of the galaxy"
           call Assert(trim(label),inclinations,[inclinationHost],absTol=1.0d-12)
        end if
        ! The misaligned torus transmits the weighted sum of the aligned torus at the inclinations returned.
        transmissionExpected=0.0d0
        do k=1,size(inclinations)
           transmissionExpected=transmissionExpected+weights(k)*torusAligned_%transmission(node,descriptors,inclinations(k))
        end do
        transmission=torusMisaligned_%transmission(node,descriptors,inclinationHost)
        write (label,'(a,i2,a)') "i=",30*i,"°: torus transmits the average over its orientations"
        call Assert(trim(label),transmission,transmissionExpected,relTol=1.0d-12)
        ! A torus whose axis is independent of its galaxy transmits the average over cos(θ) of the aligned torus. The
        ! tolerance is the accuracy of the default quadratures (8 points in cos(ψ), 16 in φ): the error was measured to
        ! be at most 9×10⁻⁵, edge-on.
        if (anglesMisalignmentMaximum(j) == 90.0d0) then
           do k=1,countFine
              transmissionExpected=torusAligned_%transmission(node,descriptors,acos(cosineInclinationFine(k)))
              integrand(k)        =transmissionExpected(1)
           end do
           write (label,'(a,i2,a)') "i=",30*i,"°: independent torus transmits the average over cos(θ)"
           call Assert(trim(label),transmission(1),trapezoid(cosineInclinationFine,integrand),relTol=2.0d-4)
        end if
     end do
     !![
     <objectDestructor name="torusMisaligned_"    />
     <objectDestructor name="galacticInclination_"/>
     !!]
     call Unit_Tests_End_Group()
  end do

  ! A galaxy with its own misalignment, seen at 60°, with a torus misaligned by ψ=60° (a deviate of 1/2, for ψ_max=90°).
  call Unit_Tests_Begin_Group("Misaligned inclination: for each galaxy")
  allocate(galacticInclination_)
  allocate(torusMisaligned_    )
  !![
  <referenceConstruct object="galacticInclination_" constructor="galacticInclinationMisaligned(90.0d0,8,16,galacticInclinationFixed_)"                                           />
  <referenceConstruct object="torusMisaligned_"     constructor="dustAttenuationAGNTorus(5.0d0,30.0d0,dustExtinctionCurve_,galacticInclination_,blackHoleTorusCoveringFraction_)"/>
  !!]
  call Assert("an inclination is available for each galaxy",galacticInclination_%isAvailable(),.true.)
  call basic%floatRank0MetaPropertySet(misalignmentDeviateID,0.5d0)
  call basic%floatRank0MetaPropertySet(misalignmentAzimuthID,0.0d0)
  call Assert("tilted towards the observer: seen at |i-ψ|"              ,galacticInclination_%inclination (node            ),                                            0.0d0 ,absTol=1.0d-6)
  call Assert("tilted towards the observer: torus transmits as at |i-ψ|",torusMisaligned_    %transmission(node,descriptors),torusAligned_%transmission(node,descriptors,0.0d0),relTol=1.0d-6)
  call basic%floatRank0MetaPropertySet(misalignmentAzimuthID,Pi   )
  call Assert("tilted away from the observer: seen at i+ψ, folded"    ,    galacticInclination_%inclination(node) ,Pi/3.0d0,absTol=1.0d-06)
  call basic%floatRank0MetaPropertySet(misalignmentAzimuthID,0.5d0*Pi)
  call Assert("tilted across the line of sight: cos θ = cos i cos ψ"  ,cos(galacticInclination_%inclination(node)),0.25d0  ,absTol=1.0d-12)
  ! Averaging over orientation uses the galaxy's own misalignment, and averages over azimuth only.
  call galacticInclination_%inclinationsMisaligned(node,0.0d0,inclinations,weights)
  call Assert("averaged face-on: the galaxy's own misalignment is used",     inclinations ,spread(Pi/3.0d0,1,size(inclinations)),absTol=1.0d-12)
  call Assert("averaged: points in azimuth only"                       ,size(inclinations),16                                                  )
  !![
  <objectDestructor name="torusMisaligned_"    />
  <objectDestructor name="galacticInclination_"/>
  !!]
  call Unit_Tests_End_Group()
  !![
  <objectDestructor name="torusAligned_"                  />
  <objectDestructor name="blackHoleTorusCoveringFraction_"/>
  <objectDestructor name="galacticInclinationFixed_"      />
  <objectDestructor name="galacticInclinationNull_"       />
  <objectDestructor name="dustExtinctionCurve_"           />
  !!]
  call Unit_Tests_Finish()

contains

  double precision function trapezoid(x,y)
    !!{RST
    Return the integral of ``y`` over ``x`` by the trapezoidal rule.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x, y

    trapezoid=sum(0.5d0*(y(2:)+y(:size(y)-1))*(x(2:)-x(:size(x)-1)))
    return
  end function trapezoid

end program Test_Dust_Attenuation_AGN_Torus_Misaligned
