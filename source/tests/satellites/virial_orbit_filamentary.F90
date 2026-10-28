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
  Contains a program to test the filamentary virial orbit class and the filaments of halos.
  !!}

program Test_Satellite_Virial_Orbit_Filamentary
  !!{RST
  Test the :galacticus-class:`virialOrbitFilamentary` class, and the :galacticus-class:`darkMatterHaloFilamentsRandomSequential`
  filaments on which it relies.

  The filaments are tested first. On a small tree, the filaments requested for a node part-way along a branch must be stored,
  identically, on every node of that branch, and must differ from those of another branch. They must be unit vectors, the
  first two must be opposite, and all must be separated by at least the minimum angle. For a non-integer mean number of
  filaments, the fraction of halos with the larger number of filaments must match its expectation, and the smallest
  separation found among many halos must approach the minimum angle (so that the separation test can not pass merely because
  the constraint is never active).

  The infall directions are then tested against the distribution from which they are drawn,

  .. math::

     p(\hat{\mathbf{n}}) = \frac{1-f}{4\pi} + \frac{f}{N} \sum_{i=1}^N \frac{\kappa}{4\pi \sinh\kappa}
     \exp\left(\kappa\, \hat{\mathbf{n}}\cdot\hat{\mathbf{f}}_i\right),

  for the filaments, :math:`\hat{\mathbf{f}}_i`, of the host. Directions are binned in :math:`\cos\theta` and :math:`\phi`, and
  the expected number in each bin is found by integrating :math:`p(\hat{\mathbf{n}})` numerically over the bin. The two are
  compared with a :math:`\chi^2` test, which must lie within five standard deviations of its expectation. This is done for a
  fixed filamentary fraction, for a zero fraction (which must give isotropic infall), and for the default, mass-dependent,
  fraction in hosts of two masses, the lesser of which has a fraction limited to zero.

  So that the test can not pass merely by being insensitive, the sample with a fixed fraction is also required to be
  inconsistent, by the same criterion, with isotropy, with a fraction 10% smaller than the true one, and with infall along
  only the first two of the filaments of the host.
  !!}
  use, intrinsic :: ISO_C_Binding             , only : c_long
  use            :: Cosmological_Density_Field, only : criticalOverdensityClass
  use            :: Dark_Matter_Halo_Filaments, only : darkMatterHaloFilamentsClass  , darkMatterHaloFilamentsRandomSequential
  use            :: Display                   , only : displayMessage                , displayVerbositySet                    , verbosityLevelStandard
  use            :: Error                     , only : Error_Handler_Register
  use            :: Events_Hooks              , only : eventsHooksInitialize
  use            :: Functions_Global_Utilities, only : Functions_Global_Set
  use            :: Galacticus_Nodes          , only : mergerTree                    , nodeClassHierarchyFinalize             , nodeClassHierarchyInitialize       , nodeComponentBasic          , &
       &                                               nodeComponentDarkMatterProfile, treeNode
  use            :: IO_HDF5                   , only : ioHDF5AccessInitialize
  use            :: Input_Parameters          , only : inputParameters
  use            :: Kepler_Orbits             , only : keplerOrbit
  use            :: Node_Components           , only : Node_Components_Initialize    , Node_Components_Thread_Initialize      , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use            :: Nodes_Operators           , only : nodeOperatorHaloFilaments
  use            :: Numerical_Constants_Math  , only : Pi
  use            :: Numerical_Random_Numbers  , only : randomNumberGeneratorGSL
  use            :: Unit_Tests                , only : Assert                        , Unit_Tests_Begin_Group                 , Unit_Tests_End_Group               , Unit_Tests_Finish
  use            :: Virial_Orbits             , only : virialOrbitClass              , virialOrbitFilamentary
  implicit none
  type            (inputParameters                        )                              :: parameters                         , parametersFilamentary
  type            (nodeOperatorHaloFilaments              ), pointer                     :: nodeOperator_
  class           (virialOrbitClass                       ), pointer                     :: virialOrbit_                       , virialOrbitWrapped_
  type            (virialOrbitFilamentary                 ), pointer                     :: virialOrbitFixed_                  , virialOrbitNone_
  class           (criticalOverdensityClass               ), pointer                     :: criticalOverdensity_
  class           (darkMatterHaloFilamentsClass           ), pointer                     :: darkMatterHaloFilaments_
  type            (darkMatterHaloFilamentsRandomSequential), pointer                     :: darkMatterHaloFilamentsMixed_
  type            (mergerTree                             )                              :: tree
  type            (treeNode                               ), pointer                     :: nodeLeaf                           , nodeMiddle                      , &
       &                                                                                    nodeRoot                           , nodeSecondary                   , &
       &                                                                                    nodeIsolated
  ! The epoch [Gyr] at which halos are placed.
  double precision                                         , parameter                   :: time                      =1.3800d1
  ! The default parameters of the filamentary class, and the concentration implied by its default filament radius.
  double precision                                         , parameter                   :: fractionCharacteristic    =0.2910d0, fractionLogarithmicSlope=0.102d0, &
       &                                                                                    radiusFilamentFractional  =0.4325d0
  ! The default number of filaments and their minimum separation [degrees].
  integer                                                  , parameter                   :: countFilamentsDefault     =3
  double precision                                         , parameter                   :: angleSeparationMinimum    =6.0d1
  ! The non-integer mean number of filaments used to test the distribution of the number of filaments, and the number of halos
  ! drawn in that test.
  double precision                                         , parameter                   :: countFilamentsMixed       =3.4d0
  integer                                                  , parameter                   :: countHalos                =20000
  ! The fixed fraction used to test the shape of the distribution.
  double precision                                         , parameter                   :: fractionFixed             =0.3000d0
  ! The number of orbits drawn in each test, the numbers of bins in cos(θ) and in φ, and the number of sub-intervals per bin in
  ! each used to integrate the expected distribution.
  integer                                                  , parameter                   :: countOrbits               =1000000 , countBinsCosine         =20     , &
       &                                                                                    countBinsAzimuth          =16      , countSubIntervals       =16
  double precision                                         , allocatable, dimension(:,:) :: axesMiddle                         , axesSecondary                   , &
       &                                                                                    axesHalo
  double precision                                         , allocatable, dimension(:  ) :: axesLeafStored                     , axesRootStored                  , &
       &                                                                                    axesIsolatedStored
  class           (nodeComponentBasic                     ), pointer                     :: basic
  double precision                                                                       :: concentration                      , massCollapsing                  , &
       &                                                                                    fractionExpected                   , cosineSeparationLargest         , &
       &                                                                                    cosineSeparationLargestAll         , fractionMixed                   , &
       &                                                                                    fractionMixedExpected
  integer                                                                                :: filamentAxesID                     , i                               , &
       &                                                                                    countHalosLarger
  logical                                                                                :: countsValid                        , separationsValid
  character       (len=256                                )                              :: message

  call displayVerbositySet(verbosityLevelStandard)
  call Error_Handler_Register()
  call ioHDF5AccessInitialize()
  call Unit_Tests_Begin_Group("Filamentary virial orbits")
  parameters=inputParameters('testSuite/parameters/satellites/virialOrbitFilamentary.xml')
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Build the filaments (with their default parameters), the node operator which assigns them, the filamentary orbit class from
  ! the parameters (with its default parameters), the class which it wraps, and the critical overdensity from which the
  ! collapsing mass is found.
  parametersFilamentary=parameters%subParameters('virialOrbit')
  !![
  <objectBuilder class="darkMatterHaloFilaments" name="darkMatterHaloFilaments_" source="parameters"           />
  <objectBuilder class="virialOrbit"             name="virialOrbit_"             source="parameters"           />
  <objectBuilder class="virialOrbit"             name="virialOrbitWrapped_"      source="parametersFilamentary"/>
  <objectBuilder class="criticalOverdensity"     name="criticalOverdensity_"     source="parameters"           />
  <addMetaProperty component="basic" name="filamentAxes" id="filamentAxesID" rank="1" isEvolvable="no"/>
  !!]
  allocate(nodeOperator_                )
  allocate(darkMatterHaloFilamentsMixed_)
  allocate(virialOrbitFixed_            )
  allocate(virialOrbitNone_             )
  nodeOperator_                =nodeOperatorHaloFilaments              (darkMatterHaloFilaments_                                                                                                                  )
  darkMatterHaloFilamentsMixed_=darkMatterHaloFilamentsRandomSequential(countFilamentsMixed     ,angleSeparationMinimum                                                                                           )
  virialOrbitFixed_            =virialOrbitFilamentary                 (fractionFixed           ,0.0d0                 ,radiusFilamentFractional,virialOrbitWrapped_,criticalOverdensity_,darkMatterHaloFilaments_)
  virialOrbitNone_             =virialOrbitFilamentary                 (0.0d0                   ,0.0d0                 ,radiusFilamentFractional,virialOrbitWrapped_,criticalOverdensity_,darkMatterHaloFilaments_)
  concentration                =1.0d0/asin(radiusFilamentFractional)**2
  ! Build a random number generator for the tree.
  allocate(randomNumberGeneratorGSL :: tree%randomNumberGenerator_)
  select type (randomNumberGenerator_ => tree%randomNumberGenerator_)
  type is (randomNumberGeneratorGSL)
     randomNumberGenerator_=randomNumberGeneratorGSL(seed_=8317_c_long)
  end select
  call tree%properties%initialize()
  ! Test the filaments on a small tree: a branch of three nodes, with a secondary progenitor joining at the top.
  call Unit_Tests_Begin_Group("Filaments")
  call buildNode(nodeRoot     ,1.0d13)
  call buildNode(nodeMiddle   ,5.0d12)
  call buildNode(nodeLeaf     ,1.0d12)
  call buildNode(nodeSecondary,1.0d12)
  tree%nodeBase            => nodeRoot
  nodeRoot     %firstChild => nodeMiddle
  nodeMiddle   %parent     => nodeRoot
  nodeMiddle   %sibling    => nodeSecondary
  nodeSecondary%parent     => nodeRoot
  nodeMiddle   %firstChild => nodeLeaf
  nodeLeaf     %parent     => nodeMiddle
  axesMiddle    =  darkMatterHaloFilaments_%axes                     (nodeMiddle    )
  basic         => nodeLeaf                %basic                    (              )
  axesLeafStored=  basic                   %floatRank1MetaPropertyGet(filamentAxesID)
  basic         => nodeRoot                %basic                    (              )
  axesRootStored=  basic                   %floatRank1MetaPropertyGet(filamentAxesID)
  axesSecondary =  darkMatterHaloFilaments_%axes                     (nodeSecondary )
  call Assert('default number of filaments'                             ,shape(axesMiddle)                                                         ,[3,countFilamentsDefault]                                   )
  call Assert('filaments are unit vectors'                              ,sqrt(sum(axesMiddle**2,dim=1))                                            ,spread(1.0d0,1,countFilamentsDefault)        ,relTol=1.0d-12)
  call Assert('first two filaments are opposite'                        ,axesMiddle(:,2)                                                           ,-axesMiddle(:,1)                                            )
  call Assert('filaments are separated by the minimum angle'            ,largestCosine(axesMiddle) <= cos(angleSeparationMinimum*Pi/180.0d0)       ,.true.                                                      )
  call Assert('filaments are stored on the earliest progenitor'         ,axesLeafStored                                                            ,reshape(axesMiddle,[3*countFilamentsDefault])               )
  call Assert('filaments are stored on the node the branch joins'       ,axesRootStored                                                            ,reshape(axesMiddle,[3*countFilamentsDefault])               )
  call Assert('filaments are unchanged when requested again'            ,reshape(darkMatterHaloFilaments_%axes(nodeRoot),[3*countFilamentsDefault]),reshape(axesMiddle,[3*countFilamentsDefault])               )
  call Assert('filaments of a secondary progenitor are drawn separately',abs(dot_product(axesSecondary(:,1),axesMiddle(:,1))) < 0.999999d0         ,.true.                                                      )
  call buildNode(nodeIsolated,1.0d12)
  call nodeOperator_%nodeTreeInitialize(nodeIsolated)
  basic             => nodeIsolated%basic                    (              )
  axesIsolatedStored=  basic       %floatRank1MetaPropertyGet(filamentAxesID)
  call Assert('node operator assigns filaments',size(axesIsolatedStored),3*countFilamentsDefault)
  call nodeIsolated%destroy()
  ! Test the distribution of the number of filaments, and their separations, over many halos.
  countHalosLarger          =0
  countsValid               =.true.
  separationsValid          =.true.
  cosineSeparationLargestAll=-1.0d0
  do i=1,countHalos
     call buildNode(nodeIsolated,1.0d12)
     axesHalo=darkMatterHaloFilamentsMixed_%axes(nodeIsolated)
     if      (size(axesHalo,dim=2) == int(countFilamentsMixed)+1) then
        countHalosLarger=countHalosLarger+1
     else if (size(axesHalo,dim=2) /= int(countFilamentsMixed)  ) then
        countsValid     =.false.
     end if
     cosineSeparationLargest   =largestCosine(axesHalo)
     cosineSeparationLargestAll=max(cosineSeparationLargestAll,cosineSeparationLargest)
     if (cosineSeparationLargest > cos(angleSeparationMinimum*Pi/180.0d0)) separationsValid=.false.
     call nodeIsolated%destroy()
  end do
  fractionMixed        =dble(countHalosLarger)/dble(countHalos)
  fractionMixedExpected=countFilamentsMixed-dble(int(countFilamentsMixed))
  write (message,'(a,f6.4,a,f6.4,a,f6.2,a)') '  fraction with more filaments = ',fractionMixed,' (expected ',fractionMixedExpected,'); smallest separation = ',acos(cosineSeparationLargestAll)*180.0d0/Pi,'°'
  call displayMessage(trim(message))
  call Assert('number of filaments brackets the mean'                    ,countsValid                                                                                                                   ,.true.)
  call Assert('fraction of halos with more filaments'                    ,abs(fractionMixed-fractionMixedExpected)    < 5.0d0*sqrt(fractionMixedExpected*(1.0d0-fractionMixedExpected)/dble(countHalos)),.true.)
  call Assert('filaments of all halos are separated by the minimum angle',separationsValid                                                                                                              ,.true.)
  call Assert('minimum separation is approached'                         ,acos(cosineSeparationLargestAll)*180.0d0/Pi < angleSeparationMinimum+1.0d0                                                    ,.true.)
  call Unit_Tests_End_Group()
  ! Test the distribution of infall directions.
  call Unit_Tests_Begin_Group("Infall directions")
  call testDirections('fixed filamentary fraction'                   ,virialOrbitFixed_,1.0d13               ,fractionFixed   ,checkSensitivity=.true. )
  call testDirections('zero filamentary fraction (isotropic)'        ,virialOrbitNone_ ,1.0d13               ,0.0d0           ,checkSensitivity=.false.)
  massCollapsing  =criticalOverdensity_%collapsingMass(time=time)
  fractionExpected=fractionCharacteristic+fractionLogarithmicSlope*log10(10.0d0)
  call testDirections('default fraction, M = 10 M*'                  ,virialOrbit_     ,1.0d+1*massCollapsing,fractionExpected,checkSensitivity=.false.)
  call testDirections('default fraction limited to zero, M = 10⁻⁴ M*',virialOrbit_     ,1.0d-4*massCollapsing,0.0d0           ,checkSensitivity=.false.)
  call Unit_Tests_End_Group()
  ! Clean up.
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()

contains

  subroutine testDirections(label,virialOrbitTest_,massHost,fraction,checkSensitivity)
    !!{RST
    Draw orbits for a satellite of a host of the given mass, and compare the distribution of the directions of infall with that
    expected for the filaments of the host and the given filamentary fraction.
    !!}
    implicit none
    character       (len=*           ), intent(in   )                                               :: label
    class           (virialOrbitClass), intent(inout)                                               :: virialOrbitTest_
    double precision                  , intent(in   )                                               :: massHost        , fraction
    logical                           , intent(in   )                                               :: checkSensitivity
    type            (treeNode        ), pointer                                                     :: host            , satellite
    double precision                  , allocatable  , dimension(:,:)                               :: axes
    integer                                          , dimension(countBinsCosine,countBinsAzimuth) :: counts
    type            (keplerOrbit     )                                                              :: orbit
    double precision                                                                                :: chiSquaredTrue  , chiSquaredIsotropic, &
         &                                                                                             chiSquaredBiased, chiSquaredPair
    integer                                                                                         :: i               , iCosine            , &
         &                                                                                             iAzimuth
    character       (len=256         )                                                              :: message

    call Unit_Tests_Begin_Group(label)
    ! Build a host and a satellite, and find the filaments of the host.
    call buildNode(host     ,       massHost)
    call buildNode(satellite,1.0d-2*massHost)
    tree%nodeBase => host
    axes=darkMatterHaloFilaments_%axes(host)
    ! Draw orbits, and accumulate the distribution of directions in bins of cos(θ) and φ.
    counts=0
    do i=1,countOrbits
       orbit   =virialOrbitTest_%orbit(satellite,host,acceptUnboundOrbits=.true.)
       iCosine =min(int((cos(orbit%theta())+1.0d0)*dble(countBinsCosine )/2.0d0   )+1,countBinsCosine )
       iAzimuth=min(int(     orbit%phi  ()        *dble(countBinsAzimuth)/2.0d0/Pi)+1,countBinsAzimuth)
       counts(iCosine,iAzimuth)=counts(iCosine,iAzimuth)+1
    end do
    chiSquaredTrue=chiSquared(counts,expectedCounts(axes,fraction))
    write (message,'(a,f8.2,a,i3,a)') '  χ² = ',chiSquaredTrue,' for ',countBinsCosine*countBinsAzimuth,' bins'
    call displayMessage(trim(message))
    call Assert('distribution of directions matches',chiSquaredTrue < chiSquaredLimit(),.true.)
    if (checkSensitivity) then
       chiSquaredIsotropic=chiSquared(counts,expectedCounts(axes      ,0.0d0         ))
       chiSquaredBiased   =chiSquared(counts,expectedCounts(axes      ,0.9d0*fraction))
       chiSquaredPair     =chiSquared(counts,expectedCounts(axes(:,1:2),     fraction))
       write (message,'(a,f12.2,a,f12.2,a,f12.2)') '  χ² against isotropy = ',chiSquaredIsotropic,', against 0.9 f = ',chiSquaredBiased,', against the first two filaments only = ',chiSquaredPair
       call displayMessage(trim(message))
       call Assert('test distinguishes isotropic infall'            ,chiSquaredIsotropic > chiSquaredLimit(),.true.)
       call Assert('test distinguishes a 10% smaller fraction'      ,chiSquaredBiased    > chiSquaredLimit(),.true.)
       call Assert('test distinguishes infall along fewer filaments',chiSquaredPair      > chiSquaredLimit(),.true.)
    end if
    call host     %destroy()
    call satellite%destroy()
    call Unit_Tests_End_Group()
    return
  end subroutine testDirections

  function expectedCounts(axes,fraction)
    !!{RST
    Return the expected number of directions in each bin of :math:`\cos\theta` and :math:`\phi`, by integrating the distribution
    of directions over each bin with the midpoint rule on a grid of sub-intervals.
    !!}
    implicit none
    double precision, dimension(countBinsCosine,countBinsAzimuth) :: expectedCounts
    double precision, dimension(:,:), intent(in   )               :: axes
    double precision                , intent(in   )               :: fraction
    double precision, dimension(3  )                              :: direction
    double precision                                              :: cosineTheta   , phi     , &
         &                                                           density
    integer                                                       :: iCosine       , iAzimuth, &
         &                                                           jCosine       , jAzimuth

    expectedCounts=0.0d0
    do iCosine=1,countBinsCosine
       do iAzimuth=1,countBinsAzimuth
          do jCosine=1,countSubIntervals
             cosineTheta=-1.0d0+2.0d0   *(dble((iCosine -1)*countSubIntervals+jCosine )-0.5d0)/dble(countBinsCosine *countSubIntervals)
             do jAzimuth=1,countSubIntervals
                phi      =      2.0d0*Pi*(dble((iAzimuth-1)*countSubIntervals+jAzimuth)-0.5d0)/dble(countBinsAzimuth*countSubIntervals)
                direction=[                                     &
                     &     sqrt(1.0d0-cosineTheta**2)*cos(phi), &
                     &     sqrt(1.0d0-cosineTheta**2)*sin(phi), &
                     &                cosineTheta               &
                     &    ]
                ! Evaluate the density, writing each von Mises-Fisher kernel in a form which can not overflow.
                density  =+(1.0d0-fraction)/4.0d0/Pi                                &
                     &    +       fraction /dble(size(axes,dim=2))                  &
                     &    *concentration/2.0d0/Pi/(1.0d0-exp(-2.0d0*concentration)) &
                     &    *sum(exp(concentration*(matmul(direction,axes)-1.0d0)))
                expectedCounts(iCosine,iAzimuth)=+expectedCounts(iCosine,iAzimuth)  &
                     &                           +density
             end do
          end do
       end do
    end do
    ! Multiply by the solid angle of each sub-interval, and by the number of orbits.
    expectedCounts=+expectedCounts                                      &
         &         *2.0d0     /dble(countBinsCosine *countSubIntervals) &
         &         *2.0d0*Pi  /dble(countBinsAzimuth*countSubIntervals) &
         &         *dble(countOrbits)
    return
  end function expectedCounts

  double precision function largestCosine(axes)
    !!{RST
    Return the largest cosine of the angle between any two of the given directions, i.e. the cosine of their smallest separation.
    !!}
    implicit none
    double precision, dimension(:,:), intent(in   ) :: axes
    integer                                         :: i   , j

    largestCosine=-1.0d0
    do i=1,size(axes,dim=2)
       do j=i+1,size(axes,dim=2)
          largestCosine=max(largestCosine,dot_product(axes(:,i),axes(:,j)))
       end do
    end do
    return
  end function largestCosine

  double precision function chiSquared(counts,expected)
    !!{RST
    Return the :math:`\chi^2` of binned counts relative to their expectations.
    !!}
    implicit none
    integer         , dimension(:,:), intent(in   ) :: counts
    double precision, dimension(:,:), intent(in   ) :: expected

    chiSquared=sum((dble(counts)-expected)**2/expected)
    return
  end function chiSquared

  double precision function chiSquaredLimit()
    !!{RST
    Return the largest acceptable :math:`\chi^2`: five standard deviations above its expectation. The constraint that the counts
    sum to the number of orbits is ignored, which makes this limit slightly more permissive.
    !!}
    implicit none
    integer :: countBins

    countBins      =countBinsCosine*countBinsAzimuth
    chiSquaredLimit=dble(countBins)+5.0d0*sqrt(2.0d0*dble(countBins))
    return
  end function chiSquaredLimit

  subroutine buildNode(node_,mass_)
    !!{RST
    Build a node of the given mass, with a dark matter profile scale length following a plausible mass-concentration relation.
    !!}
    implicit none
    type            (treeNode                      ), intent(  out), pointer :: node_
    double precision                                , intent(in   )          :: mass_
    class           (nodeComponentBasic            )               , pointer :: basic_
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile_

    node_             => treeNode(hostTree=tree)
    basic_            => node_   %basic            (autoCreate=.true.)
    darkMatterProfile_=> node_   %darkMatterProfile(autoCreate=.true.)
    call basic_            %massSet            (mass_                               )
    call basic_            %timeSet            (time                                )
    call basic_            %timeLastIsolatedSet(time                                )
    call darkMatterProfile_%scaleSet           (2.0d-2*(mass_/1.0d12)**(1.0d0/3.0d0))
    return
  end subroutine buildNode

end program Test_Satellite_Virial_Orbit_Filamentary
