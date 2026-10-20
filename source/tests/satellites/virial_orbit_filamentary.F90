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
  Contains a program to test the filamentary virial orbit class and the filament axes of halos.
  !!}

program Test_Satellite_Virial_Orbit_Filamentary
  !!{RST
  Test the :galacticus-class:`virialOrbitFilamentary` class, and the filament axes on which it relies.

  The filament axes are tested first: on a small tree, the axes requested for a node part-way along a branch must be stored,
  identically, on every node of that branch, must form an opposite pair of unit vectors, and must differ from those of another
  branch.

  The infall directions are then tested against the distribution from which they are drawn. For a host with a filament along
  :math:`\hat{\mathbf{a}}`, the cosine of the angle between the direction of infall and the filament, :math:`\mu`, has the
  distribution

  .. math::

     P(|\mu|) = (1-f) + f \kappa \frac{\cosh(\kappa|\mu|)}{\sinh\kappa},

  whose cumulative distribution, :math:`(1-f)|\mu| + f \sinh(\kappa|\mu|)/\sinh\kappa`, gives the expected number of directions
  in each bin of :math:`|\mu|`. The azimuth of the direction about the filament must be uniformly distributed. Each is compared
  with a :math:`\chi^2` test, which must lie within five standard deviations of its expectation. This is done for a fixed
  filamentary fraction, for a zero fraction (which must give isotropic infall), and for the default, mass-dependent, fraction
  in hosts of two masses, the lesser of which has a fraction limited to zero.

  So that the test can not pass merely by being insensitive, the sample with a fixed fraction is also required to be
  inconsistent, by the same criterion, both with isotropy and with a fraction 10% smaller than the true one.
  !!}
  use, intrinsic :: ISO_C_Binding             , only : c_long
  use            :: Cosmological_Density_Field, only : criticalOverdensityClass
  use            :: Dark_Matter_Halo_Filaments, only : Dark_Matter_Halo_Filament_Axes
  use            :: Display                   , only : displayMessage                , displayVerbositySet              , verbosityLevelStandard
  use            :: Error                     , only : Error_Handler_Register
  use            :: Events_Hooks              , only : eventsHooksInitialize
  use            :: Functions_Global_Utilities, only : Functions_Global_Set
  use            :: Galacticus_Nodes          , only : mergerTree                    , nodeClassHierarchyFinalize       , nodeClassHierarchyInitialize       , nodeComponentBasic          , &
       &                                               nodeComponentDarkMatterProfile, treeNode
  use            :: IO_HDF5                   , only : ioHDF5AccessInitialize
  use            :: Input_Parameters          , only : inputParameters
  use            :: Kepler_Orbits             , only : keplerOrbit
  use            :: Node_Components           , only : Node_Components_Initialize    , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use            :: Nodes_Operators           , only : nodeOperatorHaloFilaments
  use            :: Numerical_Constants_Math  , only : Pi
  use            :: Numerical_Random_Numbers  , only : randomNumberGeneratorGSL
  use            :: Unit_Tests                , only : Assert                        , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  use            :: Vectors                   , only : Vector_Product
  use            :: Virial_Orbits             , only : virialOrbitClass              , virialOrbitFilamentary
  implicit none
  type            (inputParameters          )                              :: parameters                         , parametersFilamentary
  type            (nodeOperatorHaloFilaments), pointer                     :: nodeOperator_
  class           (virialOrbitClass         ), pointer                     :: virialOrbit_                       , virialOrbitWrapped_
  type            (virialOrbitFilamentary   ), pointer                     :: virialOrbitFixed_                  , virialOrbitNone_
  class           (criticalOverdensityClass ), pointer                     :: criticalOverdensity_
  type            (mergerTree               )                              :: tree
  type            (treeNode                 ), pointer                     :: nodeLeaf                           , nodeMiddle                      , &
       &                                                                      nodeRoot                           , nodeSecondary
  ! The epoch [Gyr] at which halos are placed.
  double precision                           , parameter                   :: time                      =1.3800d1
  ! The default parameters of the filamentary class, and the concentration implied by its default filament radius.
  double precision                           , parameter                   :: fractionCharacteristic    =0.1840d0, fractionLogarithmicSlope=0.055d0, &
       &                                                                      radiusFilamentFractional  =0.4546d0
  ! The fixed fraction used to test the shape of the distribution.
  double precision                           , parameter                   :: fractionFixed             =0.3000d0
  ! The number of orbits drawn in each test, and the numbers of bins in |μ| and in azimuth.
  integer                                    , parameter                   :: countOrbits               =1000000 , countBinsCosine         =20     , &
       &                                                                      countBinsAzimuth          =16
  double precision                           , allocatable, dimension(:,:) :: axesMiddle                         , axesSecondary
  double precision                           , allocatable, dimension(:  ) :: axesLeafStored                     , axesRootStored
  class           (nodeComponentBasic       ), pointer                     :: basic
  double precision                                                         :: concentration                      , massCollapsing                  , &
       &                                                                      fractionExpected
  integer                                                                  :: filamentAxesID

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
  ! Build the node operator which creates the filament axes, the filamentary orbit class from the parameters (with its default
  ! parameters), the class which it wraps, and the critical overdensity from which the collapsing mass is found.
  allocate(nodeOperator_)
  nodeOperator_=nodeOperatorHaloFilaments()
  parametersFilamentary=parameters%subParameters('virialOrbit')
  !![
  <objectBuilder class="virialOrbit"         name="virialOrbit_"         source="parameters"           />
  <objectBuilder class="virialOrbit"         name="virialOrbitWrapped_"  source="parametersFilamentary"/>
  <objectBuilder class="criticalOverdensity" name="criticalOverdensity_" source="parameters"           />
  <addMetaProperty component="basic" name="filamentAxes" id="filamentAxesID" rank="1" isEvolvable="no"/>
  !!]
  allocate(virialOrbitFixed_)
  allocate(virialOrbitNone_ )
  virialOrbitFixed_=virialOrbitFilamentary(fractionFixed,0.0d0,radiusFilamentFractional,virialOrbitWrapped_,criticalOverdensity_)
  virialOrbitNone_ =virialOrbitFilamentary(0.0d0        ,0.0d0,radiusFilamentFractional,virialOrbitWrapped_,criticalOverdensity_)
  concentration    =1.0d0/asin(radiusFilamentFractional)**2
  ! Build a random number generator for the tree.
  allocate(randomNumberGeneratorGSL :: tree%randomNumberGenerator_)
  select type (randomNumberGenerator_ => tree%randomNumberGenerator_)
  type is (randomNumberGeneratorGSL)
     randomNumberGenerator_=randomNumberGeneratorGSL(seed_=8317_c_long)
  end select
  call tree%properties%initialize()
  ! Test the filament axes on a small tree: a branch of three nodes, with a secondary progenitor joining at the top.
  call Unit_Tests_Begin_Group("Filament axes")
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
  axesMiddle    =  Dark_Matter_Halo_Filament_Axes    (nodeMiddle    ,filamentAxesID)
  basic         => nodeLeaf%basic                    (                             )
  axesLeafStored=  basic   %floatRank1MetaPropertyGet(filamentAxesID               )
  basic         => nodeRoot%basic                    (                             )
  axesRootStored=  basic   %floatRank1MetaPropertyGet(filamentAxesID               )
  axesSecondary =  Dark_Matter_Halo_Filament_Axes    (nodeSecondary ,filamentAxesID)
  call Assert('one filament, as a pair of axes'                     ,shape(axesMiddle)                                                   ,[3,2]                                 )
  call Assert('axes are unit vectors'                               ,sqrt(sum(axesMiddle**2,dim=1))                                      ,[1.0d0,1.0d0]          ,relTol=1.0d-12)
  call Assert('axes are opposite'                                   ,axesMiddle(:,2)                                                     ,-axesMiddle(:,1)                      )
  call Assert('axes are stored on the earliest progenitor'          ,axesLeafStored                                                      ,reshape(axesMiddle,[6])               )
  call Assert('axes are stored on the node the branch joins'        ,axesRootStored                                                      ,reshape(axesMiddle,[6])               )
  call Assert('axes are unchanged when requested again'             ,reshape(Dark_Matter_Halo_Filament_Axes(nodeRoot,filamentAxesID),[6]),reshape(axesMiddle,[6])               )
  call Assert('axes of a secondary progenitor are drawn separately' ,abs(dot_product(axesSecondary(:,1),axesMiddle(:,1))) < 0.999999d0   ,.true.                                )
  call Unit_Tests_End_Group()
  ! Test the distribution of infall directions.
  call Unit_Tests_Begin_Group("Infall directions")
  call testDirections('fixed filamentary fraction'                   ,virialOrbitFixed_,1.0d13,fractionFixed,checkSensitivity=.true. )
  call testDirections('zero filamentary fraction (isotropic)'        ,virialOrbitNone_ ,1.0d13,0.0d0        ,checkSensitivity=.false.)
  massCollapsing  =criticalOverdensity_%collapsingMass(time=time)
  fractionExpected=fractionCharacteristic+fractionLogarithmicSlope*log10(10.0d0)
  call testDirections('default fraction, M = 10 M*'                  ,virialOrbit_     ,1.0d1*massCollapsing,fractionExpected,checkSensitivity=.false.)
  call testDirections('default fraction limited to zero, M = 10⁻⁴ M*',virialOrbit_     ,1.0d-4*massCollapsing,0.0d0          ,checkSensitivity=.false.)
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
    Draw orbits for a satellite of a host of the given mass, and compare the distribution of the directions of infall, relative
    to the filament axis of the host, with that expected for the given filamentary fraction.
    !!}
    implicit none
    character       (len=*           ), intent(in   )                              :: label
    class           (virialOrbitClass), intent(inout)                              :: virialOrbitTest_
    double precision                  , intent(in   )                              :: massHost           , fraction
    logical                           , intent(in   )                              :: checkSensitivity
    type            (treeNode        ), pointer                                    :: host               , satellite
    double precision                  , allocatable  , dimension(:,:)              :: axes
    double precision                                 , dimension(3  )              :: axis               , tangent          , &
         &                                                                            binormal           , direction
    integer                                          , dimension(countBinsCosine ) :: countsCosine
    integer                                          , dimension(countBinsAzimuth) :: countsAzimuth
    double precision                                 , dimension(countBinsAzimuth) :: expectedAzimuth
    type            (keplerOrbit     )                                             :: orbit
    double precision                                                               :: chiSquaredCosine   , chiSquaredAzimuth, &
         &                                                                            chiSquaredIsotropic, chiSquaredBiased , &
         &                                                                            cosine             , azimuth
    integer                                                                        :: i
    character       (len=256         )                                             :: message

    call Unit_Tests_Begin_Group(label)
    ! Build a host and a satellite, and find the filament axis of the host, with an orthonormal basis about it.
    call buildNode(host     ,       massHost)
    call buildNode(satellite,1.0d-2*massHost)
    tree%nodeBase => host
    axes=Dark_Matter_Halo_Filament_Axes(host,filamentAxesID)
    axis=axes(:,1)
    if (abs(axis(1)) < 0.9d0) then
       tangent=Vector_Product(axis,[1.0d0,0.0d0,0.0d0])
    else
       tangent=Vector_Product(axis,[0.0d0,1.0d0,0.0d0])
    end if
    tangent =tangent/sqrt(sum(tangent**2))
    binormal=Vector_Product(axis,tangent)
    ! Draw orbits, and accumulate the distributions of |μ| and of azimuth about the axis.
    countsCosine =0
    countsAzimuth=0
    do i=1,countOrbits
       orbit    =virialOrbitTest_%orbit(satellite,host,acceptUnboundOrbits=.true.)
       direction=[                                     &
            &     sin(orbit%theta())*cos(orbit%phi()), &
            &     sin(orbit%theta())*sin(orbit%phi()), &
            &     cos(orbit%theta())                   &
            &    ]
       cosine   =abs(dot_product(direction,axis))
       azimuth  =modulo(atan2(dot_product(direction,binormal),dot_product(direction,tangent)),2.0d0*Pi)
       countsCosine (min(int(cosine *dble(countBinsCosine )/1.0d0   )+1,countBinsCosine )) &
            & =countsCosine (min(int(cosine *dble(countBinsCosine )/1.0d0   )+1,countBinsCosine ))+1
       countsAzimuth(min(int(azimuth*dble(countBinsAzimuth)/2.0d0/Pi)+1,countBinsAzimuth)) &
            & =countsAzimuth(min(int(azimuth*dble(countBinsAzimuth)/2.0d0/Pi)+1,countBinsAzimuth))+1
    end do
    expectedAzimuth  =dble(countOrbits)/dble(countBinsAzimuth)
    chiSquaredCosine =chiSquared(countsCosine ,expectedCosine (fraction))
    chiSquaredAzimuth=chiSquared(countsAzimuth,expectedAzimuth          )
    write (message,'(a,f8.2,a,f8.2,a,i2,a,i2,a)') '  χ²(|μ|) = ',chiSquaredCosine,', χ²(azimuth) = ',chiSquaredAzimuth,' for ',countBinsCosine,' and ',countBinsAzimuth,' bins'
    call displayMessage(trim(message))
    call Assert('distribution of |μ| matches'        ,chiSquaredCosine  < chiSquaredLimit(countBinsCosine ),.true.)
    call Assert('distribution of azimuth is uniform' ,chiSquaredAzimuth < chiSquaredLimit(countBinsAzimuth),.true.)
    if (checkSensitivity) then
       chiSquaredIsotropic=chiSquared(countsCosine,expectedCosine(0.0d0         ))
       chiSquaredBiased   =chiSquared(countsCosine,expectedCosine(0.9d0*fraction))
       write (message,'(a,f12.2,a,f12.2)') '  χ²(|μ|) against isotropy = ',chiSquaredIsotropic,', against 0.9 f = ',chiSquaredBiased
       call displayMessage(trim(message))
       call Assert('test distinguishes isotropic infall'      ,chiSquaredIsotropic > chiSquaredLimit(countBinsCosine),.true.)
       call Assert('test distinguishes a 10% smaller fraction',chiSquaredBiased    > chiSquaredLimit(countBinsCosine),.true.)
    end if
    call host     %destroy()
    call satellite%destroy()
    call Unit_Tests_End_Group()
    return
  end subroutine testDirections

  function expectedCosine(fraction)
    !!{RST
    Return the expected number of directions in each bin of :math:`|\mu|`, from the cumulative distribution
    :math:`(1-f)|\mu| + f \sinh(\kappa|\mu|)/\sinh\kappa`.
    !!}
    implicit none
    double precision, dimension(countBinsCosine)   :: expectedCosine
    double precision, intent(in   )                :: fraction
    double precision, dimension(0:countBinsCosine) :: cumulative
    integer                                        :: i

    do i=0,countBinsCosine
       cumulative(i)=+(1.0d0-fraction)*dble(i)/dble(countBinsCosine)                                         &
            &        +       fraction *sinh(concentration*dble(i)/dble(countBinsCosine))/sinh(concentration)
    end do
    expectedCosine=dble(countOrbits)*(cumulative(1:countBinsCosine)-cumulative(0:countBinsCosine-1))
    return
  end function expectedCosine

  double precision function chiSquared(counts,expected)
    !!{RST
    Return the :math:`\chi^2` of binned counts relative to their expectations.
    !!}
    implicit none
    integer         , dimension(:), intent(in   ) :: counts
    double precision, dimension(:), intent(in   ) :: expected

    chiSquared=sum((dble(counts)-expected)**2/expected)
    return
  end function chiSquared

  double precision function chiSquaredLimit(countBins)
    !!{RST
    Return the largest acceptable :math:`\chi^2` for the given number of bins: five standard deviations above its expectation.
    The constraint that the counts sum to the number of orbits is ignored, which makes this limit slightly more permissive.
    !!}
    implicit none
    integer, intent(in   ) :: countBins

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
