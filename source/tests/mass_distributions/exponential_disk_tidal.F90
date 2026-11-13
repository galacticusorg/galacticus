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
Contains a program to test the tidally truncated and heated exponential disk mass distribution.
!!}

program Test_Mass_Distributions_Exponential_Disk_Tidal
  !!{RST
  Tests the :galacticus-class:`massDistributionExponentialDiskTidal` mass distribution.
  !!}
  use :: Coordinates                          , only : coordinateCylindrical                 , coordinateCartesian                 , assignment(=)
  use :: Display                              , only : displayVerbositySet                   , verbosityLevelStandard
  use :: Events_Hooks                         , only : eventsHooksInitialize
  use :: Galactic_Component_Truncations       , only : galacticComponentTruncationExponential
  use :: IO_HDF5                              , only : ioHDF5AccessInitialize
  use :: Mass_Distributions                   , only : massDistributionExponentialDisk       , massDistributionExponentialDiskTidal
  use :: Numerical_Constants_Math             , only : Pi
  use :: Numerical_Integration_Tabulated      , only : integrateSimpson
  use :: Tensors                              , only : assignment(=)
  use :: Test_Exponential_Disk_Tidal_Reference, only : referenceRotationCurveSquared         , referencePotential                  , referenceRadiusTruncate, referenceFieldOffPlane
  use :: Unit_Tests                           , only : Assert                                , Unit_Tests_Begin_Group              , Unit_Tests_End_Group   , Unit_Tests_Finish
  implicit none
  double precision                                      , parameter                     :: scaleHeight      =0.137d0                         , radiusDecayOverRadiusTruncate=0.1d0, &
       &                                                                                   stepFractional   =1.0d-4
  double precision                                      , dimension( 5  )               :: radii            =[0.1d0,0.5d0,1.0d0,2.0d0,5.0d0]
  double precision                                      , dimension( 3  )               :: radiiTruncated   =[0.5d0,1.0d0,2.0d0]
  double precision                                      , dimension( 3  )               :: anglesOffPlane   =[0.0d0,45.0d0,80.0d0]
  double precision                                      , dimension( 2  )               :: radiiBlend       =[1.5d0,0.4d0]
  character       (len=2)                               , dimension( 3  )               :: labelsOffPlane   =['0 ','45','80']
  double precision                                      , dimension( 5  )               :: valuesTidal                                       , valuesParent
  double precision                                      , dimension( 3  )               :: valuesReference                                   , valuesTest                         , &
       &                                                                                   accelerationPlus                                  , accelerationMinus
  double precision                                      , dimension( 3,3)               :: tidalTensor                                       , tidalTensorDifference
  double precision                                      , dimension( 3  )               :: position                                          , positionPlus                       , &
       &                                                                                   positionMinus
  type            (massDistributionExponentialDisk     )                                :: parent
  type            (massDistributionExponentialDiskTidal)                                :: disk                                              , diskPlus                           , &
       &                                                                                   diskMinus
  type            (galacticComponentTruncationExponential)                              :: truncation
  double precision                                      , dimension( 5  )               :: truncations      =[0.01d0,0.1d0,0.5d0,2.0d0,8.0d0]
  double precision                                      , dimension( 2  )               :: truncationsInner =[1.0d0,3.0d0]
  double precision                                      , dimension( 3  )               :: radiiInner       =[0.3d0,0.6d0,1.0d0]
  double precision                                      , dimension( 4  )               :: heightsInner     =[0.5d0,1.0d0,2.0d0,4.0d0]
  double precision                                      , dimension(12  )               :: potentialInner                                    , potentialInnerReference            , &
       &                                                                                   radialInner                                       , radialInnerReference               , &
       &                                                                                   verticalInner                                     , verticalInnerReference
  double precision                                      , dimension( 2  )               :: heatingsInner    =[0.1d0,0.3d0]
  double precision                                      , dimension(12,2)               :: potentialHeatedReference                         , radialHeatedReference              , &
       &                                                                                   verticalHeatedReference
  character       (len=32)                                                              :: labelInner
  integer                                               , parameter                     :: countIntegration =400
  double precision                                      , dimension(0:countIntegration) :: logarithmTruncationIntegration     , integrandIntegration
  double precision                                                                      :: radiusTruncate                                    , radius                             , &
       &                                                                                   jumpMaximum                                       , massMaximum                        , &
       &                                                                                   quadrupole                                        , angle                              , &
       &                                                                                   potentialPrevious                                 , accelerationPrevious               , &
       &                                                                                   potentialCurrent                                  , accelerationCurrent
  integer                                                                               :: i                                                 , j                                  , &
       &                                                                                   k                                                 , l

  call displayVerbositySet   (verbosityLevelStandard)
  call eventsHooksInitialize ()
  call ioHDF5AccessInitialize()
  call Unit_Tests_Begin_Group("Tidally-modified exponential disk")
  parent=massDistributionExponentialDisk(scaleHeight=scaleHeight,dimensionless=.true.)
  ! Negligible truncation and no heating should recover the exponential disk.
  call Unit_Tests_Begin_Group("Untruncated, unheated limit")
  disk=massDistributionExponentialDiskTidal(truncation=1.0d-8,heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
  do i=1,size(radii)
     valuesTidal (i)=disk  %rotationCurve(radii(i))
     valuesParent(i)=parent%rotationCurve(radii(i))
  end do
  call Assert('rotation curve',valuesTidal,valuesParent,relTol=1.0d-5)
  do i=1,size(radii)
     valuesTidal (i)=disk  %massEnclosedByCylinder(radii(i))
     valuesParent(i)=parent%massEnclosedByCylinder(radii(i))
  end do
  call Assert('cylindrical mass',valuesTidal,valuesParent,relTol=1.0d-5)
  do i=1,size(radii)
     valuesTidal (i)=disk  %surfaceDensity(coordinateCylindrical(radii(i),0.0d0,0.0d0))
     valuesParent(i)=parent%surfaceDensity(coordinateCylindrical(radii(i),0.0d0,0.0d0))
  end do
  call Assert('surface density' ,valuesTidal,valuesParent,relTol=1.0d-5)
  call Assert('half-mass radius',disk%radiusHalfMass(),parent%radiusHalfMass(),relTol=1.0d-5)
  call Unit_Tests_End_Group()
  ! A truncated disk.
  call Unit_Tests_Begin_Group("Truncated disk")
  disk=massDistributionExponentialDiskTidal(truncation=1.0d0,heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
  call Assert('total mass'                               ,disk%massTotal             (     ),exp(-1.0d0)     ,relTol=1.0d-12)
  call Assert('enclosed mass converges to the total mass',disk%massEnclosedByCylinder(1.0d3),disk%massTotal(),relTol=1.0d-06)
  ! Compare the rotation curve to an independent calculation, at radii relative to the truncation radius.
  radiusTruncate=referenceRadiusTruncate(1.0d0,radiusDecayOverRadiusTruncate)
  do i=1,size(radiiTruncated)
     valuesReference(i)=referenceRotationCurveSquared(radiiTruncated(i)*radiusTruncate,radiusTruncate,radiusDecayOverRadiusTruncate,scaleHeight)
     valuesTest     (i)=disk%rotationCurve(radiiTruncated(i)*radiusTruncate)**2
  end do
  call Assert('rotation curve vs. direct Hankel transform',valuesTest,valuesReference,relTol=2.0d-3)
  ! Compare the enclosed mass to an integral of the surface density.
  do i=1,size(radiiTruncated)
     valuesReference(i)=2.0d0*Pi*disk%surfaceDensityRadialMoment(1.0d0,radiusMaximum=radiiTruncated(i)*radiusTruncate)
     valuesTest     (i)=         disk%massEnclosedByCylinder    (                    radiiTruncated(i)*radiusTruncate)
  end do
  call Assert('enclosed mass vs. integrated surface density',valuesTest,valuesReference,relTol=1.0d-3)
  ! Compare the mid-plane potential to an independent calculation.
  do i=1,size(radiiTruncated)
     valuesReference(i)=referencePotential(radiiTruncated(i)*radiusTruncate,radiusTruncate,radiusDecayOverRadiusTruncate,scaleHeight)
     valuesTest     (i)=disk%potential   (coordinateCylindrical(radiiTruncated(i)*radiusTruncate,0.0d0,0.0d0))
  end do
  call Assert('mid-plane potential vs. direct Hankel transform',valuesTest,valuesReference,relTol=2.0d-3)
  ! Check the mid-plane tidal tensor against the rotation curve, T_RR = -(dV²/dR)/R + V²/R², and Poisson's equation. (A finite
  ! difference of the acceleration is not a useful reference here: the acceleration of the untruncated, unheated disk is
  ! interpolated from tables, and differences of that interpolant are accurate only to a few percent.)
  do i=1,size(radiiTruncated)
     radius            =radiiTruncated(i)*radiusTruncate
     valuesReference(i)=-disk%rotationCurveGradient(radius)/radius+disk%rotationCurve(radius)**2/radius**2
     tidalTensor       = disk%tidalTensor (coordinateCartesian(radius,0.0d0,0.0d0))
     valuesTest     (i)= tidalTensor(1,1)
  end do
  call Assert('mid-plane tidal tensor, radial-radial, vs. rotation curve',valuesTest,valuesReference,relTol=2.0d-3)
  do i=1,size(radiiTruncated)
     radius           =radiiTruncated(i)*radiusTruncate
     accelerationPlus =disk%acceleration(coordinateCartesian(radius,0.0d0,0.0d0))
     valuesReference(i)=accelerationPlus(1)/radius
     tidalTensor      =disk%tidalTensor (coordinateCartesian(radius,0.0d0,0.0d0))
     valuesTest     (i)=tidalTensor(2,2)
  end do
  call Assert('mid-plane tidal tensor, azimuthal-azimuthal, vs. acceleration',valuesTest,valuesReference,relTol=1.0d-6)
  do i=1,size(radiiTruncated)
     radius           =radiiTruncated(i)*radiusTruncate
     tidalTensor      =disk%tidalTensor (coordinateCartesian(radius,0.0d0,0.0d0))
     valuesReference(i)=-4.0d0*Pi*disk%density(coordinateCartesian(radius,0.0d0,0.0d0))
     valuesTest     (i)=tidalTensor(1,1)+tidalTensor(2,2)+tidalTensor(3,3)
  end do
  call Assert('tidal tensor trace vs. Poisson equation',valuesTest,valuesReference,relTol=1.0d-6)
  call Unit_Tests_End_Group()
  ! A strongly truncated disk, with a truncation radius much smaller than the scale height.
  call Unit_Tests_Begin_Group("Strongly truncated disk")
  disk=massDistributionExponentialDiskTidal(truncation=6.0d0,heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
  radiusTruncate=referenceRadiusTruncate(6.0d0,radiusDecayOverRadiusTruncate)
  do i=1,size(radiiTruncated)
     valuesReference(i)=referenceRotationCurveSquared(radiiTruncated(i)*radiusTruncate,radiusTruncate,radiusDecayOverRadiusTruncate,scaleHeight)
     valuesTest     (i)=disk%rotationCurve(radiiTruncated(i)*radiusTruncate)**2
  end do
  call Assert('rotation curve vs. direct Hankel transform',valuesTest,valuesReference,relTol=2.0d-3)
  do i=1,size(radiiTruncated)
     valuesReference(i)=referencePotential(radiiTruncated(i)*radiusTruncate,radiusTruncate,radiusDecayOverRadiusTruncate,scaleHeight)
     valuesTest     (i)=disk%potential   (coordinateCylindrical(radiiTruncated(i)*radiusTruncate,0.0d0,0.0d0))
  end do
  call Assert('mid-plane potential vs. direct Hankel transform',valuesTest,valuesReference,relTol=2.0d-3)
  ! Beyond the tabulated range (at about 300 truncation radii) the multipole expansion of the potential is used. The rotation
  ! curve must be continuous across that boundary, and so V²R (which is almost constant there) must vary smoothly. Step through
  ! the boundary in small increments and find the largest change in ln(V²R) between steps. (The remaining discontinuity, of a few
  ! times 10⁻⁴, reflects the accuracy of the tabulation at its outer edge.)
  disk=massDistributionExponentialDiskTidal(truncation=8.0d0,heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
  jumpMaximum=0.0d0
  do i=1,4000
     radius     =10.0d0**(dble(i-1)/2000.0d0)
     jumpMaximum=max(                                                                                       &
          &          jumpMaximum                                                                          , &
          &          abs(log(disk%rotationCurve(radius*1.001d0)**2*1.001d0/disk%rotationCurve(radius)**2))  &
          &         )
  end do
  call Assert('rotation curve is continuous at the edge of the table',jumpMaximum < 1.0d-3,.true.)
  jumpMaximum=0.0d0
  do i=1,4000
     radius     =10.0d0**(dble(i-1)/2000.0d0)
     jumpMaximum=max(                                                                                                                                               &
          &          jumpMaximum                                                                                                                                  , &
          &          abs(log(disk%potential(coordinateCylindrical(radius*1.001d0,0.0d0,0.0d0))*1.001d0/disk%potential(coordinateCylindrical(radius,0.0d0,0.0d0))))  &
          &         )
  end do
  call Assert('potential is continuous at the edge of the table',jumpMaximum < 1.0d-3,.true.)
  ! Beyond the table the multipole expansion is used. Check its acceleration and tidal tensor against finite differences of its
  ! potential and acceleration at a point away from the mid-plane.
  position=[20.0d0,5.0d0,10.0d0]
  do i=1,3
     positionPlus              =position
     positionMinus             =position
     positionPlus         (  i)=position(i)+stepFractional*sqrt(sum(position**2))
     positionMinus        (  i)=position(i)-stepFractional*sqrt(sum(position**2))
     valuesReference      (  i)=-(                                                                                         &
          &                       +disk%potential(coordinateCartesian(positionPlus (1),positionPlus (2),positionPlus (3))) &
          &                       -disk%potential(coordinateCartesian(positionMinus(1),positionMinus(2),positionMinus(3))) &
          &                      )                                                                                         &
          &                     /2.0d0/stepFractional/sqrt(sum(position**2))
     accelerationPlus          =disk%acceleration(coordinateCartesian(positionPlus (1),positionPlus (2),positionPlus (3)))
     accelerationMinus         =disk%acceleration(coordinateCartesian(positionMinus(1),positionMinus(2),positionMinus(3)))
     tidalTensorDifference(:,i)=(accelerationPlus-accelerationMinus)/2.0d0/stepFractional/sqrt(sum(position**2))
  end do
  valuesTest =disk%acceleration(coordinateCartesian(position(1),position(2),position(3)))
  tidalTensor=disk%tidalTensor (coordinateCartesian(position(1),position(2),position(3)))
  call Assert('multipole acceleration vs. finite difference of potential'  ,valuesTest ,valuesReference      ,relTol=1.0d-6)
  call Assert('multipole tidal tensor vs. finite difference of acceleration',tidalTensor,tidalTensorDifference,relTol=1.0d-6)
  call Unit_Tests_End_Group()
  ! Away from the mid-plane. A strongly truncated disk (with R_t ≈ 0.01 R_d, and a root-mean-square radius, ≈ 0.12 R_d, set by
  ! its thickness) is close to a point mass at radii large compared to its own extent, but still within the tabulated range (which
  ! extends to about 3 R_d). There its field must match the multipole expansion of the potential to quadrupole order (rather than
  ! retain the flattening of the untruncated disk). The points include one on the axis of symmetry.
  call Unit_Tests_Begin_Group("Away from the mid-plane")
  disk      =massDistributionExponentialDiskTidal(truncation=9.0d0,heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
  quadrupole=+disk%massTotal()*Pi**2*scaleHeight**2/12.0d0 &
       &     -Pi*disk%surfaceDensityRadialMoment(3.0d0)
  do i=1,3
     angle   =anglesOffPlane(i)*Pi/180.0d0
     position=2.0d0*[sin(angle),0.0d0,cos(angle)]
     radius  =sqrt(sum(position**2))
     valuesTest     (i)=disk%potential(coordinateCartesian(position(1),position(2),position(3)))
     valuesReference(i)=-disk%massTotal()/radius-quadrupole*(3.0d0*position(3)**2-radius**2)/2.0d0/radius**5
  end do
  call Assert('potential vs. multipole expansion',valuesTest,valuesReference,relTol=2.0d-3)
  do i=1,3
     angle             =anglesOffPlane(i)*Pi/180.0d0
     position          =2.0d0*[sin(angle),0.0d0,cos(angle)]
     radius             =sqrt(sum(position**2))
     valuesTest        =disk%acceleration(coordinateCartesian(position(1),position(2),position(3)))
     valuesReference   =-disk%massTotal()*position/radius**3+quadrupole*(1.5d0*position/radius**5-7.5d0*position(3)**2*position/radius**7)
     valuesReference(3)=valuesReference(3)+quadrupole*3.0d0*position(3)/radius**5
     call Assert('acceleration vs. multipole expansion, at '//trim(labelsOffPlane(i))//' degrees from the axis',valuesTest,valuesReference,relTol=2.0d-3,absTol=1.0d-12)
     ! The tidal tensor (whose trace is set by Poisson's equation) on the axis requires the density at R=0.
     if (i == 1) then
        tidalTensor=disk%tidalTensor(coordinateCartesian(position(1),position(2),position(3)))
        call Assert('tidal tensor trace vs. Poisson equation on the axis',tidalTensor(1,1)+tidalTensor(2,2)+tidalTensor(3,3),-4.0d0*Pi*disk%density(coordinateCartesian(position(1),position(2),position(3))),absTol=1.0d-12)
     end if
  end do
  ! The tidal tensor is the symmetrized gradient of the acceleration. Check this against finite differences, both where the field
  ! is fully that of the normalized multipole expansion, and where it is blended between that and the scaled field of the
  ! untruncated disk (where the finite differences are limited by the accuracy of the tabulated acceleration of the untruncated
  ! disk). The vertical-vertical component is set by Poisson's equation, and so is excluded.
  do j=1,2
     position=radiiBlend(j)*[sin(Pi/4.0d0),0.0d0,cos(Pi/4.0d0)]
     do i=1,3
        positionPlus              =position
        positionMinus             =position
        positionPlus         (  i)=position(i)+stepFractional*sqrt(sum(position**2))
        positionMinus        (  i)=position(i)-stepFractional*sqrt(sum(position**2))
        accelerationPlus          =disk%acceleration(coordinateCartesian(positionPlus (1),positionPlus (2),positionPlus (3)))
        accelerationMinus         =disk%acceleration(coordinateCartesian(positionMinus(1),positionMinus(2),positionMinus(3)))
        tidalTensorDifference(:,i)=(accelerationPlus-accelerationMinus)/2.0d0/stepFractional/sqrt(sum(position**2))
     end do
     tidalTensorDifference     =0.5d0*(tidalTensorDifference+transpose(tidalTensorDifference))
     tidalTensor               =disk%tidalTensor(coordinateCartesian(position(1),position(2),position(3)))
     tidalTensorDifference(3,3)=tidalTensor(3,3)
     if (j == 1) then
        call Assert('tidal tensor vs. finite difference of acceleration, normalized multipole',tidalTensor,tidalTensorDifference,relTol=1.0d-4,absTol=1.0d-4*maxval(abs(tidalTensor)))
     else
        call Assert('tidal tensor vs. finite difference of acceleration, blended field'       ,tidalTensor,tidalTensorDifference,relTol=1.0d-2,absTol=1.0d-2*maxval(abs(tidalTensor)))
     end if
  end do
  ! Away from the mid-plane, the potential and acceleration must be continuous at the edge of the table (here at about 5 R_d). The
  ! tidal tensor is also evaluated along this line, as beyond the table the multipole tidal tensor must be exactly symmetric at a
  ! general position.
  disk=massDistributionExponentialDiskTidal(truncation=8.0d0,heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
  do j=1,2
     jumpMaximum=0.0d0
     do i=1,4000
        radius  =10.0d0**(dble(i-1)/2000.0d0)
        position=radius*[sin(Pi/6.0d0),0.0d0,cos(Pi/6.0d0)]
        if (j == 1) then
           potentialCurrent   =disk%potential(coordinateCartesian(position(1),position(2),position(3)))*radius
        else
           valuesTest         =disk%acceleration(coordinateCartesian(position(1),position(2),position(3)))
           accelerationCurrent=sqrt(sum(valuesTest**2))*radius**2
           tidalTensor        =disk%tidalTensor (coordinateCartesian(position(1),position(2),position(3)))
        end if
        if (i > 1) then
           if (j == 1) then
              jumpMaximum=max(jumpMaximum,abs(log(potentialCurrent   /potentialPrevious   )))
           else
              jumpMaximum=max(jumpMaximum,abs(log(accelerationCurrent/accelerationPrevious)))
           end if
        end if
        potentialPrevious   =potentialCurrent
        accelerationPrevious=accelerationCurrent
     end do
     if (j == 1) then
        call Assert('potential is continuous at the edge of the table, away from the mid-plane'   ,jumpMaximum < 2.0d-3,.true.)
     else
        call Assert('acceleration is continuous at the edge of the table, away from the mid-plane',jumpMaximum < 2.0d-3,.true.)
     end if
  end do
  call Unit_Tests_End_Group()
  ! Away from the mid-plane, inside the truncation radius. The scaled field of the untruncated, unheated disk overestimates the
  ! vertical acceleration of a strongly truncated disk by factors of two to three here, so this checks the tabulated corrections.
  ! The reference is the exact field, from direct evaluation of the Hankel transforms.
  call Unit_Tests_Begin_Group("Inside the truncation radius, away from the mid-plane")
  do k=1,size(truncationsInner)
     disk          =massDistributionExponentialDiskTidal(truncation=truncationsInner(k),heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
     radiusTruncate=referenceRadiusTruncate(truncationsInner(k),radiusDecayOverRadiusTruncate)
     l=0
     do i=1,size(radiiInner)
        do j=1,size(heightsInner)
           l       =l+1
           position=[radiiInner(i)*radiusTruncate,0.0d0,heightsInner(j)*scaleHeight]
           call referenceFieldOffPlane(position(1),position(3),radiusTruncate,radiusDecayOverRadiusTruncate,scaleHeight,potentialInnerReference(l),radialInnerReference(l),verticalInnerReference(l))
           valuesTest       =disk%acceleration(coordinateCartesian(position(1),position(2),position(3)))
           potentialInner(l)=disk%potential   (coordinateCartesian(position(1),position(2),position(3)))
           radialInner   (l)=valuesTest(1)
           verticalInner (l)=valuesTest(3)
        end do
     end do
     write (labelInner,'(a,f3.1)') 'x = ',truncationsInner(k)
     ! The measured agreement is 1.5×10⁻³ for the potential, 8×10⁻³ for the radial acceleration, and 1.2×10⁻² for the vertical
     ! acceleration (limited by interpolation in the coarse grids on which the corrections are tabulated). Without the
     ! corrections the vertical acceleration differs by up to a factor of about three for x=3.
     call Assert('potential vs. direct Hankel transform, '            //trim(labelInner),potentialInner,potentialInnerReference,relTol=3.0d-3)
     call Assert('radial acceleration vs. direct Hankel transform, '  //trim(labelInner),radialInner   ,radialInnerReference   ,relTol=1.5d-2)
     call Assert('vertical acceleration vs. direct Hankel transform, '//trim(labelInner),verticalInner ,verticalInnerReference ,relTol=2.5d-2)
     ! The tidal tensor must be the gradient of the corrected acceleration.
     position   =[0.6d0*radiusTruncate,0.0d0,2.0d0*scaleHeight]
     tidalTensor=disk%tidalTensor(coordinateCartesian(position(1),position(2),position(3)))
     do i=1,3
        positionPlus     =position
        positionMinus    =position
        positionPlus (i) =position(i)+stepFractional*radiusTruncate
        positionMinus(i) =position(i)-stepFractional*radiusTruncate
        accelerationPlus =disk%acceleration(coordinateCartesian(positionPlus (1),positionPlus (2),positionPlus (3)))
        accelerationMinus=disk%acceleration(coordinateCartesian(positionMinus(1),positionMinus(2),positionMinus(3)))
        tidalTensorDifference(:,i)=(accelerationPlus-accelerationMinus)/2.0d0/stepFractional/radiusTruncate
     end do
     ! Only the off-diagonal and radial-radial components are compared: the vertical-vertical component is set by Poisson's
     ! equation, and the azimuthal-azimuthal component is a_R/R. The approximate field is not exactly curl-free, so the
     ! off-diagonal component is compared with the mean of the two finite differences. The measured agreement is 3×10⁻³ for the
     ! off-diagonal component, and 1.1×10⁻² for the radial-radial component (which is the same, 1.1×10⁻², without the
     ! corrections, arising from the tabulation of the field of the untruncated disk).
     call Assert('tidal tensor vs. finite difference of acceleration, '//trim(labelInner),[tidalTensor(1,1),tidalTensor(1,3)],[tidalTensorDifference(1,1),0.5d0*(tidalTensorDifference(1,3)+tidalTensorDifference(3,1))],relTol=2.0d-2)
     ! At an azimuth of 45°, where both horizontal components of the acceleration are non-zero, such that the gradients of the
     ! corrections enter every component of the tidal tensor. Here T_xy=(T_RR-T_φφ)/2 is a small difference (6% of T_xx for x=3),
     ! so the error of about 1% in T_RR (see above) gives it a relative error of up to 5%. Components are therefore also accepted
     ! if they agree to within 1% of the largest component of the tensor. The measured agreement is 0.4% of the largest component,
     ! while applying the gradients of the correction to the vertical acceleration to the y component gave errors of 3-7%.
     position   =[0.6d0*radiusTruncate/sqrt(2.0d0),0.6d0*radiusTruncate/sqrt(2.0d0),2.0d0*scaleHeight]
     tidalTensor=disk%tidalTensor(coordinateCartesian(position(1),position(2),position(3)))
     do i=1,3
        positionPlus     =position
        positionMinus    =position
        positionPlus (i) =position(i)+stepFractional*radiusTruncate
        positionMinus(i) =position(i)-stepFractional*radiusTruncate
        accelerationPlus =disk%acceleration(coordinateCartesian(positionPlus (1),positionPlus (2),positionPlus (3)))
        accelerationMinus=disk%acceleration(coordinateCartesian(positionMinus(1),positionMinus(2),positionMinus(3)))
        tidalTensorDifference(:,i)=(accelerationPlus-accelerationMinus)/2.0d0/stepFractional/radiusTruncate
     end do
     call Assert(                                                                                                                                 &
          &      'tidal tensor vs. finite difference of acceleration, at 45 degrees azimuth, '//trim(labelInner)                                , &
          &      [tidalTensor          (1,1),      tidalTensor          (1,2)                             ,      tidalTensor          (1,3)                              , &
          &       tidalTensor          (2,2),      tidalTensor          (2,3)                                                                                          ], &
          &      [tidalTensorDifference(1,1),0.5d0*(tidalTensorDifference(1,2)+tidalTensorDifference(2,1)),0.5d0*(tidalTensorDifference(1,3)+tidalTensorDifference(3,1)), &
          &       tidalTensorDifference(2,2),0.5d0*(tidalTensorDifference(2,3)+tidalTensorDifference(3,2))                                                         ], &
          &      absTol=1.0d-2*maxval(abs(tidalTensorDifference))                                                                                                  , &
          &      relTol=2.0d-2                                                                                                                                       &
          &     )
  end do
  call Unit_Tests_End_Group()
  ! A truncation stronger than that of the smallest tabulated truncation radius.
  call Unit_Tests_Begin_Group("Truncation beyond the tabulated range")
  disk=massDistributionExponentialDiskTidal(truncation=12.0d0,heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
  call Assert('total mass'                                ,disk%massTotal()                  ,exp(-12.0d0)    ,relTol=1.0d-12)
  call Assert('enclosed mass converges to the total mass' ,disk%massEnclosedByCylinder(1.0d3),disk%massTotal(),relTol=1.0d-06)
  massMaximum=0.0d0
  do i=1,4000
     massMaximum=max(massMaximum,disk%massEnclosedByCylinder(10.0d0**(-3.0d0+dble(i-1)/1000.0d0)))
  end do
  call Assert('enclosed mass never exceeds the total mass',     massMaximum                                  ,      disk%massTotal(),relTol=1.0d-6)
  call Assert('half-mass radius encloses half of the mass',disk%massEnclosedByCylinder(disk%radiusHalfMass()),0.5d0*disk%massTotal(),relTol=1.0d-4)
  call Unit_Tests_End_Group()
  ! A heated disk.
  call Unit_Tests_Begin_Group("Heated disk")
  disk=massDistributionExponentialDiskTidal(truncation=1.0d-8,heating=0.1d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
  call Assert('mass is conserved'       ,         disk%massEnclosedByCylinder    (1.0d3)                 ,                                                  parent%massTotal(),relTol=1.0d-5)
  call Assert('mean radius is conserved',2.0d0*Pi*disk%surfaceDensityRadialMoment(2.0d0)/disk%massTotal(),2.0d0*Pi*parent%surfaceDensityRadialMoment(2.0d0)/parent%massTotal(),relTol=1.0d-3)
  ! Heating spreads the disk, moving mass to large radii. (It need not increase the half-mass radius: the distribution of radii
  ! into which each ring is spread conserves the mean radius, but is skewed, such that its median is reduced.)
  call Assert('heating increases the mass beyond 5 scale lengths',disk%massTotal()-disk%massEnclosedByCylinder(5.0d0) > parent%massTotal()-parent%massEnclosedByCylinder(5.0d0),.true.)
  call Unit_Tests_End_Group()
  ! Consistency of the rotation curve gradient.
  call Unit_Tests_Begin_Group("Truncated and heated disk")
  disk=massDistributionExponentialDiskTidal(truncation=1.0d0,heating=0.1d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
  do i=1,size(radii)
     valuesTidal (i)=disk%rotationCurveGradient(radii(i))
     valuesParent(i)=+(                                                        &
          &            +disk%rotationCurve(radii(i)*(1.0d0+stepFractional))**2 &
          &            -disk%rotationCurve(radii(i)*(1.0d0-stepFractional))**2 &
          &           )                                                        &
          &          /2.0d0/stepFractional/radii(i)
  end do
  call Assert('rotation curve gradient vs. finite difference',valuesTidal,valuesParent,relTol=1.0d-4)
  call Unit_Tests_End_Group()
  ! Away from the mid-plane, inside the truncation radius, of truncated and heated disks. This checks the interpolation of the
  ! tabulated corrections in heating parameter. The reference values are from exponentialDiskTidalHeated.py (in
  ! galacticusDevTools/physicsReferences), which evaluates the Hankel transform of the heated disk as that of a superposition of
  ! rings spread in radius, and shares none of the tabulation used here. Points are at R/R_t = 0.3, 0.6, and 1.0 (outer), and
  ! z/h = 0.5, 1, 2, and 4 (inner).
  potentialHeatedReference(:,1)=[-6.07182219d-01,-5.81497759d-01,-5.15548660d-01,-4.03418393d-01,-5.50361311d-01,-5.29387953d-01, &
                            &   -4.74764306d-01,-3.79021270d-01,-4.49682152d-01,-4.35392483d-01,-3.97980018d-01,-3.30433301d-01]
  radialHeatedReference   (:,1)=[-1.70869072d-01,-1.53774609d-01,-1.15061956d-01,-6.54096624d-02,-2.50462481d-01,-2.32540744d-01, &
                            &   -1.87414720d-01,-1.16075520d-01,-3.18539208d-01,-2.96459610d-01,-2.41007352d-01,-1.53367584d-01]
  verticalHeatedReference (:,1)=[-2.82876380d-01,-4.44028342d-01,-4.79449702d-01,-3.38624535d-01,-2.30153474d-01,-3.63941392d-01, &
                            &   -4.01471893d-01,-2.94778224d-01,-1.56729649d-01,-2.48145432d-01,-2.76803181d-01,-2.14219896d-01]
  potentialHeatedReference(:,2)=[-6.31724356d-01,-6.03832430d-01,-5.32149525d-01,-4.11197937d-01,-5.72163171d-01,-5.48040663d-01, &
                            &   -4.86688578d-01,-3.83634346d-01,-4.41316678d-01,-4.28464994d-01,-3.93965515d-01,-3.29080273d-01]
  radialHeatedReference   (:,2)=[-1.56884463d-01,-1.46310469d-01,-1.18955842d-01,-7.28043921d-02,-3.07049828d-01,-2.83289075d-01, &
                            &   -2.23790931d-01,-1.32343820d-01,-3.80229656d-01,-3.49565369d-01,-2.75334770d-01,-1.66879461d-01]
  verticalHeatedReference (:,2)=[-3.06971014d-01,-4.82481420d-01,-5.20744707d-01,-3.61492931d-01,-2.66432144d-01,-4.15848874d-01, &
                            &   -4.43323911d-01,-3.10319963d-01,-1.39969078d-01,-2.24749902d-01,-2.59637112d-01,-2.09605943d-01]
  call Unit_Tests_Begin_Group("Truncated and heated disk, away from the mid-plane")
  radiusTruncate=referenceRadiusTruncate(1.0d0,radiusDecayOverRadiusTruncate)
  do k=1,size(heatingsInner)
     disk=massDistributionExponentialDiskTidal(truncation=1.0d0,heating=heatingsInner(k),radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
     l=0
     do i=1,size(radiiInner)
        do j=1,size(heightsInner)
           l                =l+1
           position         =[radiiInner(i)*radiusTruncate,0.0d0,heightsInner(j)*scaleHeight]
           valuesTest       =disk%acceleration(coordinateCartesian(position(1),position(2),position(3)))
           potentialInner(l)=disk%potential   (coordinateCartesian(position(1),position(2),position(3)))
           radialInner   (l)=valuesTest(1)
           verticalInner (l)=valuesTest(3)
        end do
     end do
     write (labelInner,'(a,f3.1)') 'x = 1.0, s = ',heatingsInner(k)
     ! The measured agreement is 1.6×10⁻³ for the potential, 9.6×10⁻³ for the radial acceleration, and 1.5×10⁻² for the vertical
     ! acceleration (limited by interpolation in the coarse grids on which the corrections are tabulated). The field of these disks
     ! differs from that of the unheated disk by up to 26%.
     call Assert('potential vs. superposition of rings, '            //trim(labelInner),potentialInner,potentialHeatedReference(:,k),relTol=3.0d-3)
     call Assert('radial acceleration vs. superposition of rings, '  //trim(labelInner),radialInner   ,radialHeatedReference   (:,k),relTol=1.5d-2)
     call Assert('vertical acceleration vs. superposition of rings, '//trim(labelInner),verticalInner ,verticalHeatedReference (:,k),relTol=2.5d-2)
  end do
  call Unit_Tests_End_Group()
  ! The specific angular momentum of material removed by further truncation, j/⟨j⟩ = 1 - d log⟨R⟩/dx.
  call Unit_Tests_Begin_Group("Specific angular momentum of stripped material")
  truncation=galacticComponentTruncationExponential(radiusDecayOverRadiusTruncate)
  ! Compare with values computed independently (using incomplete gamma functions in Python's scipy.special) for truncation radii
  ! of 1 and 0.1 scale lengths.
  call Assert('ratio vs. reference, R_t = 1.0',truncation%ratioAngularMomentumCylindrical(exp(-8.520917d-1)),1.710221d0,relTol=1.0d-5)
  call Assert('ratio vs. reference, R_t = 0.1',truncation%ratioAngularMomentumCylindrical(exp(-4.602535d+0)),1.517677d0,relTol=1.0d-5)
  ! Compare with the change in the mean radius of the tabulated truncated disk between successive degrees of truncation. Since
  ! j/⟨j⟩ = 1 - d log⟨R⟩/dx, the integral of j/⟨j⟩-1 over x between two truncations should equal the change in log⟨R⟩. (The
  ! derivative of the mean radius of the tabulated disk is not compared directly, as it has small discontinuities at the
  ! tabulated truncation radii.) The integral is evaluated using Simpson's rule in log x.
  do i=1,size(truncations)-1
     diskPlus       =massDistributionExponentialDiskTidal(truncation=truncations(i+1),heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
     diskMinus      =massDistributionExponentialDiskTidal(truncation=truncations(i  ),heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
     valuesParent(i)=+log(diskMinus%surfaceDensityRadialMoment(2.0d0)/diskMinus%surfaceDensityRadialMoment(1.0d0)) &
          &          -log(diskPlus %surfaceDensityRadialMoment(2.0d0)/diskPlus %surfaceDensityRadialMoment(1.0d0))
     do j=0,countIntegration
        logarithmTruncationIntegration(j)=+log(truncations(i))+log(truncations(i+1)/truncations(i))*dble(j)/dble(countIntegration)
        integrandIntegration          (j)=+exp(logarithmTruncationIntegration(j))                                                           &
             &                            *(truncation%ratioAngularMomentumCylindrical(exp(-exp(logarithmTruncationIntegration(j))))-1.0d0)
     end do
     valuesTidal(i)=integrateSimpson(logarithmTruncationIntegration,integrandIntegration)
  end do
  call Assert('integrated ratio vs. change in mean radius of the truncated disk',valuesTidal(1:size(truncations)-1),valuesParent(1:size(truncations)-1),relTol=1.0d-3)
  ! The shape of the disk is unchanged for negligible truncations, and for truncations beyond the tabulated range.
  call Assert('ratio for negligible truncation' ,truncation%ratioAngularMomentumCylindrical(exp(-1.0d-5)),1.0d0,relTol=1.0d-12)
  call Assert('ratio beyond the tabulated range',truncation%ratioAngularMomentumCylindrical(exp(-12.0d0)),1.0d0,relTol=1.0d-12)
  ! The shape of the disk is held fixed for truncations beyond that of the smallest tabulated truncation radius, 0.01 scale lengths,
  ! for which x = 9.114112 (computed using Python's scipy.special). Check that the ratio becomes unity at this point, and not before.
  diskPlus =massDistributionExponentialDiskTidal(truncation=9.10d0,heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
  diskMinus=massDistributionExponentialDiskTidal(truncation=9.05d0,heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
  call Assert('disk shape changes below the limit',diskPlus%surfaceDensityRadialMoment(2.0d0)/diskPlus%surfaceDensityRadialMoment(1.0d0) < 0.9999d0*diskMinus%surfaceDensityRadialMoment(2.0d0)/diskMinus%surfaceDensityRadialMoment(1.0d0),.true.)
  diskPlus =massDistributionExponentialDiskTidal(truncation=9.20d0,heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
  diskMinus=massDistributionExponentialDiskTidal(truncation=9.13d0,heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
  call Assert('disk shape is fixed beyond the limit',diskPlus%surfaceDensityRadialMoment(2.0d0)/diskPlus%surfaceDensityRadialMoment(1.0d0),diskMinus%surfaceDensityRadialMoment(2.0d0)/diskMinus%surfaceDensityRadialMoment(1.0d0),relTol=1.0d-12)
  call Assert('ratio is not unity below the limit',truncation%ratioAngularMomentumCylindrical(exp(-9.10d0)) > 1.4d0,.true.               )
  call Assert('ratio is unity beyond the limit'   ,truncation%ratioAngularMomentumCylindrical(exp(-9.13d0))         ,1.0d0,relTol=1.0d-12)
  call Unit_Tests_End_Group()
  ! The mean radius of the truncated disk, and the relation between truncation radius and the fraction of mass retained.
  call Unit_Tests_Begin_Group("Mean radius and truncation radius")
  ! Compare with values computed independently (using incomplete gamma functions and numerical integration in Python's scipy) for
  ! truncation radii of 3, 1, 0.1, and 0.02 scale lengths.
  call Assert('mean radius vs. reference, R_t = 3.0' ,truncation%ratioRadiusMeanCylindrical     (exp(-8.2510504329d-2)),8.6293985693d-1,relTol=1.0d-5)
  call Assert('mean radius vs. reference, R_t = 1.0' ,truncation%ratioRadiusMeanCylindrical     (exp(-8.5209171170d-1)),4.3500690895d-1,relTol=1.0d-5)
  call Assert('mean radius vs. reference, R_t = 0.1' ,truncation%ratioRadiusMeanCylindrical     (exp(-4.6025345784d+0)),5.1576214022d-2,relTol=1.0d-5)
  call Assert('mean radius vs. reference, R_t = 0.02',truncation%ratioRadiusMeanCylindrical     (exp(-7.7382916764d+0)),1.0464576998d-2,relTol=1.0d-5)
  call Assert('fraction retained, R_t = 3.0'         ,truncation%fractionMassRetainedCylindrical(     3.0000000000d+0 ),9.2080176540d-1,relTol=1.0d-5)
  call Assert('fraction retained, R_t = 1.0'         ,truncation%fractionMassRetainedCylindrical(     1.0000000000d+0 ),4.2652183751d-1,relTol=1.0d-5)
  call Assert('fraction retained, R_t = 0.1'         ,truncation%fractionMassRetainedCylindrical(     1.0000000000d-1 ),1.0026390838d-2,relTol=1.0d-5)
  call Assert('fraction retained, R_t = 0.02'        ,truncation%fractionMassRetainedCylindrical(     2.0000000000d-2 ),4.3581545332d-4,relTol=1.0d-5)
  ! Compare the mean radius with that of the tabulated truncated disk.
  do i=1,size(truncations)
     diskPlus       =massDistributionExponentialDiskTidal(truncation=truncations(i),heating=0.0d0,radiusDecayOverRadiusTruncate=radiusDecayOverRadiusTruncate,massDistribution_=parent)
     valuesParent(i)=+diskPlus%surfaceDensityRadialMoment(2.0d0)/diskPlus%surfaceDensityRadialMoment(1.0d0) &
          &          /(parent %surfaceDensityRadialMoment(2.0d0)/parent  %surfaceDensityRadialMoment(1.0d0))
     valuesTidal (i)=truncation%ratioRadiusMeanCylindrical(exp(-truncations(i)))
  end do
  call Assert('mean radius vs. tabulated truncated disk',valuesTidal(1:size(truncations)),valuesParent(1:size(truncations)),relTol=1.0d-3)
  ! Limits: negligible truncation, truncation radii beyond which truncation is negligible, and within the smallest modeled.
  call Assert('mean radius for negligible truncation'     ,truncation%ratioRadiusMeanCylindrical     (exp( -1.0d-5)),1.0d0                                              ,relTol=1.0d-12)
  call Assert('mean radius fixed beyond the limit'        ,truncation%ratioRadiusMeanCylindrical     (exp(-12.0d+0)),truncation%ratioRadiusMeanCylindrical(exp(-9.13d0)),relTol=1.0d-12)
  call Assert('fraction retained for large R_t'           ,truncation%fractionMassRetainedCylindrical(      1.0d+2 ),1.0d0                                              ,relTol=1.0d-12)
  call Assert('fraction retained within the smallest R_t' ,truncation%fractionMassRetainedCylindrical(      5.0d-3 ),0.0d0                                              ,absTol=1.0d-30)
  call Unit_Tests_End_Group()
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_Mass_Distributions_Exponential_Disk_Tidal
