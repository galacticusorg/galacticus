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
  use :: Coordinates                          , only : coordinateCylindrical          , coordinateCartesian                 , assignment(=)
  use :: Display                              , only : displayVerbositySet            , verbosityLevelStandard
  use :: Events_Hooks                         , only : eventsHooksInitialize
  use :: IO_HDF5                              , only : ioHDF5AccessInitialize
  use :: Mass_Distributions                   , only : massDistributionExponentialDisk, massDistributionExponentialDiskTidal
  use :: Numerical_Constants_Math             , only : Pi
  use :: Tensors                              , only : assignment(=)
  use :: Test_Exponential_Disk_Tidal_Reference, only : referenceRotationCurveSquared  , referencePotential                  , referenceRadiusTruncate
  use :: Unit_Tests                           , only : Assert                         , Unit_Tests_Begin_Group              , Unit_Tests_End_Group   , Unit_Tests_Finish
  implicit none
  double precision                                      , parameter      :: scaleHeight      =0.137d0                        , radiusDecayOverRadiusTruncate=0.1d0, &
       &                                                                    stepFractional   =1.0d-4
  double precision                                      , dimension(5  ) :: radii            =[0.1d0,0.5d0,1.0d0,2.0d0,5.0d0]
  double precision                                      , dimension(3  ) :: radiiTruncated   =[0.5d0,1.0d0,2.0d0]
  double precision                                      , dimension(3  ) :: anglesOffPlane   =[0.0d0,45.0d0,80.0d0]
  double precision                                      , dimension(2  ) :: radiiBlend       =[1.5d0,0.4d0]
  character       (len=2)                               , dimension(3  ) :: labelsOffPlane   =['0 ','45','80']
  double precision                                      , dimension(5  ) :: valuesTidal                                      , valuesParent
  double precision                                      , dimension(3  ) :: valuesReference                                  , valuesTest                         , &
       &                                                                    accelerationPlus                                 , accelerationMinus
  double precision                                      , dimension(3,3) :: tidalTensor                                      , tidalTensorDifference
  double precision                                      , dimension(3  ) :: position                                         , positionPlus                       , &
       &                                                                    positionMinus
  type            (massDistributionExponentialDisk     )                 :: parent
  type            (massDistributionExponentialDiskTidal)                 :: disk
  double precision                                                       :: radiusTruncate                                   , radius                             , &
       &                                                                    jumpMaximum                                      , massMaximum                        , &
       &                                                                    quadrupole                                       , angle                              , &
       &                                                                    potentialPrevious                                , accelerationPrevious               , &
       &                                                                    potentialCurrent                                 , accelerationCurrent
  integer                                                                :: i                                                , j

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
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_Mass_Distributions_Exponential_Disk_Tidal
