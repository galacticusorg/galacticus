!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the exponential disk node component.

module Node_Component_Disk_Exponential_Data
  use Kind_Numbers
  use Tables
  implicit none
  public

  ! Record of unique ID of node which we last computed results for.
  integer(kind=kind_int8)                     :: lastUniqueID=-1
  !$omp threadprivate(lastUniqueID)
  
  ! Records of previously computed and stored quantities.
  logical                                     :: surfaceDensityCentralGasComputed,surfaceDensityCentralStellarComputed&
       &,surfaceDensityCentralTotalComputed 
  !$omp threadprivate(surfaceDensityCentralGasComputed,surfaceDensityCentralStellarComputed,surfaceDensityCentralTotalComputed)
  double precision                            :: surfaceDensityCentralGas,surfaceDensityCentralStellar,surfaceDensityCentralTotal
  !$omp threadprivate(surfaceDensityCentralGas,surfaceDensityCentralStellar,surfaceDensityCentralTotal)
  logical                                     :: radiusScaleDiskComputed
  !$omp threadprivate(radiusScaleDiskComputed)
  double precision                            :: radiusScaleDisk
  !$omp threadprivate(radiusScaleDisk)

  ! Luminosity work arrays.
  double precision, allocatable, dimension(:) :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)

  ! Tabulation of the exponential disk rotation curve.
  integer,          parameter                 :: rotationCurvePointsPerDecade=10
  integer                                     :: rotationCurvePointsCount
  logical                                     :: rotationCurveInitialized=.false.
  double precision, parameter                 :: rotationCurveHalfRadiusMinimumDefault=1.0d-6,rotationCurveHalfRadiusMaximumDefault=10.0d0
  double precision                            :: rotationCurveHalfRadiusMinimum=rotationCurveHalfRadiusMinimumDefault&
       &,rotationCurveHalfRadiusMaximum=rotationCurveHalfRadiusMaximumDefault
  double precision                            :: scaleLengthFactor,diskStructureSolverSpecificAngularMomentum,diskRadiusSolverFlatVsSphericalFactor
  logical                                     :: scaleLengthFactorSet=.false.
  type(table1DLogarithmicLinear)              :: rotationCurveTable

  ! Tabulation of the exponential disk rotation curve gradient.
  integer,          parameter                 :: rotationCurveGradientPointsPerDecade=10
  integer                                     :: rotationCurveGradientPointsCount
  logical                                     :: rotationCurveGradientInitialized=.false.
  double precision, parameter                 :: rotationCurveGradientHalfRadiusMinimumDefault=1.0d-6,rotationCurveGradientHalfRadiusMaximumDefault=10.0d0
  double precision                            :: rotationCurveGradientHalfRadiusMinimum=rotationCurveGradientHalfRadiusMinimumDefault&
       &,rotationCurveGradientHalfRadiusMaximum=rotationCurveGradientHalfRadiusMaximumDefault
  type(table1DLogarithmicLinear)              :: rotationCurveGradientTable

  ! The radius (in units of the disk scale length) beyond which the disk is treated as a point mass for the purposes of computing
  ! rotation curves.
  double precision, parameter                 :: fractionalRadiusMaximum=30.0d0

  ! Parameter controlling the scale height of the disk.
  double precision                            :: heightToRadialScaleDisk
contains

  subroutine Node_Component_Disk_Exponential_Reset(uniqueID)
    !% Reset calculations for the exponential disk component.
    use Kind_Numbers
    implicit none
    integer(kind=kind_int8), intent(in   ) :: uniqueID

    radiusScaleDiskComputed             =.false.
    surfaceDensityCentralGasComputed    =.false.
    surfaceDensityCentralStellarComputed=.false.
    surfaceDensityCentralTotalComputed  =.false.
    lastUniqueID                        =uniqueID
    return
  end subroutine Node_Component_Disk_Exponential_Reset
  
  double precision function Node_Component_Disk_Exponential_Enclosed_Mass_Dimensionless(radius)
    !% Returns the fractional mass enclosed within {\tt radius} in a dimensionless exponential disk.
    implicit none
    double precision, intent(in) :: radius
    
    Node_Component_Disk_Exponential_Enclosed_Mass_Dimensionless=1.0d0-(1.0d0+radius)*exp(-radius)
    return
  end function Node_Component_Disk_Exponential_Enclosed_Mass_Dimensionless
   
   double precision function Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors(halfRadius)
     !% Compute Bessel function factors appearing in the expression for an razor-thin exponential disk rotation curve.
     use Memory_Management
     use Numerical_Constants_Math
     use Numerical_Interpolation
     use Bessel_Functions
     implicit none
     double precision, intent(in) :: halfRadius
     double precision, parameter  :: halfRadiusSmall=1.0d-3
     integer                      :: iPoint
     double precision             :: x


     ! For small half-radii, use a series expansion for a more accurate result.
     if (halfRadius <= 0.0d0) then
        Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors=0.0d0
        return
     else if (halfRadius < halfRadiusSmall) then
        Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors=(ln2-eulersConstant-0.5d0-log(halfRadius))*halfRadius**2
        return
     end if
     
     !$omp critical(Exponential_Disk_Rotation_Curve_Tabulate)
     if (.not.rotationCurveInitialized .or. halfRadius <  rotationCurveHalfRadiusMinimum .or. halfRadius > rotationCurveHalfRadiusMaximum) then
        ! Find the minimum and maximum half-radii to tabulate.
        rotationCurveHalfRadiusMinimum=min(rotationCurveHalfRadiusMinimum,0.5d0*halfRadius)
        rotationCurveHalfRadiusMaximum=max(rotationCurveHalfRadiusMaximum,2.0d0*halfRadius)
        
        ! Determine how many points to tabulate.
        rotationCurvePointsCount=int(dlog10(rotationCurveHalfRadiusMaximum/rotationCurveHalfRadiusMinimum)*dble(rotationCurvePointsPerDecade))+1
        
        ! Allocate table arrays.
        call rotationCurveTable%destroy()
        call rotationCurveTable%create(rotationCurveHalfRadiusMinimum,rotationCurveHalfRadiusMaximum,rotationCurvePointsCount)
        
        ! Compute Bessel factors.
        do iPoint=1,rotationCurvePointsCount
           x=rotationCurveTable%x(iPoint)
           call rotationCurveTable%populate(                                               &
                &                            x**2                                          &
                &                           *(                                             &
                &                              Bessel_Function_I0(x)*Bessel_Function_K0(x) &
                &                             -Bessel_Function_I1(x)*Bessel_Function_K1(x) &
                &                            ),                                            &
                &                           iPoint                                         &
                &                          )
        end do
        
        ! Flag that the rotation curve is now initialized.
        rotationCurveInitialized=.true.
     end if
     
     ! Interpolate in the tabulated function.
     Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors=rotationCurveTable%interpolate(halfRadius)
     !$omp end critical(Exponential_Disk_Rotation_Curve_Tabulate)
     
     return
   end function Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors

  double precision function Node_Component_Disk_Exponential_Rttn_Crv_Grdnt_Bssl_Fctrs(halfRadius)
    !% Compute Bessel function factors appearing in the expression for a razor-thin exponential disk rotation curve gradient.
    use Memory_Management
    use Numerical_Constants_Math
    use Numerical_Interpolation
    use Bessel_Functions
    implicit none
    double precision, intent(in) :: halfRadius
    double precision, parameter  :: halfRadiusSmall=1.0d-3
    double precision, parameter  :: halfRadiusLarge=1.0d+2
    integer                      :: iPoint
    double precision             :: x

    ! For small and large half-radii, use a series expansion for a more accurate result.
    if (halfRadius == 0.0d0) then
       Node_Component_Disk_Exponential_Rttn_Crv_Grdnt_Bssl_Fctrs=0.0d0
       return
    else if (halfRadius < halfRadiusSmall) then
       Node_Component_Disk_Exponential_Rttn_Crv_Grdnt_Bssl_Fctrs=(ln2-log(halfRadius)-eulersConstant-1.0d0)*halfRadius &
            & +(1.5d0*ln2-1.5d0*log(halfRadius)+0.25d0-1.5d0*eulersConstant)*halfRadius**3
       return
    else if (halfRadius > halfRadiusLarge) then
       Node_Component_Disk_Exponential_Rttn_Crv_Grdnt_Bssl_Fctrs=-0.125d0-27.0d0/64.0d0/halfRadius**2
       return
    end if

    !$omp critical(Exponential_Disk_Rotation_Curve_Gradient_Tabulate)
    if (.not.rotationCurveGradientInitialized .or. halfRadius <  rotationCurveHalfRadiusMinimum .or. halfRadius > rotationCurveGradientHalfRadiusMaximum) then
       ! Find the minimum and maximum half-radii to tabulate.
       rotationCurveGradientHalfRadiusMinimum=min(rotationCurveGradientHalfRadiusMinimum,0.5d0*halfRadius)
       rotationCurveGradientHalfRadiusMaximum=max(rotationCurveGradientHalfRadiusMaximum,2.0d0*halfRadius)

       ! Determine how many points to tabulate.
       rotationCurveGradientPointsCount=int(dlog10(rotationCurveGradientHalfRadiusMaximum/rotationCurveGradientHalfRadiusMinimum)*dble(rotationCurveGradientPointsPerDecade))+1

       ! Allocate table arrays.
       call rotationCurveGradientTable%destroy()
       call rotationCurveGradientTable%create(rotationCurveGradientHalfRadiusMinimum,rotationCurveGradientHalfRadiusMaximum,rotationCurveGradientPointsCount)

       ! Compute Bessel factors.
       do iPoint=1,rotationCurveGradientPointsCount
          x=rotationCurveGradientTable%x(iPoint)
          call rotationCurveGradientTable%populate(                                                 &
               &                                    x**2                                            &
               &                                   *( x                                             &
               &                                     *  Bessel_Function_I0(x)*Bessel_Function_K0(x) &
               &                                     +x**2                                          &
               &                                     *( Bessel_Function_I1(x)*Bessel_Function_K0(x) &
               &                                       -Bessel_Function_I0(x)*Bessel_Function_K1(x) &
               &                                      )                                             &
               &                                    ),                                              &
               &                                   iPoint                                           &
               &                                  )
       end do
       
       ! Flag that the rotation curve is now initialized.
       rotationCurveGradientInitialized=.true.
    end if
    
    ! Interpolate in the tabulated function.
    Node_Component_Disk_Exponential_Rttn_Crv_Grdnt_Bssl_Fctrs=rotationCurveGradientTable%interpolate(halfRadius)
    !$omp end critical(Exponential_Disk_Rotation_Curve_Gradient_Tabulate)
    return
  end function Node_Component_Disk_Exponential_Rttn_Crv_Grdnt_Bssl_Fctrs

  !# <galacticusStateStoreTask>
  !#  <unitName>Node_Component_Disk_Exponential_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Node_Component_Disk_Exponential_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) rotationCurveHalfRadiusMinimum,rotationCurveHalfRadiusMaximum,scaleLengthFactor,scaleLengthFactorSet,&
                    & rotationCurveGradientHalfRadiusMinimum,rotationCurveGradientHalfRadiusMaximum
    return
  end subroutine Node_Component_Disk_Exponential_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Node_Component_Disk_Exponential_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Node_Component_Disk_Exponential_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) rotationCurveHalfRadiusMinimum,rotationCurveHalfRadiusMaximum,scaleLengthFactor,scaleLengthFactorSet,&
                   & rotationCurveGradientHalfRadiusMinimum,rotationCurveGradientHalfRadiusMaximum
    ! Flag that the table is now uninitialized.
    rotationCurveInitialized=.false.
    return
  end subroutine Node_Component_Disk_Exponential_State_Retrieve
  
end module Node_Component_Disk_Exponential_Data
