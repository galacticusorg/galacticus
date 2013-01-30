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

! S\'ersic functions.
subroutine Mass_Distribution_Sersic_Initialize(self,index,halfMassRadius,densityNormalization,mass,isDimensionless)
  !% Initialize a S\'ersic mass distribution.
  use Numerical_Constants_Math
  use Numerical_Comparison
  use Galacticus_Error
  implicit none
  class           (massDistributionSersic), intent(  out)           :: self
  double precision                        , intent(in   )           :: index
  double precision                        , intent(in   ), optional :: halfMassRadius,densityNormalization,mass
  logical                                 , intent(in   ), optional :: isDimensionless

  ! Determine if profile is dimensionless.
  if      (present(isDimensionless     )) then
     self%dimensionless=isDimensionless
  else
     self%dimensionless=.false.
  end if
  ! Initialize profile or check for consistency with request index.
  !$omp critical (Mass_Distribution_Sersic_Initialize)
  if (sersicProfileTabulated) then
     ! Profile has already been tabulated. Check that it has the same index as that requested.
     if (Values_Differ(index,sersicProfileIndex,absTol=1.0d-2)) call Galacticus_Error_Report('Mass_Distribution_Sersic_Initialize','requested Sersic index differs from that tabulated')
  else
     ! Store the index.
     sersicProfileIndex=index
     ! Tabulate the profile.
     call Sersic_Profile_Tabulate()
     sersicProfileTabulated=.true.
  end if
  !$omp end critical (Mass_Distribution_Sersic_Initialize)
  ! If dimensionless, then set scale length and mass to unity.
  if (self%dimensionless) then
     if (present(halfMassRadius      )) then
        if (Values_Differ(halfMassRadius      ,1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('Mass_Distribution_Sersic_Initialize','halfMassRadius should be unity for a dimensionless profile (or simply do not specify a half mass radius)')
     end if
     if (present(mass                )) then
        if (Values_Differ(mass                ,1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('Mass_Distribution_Sersic_Initialize','mass should be unity for a dimensionless profile (or simply do not specify a mass)')
     end if
     self%radiusHalfMass=1.0d0
     self%mass          =1.0d0
  else
     if (present(halfMassRadius)) then
        self%radiusHalfMass=halfMassRadius
     else
        call Galacticus_Error_Report('Mass_Distribution_Sersic_Initialize','"halfMassRadius" must be specified')
     end if
     if (present(mass)) then
        self%mass          =mass
     else
        call Galacticus_Error_Report('Mass_Distribution_Sersic_Initialize','"mass" must be specified'       )
     end if
  end if
  return
end subroutine Mass_Distribution_Sersic_Initialize

double precision function Mass_Distribution_Sersic_Density(self,coordinates)
  !% Return the density at the specified {\tt coordinates} in a Sersic mass distribution.
  use Coordinates
  use Numerical_Interpolation
  implicit none
  class           (massDistributionSersic), intent(in   ) :: self
  class           (coordinate            ), intent(in   ) :: coordinates
  type            (coordinateSpherical   )                :: position
  double precision                                        :: r

   ! Get position in spherical coordinate system.
   position=coordinates
   ! Compute the density at this position.
   r=position%r()/self%radiusHalfMass
   call Sersic_Profile_Tabulate(r)
  !$omp critical (Sersic_Profile_Tabulation)
   Mass_Distribution_Sersic_Density=(self%mass/self%radiusHalfMass**3)*Interpolate(sersicTableCount,sersicTableRadius&
        &,sersicTableDensity,sersicTableInterpolationObject,sersicTableInterpolationAccelerator,r,reset &
        &=sersicTableInterpolationReset,extrapolationType=extrapolationTypeLinear)
  !$omp end critical (Sersic_Profile_Tabulation)
  return
end function Mass_Distribution_Sersic_Density

double precision function Mass_Distribution_Sersic_Density_Radial_Moment(self,moment,isInfinite)
  !% Returns a radial density moment for the Sersic mass distribution.
  use Numerical_Constants_Math
  use Numerical_Comparison
  use Galacticus_Error
  implicit none
  class(massDistributionSersic), intent(in   )           :: self
  double precision             , intent(in   )           :: moment
  logical                      , intent(  out), optional :: isInfinite
  integer                                                :: iRadius
  double precision                                       :: integrand,previousIntegrand,deltaRadius

  isInfinite=.false.
  Mass_Distribution_Sersic_Density_Radial_Moment=0.0d0
  !$omp critical (Sersic_Profile_Tabulation)
  do iRadius=1,sersicTableCount
     if (iRadius == 1) then
        deltaRadius      =sersicTableRadius(iRadius)
        previousIntegrand=0.0d0
     else
        deltaRadius      =sersicTableRadius(iRadius)-sersicTableRadius(iRadius-1)
        previousIntegrand=sersicTableRadius(iRadius-1)**moment*sersicTableDensity(iRadius-1)
     end if
     integrand=sersicTableRadius(iRadius)**moment*sersicTableDensity(iRadius)
     Mass_Distribution_Sersic_Density_Radial_Moment=Mass_Distribution_Sersic_Density_Radial_Moment+0.5d0*(previousIntegrand&
          &+integrand)*deltaRadius
  end do
  !$omp end critical (Sersic_Profile_Tabulation)
  Mass_Distribution_Sersic_Density_Radial_Moment=Mass_Distribution_Sersic_Density_Radial_Moment*self%mass&
       &*self%radiusHalfMass**(moment-3.0d0)
  return
end function Mass_Distribution_Sersic_Density_Radial_Moment

double precision function Mass_Distribution_Sersic_Mass_Enc_By_Sphere(self,radius)
  !% Computes the mass enclosed within a sphere of given {\tt radius} for Sersic mass distributions.
  use Numerical_Constants_Math
  use Numerical_Interpolation
  implicit none
  class(massDistributionSersic), intent(in   ), target :: self
  double precision             , intent(in   )         :: radius
  double precision                                     :: fractionalRadius

  fractionalRadius=radius/self%radiusHalfMass
  call Sersic_Profile_Tabulate(fractionalRadius)
  if (fractionalRadius < sersicTableRadius(sersicTableCount)) then
     !$omp critical (Sersic_Profile_Tabulation)
     Mass_Distribution_Sersic_Mass_Enc_By_Sphere=self%mass*Interpolate(sersicTableCount,sersicTableRadius,sersicTableEnclosedMass&
          & ,sersicTableInterpolationObject,sersicTableInterpolationAccelerator,fractionalRadius,reset &
          &=sersicTableInterpolationReset)
     !$omp end critical (Sersic_Profile_Tabulation)
  else
     Mass_Distribution_Sersic_Mass_Enc_By_Sphere=self%mass
  end if
  return
end function Mass_Distribution_Sersic_Mass_Enc_By_Sphere

double precision function Mass_Distribution_Sersic_Potential(self,coordinates)
  !% Return the potential at the specified {\tt coordinates} in a Sersic mass distribution.
  use Numerical_Constants_Physical
  use Coordinates
  use Numerical_Interpolation
  implicit none
  class           (massDistributionSersic), intent(in   ) :: self
  class           (coordinate            ), intent(in   ) :: coordinates
  type            (coordinateSpherical   )                :: position
  double precision                                        :: r

  ! Get position in spherical coordinate system.
  position=coordinates
  ! Compute the potential at this position.
  r=position%r()/self%radiusHalfMass
  call Sersic_Profile_Tabulate(r)
  if (r < sersicTableRadius(sersicTableCount)) then
     !$omp critical (Sersic_Profile_Tabulation)
     Mass_Distribution_Sersic_Potential=(self%mass/self%radiusHalfMass)*Interpolate(sersicTableCount,sersicTableRadius&
          &,sersicTablePotential ,sersicTableInterpolationObject,sersicTableInterpolationAccelerator,r,reset &
          &=sersicTableInterpolationReset)
     !$omp end critical (Sersic_Profile_Tabulation)
  else
     Mass_Distribution_Sersic_Potential=0.0d0
  end if
  if (.not.self%isDimensionless()) Mass_Distribution_Sersic_Potential=gravitationalConstantGalacticus&
       &*Mass_Distribution_Sersic_Potential
  return
end function Mass_Distribution_Sersic_Potential

double precision function Mass_Distribution_Sersic_Half_Mass_Radius(self)
  !% Return the half-mass radius of a Sersic mass distribution.
  implicit none
  class(massDistributionSersic), intent(in   ) :: self

  !$omp critical (Sersic_Profile_Tabulation)
  Mass_Distribution_Sersic_Half_Mass_Radius=self%radiusHalfMass
  !$omp end critical (Sersic_Profile_Tabulation)
  return
end function Mass_Distribution_Sersic_Half_Mass_Radius

double precision function Mass_Distribution_Sersic_Half_Mass_Radius_Projected(self)
  !% Return the half-mass radius in projection of a Sersic mass distribution.
  implicit none
  class(massDistributionSersic), intent(in   ) :: self

  !$omp critical (Sersic_Profile_Tabulation)
  Mass_Distribution_Sersic_Half_Mass_Radius_Projected=self%radiusHalfMass*sersicTable2dHalfMassRadius
  !$omp end critical (Sersic_Profile_Tabulation)
  return
end function Mass_Distribution_Sersic_Half_Mass_Radius_Projected

subroutine Sersic_Profile_Tabulate(radius)
  !% Tabulate the density and enclosed mass in a dimensionless S\'ersic profile.
  use, intrinsic :: ISO_C_Binding                             
  use Memory_Management
  use Numerical_Ranges
  use Numerical_Integration
  use Numerical_Interpolation
  use Numerical_Constants_Math         
  use FGSL
  use Root_Finder
  implicit none
  double precision                            , intent(in), optional :: radius
  double precision                            , parameter            :: radiusMaximumToTabulate   =1000.0d0
  double precision                            , parameter            :: sersicCoefficientTolerance=1.0d-6
  double precision                            , parameter            :: sersicCoefficientGuess    =7.67d0
  logical                                                            :: rebuildTable,tableHasSufficientExtent
  integer                                                            :: iRadius
  double precision                                                   :: deltaRadius,massPrevious,radiusInfinity,integrand&
       &,previousIntegrand,radiusActual,sersicCoefficientMinimum,sersicCoefficientMaximum
  type            (c_ptr                     )                       :: parameterPointer
  type            (fgsl_function             )                       :: integrandFunction
  type            (fgsl_integration_workspace)                       :: integrationWorkspace
  type            (fgsl_function             )                       :: rootFunction
  type            (fgsl_root_fsolver         )                       :: rootFunctionSolver
  
  ! Check if a radius was specified. Use it if so, otherwise use a midpoint radius.
  if (present(radius)) then
     radiusActual=min(radius,radiusMaximumToTabulate)
  else
     radiusActual=dsqrt(sersicTableRadiusMinimum*sersicTableRadiusMaximum)
  end if

  !$omp critical (Sersic_Profile_Tabulation)
  ! Determine if the table must be rebuilt.
  if (sersicTableInitialized) then
     rebuildTable=(radiusActual*sersicTable3dHalfMassRadius < sersicTableRadiusMinimum) &
          &        .or.                               &
          &       (radiusActual*sersicTable3dHalfMassRadius > sersicTableRadiusMaximum)
  else
     rebuildTable=.true.
  end if

  ! Rebuild the table if necessary.
  if (rebuildTable) then
     ! Try building the table until it has sufficient extent to encompass the requested radius.
     tableHasSufficientExtent=.false.
     do while (.not.tableHasSufficientExtent)
        ! Find suitable radius limits.
        sersicTableRadiusMinimum=min(sersicTableRadiusMinimum,0.5d0*radiusActual*sersicTable3dHalfMassRadius)
        sersicTableRadiusMaximum=max(sersicTableRadiusMaximum,2.0d0*radiusActual*sersicTable3dHalfMassRadius)
        ! Determine the number of points at which to tabulate the profile.
        sersicTableCount=int(dlog10(sersicTableRadiusMaximum/sersicTableRadiusMinimum)*dble(sersicTablePointsPerDecade))+1
        ! Allocate arrays for storing the tables.
        if (allocated(sersicTableRadius)) then
           call Dealloc_Array(sersicTableRadius      )
           call Dealloc_Array(sersicTableDensity     )
           call Dealloc_Array(sersicTableEnclosedMass)
           call Dealloc_Array(sersicTablePotential   )
        end if
        call Alloc_Array(sersicTableRadius      ,[sersicTableCount])
        call Alloc_Array(sersicTableDensity     ,[sersicTableCount])
        call Alloc_Array(sersicTableEnclosedMass,[sersicTableCount])
        call Alloc_Array(sersicTablePotential   ,[sersicTableCount])
        ! Create an array of logarithmically distributed radii.
        sersicTableRadius=Make_Range(sersicTableRadiusMinimum,sersicTableRadiusMaximum,sersicTableCount,rangeType=rangeTypeLogarithmic)
        ! Compute the coefficient appearing in the Sérsic profile.
        sersicCoefficientMinimum=sersicCoefficientGuess
        sersicCoefficientMaximum=sersicCoefficientGuess
        do while (Sersic_Coefficient_Root(sersicCoefficientMinimum,parameterPointer) <= 0.0d0)
           sersicCoefficientMinimum=0.5d0*sersicCoefficientMinimum
        end do
        do while (Sersic_Coefficient_Root(sersicCoefficientMaximum,parameterPointer) >= 0.0d0)
           sersicCoefficientMaximum=2.0d0*sersicCoefficientMaximum
        end do
        sersicCoefficient=Root_Find(sersicCoefficientMinimum,sersicCoefficientMaximum,Sersic_Coefficient_Root,parameterPointer&
             &,rootFunction,rootFunctionSolver,toleranceRelative=sersicCoefficientTolerance)
        ! Compute a suitably large approximation to infinite radius for use in integration.
        radiusInfinity=10.0d0*sersicTableRadiusMaximum
        ! Loop over radii and compute the inverse Abel integral required to get the 3D Sérsic profile.
        do iRadius=1,sersicTableCount
           sersicRadiusStart=sersicTableRadius(iRadius)
           sersicTableDensity(iRadius)=Integrate( sersicRadiusStart,radiusInfinity                       &
                &                                ,Sersic_Abel_Integrand,parameterPointer           &
                &                                ,integrandFunction,integrationWorkspace           &
                &                                ,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3 &
                &                               )
           call Integrate_Done(integrandFunction,integrationWorkspace)       
           ! Accumulate the enclosed mass using a simple trapezoidal integration.
           if (iRadius == 1) then
              deltaRadius      =sersicTableRadius(iRadius)
              massPrevious     =0.0d0
              previousIntegrand=0.0d0
           else
              deltaRadius      =sersicTableRadius(iRadius)-sersicTableRadius      (iRadius-1)
              massPrevious     =                           sersicTableEnclosedMass(iRadius-1)
              previousIntegrand=4.0d0*Pi*(sersicTableRadius(iRadius-1)**2)*sersicTableDensity(iRadius-1)
           end if
           integrand           =4.0d0*Pi*(sersicTableRadius(iRadius  )**2)*sersicTableDensity(iRadius  )
           sersicTableEnclosedMass(iRadius)=0.5d0*(previousIntegrand+integrand)*deltaRadius+massPrevious
        end do
        ! Normalize the mass and density to unit total mass.
        sersicTableDensity     =sersicTableDensity     /sersicTableEnclosedMass(sersicTableCount)
        sersicTableEnclosedMass=sersicTableEnclosedMass/sersicTableEnclosedMass(sersicTableCount)
        ! Find the half mass radius.
        call Interpolate_Done(sersicTableInterpolationObject,sersicTableInterpolationAccelerator,sersicTableInterpolationReset)
        sersicTableInterpolationReset=.true.
        sersicTable3dHalfMassRadius=Interpolate(                                                                  &
             &                                                            maxloc(sersicTableEnclosedMass,dim=1),  &
             &                                  sersicTableEnclosedMass(1:maxloc(sersicTableEnclosedMass,dim=1)), &
             &                                  sersicTableRadius      (1:maxloc(sersicTableEnclosedMass,dim=1)), &
             &                                  sersicTableInterpolationObject,                                   &
             &                                  sersicTableInterpolationAccelerator,                              &
             &                                  0.5d0,                                                            &
             &                                  reset=sersicTableInterpolationReset                               &
             &                                 )
        call Interpolate_Done(sersicTableInterpolationObject,sersicTableInterpolationAccelerator,sersicTableInterpolationReset)
        sersicTableInterpolationReset=.true.
        ! Scale radii and densities to be in units of the 3D half mass radius.
        sersicTableRadius =sersicTableRadius /sersicTable3dHalfMassRadius
        sersicTableDensity=sersicTableDensity*sersicTable3dHalfMassRadius**3

        ! Store the 2d half mass radius.
        sersicTable2dHalfMassRadius=1.0d0/sersicTable3dHalfMassRadius

        ! Compute the gravitational potential at each radius using a simple trapezoidal rule integration.
        sersicTablePotential(sersicTableCount)=0.0d0 ! Assume zero potential at effective infinity.
        do iRadius=sersicTableCount-1,1,-1
           sersicTablePotential(iRadius)= sersicTablePotential(iRadius+1)                                      &
                &                        -(                                                                    &
                &                           sersicTableEnclosedMass(iRadius+1)/sersicTableRadius(iRadius+1)**2 &
                &                          +sersicTableEnclosedMass(iRadius  )/sersicTableRadius(iRadius  )**2 &
                &                         )                                                                    &
                &                        *0.5d0                                                                &
                &                        *(                                                                    &
                &                           sersicTableRadius(iRadius+1)                                       &
                &                          -sersicTableRadius(iRadius  )                                       &
                &                         )
        end do

        ! Test that the table has sufficient extent for the requested radius.
        tableHasSufficientExtent=(radiusActual >= sersicTableRadiusMinimum) & 
             &                    .and. &
             &                   (radiusActual <= sersicTableRadiusMaximum)
     end do
     ! Flag that the table is initialized.
     sersicTableInitialized=.true.
  end if
  !$omp end critical (Sersic_Profile_Tabulation)
  return
end subroutine Sersic_Profile_Tabulate

function Sersic_Coefficient_Root(coefficient,parameterPointer) bind(c)
  use Gamma_Functions
  use, intrinsic :: ISO_C_Binding                             
  implicit none
  real(c_double), value   :: coefficient
  type(c_ptr),    value   :: parameterPointer
  real(c_double)          :: Sersic_Coefficient_Root
  
  Sersic_Coefficient_Root=Gamma_Function_Incomplete(2.0d0*sersicProfileIndex,coefficient)-0.5d0
end function Sersic_Coefficient_Root

function Sersic_Abel_Integrand(radius,parameterPointer) bind(c)
  !% The integrand in the Abel integral used to invert the S\'ersic profile to get the corresponding 3-D profile.
  use, intrinsic :: ISO_C_Binding
  use Numerical_Constants_Math
  implicit none
  real(c_double)        :: Sersic_Abel_Integrand
  real(c_double), value :: radius
  type(c_ptr   ), value :: parameterPointer

  if (radius > sersicRadiusStart) then
     Sersic_Abel_Integrand= (1.0d0/Pi)                                                                &
          &                *      sersicCoefficient*(radius**(1.0d0/dble(sersicProfileIndex) -1.0d0)) &
          &                *dexp(-sersicCoefficient*(radius**(1.0d0/dble(sersicProfileIndex))-1.0d0)) &
          &                /dble(sersicProfileIndex)                                                  &
          &                /dsqrt(radius**2-sersicRadiusStart**2)
  else
     Sersic_Abel_Integrand=0.0d0
  end if
  return
end function Sersic_Abel_Integrand

!# <galacticusStateStoreTask>
!#  <unitName>Sersic_Profile_Tabulate_State_Store</unitName>
!# </galacticusStateStoreTask>
subroutine Sersic_Profile_Tabulate_State_Store(stateFile,fgslStateFile)
  !% Write the history state to file.
  use FGSL
  implicit none
  integer           , intent(in   ) :: stateFile
  type   (fgsl_file), intent(in   ) :: fgslStateFile
  
  write (stateFile) sersicTableRadiusMinimum,sersicTableRadiusMaximum,sersicTable3dHalfMassRadius
  return
end subroutine Sersic_Profile_Tabulate_State_Store

!# <galacticusStateRetrieveTask>
!#  <unitName>Sersic_Profile_Tabulate_State_Retrieve</unitName>
!# </galacticusStateRetrieveTask>
subroutine Sersic_Profile_Tabulate_State_Retrieve(stateFile,fgslStateFile)
  !% Retrieve the history state from the file.
  use FGSL
  implicit none
  integer           , intent(in   ) :: stateFile
  type   (fgsl_file), intent(in   ) :: fgslStateFile
  
  read (stateFile) sersicTableRadiusMinimum,sersicTableRadiusMaximum,sersicTable3dHalfMassRadius
  sersicTableInitialized=.false.
  return
end subroutine Sersic_Profile_Tabulate_State_Retrieve
