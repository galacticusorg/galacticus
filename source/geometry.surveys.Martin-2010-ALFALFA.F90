!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Implements the survey geometry used by \cite{martin_arecibo_2010}.

  !# <surveyGeometry name="surveyGeometryMartin2010ALFALFA">
  !#  <description>Implements the survey geometry of the SDSS sample used by \cite{martin_arecibo_2010}.</description>
  !# </surveyGeometry>

  type, extends(surveyGeometryClass) :: surveyGeometryMartin2010ALFALFA
   contains
     procedure :: distanceMaximum => martin2010ALFALFADistanceMaximum
     procedure :: solidAngle      => martin2010ALFALFASolidAngle
     procedure :: windowFunctions => martin2010ALFALFAWindowFunctions
  end type surveyGeometryMartin2010ALFALFA

  interface surveyGeometryMartin2010ALFALFA
     !% Constructors for the \cite{martin_arecibo_2010} survey geometry class.
     module procedure martin2010ALFALFADefaultConstructor
  end interface surveyGeometryMartin2010ALFALFA

  ! Minimum and maximum recession velocities for a galaxy to be admitted to the sample.
  double precision, parameter :: martin2010ALFALFASampleVelocityMinimum=0.0d0, martin2010ALFALFASampleVelocityMaximum=15.0d3

contains

  function martin2010ALFALFADefaultConstructor()
    !% Default constructor for the \cite{martin_arecibo_2010} conditional mass function class.
    use Input_Parameters
    implicit none
    type(surveyGeometryMartin2010ALFALFA) :: martin2010ALFALFADefaultConstructor

    return
  end function martin2010ALFALFADefaultConstructor

  double precision function martin2010ALFALFADistanceMaximum(self,mass,field)
    !% Compute the maximum distance at which a galaxy is visible.
    use Galacticus_Error
    use Cosmology_Parameters
    implicit none
    class           (surveyGeometryMartin2010ALFALFA), intent(inout)           :: self
    double precision                                 , intent(in   )           :: mass
    integer                                          , intent(in   ), optional :: field
    ! The signal-to-noise limit used by Martin et al. (2010).
    double precision                                 , parameter               :: signalToNoise                    = 6.5d0
    ! Coefficients of the polynomial approximation for log10(lineWidth) vs. log10(HI mass).
    double precision                                 , parameter               :: lineWidthCoefficient0            =-0.769635671616885d0, lineWidthCoefficient1=0.314983275066432d0
    ! Line width characteristic scale.
    double precision                                 , parameter               :: lineWidthCharacteristic         =200.0d0
    ! Normalization of the flux limit for unit signal-to-noise at characteristic line width.
    double precision                                 , parameter               :: integratedFluxLimitNormalization=0.15d0
    ! Normalization of the mass-integrated flux-distance relation.
    double precision                                 , parameter               :: massNormalization               =2.356d5
    double precision                                                           :: logarithmicMass                                       , lineWidth                                , &
         &                                                                        integratedFluxLimit
    class    (cosmologyParametersClass              ), pointer                 :: cosmologyParameters_

    ! Validate field.
    if (present(field).and.field /= 1) call Galacticus_Error_Report('martin2010ALFALFADistanceMaximum','field = 1 required')
    ! Get the logarithm of the mass.
    logarithmicMass=log10(mass)
    ! Find the median line width for this mass. (See
    ! constraints/dataAnalysis/hiMassFunction_ALFALFA_z0.00/lineWidthMassRelation.pl for details.)
    lineWidth=10.0d0**(lineWidthCoefficient0+lineWidthCoefficient1*logarithmicMass)
    ! Compute the limiting integrated flux using equation (A1) of Martin et al. (2010).
    if (lineWidth < lineWidthCharacteristic) then
       integratedFluxLimit=integratedFluxLimitNormalization*signalToNoise*sqrt(lineWidth/lineWidthCharacteristic)
    else
       integratedFluxLimit=integratedFluxLimitNormalization*signalToNoise*    (lineWidth/lineWidthCharacteristic)
    end if
    ! Get the default cosmology.
    cosmologyParameters_ => cosmologyParameters()
    ! Convert from mass and limiting integrated flux to maximum distance using relation given in text of section 2.2 of Martin et
    ! al. (2010). Limit by the maximum velocity allowed for galaxies to make it into the sample.
    martin2010ALFALFADistanceMaximum=min(sqrt(mass/massNormalization/integratedFluxLimit),martin2010ALFALFASampleVelocityMaximum/cosmologyParameters_%hubbleConstant())
    return
  end function martin2010ALFALFADistanceMaximum

  double precision function martin2010ALFALFASolidAngle(self,field)
    !% Return the solid angle of the \cite{martin_arecibo_2010} sample.
    use Galacticus_Error
    implicit none
    class           (surveyGeometryMartin2010ALFALFA), intent(inout)           :: self
    integer                                          , intent(in   ), optional :: field
    double precision                                 , parameter               :: solidAngleSurvey=0.79415674617213461d0 ! Computed from survey bounds in Martin et al. (2010; ApJ; 723; 1359)
    
    ! Validate field.
    if (present(field).and.field /= 1) call Galacticus_Error_Report('martin2010ALFALFASolidAngle','field = 1 required')
    martin2010ALFALFASolidAngle=solidAngleSurvey
    return
  end function martin2010ALFALFASolidAngle
  
  subroutine martin2010ALFALFAWindowFunctions(self,mass1,mass2,gridCount,boxLength,windowFunction1,windowFunction2)
    !% Compute the window function for the survey.
    use FFTW3
    use Vectors
    use Pseudo_Random
    use File_Utilities
    use FGSL
    use Meshes
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Math
    use Numerical_Constants_Astronomical
    use Memory_Management
    use Galacticus_Display
    use ISO_Varying_String
    use Cosmology_Functions
    use Cosmology_Parameters
    use String_Handling
    use File_Utilities
    use Galacticus_Input_Paths
    use System_Command
    use Galacticus_Error
    use Pseudo_Random
    implicit none
    class           (surveyGeometryMartin2010ALFALFA), intent(inout) :: self
    double precision                                 , intent(in   )                                               :: mass1,mass2
    integer                                          , intent(in   )                                               :: gridCount
    double precision                                 , intent(  out)                                               :: boxLength
    complex         (c_double_complex               ), intent(  out),     dimension(gridCount,gridCount,gridCount) :: windowFunction1,windowFunction2
    double precision                                 ,                    dimension(3                            ) :: origin,position1,position2
    double precision                                 ,                    dimension(2                            ) :: phiRange,thetaRange
    type            (c_ptr                          )                                                              :: plan
    complex         (c_double_complex               ),                    dimension(gridCount,gridCount,gridCount) :: selectionFunction1,selectionFunction2
    complex         (c_double_complex               )                                                              :: normalization
    logical                                          , save                                                        :: geometryInitialized=.false.
    double precision                                 , save, allocatable, dimension(:                            ) :: randomTheta,randomPhi
    integer                                          , parameter                                                   :: randomsCount=1000000
    type            (fgsl_rng                       ), save                                                        :: pseudoSequenceObject
    logical                                          , save                                                        :: reset=.true.
    type            (varying_string                 )                                                              :: message
    double precision                                 , save                                                        :: surveyDistanceMinimum&
         &,surveyDistanceMaximum
    integer                                          , parameter                                                   :: regionCount=3
    ! Survey geometry from Haynes et al. (2011; http://adsabs.harvard.edu/abs/2011AJ....142..170H).
    double precision                                 , parameter        , dimension(2,regionCount                ) :: regionRightAscensionRange=reshape([22.0d0,03.0d0,07.5d0,16.5d0,07.5d0,16.5d0],[2,regionCount]), &
         &                                                                                                            regionDeclinationRange   =reshape([24.0d0,32.0d0,04.0d0,16.0d0,24.0d0,28.0d0],[2,regionCount])
    double precision                                                    , dimension(2,regionCount                ) :: regionPhiRange,regionThetaRange
    double precision                                                    , dimension(  regionCount                ) :: regionSolidAngle
    integer                                                                                                        :: i,j,iRegion,iRandom
    double precision                                                                                               :: comovingDistanceMaximum1, comovingDistanceMaximum2, &
         &                                                                                                            comovingDistanceMinimum1, comovingDistanceMinimum2, &
         &                                                                                                            rightAscension          , declination             , &
         &                                                                                                            distance1               , distance2               , &
         &                                                                                                            uniformRandom
    class           (cosmologyParametersClass       ), pointer                                                      :: cosmologyParameters_

    ! Initialize geometry if necessary.
    if (.not.geometryInitialized) then       
       ! Determine the solid angles of the different survey regions.
       regionPhiRange  =        +regionRightAscensionRange*  hoursToRadians
       regionThetaRange=0.5d0*Pi-   regionDeclinationRange*degreesToRadians
       do iRegion=1,regionCount
          if (regionPhiRange(2,iRegion) < regionPhiRange(1,iRegion)) regionPhiRange(2,iRegion)=regionPhiRange(2,iRegion)+2.0d0*Pi
       end do
       regionSolidAngle=(regionPhiRange(2,:)-regionPhiRange(1,:))*(cos(regionThetaRange(2,:))-cos(regionThetaRange(1,:)))
       ! Cumulate and normalize the region solid angles.
       do iRegion=2,regionCount
          regionSolidAngle(iRegion)=regionSolidAngle(iRegion)+regionSolidAngle(iRegion-1)
       end do
       regionSolidAngle=regionSolidAngle/regionSolidAngle(regionCount)
       ! Generate random points.
       call Alloc_Array(randomTheta,[randomsCount])
       call Alloc_Array(randomPhi  ,[randomsCount])
       do iRandom=1,randomsCount
          ! Select a region at random.
          uniformRandom=Pseudo_Random_Get(pseudoSequenceObject,reset)
          iRegion=1
          do while (uniformRandom > regionSolidAngle(iRegion))
             iRegion=iRegion+1
          end do
          ! Select coordinates at random within this region.
          uniformRandom=Pseudo_Random_Get(pseudoSequenceObject,reset)
          randomPhi  (iRandom)=     uniformRandom*(    regionPhiRange  (2,iRegion) -    regionPhiRange  (1,iRegion)) +    regionPhiRange  (1,iRegion)
          uniformRandom=Pseudo_Random_Get(pseudoSequenceObject,reset)
          randomTheta(iRandom)=acos(uniformRandom*(cos(regionThetaRange(2,iRegion))-cos(regionThetaRange(1,iRegion)))+cos(regionThetaRange(1,iRegion)))
          uniformRandom=Pseudo_Random_Get(pseudoSequenceObject,reset)
       end do
       ! Get the default cosmology.
       cosmologyParameters_ => cosmologyParameters()
      ! Compute the distances corresponding to the minimum and maximum redshifts.
       surveyDistanceMinimum=martin2010ALFALFASampleVelocityMinimum/cosmologyParameters_%hubbleConstant()
       surveyDistanceMaximum=martin2010ALFALFASampleVelocityMaximum/cosmologyParameters_%hubbleConstant()
       geometryInitialized=.true.
    end if
    ! Find the comoving distance corresponding to this distance module.
    comovingDistanceMaximum1=min(self%distanceMaximum(mass1),surveyDistanceMaximum)
    comovingDistanceMaximum2=min(self%distanceMaximum(mass2),surveyDistanceMaximum)
    comovingDistanceMinimum1=surveyDistanceMinimum
    comovingDistanceMinimum2=surveyDistanceMinimum
    ! Determine a suitable box length for the window function calculation.
    boxLength=3.0d0*max(comovingDistanceMaximum1,comovingDistanceMaximum2)
    ! Set up origin.
    origin=([0.5d0,0.5d0,0.5d0]/dble(gridCount)+[0.5d0,0.5d0,0.5d0])*dble(boxLength)
    ! Populate the cube with the survey selection function.
    selectionFunction1=cmplx(0.0d0,0.0d0,kind=c_double_complex)
    selectionFunction2=cmplx(0.0d0,0.0d0,kind=c_double_complex)
    ! Loop over randoms.
    do i=1,randomsCount
       ! Choose random distances.
       distance1=(Pseudo_Random_Get(pseudoSequenceObject,reset)*(comovingDistanceMaximum1**3-comovingDistanceMinimum1**3)+comovingDistanceMinimum1**3)**(1.0d0/3.0d0)
       distance2=(Pseudo_Random_Get(pseudoSequenceObject,reset)*(comovingDistanceMaximum2**3-comovingDistanceMinimum2**3)+comovingDistanceMinimum2**3)**(1.0d0/3.0d0)
       ! Convert to Cartesian coordinates.
       position1=distance1*[sin(randomTheta(i))*cos(randomPhi(i)),sin(randomTheta(i))*sin(randomPhi(i)),cos(randomTheta(i))]+origin
       position2=distance2*[sin(randomTheta(i))*cos(randomPhi(i)),sin(randomTheta(i))*sin(randomPhi(i)),cos(randomTheta(i))]+origin
       ! Apply to the mesh.
       call Meshes_Apply_Point(selectionFunction1,boxLength,position1,pointWeight=cmplx(1.0d0,0.0d0,kind=c_double_complex),cloudType=cloudTypeTriangular)
       call Meshes_Apply_Point(selectionFunction2,boxLength,position2,pointWeight=cmplx(1.0d0,0.0d0,kind=c_double_complex),cloudType=cloudTypeTriangular)
    end do
    ! Take the Fourier transform of the selection function.
    plan=fftw_plan_dft_3d(gridCount,gridCount,gridCount,selectionFunction1,windowFunction1,FFTW_FORWARD,FFTW_ESTIMATE)
    call fftw_execute_dft(plan,selectionFunction1,windowFunction1)
    call fftw_execute_dft(plan,selectionFunction2,windowFunction2)
    call fftw_destroy_plan(plan)
    ! Normalize the window function.
    normalization=windowFunction1(1,1,1)
    if (real(normalization) > 0.0d0) windowFunction1=windowFunction1/normalization
    normalization=windowFunction2(1,1,1)
    if (real(normalization) > 0.0d0) windowFunction2=windowFunction2/normalization
    return
  end subroutine martin2010ALFALFAWindowFunctions
