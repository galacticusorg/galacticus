!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  !!{
  Implementation of a kinematic distribution class for finite-resolution NFW mass distributions.
  !!}

  use :: Numerical_Interpolation, only : interpolator

  !![
  <kinematicsDistribution name="kinematicsDistributionFiniteResolutionNFW">
   <description>A kinematic distribution class for finite-resolution NFW mass distributions.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionFiniteResolutionNFW
     !!{
     A kinematics distribution for finite-resolution NFW distributions.
     !!}
     double precision                                            :: velocityDispersion1DRadiusPrevious                   , velocityDispersion1DPrevious
     ! Velocity dispersion tabulation.
     logical                                                     :: velocityDispersion1DTableInitialized
     integer                                                     :: velocityDispersion1DTableLengthResolutionCount       , velocityDispersion1DTableRadiusCount
     double precision              , allocatable, dimension(:  ) :: velocityDispersion1DTableLengthResolution            , velocityDispersion1DTableRadius
     double precision              , allocatable, dimension(:,:) :: velocityDispersion1DTable
     type            (interpolator), allocatable                 :: velocityDispersion1DTableLengthResolutionInterpolator, velocityDispersion1DTableRadiusInterpolator
     double precision                                            :: velocityDispersion1DRadiusMinimum                    , velocityDispersion1DRadiusMaximum          , &
          &                                                         velocityDispersion1DLengthResolutionMinimum          , velocityDispersion1DLengthResolutionMaximum
   contains
     !![
     <methods>
       <method method="velocityDispersion1DTabulate"   description="Tabulate the enclosed mass as a function of radius and core radius."                          />
       <method method="storeVelocityDispersionTable"   description="Store the tabulated velocity dispersion to file."                                             />
       <method method="restoreVelocityDispersionTable" description="Attempt to restore the tabulated velocity dispersion from file, returning true if successful."/>
     </methods>
     !!]
     procedure :: isCollisional                  => finiteResolutionNFWIsCollisional
     procedure :: velocityDispersion1D           => finiteResolutionNFWVelocityDispersion1D
     procedure :: velocityDispersion1DTabulate   => finiteResolutionNFWVelocityDispersionRadialTabulate
     procedure :: storeVelocityDispersionTable   => finiteResolutionNFWStoreVelocityDispersionTable
     procedure :: restoreVelocityDispersionTable => finiteResolutionNFWRestoreVelocityDispersionTable
  end type kinematicsDistributionFiniteResolutionNFW

  interface kinematicsDistributionFiniteResolutionNFW
     !!{
     Constructors for the {\normalfont \ttfamily nfw} kinematic distribution class.
     !!}
     module procedure finiteResolutionNFWConstructorParameters
     module procedure finiteResolutionNFWConstructorInternal
  end interface kinematicsDistributionFiniteResolutionNFW

  ! Tabulation resolution parameters.
  integer, parameter :: velocityDispersion1DTableRadiusPointsPerDecade          =100
  integer, parameter :: velocityDispersion1DTableLengthResolutionPointsPerDecade=100

  class(kinematicsDistributionFiniteResolutionNFW   ), pointer :: self_
  class(massDistributionSphericalFiniteResolutionNFW), pointer :: massDistributionEmbedding_
  integer                                                      :: iLengthResolution_
  !$omp threadprivate(iLengthResolution_,self_,massDistributionEmbedding_)

contains

  function finiteResolutionNFWConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily finiteResolutionNFW} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(kinematicsDistributionFiniteResolutionNFW)                :: self
    type(inputParameters                          ), intent(inout) :: parameters

    self=kinematicsDistributionFiniteResolutionNFW()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function finiteResolutionNFWConstructorParameters

  function finiteResolutionNFWConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily finiteResolutionNFW} kinematic distribution class.
    !!}
    implicit none
    type(kinematicsDistributionFiniteResolutionNFW) :: self

    self%velocityDispersion1DRadiusMinimum          =+huge(0.0d0)
    self%velocityDispersion1DRadiusMaximum          =-huge(0.0d0)
    self%velocityDispersion1DLengthResolutionMinimum=+huge(0.0d0)
    self%velocityDispersion1DLengthResolutionMaximum=-huge(0.0d0)
    self%velocityDispersion1DRadiusPrevious         =-huge(0.0d0)
    self%velocityDispersion1DPrevious               =-huge(0.0d0)
    return
  end function finiteResolutionNFWConstructorInternal
  
  logical function finiteResolutionNFWIsCollisional(self)
    !!{
    Return true indicating that the finiteResolutionNFW kinematic distribution represents collisional particles.
    !!}
    implicit none
    class(kinematicsDistributionFiniteResolutionNFW), intent(inout) :: self
    
    finiteResolutionNFWIsCollisional=.false.
    return
  end function finiteResolutionNFWIsCollisional

  double precision function finiteResolutionNFWVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in an finite-resolution NFW kinematic distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (kinematicsDistributionFiniteResolutionNFW), intent(inout)           :: self
    class           (coordinate                               ), intent(in   )           :: coordinates
    class           (massDistributionClass                    ), intent(inout) , target  :: massDistribution_                    , massDistributionEmbedding
    class           (massDistributionClass                    )                , pointer :: massDistribution__
    double precision                                           , parameter               :: lengthResolutionScaleFreeSmall=1.0d-3
    integer         (c_size_t                                 ), dimension(0:1)          :: jLengthResolution
    double precision                                           , dimension(0:1)          :: hLengthResolution
    integer                                                                              :: iLengthResolution
    double precision                                                                     :: radiusScaleFree                      , radiusScaleFreeEffective

    massDistribution__ => massDistribution_
    if (associated(massDistribution__,massDistributionEmbedding)) then
       ! For the case of a self-gravitating finite-resolution NFW distribution we have a tabulated solution for the velocity dispersion.
       select type (massDistributionEmbedding)
       class is (massDistributionSphericalFiniteResolutionNFW)
          if (coordinates%rSpherical() /= self%velocityDispersion1DRadiusPrevious) then
             self%velocityDispersion1DRadiusPrevious=coordinates%rSpherical()
             ! Compute the effective radius. In the core of the profile the velocity dispersion must become constant. Therefore, we
             ! limit the smallest radius we consider to a small fraction of the core radius. Below this radius a constant velocity
             ! dispersion is assumed.
             radiusScaleFree         =coordinates%rSpherical()/massDistributionEmbedding%radiusScale
             radiusScaleFreeEffective=max(radiusScaleFree,lengthResolutionScaleFreeSmall*massDistributionEmbedding%lengthResolutionScaleFree)
             ! Ensure table is sufficiently extensive.
             call self%velocityDispersion1DTabulate(massDistributionEmbedding,radiusScaleFreeEffective,massDistributionEmbedding%lengthResolutionScaleFree)
             ! Interpolate to get the velocity dispersion.
             call self%velocityDispersion1DTableLengthResolutionInterpolator%linearFactors(massDistributionEmbedding%lengthResolutionScaleFree,jLengthResolution(0),hLengthResolution)
             jLengthResolution(1)=jLengthResolution(0)+1
             self%velocityDispersion1DPrevious=0.0d0
             do iLengthResolution=0,1
                self%velocityDispersion1DPrevious=+self%velocityDispersion1DPrevious                                                                                                                             &
                     &                            +self%velocityDispersion1DTableRadiusInterpolator%interpolate(radiusScaleFreeEffective,self%velocityDispersion1DTable(:,jLengthResolution(iLengthResolution))) &
                     &                            *                                                                                                                       hLengthResolution(iLengthResolution)
             end do
             self%velocityDispersion1DPrevious=+self%velocityDispersion1DPrevious                       &
                  &                            *sqrt(                                                   &
                  &                                  +gravitationalConstant_internal                    &
                  &                                  *massDistributionEmbedding%densityNormalization    &
                  &                                  *massDistributionEmbedding%radiusScale         **2 &
                  &                                 )
          end if
          velocityDispersion=self%velocityDispersion1DPrevious
       class default
          velocityDispersion=0.0d0
          call Error_Report('expecting a finite-resolution NFW mass distribution, but received '//char(massDistributionEmbedding%objectType())//{introspection:location})
       end select
    else
       ! Our finite resolution NFW distribution is embedded in another distribution. We must compute the velocity dispersion numerically.
       velocityDispersion=self%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
    end if
    return
  end function finiteResolutionNFWVelocityDispersion1D
  
  subroutine finiteResolutionNFWVelocityDispersionRadialTabulate(self,massDistributionEmbedding,radius,radiusCore)
    !!{
    Tabulates the mass enclosed within a given radius for finite resolution NFW density profiles.
    !!}
    use :: Numerical_Ranges     , only : Make_Range, rangeTypeLogarithmic
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (kinematicsDistributionFiniteResolutionNFW   ), intent(inout), target :: self
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout), target :: massDistributionEmbedding
    double precision                                              , intent(in   )         :: radius                           , radiusCore
    double precision                                              , parameter             :: radiusTiny               =1.0d-2
    type            (integrator                                  ), save                  :: integrator_
    logical                                                       , save                  :: initialized              =.false.
    !$omp threadprivate(integrator_,initialized)
    logical                                                                               :: retabulate
    integer                                                                               :: iLengthResolution                , iRadius              , &
         &                                                                                   i
    double precision                                                                      :: jeansIntegral                    , jeansIntegralPrevious, &
         &                                                                                   radiusLower                      , radiusUpper          , &
         &                                                                                   radiusOuter                      , density

    do i=1,2
       retabulate=.false.
       if (.not.self%velocityDispersion1DTableInitialized) then
          retabulate=.true.
       else if (                                                               &
            &    radius     < self%velocityDispersion1DRadiusMinimum           &
            &   .or.                                                           &
            &    radius     > self%velocityDispersion1DRadiusMaximum           &
            &   .or.                                                           &
            &    radiusCore < self%velocityDispersion1DLengthResolutionMinimum &
            &   .or.                                                           &
            &    radiusCore > self%velocityDispersion1DLengthResolutionMaximum &
            &  ) then
          retabulate=.true.
       end if
       if (retabulate     .and.i==1) call self%restoreVelocityDispersionTable()
       if (.not.retabulate         ) exit
    end do
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%velocityDispersion1DRadiusMinimum             =min(0.5d0*radius    ,self%velocityDispersion1DRadiusMinimum          )
       self%velocityDispersion1DRadiusMaximum             =max(2.0d0*radius    ,self%velocityDispersion1DRadiusMaximum          )
       self%velocityDispersion1DLengthResolutionMinimum   =min(0.5d0*radiusCore,self%velocityDispersion1DLengthResolutionMinimum)
       self%velocityDispersion1DLengthResolutionMaximum   =max(2.0d0*radiusCore,self%velocityDispersion1DLengthResolutionMaximum)
       self%velocityDispersion1DTableRadiusCount          =int(log10(self%velocityDispersion1DRadiusMaximum          /self%velocityDispersion1DRadiusMinimum          )*dble(velocityDispersion1DTableRadiusPointsPerDecade          ))+1
       self%velocityDispersion1DTableLengthResolutionCount=int(log10(self%velocityDispersion1DLengthResolutionMaximum/self%velocityDispersion1DLengthResolutionMinimum)*dble(velocityDispersion1DTableLengthResolutionPointsPerDecade))+1
       if (allocated(self%velocityDispersion1DTableRadius)) then
          deallocate(self%velocityDispersion1DTableLengthResolution)
          deallocate(self%velocityDispersion1DTableRadius    )
          deallocate(self%velocityDispersion1DTable          )
       end if
       allocate(self%velocityDispersion1DTableLengthResolution(                                          self%velocityDispersion1DTableLengthResolutionCount))
       allocate(self%velocityDispersion1DTableRadius          (self%velocityDispersion1DTableRadiusCount                                                    ))
       allocate(self%velocityDispersion1DTable                (self%velocityDispersion1DTableRadiusCount,self%velocityDispersion1DTableLengthResolutionCount))
       ! Create a range of radii and core radii.
       self%velocityDispersion1DTableRadius          =Make_Range(self%velocityDispersion1DRadiusMinimum          ,self%velocityDispersion1DRadiusMaximum          ,self%velocityDispersion1DTableRadiusCount          ,rangeType=rangeTypeLogarithmic)
       self%velocityDispersion1DTableLengthResolution=Make_Range(self%velocityDispersion1DLengthResolutionMinimum,self%velocityDispersion1DLengthResolutionMaximum,self%velocityDispersion1DTableLengthResolutionCount,rangeType=rangeTypeLogarithmic)
       ! Initialize integrator if necessary.
       if (.not.initialized) then
          integrator_=integrator(jeansEquationIntegrand,toleranceRelative=1.0d-2)
          initialized=.true.
       end if
       ! Loop over radii and α and populate tables.
       self_                      => self
       massDistributionEmbedding_ => massDistributionEmbedding
       radiusOuter                =  max(10.0d0*self%velocityDispersion1DRadiusMaximum,1000.0d0)
       do iLengthResolution=1,self%velocityDispersion1DTableLengthResolutionCount
          iLengthResolution_   =iLengthResolution
          jeansIntegralPrevious=0.0d0
          do iRadius=self%velocityDispersion1DTableRadiusCount,1,-1
             ! For radii that are tiny compared to the core radius the velocity dispersion become almost constant. Simply assume this to avoid floating point errors.
             if     (                                                                                                                              &
                  &   self%velocityDispersion1DTableRadius(iRadius) < radiusTiny                                                                   &
                  &  .and.                                                                                                                         &
                  &   self%velocityDispersion1DTableRadius(iRadius) < radiusTiny*self%velocityDispersion1DTableLengthResolution(iLengthResolution) &
                  &  .and.                                                                                                                         &
                  &   iRadius                                       < self%velocityDispersion1DTableRadiusCount                                    &
                  & ) then
                self%velocityDispersion1DTable(iRadius,iLengthResolution)=self%velocityDispersion1DTable(iRadius+1,iLengthResolution)
             else
                ! Find the limits for the integral.
                if (iRadius == self%velocityDispersion1DTableRadiusCount) then
                   radiusUpper=radiusOuter
                else
                   radiusUpper=self%velocityDispersion1DTableRadius(iRadius+1)
                end if
                radiusLower                                              =self                     %velocityDispersion1DTableRadius(                                                           iRadius           )
                density                                                  =massDistributionEmbedding%densityScaleFree               (radiusLower,self%velocityDispersion1DTableLengthResolution(iLengthResolution))
                jeansIntegral                                            =integrator_              %integrate                      (radiusLower,     radiusUpper                                                 )
                self%velocityDispersion1DTable(iRadius,iLengthResolution)=+sqrt(                         &
                     &                                                          +(                       &
                     &                                                            +jeansIntegral         &
                     &                                                            +jeansIntegralPrevious &
                     &                                                           )                       &
                     &                                                          /density                 &
                     &                                                         )
                jeansIntegralPrevious                                    =+jeansIntegralPrevious &
                     &                                                    +jeansIntegral
             end if
          end do
       end do
       ! Build interpolators.
       if (allocated(self%velocityDispersion1DTableLengthResolutionInterpolator)) deallocate(self%velocityDispersion1DTableLengthResolutionInterpolator)
       if (allocated(self%velocityDispersion1DTableRadiusInterpolator          )) deallocate(self%velocityDispersion1DTableRadiusInterpolator          )
       allocate(self%velocityDispersion1DTableLengthResolutionInterpolator)
       allocate(self%velocityDispersion1DTableRadiusInterpolator          )
       self%velocityDispersion1DTableLengthResolutionInterpolator=interpolator(self%velocityDispersion1DTableLengthResolution)
       self%velocityDispersion1DTableRadiusInterpolator          =interpolator(self%velocityDispersion1DTableRadius          )
       ! Specify that tabulation has been made.
       self%velocityDispersion1DTableInitialized=.true.
       call self%storeVelocityDispersionTable()
    end if
    return
  end subroutine finiteResolutionNFWVelocityDispersionRadialTabulate
  
  double precision function jeansEquationIntegrand(radius)
    !!{
    Integrand for dark matter profile Jeans equation.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    
    if (radius > 0.0d0) then
       jeansEquationIntegrand=+massDistributionEmbedding_%massEnclosedScaleFree(radius,self_%velocityDispersion1DTableLengthResolution(iLengthResolution_))    &
            &                 *massDistributionEmbedding_%densityScaleFree     (radius,self_%velocityDispersion1DTableLengthResolution(iLengthResolution_))    &
            &                 /                                                 radius                                                                     **2
    else
       jeansEquationIntegrand=0.0d0
    end if
    return
  end function jeansEquationIntegrand

  subroutine finiteResolutionNFWStoreVelocityDispersionTable(self)
    !!{
    Store the tabulated velocity dispersion data to file.
    !!}
    use :: File_Utilities    , only : File_Lock     , File_Unlock        , lockDescriptor, Directory_Make, &
         &                            File_Path
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)       , char
    implicit none
    class(kinematicsDistributionFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                           )                :: fileLock
    type (hdf5Object                               )                :: file
    type (varying_string                           )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'VelocityDispersion_'                            // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    call Directory_Make(char(File_Path(char(fileName))))
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(char(fileName),fileLock,lockIsShared=.false.)
    !$ call hdf5Access%set()
    call file%openFile(char(fileName),overWrite=.true.,objectsOverwritable=.true.,readOnly=.false.)
    call file%writeDataset(self%velocityDispersion1DTableLengthResolution,'radiusCore'        )
    call file%writeDataset(self%velocityDispersion1DTableRadius          ,'radius'            )
    call file%writeDataset(self%velocityDispersion1DTable                ,'velocityDispersion')
    call file%close()
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    return
  end subroutine finiteResolutionNFWStoreVelocityDispersionTable

  subroutine finiteResolutionNFWRestoreVelocityDispersionTable(self)
    !!{
    Restore the tabulated velocity dispersion data from file.
    !!}
    use :: File_Utilities    , only : File_Exists   , File_Lock          , File_Unlock, lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)
    implicit none
    class(kinematicsDistributionFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                           )                :: fileLock
    type (hdf5Object                               )                :: file
    type (varying_string                           )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'VelocityDispersion_'                            // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    if (File_Exists(fileName)) then
       if (allocated(self%velocityDispersion1DTableRadius)) then
          deallocate(self%velocityDispersion1DTableLengthResolution)
          deallocate(self%velocityDispersion1DTableRadius          )
          deallocate(self%velocityDispersion1DTable                )
       end if
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
       !$ call hdf5Access%set()
       call file%openFile(char(fileName))
       call file%readDataset('radiusCore'        ,self%velocityDispersion1DTableLengthResolution)
       call file%readDataset('radius'            ,self%velocityDispersion1DTableRadius          )
       call file%readDataset('velocityDispersion',self%velocityDispersion1DTable                )
       call file%close()
       !$ call hdf5Access%unset()
       call File_Unlock(fileLock)
       self%velocityDispersion1DTableRadiusCount          =size(self%velocityDispersion1DTableRadius          )
       self%velocityDispersion1DTableLengthResolutionCount=size(self%velocityDispersion1DTableLengthResolution)
       self%velocityDispersion1DRadiusMinimum             =self%velocityDispersion1DTableRadius          (                                                  1)
       self%velocityDispersion1DRadiusMaximum             =self%velocityDispersion1DTableRadius          (self%velocityDispersion1DTableRadiusCount          )
       self%velocityDispersion1DLengthResolutionMinimum   =self%velocityDispersion1DTableLengthResolution(                                                  1)
       self%velocityDispersion1DLengthResolutionMaximum   =self%velocityDispersion1DTableLengthResolution(self%velocityDispersion1DTableLengthResolutionCount)
       if (allocated(self%velocityDispersion1DTableLengthResolutionInterpolator)) deallocate(self%velocityDispersion1DTableLengthResolutionInterpolator)
       if (allocated(self%velocityDispersion1DTableRadiusInterpolator          )) deallocate(self%velocityDispersion1DTableRadiusInterpolator          )
       allocate(self%velocityDispersion1DTableLengthResolutionInterpolator)
       allocate(self%velocityDispersion1DTableRadiusInterpolator          )
       self%velocityDispersion1DTableLengthResolutionInterpolator=interpolator(self%velocityDispersion1DTableLengthResolution)
       self%velocityDispersion1DTableRadiusInterpolator          =interpolator(self%velocityDispersion1DTableRadius          )
       self%velocityDispersion1DTableInitialized                 =.true.
    end if    
    return
  end subroutine finiteResolutionNFWRestoreVelocityDispersionTable
