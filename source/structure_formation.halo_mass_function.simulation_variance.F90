!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Implements a dark matter halo mass function class which modifies another mass function by accounting for cosmic variance in a
simulation cube.
!!}

  use :: Dark_Matter_Halo_Biases, only : darkMatterHaloBiasClass
  use :: Power_Spectra          , only : powerSpectrumClass
  use :: Numerical_Integration  , only : integrator

  !![
  <haloMassFunction name="haloMassFunctionSimulationVariance">
   <description>
    The halo mass function is computed by modifying another mass function to account for cosmic variance in a simulation cube.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionSimulationVariance
     !!{
     A halo mass function class that modifies another mass function by accounting for cosmic variance in a
     simulation cube.
     !!}
     private
     double precision                                   :: lengthSimulationCube          , varianceSimulation    , &
          &                                                timePrevious                  , perturbationFractional
     class           (haloMassFunctionClass  ), pointer :: massFunction_        => null()
     class           (darkMatterHaloBiasClass), pointer :: darkMatterHaloBias_  => null()
     class           (powerSpectrumClass     ), pointer :: powerSpectrum_       => null()
   contains
     final     ::                 simulationVarianceDestructor
     procedure :: differential => simulationVarianceDifferential
  end type haloMassFunctionSimulationVariance

  interface haloMassFunctionSimulationVariance
     !!{
     Constructors for the {\normalfont \ttfamily simulationVariance} halo mass function class.
     !!}
     module procedure simulationVarianceConstructorParameters
     module procedure simulationVarianceConstructorInternal
  end interface haloMassFunctionSimulationVariance

  ! Sub-module scope variables used in compute the variance in the simulation cube.
  class           (powerSpectrumClass), pointer :: powerSpectrum_
  type            (integrator        )          :: integratorX            , integratorY      , &
       &                                           integratorZ
  double precision                              :: wavenumberX            , windowFunctionX  , &
       &                                           wavenumberY            , windowFunctionY  , &
       &                                           wavenumberZ            , windowFunctionZ  , &
       &                                           wavenumberMinimum      , wavenumberMaximum, &
       &                                           lengthSimulationCube_  , time_
  !$omp threadprivate(powerSpectrum_,wavenumberX,windowFunctionX,wavenumberY,windowFunctionY,wavenumberZ,windowFunctionZ,lengthSimulationCube_,integratorX,integratorY,integratorZ,wavenumberMinimum,wavenumberMaximum,time_)

  ! Cached copies of computed variances. These are used to avoid re-reading from file if the same solution is requested multiple times.
  type :: cachedVariance
     type            (varying_string) :: fileName
     double precision                 :: variance
  end type cachedVariance

  integer                , parameter            :: sizeCache      =25
  integer                                       :: countCache     = 0, lastCache=0
  type   (cachedVariance), dimension(sizeCache) :: cachedVariances
  
contains

  function simulationVarianceConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily simulationVariance} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloMassFunctionSimulationVariance)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (haloMassFunctionClass             ), pointer       :: massFunction_
    class           (cosmologyParametersClass          ), pointer       :: cosmologyParameters_
    class           (darkMatterHaloBiasClass           ), pointer       :: darkMatterHaloBias_
    class           (powerSpectrumClass                ), pointer       :: powerSpectrum_
    double precision                                                    :: lengthSimulationCube, perturbationFractional
    
    !![
    <inputParameter>
      <name>lengthSimulationCube</name>
      <description>The length of the simulation cube from which the target mass function was derived.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>perturbationFractional</name>
      <description>The fractional perturbation (in units of the simulation cube root-variance) to apply.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="haloMassFunction"    name="massFunction_"        source="parameters"/>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="darkMatterHaloBias"  name="darkMatterHaloBias_"  source="parameters"/>
    <objectBuilder class="powerSpectrum"       name="powerSpectrum_"       source="parameters"/>
    !!]
    self=haloMassFunctionSimulationVariance(lengthSimulationCube,perturbationFractional,massFunction_,cosmologyParameters_,powerSpectrum_,darkMatterHaloBias_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massFunction_"       />
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="darkMatterHaloBias_" />
    <objectDestructor name="powerSpectrum_"      />
    !!]
    return
  end function simulationVarianceConstructorParameters

  function simulationVarianceConstructorInternal(lengthSimulationCube,perturbationFractional,massFunction_,cosmologyParameters_,powerSpectrum_,darkMatterHaloBias_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily simulationVariance} halo mass function class.
    !!}
    implicit none
    type            (haloMassFunctionSimulationVariance)                        :: self
    class           (haloMassFunctionClass             ), intent(in   ), target :: massFunction_
    class           (cosmologyParametersClass          ), intent(in   ), target :: cosmologyParameters_
    class           (darkMatterHaloBiasClass           ), intent(in   ), target :: darkMatterHaloBias_
    class           (powerSpectrumClass                ), intent(in   ), target :: powerSpectrum_
    double precision                                    , intent(in   )         :: lengthSimulationCube, perturbationFractional
    !![
    <constructorAssign variables="lengthSimulationCube, perturbationFractional, *cosmologyParameters_, *massFunction_, *powerSpectrum_, *darkMatterHaloBias_"/>
    !!]

    self%timePrevious=-huge(0.0d0)
    return
  end function simulationVarianceConstructorInternal

  subroutine simulationVarianceDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily simulationVariance} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionSimulationVariance), intent(inout) :: self

    !![
    <objectDestructor name="self%massFunction_"       />
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%darkMatterHaloBias_" />
    <objectDestructor name="self%powerSpectrum_"      />
    !!]
    return
  end subroutine simulationVarianceDestructor

  double precision function simulationVarianceDifferential(self,time,mass,node) result(massFunction)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    use :: Display       , only : displayIndent, displayUnindent    , displayMessage, verbosityLevelWorking
    use :: HDF5_Access   , only : hdf5Access
    use :: IO_HDF5       , only : hdf5Object
    use :: File_Utilities, only : File_Exists  , File_Lock          , File_Unlock   , lockDescriptor
    use :: Input_Paths   , only : inputPath    , pathTypeDataDynamic
    implicit none
    class           (haloMassFunctionSimulationVariance), intent(inout), target   :: self
    double precision                                    , intent(in   )           :: time                             , mass
    type            (treeNode                          ), intent(inout), optional :: node
    double precision                                    , parameter               :: wavenumberMaximumFractional=1.0d2
 
    ! Find the variance of a cosmological box simulation if the time differs from the previous time for which this was calculated.
    if (time /= self%timePrevious) then
       block
         type     (varying_string) :: fileNameVariance
         type     (lockDescriptor) :: fileLock
         type     (hdf5Object    ) :: varianceFile
         character(len=12        ) :: timeLabel       , lengthCubeLabel
         integer                   :: useCache        , i
         
         self%timePrevious=time
         write (lengthCubeLabel,'(e12.6)') self%lengthSimulationCube
         write (      timeLabel,'(e12.6)')      time
         fileNameVariance=inputPath(pathTypeDataDynamic)                                  // &
              &           'largeScaleStructure/'                                          // &
              &           self%objectType()                                               // &
              &           'CubeVariance_'                                                 // &
              &           'lengthSimulationCube_'//trim(lengthCubeLabel)//'_'             // &
              &           'time_'                //trim(      timeLabel)//'_'             // &
              &           self%powerSpectrum_%hashedDescriptor(includeSourceDigest=.true.)// &
              &           '.hdf5'
         !$omp critical(haloMassFunctionSimulationVarianceCache)
         useCache=0
         if (countCache > 0) then
            do i=1,countCache
               if (cachedVariances(i)%fileName == fileNameVariance) then
                  useCache=i
                  exit
               end if
            end do
         end if
         if (useCache /= 0) self%varianceSimulation=cachedVariances(useCache)%variance
         !$omp end critical(haloMassFunctionSimulationVarianceCache)
         if (useCache ==0) then
            if (File_Exists(fileNameVariance)) then
               call File_Lock(char(fileNameVariance),fileLock,lockIsShared=.true.)
               !$ call hdf5Access%set()
               call varianceFile%openFile     (char(fileNameVariance),readOnly=.true.                 )
               call varianceFile%readAttribute('variance'             ,        self%varianceSimulation)
               call varianceFile%close        (                                                       )
               !$ call hdf5Access%unset()
               call File_Unlock(fileLock)
            else
               call displayIndent('computing simulation variance',verbosityLevelWorking)
               call displayMessage('                time = '//timeLabel      //' Gyr',verbosityLevelWorking)
               call displayMessage('lengthSimulationCube = '//lengthCubeLabel//' Mpc',verbosityLevelWorking)
               powerSpectrum_         => self%powerSpectrum_
               time_                   =       time
               lengthSimulationCube_   =  self%lengthSimulationCube
               integratorX             =  integrator(integrandVarianceSimulationX,toleranceRelative=1.0d-6)
               integratorY             =  integrator(integrandVarianceSimulationY,toleranceRelative=1.0d-6)
               integratorZ             =  integrator(integrandVarianceSimulationZ,toleranceRelative=1.0d-6)
               wavenumberMinimum       =  +0.0d0
               wavenumberMaximum       =  +wavenumberMaximumFractional/self%lengthSimulationCube
               self%varianceSimulation =  integratorX%integrate(wavenumberMinimum,wavenumberMaximum)
               call File_Lock(char(fileNameVariance),fileLock,lockIsShared=.false.)
               !$ call hdf5Access%set()
               call varianceFile%openFile      (char(fileNameVariance) ,readOnly=.false.   ,overWrite=.true.)
               call varianceFile%writeAttribute(self%varianceSimulation,         'variance'                 )
               call varianceFile%close         (                                                            )
               !$ call hdf5Access%unset()
               call File_Unlock(fileLock)
               call displayUnindent('done',verbosityLevelWorking)
            end if
            !$omp critical(haloMassFunctionSimulationVarianceCache)
            lastCache=lastCache+1
            if (lastCache > sizeCache) lastCache=1
            countCache=max(countCache,lastCache)
            cachedVariances(lastCache)%fileName=fileNameVariance
            cachedVariances(lastCache)%variance=self%varianceSimulation
            !$omp end critical(haloMassFunctionSimulationVarianceCache)
         end if
      end block
    end if
    ! Modify the mass function by the perturbation.
    massFunction=+self%massFunction_%differential(time,mass,node)                        &
         &       *max(                                                                   &
         &            +0.0d0                                                           , &
         &            +1.0d0                                                             &
         &            +     self                    %perturbationFractional              &
         &            *     self%darkMatterHaloBias_%bias                  (mass,time)   &
         &            *sqrt(self                    %varianceSimulation               )  &
         &           )
    return
  end function simulationVarianceDifferential

  double precision function integrandVarianceSimulationX(wavenumber)
    !!{
    Integrand function for cosmic variance in a simulation cube: $x$-axis.
    !!}
    implicit none
    double precision, intent(in   ) :: wavenumber

    wavenumberX                 =wavenumber
    windowFunctionX             =windowFunctionCube(wavenumberX*lengthSimulationCube_)**2
    integrandVarianceSimulationX=integratorY%integrate(wavenumberMinimum,wavenumberMaximum)
    return
  end function integrandVarianceSimulationX
  
  double precision function integrandVarianceSimulationY(wavenumber)
    !!{
    Integrand function for cosmic variance in a simulation cube: $y$-axis.
    !!}
    implicit none
    double precision, intent(in   ) :: wavenumber

    wavenumberY                 =wavenumber
    windowFunctionY             =windowFunctionCube(wavenumberY*lengthSimulationCube_)**2
    integrandVarianceSimulationY=integratorZ%integrate(wavenumberMinimum,wavenumberMaximum)
    return
  end function integrandVarianceSimulationY
  
  double precision function integrandVarianceSimulationZ(wavenumber)
    !!{
    Integrand function for cosmic variance in a simulation cube: $z$-axis.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: wavenumber
    double precision                :: wavenumberMagnitude
    
    wavenumberZ                 =wavenumber
    windowFunctionZ             =+windowFunctionCube(wavenumberZ*lengthSimulationCube_)**2
    wavenumberMagnitude         =+sqrt(+wavenumberX**2+wavenumberY**2+wavenumberZ**2)    
    if (wavenumberMagnitude == 0.0d0) then
       integrandVarianceSimulationZ=+0.0d0
    else
       ! We integrate only over positive wavenumbers, so we include a factor 8 here to account for the other octants.
       integrandVarianceSimulationZ=+8.0d0                                           &
            &                       *powerSpectrum_%power(wavenumberMagnitude,time_) &
            &                       *windowFunctionX                                 &
            &                       *windowFunctionY                                 &
            &                       *windowFunctionZ
    end if
    return
  end function integrandVarianceSimulationZ

  double precision function windowFunctionCube(x)
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: x

    if (x == 0.0d0) then
       windowFunctionCube=+1.0d0               /sqrt(2.0d0*Pi)
    else
       windowFunctionCube=+2.0d0*sin(0.5d0*x)/x/sqrt(2.0d0*Pi)
    end if
    return
  end function windowFunctionCube
  
