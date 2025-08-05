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
  Implements a stellar tracks class in which the tracks are read from file and interpolated.
  !!}

  use :: Numerical_Interpolation, only : interpolator
  
  !![
  <stellarTracks name="stellarTracksFile">
   <description>
    A stellar tracks class in which luminosities and effective temperatures of stars are computed from a tabulated set of
    stellar tracks, read from file and interpolated. The file containing the tracks to use is specified via the {\normalfont
    \ttfamily stellarTracksFile} parameter. The file specified must be an HDF5 file with the following structure:
    \begin{verbatim}
     stellarTracksFile
      |
      +-> metallicity1
      |    |
      |    +-> metallicity
      |    |
      |    +-> mass1
      |    |    |
      |    |    +-> mass
      |    |    |
      |    |    +-> age
      |    |    |
      |    |    +-> luminosity
      |    |    |
      |    |    +-> effectiveTemperature
      |    |
      |    x-> massN
      |
      x-> metallicityN
    \end{verbatim}
    Each {\normalfont \ttfamily metallicityN} group tabulates tracks for a given metallicity (the value of which is stored in
    the {\normalfont \ttfamily metallicity} dataset within each group), and may contain an arbitrary number of {\normalfont
    \ttfamily massN} groups. Each {\normalfont \ttfamily massN} group should contain a track for a star of some mass (the value
    of which is given in the {\normalfont \ttfamily mass} dataset). Within each track three datasets specify the {\normalfont
    \ttfamily age} (in Gyr), {\normalfont \ttfamily luminosity} (in $L_\odot$) and {\normalfont \ttfamily effectiveTemperature}
    (in Kelvin) along the track.
   </description>
   <runTimeFileDependencies paths="fileName"/>
  </stellarTracks>
  !!]
  type, extends(stellarTracksClass) :: stellarTracksFile
     !!{
     A stellar tracks class in which the tracks are read from file and interpolated.
     !!}
     private
     type            (varying_string)                                :: fileName
     double precision                , allocatable, dimension(:    ) :: metallicityLogarithmic
     double precision                , allocatable, dimension(:,:  ) :: massInitial
     double precision                , allocatable, dimension(:,:,:) :: age                    , luminosityTrack, &
          &                                                             temperatureTrack
     integer                         , allocatable, dimension(:    ) :: countMassInitial
     integer                         , allocatable, dimension(:,:  ) :: countAge
     integer                                                         :: countMetallicity
     type            (interpolator  )                                :: interpolatorMetallicity
     type            (interpolator  ), allocatable, dimension(:,:  ) :: interpolatorAge
     type            (interpolator  ), allocatable, dimension(:    ) :: interpolatorMass
     logical                                                         :: initialized
   contains
     !![
     <methods>
       <method description="\textcolor{red}{\textless integer(c\_size\_t)(2)\textgreater} interpolationIndicesMetallicity\argout, {\textless integer(c\_size\_t)(2,2)\textgreater} interpolationIndicesMass\argout, {\textless integer(c\_size\_t)(2,2,2)\textgreater} interpolationIndicesAge\argout, {\textless double(2)\textgreater} interpolationFactorsMetallicity\argout, {\textless double(2,2)\textgreater} interpolationFactorsMass\argout, {\textless double(2,2,2)\textgreater} interpolationFactorsAge\argout, \logicalzero\ metallicityOutOfRange\argout, \logicalzero\ massOutOfRange\argout, \logicalzero\ ageOutOfRange\argout" method="interpolationCompute" />
       <method description="\textcolor{red}{\textless integer(c\_size\_t)(2)\textgreater} interpolationIndicesMetallicity\argin, {\textless integer(c\_size\_t)(2,2)\textgreater} interpolationIndicesMass\argin, {\textless integer(c\_size\_t)(2,2,2)\textgreater} interpolationIndicesAge\argin, {\textless double(2)\textgreater} interpolationFactorsMetallicity\argin, {\textless double(2,2)\textgreater} interpolationFactorsMass\argin, {\textless double(2,2,2)\textgreater} interpolationFactorsAge\argin, {\textless double(:,:,:)\textgreater} stellarTracks\argin" method="interpolate" />
       <method method="initialize" description="Initialize stellar data."/>
     </methods>
     !!]
     procedure :: luminosity           => fileLuminosity
     procedure :: temperatureEffective => fileTemperatureEffective
     procedure :: interpolationCompute => fileInterpolationCompute
     procedure :: interpolate          => fileInterpolate
     procedure :: initialize           => fileInitialize
  end type stellarTracksFile

  interface stellarTracksFile
     !!{
     Constructors for the \refClass{stellarTracksFile} stellar tracks class.
     !!}
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface stellarTracksFile

  ! The current file format version.
  integer, parameter :: fileFormatVersionCurrent=1

contains

  function fileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarTracksFile} stellar tracks class which takes a parameter list as input.
    !!}
    use :: Input_Paths     , only : inputPath     , pathTypeDataStatic
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(stellarTracksFile)                :: self
    type(inputParameters  ), intent(inout) :: parameters
    type(varying_string   )                :: fileName

    !![
    <inputParameter>
      <name>fileName</name>
      <defaultValue>inputPath(pathTypeDataStatic)//'stellarAstrophysics/Stellar_Tracks_Padova.hdf5'</defaultValue>
      <description>The name of the HDF5 file from which to read stellar tracks.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=stellarTracksFile(char(fileName))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fileConstructorParameters

  function fileConstructorInternal(fileName) result(self)
    !!{
    Internal constructor for the \refClass{stellarTracksFile} stellar tracks class.
    !!}
    implicit none
    type     (stellarTracksFile)                :: self
    character(len=*            ), intent(in   ) :: fileName
    !![
    <constructorAssign variables="fileName"/>
    !!]

    self%initialized=.false.
    return
  end function fileConstructorInternal

  subroutine fileInitialize(self)
    !!{
    Read data for the {\normalfont \ttfamily file} stellar tracks class.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : assignment(=), operator(//), varying_string
    use :: String_Handling   , only : operator(//)
    implicit none
    class(stellarTracksFile), intent(inout) :: self

    if (self%initialized) return
    block
      type     (varying_string) :: groupName
      integer                   :: ageCountMaximum        , fileFormatVersion      , &
           &                       initialMassCount       , initialMassCountMaximum, &
           &                       metallicityCountMaximum, metallicityCount
      type     (hdf5Object    ) :: ageDataset             , massGroup              , &
           &                       metallicityGroup       , stellarTracks
      logical                   :: foundMassGroup         , foundMetallicityGroup

      ! Open the HDF5 file.
      !$ call hdf5Access%set()
      stellarTracks=hdf5Object(char(self%fileName),readOnly=.true.)
      ! Check that this file has the correct format.
      call stellarTracks%readAttribute('fileFormat',fileFormatVersion,allowPseudoScalar=.true.)
      if (fileFormatVersion /= fileFormatVersionCurrent) call Error_Report('format of stellar tracks file is out of date'//{introspection:location})
      ! Count up number of metallicities present, the number of stellar masses tabulated and the number of ages tabulated.
      metallicityCountMaximum=0
      initialMassCountMaximum=0
      ageCountMaximum        =0
      ! Count metallicity groups.
      foundMetallicityGroup=.true.
      do while (foundMetallicityGroup)
         groupName="metallicity"
         groupName=groupName//(metallicityCountMaximum+1)
         foundMetallicityGroup=stellarTracks%hasGroup(char(groupName))
         if (foundMetallicityGroup) then
            metallicityCountMaximum=metallicityCountMaximum+1
            metallicityGroup=stellarTracks%openGroup(char(groupName))
            ! Find mass groups.
            foundMassGroup=.true.
            do while (foundMassGroup)
               groupName="mass"
               groupName=groupName//(initialMassCountMaximum+1)
               foundMassGroup=metallicityGroup%hasGroup(char(groupName))
               if (foundMassGroup) then
                  initialMassCountMaximum=initialMassCountMaximum+1
                  massGroup=metallicityGroup%openGroup(char(groupName))
                  ageDataset=massGroup%openDataset('age')
                  ageCountMaximum=max(ageCountMaximum,int(ageDataset%size(1)))
               end if
            end do
         end if
      end do
      ! Allocate storage space for data.
      allocate(self%metallicityLogarithmic(metallicityCountMaximum))
      allocate(self%countMassInitial      (metallicityCountMaximum))
      allocate(self%massInitial           (initialMassCountMaximum,metallicityCountMaximum))
      allocate(self%countAge              (initialMassCountMaximum,metallicityCountMaximum))
      allocate(self%age                   (ageCountMaximum,initialMassCountMaximum,metallicityCountMaximum))
      allocate(self%luminosityTrack       (ageCountMaximum,initialMassCountMaximum,metallicityCountMaximum))
      allocate(self%temperatureTrack      (ageCountMaximum,initialMassCountMaximum,metallicityCountMaximum))
      ! Read in all data.
      do metallicityCount=1,metallicityCountMaximum
         ! Open the metallicity group.
         groupName="metallicity"
         groupName=groupName//metallicityCount
         metallicityGroup=stellarTracks%openGroup(char(groupName))
         ! Get the metallicity.
         call metallicityGroup%readDatasetStatic('metallicity', self%metallicityLogarithmic(metallicityCount:metallicityCount))
         ! Count how many masses are tabulated at this metallicity.
         initialMassCount=0
         foundMassGroup=.true.
         do while (foundMassGroup)
            groupName="mass"
            groupName=groupName//(initialMassCount+1)
            foundMassGroup=metallicityGroup%hasGroup(char(groupName))
            if (foundMassGroup) initialMassCount=initialMassCount+1
         end do
         self%countMassInitial(metallicityCount)=initialMassCount
         ! Loop through all tabulated masses.
         do initialMassCount=1,self%countMassInitial(metallicityCount)
            ! Open the mass group.
            groupName="mass"
            groupName=groupName//initialMassCount
            massGroup=metallicityGroup%openGroup(char(groupName))
            ! Get initial mass.
            call massGroup%readDatasetStatic('mass',self%massInitial(initialMassCount:initialMassCount,metallicityCount))
            ! Read tracks.
            ageDataset=massGroup%openDataset('age')
            self%countAge(initialMassCount,metallicityCount)=int(ageDataset%size(1))
            call massGroup%readDatasetStatic('age'                 ,self%age             (1:self%countAge(initialMassCount,metallicityCount),initialMassCount,metallicityCount))
            call massGroup%readDatasetStatic('luminosity'          ,self%luminosityTrack (1:self%countAge(initialMassCount,metallicityCount),initialMassCount,metallicityCount))
            call massGroup%readDatasetStatic('effectiveTemperature',self%temperatureTrack(1:self%countAge(initialMassCount,metallicityCount),initialMassCount,metallicityCount))
         end do
      end do
      ! Convert metallicities to logarithmic scale.
      self%metallicityLogarithmic=log(self%metallicityLogarithmic)
      self%countMetallicity=metallicityCountMaximum
      !$ call hdf5Access%unset()
      ! Initialize interpolators.
      allocate(self%interpolatorMass(                        metallicityCountMaximum))
      allocate(self%interpolatorAge (initialMassCountMaximum,metallicityCountMaximum))
      self%interpolatorMetallicity=interpolator(self%metallicityLogarithmic)
      do metallicityCount=1,metallicityCountMaximum
         self%interpolatorMass(metallicityCount)=interpolator(self%massInitial(1:self%countMassInitial(metallicityCount),metallicityCount))
         do initialMassCount=1,self%countMassInitial(metallicityCount)
            self%interpolatorAge(initialMassCount,metallicityCount)=interpolator(self%age(1:self%countAge(initialMassCount,metallicityCount),initialMassCount,metallicityCount))
         end do
      end do
      self%initialized=.true.
    end block
    return
  end subroutine fileInitialize

  double precision function fileLuminosity(self,initialMass,metallicity,age)
    !!{
    Return the bolometric luminosity (in $L_\odot$) for a star of given {\normalfont \ttfamily initialMass}, {\normalfont \ttfamily metallicity} and {\normalfont \ttfamily age}.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (stellarTracksFile), intent(inout)    :: self
    double precision                   , intent(in   )    :: age                            , initialMass   , &
         &                                                   metallicity
    integer         (c_size_t         ), dimension(2,2,2) :: interpolationIndicesAge
    integer         (c_size_t         ), dimension(2,2  ) :: interpolationIndicesMass
    integer         (c_size_t         ), dimension(2    ) :: interpolationIndicesMetallicity
    double precision                   , dimension(2,2,2) :: interpolationFactorsAge
    double precision                   , dimension(2,2  ) :: interpolationFactorsMass
    double precision                    ,dimension(2    ) :: interpolationFactorsMetallicity
    logical                                               :: ageOutOfRange                  , massOutOfRange, &
         &                                                   metallicityOutOfRange

    ! Get the interpolating factors.
    call self%interpolationCompute(                                 &
         &                         initialMass                    , &
         &                         metallicity                    , &
         &                         age                            , &
         &                         interpolationIndicesMetallicity, &
         &                         interpolationIndicesMass       , &
         &                         interpolationIndicesAge        , &
         &                         interpolationFactorsMetallicity, &
         &                         interpolationFactorsMass       , &
         &                         interpolationFactorsAge        , &
         &                         metallicityOutOfRange          , &
         &                         massOutOfRange                 , &
         &                         ageOutOfRange                    &
         &                        )
    ! Do the interpolation.
    if (massOutOfRange .or. ageOutOfRange) then
       fileLuminosity=0.0d0
    else
       fileLuminosity=self%interpolate(                                 &
            &                          interpolationIndicesMetallicity, &
            &                          interpolationIndicesMass       , &
            &                          interpolationIndicesAge        , &
            &                          interpolationFactorsMetallicity, &
            &                          interpolationFactorsMass       , &
            &                          interpolationFactorsAge        , &
            &                          self%luminosityTrack             &
            &                         )
    end if
    return
  end function fileLuminosity

  double precision function fileTemperatureEffective(self,initialMass,metallicity,age)
    !!{
    Return the effective temperature (in Kelvin) for a star of given {\normalfont \ttfamily initialMass}, {\normalfont \ttfamily metallicity} and {\normalfont \ttfamily age}.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
     class           (stellarTracksFile), intent(inout)    :: self
    double precision                    , intent(in   )    :: age                            , initialMass   , &
         &                                                    metallicity
    integer         (c_size_t          ), dimension(2,2,2) :: interpolationIndicesAge
    integer         (c_size_t          ), dimension(2,2  ) :: interpolationIndicesMass
    integer         (c_size_t          ), dimension(2    ) :: interpolationIndicesMetallicity
    double precision                    , dimension(2,2,2) :: interpolationFactorsAge
    double precision                    , dimension(2,2  ) :: interpolationFactorsMass
    double precision                     ,dimension(2    ) :: interpolationFactorsMetallicity
    logical                                                :: ageOutOfRange                  , massOutOfRange, &
         &                                                    metallicityOutOfRange

    ! Get the interpolating factors.
    call self%interpolationCompute(                                 &
         &                         initialMass                    , &
         &                         metallicity                    , &
         &                         age                            , &
         &                         interpolationIndicesMetallicity, &
         &                         interpolationIndicesMass       , &
         &                         interpolationIndicesAge        , &
         &                         interpolationFactorsMetallicity, &
         &                         interpolationFactorsMass       , &
         &                         interpolationFactorsAge        , &
         &                         metallicityOutOfRange          , &
         &                         massOutOfRange                 , &
         &                         ageOutOfRange                    &
         &                        )
    ! Do the interpolation.
    if (massOutOfRange.or.ageOutOfRange) then
       fileTemperatureEffective=0.0d0
    else
       fileTemperatureEffective=self%interpolate(                                 &
            &                                    interpolationIndicesMetallicity, &
            &                                    interpolationIndicesMass       , &
            &                                    interpolationIndicesAge        , &
            &                                    interpolationFactorsMetallicity, &
            &                                    interpolationFactorsMass       , &
            &                                    interpolationFactorsAge        , &
            &                                    self%temperatureTrack            &
            &                                   )
    end if
    return
  end function fileTemperatureEffective

  double precision function fileInterpolate(self,interpolationIndicesMetallicity,interpolationIndicesMass,interpolationIndicesAge,interpolationFactorsMetallicity,interpolationFactorsMass,interpolationFactorsAge,stellarTracks)
    !!{
    Using precomputed factors, interpolate in metallicity, mass and age in the given {\normalfont \ttfamily stellarTracks}.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: Error        , only : Error_Report
    implicit none
    class           (stellarTracksFile), intent(inout)                   :: self
    integer         (c_size_t         ), intent(in   ), dimension(2,2,2) :: interpolationIndicesAge
    integer         (c_size_t         ), intent(in   ), dimension(2,2  ) :: interpolationIndicesMass
    integer         (c_size_t         ), intent(in   ), dimension(2    ) :: interpolationIndicesMetallicity
    double precision                   , intent(in   ), dimension(2,2,2) :: interpolationFactorsAge
    double precision                   , intent(in   ), dimension(2,2  ) :: interpolationFactorsMass
    double precision                   , intent(in   ), dimension(2    ) :: interpolationFactorsMetallicity
    double precision                   , intent(in   ), dimension(:,:,:) :: stellarTracks
    integer                                                              :: iAge                           , iMass       , &
         &                                                                  iMetallicity
    integer         (c_size_t         )                                  :: jAge                           , jMetallicity, &
         &                                                                  jMass
    !$GLC attributes unused :: self

    fileInterpolate=0.0d0
    do iMetallicity=1,2
       jMetallicity=interpolationIndicesMetallicity(iMetallicity)
       if       (jMetallicity < 1 .or. jMetallicity > size(stellarTracks,dim=3)                ) call Error_Report('metallicity index out of range'//{introspection:location})
       do iMass=1,2
          jMass=interpolationIndicesMass(iMetallicity,iMass)
          if    (jMass        < 1 .or. jMass        > self%countMassInitial(      jMetallicity)) call Error_Report('mass index out of range'       //{introspection:location})
          do iAge=1,2
             jAge=interpolationIndicesAge(iMetallicity,iMass,iAge)
             if (jAge         < 1 .or. jAge         > self%countAge        (jMass,jMetallicity)) call Error_Report('age index out of range'        //{introspection:location})
             fileInterpolate=+fileInterpolate                                          &
                  &          +stellarTracks                  (jAge,jMass,jMetallicity) &
                  &          *interpolationFactorsMetallicity(iMetallicity           ) &
                  &          *interpolationFactorsMass       (iMetallicity,iMass     ) &
                  &          *interpolationFactorsAge        (iMetallicity,iMass,iAge)
          end do
       end do
    end do
    return
  end function fileInterpolate

  subroutine fileInterpolationCompute(self,initialMass,metallicity,age,interpolationIndicesMetallicity,interpolationIndicesMass,interpolationIndicesAge,interpolationFactorsMetallicity,interpolationFactorsMass,interpolationFactorsAge,metallicityOutOfRange,massOutOfRange,ageOutOfRange)
    !!{
    Get interpolating factors for stellar tracks.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (stellarTracksFile), intent(inout)                   :: self
    double precision                   , intent(in   )                   :: age                            , initialMass   , &
         &                                                                  metallicity
    integer         (c_size_t         ), intent(  out), dimension(2,2,2) :: interpolationIndicesAge
    integer         (c_size_t         ), intent(  out), dimension(2,2  ) :: interpolationIndicesMass
    integer         (c_size_t         ), intent(  out), dimension(2    ) :: interpolationIndicesMetallicity
    double precision                   , intent(  out), dimension(2,2,2) :: interpolationFactorsAge
    double precision                   , intent(  out), dimension(2,2  ) :: interpolationFactorsMass
    double precision                   , intent(  out), dimension(2    ) :: interpolationFactorsMetallicity
    logical                            , intent(  out)                   :: ageOutOfRange                  , massOutOfRange, &
         &                                                                  metallicityOutOfRange
    integer                                                              :: iMass                          , iMetallicity
    integer         (c_size_t         )                                  :: jMass                          , jMetallicity
    double precision                                                     :: logMetallicity

    call self%initialize()
    ! Assume everything is in range initially.
    metallicityOutOfRange=.false.
    massOutOfRange       =.false.
    ageOutOfRange        =.false.
    ! Interpolate in metallicity.
    if (metallicity <= 0.0d0) then
       logMetallicity=self%metallicityLogarithmic(1)
    else
       logMetallicity=log(metallicity)
    end if
    if (logMetallicity < self%metallicityLogarithmic(1)) then
       interpolationIndicesMetallicity=[    1,    2]
       interpolationFactorsMetallicity=[1.0d0,0.0d0]
       metallicityOutOfRange=.true.
    else if (logMetallicity > self%metallicityLogarithmic(self%countMetallicity)) then
       interpolationIndicesMetallicity=[self%countMetallicity-1,self%countMetallicity]
       interpolationFactorsMetallicity=[0.0d0,1.0d0]
       metallicityOutOfRange=.true.
    else
       call self%interpolatorMetallicity%linearFactors(logMetallicity,interpolationIndicesMetallicity(1),interpolationFactorsMetallicity)
       interpolationIndicesMetallicity(2)=interpolationIndicesMetallicity(1)+1
    end if
    ! Loop over metallicities.
    do iMetallicity=1,2
       jMetallicity=interpolationIndicesMetallicity(iMetallicity)
       ! Interpolate in mass at each metallicity.
       if (initialMass < self%massInitial(1,jMetallicity)) then
          interpolationIndicesMass(iMetallicity,:)=[    1,    2]
          interpolationFactorsMass(iMetallicity,:)=[1.0d0,0.0d0]
          massOutOfRange=.true.
       else if (initialMass > self%massInitial(self%countMassInitial(jMetallicity),jMetallicity)) then
          interpolationIndicesMass(iMetallicity,:)=[self%countMassInitial(jMetallicity)-1,self%countMassInitial(jMetallicity)]
          interpolationFactorsMass(iMetallicity,:)=[0.0d0,1.0d0]
          massOutOfRange=.true.
       else
          call self%interpolatorMass(jMetallicity)%linearFactors(initialMass,interpolationIndicesMass(iMetallicity,1),interpolationFactorsMass(iMetallicity,:))
          interpolationIndicesMass(iMetallicity,2)=interpolationIndicesMass(iMetallicity,1)+1
       end if
       ! Loop over masses.
       do iMass=1,2
          jMass=interpolationIndicesMass(iMetallicity,iMass)
          ! Interpolate in age at each mass.
          if (age < self%age(1,jMass,jMetallicity)) then
             interpolationIndicesAge(iMetallicity,iMass,:)=[    1,    2]
             interpolationFactorsAge(iMetallicity,iMass,:)=[1.0d0,0.0d0]
             ageOutOfRange=.true.
          else if (age > self%age(self%countAge(jMass,jMetallicity),jMass,jMetallicity)) then
             interpolationIndicesAge(iMetallicity,iMass,:)=[self%countAge(jMass,jMetallicity)-1,self%countAge(jMass,jMetallicity)]
             interpolationFactorsAge(iMetallicity,iMass,:)=[0.0d0,1.0d0]
             ageOutOfRange=.true.
          else
             call self%interpolatorAge(jMass,jMetallicity)%linearFactors(age,interpolationIndicesAge(iMetallicity,iMass,1),interpolationFactorsAge(iMetallicity,iMass,:))
             interpolationIndicesAge(iMetallicity,iMass,2)=interpolationIndicesAge(iMetallicity,iMass,1)+1
          end if
       end do
    end do
    return
  end subroutine fileInterpolationCompute
