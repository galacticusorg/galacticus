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
Implements a property extractor class for broadband stellar luminosities computed directly from the star formation history.
!!}

  use :: Cosmology_Functions                       , only : cosmologyFunctionsClass
  use :: Galactic_Structure_Options                , only : enumerationComponentTypeType
  use :: ISO_Varying_String                        , only : varying_string
  use :: Numerical_Interpolation                   , only : interpolator
  use :: Star_Formation_Histories                  , only : starFormationHistoryClass
  use :: Stellar_Luminosities_Structure            , only : enumerationFrameType
  use :: Stellar_Population_Broad_Band_Luminosities, only : stellarPopulationBroadBandLuminositiesClass
  use :: Stellar_Population_Spectra_Postprocess    , only : stellarPopulationSpectraPostprocessorClass , &
  &                                                         stellarPopulationSpectraPostprocessorList
  use :: Stellar_Populations                       , only : stellarPopulationClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorLuminosityStellarFromStarFormationHistory" docformat="rst">
   <description>
   A property extractor that computes broadband stellar luminosities in AB units directly from the star formation history, by contracting it against a template of filter-convolved simple stellar population luminosities. For each filter named in ``filterNames`` the luminosity is

   .. math::

      L = \sum_{i,j} M_{ij} \bar{L}_{ij},

   where :math:`M_{ij}` is the mass formed in bin :math:`i` of age and bin :math:`j` of metallicity in the galaxy's star formation history, and :math:`\bar{L}_{ij}` is the mean luminosity per unit mass of a simple stellar population over that bin, convolved with the filter response. The cost per galaxy is therefore a single contraction over (age :math:`\times` metallicity) bins per filter, with no wavelength integration---in contrast to :galacticus-class:`nodePropertyExtractorLuminosityStellarFromSED`, which integrates a filter response over an SED. That makes this the efficient route when many galaxies and few filters dominate the cost, which is the usual case for lightcone output. It is also free of the resolution of any SED: the filter integral is performed when the simple stellar population luminosities are tabulated, not over an SED sampled at finite ``resolution``.

   The template depends on the bin structure of the star formation history, which must therefore be the same for every galaxy: the star formation history class must report fixed ages (see :galacticus-class:`starFormationHistoryFixedAges`). The template is built on first use---not at construction---because the bin structure is obtained from a galaxy's own star formation history.

   For ``frame`` = ``observed`` the simple stellar population luminosities depend on redshift, through both the shift of the filter relative to the spectrum and the redshift-dependent postprocessing of the spectrum (e.g. attenuation by the :term:`IGM`). Since a lightcone galaxy is observed at its own crossing redshift, and not at an output redshift, the template is tabulated on a grid spanning ``redshiftMinimum`` to ``redshiftMaximum``, spaced by ``redshiftSpacingLogarithmic`` in :math:`\log(1+z)`, and interpolated to the redshift of each galaxy by a local cubic in that same variable. The rest-frame position of a filter scales as :math:`1/(1+z)`, so equal steps in :math:`\log(1+z)` sweep the filter across the spectrum by equal fractions of wavelength, which is what sets how fast the tabulated luminosity varies: a single spacing in this variable is therefore the right criterion at every redshift. A galaxy falling outside the grid is reported as an error rather than extrapolated.

   The spacing which is needed depends on the filter, and is set by how fast the filter sweeps across structure in the spectrum rather than by the smoothness of the postprocessing alone. Measured against :galacticus-class:`nodePropertyExtractorLuminosityStellarFromSED`, which needs no grid and is exact, for galaxies at :math:`z \simeq 0.44`: a broad filter (``Galex_NUV``, ``SDSS_r``) agrees to :math:`3\times10^{-4}` at a spacing of :math:`\Delta z = 0.05`, and to :math:`\sim10^{-3}` even at :math:`\Delta z = 0.1`; but a narrow top-hat at :math:`1300\,\AA`, which at that redshift samples the spectrum across the Lyman break, is in error by :math:`2\times10^{-2}` at :math:`\Delta z = 0.05`, :math:`8\times10^{-3}` at :math:`0.025`, and :math:`10^{-3}` only by :math:`\Delta z = 0.01`. A filter which straddles a break at the redshifts of interest is therefore the case which sets the spacing, and the default is chosen to satisfy it. Where the two routes are both available, comparing them is the way to confirm the choice; and note that a residual difference of :math:`\sim10^{-3}` remains however fine the grid, since the from-SED route integrates the filter over an SED of finite ``resolution`` while this route uses the exact filter integral.

   For ``frame`` = ``rest`` no grid is needed and none is built: rest-frame luminosities are tabulated at zero redshift, following the convention of the broadband luminosity tabulation, which is passed the redshift to which the *filter* is shifted (zero for a rest-frame filter). Note that this differs from :galacticus-class:`nodePropertyExtractorLuminosityStellarFromSED`, which for a rest-frame luminosity applies the postprocessing at the redshift of the galaxy. The two routes therefore agree in the rest frame only for redshift-independent postprocessing.

   The template is the mean over each bin of the same bilinear (in age, and in :math:`\log_{10}` metallicity) interpolant of the tabulated simple stellar population luminosities that the broadband luminosity tabulation itself uses, and it is computed in closed form rather than by quadrature: on each cell of the tabulation the interpolant is linear in age and in :math:`\log_{10} Z`, so both integrals are elementary. There is consequently no integration tolerance here, and no possibility of an unmet one. Outside the tabulated range of age or of metallicity the interpolant is held constant at the edge of the tabulation, which matches the tabulation's own treatment of metallicity, and which for age is bounded where a linear extrapolation of a steeply rising young-population luminosity need not be.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorLuminosityStellarFromStarFormationHistory
     !!{RST
     A stellar luminosity property extractor class which works directly from the star formation history.
     !!}
     private
     type            (varying_string                             ), allocatable, dimension(:)       :: filterNames
     integer                                                      , allocatable, dimension(:)       :: filterIndex                                      , &
          &                                                                                            luminosityIndex
     double precision                                             , allocatable, dimension(:)       :: redshiftGrid                                     , &
          &                                                                                            redshiftTabulation                               , &
          &                                                                                            metallicityBoundaries
     type            (stellarPopulationSpectraPostprocessorList  ), allocatable, dimension(:)       :: postprocessor
     double precision                                             , allocatable, dimension(:,:,:,:) :: template
     type            (interpolator                               ), allocatable                     :: interpolatorRedshift
     class           (starFormationHistoryClass                  ), pointer                         :: starFormationHistory_                   => null()
     class           (stellarPopulationClass                     ), pointer                         :: stellarPopulation_                      => null()
     class           (stellarPopulationBroadBandLuminositiesClass), pointer                         :: stellarPopulationBroadBandLuminosities_ => null()
     class           (stellarPopulationSpectraPostprocessorClass ), pointer                         :: stellarPopulationSpectraPostprocessor_  => null()
     class           (cosmologyFunctionsClass                    ), pointer                         :: cosmologyFunctions_                     => null()
     type            (enumerationComponentTypeType               )                                  :: component
     type            (enumerationFrameType                       )                                  :: frame
     double precision                                                                               :: redshiftMinimum                                  , &
          &                                                                                            redshiftMaximum                                  , &
          &                                                                                            redshiftSpacingLogarithmic
     integer                                                                                        :: countRedshifts
     logical                                                                                        :: templateBuilt
   contains
     final              ::                  luminosityStellarFromSFHDestructor
     procedure          :: elementCount  => luminosityStellarFromSFHElementCount
     procedure          :: extract       => luminosityStellarFromSFHExtract
     procedure          :: quantity      => luminosityStellarFromSFHQuantity
     procedure          :: names         => luminosityStellarFromSFHNames
     procedure          :: descriptions  => luminosityStellarFromSFHDescriptions
     procedure          :: unitsInSI     => luminosityStellarFromSFHUnitsInSI
     procedure          :: units         => luminosityStellarFromSFHUnits
     procedure, private :: templateBuild => luminosityStellarFromSFHTemplateBuild
  end type nodePropertyExtractorLuminosityStellarFromStarFormationHistory

  interface nodePropertyExtractorLuminosityStellarFromStarFormationHistory
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorLuminosityStellarFromStarFormationHistory` property extractor class.
     !!}
     module procedure luminosityStellarFromSFHConstructorParameters
     module procedure luminosityStellarFromSFHConstructorInternal
  end interface nodePropertyExtractorLuminosityStellarFromStarFormationHistory

  ! Smallest separation in redshift which the broadband luminosity tabulation can distinguish. It names the dataset holding each
  ! tabulation `redshift<z>`, with `z` written in `f7.4` format, so two grid redshifts closer than this would silently share a
  ! tabulation.
  double precision, parameter :: redshiftSeparationMinimum=1.0d-4

  ! Counter from which luminosity indices are allocated. A luminosity index addresses one tabulation held by the
  ! `stellarPopulationBroadBandLuminosities` object, and that object records only whether an index has been tabulated---not the
  ! filter or the redshift which that tabulation is for. The object is shared with every other user of it which does not
  ! declare its own: the evolved luminosities of `Stellar_Luminosities_Structure`, which use indices 1 to
  ! `unitStellarLuminosities%luminosityCount()`, and any other object of this class. Two users given the same index would
  ! silently share one tabulation, so indices are allocated here from a single counter which starts above the block reserved
  ! for the evolved luminosities.
  integer                     :: luminosityIndexAllocated =-1

contains

  function luminosityStellarFromSFHConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorLuminosityStellarFromStarFormationHistory` property extractor class which takes a parameter set as input.
    !!}
    use :: Galactic_Structure_Options    , only : enumerationComponentTypeEncode
    use :: Input_Parameters              , only : inputParameter                , inputParameters
    use :: ISO_Varying_String            , only : char                          , var_str
    use :: Stellar_Luminosities_Structure, only : enumerationFrameEncode
    implicit none
    type            (nodePropertyExtractorLuminosityStellarFromStarFormationHistory)                              :: self
    type            (inputParameters                                               ), intent(inout)               :: parameters
    class           (starFormationHistoryClass                                     ), pointer                     :: starFormationHistory_
    class           (stellarPopulationClass                                        ), pointer                     :: stellarPopulation_
    class           (stellarPopulationBroadBandLuminositiesClass                   ), pointer                     :: stellarPopulationBroadBandLuminosities_
    class           (stellarPopulationSpectraPostprocessorClass                    ), pointer                     :: stellarPopulationSpectraPostprocessor_
    class           (cosmologyFunctionsClass                                       ), pointer                     :: cosmologyFunctions_
    type            (varying_string                                                ), allocatable  , dimension(:) :: filterNames
    type            (varying_string                                                )                              :: component                              , &
         &                                                                                                           frame
    double precision                                                                                              :: redshiftMinimum                        , &
         &                                                                                                           redshiftMaximum                        , &
         &                                                                                                           redshiftSpacingLogarithmic

    allocate(filterNames(parameters%count('filterNames')))
    !![
    <inputParameter docformat="rst">
      <name>filterNames</name>
      <source>parameters</source>
      <description>
      The filters for which luminosities are to be computed.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>component</name>
      <source>parameters</source>
      <defaultValue>var_str('disk')</defaultValue>
      <description>
      The component (``disk``, ``spheroid``, or ``nuclearStarCluster``) whose star formation history is to be used.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>frame</name>
      <source>parameters</source>
      <defaultValue>var_str('rest')</defaultValue>
      <description>
      The frame (``rest`` or ``observed``) for which to compute luminosities.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>redshiftMinimum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The minimum redshift of the grid on which observed-frame luminosities are tabulated. Ignored for rest-frame luminosities, which need no grid.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>redshiftMaximum</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The maximum redshift of the grid on which observed-frame luminosities are tabulated. Ignored for rest-frame luminosities, which need no grid.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>redshiftSpacingLogarithmic</name>
      <source>parameters</source>
      <defaultValue>5.0d-3</defaultValue>
      <description>
      The target spacing, in :math:`\log(1+z)`, of the grid on which observed-frame luminosities are tabulated. The number of grid points follows from this and the range, and is at least four. Ignored for rest-frame luminosities, which need no grid.
      </description>
    </inputParameter>
    <objectBuilder class="starFormationHistory"                   name="starFormationHistory_"                   source="parameters"/>
    <objectBuilder class="stellarPopulation"                      name="stellarPopulation_"                      source="parameters"/>
    <objectBuilder class="stellarPopulationBroadBandLuminosities" name="stellarPopulationBroadBandLuminosities_" source="parameters"/>
    <objectBuilder class="stellarPopulationSpectraPostprocessor"  name="stellarPopulationSpectraPostprocessor_"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"                     name="cosmologyFunctions_"                     source="parameters"/>
    !!]
    self=nodePropertyExtractorLuminosityStellarFromStarFormationHistory(                                                                              &
         &                                                              filterNames                                                                 , &
         &                                                              enumerationComponentTypeEncode(char(component),includesPrefix=.false.)       , &
         &                                                              enumerationFrameEncode        (char(frame    ),includesPrefix=.false.)       , &
         &                                                              redshiftMinimum                                                             , &
         &                                                              redshiftMaximum                                                             , &
         &                                                              redshiftSpacingLogarithmic                                                  , &
         &                                                              starFormationHistory_                                                       , &
         &                                                              stellarPopulation_                                                          , &
         &                                                              stellarPopulationBroadBandLuminosities_                                     , &
         &                                                              stellarPopulationSpectraPostprocessor_                                      , &
         &                                                              cosmologyFunctions_                                                           &
         &                                                             )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationHistory_"                  />
    <objectDestructor name="stellarPopulation_"                     />
    <objectDestructor name="stellarPopulationBroadBandLuminosities_"/>
    <objectDestructor name="stellarPopulationSpectraPostprocessor_" />
    <objectDestructor name="cosmologyFunctions_"                    />
    !!]
    return
  end function luminosityStellarFromSFHConstructorParameters

  function luminosityStellarFromSFHConstructorInternal(filterNames,component,frame,redshiftMinimum,redshiftMaximum,redshiftSpacingLogarithmic,starFormationHistory_,stellarPopulation_,stellarPopulationBroadBandLuminosities_,stellarPopulationSpectraPostprocessor_,cosmologyFunctions_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorLuminosityStellarFromStarFormationHistory` property extractor class.
    !!}
    use :: Error                         , only : Error_Report
    use :: Galactic_Structure_Options    , only : componentTypeDisk            , componentTypeNuclearStarCluster, &
    &                                             componentTypeSpheroid        , enumerationComponentTypeDecode
    use :: Instruments_Filters           , only : Filter_Get_Index
    use :: ISO_Varying_String            , only : operator(//)                 , var_str
    use :: Star_Formation_Histories      , only : starFormationHistoryAgesFixed
    use :: Stellar_Luminosities_Structure, only : enumerationFrameDecode       , frameObserved                  , &
    &                                             frameRest
    use :: String_Handling               , only : operator(//)
    implicit none
    type            (nodePropertyExtractorLuminosityStellarFromStarFormationHistory)                              :: self
    type            (varying_string                                                ), intent(in   ), dimension(:) :: filterNames
    type            (enumerationComponentTypeType                                  ), intent(in   )               :: component
    type            (enumerationFrameType                                          ), intent(in   )               :: frame
    double precision                                                                , intent(in   )               :: redshiftMinimum                        , &
         &                                                                                                           redshiftMaximum                        , &
         &                                                                                                           redshiftSpacingLogarithmic
    class           (starFormationHistoryClass                                     ), intent(in   ), target       :: starFormationHistory_
    class           (stellarPopulationClass                                        ), intent(in   ), target       :: stellarPopulation_
    class           (stellarPopulationBroadBandLuminositiesClass                   ), intent(in   ), target       :: stellarPopulationBroadBandLuminosities_
    class           (stellarPopulationSpectraPostprocessorClass                    ), intent(in   ), target       :: stellarPopulationSpectraPostprocessor_
    class           (cosmologyFunctionsClass                                       ), intent(in   ), target       :: cosmologyFunctions_
    integer                                                                                                       :: i                                      , &
         &                                                                                                           j                                      , &
         &                                                                                                           k
    !![
    <constructorAssign variables="filterNames, component, frame, redshiftMinimum, redshiftMaximum, redshiftSpacingLogarithmic, *starFormationHistory_, *stellarPopulation_, *stellarPopulationBroadBandLuminosities_, *stellarPopulationSpectraPostprocessor_, *cosmologyFunctions_"/>
    !!]

    ! The template is built on first use, from a galaxy's own star formation history: the bin structure is not available until
    ! then.
    self%templateBuilt=.false.
    ! Only the components which carry a star formation history can be used.
    if     (                                                   &
         &   self%component /= componentTypeDisk               &
         &  .and.                                              &
         &   self%component /= componentTypeSpheroid           &
         &  .and.                                              &
         &   self%component /= componentTypeNuclearStarCluster &
         & ) call Error_Report(                                                                                        &
         &                     'component "'                                                                        // &
         &                     enumerationComponentTypeDecode(self%component,includePrefix=.false.)                 // &
         &                     '" carries no star formation history: use disk, spheroid, or nuclearStarCluster'      // &
         &                     {introspection:location}                                                                &
         &                    )
    ! The template is shared by every galaxy, so it is usable only if every galaxy's star formation history has the same bin
    ! structure, which is what fixed ages mean.
    if (self%starFormationHistory_%ageDistribution() /= starFormationHistoryAgesFixed)                                      &
         & call Error_Report(                                                                                               &
         &                   'the star formation history must use fixed ages, so that one template serves every galaxy: ' // &
         &                   'use `starFormationHistory` of class `fixedAges`'                                            // &
         &                   {introspection:location}                                                                       &
         &                  )
    ! The tabulation is passed a single redshift, which serves both the shift of the filter relative to the spectrum and the
    ! postprocessing multiplier, so a multiplier which couples age and redshift can not be represented by it.
    if     (                                                                   &
         &   self%frame == frameObserved                                       &
         &  .and.                                                              &
         &   self%stellarPopulationSpectraPostprocessor_%isAgeDependent     () &
         &  .and.                                                              &
         &   self%stellarPopulationSpectraPostprocessor_%isRedshiftDependent() &
         & ) call Error_Report(                                                                                               &
         &                     'the spectra postprocessor depends on both age and redshift, which the broadband tabulation '// &
         &                     'can not represent: use `nodePropertyExtractorLuminosityStellarFromSED` for this case'       // &
         &                     {introspection:location}                                                                       &
         &                    )
    ! Build the redshift grid. Rest-frame luminosities are tabulated at zero redshift---the tabulation is passed the redshift to
    ! which the filter is shifted, which is zero for a rest-frame filter---so a single grid point serves every galaxy.
    select case (self%frame%ID)
    case (frameRest    %ID)
       self%countRedshifts=1
       allocate(self%redshiftGrid(1))
       self%redshiftGrid=0.0d0
    case (frameObserved%ID)
       if (self%redshiftSpacingLogarithmic <= 0.0d0)                                                      &
            & call Error_Report(                                                                          &
            &                   'redshiftSpacingLogarithmic must be positive'                           // &
            &                   {introspection:location}                                                  &
            &                  )
       if (self%redshiftMaximum <= self%redshiftMinimum)                                                  &
            & call Error_Report(                                                                          &
            &                   'redshiftMaximum must exceed redshiftMinimum'                           // &
            &                   {introspection:location}                                                  &
            &                  )
       ! The number of grid points follows from the requested spacing. At least four are used, so that the interpolation in
       ! redshift is cubic.
       self%countRedshifts=max(                                                                        &
            &                  4                                                                     , &
            &                  1+int(                                                                 &
            &                        ceiling(                                                         &
            &                                 log((1.0d0+self%redshiftMaximum)/(1.0d0+self%redshiftMinimum)) &
            &                                /self%redshiftSpacingLogarithmic                          &
            &                               )                                                          &
            &                       )                                                                  &
            &                 )
       allocate(self%redshiftGrid(self%countRedshifts))
       ! Space the grid logarithmically in 1+z. The rest-frame position of the filter scales as 1/(1+z), so equal steps in
       ! log(1+z) sweep the filter across the spectrum by equal fractions of wavelength---which is what sets how fast the
       ! tabulated luminosity varies, and so makes a single spacing in this variable the right criterion at every redshift.
       do i=1,self%countRedshifts
          self%redshiftGrid(i)=+exp(                                                      &
               &                    +log(1.0d0+self%redshiftMinimum)                       &
               &                    +log(                                                  &
               &                         +(1.0d0+self%redshiftMaximum)                     &
               &                         /(1.0d0+self%redshiftMinimum)                     &
               &                        )                                                  &
               &                    *dble(i-1)                                             &
               &                    /dble(self%countRedshifts-1)                           &
               &                   )                                                       &
               &               -1.0d0
       end do
       ! The tabulation names the dataset holding each redshift in `f7.4` format, so grid points closer than this would silently
       ! share a tabulation.
       do i=2,self%countRedshifts
          if (self%redshiftGrid(i)-self%redshiftGrid(i-1) < redshiftSeparationMinimum)                                       &
               & call Error_Report(                                                                                          &
               &                   var_str('the redshift grid is too fine: the broadband luminosity tabulation resolves ')// &
               &                   'redshifts only to 1.0e-4, so these grid points would share a tabulation. Increase '   // &
               &                   'redshiftSpacingLogarithmic, or widen the range'                                       // &
               &                   {introspection:location}                                                                  &
               &                  )
       end do
    case default
       call Error_Report('unknown frame'//{introspection:location})
    end select
    ! Build the private index space for the tabulation: one entry per (filter, grid redshift). These indices address this
    ! object's own tabulation, and are unrelated to the indices of the evolved luminosities registered in
    ! `Stellar_Luminosities_Structure`, so no luminosities need be registered there and none are evolved.
    allocate(self%luminosityIndex  (size(filterNames)*self%countRedshifts))
    allocate(self%filterIndex      (size(filterNames)*self%countRedshifts))
    allocate(self%redshiftTabulation(size(filterNames)*self%countRedshifts))
    allocate(self%postprocessor    (size(filterNames)*self%countRedshifts))
    k=0
    do i=1,self%countRedshifts
       do j=1,size(filterNames)
          k                          =  k+1
          self%filterIndex       (k) =  Filter_Get_Index(filterNames(j))
          self%redshiftTabulation(k) =  self%redshiftGrid(i)
          self%postprocessor     (k)%stellarPopulationSpectraPostprocessor_ => self%stellarPopulationSpectraPostprocessor_
       end do
    end do
    ! The indices themselves are allocated when the template is built: the block reserved for the evolved luminosities is not
    ! known until the stellar luminosities registry has been initialized, which happens after property extractors are built.
    self%luminosityIndex=-1
    return
  end function luminosityStellarFromSFHConstructorInternal

  subroutine luminosityStellarFromSFHDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodePropertyExtractorLuminosityStellarFromStarFormationHistory` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorLuminosityStellarFromStarFormationHistory), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationHistory_"                  />
    <objectDestructor name="self%stellarPopulation_"                     />
    <objectDestructor name="self%stellarPopulationBroadBandLuminosities_"/>
    <objectDestructor name="self%stellarPopulationSpectraPostprocessor_" />
    <objectDestructor name="self%cosmologyFunctions_"                    />
    !!]
    return
  end subroutine luminosityStellarFromSFHDestructor

  integer function luminosityStellarFromSFHElementCount(self,time)
    !!{RST
    Return the number of elements in the ``luminosityStellarFromStarFormationHistory`` property extractor.
    !!}
    implicit none
    class           (nodePropertyExtractorLuminosityStellarFromStarFormationHistory), intent(inout) :: self
    double precision                                                                , intent(in   ) :: time
    !$GLC attributes unused :: time

    luminosityStellarFromSFHElementCount=size(self%filterNames)
    return
  end function luminosityStellarFromSFHElementCount

  subroutine luminosityStellarFromSFHTemplateBuild(self,node,starFormationHistory,time,countAgeBins,countMetallicityBins)
    !!{RST
    Build the template of mean filter-convolved luminosities per unit mass, over each bin of age and metallicity of the star
    formation history, for each filter and each grid redshift.

    The template is the mean over each bin of the bilinear interpolant---linear in age, and linear in :math:`\log_{10}` of
    metallicity---of the tabulated simple stellar population luminosities, which is the interpolation that the broadband
    luminosity tabulation itself performs. On each cell of that tabulation the interpolant is linear in each variable, so the
    mean over a bin is elementary and is evaluated in closed form: there is no quadrature here, and so no tolerance to be met.
    !!}
    use :: Abundances_Structure            , only : abundances                   , metallicityTypeLinearByMassSolar
    use :: Error                           , only : Error_Report
    use :: Histories                       , only : history
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    use :: Stellar_Luminosities_Structure  , only : unitStellarLuminosities
    use :: Stellar_Population_Spectra      , only : stellarPopulationSpectraClass
    use :: Table_Labels                    , only : extrapolationTypeAbort
    implicit none
    class           (nodePropertyExtractorLuminosityStellarFromStarFormationHistory), intent(inout)                   :: self
    type            (treeNode                                                      ), intent(inout)                   :: node
    type            (history                                                       ), intent(in   )                   :: starFormationHistory
    double precision                                                                , intent(in   )                   :: time
    integer                                                                         , intent(in   )                   :: countAgeBins             , &
         &                                                                                                               countMetallicityBins
    double precision                                                                , allocatable  , dimension(:)     :: times                    , &
         &                                                                                                               ages                     , &
         &                                                                                                               ageTable                 , &
         &                                                                                                               metallicityTable         , &
         &                                                                                                               metallicityTableLog
    double precision                                                                , allocatable  , dimension(:,:)   :: luminositiesNode
    double precision                                                                , allocatable  , dimension(:,:,:) :: tableLuminosity
    type            (abundances                                                    )                                  :: abundancesStellar
    integer                                                                                                           :: countAgesTable           , &
         &                                                                                                               countMetallicitiesTable  , &
         &                                                                                                               iAgeBin                  , &
         &                                                                                                               iMetallicityBin          , &
         &                                                                                                               iFilter                  , &
         &                                                                                                               iRedshift                , &
         &                                                                                                               iTable                   , &
         &                                                                                                               iIndex
    double precision                                                                                                  :: timeStart                , &
         &                                                                                                               ageLow                   , &
         &                                                                                                               ageHigh                  , &
         &                                                                                                               metallicityLow           , &
         &                                                                                                               metallicityHigh          , &
         &                                                                                                               metallicityPopulationLow , &
         &                                                                                                               metallicityPopulationHigh

    ! Allocate the luminosity indices, now that the stellar luminosities registry is initialized and the block it reserves is
    ! therefore known.
    !$omp critical(luminosityStellarFromSFHIndexAllocate)
    if (luminosityIndexAllocated < 0)                                             &
         & luminosityIndexAllocated=max(                                          &
         &                              unitStellarLuminosities%luminosityCount(               ), &
         &                              unitStellarLuminosities%luminosityCount(unmapped=.true.)  &
         &                             )
    do iIndex=1,size(self%luminosityIndex)
       self%luminosityIndex(iIndex)=luminosityIndexAllocated+iIndex
    end do
    luminosityIndexAllocated=luminosityIndexAllocated+size(self%luminosityIndex)
    !$omp end critical(luminosityStellarFromSFHIndexAllocate)
    ! Get the bin structure of the star formation history. The ages are fixed---the constructor refuses any other case---so this
    ! structure is the same for every galaxy, and the template built here serves all of them.
    times                     =self%starFormationHistory_%times                (node=node,starFormationHistory=starFormationHistory,allowTruncation=.false.,timeStart=timeStart)
    self%metallicityBoundaries=self%starFormationHistory_%metallicityBoundaries()
    ! The bin counts are those of the mass array which will be contracted against the template, so that the two cannot disagree
    ! about the bin structure.
    if (size(times) /= countAgeBins)                                                                                        &
         & call Error_Report(                                                                                               &
         &                   'the star formation history reports a different number of age bins from its mass array'      // &
         &                   {introspection:location}                                                                       &
         &                  )
    ! Get the tabulation of the stellar population spectra: the grid on which the broadband luminosities are tabulated, and so
    ! the grid of the interpolant whose bin means are wanted here.
    block
      class(stellarPopulationSpectraClass), pointer :: stellarPopulationSpectra_

      stellarPopulationSpectra_ => self%stellarPopulation_%spectra()
      call stellarPopulationSpectra_%tabulation(countAgesTable,countMetallicitiesTable,ageTable,metallicityTable)
    end block
    allocate(metallicityTableLog(countMetallicitiesTable))
    metallicityTableLog     =log10(metallicityTable/metallicitySolar)
    metallicityPopulationLow =     metallicityTable(                      1)/metallicitySolar
    metallicityPopulationHigh=     metallicityTable(countMetallicitiesTable)/metallicitySolar
    ! Collect the tabulated luminosities at each tabulated metallicity. At a tabulated metallicity the tabulation's own
    ! interpolation returns the tabulated value exactly, so these calls recover the table itself, for every filter and grid
    ! redshift at once.
    allocate(tableLuminosity(size(self%luminosityIndex),countAgesTable,countMetallicitiesTable))
    do iTable=1,countMetallicitiesTable
       call abundancesStellar%metallicitySet(metallicityTable(iTable)/metallicitySolar,metallicityType=metallicityTypeLinearByMassSolar)
       call self%stellarPopulationBroadBandLuminosities_%luminosityTracks(      &
            &                                                            self%luminosityIndex   , &
            &                                                            self%filterIndex       , &
            &                                                            self%postprocessor     , &
            &                                                            self%stellarPopulation_, &
            &                                                            abundancesStellar      , &
            &                                                            self%redshiftTabulation, &
            &                                                            ages                   , &
            &                                                            luminositiesNode         &
            &                                                           )
       if (size(ages) /= countAgesTable) call Error_Report('the luminosity tabulation is not on the spectra tabulation grid'//{introspection:location})
       do iIndex=1,size(self%luminosityIndex)
          tableLuminosity(iIndex,:,iTable)=luminositiesNode(:,iIndex)
       end do
    end do
    ! Form the bin means.
    allocate(self%template(size(self%filterNames),countAgeBins,countMetallicityBins,self%countRedshifts))
    self%template=0.0d0
    do iRedshift=1,self%countRedshifts
       do iFilter=1,size(self%filterNames)
          iIndex=(iRedshift-1)*size(self%filterNames)+iFilter
          do iAgeBin=1,countAgeBins
             ! The bin spans times [timeStart or times(iAgeBin-1), times(iAgeBin)], and so ages measured back from the time of
             ! observation. Fixed ages make these the same for every galaxy.
             if (iAgeBin == 1) then
                ageHigh=time-     timeStart
             else
                ageHigh=time-     times(iAgeBin-1)
             end if
             ageLow    =time-min(times(iAgeBin  ),time)
             if (ageHigh <= ageLow) cycle
             do iMetallicityBin=1,countMetallicityBins
                ! Determine the metallicity range of the bin, clamped to the tabulated range, following the convention of
                ! `nodePropertyExtractorSED`.
                if (iMetallicityBin == 1) then
                   metallicityLow =                                                                   metallicityPopulationLow
                else
                   metallicityLow =min(max(self%metallicityBoundaries(iMetallicityBin-1),metallicityPopulationLow),metallicityPopulationHigh)
                end if
                metallicityHigh   =min(max(self%metallicityBoundaries(iMetallicityBin  ),metallicityPopulationLow),metallicityPopulationHigh)
                ! A bin which lies wholly outside the tabulated range of metallicity clamps to zero width. Its mass is not
                ! discarded: the luminosity is taken at the edge of the tabulation, which is what `nodePropertyExtractorSED`
                ! does for such a bin, and is the same treatment of metallicity which the tabulation's own lookup applies.
                self%template(iFilter,iAgeBin,iMetallicityBin,iRedshift)=                                       &
                     & +luminosityStellarFromSFHBinMean(                                                                                &
                     &          tableLuminosity    (iIndex,:,:),                                                &
                     &          ageTable                       ,                                                &
                     &          metallicityTable  /metallicitySolar,                                            &
                     &          metallicityTableLog            ,                                                &
                     &          ageLow                         ,                                                &
                     &          ageHigh                        ,                                                &
                     &          metallicityLow                 ,                                                &
                     &          metallicityHigh                                                                 &
                     &         )
             end do
          end do
       end do
    end do
    ! Build the interpolator used to interpolate the template to the redshift of a galaxy. It is held by value, so that it is
    ! destroyed with this object and a deep copy of this object is given an interpolator of its own.
    if (self%countRedshifts > 1) then
       allocate(self%interpolatorRedshift)
       self%interpolatorRedshift=interpolator(log(1.0d0+self%redshiftGrid),extrapolationType=extrapolationTypeAbort)
    end if
    self%templateBuilt=.true.
    return
  end subroutine luminosityStellarFromSFHTemplateBuild

  double precision function luminosityStellarFromSFHBinMean(table,ageTable,metallicityTable,metallicityTableLog,ageLow,ageHigh,metallicityLow,metallicityHigh)
    !!{RST
    Return the mean, over a bin of age and of metallicity, of the bilinear interpolant of a tabulated luminosity---linear in
    age, and linear in :math:`\log_{10}` of metallicity, which is the interpolation performed by the broadband luminosity
    tabulation. The integral is evaluated in closed form, cell by cell of the tabulation, so no quadrature is involved. Outside
    the tabulated range the interpolant is held constant at the edge of the tabulation.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:,:) :: table
    double precision, intent(in   ), dimension(:)   :: ageTable            , metallicityTable   , &
         &                                             metallicityTableLog
    double precision, intent(in   )                 :: ageLow              , ageHigh            , &
         &                                             metallicityLow      , metallicityHigh
    integer                                         :: iAge                , iMetallicity
    double precision                                :: ageCellLow          , ageCellHigh        , &
         &                                             metallicityCellLow  , metallicityCellHigh, &
         &                                             weightAgeLow        , weightAgeHigh      , &
         &                                             integralLow         , integralHigh
    logical                                         :: metallicityCollapsed

    ! A bin which lies wholly outside the tabulated range of metallicity has clamped to zero width. The luminosity is then
    ! taken at that edge of the tabulation and averaged over age alone, rather than over a range of zero width.
    metallicityCollapsed           =metallicityHigh <= metallicityLow
    luminosityStellarFromSFHBinMean=0.0d0
    ! Sum over the cells of the tabulation which overlap the bin, in both age and metallicity. Each cell contributes the exact
    ! integral of the interpolant over the overlap.
    do iAge=0,size(ageTable)
       ! Cells `0` and `size(ageTable)` are the regions below and above the tabulation, where the interpolant is held constant.
       if      (iAge == 0              ) then
          ageCellLow =-huge(0.0d0)
          ageCellHigh=ageTable(1)
       else if (iAge == size(ageTable) ) then
          ageCellLow = ageTable(size(ageTable))
          ageCellHigh= huge(0.0d0)
       else
          ageCellLow = ageTable(iAge  )
          ageCellHigh= ageTable(iAge+1)
       end if
       if (min(ageCellHigh,ageHigh) <= max(ageCellLow,ageLow)) cycle
       call luminosityStellarFromSFHAgeWeights(ageTable,iAge,max(ageCellLow,ageLow),min(ageCellHigh,ageHigh),weightAgeLow,weightAgeHigh)
       if (metallicityCollapsed) then
          luminosityStellarFromSFHBinMean=+luminosityStellarFromSFHBinMean                                                                                      &
               &                          +weightAgeLow *luminosityStellarFromSFHLuminosity(table,metallicityTable,metallicityTableLog,iAge,.false.,metallicityLow) &
               &                          +weightAgeHigh*luminosityStellarFromSFHLuminosity(table,metallicityTable,metallicityTableLog,iAge,.true. ,metallicityLow)
          cycle
       end if
       do iMetallicity=0,size(metallicityTable)
          if      (iMetallicity == 0                      ) then
             metallicityCellLow =-huge(0.0d0)
             metallicityCellHigh= metallicityTable(1)
          else if (iMetallicity == size(metallicityTable) ) then
             metallicityCellLow = metallicityTable(size(metallicityTable))
             metallicityCellHigh= huge(0.0d0)
          else
             metallicityCellLow = metallicityTable(iMetallicity  )
             metallicityCellHigh= metallicityTable(iMetallicity+1)
          end if
          if (min(metallicityCellHigh,metallicityHigh) <= max(metallicityCellLow,metallicityLow)) cycle
          ! Integrate over metallicity at each of the two bounding ages of the cell, then combine with the age weights: the
          ! interpolant is a product of a function linear in age and one linear in log-metallicity, so the double integral
          ! factorizes cell by cell.
          integralLow =luminosityStellarFromSFHMetallicityIntegral(table,metallicityTable,metallicityTableLog,iAge,iMetallicity,.false.,max(metallicityCellLow,metallicityLow),min(metallicityCellHigh,metallicityHigh))
          integralHigh=luminosityStellarFromSFHMetallicityIntegral(table,metallicityTable,metallicityTableLog,iAge,iMetallicity,.true. ,max(metallicityCellLow,metallicityLow),min(metallicityCellHigh,metallicityHigh))
          luminosityStellarFromSFHBinMean=+luminosityStellarFromSFHBinMean                   &
               &  +weightAgeLow *integralLow  &
               &  +weightAgeHigh*integralHigh
       end do
    end do
    ! Normalize to a mean.
    luminosityStellarFromSFHBinMean=+luminosityStellarFromSFHBinMean &
         &                          /(ageHigh-ageLow)
    if (.not.metallicityCollapsed)                                      &
         & luminosityStellarFromSFHBinMean=+luminosityStellarFromSFHBinMean   &
         &                                 /(metallicityHigh-metallicityLow)
    return
  end function luminosityStellarFromSFHBinMean

  subroutine luminosityStellarFromSFHAgeWeights(ageTable,iAge,ageLow,ageHigh,weightAgeLow,weightAgeHigh)
    !!{RST
    Return the integrals, over the age interval [``ageLow``,``ageHigh``], of the two linear basis functions of the age
    interpolation in cell ``iAge`` of the tabulation---that is, the weights with which the luminosities at the two bounding ages
    of the cell enter the integral of the interpolant. Outside the tabulation the interpolant is held constant, so all of the
    weight falls on the single bounding age.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: ageTable
    integer         , intent(in   )               :: iAge
    double precision, intent(in   )               :: ageLow      , ageHigh
    double precision, intent(  out)               :: weightAgeLow, weightAgeHigh
    double precision                              :: ageCellLow  , ageCellHigh

    if      (iAge == 0             ) then
       weightAgeLow =ageHigh-ageLow
       weightAgeHigh=0.0d0
    else if (iAge == size(ageTable)) then
       weightAgeLow =0.0d0
       weightAgeHigh=ageHigh-ageLow
    else
       ageCellLow   =ageTable(iAge  )
       ageCellHigh  =ageTable(iAge+1)
       ! The integral over [ageLow,ageHigh] of (ageCellHigh-age)/(ageCellHigh-ageCellLow), and of its complement.
       weightAgeLow =+(                                                         &
            &          +(ageCellHigh-ageLow )**2                                &
            &          -(ageCellHigh-ageHigh)**2                                &
            &         )                                                         &
            &        /2.0d0                                                     &
            &        /(ageCellHigh-ageCellLow)
       weightAgeHigh=+(ageHigh-ageLow)-weightAgeLow
    end if
    return
  end subroutine luminosityStellarFromSFHAgeWeights

  double precision function luminosityStellarFromSFHMetallicityIntegral(table,metallicityTable,metallicityTableLog,iAge,iMetallicity,useHighAge,metallicityLow,metallicityHigh)
    !!{RST
    Return the integral, over the metallicity interval [``metallicityLow``,``metallicityHigh``], of the luminosity at one
    bounding age of cell ``iAge`` of the tabulation, interpolated linearly in :math:`\log_{10}` of metallicity across cell
    ``iMetallicity``. Outside the tabulated range of metallicity the luminosity is held constant, as the tabulation's own
    lookup does, which also makes zero metallicity well defined.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:,:) :: table
    double precision, intent(in   ), dimension(:)   :: metallicityTable, metallicityTableLog
    integer         , intent(in   )                 :: iAge            , iMetallicity
    logical         , intent(in   )                 :: useHighAge
    double precision, intent(in   )                 :: metallicityLow  , metallicityHigh
    integer                                         :: indexAge
    double precision                                :: gradient        , intercept          , &
         &                                             luminosityLow   , luminosityHigh

    ! Select the bounding age of the cell at which the luminosity is taken.
    indexAge=luminosityStellarFromSFHIndexAge(size(table,dim=1),iAge,useHighAge)
    ! Below and above the tabulated range of metallicity the luminosity is constant.
    if      (iMetallicity == 0                      ) then
       luminosityStellarFromSFHMetallicityIntegral=+table(indexAge,1                     )*(metallicityHigh-metallicityLow)
    else if (iMetallicity == size(metallicityTable) ) then
       luminosityStellarFromSFHMetallicityIntegral=+table(indexAge,size(metallicityTable))*(metallicityHigh-metallicityLow)
    else
       ! Within the tabulation the luminosity is linear in log-metallicity, so its integral over linear metallicity is
       ! elementary: with L = intercept + gradient log₁₀Z, ∫L dZ = intercept Z + gradient (Z log₁₀Z - Z/ln10).
       luminosityLow      =table(indexAge,iMetallicity  )
       luminosityHigh     =table(indexAge,iMetallicity+1)
       gradient           =+(luminosityHigh        -luminosityLow        ) &
            &              /(metallicityTableLog(iMetallicity+1)-metallicityTableLog(iMetallicity))
       intercept          =+luminosityLow                                  &
            &              -gradient*metallicityTableLog(iMetallicity)
       luminosityStellarFromSFHMetallicityIntegral=+intercept*(metallicityHigh-metallicityLow)                                   &
            &              +gradient *(                                                                  &
            &                          +metallicityHigh*log10(metallicityHigh)-metallicityHigh/log(10.0d0) &
            &                          -metallicityLow *log10(metallicityLow )+metallicityLow /log(10.0d0) &
            &                         )
    end if
    ! The tabulation's lookup clamps interpolated luminosities to be non-negative; do the same here.
    luminosityStellarFromSFHMetallicityIntegral=max(luminosityStellarFromSFHMetallicityIntegral,0.0d0)
    return
  end function luminosityStellarFromSFHMetallicityIntegral

  double precision function luminosityStellarFromSFHLuminosity(table,metallicityTable,metallicityTableLog,iAge,useHighAge,metallicity)
    !!{RST
    Return the tabulated luminosity at one bounding age of cell ``iAge`` of the tabulation, interpolated linearly in
    :math:`\log_{10}` of metallicity to a single ``metallicity`` and clamped to the tabulated range, as the tabulation's own
    lookup is.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:,:) :: table
    double precision, intent(in   ), dimension(:)   :: metallicityTable, metallicityTableLog
    integer         , intent(in   )                 :: iAge
    logical         , intent(in   )                 :: useHighAge
    double precision, intent(in   )                 :: metallicity
    integer                                         :: indexAge        , iMetallicity
    double precision                                :: h

    indexAge=luminosityStellarFromSFHIndexAge(size(table,dim=1),iAge,useHighAge)
    if      (metallicity <= metallicityTable(1                      )) then
       luminosityStellarFromSFHLuminosity=table(indexAge,1                      )
    else if (metallicity >= metallicityTable(size(metallicityTable)) ) then
       luminosityStellarFromSFHLuminosity=table(indexAge,size(metallicityTable))
    else
       iMetallicity=1
       do while (iMetallicity < size(metallicityTable)-1 .and. metallicityTable(iMetallicity+1) < metallicity)
          iMetallicity=iMetallicity+1
       end do
       h                                 =+(log10(metallicity                  )-metallicityTableLog(iMetallicity)) &
            &                             /(metallicityTableLog(iMetallicity+1)-metallicityTableLog(iMetallicity))
       luminosityStellarFromSFHLuminosity=+table(indexAge,iMetallicity  )*(1.0d0-h) &
            &                             +table(indexAge,iMetallicity+1)*       h
    end if
    luminosityStellarFromSFHLuminosity=max(luminosityStellarFromSFHLuminosity,0.0d0)
    return
  end function luminosityStellarFromSFHLuminosity

  integer function luminosityStellarFromSFHIndexAge(countAgesTable,iAge,useHighAge)
    !!{RST
    Return the index, in the tabulation, of one bounding age of cell ``iAge``. Cells below and above the tabulation have only
    the single bounding age at that edge of it.
    !!}
    implicit none
    integer, intent(in   ) :: countAgesTable, iAge
    logical, intent(in   ) :: useHighAge

    if      (iAge == 0             ) then
       luminosityStellarFromSFHIndexAge=1
    else if (iAge == countAgesTable) then
       luminosityStellarFromSFHIndexAge=countAgesTable
    else if (useHighAge            ) then
       luminosityStellarFromSFHIndexAge=iAge+1
    else
       luminosityStellarFromSFHIndexAge=iAge
    end if
    return
  end function luminosityStellarFromSFHIndexAge

  subroutine luminosityStellarFromSFHWeightsCubic(x,x_,weight)
    !!{RST
    Return the weights of a cubic Lagrange interpolation, at ``x_``, through the four abscissae ``x``. The interpolation is
    local: only the four grid points bracketing the point of interest contribute.
    !!}
    implicit none
    double precision, intent(in   ), dimension(4) :: x
    double precision, intent(in   )               :: x_
    double precision, intent(  out), dimension(4) :: weight
    integer                                       :: i     , j

    weight=1.0d0
    do i=1,4
       do j=1,4
          if (j == i) cycle
          weight(i)=+weight(i)       &
               &    *(x_   -x    (j)) &
               &    /(x (i)-x    (j))
       end do
    end do
    return
  end subroutine luminosityStellarFromSFHWeightsCubic

  function luminosityStellarFromSFHExtract(self,node,time,instance) result(luminosities)
    !!{RST
    Implement a stellar luminosity property extractor which works directly from the star formation history.
    !!}
    use, intrinsic :: ISO_C_Binding                 , only : c_size_t
    use            :: Error                         , only : Error_Report
    use            :: Galactic_Structure_Options    , only : componentTypeDisk    , componentTypeNuclearStarCluster, &
    &                                                        componentTypeSpheroid
    use            :: Galacticus_Nodes              , only : nodeComponentDisk    , nodeComponentNSC               , &
    &                                                        nodeComponentSpheroid
    use            :: Histories                     , only : history
    use            :: ISO_Varying_String            , only : operator(//)         , var_str
    use            :: Stellar_Luminosities_Structure, only : frameObserved
    use            :: String_Handling               , only : operator(//)
    implicit none
    double precision                                                                , allocatable   , dimension(:)   :: luminosities
    class           (nodePropertyExtractorLuminosityStellarFromStarFormationHistory), intent(inout) , target         :: self
    type            (treeNode                                                      ), intent(inout) , target         :: node
    double precision                                                                , intent(in   )                  :: time
    type            (multiCounter                                                  ), intent(inout) , optional       :: instance
    class           (nodeComponentDisk                                             ), pointer                        :: disk
    class           (nodeComponentSpheroid                                         ), pointer                        :: spheroid
    class           (nodeComponentNSC                                              ), pointer                        :: nuclearStarCluster
    double precision                                                                , allocatable   , dimension(:,:) :: masses
    type            (history                                                       )                                 :: starFormationHistory
    integer                                                                                                          :: i                   , &
         &                                                                                                              j                   , &
         &                                                                                                              iRedshift           , &
         &                                                                                                              countWeights
    integer                                                                         , dimension(4)                   :: indexWeight
    double precision                                                                , dimension(4)                   :: weight
    double precision                                                                , dimension(0:1)                 :: hRedshift
    double precision                                                                                                 :: expansionFactor     , &
         &                                                                                                              redshift
    integer         (c_size_t                                                      )                                 :: iRedshift_
    !$GLC attributes unused :: instance

    allocate(luminosities(size(self%filterNames)))
    luminosities=0.0d0
    ! Get the star formation history of the relevant component.
    select case (self%component%ID)
    case (componentTypeDisk              %ID)
       disk                 => node              %disk                ()
       starFormationHistory =  disk              %starFormationHistory()
    case (componentTypeSpheroid          %ID)
       spheroid             => node              %spheroid            ()
       starFormationHistory =  spheroid          %starFormationHistory()
    case (componentTypeNuclearStarCluster%ID)
       nuclearStarCluster   => node              %NSC                 ()
       starFormationHistory =  nuclearStarCluster%starFormationHistory()
    end select
    ! A galaxy with no star formation history has zero luminosity.
    if (.not.starFormationHistory%exists()) return
    ! Get the mass formed in each bin of age and metallicity.
    masses=self%starFormationHistory_%masses(node,starFormationHistory,allowTruncation=.false.)
    ! Build the template if this is the first galaxy. The bin structure of the star formation history is needed to do so, and the
    ! constructor has established that it is the same for every galaxy.
    if (.not.self%templateBuilt) call self%templateBuild(node,starFormationHistory,time,size(masses,dim=1),size(masses,dim=2))
    ! The template is shared by every galaxy, so a galaxy whose star formation history has a different bin structure would be
    ! silently wrong.
    if     (                                                      &
         &   size(masses,dim=1) /= size(self%template,dim=2)      &
         &  .or.                                                  &
         &   size(masses,dim=2) /= size(self%template,dim=3)      &
         & ) call Error_Report(                                                                                             &
         &                     'the star formation history of this galaxy has a different bin structure from the first: ' // &
         &                     'the template can not be shared'                                                           // &
         &                     {introspection:location}                                                                     &
         &                    )
    ! Locate the galaxy in the redshift grid.
    if (self%frame == frameObserved) then
       expansionFactor=self%cosmologyFunctions_%expansionFactor            (time           )
       redshift       =self%cosmologyFunctions_%redshiftFromExpansionFactor(expansionFactor)
       if (redshift < self%redshiftMinimum .or. redshift > self%redshiftMaximum) call reportRangeError(redshift)
       call self%interpolatorRedshift%linearFactors(log(1.0d0+redshift),iRedshift_,hRedshift)
       iRedshift=int(iRedshift_)
       ! Interpolate in log(1+z). A cubic interpolation is used where the grid allows it: the error of a linear interpolation
       ! falls only as the square of the grid spacing, which for a narrow filter forces an impractically dense grid, while a
       ! cubic costs merely two further contractions of the star formation history---which is negligible beside the tabulation
       ! which each further grid point would otherwise require.
       if (self%countRedshifts >= 4) then
          ! Take the four grid points bracketing the galaxy, shifted inward at the ends of the grid.
          countWeights  =4
          indexWeight(1)=min(max(iRedshift-1,1),self%countRedshifts-3)
          do j=2,4
             indexWeight(j)=indexWeight(1)+j-1
          end do
          call luminosityStellarFromSFHWeightsCubic(log(1.0d0+self%redshiftGrid(indexWeight(1:4))),log(1.0d0+redshift),weight)
       else
          countWeights     =2
          indexWeight(1:2) =[iRedshift,iRedshift+1]
          weight     (1:2) =[hRedshift(0),hRedshift(1)]
       end if
    else
       countWeights   =1
       indexWeight(1) =1
       weight     (1) =1.0d0
    end if
    ! Contract the star formation history with the template at each grid redshift which carries weight, and combine the results.
    ! This is identical to interpolating the template itself, and avoids building an interpolated template for every galaxy.
    do i=1,size(self%filterNames)
       do j=1,countWeights
          luminosities(i)=+luminosities(i)                                             &
               &          +sum(self%template(i,:,:,indexWeight(j))*masses)*weight(j)
       end do
    end do
    return

  contains

    subroutine reportRangeError(redshift)
      !!{RST
      Report that a galaxy lies outside the tabulated range of redshift.
      !!}
      use :: Display           , only : displayGreen   , displayReset
      use :: ISO_Varying_String, only : operator(//)   , var_str
      use :: String_Handling   , only : stringXMLFormat
      implicit none
      double precision                , intent(in   ) :: redshift
      character       (len=12        )                :: label       , labelMinimum, &
           &                                             labelMaximum
      type            (varying_string)                :: message

      !$omp critical(gfortranInternalIO)
      write (label       ,'(f12.6)') redshift
      write (labelMinimum,'(f12.6)') self%redshiftMinimum
      write (labelMaximum,'(f12.6)') self%redshiftMaximum
      !$omp end critical(gfortranInternalIO)
      message=var_str('a galaxy at redshift ')//trim(adjustl(label))//' lies outside the tabulated range ['// &
           &  trim(adjustl(labelMinimum))//','//trim(adjustl(labelMaximum))//']'//char(10)
      message=message//displayGreen()//'HELP:'//displayReset()//' the luminosities are tabulated on a grid of redshifts and'//char(10)
      message=message//'      interpolated to the redshift of each galaxy, and are not extrapolated beyond the grid, since'//char(10)
      message=message//'      the postprocessing of the spectrum varies sharply with redshift. Widen the range in your'//char(10)
      message=message//'      parameter file, for example to:'//char(10)//char(10)
      !$omp critical(gfortranInternalIO)
      write (labelMinimum,'(f12.6)') min(self%redshiftMinimum,0.9d0*redshift)
      write (labelMaximum,'(f12.6)') max(self%redshiftMaximum,1.1d0*redshift)
      !$omp end critical(gfortranInternalIO)
      message=message//stringXMLFormat('<nodePropertyExtractor value="luminosityStellarFromStarFormationHistory">**B<redshiftMinimum value="'//trim(adjustl(labelMinimum))//'"/>**C**B<redshiftMaximum value="'//trim(adjustl(labelMaximum))//'"/>**C</nodePropertyExtractor>',indentInitial=6)//char(10)
      call Error_Report(message//{introspection:location})
      return
    end subroutine reportRangeError

  end function luminosityStellarFromSFHExtract

  function luminosityStellarFromSFHQuantity(self)
    !!{RST
    Return the class of the stellar luminosity property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyQuantityLuminosity
    implicit none
    type (enumerationOutputAnalysisPropertyQuantityType                 )                :: luminosityStellarFromSFHQuantity
    class(nodePropertyExtractorLuminosityStellarFromStarFormationHistory), intent(inout) :: self
    !$GLC attributes unused :: self

    luminosityStellarFromSFHQuantity=outputAnalysisPropertyQuantityLuminosity
    return
  end function luminosityStellarFromSFHQuantity

  subroutine luminosityStellarFromSFHNames(self,time,names)
    !!{RST
    Return the names of the ``luminosityStellarFromStarFormationHistory`` properties.
    !!}
    use :: Galactic_Structure_Options    , only : enumerationComponentTypeDecode
    use :: ISO_Varying_String            , only : operator(//)
    use :: Stellar_Luminosities_Structure, only : enumerationFrameDecode
    implicit none
    class           (nodePropertyExtractorLuminosityStellarFromStarFormationHistory), intent(inout)                            :: self
    double precision                                                                , intent(in   )                            :: time
    type            (varying_string                                                ), intent(inout), allocatable, dimension(:) :: names
    integer                                                                                                                    :: i
    !$GLC attributes unused :: time

    allocate(names(size(self%filterNames)))
    do i=1,size(names)
       names(i)=enumerationComponentTypeDecode(self%component,includePrefix=.false.)                  // &
            &   "StellarSFH:"                                                                         // &
            &   self%stellarPopulationSpectraPostprocessor_%objectType(short=.true.)                  // &
            &   ":"                                                                                   // &
            &   self%filterNames(i)                                                                   // &
            &   ":"                                                                                   // &
            &   enumerationFrameDecode(self%frame,includePrefix=.false.)
    end do
    return
  end subroutine luminosityStellarFromSFHNames

  subroutine luminosityStellarFromSFHDescriptions(self,time,descriptions)
    !!{RST
    Return descriptions of the ``luminosityStellarFromStarFormationHistory`` properties.
    !!}
    use :: Galactic_Structure_Options    , only : enumerationComponentTypeDecode
    use :: ISO_Varying_String            , only : operator(//)
    use :: Stellar_Luminosities_Structure, only : enumerationFrameDecode
    implicit none
    class           (nodePropertyExtractorLuminosityStellarFromStarFormationHistory), intent(inout)                            :: self
    double precision                                                                , intent(in   )                            :: time
    type            (varying_string                                                ), intent(inout), allocatable, dimension(:) :: descriptions
    integer                                                                                                                    :: i
    !$GLC attributes unused :: time

    allocate(descriptions(size(self%filterNames)))
    do i=1,size(descriptions)
       descriptions(i)="Luminosity in "                                                               // &
            &          self%filterNames(i)                                                           // &
            &          " filter in the "                                                              // &
            &          enumerationFrameDecode(self%frame,includePrefix=.false.)                       // &
            &          " frame in units of the AB-magnitude system zero-point; computed directly from the "// &
            &          enumerationComponentTypeDecode(self%component,includePrefix=.false.)           // &
            &          " star formation history."
    end do
    return
  end subroutine luminosityStellarFromSFHDescriptions

  function luminosityStellarFromSFHUnitsInSI(self,time) result(unitsInSI)
    !!{RST
    Return the units of the ``luminosityStellarFromStarFormationHistory`` properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : luminosityZeroPointAB
    implicit none
    double precision                                                                , allocatable  , dimension(:) :: unitsInSI
    class           (nodePropertyExtractorLuminosityStellarFromStarFormationHistory), intent(inout)               :: self
    double precision                                                                , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(unitsInSI(size(self%filterNames)))
    unitsInSI=luminosityZeroPointAB
    return
  end function luminosityStellarFromSFHUnitsInSI

  function luminosityStellarFromSFHUnits(self,time) result(units)
    !!{RST
    Return the units of the ``luminosityStellarFromStarFormationHistory`` properties.
    !!}
    implicit none
    type            (unitType                                                      ), allocatable  , dimension(:) :: units
    class           (nodePropertyExtractorLuminosityStellarFromStarFormationHistory), intent(inout)               :: self
    double precision                                                                , intent(in   )               :: time
    double precision                                                                , allocatable  , dimension(:) :: siValues
    integer                                                                                                       :: i

    siValues=self%unitsInSI(time)
    allocate(units(size(siValues)))
    do i=1,size(siValues)
       units(i)=unitType(siValues(i),description='AB-magnitude zero point',quantity='4.465920e17 W/Hz')
    end do
    return
  end function luminosityStellarFromSFHUnits
