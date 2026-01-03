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

!!{
Implements the survey geometry used by \cite{caputi_stellar_2011}.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <surveyGeometry name="surveyGeometryCaputi2011UKIDSSUDS">
   <description>
    A survey geometry class which implements the UKIDSS UDS survey used by \cite{caputi_stellar_2011}. The survey window function is
    determined from the set of galaxy positions provided by Caputi (private communication), by finding a suitable bounding box and
    then cutting out empty regions (corresponding to regions that were removed around bright stars). A set of random points are then
    found within this mask and are used to find the Fourier transform of the survey volume.
    
    To estimate the depth of the \cite{caputi_stellar_2011} sample as a function of galaxy stellar mass we make use of semi-analytic
    models in the Millennium Database. Specifically, we use the \glspl{sam} of \cite{guo_dwarf_2011} and
    \cite{henriques_confronting_2012} specifically the {\normalfont \ttfamily Guo2010a..MR} and {\normalfont \ttfamily
      Henriques2012a.wmap1.BC03\_001} tables in the Millennium Database. For each snapshot in the database, we extract the stellar
    masses and observed-frame IRAC 4.5$\mu$m apparent magnitudes (including dust extinction), and determine the median apparent
    magnitude as a function of stellar mass. Using the limiting apparent magnitude of the \cite{caputi_stellar_2011} sample,
    $i_{4.5}=24$, we infer the corresponding absolute magnitude at each redshift and, using our derived apparent magnitude--stellar
    mass relation, infer the corresponding stellar mass.
    
    The end result of this procedure is the limiting stellar mass as a function of redshift, accounting for k-corrections, evolution,
    and the effects of dust. Figure~\ref{fig:UKIDSSUDSMassRedshift} shows the resulting relation between stellar mass and the maximum
    redshift at which such a galaxy would be included in the sample. Points indicate measurements from the \gls{sam}, while the line
    shows a polynomial fit:
    \begin{equation}
     z(M_\star) = -56.247 + 5.881 m,
     \label{eq:UKIDSSUDSDepthPolynomial}
    \end{equation}
    where $m= \log_{10}(M_\star/M_\odot)$. We use this polynomial fit to determine the depth of the sample as a function of stellar mass.
    
    \begin{figure}
     \begin{center}
     \includegraphics[width=85mm,trim=0mm 0mm 0mm 4mm,clip]{Plots/DataAnalysis/UKIDSSUDSMassLuminosityRelation.pdf}
     \caption{The maximum redshift at which a galaxy of given stellar mass can be detected in the sample of
       \protect\cite{caputi_stellar_2011}. Points show the results obtained using the \protect\cite{henriques_confronting_2012} model
       from the Millennium Database, while the lines shows a polynomial fit to these results (given in
       eqn.~\ref{eq:UKIDSSUDSDepthPolynomial}).}
     \end{center}
     \label{fig:UKIDSSUDSMassRedshift}
    \end{figure}
   </description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryRandomPoints) :: surveyGeometryCaputi2011UKIDSSUDS
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     integer                                            :: redshiftBin
     double precision                                   :: binDistanceMinimum           , binDistanceMaximum, &
          &                                                redshiftMinimum              , redshiftMaximum
   contains
     final     ::                      caputi2011UKIDSSUDSDestructor
     procedure :: distanceMinimum   => caputi2011UKIDSSUDSDistanceMinimum
     procedure :: distanceMaximum   => caputi2011UKIDSSUDSDistanceMaximum
     procedure :: volumeMaximum     => caputi2011UKIDSSUDSVolumeMaximum
     procedure :: solidAngle        => caputi2011UKIDSSUDSSolidAngle
     procedure :: randomsInitialize => caputi2011UKIDSSUDSRandomsInitialize
  end type surveyGeometryCaputi2011UKIDSSUDS

  interface surveyGeometryCaputi2011UKIDSSUDS
     !!{
     Constructors for the \cite{caputi_stellar_2011} survey geometry class.
     !!}
     module procedure caputi2011UKIDSSUDSConstructorParameters
     module procedure caputi2011UKIDSSUDSConstructorInternal
  end interface surveyGeometryCaputi2011UKIDSSUDS

contains

  function caputi2011UKIDSSUDSConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the \cite{caputi_stellar_2011} conditional mass function class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (surveyGeometryCaputi2011UKIDSSUDS)                :: self
    type   (inputParameters                  ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class  (randomNumberGeneratorClass       ), pointer       :: randomNumberGenerator_
    integer                                                   :: redshiftBin

    !![
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    <inputParameter>
      <name>redshiftBin</name>
      <source>parameters</source>
      <description>The redshift bin (0, 1, or 2) of the \cite{caputi_stellar_2011} to use.</description>
    </inputParameter>
    !!]
    self=surveyGeometryCaputi2011UKIDSSUDS(redshiftBin,cosmologyFunctions_,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function caputi2011UKIDSSUDSConstructorParameters

  function caputi2011UKIDSSUDSConstructorInternal(redshiftBin,cosmologyFunctions_,randomNumberGenerator_) result(self)
    !!{
    Internal constructor for the \cite{caputi_stellar_2011} conditional mass function class.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Error                      , only : Error_Report
    implicit none
    type            (surveyGeometryCaputi2011UKIDSSUDS)                                  :: self
    integer                                            , intent(in   )                   :: redshiftBin
    class           (cosmologyFunctionsClass          ), intent(in   ), target           :: cosmologyFunctions_
    class           (randomNumberGeneratorClass       ), intent(in   ), target, optional :: randomNumberGenerator_
    !![
    <constructorAssign variables="redshiftBin, *cosmologyFunctions_, *randomNumberGenerator_"/>
    !!]

    ! Find distance limits for this redshift bin.
    select case (redshiftBin)
    case(0)
       self%redshiftMinimum=3.00d0
       self%redshiftMaximum=3.50d0
    case(1)
       self%redshiftMinimum=3.50d0
       self%redshiftMaximum=4.25d0
    case(2)
       self%redshiftMinimum=4.25d0
       self%redshiftMaximum=5.00d0
    case default
       call Error_Report('0≤redshiftBin≤2 is required'//{introspection:location})
    end select
    self%binDistanceMinimum                                                                 &
         & =self%cosmologyFunctions_%distanceComovingConvert(                               &
         &                                                   output  =distanceTypeComoving, &
         &                                                   redshift=self%redshiftMinimum  &
         &                                                  )
    self%binDistanceMaximum                                                                 &
         & =self%cosmologyFunctions_%distanceComovingConvert(                               &
         &                                                   output  =distanceTypeComoving, &
         &                                                   redshift=self%redshiftMaximum  &
         &                                                  )
    self%geometryInitialized=.false.
    return
  end function caputi2011UKIDSSUDSConstructorInternal

  subroutine caputi2011UKIDSSUDSDestructor(self)
    !!{
    Destructor for the \refClass{surveyGeometryCaputi2011UKIDSSUDS} survey geometry class.
    !!}
    implicit none
    type(surveyGeometryCaputi2011UKIDSSUDS), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine caputi2011UKIDSSUDSDestructor

  double precision function caputi2011UKIDSSUDSDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the minimum distance at which a galaxy is included.
    !!}
    implicit none
    class           (surveyGeometryCaputi2011UKIDSSUDS), intent(inout)           :: self
    double precision                                   , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                          luminosity, starFormationRate
    integer                                            , intent(in   ), optional :: field
    !$GLC attributes unused :: mass, field, magnitudeAbsolute, luminosity, starFormationRate

    caputi2011UKIDSSUDSDistanceMinimum=self%binDistanceMinimum
    return
  end function caputi2011UKIDSSUDSDistanceMinimum

  double precision function caputi2011UKIDSSUDSDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Error                      , only : Error_Report
    implicit none
    class           (surveyGeometryCaputi2011UKIDSSUDS), intent(inout)           :: self
    double precision                                   , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                          luminosity, starFormationRate
    integer                                            , intent(in   ), optional :: field
    double precision                                                             :: redshift  , logarithmicMass

    ! Validate field.
    if (present(field).and.field /= 1) call Error_Report('field = 1 required'//{introspection:location})
    ! Validate arguments.
    if (present(magnitudeAbsolute)) call Error_Report('`magnitudeAbsolute` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! Find the limiting redshift for this mass using a fit derived from Millennium Simulation SAMs. (See
    ! constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/massLuminosityRelation.pl for details.)
    if (present(mass)) then
       logarithmicMass=log10(mass)
       redshift=-56.247426278132d0+logarithmicMass*(5.88091022342758d0)
       ! Convert from redshift to comoving distance.
       caputi2011UKIDSSUDSDistanceMaximum                                                                                             &
            &=self%cosmologyFunctions_%distanceComovingConvert(                                                                       &
            &                                                  output  =distanceTypeComoving                                        , &
            &                                                  redshift=min(max(redshift,self%redshiftMinimum),self%redshiftMaximum)  &
            &                                                 )
       ! Limit the maximum distance.
       caputi2011UKIDSSUDSDistanceMaximum=min(caputi2011UKIDSSUDSDistanceMaximum,self%binDistanceMaximum)
    else
       caputi2011UKIDSSUDSDistanceMaximum=                                       self%binDistanceMaximum
    end if
       return
  end function caputi2011UKIDSSUDSDistanceMaximum

  double precision function caputi2011UKIDSSUDSVolumeMaximum(self,mass,field)
    !!{
    Compute the maximum volume within which a galaxy is visible.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (surveyGeometryCaputi2011UKIDSSUDS), intent(inout)           :: self
    double precision                                   , intent(in   )           :: mass
    integer                                            , intent(in   ), optional :: field

    ! Validate field.
    if (present(field).and.field /= 1) call Error_Report('field = 1 required'//{introspection:location})
    ! Compute the volume.
    caputi2011UKIDSSUDSVolumeMaximum                             &
         & =max(                                                 &
         &       0.0d0                                         , &
         &       self%solidAngle(field)                          &
         &      *(                                               &
         &        +self%distanceMaximum   (mass,field=field)**3  &
         &        -self%binDistanceMinimum                  **3  &
         &       )                                               &
         &      /3.0d0                                           &
         &     )
    return
  end function caputi2011UKIDSSUDSVolumeMaximum

  double precision function caputi2011UKIDSSUDSSolidAngle(self,field)
    !!{
    Return the solid angle of the \cite{caputi_stellar_2011} sample. Computed from survey mask (see {\normalfont \ttfamily
    constraints/dataAnalysis/stellarMassFunctions\_UKIDSS\_UDS\_z3\_5/surveyGeometryRandoms.pl}).
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (surveyGeometryCaputi2011UKIDSSUDS), intent(inout)           :: self
    integer                                            , intent(in   ), optional :: field
    double precision                                   , parameter               :: solidAngleSurvey=1.59233703487973d-4
    !$GLC attributes unused :: self

    ! Validate field.
    if (present(field).and.field /= 1) call Error_Report('field = 1 required'//{introspection:location})
    caputi2011UKIDSSUDSSolidAngle=solidAngleSurvey
    return
  end function caputi2011UKIDSSUDSSolidAngle

  subroutine caputi2011UKIDSSUDSRandomsInitialize(self)
    !!{
    Load random points for the survey.
    !!}
    use :: File_Utilities , only : File_Exists
    use :: Error          , only : Error_Report
    use :: Input_Paths    , only : inputPath        , pathTypeDataDynamic
    use :: HDF5_Access    , only : hdf5Access
    use :: IO_HDF5        , only : hdf5Object
    use :: String_Handling, only : operator(//)
    use :: System_Command , only : System_Command_Do
    implicit none
    class(surveyGeometryCaputi2011UKIDSSUDS), intent(inout) :: self
    type (hdf5Object                       )                :: surveyGeometryRandomsFile

    ! Generate the randoms file if necessary.
    if (.not.File_Exists(inputPath(pathTypeDataDynamic)//&
         & "surveys/UKIDSS_UDS/data/surveyGeometryRandoms.hdf5")) then
       call System_Command_Do(inputPath(pathTypeDataDynamic)//"surveyGeometry/UKIDDS_UDS/surveyGeometryRandoms.pl")
       if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"surveys/UKIDSS_UDS/surveyGeometryRandoms.hdf5")) call Error_Report('unable to create survey geometry randoms file'//{introspection:location})
    end if
    ! Read the distribution of random points from file.
    !$ call hdf5Access%set()
    call surveyGeometryRandomsFile%openFile(char(inputPath(pathTypeDataDynamic)//&
         &'surveys/UKIDSS_UDS/surveyGeometryRandoms.hdf5')&
         &,readOnly=.true.)
    call surveyGeometryRandomsFile%readDataset('theta',self%randomTheta)
    call surveyGeometryRandomsFile%readDataset('phi'  ,self%randomPhi  )
    call surveyGeometryRandomsFile%close()
    !$ call hdf5Access%unset()
    return
  end subroutine caputi2011UKIDSSUDSRandomsInitialize
