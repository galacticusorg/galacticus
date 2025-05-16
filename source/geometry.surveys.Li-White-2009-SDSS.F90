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
Implements the survey geometry of the SDSS sample used by \cite{li_distribution_2009}.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <surveyGeometry name="surveyGeometryLiWhite2009SDSS">
   <description>
    A survey geometry class that describes the survey geometry of \cite{li_distribution_2009}. 
    
    For the angular mask, we make use of the catalog of random points within the survey footprint provided by the
    NYU-VAGC\footnote{Specifically, \href{https://zenodo.org/records/10257229/files/lss_random-0.dr72.dat}{https://zenodo.org/records/10257229/files/lss\_random-0.dr72.dat}
     (which is a copy of the dataset originally found at the, now defunct, URL {\normalfont \ttfamily http://sdss.physics.nyu.edu/lss/dr72/random/lss\_random-0.dr72.dat}).}
    (\citealt{blanton_new_2005}; see also
    \citealt{adelman-mccarthy_sixth_2008,padmanabhan_improved_2008}). \cite{li_distribution_2009} consider only the main,
    contiguous region and so we keep only those points which satisfy RA$>100^\circ$, RA$&lt;300^\circ$, and RA$&lt;247^\circ$ or
    $\delta&lt; 51^\circ$. When the survey window function is needed, these points are used to determine which elements of a 3D
    grid fall within the window function.
    
    To estimate the depth of the \cite{li_distribution_2009} sample as a function of galaxy stellar mass we make use of
    semi-analytic models in the Millennium Database. Specifically, we use the \gls{sam} of
    \citeauthor{de_lucia_hierarchical_2007}~(\citeyear{de_lucia_hierarchical_2007}; specifically the {\normalfont \ttfamily
    millimil..DeLucia2006a} and {\normalfont \ttfamily millimil..DeLucia2006a\_sdss2mass} tables in the Millennium
    Database). For each snapshot in the database, we extract the stellar masses and observed-frame SDSS r-band absolute
    magnitudes (including dust extinction), and determine the median absolute magnitude as a function of stellar mass. Using
    the limiting apparent magnitude of the \cite{li_distribution_2009} sample, $r=17.6$, we infer the corresponding absolute
    magnitude at each redshift and, using our derived absolute magnitude--stellar mass relation, infer the corresponding
    stellar mass.
    
    The end result of this procedure is the limiting stellar mass as a function of redshift, accounting for k-corrections,
    evolution, and the effects of dust. Figure~\ref{fig:SDSSDepthFit} shows the resulting relation between stellar mass and the
    maximum redshift at which such a galaxy would be included in the sample. Points indicate measurements from the \gls{sam},
    while the line shows a polynomial fit:
    \begin{eqnarray}
     z(M_\star) &amp;=&amp; -5.950 + 2.638 m - 0.4211 m^2 \nonumber \\
                &amp; &amp; + 2.852\times 10^{-2} m^3 - 6.783 \times 10^{-4} m^4,
     \label{eq:DepthPolynomial}
    \end{eqnarray}
    where $m= \log_{10}(M_\star/M_\odot)$. We use this polynomial fit to determine the depth of the sample as a function of
    stellar mass. We adopt a solid angle of $2.1901993$~sr \citep{percival_shape_2007} for the sample.
    
    \begin{figure}
     \begin{center}
     \includegraphics[width=85mm,trim=0mm 0mm 0mm 4mm,clip]{Plots/DataAnalysis/SDSSMassLuminosityRelation.pdf}
     \end{center}
     \caption{The maximum redshift at which a galaxy of given stellar mass can be detected in the sample of
     \protect\cite{li_distribution_2009}. Points show the results obtained using the \protect\cite{de_lucia_hierarchical_2007}
     model from the Millennium Database, while the lines shows a polynomial fit to these results (given in
     eqn.~\ref{eq:DepthPolynomial}).}
     \label{fig:SDSSDepthFit}
    \end{figure}
   </description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryRandomPoints) :: surveyGeometryLiWhite2009SDSS
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_  => null()
     double precision                                   :: redshiftMinimum               , redshiftMaximum
     double precision                                   :: limitDistanceMinimum          , limitDistanceMaximum
   contains
     final     ::                      liWhite2009SDSSDestructor
     procedure :: distanceMinimum   => liWhite2009SDSSDistanceMinimum
     procedure :: distanceMaximum   => liWhite2009SDSSDistanceMaximum
     procedure :: solidAngle        => liWhite2009SDSSSolidAngle
     procedure :: randomsInitialize => liWhite2009SDSSRandomsInitialize
  end type surveyGeometryLiWhite2009SDSS

  interface surveyGeometryLiWhite2009SDSS
     !!{
     Constructors for the \cite{li_distribution_2009} survey geometry class.
     !!}
     module procedure liWhite2009SDSSConstructorParameters
     module procedure liWhite2009SDSSConstructorInternal
  end interface surveyGeometryLiWhite2009SDSS

contains

  function liWhite2009SDSSConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{li_distribution_2009} survey geometry class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (surveyGeometryLiWhite2009SDSS)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class           (randomNumberGeneratorClass   ), pointer       :: randomNumberGenerator_
    double precision                                               :: redshiftMinimum       , redshiftMaximum

    !![
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    <inputParameter>
      <name>redshiftMinimum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum redshift for the survey.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftMaximum</name>
      <defaultValue>huge(1.0d0)</defaultValue>
      <source>parameters</source>
      <description>The maximum redshift for the survey.</description>
    </inputParameter>
    !!]
    ! Build the object.
    self=surveyGeometryLiWhite2009SDSS(redshiftMinimum,redshiftMaximum,cosmologyFunctions_,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function liWhite2009SDSSConstructorParameters

  function liWhite2009SDSSConstructorInternal(redshiftMinimum,redshiftMaximum,cosmologyFunctions_,randomNumberGenerator_) result(self)
    !!{
    Constructor for the \cite{li_distribution_2009} survey geometry class which allows specification of minimum and maximum redshifts.
    !!}
    use :: Cosmology_Functions        , only : cosmologyFunctionsClass
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    implicit none
    type            (surveyGeometryLiWhite2009SDSS)                                  :: self
    double precision                               , intent(in   )                   :: redshiftMinimum    , redshiftMaximum
    class           (cosmologyFunctionsClass      ), intent(in   ), target           :: cosmologyFunctions_
    class           (randomNumberGeneratorClass   ), intent(in   ), target, optional :: randomNumberGenerator_ 
    !![
    <constructorAssign variables="*cosmologyFunctions_, *randomNumberGenerator_, redshiftMinimum, redshiftMaximum"/>
    !!]

    self   %geometryInitialized =.false.
    self   %limitDistanceMinimum=self%cosmologyFunctions_%distanceComovingConvert(                               &
         &                                                                        output  =distanceTypeComoving, &
         &                                                                        redshift=redshiftMinimum       &
         &                                                                       )
    if (redshiftMaximum < huge(1.0d0)) then
       self%limitDistanceMaximum=self%cosmologyFunctions_%distanceComovingConvert(                               &
            &                                                                     output  =distanceTypeComoving, &
            &                                                                     redshift=redshiftMaximum       &
            &                                                                    )
    else
       self%limitDistanceMaximum=huge(1.0d0)
    end if
    return
  end function liWhite2009SDSSConstructorInternal

  subroutine liWhite2009SDSSDestructor(self)
    !!{
    Destructor for the \refClass{surveyGeometryLiWhite2009SDSS} survey geometry class.
    !!}
    implicit none
    type(surveyGeometryLiWhite2009SDSS), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine liWhite2009SDSSDestructor

  double precision function liWhite2009SDSSDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the minimum distance at which a galaxy is visible.
    !!}
    implicit none
    class           (surveyGeometryLiWhite2009SDSS), intent(inout)           :: self
    double precision                               , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                      luminosity, starFormationRate
    integer                                        , intent(in   ), optional :: field
    !$GLC attributes unused :: mass, field, magnitudeAbsolute, luminosity, starFormationRate

    liWhite2009SDSSDistanceMinimum=self%limitDistanceMinimum
    return
  end function liWhite2009SDSSDistanceMinimum

  double precision function liWhite2009SDSSDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Error                      , only : Error_Report
    implicit none
    class           (surveyGeometryLiWhite2009SDSS), intent(inout)           :: self
    double precision                               , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                      luminosity, starFormationRate
    integer                                        , intent(in   ), optional :: field
    double precision                                                         :: redshift, logarithmicMass

    ! Validate field.
    if (present(field).and.field /= 1) call Error_Report('field = 1 required'//{introspection:location})
    ! Validate arguments.
    if (present(magnitudeAbsolute)) call Error_Report('`magnitudeAbsolute` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! Find the limiting redshift for this mass using a fit derived from Millennium Simulation SAMs. (See
    ! constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07/massLuminosityRelation.pl for details.)
    if (present(mass)) then
       logarithmicMass=log10(mass)
       redshift=                                   &
            & max(                                 &
            &     -5.9502006195004d0               &
            &     +logarithmicMass                 &
            &     *(                               &
            &       +2.63793788603951d0            &
            &       +logarithmicMass               &
            &       *(                             &
            &         -0.421075858899237d0         &
            &         +logarithmicMass             &
            &         *(                           &
            &           +0.0285198776926787d0      &
            &           +logarithmicMass           &
            &           *(                         &
            &             -0.000678327494720407d0  &
            &            )                         &
            &          )                           &
            &        )                             &
            &      )                             , &
            &     +0.0d0                           &
            &    )
       ! Convert from redshift to comoving distance.
       liWhite2009SDSSDistanceMaximum=min(                                                                                &
            &                             self%limitDistanceMaximum                                                     , &
            &                             self%cosmologyFunctions_%distanceComovingConvert(                               &
            &                                                                              output  =distanceTypeComoving, &
            &                                                                              redshift=redshift              &
            &                                                                             )                               &
            &                            )
    else
       liWhite2009SDSSDistanceMaximum=    self%limitDistanceMaximum
    end if
    return
  end function liWhite2009SDSSDistanceMaximum

  double precision function liWhite2009SDSSSolidAngle(self,field)
    !!{
    Return the solid angle of the \cite{li_distribution_2009} sample.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (surveyGeometryLiWhite2009SDSS), intent(inout)           :: self
    integer                                        , intent(in   ), optional :: field
    double precision                               , parameter               :: solidAngleSurvey=2.1901993d0 ! From Percival et al. (2010; MNRAS; 401; 2148)
    !$GLC attributes unused :: self

    ! Validate field.
    if (present(field).and.field /= 1) call Error_Report('field = 1 required'//{introspection:location})
    liWhite2009SDSSSolidAngle=solidAngleSurvey
    return
  end function liWhite2009SDSSSolidAngle

  subroutine liWhite2009SDSSRandomsInitialize(self)
    !!{
    Compute the window function for the survey.
    !!}
    use :: Display                 , only : displayMessage
    use :: File_Utilities          , only : Count_Lines_In_File, Directory_Make     , File_Exists, File_Lock, &
         &                                  File_Unlock        , lockDescriptor
    use :: Error                   , only : Error_Report
    use :: Input_Paths             , only : inputPath          , pathTypeDataDynamic
    use :: ISO_Varying_String      , only : varying_string
    use :: Numerical_Constants_Math, only : Pi
    use :: String_Handling         , only : operator(//)
    use :: System_Download         , only : download
    implicit none
    class           (surveyGeometryLiWhite2009SDSS), intent(inout)             :: self
    double precision                               , allocatable, dimension(:) :: angleTmp
    integer                                                                    :: randomsCount  , j          , &
         &                                                                        i             , randomUnit
    double precision                                                           :: rightAscension, declination
    type            (varying_string               )                            :: message
    type            (lockDescriptor               )                            :: lock
    integer                                                                    :: status

    ! Randoms file obtained from:  http://sdss.physics.nyu.edu/lss/dr72/random/
    if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"surveyGeometry/lss_random-0.dr72.dat")) then
       call Directory_Make(inputPath(pathTypeDataDynamic)//"surveyGeometry")
       call File_Lock  (char(inputPath(pathTypeDataDynamic)//"surveyGeometry/lss_random-0.dr72.dat"),lock,lockIsShared=.false.)
       if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"surveyGeometry/lss_random-0.dr72.dat")) then
          call download("https://zenodo.org/records/10257229/files/lss_random-0.dr72.dat",char(inputPath(pathTypeDataDynamic))//"surveyGeometry/lss_random-0.dr72.dat",status=status)
          if (status /= 0 .or. .not.File_Exists(inputPath(pathTypeDataDynamic)//"surveyGeometry/lss_random-0.dr72.dat")) call Error_Report('unable to download SDSS survey geometry randoms file'//{introspection:location})
       end if
       call File_Unlock(                     lock                     )
    end if
    randomsCount=Count_Lines_In_File(inputPath(pathTypeDataDynamic)//"surveyGeometry/lss_random-0.dr72.dat")
    allocate(self%randomTheta(randomsCount))
    allocate(self%randomPhi  (randomsCount))
    open(newUnit=randomUnit,file=char(inputPath(pathTypeDataDynamic)//"surveyGeometry/lss_random-0.dr72.dat"),status="old",form="formatted")
    j=0
    do i=1,randomsCount
       read (randomUnit,*) rightAscension,declination
       if     (                                &
            &          rightAscension  > 100.0 &
            &  .and.   rightAscension  < 300.0 &
            &  .and. (                         &
            &          rightAscension  < 247.0 &
            &         .or. declination <  51.0 &
            &        )                         &
            & ) then
          j=j+1
          self%randomTheta(j)=Pi*(90.0d0-declination   )/180.0d0
          self%randomPhi  (j)=Pi*        rightAscension /180.0d0
       end if
    end do
    close(randomUnit)
    randomsCount=j
    call Move_Alloc   (self%randomTheta,angleTmp      )
    allocate(self%randomTheta(randomsCount))
    self%randomTheta=angleTmp(1:randomsCount)
    deallocate(angleTmp                       )
    call Move_Alloc   (self%randomPhi  ,angleTmp      )
    allocate(self%randomPhi  (randomsCount))
    self%randomPhi  =angleTmp(1:randomsCount)
    deallocate(angleTmp                       )
    message="Read "
    message=message//randomsCount//" random points and kept "//randomsCount//" of them"
    call displayMessage(message)
    return
  end subroutine liWhite2009SDSSRandomsInitialize
