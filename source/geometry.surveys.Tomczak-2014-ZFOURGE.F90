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
Implements the geometry of the ZFOURGE survey used by \cite{tomczak_galaxy_2014}.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <surveyGeometry name="surveyGeometryTomczak2014ZFOURGE">
   <description>
    A survey geometry class that describes the survey geometry of \cite{tomczak_galaxy_2014}. 
    
    For the angular mask, we make use of {\normalfont \ttfamily mangle} polygon files constructed by hand using vertices
    matched approximately to the distribution of galaxies in the survey (positions of which were provided by R.~Quadri; private
    communication). The solid angle of each mask is computed using the {\normalfont \ttfamily mangle} {\normalfont \ttfamily
    harmonize} command.
    
    To determine the depth as a function of stellar mass, we make use of the tabulated mass completeness limits as a function
    of redshift for ZFOURGE and NMBS fields provided by R.~Quadri (private communication). These are fit with fourth-order
    polynomials. Figure~\ref{fig:Tomczak2014DepthFit} shows the resulting relation between stellar mass and the maximum
    redshift at which such a galaxy would be included in the sample. Dotted lines indicate the tabulated result from ZFOURGE,
    while the lines show polynomial fits:
    \begin{equation}
     z_\mathrm{max}(M_\star) = \left\{ \begin{array}{ll} -114.66+m*(45.901+m*(-6.1617+m*(0.27822))) &amp; \hbox{ZFOURGE fields} \\
     -58.483+m*(20.250+m*(-2.3563+m*(0.092705))) &amp; \hbox{NMBS fields} \end{array} \right.
     \label{eq:TomczakDepthPolynomial}
    \end{equation}
    where $m= \log_{10}(M_\star/M_\odot)$. We use this polynomial fit to determine the depth of the sample as a function of
    stellar mass.
    
    \begin{figure}
     \begin{center}
     \includegraphics[width=85mm,trim=0mm 0mm 0mm 4mm,clip]{Plots/DataAnalysis/TomczakZFOURGEMassRedshiftRelation.pdf}
     \end{center}
     \caption{The maximum redshift at which a galaxy of given stellar mass can be detected in the sample of
     \protect\cite{tomczak_galaxy_2014}. Points show the results obtained from data provided by Davidzon, while the lines shows
     a polynomial fit to these results (given in eqn.~\ref{eq:TomczakDepthPolynomial}).}
     \label{fig:Tomczak2014DepthFit}
    \end{figure}
   </description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryMangle) :: surveyGeometryTomczak2014ZFOURGE
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     integer                                            :: redshiftBin
     double precision                                   :: binDistanceMinimum           , binDistanceMaximum, &
          &                                                redshiftMinimum              , redshiftMaximum
   contains
     final     ::                              tomczak2014ZFOURGEDestructor
     procedure :: fieldCount                => tomczak2014ZFOURGEFieldCount
     procedure :: distanceMinimum           => tomczak2014ZFOURGEDistanceMinimum
     procedure :: distanceMaximum           => tomczak2014ZFOURGEDistanceMaximum
     procedure :: volumeMaximum             => tomczak2014ZFOURGEVolumeMaximum
     procedure :: angularPowerMaximumDegree => tomczak2014ZFOURGEAngularPowerMaximumDegree
     procedure :: mangleDirectory           => tomczak2014ZFOURGEMangleDirectory
     procedure :: mangleFiles               => tomczak2014ZFOURGEMangleFiles
   end type surveyGeometryTomczak2014ZFOURGE

  interface surveyGeometryTomczak2014ZFOURGE
     !!{
     Constructors for the \cite{tomczak_galaxy_2014} survey geometry class.
     !!}
     module procedure tomczak2014ZFOURGEConstructorParameters
     module procedure tomczak2014ZFOURGEConstructorInternal
  end interface surveyGeometryTomczak2014ZFOURGE

  ! Paths and file names for mangle polygon files.
  integer, parameter :: countFields         =   2

  ! Angular power spectra.
  integer, parameter :: angularPowerMaximumL=5000

contains

  function tomczak2014ZFOURGEConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the \cite{tomczak_galaxy_2014} conditional mass function class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (surveyGeometryTomczak2014ZFOURGE)                :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    integer                                                  :: redshiftBin

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <inputParameter>
      <name>redshiftBin</name>
      <source>parameters</source>
      <description>The redshift bin (0, 1, 2, 3, 4, 5, 6, or 7) of the \cite{tomczak_galaxy_2014} mass function to use.</description>
    </inputParameter>
    !!]
    self=surveyGeometryTomczak2014ZFOURGE(redshiftBin,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function tomczak2014ZFOURGEConstructorParameters

  function tomczak2014ZFOURGEConstructorInternal(redshiftBin,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \cite{tomczak_galaxy_2014} mass function class.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Error                      , only : Error_Report
    implicit none
    type   (surveyGeometryTomczak2014ZFOURGE)                        :: self
    integer                                  , intent(in   )         :: redshiftBin
    class  (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="redshiftBin, *cosmologyFunctions_"/>
    !!]

    ! Find distance limits for this redshift bin.
    select case (redshiftBin)
    case(0)
       self%redshiftMinimum=0.20d0
       self%redshiftMaximum=0.50d0
    case(1)
       self%redshiftMinimum=0.50d0
       self%redshiftMaximum=0.75d0
    case(2)
       self%redshiftMinimum=0.75d0
       self%redshiftMaximum=1.00d0
    case(3)
       self%redshiftMinimum=1.00d0
       self%redshiftMaximum=1.25d0
    case(4)
       self%redshiftMinimum=1.25d0
       self%redshiftMaximum=1.50d0
    case(5)
       self%redshiftMinimum=1.50d0
       self%redshiftMaximum=2.00d0
    case(6)
       self%redshiftMinimum=2.00d0
       self%redshiftMaximum=2.50d0
    case(7)
       self%redshiftMinimum=2.50d0
       self%redshiftMaximum=3.00d0
    case default
       call Error_Report('0≤redshiftBin≤7 is required'//{introspection:location})
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
    call self%initialize()
    return
  end function tomczak2014ZFOURGEConstructorInternal

  subroutine tomczak2014ZFOURGEDestructor(self)
    !!{
    Destructor for the \refClass{surveyGeometryTomczak2014ZFOURGE} survey geometry class.
    !!}
    implicit none
    type(surveyGeometryTomczak2014ZFOURGE), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine tomczak2014ZFOURGEDestructor

  integer function tomczak2014ZFOURGEFieldCount(self)
    !!{
    Return the number of fields in this sample.
    !!}
    implicit none
    class(surveyGeometryTomczak2014ZFOURGE), intent(inout) :: self
    !$GLC attributes unused :: self

    tomczak2014ZFOURGEFieldCount=countFields
    return
  end function tomczak2014ZFOURGEFieldCount

  double precision function tomczak2014ZFOURGEDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the minimum distance at which a galaxy is included.
    !!}
    implicit none
    class           (surveyGeometryTomczak2014ZFOURGE), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                         luminosity, starFormationRate
    integer                                           , intent(in   ), optional :: field
    !$GLC attributes unused :: mass, field, magnitudeAbsolute, luminosity, starFormationRate

    tomczak2014ZFOURGEDistanceMinimum=self%binDistanceMinimum
    return
  end function tomczak2014ZFOURGEDistanceMinimum

  double precision function tomczak2014ZFOURGEDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Error                      , only : Error_Report
    implicit none
    class           (surveyGeometryTomczak2014ZFOURGE), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                         luminosity, starFormationRate
    integer                                           , intent(in   ), optional :: field
    double precision                                                            :: redshift  , logarithmicMass

    ! Validate field.
    if (.not.present(field)) call Error_Report('field must be specified'//{introspection:location})
    ! Validate arguments.
    if (present(magnitudeAbsolute)) call Error_Report('`magnitudeAbsolute` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! Find the limiting redshift for this mass. (See
    ! constraints/dataAnalysis/stellarMassFunctions_ZFOURGE_z0.2_2.5/massRedshiftRelation.pl for details.)
    if (present(mass)) then
       logarithmicMass=log10(mass)
       select case (field)
       case (1)
          redshift=-114.659703302477d0+logarithmicMass*(45.9008293873578d0+logarithmicMass*(-6.16172903321726d0+logarithmicMass*(0.278223082977791d0)))
       case (2)
          redshift=-58.4827675323933d0+logarithmicMass*(20.2501133330335d0+logarithmicMass*(-2.35628286307041d0+logarithmicMass*(0.0927047000361006d0)))
       case default
          redshift=0.0d0
          call Error_Report('1 ≤ field ≤ 2 required'//{introspection:location})
       end select
       ! Convert from redshift to comoving distance.
       tomczak2014ZFOURGEDistanceMaximum                                                                                              &
            &=self%cosmologyFunctions_%distanceComovingConvert(                                                                       &
            &                                                  output  =distanceTypeComoving                                        , &
            &                                                  redshift=min(max(redshift,self%redshiftMinimum),self%redshiftMaximum)  &
            &                                                 )
       ! Limit the maximum distance.
       tomczak2014ZFOURGEDistanceMaximum=min(tomczak2014ZFOURGEDistanceMaximum,self%binDistanceMaximum)
    else
       tomczak2014ZFOURGEDistanceMaximum=self%binDistanceMaximum
    end if
    return
  end function tomczak2014ZFOURGEDistanceMaximum

  double precision function tomczak2014ZFOURGEVolumeMaximum(self,mass,field)
    !!{
    Compute the maximum volume within which a galaxy is visible.
    !!}
    implicit none
    class           (surveyGeometryTomczak2014ZFOURGE), intent(inout)           :: self
    double precision                                  , intent(in   )           :: mass
    integer                                           , intent(in   ), optional :: field

     ! Compute the volume.
    tomczak2014ZFOURGEVolumeMaximum                              &
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
  end function tomczak2014ZFOURGEVolumeMaximum

  function tomczak2014ZFOURGEMangleDirectory(self)
    !!{
    Return the path to the directory containing \gls{mangle} files.
    !!}
    use :: Input_Paths, only : inputPath, pathTypeDataStatic
    implicit none
    class(surveyGeometryTomczak2014ZFOURGE), intent(inout) :: self
    type (varying_string                  )                :: tomczak2014ZFOURGEMangleDirectory
    !$GLC attributes unused :: self

    tomczak2014ZFOURGEMangleDirectory=inputPath(pathTypeDataStatic)//"surveyGeometry/ZFOURGE/"
    return
  end function tomczak2014ZFOURGEMangleDirectory

  subroutine tomczak2014ZFOURGEMangleFiles(self,mangleFiles)
    !!{
    Return a list of \gls{mangle} files.
    !!}
    implicit none
    class(surveyGeometryTomczak2014ZFOURGE)                           , intent(inout) :: self
    type (varying_string                  ), allocatable, dimension(:), intent(inout) :: mangleFiles

    allocate(mangleFiles(2))
    mangleFiles(1)="+"//self%mangleDirectory()//"/ZFOURGE-CDFS.ply:"  // &
         &         "+"//self%mangleDirectory()//"/ZFOURGE-COSMOS.ply:"// &
         &         "+"//self%mangleDirectory()//"/ZFOURGE-UDS.ply"
    mangleFiles(2)="+"//self%mangleDirectory()//"/NMBS-COSMOS.ply:"   // &
         &         "+"//self%mangleDirectory()//"/NMBS-AEGIS.ply"
    return
  end subroutine tomczak2014ZFOURGEMangleFiles

  integer function tomczak2014ZFOURGEAngularPowerMaximumDegree(self)
    !!{
    Return the maximum degree for which angular power is computed for the \cite{tomczak_galaxy_2014} survey.
    !!}
    implicit none
    class(surveyGeometryTomczak2014ZFOURGE), intent(inout) :: self
    !$GLC attributes unused :: self

    tomczak2014ZFOURGEAngularPowerMaximumDegree=angularPowerMaximumL
    return
  end function tomczak2014ZFOURGEAngularPowerMaximumDegree

