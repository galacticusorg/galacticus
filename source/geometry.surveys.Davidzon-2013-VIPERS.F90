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
Implements the geometry of the VIPERS survey used by \cite{davidzon_vimos_2013}.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <surveyGeometry name="surveyGeometryDavidzon2013VIPERS">
   <description>
    A survey geometry class that describes the survey geometry of \cite{davidzon_vimos_2013}. 
    
    For the angular mask, we make use of {\normalfont \ttfamily mangle} polygon files provided by I.~Davidzon (private
    communication) corresponding to the VIPERS fields. The solid angle of each mask is computed using the {\normalfont
    \ttfamily mangle} {\normalfont \ttfamily harmonize} command.
    
    To determine the depth as a function of stellar mass, we make use of the tabulated mass function, $\phi$, and number of
    galaxies per bin, $N$, supplied by I.~Davidzon (private communication). The effective volume of each bin is found as $V_i =
    N_i/f_\mathrm{complete}\phi_i\Delta\log_{10}M_\star$, where $\Delta\log_{10}M_\star$ is the width of the bin, and
    $f_\mathrm{complete}$ is the completeness of the survey, estimated to be approximately 40\% \citep{guzzo_vimos_2013}. These
    volumes are converted to maximum distances in each field using the survey solid angle. The resulting mass vs. distance
    relation in each field is fit with a $1^\mathrm{st}$-order polynomial in log-log space over the range where the maximum
    volume is limited by the survey depth and not by the imposed upper limit to redshift. Figure~\ref{fig:Davidzon2013DepthFit}
    shows the resulting relation between stellar mass and the maximum distance at which such a galaxy would be included in the
    sample. Points indicate results from VIPERS, while the lines show polynomial fits:
    \begin{equation}
     \log_{10} \left[ {D_\mathrm{max}(M_\star) \over \hbox{Mpc}}\right] = \left\{ \begin{array}{ll} 3.207 + 0.0124m &amp; 0.5 &lt; z &lt;
     0.6 \\ 3.148 + 0.0268m &amp; 0.6 &lt; z &lt; 0.8 \\ 3.207 + 0.0273m &amp; 0.8 &lt; z &lt; 1.0 \end{array} \right.
     \label{eq:DavidzonDepthPolynomial}
    \end{equation}
    where $m= \log_{10}(M_\star/M_\odot)$. We use this polynomial fit to determine the depth of the sample as a function of
    stellar mass.
    
    \begin{figure}
     \begin{center}
     \includegraphics[width=85mm,trim=0mm 0mm 0mm 4mm,clip]{Plots/DataAnalysis/DavidzonVIPERSMassDistanceRelation.pdf}
     \end{center}
     \caption{The maximum distance at which a galaxy of given stellar mass can be detected in the sample of
     \protect\cite{davidzon_vimos_2013}. Points show the results obtained from data provided by Davidzon, while the lines shows
     a polynomial fit to these results (given in eqn.~\ref{eq:DavidzonDepthPolynomial}). Note that at high masses the distance
     is limited by the imposed upper limit---the polynomial fit does not consider these points.}
     \label{fig:Davidzon2013DepthFit}
    \end{figure}
   </description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryMangle) :: surveyGeometryDavidzon2013VIPERS
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     integer                                            :: redshiftBin
     double precision                                   :: binDistanceMinimum           , binDistanceMaximum
   contains
     final     ::                              davidzon2013VIPERSDestructor
     procedure :: fieldCount                => davidzon2013VIPERSFieldCount
     procedure :: distanceMinimum           => davidzon2013VIPERSDistanceMinimum
     procedure :: distanceMaximum           => davidzon2013VIPERSDistanceMaximum
     procedure :: volumeMaximum             => davidzon2013VIPERSVolumeMaximum
     procedure :: angularPowerMaximumDegree => davidzon2013VIPERSAngularPowerMaximumDegree
     procedure :: mangleDirectory           => davidzon2013VIPERSMangleDirectory
     procedure :: mangleFiles               => davidzon2013VIPERSMangleFiles
   end type surveyGeometryDavidzon2013VIPERS

  interface surveyGeometryDavidzon2013VIPERS
     !!{
     Constructors for the \cite{davidzon_vimos_2013} survey geometry class.
     !!}
     module procedure davidzon2013VIPERSConstructorParameters
     module procedure davidzon2013VIPERSConstructorInternal
  end interface surveyGeometryDavidzon2013VIPERS

  ! Number of fields.
  integer, parameter :: countFields         =  1

  ! Angular power spectra.
  integer, parameter :: angularPowerMaximumL=720

contains

  function davidzon2013VIPERSConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{davidzon_vimos_2013} conditional mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (surveyGeometryDavidzon2013VIPERS)                :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    integer                                                  :: redshiftBin

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <inputParameter>
      <name>redshiftBin</name>
      <source>parameters</source>
      <description>The redshift bin (0, 1, 2) of the \cite{davidzon_vimos_2013} mass function to use.</description>
    </inputParameter>
    !!]
    self=surveyGeometryDavidzon2013VIPERS(redshiftBin,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function davidzon2013VIPERSConstructorParameters

  function davidzon2013VIPERSConstructorInternal(redshiftBin,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \cite{davidzon_vimos_2013} mass function class.
    !!}
    use :: Cosmology_Functions        , only : cosmologyFunctionsClass
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Error                      , only : Error_Report
    implicit none
    type            (surveyGeometryDavidzon2013VIPERS)                        :: self
    integer                                           , intent(in   )         :: redshiftBin
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    double precision                                                          :: redshiftMinimum    , redshiftMaximum
    !![
    <constructorAssign variables="redshiftBin, *cosmologyFunctions_"/>
    !!]

    ! Find distance limits for this redshift bin.
    select case (redshiftBin)
    case(0)
       redshiftMinimum=0.50d0
       redshiftMaximum=0.60d0
    case(1)
       redshiftMinimum=0.60d0
       redshiftMaximum=0.80d0
    case(2)
       redshiftMinimum=0.80d0
       redshiftMaximum=1.00d0
    case default
       call Error_Report('0≤redshiftBin≤3 is required'//{introspection:location})
    end select
    self%binDistanceMinimum                                                                 &
         & =self%cosmologyFunctions_%distanceComovingConvert(                               &
         &                                                   output  =distanceTypeComoving, &
         &                                                   redshift=redshiftMinimum       &
         &                                                  )
    self%binDistanceMaximum                                                                 &
         & =self%cosmologyFunctions_%distanceComovingConvert(                               &
         &                                                   output  =distanceTypeComoving, &
         &                                                   redshift=redshiftMaximum       &
         &                                                  )
    call self%initialize()
    return
  end function davidzon2013VIPERSConstructorInternal

  subroutine davidzon2013VIPERSDestructor(self)
    !!{
    Destructor for the \refClass{surveyGeometryDavidzon2013VIPERS} survey geometry class.
    !!}
    implicit none
    type(surveyGeometryDavidzon2013VIPERS), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine davidzon2013VIPERSDestructor

  integer function davidzon2013VIPERSFieldCount(self)
    !!{
    Return the number of fields in this sample.
    !!}
    implicit none
    class(surveyGeometryDavidzon2013VIPERS), intent(inout) :: self
    !$GLC attributes unused :: self

    davidzon2013VIPERSFieldCount=countFields
    return
  end function davidzon2013VIPERSFieldCount

  double precision function davidzon2013VIPERSDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the minimum distance at which a galaxy is included.
    !!}
    implicit none
    class           (surveyGeometryDavidzon2013VIPERS), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                         luminosity, starFormationRate
    integer                                           , intent(in   ), optional :: field
    !$GLC attributes unused :: mass, field, magnitudeAbsolute, luminosity, starFormationRate

    davidzon2013VIPERSDistanceMinimum=self%binDistanceMinimum
    return
  end function davidzon2013VIPERSDistanceMinimum

  double precision function davidzon2013VIPERSDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (surveyGeometryDavidzon2013VIPERS), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: mass           , magnitudeAbsolute, &
         &                                                                         luminosity     , starFormationRate
    integer                                           , intent(in   ), optional :: field
    double precision                                                            :: logarithmicMass
    !$GLC attributes unused :: field

    ! Validate arguments.
    if (present(magnitudeAbsolute)) call Error_Report('`magnitudeAbsolute` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! Find the limiting distance for this mass. (See
    ! constraints/dataAnalysis/stellarMassFunctions_VIPERS_z0_1/massDistanceRelation.pl for details.)
    if (present(mass)) then
       logarithmicMass=log10(mass)
       select case (self%redshiftBin)
       case (0)
          davidzon2013VIPERSDistanceMaximum=3.20663737335189d0+logarithmicMass*(0.0124101908903665d0)
       case (1)
          davidzon2013VIPERSDistanceMaximum=3.14840402683405d0+logarithmicMass*(0.0268494389098537d0)
       case (2)
          davidzon2013VIPERSDistanceMaximum=3.20688538211015d0+logarithmicMass*(0.0273132827274515d0)
       case default
          davidzon2013VIPERSDistanceMaximum=0.0d0
          call Error_Report('invalid redshift bin'//{introspection:location})
       end select
       ! Limit the maximum distance.
       davidzon2013VIPERSDistanceMaximum=min(10.0d0**davidzon2013VIPERSDistanceMaximum,self%binDistanceMaximum)
    else
        davidzon2013VIPERSDistanceMaximum=                                             self%binDistanceMaximum
    end if
    return
  end function davidzon2013VIPERSDistanceMaximum

  double precision function davidzon2013VIPERSVolumeMaximum(self,mass,field)
    !!{
    Compute the maximum volume within which a galaxy is visible.
    !!}
    implicit none
    class           (surveyGeometryDavidzon2013VIPERS), intent(inout)           :: self
    double precision                                  , intent(in   )           :: mass
    integer                                           , intent(in   ), optional :: field

     ! Compute the volume.
    davidzon2013VIPERSVolumeMaximum                              &
         & =max(                                                 &
         &       0.0d0                                         , &
         &       self%solidAngle()                               &
         &      *(                                               &
         &        +self%distanceMaximum   (mass,field=field)**3  &
         &        -self%binDistanceMinimum                  **3  &
         &       )                                               &
         &      /3.0d0                                           &
         &     )
    return
  end function davidzon2013VIPERSVolumeMaximum

  function davidzon2013VIPERSMangleDirectory(self)
    !!{
    Return the path to the directory containing \gls{mangle} files.
    !!}
    use :: Input_Paths, only : inputPath, pathTypeDataStatic
    implicit none
    class(surveyGeometryDavidzon2013VIPERS), intent(inout) :: self
    type (varying_string                  )                :: davidzon2013VIPERSMangleDirectory
    !$GLC attributes unused :: self

    davidzon2013VIPERSMangleDirectory=inputPath(pathTypeDataStatic)//"surveyGeometry/VIPERS/"
    return
  end function davidzon2013VIPERSMangleDirectory

  subroutine davidzon2013VIPERSMangleFiles(self,mangleFiles)
    !!{
    Return a list of \gls{mangle} files.
    !!}
    implicit none
    class(surveyGeometryDavidzon2013VIPERS)                           , intent(inout) :: self
    type (varying_string                  ), allocatable, dimension(:), intent(inout) :: mangleFiles

    allocate(mangleFiles(1))
    mangleFiles(1)="+"//self%mangleDirectory()//"/maskCombinedBU.ply:"// &
         &         "-"//self%mangleDirectory()//"/photoW1.BU2.ply:"   // &
         &         "-"//self%mangleDirectory()//"/photoW4.BU2.ply"
    return
  end subroutine davidzon2013VIPERSMangleFiles

  integer function davidzon2013VIPERSAngularPowerMaximumDegree(self)
    !!{
    Return the maximum degree for which angular power is computed for the \cite{davidzon_vimos_2013} survey.
    !!}
    implicit none
    class(surveyGeometryDavidzon2013VIPERS), intent(inout) :: self
    !$GLC attributes unused :: self

    davidzon2013VIPERSAngularPowerMaximumDegree=angularPowerMaximumL
    return
  end function davidzon2013VIPERSAngularPowerMaximumDegree

