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
Implements the geometry of the ULTRAVISTA survey used by \cite{muzzin_evolution_2013}.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <surveyGeometry name="surveyGeometryMuzzin2013ULTRAVISTA">
   <description>
    A survey geometry class that describes the survey geometry of \cite{muzzin_evolution_2013}. 
    
    For the angular mask, we generate a \gls{mangle} polygon file, by first defining a rectangle encompassing the bounds of the
    ULTRAVISTA field ($149.373^\circ &lt; \alpha &lt; 150.779^\circ$ and $1.604^\circ &lt; \delta &lt; 2.81^\circ$). From this rectangle, we
    then remove circles of radii $75^{\prime\prime}$ around bright stars (i.e. those bright than 10$^\mathrm{th}$ and
    $8^\mathrm{th}$ magnitudes in the USNO and 2MASS star lists respectively) and radii $30^{\prime\prime}$ around medium stars
    (i.e. those bright than $13^\mathrm{th}$ and $10.5^\mathrm{th}$ magnitudes in the USNO and 2MASS star lists
    respectively). Finally, we mask regions of one detector for which 75\% of pixels are dead by clipping pixels with weights
    below $0.02$ in the K$_\mathrm{s}$-band weight map. These choices match those made in the ULTRAVISTA survey (A.~Muzzin,
    private communication). The solid angle of each mask is computed using the \gls{mangle} {\normalfont \ttfamily harmonize}
    command.
    
    To determine the depth as a function of stellar mass, we simply fit the
    \href{https://github.com/galacticusorg/datasets/blob/master/static/surveyGeometry/ULTRAVISTA/Mstar_redshift_completeness_emp_uvista_v4.1_100.dat}{tabulated
    relations} provided by the ULTRAVISTA survey:
    \begin{equation}
    z_\mathrm{max}(M_\star) = {-8364.45 + m (4331.82 + m (-896.596 + m (92.6999 + m (-4.78750 + m (0.0988215))))) \over 1 -
    \exp[(m-11.24)/0.02] }
     \label{eq:MuzzinDepthPolynomial}
    \end{equation}
    where $m= \log_{10}(M_\star/M_\odot)$.
    
    \begin{figure}
     \begin{center}
     \includegraphics[width=85mm,trim=0mm 0mm 0mm 4mm,clip]{Plots/DataAnalysis/MuzzinULTRAVISTAMassRedshiftRelation.pdf}
     \end{center}
     \caption{The maximum distance at which a galaxy of given stellar mass can be detected in the sample of
     \protect\cite{muzzin_evolution_2013}. The dotted line shows the results obtained from the ULTRAVISTA survey
     \protect\citep{muzzin_evolution_2013}, while the solid line shows the polynomial fit to these results (given in
     eqn.~\ref{eq:MuzzinDepthPolynomial}).}
     \label{fig:MuzzinULTRAVISTADepthFit}
    \end{figure}
   </description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryMangle) :: surveyGeometryMuzzin2013ULTRAVISTA
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     integer                                            :: redshiftBin
     double precision                                   :: binDistanceMinimum , binDistanceMaximum
   contains
     final     ::                              muzzin2013ULTRAVISTADestructor
     procedure :: fieldCount                => muzzin2013ULTRAVISTAFieldCount
     procedure :: distanceMinimum           => muzzin2013ULTRAVISTADistanceMinimum
     procedure :: distanceMaximum           => muzzin2013ULTRAVISTADistanceMaximum
     procedure :: volumeMaximum             => muzzin2013ULTRAVISTAVolumeMaximum
     procedure :: angularPowerMaximumDegree => muzzin2013ULTRAVISTAAngularPowerMaximumDegree
     procedure :: mangleDirectory           => muzzin2013ULTRAVISTAMangleDirectory
     procedure :: mangleFiles               => muzzin2013ULTRAVISTAMangleFiles
   end type surveyGeometryMuzzin2013ULTRAVISTA

  interface surveyGeometryMuzzin2013ULTRAVISTA
     !!{
     Constructors for the \cite{muzzin_evolution_2013} survey geometry class.
     !!}
     module procedure muzzin2013ULTRAVISTAConstructorParameters
     module procedure muzzin2013ULTRAVISTAConstructorInternal
  end interface surveyGeometryMuzzin2013ULTRAVISTA

  ! Paths and file names for mangle polygon files.
  integer, parameter :: countFields         =   1

  ! Angular power spectra.
  integer, parameter :: angularPowerMaximumL=3600

contains

  function muzzin2013ULTRAVISTAConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the \cite{muzzin_evolution_2013} conditional mass function class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(surveyGeometryMuzzin2013ULTRAVISTA) :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    integer                                                  :: redshiftBin

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <inputParameter>
      <name>redshiftBin</name>
      <source>parameters</source>
      <description>The redshift bin (0, 1, 2, 3, 4, 5, or 6) of the \cite{muzzin_evolution_2013} mass function to use.</description>
    </inputParameter>
    !!]
    self=surveyGeometryMuzzin2013ULTRAVISTA(redshiftBin,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function muzzin2013ULTRAVISTAConstructorParameters

  function muzzin2013ULTRAVISTAConstructorInternal(redshiftBin,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \cite{muzzin_evolution_2013} mass function class.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Error                      , only : Error_Report
    implicit none
    type            (surveyGeometryMuzzin2013ULTRAVISTA)                        :: self
    integer                                             , intent(in   )         :: redshiftBin
    class           (cosmologyFunctionsClass           ), intent(in   ), target :: cosmologyFunctions_
    double precision                                                            :: redshiftMinimum    , redshiftMaximum
    !![
    <constructorAssign variables="redshiftBin, *cosmologyFunctions_"/>
    !!]

    ! Find distance limits for this redshift bin.
    select case (redshiftBin)
    case(0)
       redshiftMinimum=0.20d0
       redshiftMaximum=0.50d0
    case(1)
       redshiftMinimum=0.50d0
       redshiftMaximum=1.00d0
    case(2)
       redshiftMinimum=1.00d0
       redshiftMaximum=1.50d0
    case(3)
       redshiftMinimum=1.50d0
       redshiftMaximum=2.00d0
    case(4)
       redshiftMinimum=2.00d0
       redshiftMaximum=2.50d0
    case(5)
       redshiftMinimum=2.50d0
       redshiftMaximum=3.00d0
    case(6)
       redshiftMinimum=3.00d0
       redshiftMaximum=4.00d0
    case default
       call Error_Report('0≤redshiftBin≤6 is required'//{introspection:location})
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
  end function muzzin2013ULTRAVISTAConstructorInternal

  subroutine muzzin2013ULTRAVISTADestructor(self)
    !!{
    Destructor for the \refClass{surveyGeometryMuzzin2013ULTRAVISTA} survey geometry class.
    !!}
    implicit none
    type(surveyGeometryMuzzin2013ULTRAVISTA), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine muzzin2013ULTRAVISTADestructor

  integer function muzzin2013ULTRAVISTAFieldCount(self)
    !!{
    Return the number of fields in this sample.
    !!}
    implicit none
    class(surveyGeometryMuzzin2013ULTRAVISTA), intent(inout) :: self
    !$GLC attributes unused :: self

    muzzin2013ULTRAVISTAFieldCount=countFields
    return
  end function muzzin2013ULTRAVISTAFieldCount

  double precision function muzzin2013ULTRAVISTADistanceMinimum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the minimum distance at which a galaxy is included.
    !!}
    implicit none
    class           (surveyGeometryMuzzin2013ULTRAVISTA), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                           luminosity, starFormationRate
    integer                                             , intent(in   ), optional :: field
    !$GLC attributes unused :: mass, field, magnitudeAbsolute, luminosity, starFormationRate

    muzzin2013ULTRAVISTADistanceMinimum=self%binDistanceMinimum
    return
  end function muzzin2013ULTRAVISTADistanceMinimum

  double precision function muzzin2013ULTRAVISTADistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    implicit none
    class           (surveyGeometryMuzzin2013ULTRAVISTA), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                           luminosity, starFormationRate
    integer                                             , intent(in   ), optional :: field
    double precision                                                              :: redshift  , logarithmicMass
    !$GLC attributes unused :: fieldK magnitudeAbsolute, luminosity, starFormationRate

    ! Validate arguments.
    if (present(magnitudeAbsolute)) call Error_Report('`magnitudeAbsolute` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! Find the limiting redshift for this mass. (See
    ! constraints/dataAnalysis/stellarMassFunctions_ULTRAVISTA_z0.2_4.0/massRedshiftRelation.pl for details.)
    if (present(mass)) then
       logarithmicMass=log10(mass)
       redshift=-6076.22869161192d0+logarithmicMass*(3231.43806947672d0+logarithmicMass*(-686.81594922437d0+logarithmicMass*(72.9147764759627d0+logarithmicMass*(-3.86638121388152d0+logarithmicMass*(0.0819398410572916d0)))))
       if (logarithmicMass < 11.24d0) then
          redshift=min(4.0d0,redshift/(1.0d0-exp((logarithmicMass-11.24d0)/0.02d0)))
       else
          redshift=    4.0d0
       end if
       if (redshift < 0.0d0) then
          muzzin2013ULTRAVISTADistanceMaximum=0.0d0
       else
          ! Convert from redshift to comoving distance.
          muzzin2013ULTRAVISTADistanceMaximum                                                    &
               &=self%cosmologyFunctions_%distanceComovingConvert(                               &
               &                                                  output  =distanceTypeComoving, &
               &                                                  redshift=redshift              &
               &                                                 )
          ! Limit the maximum distance.
          muzzin2013ULTRAVISTADistanceMaximum=min(muzzin2013ULTRAVISTADistanceMaximum,self%binDistanceMaximum)
       end if
    else
       muzzin2013ULTRAVISTADistanceMaximum   =                                        self%binDistanceMaximum
    end if
    return
  end function muzzin2013ULTRAVISTADistanceMaximum

  double precision function muzzin2013ULTRAVISTAVolumeMaximum(self,mass,field)
    !!{
    Compute the maximum volume within which a galaxy is visible.
    !!}
    implicit none
    class           (surveyGeometryMuzzin2013ULTRAVISTA), intent(inout)           :: self
    double precision                                    , intent(in   )           :: mass
    integer                                             , intent(in   ), optional :: field

    ! Compute the volume.
    muzzin2013ULTRAVISTAVolumeMaximum                            &
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
  end function muzzin2013ULTRAVISTAVolumeMaximum

  function muzzin2013ULTRAVISTAMangleDirectory(self)
    !!{
    Return the path to the directory containing \gls{mangle} files.
    !!}
    use :: Input_Paths, only : inputPath, pathTypeDataStatic
    implicit none
    class(surveyGeometryMuzzin2013ULTRAVISTA), intent(inout) :: self
    type (varying_string                    )                :: muzzin2013ULTRAVISTAMangleDirectory
    !$GLC attributes unused :: self

    muzzin2013ULTRAVISTAMangleDirectory=inputPath(pathTypeDataStatic)//"surveyGeometry/ULTRAVISTA/"
    return
  end function muzzin2013ULTRAVISTAMangleDirectory

  subroutine muzzin2013ULTRAVISTAMangleFiles(self,mangleFiles)
    !!{
    Return a list of \gls{mangle} files.
    !!}
    implicit none
    class(surveyGeometryMuzzin2013ULTRAVISTA)                           , intent(inout) :: self
    type (varying_string                    ), allocatable, dimension(:), intent(inout) :: mangleFiles
    !$GLC attributes unused :: self

    allocate(mangleFiles(1))
    mangleFiles(1)='surveyMask.ply'
    return
  end subroutine muzzin2013ULTRAVISTAMangleFiles

  integer function muzzin2013ULTRAVISTAAngularPowerMaximumDegree(self)
    !!{
    Return the maximum degree for which angular power is computed for the \cite{muzzin_evolution_2013} survey.
    !!}
    implicit none
    class(surveyGeometryMuzzin2013ULTRAVISTA), intent(inout) :: self
    !$GLC attributes unused :: self

    muzzin2013ULTRAVISTAAngularPowerMaximumDegree=angularPowerMaximumL
    return
  end function muzzin2013ULTRAVISTAAngularPowerMaximumDegree

