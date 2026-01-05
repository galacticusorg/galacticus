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
Implements the geometry of the PRIMUS survey used by \cite{moustakas_primus:_2013}.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <surveyGeometry name="surveyGeometryMoustakas2013PRIMUS">
   <description>
    A survey geometry class that describes the survey geometry of \cite{moustakas_primus:_2013}. 
    
    For the angular mask, we make use of \gls{mangle} polygon files provided by J.~Moustakas (private communication)
    corresponding to the give PRIMUS fields. The solid angle of each mask is computed using the \gls{mangle} {\normalfont
    \ttfamily harmonize} command.
    
    To determine the depth as a function of stellar mass, we make use of completeness limits for ``All'' galaxies given in
    Table~2 of \cite{moustakas_primus:_2013}. These are fit, for each field with a second order polynomial to give the limiting
    redshift as a function of stellar mass. Figure~\ref{fig:MoustakasPRIMUSDepthFit} shows the resulting relation between
    stellar mass and the maximum redshift at which such a galaxy would be included in the sample. Points indicate results from
    \cite{moustakas_primus:_2013}, while the line shows a polynomial fits:
    \begin{eqnarray}
    z_\mathrm{max}(M_\star) = +3.51+m(-0.941+m(+0.0651)) &amp; &amp; \hbox{COSMOS} \\
    z_\mathrm{max}(M_\star) = +2.46+m(-0.730+m(+0.0542)) &amp; &amp; \hbox{XMM-SXDS} \\
    z_\mathrm{max}(M_\star) = -3.60+m(+0.500+m(-0.0078)) &amp; &amp; \hbox{XMM-CFHTLS} \\
    z_\mathrm{max}(M_\star) = +5.87+m(-1.528+m(+0.0982)) &amp; &amp; \hbox{CDFS} \\
    z_\mathrm{max}(M_\star) = +6.87+m(-1.656+m(+0.1003)) &amp; &amp; \hbox{ELAIS-S1}
     \label{eq:MoustakasDepthPolynomial}
    \end{eqnarray}
    where $m= \log_{10}(M_\star/M_\odot)$. We use this polynomial fit to determine the depth of the sample as a function of
    stellar mass.
    
    \begin{figure}
     \begin{center}
     \includegraphics[width=85mm,trim=0mm 0mm 0mm 4mm,clip]{Plots/DataAnalysis/MoustakasPRIMUSMassRedshiftRelation.pdf}
     \end{center}
     \caption{The maximum distance at which a galaxy of given stellar mass can be detected in the sample of
     \protect\cite{moustakas_primus:_2013}. Points show the results obtained from completeness limit data taken from Table~2 of
     \protect\cite{moustakas_primus:_2013}, while the lines shows a polynomial fit to these results (given in
     eqn.~\ref{eq:MoustakasDepthPolynomial}).}
     \label{fig:MoustakasPRIMUSDepthFit}
    \end{figure}
   </description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryMangle) :: surveyGeometryMoustakas2013PRIMUS
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     integer                                            :: redshiftBin
     double precision                                   :: binDistanceMinimum           , binDistanceMaximum, &
          &                                                redshiftMinimum              , redshiftMaximum
   contains
     final     ::                              moustakas2013PRIMUSDestructor
     procedure :: fieldCount                => moustakas2013PRIMUSFieldCount
     procedure :: distanceMinimum           => moustakas2013PRIMUSDistanceMinimum
     procedure :: distanceMaximum           => moustakas2013PRIMUSDistanceMaximum
     procedure :: volumeMaximum             => moustakas2013PRIMUSVolumeMaximum
     procedure :: angularPowerMaximumDegree => moustakas2013PRIMUSAngularPowerMaximumDegree
     procedure :: mangleDirectory           => moustakas2013PRIMUSMangleDirectory
     procedure :: mangleFiles               => moustakas2013PRIMUSMangleFiles
  end type surveyGeometryMoustakas2013PRIMUS

  interface surveyGeometryMoustakas2013PRIMUS
     !!{
     Constructors for the \cite{moustakas_primus:_2013} survey geometry class.
     !!}
     module procedure moustakas2013PRIMUSConstructorParameters
     module procedure moustakas2013PRIMUSConstructorInternal
  end interface surveyGeometryMoustakas2013PRIMUS

  ! Number of fields.
  integer, parameter :: countFields         =   5

  ! Angular power spectra.
  integer, parameter :: angularPowerMaximumL=3000

contains

  function moustakas2013PRIMUSConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the \cite{moustakas_primus:_2013} conditional mass function class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (surveyGeometryMoustakas2013PRIMUS)                :: self
    type   (inputParameters                  ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    integer                                                   :: redshiftBin

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <inputParameter>
      <name>redshiftBin</name>
      <source>parameters</source>
      <description>The redshift bin (0, 1, 2, 3, 4, 5, or 5) of the \cite{moustakas_primus:_2013} mass function to use.</description>
    </inputParameter>
    !!]
    self=surveyGeometryMoustakas2013PRIMUS(redshiftBin,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function moustakas2013PRIMUSConstructorParameters

  function moustakas2013PRIMUSConstructorInternal(redshiftBin,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \cite{moustakas_primus:_2013} mass function class.
    !!}
    use :: Cosmology_Functions        , only : cosmologyFunctionsClass
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Error                      , only : Error_Report
    implicit none
    type   (surveyGeometryMoustakas2013PRIMUS)                        :: self
    integer                                   , intent(in   )         :: redshiftBin
    class  (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="redshiftBin, *cosmologyFunctions_"/>
    !!]

    ! Find distance limits for this redshift bin.
    select case (redshiftBin)
    case(0)
       self%redshiftMinimum=0.00d0
       self%redshiftMaximum=0.10d0
    case(1)
       self%redshiftMinimum=0.20d0
       self%redshiftMaximum=0.30d0
    case(2)
       self%redshiftMinimum=0.30d0
       self%redshiftMaximum=0.40d0
    case(3)
       self%redshiftMinimum=0.40d0
       self%redshiftMaximum=0.50d0
    case(4)
       self%redshiftMinimum=0.50d0
       self%redshiftMaximum=0.65d0
    case(5)
       self%redshiftMinimum=0.65d0
       self%redshiftMaximum=0.80d0
    case(6)
       self%redshiftMinimum=0.80d0
       self%redshiftMaximum=1.00d0
    case default
       call Error_Report('0≤redshiftBin≤6 is required'//{introspection:location})
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
  end function moustakas2013PRIMUSConstructorInternal

  subroutine moustakas2013PRIMUSDestructor(self)
    !!{
    Destructor for the \refClass{surveyGeometryMoustakas2013PRIMUS} survey geometry class.
    !!}
    implicit none
    type(surveyGeometryMoustakas2013PRIMUS), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine moustakas2013PRIMUSDestructor

  integer function moustakas2013PRIMUSFieldCount(self)
    !!{
    Return the number of fields in this sample.
    !!}
    implicit none
    class(surveyGeometryMoustakas2013PRIMUS), intent(inout) :: self
    !$GLC attributes unused :: self

    moustakas2013PRIMUSFieldCount=countFields
    return
  end function moustakas2013PRIMUSFieldCount

  double precision function moustakas2013PRIMUSDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the minimum distance at which a galaxy is included.
    !!}
    implicit none
    class           (surveyGeometryMoustakas2013PRIMUS), intent(inout)           :: self
    double precision                                   , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                          luminosity, starFormationRate
    integer                                            , intent(in   ), optional :: field
    !$GLC attributes unused :: mass, field, magnitudeAbsolute, luminosity, starFormationRate

    moustakas2013PRIMUSDistanceMinimum=self%binDistanceMinimum
    return
  end function moustakas2013PRIMUSDistanceMinimum

  double precision function moustakas2013PRIMUSDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Error                      , only : Error_Report
    implicit none
    class           (surveyGeometryMoustakas2013PRIMUS), intent(inout)           :: self
    double precision                                   , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                          luminosity, starFormationRate
    integer                                            , intent(in   ), optional :: field
    double precision                                                             :: redshift  , logarithmicMass

    ! Validate field.
    if (.not.present(field)) call Error_Report('field must be specified'//{introspection:location})
    ! Validate arguments.
    if (present(magnitudeAbsolute)) call Error_Report('`magnitudeAbsolute` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! Find the limiting redshift for this mass completeness limits from Moustakas et al. (2013; Table 2). (See
    ! constraints/dataAnalysis/stellarMassFunctions_PRIMUS_z0_1/massRedshiftRelation.pl for details.)
    if (present(mass)) then
       logarithmicMass=log10(mass)
       select case (field)
       case (1) ! COSMOS
          redshift=                  +3.51240871481968000d0  &
               &   +logarithmicMass*(-0.94131511297034200d0  &
               &   +logarithmicMass*(+0.06507208866075860d0) &
               &                                           )
       case (2) ! XMM-SXDS
          redshift=                  +2.46068289817352000d0  &
               &   +logarithmicMass*(-0.72960045705258400d0  &
               &   +logarithmicMass*(+0.05422457500058130d0) &
               &                                           )
       case (3) ! XMM-CFHTLS
          redshift=                  -3.60001396783385000d0  &
               &   +logarithmicMass*(+0.50007933123305700d0  &
               &   +logarithmicMass*(-0.00781044013508246d0) &
               &                                           )
       case (4) ! CDFS
          redshift=                  +5.86929723910240000d0  &
               &   +logarithmicMass*(-1.52816338306828000d0  &
               &   +logarithmicMass*(+0.09815249638377170d0) &
               &                                           )
       case (5) ! ELAIS-S1
          redshift=                  +6.87489619768651000d0  &
               &   +logarithmicMass*(-1.65556365363183000d0  &
               &   +logarithmicMass*(+0.10030052053225000d0) &
               &                                           )
       case default
          redshift=0.0d0
          call Error_Report('1 ≤ field ≤ 5 required'//{introspection:location})
       end select
       ! Convert from redshift to comoving distance.
       moustakas2013PRIMUSDistanceMaximum                                                                                             &
            &=self%cosmologyFunctions_%distanceComovingConvert(                                                                       &
            &                                                  output  =distanceTypeComoving                                        , &
            &                                                  redshift=max(min(redshift,self%redshiftMaximum),self%redshiftMinimum)  &
            &                                                 )
       ! Limit the maximum distance.
       moustakas2013PRIMUSDistanceMaximum=min(moustakas2013PRIMUSDistanceMaximum,self%binDistanceMaximum)
    else
       moustakas2013PRIMUSDistanceMaximum=                                       self%binDistanceMaximum
    end if
    return
  end function moustakas2013PRIMUSDistanceMaximum

  double precision function moustakas2013PRIMUSVolumeMaximum(self,mass,field)
    !!{
    Compute the maximum volume within which a galaxy is visible.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (surveyGeometryMoustakas2013PRIMUS), intent(inout)           :: self
    double precision                                   , intent(in   )           :: mass
    integer                                            , intent(in   ), optional :: field

    ! Validate field.
    if (.not.present(field)) call Error_Report('field must be specified'//{introspection:location})
    if (field < 1 .or. field > 5) call Error_Report('1 ≤ field ≤ 5 required'//{introspection:location})
    ! Compute the volume.
    moustakas2013PRIMUSVolumeMaximum                             &
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
  end function moustakas2013PRIMUSVolumeMaximum

  function moustakas2013PRIMUSMangleDirectory(self)
    !!{
    Return the path to the directory containing \gls{mangle} files.
    !!}
    use :: Input_Paths, only : inputPath, pathTypeDataStatic
    implicit none
    class(surveyGeometryMoustakas2013PRIMUS), intent(inout) :: self
    type (varying_string                   )                :: moustakas2013PRIMUSMangleDirectory
    !$GLC attributes unused :: self

    moustakas2013PRIMUSMangleDirectory=inputPath(pathTypeDataStatic)//"surveyGeometry/PRIMUS/"
    return
  end function moustakas2013PRIMUSMangleDirectory

  subroutine moustakas2013PRIMUSMangleFiles(self,mangleFiles)
    !!{
    Return a list of \gls{mangle} files.
    !!}
    implicit none
    class(surveyGeometryMoustakas2013PRIMUS)                           , intent(inout) :: self
    type (varying_string                   ), allocatable, dimension(:), intent(inout) :: mangleFiles

    allocate(mangleFiles(5))
    mangleFiles(1)=self%mangleDirectory()//"cosmos_field_galex_window_2mask.ply"
    mangleFiles(2)=self%mangleDirectory()//"xmm_swire_field_galex_window_2mask.ply"
    mangleFiles(3)=self%mangleDirectory()//"cfhtls_xmm_field_galex_window_2mask.ply"
    mangleFiles(4)=self%mangleDirectory()//"cdfs_field_galex_window_2mask.ply"
    mangleFiles(5)=self%mangleDirectory()//"es1_field_galex_window_2mask.ply"
    return
  end subroutine moustakas2013PRIMUSMangleFiles

  integer function moustakas2013PRIMUSAngularPowerMaximumDegree(self)
    !!{
    Return the maximum degree for which angular power is computed for the \cite{moustakas_primus:_2013} survey.
    !!}
    implicit none
    class(surveyGeometryMoustakas2013PRIMUS), intent(inout) :: self
    !$GLC attributes unused :: self

    moustakas2013PRIMUSAngularPowerMaximumDegree=angularPowerMaximumL
    return
  end function moustakas2013PRIMUSAngularPowerMaximumDegree

  integer function moustakas2013PRIMUSFieldPairIndex(i,j)
    !!{
    Compute the index of a pair of fields in the \cite{moustakas_primus:_2013} survey.
    !!}
    implicit none
    integer, intent(in   ) ::  i,  j
    integer                :: ii, jj

    ii=min(i,j)
    jj=max(i,j)
    moustakas2013PRIMUSFieldPairIndex=(ii-1)*(2*countFields-ii+2)/2+(jj-ii+1)
    return
  end function moustakas2013PRIMUSFieldPairIndex

