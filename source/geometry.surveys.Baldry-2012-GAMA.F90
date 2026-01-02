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
Implements the geometry of the GAMA survey used by \cite{baldry_galaxy_2012}.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <surveyGeometry name="surveyGeometryBaldry2012GAMA">
   <description>
    A survey geometry class that describes the survey geometry of \cite{baldry_galaxy_2012}. 
    
    For the angular mask we use the specifications of the G09, G12, and G15 fields given by \cite{driver_galaxy_2011} to
    construct {\normalfont \scshape mangle} polygon files.
    
    To determine the depth as a function of stellar mass, we make use of the publicly available tabulated mass function,
    $\phi$, and number of galaxies per bin, $N$. The effective volume of each bin is found as $V_i =
    N_i/\phi_i\Delta\log_{10}M_\star$, where $\Delta\log_{10}M_\star$ is the width of the bin. The GAMA survey consists of
    three fields, each of the same solid angle, but with differing depths. We assume that the relative depths in terms of
    stellar mass scale with the depth in terms of flux. Given this assumption, these volumes are converted to maximum distances
    in each field using the solid angle quoted above. The resulting mass vs. distance relation in each field is fit with a
    $1^\mathrm{st}$-order polynomial in log-log space over the range where the maximum volume is limited by the survey depth
    and not by the imposed $z=0.06$ upper limit to redshift. Figure~\ref{fig:BaldryGAMADepthFit} shows the resulting relation
    between stellar mass and the maximum distance at which such a galaxy would be included in the sample. Points indicate
    results from GAMA, while the line shows a polynomial fit:
    \begin{equation}
     \log_{10} \left[ {D_\mathrm{max}(M_\star) \over \hbox{Mpc}}\right] = \left\{ \begin{array}{ll} -0.521 + 0.319m &amp;
     \hbox{fields G09/G15} \\ -0.361 + 0.319m &amp; \hbox{field G12} \end{array} \right.
     \label{eq:BaldryDepthPolynomial}
    \end{equation}
    where $m= \log_{10}(M_\star/M_\odot)$. We use this polynomial fit to determine the depth of the sample as a function of
    stellar mass.
    
    \begin{figure}
     \begin{center}
     \includegraphics[width=85mm,trim=0mm 0mm 0mm 4mm,clip]{Plots/DataAnalysis/BaldryGAMAMassDistanceRelation.pdf}
     \end{center}
     \caption{The maximum distance at which a galaxy of given stellar mass can be detected in the sample of
     \protect\cite{baldry_galaxy_2012}. Points show the results obtained from data provided by Baldry, while the lines shows a
     polynomial fit to these results (given in eqn.~\ref{eq:BaldryDepthPolynomial}). Note that above $10^9M_\odot$ the distance
     is limited by the imposed upper limit of $z=0.06$ in the GAMA sample---the polynomial fit does not consider these points.}
     \label{fig:BaldryGAMADepthFit}
    \end{figure}
   </description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryMangle) :: surveyGeometryBaldry2012GAMA
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_   => null()
     double precision                                   :: distanceMaximumSurvey
   contains
     final     ::                              baldry2012GAMADestructor
     procedure :: fieldCount                => baldry2012GAMAFieldCount
     procedure :: distanceMaximum           => baldry2012GAMADistanceMaximum
     procedure :: angularPowerMaximumDegree => baldry2012GAMAAngularPowerMaximumDegree
     procedure :: mangleDirectory           => baldry2012GAMAMangleDirectory
     procedure :: mangleFiles               => baldry2012GAMAMangleFiles
  end type surveyGeometryBaldry2012GAMA

  interface surveyGeometryBaldry2012GAMA
     !!{
     Constructors for the \cite{baldry_galaxy_2012} survey geometry class.
     !!}
     module procedure baldry2012GAMAConstructorParameters
     module procedure baldry2012GAMAConstructorInternal
  end interface surveyGeometryBaldry2012GAMA

  ! Number of fields.
  integer, parameter :: countFields         =  3

  ! Maximum degree for angular power spectrum
  integer, parameter :: angularPowerMaximumL=360

contains

  function baldry2012GAMAConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \cite{baldry_galaxy_2012} conditional mass function class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass
    use :: Input_Parameters   , only : inputParameter    , inputParameters
    implicit none
    type (surveyGeometryBaldry2012GAMA)                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(cosmologyFunctionsClass     ), pointer       :: cosmologyFunctions_

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    ! Build the object.
    self=surveyGeometryBaldry2012GAMA(cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function baldry2012GAMAConstructorParameters

  function baldry2012GAMAConstructorInternal(cosmologyFunctions_) result (self)
    !!{
    Internal constructor for the \cite{baldry_galaxy_2012} conditional mass function class.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    implicit none
    type            (surveyGeometryBaldry2012GAMA)                        :: self
    class           (cosmologyFunctionsClass     ), intent(in   ), target :: cosmologyFunctions_
    double precision                              , parameter             :: redshiftMaximum    =0.06d0
    !![
    <constructorAssign variables="*cosmologyFunctions_"/>
    !!]

    call self%initialize()
    self%distanceMaximumSurvey=self%cosmologyFunctions_%distanceComovingConvert(distanceTypeComoving,redshift=redshiftMaximum)
   return
  end function baldry2012GAMAConstructorInternal

  subroutine baldry2012GAMADestructor(self)
    !!{
    Destructor for the \refClass{surveyGeometryBaldry2012GAMA} geometry class.
    !!}
    implicit none
    type(surveyGeometryBaldry2012GAMA), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine baldry2012GAMADestructor

  integer function baldry2012GAMAFieldCount(self)
    !!{
    Return the number of fields in this sample.
    !!}
    implicit none
    class(surveyGeometryBaldry2012GAMA), intent(inout) :: self
    !$GLC attributes unused :: self

    baldry2012GAMAFieldCount=countFields
    return
  end function baldry2012GAMAFieldCount

  double precision function baldry2012GAMADistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (surveyGeometryBaldry2012GAMA), intent(inout)           :: self
    double precision                              , intent(in   ), optional :: mass           , magnitudeAbsolute, &
         &                                                                     luminosity     , starFormationRate
    integer                                       , intent(in   ), optional :: field
    double precision                                                        :: logarithmicMass

    ! Validate field.
    if (.not.present(field)) call Error_Report('field must be specified'//{introspection:location})
    ! Validate arguments.
    if (present(magnitudeAbsolute)) call Error_Report('`magnitudeAbsolute` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! Compute the limiting distance.
    if (present(mass)) then
       logarithmicMass=log10(mass)
       select case (field)
       case (1,3) ! Fields G09 and G15.
          baldry2012GAMADistanceMaximum                          &
               & =10.0d0**(                                      &
               &           -0.521147071716417d0                  &
               &           +0.318557607893107d0*logarithmicMass  &
               &          )
       case (2)
          baldry2012GAMADistanceMaximum                          &
               & =10.0d0**(                                      &
               &           -0.361147071716369d0                  &
               &           +0.318557607893101d0*logarithmicMass  &
               &          )
       case default
          baldry2012GAMADistanceMaximum=0.0d0
          call Error_Report('1 ≤ field ≤ 3 required'//{introspection:location})
       end select
       baldry2012GAMADistanceMaximum=min(baldry2012GAMADistanceMaximum,self%distanceMaximumSurvey)
    else
       baldry2012GAMADistanceMaximum=self%distanceMaximumSurvey
    end if
    return
  end function baldry2012GAMADistanceMaximum

  function baldry2012GAMAMangleDirectory(self)
    !!{
    Return the path to the directory containing \gls{mangle} files.
    !!}
    use :: Input_Paths, only : inputPath, pathTypeDataStatic
    implicit none
    class(surveyGeometryBaldry2012GAMA), intent(inout) :: self
    type (varying_string              )                :: baldry2012GAMAMangleDirectory
    !$GLC attributes unused :: self

    baldry2012GAMAMangleDirectory=inputPath(pathTypeDataStatic)//"surveyGeometry/GAMA/"
    return
  end function baldry2012GAMAMangleDirectory

  subroutine baldry2012GAMAMangleFiles(self,mangleFiles)
    !!{
    Return a list of \gls{mangle} files.
    !!}
    implicit none
    class(surveyGeometryBaldry2012GAMA)                           , intent(inout) :: self
    type (varying_string              ), allocatable, dimension(:), intent(inout) :: mangleFiles

    allocate(mangleFiles(3))
    mangleFiles=                                                   &
         &      [                                                  &
         &       self%mangleDirectory()//"angularGeometryG09.ply", &
         &       self%mangleDirectory()//"angularGeometryG12.ply", &
         &       self%mangleDirectory()//"angularGeometryG15.ply"  &
         &      ]
    return
  end subroutine baldry2012GAMAMangleFiles

  integer function baldry2012GAMAAngularPowerMaximumDegree(self)
    !!{
    Return the maximum degree for which angular power is computed for the \cite{bernardi_massive_2013} survey.
    !!}
    implicit none
    class(surveyGeometryBaldry2012GAMA), intent(inout) :: self
    !$GLC attributes unused :: self

    baldry2012GAMAAngularPowerMaximumDegree=angularPowerMaximumL
    return
  end function baldry2012GAMAAngularPowerMaximumDegree

