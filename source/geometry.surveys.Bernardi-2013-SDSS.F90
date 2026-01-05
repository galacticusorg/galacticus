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
Implements the geometry of the SDSS survey used by \cite{bernardi_massive_2013}.
!!}


  !![
  <surveyGeometry name="surveyGeometryBernardi2013SDSS">
   <description>
    A survey geometry class that describes the survey geometry of \cite{bernardi_massive_2013}. 
    
    For the angular mask, we make use of the \gls{mangle} polygon file provided by the \gls{mangle}
    project\footnote{Specifically,
    \href{https://zenodo.org/records/10998446/files/sdss_dr72safe0_res6d.pol.gz}{https://zenodo.org/records/10998446/files/sdss\_dr72safe0\_res6d.pol.gz}.}
    The solid angle of this mask, computed using the \gls{mangle} {\normalfont \ttfamily harmonize} command is
    2.232262776405~sr.
    
    To determine the depth as a function of stellar mass, we make use of results provided by M. Bernardi (private
    communication), giving the mean maximum volume, $V_\mathrm{max}$, as a function of stellar mass for galaxies in this
    sample. These maximum volumes are converted to maximum distances using the solid angle quoted above. The results mass
    vs. distance relation is fit with a $5^\mathrm{th}$-order polynomial. Figure~\ref{fig:BernardiSDSSDepthFit} shows the
    resulting relation between stellar mass and the maximum distance at which such a galaxy would be included in the
    sample. Points indicate results from Bernardi, while the line shows a polynomial fit:
    \begin{equation}
     \log_{10} \left[ {D_\mathrm{max}(M_\star) \over \hbox{Mpc}}\right] = 1282.11+m (-626.644+m (122.091+m (-11.8431+m
     (0.572399+m (-0.0110301)))))
     \label{eq:BernardiDepthPolynomial}
    \end{equation}
    where $m= \log_{10}(M_\star/M_\odot)$. We use this polynomial fit to determine the depth of the sample as a function of
    stellar mass.
    
    \begin{figure}
     \begin{center}
     \includegraphics[width=85mm,trim=0mm 0mm 0mm 4mm,clip]{Plots/DataAnalysis/BernardiSDSSMassLuminosityRelation.pdf}
     \end{center}
     \caption{The maximum distance at which a galaxy of given stellar mass can be detected in the sample of
     \protect\cite{bernardi_massive_2013}. Points show the results obtained from data provided by Bernardi, while the lines
     shows a polynomial fit to these results (given in eqn.~\ref{eq:BernardiDepthPolynomial}).}
     \label{fig:BernardiSDSSDepthFit}
    \end{figure}
   </description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryMangle) :: surveyGeometryBernardi2013SDSS
     private
   contains
     procedure :: fieldCount                => bernardi2013SDSSFieldCount
     procedure :: distanceMaximum           => bernardi2013SDSSDistanceMaximum
     procedure :: angularPowerMaximumDegree => bernardi2013SDSSAngularPowerMaximumDegree
     procedure :: mangleDirectory           => bernardi2013SDSSMangleDirectory
     procedure :: mangleFiles               => bernardi2013SDSSMangleFiles
     procedure :: pointIncluded             => bernardi2013SDSSPointIncluded
  end type surveyGeometryBernardi2013SDSS

  interface surveyGeometryBernardi2013SDSS
     !!{
     Constructors for the \cite{bernardi_massive_2013} survey geometry class.
     !!}
     module procedure bernardi2013SDSSConstructorParameters
     module procedure bernardi2013SDSSConstructorInternal
  end interface surveyGeometryBernardi2013SDSS

  ! Angular power spectra.
  integer, parameter :: angularPowerMaximumL=360

contains

  function bernardi2013SDSSConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \cite{bernardi_massive_2013} conditional mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(surveyGeometryBernardi2013SDSS)                :: self
    type(inputParameters               ), intent(inout) :: parameters

    self=surveyGeometryBernardi2013SDSS()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function bernardi2013SDSSConstructorParameters

  function bernardi2013SDSSConstructorInternal() result (self)
    !!{
    Default constructor for the \cite{bernardi_massive_2013} conditional mass function class.
    !!}
    implicit none
    type(surveyGeometryBernardi2013SDSS) :: self

    call self%initialize()
    return
  end function bernardi2013SDSSConstructorInternal

  integer function bernardi2013SDSSFieldCount(self)
    !!{
    Return the number of fields in this sample.
    !!}
    implicit none
    class(surveyGeometryBernardi2013SDSS), intent(inout) :: self
    !$GLC attributes unused :: self

    bernardi2013SDSSFieldCount=1
    return
  end function bernardi2013SDSSFieldCount

  double precision function bernardi2013SDSSDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    implicit none
    class           (surveyGeometryBernardi2013SDSS), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: mass           , magnitudeAbsolute, &
         &                                                                       luminosity     , starFormationRate
    integer                                         , intent(in   ), optional :: field
    double precision                                                          :: logarithmicMass
    !$GLC attributes unused :: self, field

    ! Validate arguments.
    if (present(magnitudeAbsolute)) call Error_Report('`magnitudeAbsolute` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! Find the limiting distance for this mass completeness limits. (See
    ! constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07_Bernardi/massDistanceRelation.pl for details.)
    if (present(mass)) then
       ! Limit the mass to the range for which our empirical mass-distance relation is calibrated.
       logarithmicMass=max(8.0d0,min(12.5d0,log10(mass)))
    else
       ! If no mass is given, use the largest calibrated mass (which will give the largest distance).
       logarithmicMass=              12.5d0
    end if
    bernardi2013SDSSDistanceMaximum                                                         &
         &                         =10.0d0**(                                               &
         &                                                         +1282.1065495948200000d0 &
         &                                   +logarithmicMass*(    - 626.6442739444630000d0 &
         &                                   +logarithmicMass* (   + 122.0914916099620000d0 &
         &                                   +logarithmicMass*  (  -  11.8431000301984000d0 &
         &                                   +logarithmicMass*   ( +   0.5723990783953920d0 &
         &                                   +logarithmicMass*    (-   0.0110301089727899d0 &
         &                                                        )                         &
         &                                                       )                          &
         &                                                      )                           &
         &                                                     )                            &
         &                                                    )                             &
         &                                  )
    return
  end function bernardi2013SDSSDistanceMaximum

  integer function bernardi2013SDSSAngularPowerMaximumDegree(self)
    !!{
    Return the maximum degree for which angular power is computed for the \cite{bernardi_massive_2013} survey.
    !!}
    implicit none
    class(surveyGeometryBernardi2013SDSS), intent(inout) :: self
    !$GLC attributes unused :: self

    bernardi2013SDSSAngularPowerMaximumDegree=angularPowerMaximumL
    return
  end function bernardi2013SDSSAngularPowerMaximumDegree

  function bernardi2013SDSSMangleDirectory(self)
    !!{
    Return the path to the directory containing \gls{mangle} files.
    !!}
    use :: Input_Paths, only : inputPath, pathTypeDataDynamic
    implicit none
    class(surveyGeometryBernardi2013SDSS), intent(inout) :: self
    type (varying_string                )                :: bernardi2013SDSSMangleDirectory
    !$GLC attributes unused :: self

    bernardi2013SDSSMangleDirectory=inputPath(pathTypeDataDynamic)//"surveyGeometry/SDSS/"
    return
  end function bernardi2013SDSSMangleDirectory

  subroutine bernardi2013SDSSMangleFiles(self,mangleFiles)
    !!{
    Return a list of \gls{mangle} files.
    !!}
    use :: File_Utilities , only : File_Exists      , File_Lock, File_Unlock, lockDescriptor, &
         &                         Directory_Make
    use :: Error          , only : Error_Report
    use :: System_Download, only : download
    use :: System_Command , only : System_Command_Do
    implicit none
    class  (surveyGeometryBernardi2013SDSS)                           , intent(inout) :: self
    type   (varying_string                ), allocatable, dimension(:), intent(inout) :: mangleFiles
    type   (lockDescriptor                )                                           :: lock
    integer                                                                           :: status

    allocate(mangleFiles(1))
    mangleFiles(1)=self%mangleDirectory()//"sdss_dr72safe0_res6d.pol"
    if (.not.File_Exists(mangleFiles(1))) then
       call Directory_Make(self%mangleDirectory())
       call File_Lock  (char(mangleFiles(1)),lock,lockIsShared=.false.)
       if (.not.File_Exists(mangleFiles(1))) then
          call download("https://zenodo.org/records/10998446/files/sdss_dr72safe0_res6d.pol.gz",char(mangleFiles(1))//".gz",status=status)
          if (status /= 0 .or. .not.File_Exists(mangleFiles(1)//".gz")) &
               & call Error_Report('failed to download mangle polygon file'//{introspection:location})
          call System_Command_Do("gunzip "//mangleFiles(1)//".gz",status)
          if (status /= 0 .or. .not.File_Exists(mangleFiles(1))) &
               & call Error_Report('failed to decompress mangle polygon file'//{introspection:location})
       end if
       call File_Unlock(                     lock                     )
    end if
    return
  end subroutine bernardi2013SDSSMangleFiles

  logical function bernardi2013SDSSPointIncluded(self,point,mass)
    !!{
    Return true if a point is included in the survey geometry.
    !!}
    use :: Numerical_Constants_Units, only : degree
    use :: Vectors                  , only : Vector_Magnitude
    implicit none
    class           (surveyGeometryBernardi2013SDSS), intent(inout)               :: self
    double precision                                , intent(in   ), dimension(3) :: point
    double precision                                , intent(in   )               :: mass
    double precision                                                              :: pointDistance, rightAscension, &
         &                                                                           declination

    ! Get the distance to the point.
    pointDistance=Vector_Magnitude(point)
    ! Compute if point lies within survey bounds.
    bernardi2013SDSSPointIncluded=                    &
         & pointDistance > self%distanceMinimum(mass) &
         & .and.                                      &
         & pointDistance < self%distanceMaximum(mass)
    if (.not.bernardi2013SDSSPointIncluded) return
    ! Exclude regions with no SDSS coverage.
    bernardi2013SDSSPointIncluded=.false.
    declination   =+90.0d0-acos (point(3)/pointDistance)/degree
    if     (                       &
         &   declination < -15.0d0 &
         &  .or.                   &
         &   declination > +75.0d0 &
         & ) return
    if     ( declination > +30.0d0) then
       rightAscension=    +atan2(point(2),point(1)     )/degree
       if     (                                                           &
            &   (rightAscension >   0.0d0 .and. rightAscension <  90.0d0) &
            &  .or.                                                       &
            &   (rightAscension > 270.0d0 .and. rightAscension < 360.0d0) &
            & ) return
    end if
    ! If point has not been excluded, perform the full test.
    bernardi2013SDSSPointIncluded=manglePointIncluded(self,point,mass)
    return
  end function bernardi2013SDSSPointIncluded
