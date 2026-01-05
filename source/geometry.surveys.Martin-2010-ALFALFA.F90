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
Implements the survey geometry used by \cite{martin_arecibo_2010}.
!!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <surveyGeometry name="surveyGeometryMartin2010ALFALFA">
   <description>
    A survey geometry class that describes the survey geometry of \cite{martin_arecibo_2010}. 
    
    For the angular mask we use the three disjoint regions defined by 07$^\mathrm{h}$30$^\mathrm{m}$ $&lt;$ R.A. $&lt;$
    16$^\mathrm{h}$30$^\mathrm{m}$, +04$^\circ$ $&lt;$ decl. $&lt;$ +16$^\circ$, and +24$^\circ$ $&lt;$ decl. $&lt;$ +28$^\circ$ and
    22$^\mathrm{h}$ $&lt;$ R.A. $&lt;$ 03$^\mathrm{h}$, +14$^\circ$ $&lt;$ decl. $&lt;$ +16$^\circ$, and +24$^\circ$ $&lt;$ decl. $&lt;$
    +32$^\circ$ corresponding to the sample of \cite{martin_arecibo_2010}. When the survey window function is needed we
    generate randomly distributed points within this angular mask and out to the survey depth. These points are used to
    determine which elements of a 3D grid fall within the window function.
    
    To estimate the depth of the \cite{martin_arecibo_2010} sample as a function of galaxy HI mass we first infer the median
    line width corresponding to that mass. To do so, we have fit the median line width-mass relation from the $\alpha.40$
    sample with power-law function as shown in Fig.~\ref{fig:ALFALFALineWidthMassRelation}. We find that the median line width
    can be approximated by
    \begin{equation}
     \log_{10} (W_\mathrm{50}/\hbox{km s}^{-1}) = c_0 + c_1 \log_10(M_\mathrm{HI}/M_\odot),
     \label{eq:ALFALFALineWidthMassRelation}
    \end{equation}
    with $c_0=-0.770$ and $c_1=0.315$. Given the line width, the corresponding integrated flux limit, $S_\mathrm{int}$, for a
    signal-to-noise of $6.5$ is inferred using equation~(A1) of \cite{haynes_arecibo_2011}. Finally, this integrated flux limit
    is converted to maximum distance at which the source could be detected using the expression given in the text of
    section~2.2 of \cite{martin_arecibo_2010}:
    \begin{equation}
     M_\mathrm{HI} = 2.356\times10^5 \left({D\over \hbox{Mpc}}\right)^2 \left({S_\mathrm{int}\over\hbox{Jy km s}^{-1}}\right).
    \end{equation}
    
    \begin{figure}
     \begin{center}
      \includegraphics[width=85mm,trim=0mm 0mm 0mm 4mm,clip]{Plots/DataAnalysis/alfalfaHILineWidthMassRelation.pdf}
     \end{center}
     \caption{HI line width vs. HI mass as measured from the $\alpha.40$ survey of \protect\cite{martin_arecibo_2010}. Red
     points with error bars show individual measurements, while the larger circles indicate the running median of these
     data. The green line is a power-law fit to the running median as described in
     eqn.~(\protect\ref{eq:ALFALFALineWidthMassRelation}).}
     \label{fig:ALFALFALineWidthMassRelation}
    \end{figure}
   </description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryRandomPoints) :: surveyGeometryMartin2010ALFALFA
     private
     class(cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
   contains
     final     ::                      martin2010ALFALFADestructor
     procedure :: distanceMaximum   => martin2010ALFALFADistanceMaximum
     procedure :: solidAngle        => martin2010ALFALFASolidAngle
     procedure :: randomsInitialize => martin2010ALFALFARandomsInitialize
  end type surveyGeometryMartin2010ALFALFA

  interface surveyGeometryMartin2010ALFALFA
     !!{
     Constructors for the \cite{martin_arecibo_2010} survey geometry class.
     !!}
     module procedure martin2010ALFALFAConstructorParameters
     module procedure martin2010ALFALFAConstructorInternal
  end interface surveyGeometryMartin2010ALFALFA

  ! Minimum and maximum recession velocities for a galaxy to be admitted to the sample.
  double precision, parameter :: sampleVelocityMinimum=0.0d0, sampleVelocityMaximum=15.0d3

contains

  function martin2010ALFALFAConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{martin_arecibo_2010} conditional mass function class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (surveyGeometryMartin2010ALFALFA)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(cosmologyParametersClass       ), pointer       :: cosmologyParameters_
    class(randomNumberGeneratorClass     ), pointer       :: randomNumberGenerator_

    !![
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=surveyGeometryMartin2010ALFALFA(cosmologyParameters_,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function martin2010ALFALFAConstructorParameters

  function martin2010ALFALFAConstructorInternal(cosmologyParameters_,randomNumberGenerator_) result(self)
    !!{
    Internal constructor for the \cite{martin_arecibo_2010} conditional mass function class.
    !!}
    implicit none
    type (surveyGeometryMartin2010ALFALFA)                                  :: self
    class(cosmologyParametersClass       ), intent(in   ), target           :: cosmologyParameters_
    class(randomNumberGeneratorClass     ), intent(in   ), target, optional :: randomNumberGenerator_
    !![
    <constructorAssign variables="*cosmologyParameters_ ,*randomNumberGenerator_"/>
    !!]

    self%geometryInitialized=.false.
    return
  end function martin2010ALFALFAConstructorInternal

  subroutine martin2010ALFALFADestructor(self)
    !!{
    Destructor for the \refClass{surveyGeometryMartin2010ALFALFA} survey geometry class.
    !!}
    implicit none
    type(surveyGeometryMartin2010ALFALFA), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine martin2010ALFALFADestructor

  double precision function martin2010ALFALFADistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (surveyGeometryMartin2010ALFALFA), intent(inout)           :: self
    double precision                                 , intent(in   ), optional :: mass                                                  , magnitudeAbsolute                        , &
         &                                                                        luminosity                                            , starFormationRate
    integer                                          , intent(in   ), optional :: field
    ! The signal-to-noise limit used by Martin et al. (2010).
    double precision                                 , parameter               :: signalToNoise                    = 6.5d0
    ! Coefficients of the polynomial approximation for log10(lineWidth) vs. log10(HI mass).
    double precision                                 , parameter               :: lineWidthCoefficient0            =-0.769635671616885d0, lineWidthCoefficient1=0.314983275066432d0
    ! Line width characteristic scale.
    double precision                                 , parameter               :: lineWidthCharacteristic          =200.0d0
    ! Normalization of the flux limit for unit signal-to-noise at characteristic line width.
    double precision                                 , parameter               :: integratedFluxLimitNormalization=0.15d0
    ! Normalization of the mass-integrated flux-distance relation.
    double precision                                 , parameter               :: massNormalization               =2.356d5
    double precision                                                           :: logarithmicMass                                       , lineWidth                                , &
         &                                                                        integratedFluxLimit

    ! Validate field.
    if (present(field).and.field /= 1) call Error_Report('field = 1 required'//{introspection:location})
    ! Validate arguments.
    if (present(magnitudeAbsolute)) call Error_Report('`magnitudeAbsolute` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! Find the maximum distance.
    if (present(mass)) then
       ! Get the logarithm of the mass.
       logarithmicMass=log10(mass)
       ! Find the median line width for this mass. (See
       ! constraints/dataAnalysis/hiMassFunction_ALFALFA_z0.00/lineWidthMassRelation.pl for details.)
       lineWidth=10.0d0**(lineWidthCoefficient0+lineWidthCoefficient1*logarithmicMass)
       ! Compute the limiting integrated flux using equation (A1) of Martin et al. (2010).
       if (lineWidth < lineWidthCharacteristic) then
          integratedFluxLimit=integratedFluxLimitNormalization*signalToNoise*sqrt(lineWidth/lineWidthCharacteristic)
       else
          integratedFluxLimit=integratedFluxLimitNormalization*signalToNoise*    (lineWidth/lineWidthCharacteristic)
       end if
       ! Convert from mass and limiting integrated flux to maximum distance using relation given in text of section 2.2 of Martin et
       ! al. (2010). Limit by the maximum velocity allowed for galaxies to make it into the sample.
       martin2010ALFALFADistanceMaximum=min(sqrt(mass/massNormalization/integratedFluxLimit),sampleVelocityMaximum/self%cosmologyParameters_%hubbleConstant())
    else
       martin2010ALFALFADistanceMaximum=                                                     sampleVelocityMaximum/self%cosmologyParameters_%hubbleConstant()
    end if
    return
  end function martin2010ALFALFADistanceMaximum

  double precision function martin2010ALFALFASolidAngle(self,field)
    !!{
    Return the solid angle of the \cite{martin_arecibo_2010} sample.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (surveyGeometryMartin2010ALFALFA), intent(inout)           :: self
    integer                                          , intent(in   ), optional :: field
    double precision                                 , parameter               :: solidAngleSurvey=0.79415674617213461d0 ! Computed from survey bounds in Martin et al. (2010; ApJ; 723; 1359)
    !$GLC attributes unused :: self

    ! Validate field.
    if (present(field).and.field /= 1) call Error_Report('field = 1 required'//{introspection:location})
    martin2010ALFALFASolidAngle=solidAngleSurvey
    return
  end function martin2010ALFALFASolidAngle

  subroutine martin2010ALFALFARandomsInitialize(self)
    !!{
    Initialize random points for the survey.
    !!}
    use :: Numerical_Constants_Astronomical, only : degreesToRadians       , hoursToRadians
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (surveyGeometryMartin2010ALFALFA), intent(inout)                           :: self
    integer                                          , parameter                               :: randomsCount=1000000
    integer                                          , parameter                               :: regionCount =      3
    ! Survey geometry from Haynes et al. (2011; http://adsabs.harvard.edu/abs/2011AJ....142..170H).
    double precision                                 , parameter    , dimension(2,regionCount) :: regionRightAscensionRange=reshape([22.0d0,03.0d0,07.5d0,16.5d0,07.5d0,16.5d0],[2,regionCount]),                  &
         &                                                                                        regionDeclinationRange   =reshape([24.0d0,32.0d0,04.0d0,16.0d0,24.0d0,28.0d0],[2,regionCount])
    double precision                                                , dimension(2,regionCount) :: regionPhiRange                                                                                , regionThetaRange
    double precision                                                , dimension(  regionCount) :: regionSolidAngle
    integer                                                                                    :: iRegion                                                                                       , iRandom
    double precision                                                                           :: uniformRandom

    ! Determine the solid angles of the different survey regions.
    regionPhiRange  =        +regionRightAscensionRange*  hoursToRadians
    regionThetaRange=0.5d0*Pi-   regionDeclinationRange*degreesToRadians
    do iRegion=1,regionCount
       if (regionPhiRange(2,iRegion) < regionPhiRange(1,iRegion)) regionPhiRange(2,iRegion)=regionPhiRange(2,iRegion)+2.0d0*Pi
    end do
    regionSolidAngle=(regionPhiRange(2,:)-regionPhiRange(1,:))*(cos(regionThetaRange(2,:))-cos(regionThetaRange(1,:)))
    ! Cumulate and normalize the region solid angles.
    do iRegion=2,regionCount
       regionSolidAngle(iRegion)=regionSolidAngle(iRegion)+regionSolidAngle(iRegion-1)
    end do
    regionSolidAngle=regionSolidAngle/regionSolidAngle(regionCount)
    ! Generate random points.
    allocate(self%randomTheta(randomsCount))
    allocate(self%randomPhi  (randomsCount))
    do iRandom=1,randomsCount
       ! Select a region at random.
       uniformRandom=self%randomNumberGenerator_%uniformSample()
       iRegion=1
       do while (uniformRandom > regionSolidAngle(iRegion))
          iRegion=iRegion+1
       end do
       ! Select coordinates at random within this region.
       uniformRandom            =self%randomNumberGenerator_%uniformSample()
       self%randomPhi  (iRandom)=     uniformRandom*(    regionPhiRange  (2,iRegion) -    regionPhiRange  (1,iRegion)) +    regionPhiRange  (1,iRegion)
       uniformRandom            =self%randomNumberGenerator_%uniformSample()
       self%randomTheta(iRandom)=acos(uniformRandom*(cos(regionThetaRange(2,iRegion))-cos(regionThetaRange(1,iRegion)))+cos(regionThetaRange(1,iRegion)))
    end do
    return
  end subroutine martin2010ALFALFARandomsInitialize
