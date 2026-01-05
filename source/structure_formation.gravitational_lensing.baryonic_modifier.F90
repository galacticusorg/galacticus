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
Implements the gravitational lensing distribution by modifying another distribution for the effects of baryons.
!!}

  !![
  <gravitationalLensing name="gravitationalLensingBaryonicModifier">
   <description>
    A gravitational lensing distribution class which (approximately) modifies another distribution for the effects of
    baryons. The distribution to modify is specified via the {\normalfont \ttfamily
    [gravitationalLensingBaryonicModifierOriginalDistribution]} parameter. The modification takes the form:
    \begin{equation}
    P(\mu) \rightarrow P(\mu) + \hbox{min}[\alpha,\beta P(\mu)]
    \end{equation}
    where $\alpha=${\normalfont \ttfamily [gravitationalLensingBaryonicModifierAlpha]} and $\beta=${\normalfont \ttfamily
    [gravitationalLensingBaryonicModifierBeta]}. The distribution is then renormalized to ensure that the cumulative
    probability reaches unity for infinite magnification. As an example, values of $\alpha=2.05\times 10^{-3}$ and $\beta=0.62$
    approximately reproduce the results of \cite[][their Fig.~1]{hilbert_strong-lensing_2008}.
   </description>
  </gravitationalLensing>
  !!]
  type, extends(gravitationalLensingClass) :: gravitationalLensingBaryonicModifier
     class           (gravitationalLensingClass), pointer :: gravitationalLensing_   => null()
     double precision                                     :: alpha                            , beta               , &
          &                                                  transitionMagnification          , renormalization    , &
          &                                                  redshiftPrevious                 , scaleSourcePrevious
   contains
     !![
     <methods>
       <method description="Renormalize the gravitational lensing magnification distribution function." method="renormalize" />
     </methods>
     !!]
     final     ::                     baryonicModifierDestructor
     procedure :: magnificationPDF => baryonicModifierMagnificationPDF
     procedure :: magnificationCDF => baryonicModifierMagnificationCDF
     procedure :: renormalize      => baryonicModifierRenormalize
  end type gravitationalLensingBaryonicModifier

  interface gravitationalLensingBaryonicModifier
     !!{
     Constructors for the \refClass{gravitationalLensingBaryonicModifier} gravitational lensing class.
     !!}
     module procedure baryonicModifierConstructorParameters
     module procedure baryonicModifierConstructorInternal
  end interface gravitationalLensingBaryonicModifier

contains

  function baryonicModifierConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily baryonicModifier} gravitational lensing class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (gravitationalLensingBaryonicModifier)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (gravitationalLensingClass           ), pointer       :: gravitationalLensing_
    double precision                                                      :: alpha                , beta

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>alpha</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>Parameter $\alpha$ in the modified gravitational lensing \gls{pdf}, $P(\mu) \rightarrow P(\mu) + \hbox{min}[\alpha,\beta P(\mu)]$.</description>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>Parameter $\beta$ in the modified gravitational lensing \gls{pdf}, $P(\mu) \rightarrow P(\mu) + \hbox{min}[\alpha,\beta P(\mu)]$.</description>
    </inputParameter>
    <objectBuilder class="gravitationalLensing" name="gravitationalLensing_" source="parameters"/>
    !!]
    ! Build the object.
    self=gravitationalLensingBaryonicModifier(gravitationalLensing_,alpha,beta)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="gravitationalLensing_"/>
    !!]
    return
  end function baryonicModifierConstructorParameters

  function baryonicModifierConstructorInternal(gravitationalLensing_,alpha,beta) result(self)
    !!{
    Internal constructor for the \refClass{gravitationalLensingBaryonicModifier} gravitational lensing class.
    !!}
    implicit none
    type            (gravitationalLensingBaryonicModifier)                        :: self
    class           (gravitationalLensingClass           ), intent(in   ), target :: gravitationalLensing_
    double precision                                      , intent(in   )         :: alpha                , beta
    !![
    <constructorAssign variables="*gravitationalLensing_, alpha, beta"/>
    !!]

    ! Initialize
    self%redshiftPrevious   =-1.0d0
    self%scaleSourcePrevious=-1.0d0
    return
  end function baryonicModifierConstructorInternal

  subroutine baryonicModifierDestructor(self)
    !!{
    Destructor for the \refClass{gravitationalLensingBaryonicModifier} gravitational lensing class.
    !!}
    implicit none
    type(gravitationalLensingBaryonicModifier), intent(inout) :: self

    !![
    <objectDestructor name="self%gravitationalLensing_"/>
    !!]
    return
  end subroutine baryonicModifierDestructor

  subroutine baryonicModifierRenormalize(self,redshift,scaleSource)
    !!{
    Renormalize for \gls{pdf} for baryonic modification.
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (gravitationalLensingBaryonicModifier), intent(inout) :: self
    double precision                                      , intent(in   ) :: redshift                 , scaleSource
    double precision                                      , parameter     :: toleranceAbsolute =0.0d0 , toleranceRelative=1.0d-6
    type            (rootFinder                          ), save          :: finder
    logical                                               , save          :: finderConstructed=.false.
    !$omp threadprivate(finder,finderConstructed)

    ! Exit if nothing has changed since the previous call.
    if (redshift == self%redshiftPrevious .and. scaleSource == self%scaleSourcePrevious) return
    ! Trap case of no modification.
    if (self%beta == 0.0d0) then
       self%transitionMagnification=1.0d0
       self%renormalization        =1.0d0
       return
    end if
    ! Check for δ-function magnification PDF.
    if (self%gravitationalLensing_%magnificationCDF(1.0d0,redshift,scaleSource) >= 1.0d0) then
       self%transitionMagnification=1.0d0
       self%renormalization        =1.0d0
       return
    end if
    ! Find the magnification at which we transition from additive to multiplicative correction.
    if (.not.finderConstructed) then
       finder=rootFinder(                                                             &
            &            rootFunction                 =magnificationTransition      , &
            &            toleranceAbsolute            =toleranceAbsolute            , &
            &            toleranceRelative            =toleranceRelative            , &
            &            rangeExpandUpward            =2.0d0                        , &
            &            rangeExpandDownward          =0.5d0                        , &
            &            rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
            &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &            rangeExpandType              =rangeExpandMultiplicative      &
            &           )
       finderConstructed=.true.
    end if
    self%transitionMagnification=finder%find(rootGuess=2.0d0)
    ! Find renormalization.
    self%renormalization=+1.0d0                                                                                              &
         &               /(                                                                                                  &
         &                 +1.0d0                                                                                            &
         &                 +self%alpha                                                                                       &
         &                 *(                                                                                                &
         &                   +                                           self%transitionMagnification                        &
         &                   -1.0d0                                                                                          &
         &                  )                                                                                                &
         &                 +self%beta                                                                                        &
         &                 *(                                                                                                &
         &                   +1.0d0                                                                                          &
         &                   -self%gravitationalLensing_%magnificationCDF(self%transitionMagnification,redshift,scaleSource) &
         &                  )                                                                                                &
         &                )
    ! Record arguments.
    self%redshiftPrevious   =redshift
    self%scaleSourcePrevious=scaleSource

  contains

    double precision function magnificationTransition(magnification)
      !!{
      Root finding function used in the \refClass{gravitationalLensingBaryonicModifier} gravitational lensing class.
      !!}
      implicit none
      double precision, intent(in   ) :: magnification

      magnificationTransition=self%gravitationalLensing_%magnificationPDF(magnification,redshift,scaleSource)-self%alpha/self%beta
      return
    end function magnificationTransition

  end subroutine baryonicModifierRenormalize

  double precision function baryonicModifierMagnificationPDF(self,magnification,redshift,scaleSource)
    !!{
    Compute the magnification probability density function at the given {\normalfont \ttfamily magnification} and {\normalfont \ttfamily redshift} by modifying
    another distribution for the effects of baryons.
    !!}
    implicit none
    class           (gravitationalLensingBaryonicModifier), intent(inout) :: self
    double precision                                      , intent(in   ) :: magnification, redshift, &
         &                                                                   scaleSource

    call self%renormalize(redshift,scaleSource)
    baryonicModifierMagnificationPDF=self%gravitationalLensing_%magnificationPDF(magnification,redshift,scaleSource)
    if      (magnification > self%transitionMagnification) then
       baryonicModifierMagnificationPDF=baryonicModifierMagnificationPDF*(1.0d0+self%beta )
    else if (magnification > 1.0d0                       ) then
       baryonicModifierMagnificationPDF=baryonicModifierMagnificationPDF       +self%alpha
    end if
    baryonicModifierMagnificationPDF=baryonicModifierMagnificationPDF*self%renormalization
    return
  end function baryonicModifierMagnificationPDF

  double precision function baryonicModifierMagnificationCDF(self,magnification,redshift,scaleSource)
    !!{
    Compute the magnification probability density function at the given {\normalfont \ttfamily magnification} and {\normalfont \ttfamily redshift} by modifying
    another distribution for the effects of baryons.
    !!}
    implicit none
    class           (gravitationalLensingBaryonicModifier), intent(inout) :: self
    double precision                                      , intent(in   ) :: magnification, redshift, &
         &                                                                   scaleSource

    call self%renormalize(redshift,scaleSource)
    baryonicModifierMagnificationCDF=self%gravitationalLensing_%magnificationCDF(magnification,redshift,scaleSource)
    if      (magnification > self%transitionMagnification) then
       baryonicModifierMagnificationCDF=+baryonicModifierMagnificationCDF                                                                 &
            &                           +self%alpha                                                                                       &
            &                           *(self%transitionMagnification-1.0d0)                                                             &
            &                           +self%beta                                                                                        &
            &                           *(                                                                                                &
            &                             +baryonicModifierMagnificationCDF                                                               &
            &                             -self%gravitationalLensing_%magnificationCDF(self%transitionMagnification,redshift,scaleSource) &
            &                           )
    else if (magnification > 1.0d0                       ) then
       baryonicModifierMagnificationCDF=+baryonicModifierMagnificationCDF     &
            &                           +self%alpha                           &
            &                           *(               magnification-1.0d0)
    end if
    baryonicModifierMagnificationCDF=baryonicModifierMagnificationCDF*self%renormalization
    return
  end function baryonicModifierMagnificationCDF

