!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Implements the gravitational lensing distribution by modifying another distribution for the effects of baryons.

  !# <gravitationalLensing name="gravitationalLensingBaryonicModifier">
  !#  <description>Implements the gravitational lensing distribution by modifying another distribution for the effects of baryons.</description>
  !# </gravitationalLensing>

  type, extends(gravitationalLensingClass) :: gravitationalLensingBaryonicModifier
     class           (gravitationalLensingClass), pointer :: originalDistribution
     double precision                                     :: alpha                  , beta               , &
          &                                                  transitionMagnification, renormalization    , &
          &                                                  redshiftPrevious       , scaleSourcePrevious
   contains
     !@ <objectMethods>
     !@   <object>gravitationalLensingBaryonicModifier</object>
     !@   <objectMethod>
     !@     <method>renormalize</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ redshift\argin, \doublezero\ scaleSource\argin</arguments>
     !@     <description>Renormalize the gravitational lensing magnification distribution function.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                     baryonicModifierDestructor
     procedure :: magnificationPDF => baryonicModifierMagnificationPDF
     procedure :: magnificationCDF => baryonicModifierMagnificationCDF
     procedure :: renormalize      => baryonicModifierRenormalize
  end type gravitationalLensingBaryonicModifier
  
  interface gravitationalLensingBaryonicModifier
     !% Constructors for the ``baryonic modifier'' gravitational lensing class.
     module procedure baryonicModifierDefaultConstructor
  end interface gravitationalLensingBaryonicModifier

  ! Initialization state.
  logical                          :: baryonicModifierInitialized                             =.false.

  ! Default parameters.
  type            (varying_string) :: gravitationalLensingBaryonicModifierOriginalDistribution
  double precision                 :: gravitationalLensingBaryonicModifierAlpha
  double precision                 :: gravitationalLensingBaryonicModifierBeta

contains

  function baryonicModifierDefaultConstructor()
    !% Default constructor for the ``baryonic modifier'' gravitational lensing class.
    use Input_Parameters
    implicit none
    type(gravitationalLensingBaryonicModifier) :: baryonicModifierDefaultConstructor

    if (.not.baryonicModifierInitialized) then
       !$omp critical (baryonicModifierInitialize)
       if (.not.baryonicModifierInitialized) then
          !@ <inputParameter>
          !@   <name>gravitationalLensingBaryonicModifierOriginalDistribution</name>
          !@   <defaultValue>takahashi2011</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Name of the original gravitational lensing magnification distribution method to which to apply baryonic modifiers.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('gravitationalLensingBaryonicModifierOriginalDistribution',gravitationalLensingBaryonicModifierOriginalDistribution,defaultValue='takahashi2011')
          !@ <inputParameter>
          !@   <name>gravitationalLensingBaryonicModifierAlpha</name>
          !@   <defaultValue>$0$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Parameter $\alpha$ in the modified gravitational lensing \gls{pdf}, $P(\mu) \rightarrow P(\mu) + \hbox{min}[\alpha,\beta P(\mu)]$.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('gravitationalLensingBaryonicModifierAlpha',gravitationalLensingBaryonicModifierAlpha,defaultValue=0.0d0)
          !@ <inputParameter>
          !@   <name>gravitationalLensingBaryonicModifierBeta</name>
          !@   <defaultValue>$0$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Parameter $\beta in the modified gravitational lensing \gls{pdf}, $P(\mu) \rightarrow P(\mu) + \hbox{min}[\alpha,\beta P(\mu)]$.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('gravitationalLensingBaryonicModifierBeta',gravitationalLensingBaryonicModifierBeta,defaultValue=0.0d0)
          ! Record that we are now initialized.
          baryonicModifierInitialized=.true.
       end if
       !$omp end critical (baryonicModifierInitialize)
    end if
    baryonicModifierDefaultConstructor                                                            &
         & =baryonicModifierConstructor(                                                          &
         &                              gravitationalLensingBaryonicModifierOriginalDistribution, &
         &                              gravitationalLensingBaryonicModifierAlpha               , &
         &                              gravitationalLensingBaryonicModifierBeta                  &
         &                             )
    return
  end function baryonicModifierDefaultConstructor

  function baryonicModifierConstructor(originalDistributionName,alpha,beta)
    !% Generic constructor for the ``baryonic modifier'' gravitational lensing class.
    implicit none
    type            (gravitationalLensingBaryonicModifier)                :: baryonicModifierConstructor
    type            (varying_string                      ), intent(in   ) :: originalDistributionName
    double precision                                      , intent(in   ) :: alpha                            , beta

    ! Set parameters of this object.
    baryonicModifierConstructor%originalDistribution => gravitationalLensing(char(originalDistributionName))
    baryonicModifierConstructor%alpha                =  alpha
    baryonicModifierConstructor%beta                 =  beta
    baryonicModifierConstructor%redshiftPrevious     =  -1.0d0
    baryonicModifierConstructor%scaleSourcePrevious  =  -1.0d0
    return
  end function baryonicModifierConstructor

  subroutine baryonicModifierDestructor(self)
    !% Destructor for the ``baryonic modifier'' gravitational lensing class.
    implicit none
    type(gravitationalLensingBaryonicModifier), intent(inout) :: self

    if (self%originalDistribution%isFinalizable()) deallocate(self%originalDistribution)
    return
  end subroutine baryonicModifierDestructor

  subroutine baryonicModifierRenormalize(self,redshift,scaleSource)
    !% Renormlize for \gls{pdf} for baryonic modification.
    use Root_Finder
    implicit none
    class           (gravitationalLensingBaryonicModifier), intent(inout) :: self
    double precision                                      , intent(in   ) :: redshift               , scaleSource
    double precision                                      , parameter     :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-6
    type            (rootFinder                          ), save          :: finder
    !$omp threadprivate(finder)
    
    ! Exit if nothing has changed since the previous call.
    if (redshift == self%redshiftPrevious .and. scaleSource == self%scaleSourcePrevious) return
    ! Trap case of no modification.
    if (self%beta == 0.0d0) then
       self%transitionMagnification=1.0d0
       self%renormalization        =1.0d0
       return
    end if
    ! Check for delta-function magnification PDF.
    if (self%originalDistribution%magnificationCDF(1.0d0,redshift,scaleSource) >= 1.0d0) then
       self%transitionMagnification=1.0d0
       self%renormalization        =1.0d0
       return
    end if
    ! Find the magnification at which we transition from additive to multiplicative correction.
    if (.not.finder%isInitialized()) then
       call finder%rootFunction(magnificationTransition            )
       call finder%tolerance   (toleranceAbsolute,toleranceRelative)
       call finder%rangeExpand (                                                             &
            &                   rangeExpandUpward            =2.0d0                        , &
            &                   rangeExpandDownward          =0.5d0                        , &
            &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
            &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &                   rangeExpandType              =rangeExpandMultiplicative      &
            &                  )
    end if
    self%transitionMagnification=finder%find(rootGuess=2.0d0)
    ! Find renormalization.
    self%renormalization=+1.0d0                                                                                             &
         &               /(                                                                                                 &
         &                 +1.0d0                                                                                           &
         &                 +self%alpha                                                                                      &
         &                 *(                                                                                               &
         &                   +                                           self%transitionMagnification                       &
         &                   -1.0d0                                                                                         &
         &                  )                                                                                               &
         &                 +self%beta                                                                                       &
         &                 *(                                                                                               &
         &                   +1.0d0                                                                                         &
         &                   -self%originalDistribution%magnificationCDF(self%transitionMagnification,redshift,scaleSource) &
         &                  )                                                                                               &
         &                )
    ! Record arguments.
    self%redshiftPrevious   =redshift
    self%scaleSourcePrevious=scaleSource

  contains
    
    double precision function magnificationTransition(magnification)
      !% Root finding function used in the ``baryonic modifier'' gravitational lensing class.
      implicit none
      double precision, intent(in   ) :: magnification

      magnificationTransition=self%originalDistribution%magnificationPDF(magnification,redshift,scaleSource)-self%alpha/self%beta
      return
    end function magnificationTransition

  end subroutine baryonicModifierRenormalize

  double precision function baryonicModifierMagnificationPDF(self,magnification,redshift,scaleSource)
    !% Compute the magnification probability density function at the given {\tt magnification} and {\tt redshift} by modifying
    !% another distribution for the effects of baryons.
    implicit none
    class           (gravitationalLensingBaryonicModifier), intent(inout) :: self
    double precision                                      , intent(in   ) :: magnification, redshift, &
         &                                                                   scaleSource

    call self%renormalize(redshift,scaleSource)
    baryonicModifierMagnificationPDF=self%originalDistribution%magnificationPDF(magnification,redshift,scaleSource)
    if      (magnification > self%transitionMagnification) then
       baryonicModifierMagnificationPDF=baryonicModifierMagnificationPDF*(1.0d0+self%beta )
    else if (magnification > 1.0d0                       ) then
       baryonicModifierMagnificationPDF=baryonicModifierMagnificationPDF       +self%alpha
    end if
    baryonicModifierMagnificationPDF=baryonicModifierMagnificationPDF*self%renormalization
    return
  end function baryonicModifierMagnificationPDF
  
  double precision function baryonicModifierMagnificationCDF(self,magnification,redshift,scaleSource)
    !% Compute the magnification probability density function at the given {\tt magnification} and {\tt redshift} by modifying
    !% another distribution for the effects of baryons. 
    implicit none
    class           (gravitationalLensingBaryonicModifier), intent(inout) :: self
    double precision                                      , intent(in   ) :: magnification, redshift, &
         &                                                                   scaleSource

    call self%renormalize(redshift,scaleSource)
    baryonicModifierMagnificationCDF=self%originalDistribution%magnificationCDF(magnification,redshift,scaleSource)
    if      (magnification > self%transitionMagnification) then
       baryonicModifierMagnificationCDF=+baryonicModifierMagnificationCDF                                                                &
            &                           +self%alpha                                                                                      &
            &                           *(self%transitionMagnification-1.0d0)                                                            &
            &                           +self%beta                                                                                       &
            &                           *(                                                                                               &
            &                             +baryonicModifierMagnificationCDF                                                              &
            &                             -self%originalDistribution%magnificationCDF(self%transitionMagnification,redshift,scaleSource) &
            &                           )
    else if (magnification > 1.0d0                       ) then
       baryonicModifierMagnificationCDF=+baryonicModifierMagnificationCDF     &
            &                           +self%alpha                           &
            &                           *(               magnification-1.0d0)
    end if
    baryonicModifierMagnificationCDF=baryonicModifierMagnificationCDF*self%renormalization
    return    
  end function baryonicModifierMagnificationCDF
  
