!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a critical overdensity for collapse the \gls{wdm} modifier of \cite{barkana_constraints_2001}.
  use Cosmology_Parameters
  
  !# <criticalOverdensity name="criticalOverdensityBarkana2001WDM" defaultThreaedPrivate="yes">
  !#  <description>Provides a critical overdensity for collapse based on the \gls{wdm} modifier of \cite{barkana_constraints_2001} applied to some other critical overdensity class.</description>
  !# </criticalOverdensity>
  type, extends(criticalOverdensityClass) :: criticalOverdensityBarkana2001WDM
     !% A critical overdensity for collapse class which modifies another transfer function using the \gls{wdm} modifier of \cite{barkana_constraints_2001}.
     private
     class           (criticalOverdensityClass), pointer                   :: criticalOverdensityCDM
     class           (cosmologyParametersClass), pointer                   :: cosmologyParameters_
     logical                                                               :: useFittingFunction
     double precision                                                      :: mx                             , gX            , &
          &                                                                   jeansMass
     integer                                                               :: deltaTableCount
     double precision                          , allocatable, dimension(:) :: deltaTableDelta                , deltaTableMass
     logical                                                               :: interpolationReset      =.true.
     type            (fgsl_interp_accel       )                            :: interpolationAccelerator
     type            (fgsl_interp             )                            :: interpolationObject
   contains
     final     ::               barkana2001WDMDestructor
     procedure :: value      => barkana2001WDMValue
     procedure :: descriptor => barkana2001WDMDescriptor
  end type criticalOverdensityBarkana2001WDM

  interface criticalOverdensityBarkana2001WDM
     !% Constructors for the ``{\normalfont \ttfamily barkana2001WDM}'' critical overdensity for collapse class.
     module procedure barkana2001WDMConstructorParameters
     module procedure barkana2001WDMConstructorInternal
  end interface criticalOverdensityBarkana2001WDM

  ! Parameters of fitting functions.
  double precision, parameter :: barkana2001WDMFitParameterA            =+2.40000d0
  double precision, parameter :: barkana2001WDMFitParameterB            =+0.10000d0
  double precision, parameter :: barkana2001WDMFitParameterC            =+0.04000d0
  double precision, parameter :: barkana2001WDMFitParameterD            =+2.30000d0
  double precision, parameter :: barkana2001WDMFitParameterE            =+0.31687d0
  double precision, parameter :: barkana2001WDMFitParameterF            =+0.80900d0
  double precision, parameter :: barkana2001WDMSmallMassLogarithmicSlope=-1.93400d0

contains

  function barkana2001WDMConstructorParameters(parameters)
    !% Constructor for the ``{\normalfont \ttfamily barkana2001WDM}'' critical overdensity for collapse class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type            (criticalOverdensityBarkana2001WDM)                :: barkana2001WDMConstructorParameters
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (criticalOverdensityClass         ), pointer       :: criticalOverdensityCDM
    class           (cosmologyParametersClass         ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass    ), pointer       :: cosmologicalMassVariance_
    double precision                                                   :: gX                                 , mX
    logical                                                            :: useFittingFunction
    !# <inputParameterList label="allowedParameterNames" />
    
    !# <inputParameter>
    !#   <name>gX</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.5d0</defaultValue>
    !#   <description>The effective number of degrees of freedom for the warm dark matter particle.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>mX</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The mass (in keV) of the warm dark matter particle.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>useFittingFunction</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>Specifies whether the warm dark matter critical overdensity mass scaling should be computed from a fitting function or from tabulated data.</description>
    !#   <type>boolean</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="criticalOverdensity"      name="criticalOverdensityCDM"    source="parameters"/>
    !# <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    ! Call the internal constructor
    barkana2001WDMConstructorParameters=barkana2001WDMConstructorInternal(criticalOverdensityCDM,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,gX,mX,useFittingFunction)
    return
  end function barkana2001WDMConstructorParameters

  function barkana2001WDMConstructorInternal(criticalOverdensityCDM,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,gX,mX,useFittingFunction)
    !% Internal constructor for the ``{\normalfont \ttfamily barkana2001WDM}'' critical overdensity for collapse class.
    use FoX_DOM
    use IO_XML
    use Galacticus_Input_Paths
    use Galacticus_Error
    use ISO_Varying_String
    implicit none
    type            (criticalOverdensityBarkana2001WDM)                        :: barkana2001WDMConstructorInternal
    class           (criticalOverdensityClass         ), target, intent(in   ) :: criticalOverdensityCDM
    class           (cosmologyParametersClass         ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), target, intent(in   ) :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass    ), target, intent(in   ) :: cosmologicalMassVariance_
    double precision                                           , intent(in   ) :: gX                               , mX
    logical                                                    , intent(in   ) :: useFittingFunction
    type            (node                             ), pointer               :: doc                              , element
    double precision                                                           :: matterRadiationEqualityRedshift
    integer                                                                    :: ioStatus

    barkana2001WDMConstructorInternal%criticalOverdensityCDM    => criticalOverdensityCDM
    barkana2001WDMConstructorInternal%cosmologyParameters_      => cosmologyParameters_
    barkana2001WDMConstructorInternal%cosmologyFunctions_       => cosmologyFunctions_
    barkana2001WDMConstructorInternal%cosmologicalMassVariance_ => cosmologicalMassVariance_
    barkana2001WDMConstructorInternal%useFittingFunction        =  useFittingFunction
    barkana2001WDMConstructorInternal%gX                        =  mX
    barkana2001WDMConstructorInternal%mX                        =  gX
    ! Compute corresponding Jeans mass.
    matterRadiationEqualityRedshift            =+3600.0d0                                                                                       &
         &                                      *(                                                                                              &
         &                                        +barkana2001WDMConstructorInternal%cosmologyParameters_%OmegaMatter   (                  )    &
         &                                        *barkana2001WDMConstructorInternal%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                                        /0.15d0                                                                                       &
         &                                      )                                                                                               &
         &                                      -1.0d0
    barkana2001WDMConstructorInternal%jeansMass=+3.06d8                                                                                                 &
         &                                      *((1.0d0+matterRadiationEqualityRedshift)/3000.0d0)                                             **1.5d0 &
         &                                      *sqrt(                                                                                                  &
         &                                            +barkana2001WDMConstructorInternal%cosmologyParameters_%OmegaMatter   (                  )        &
         &                                            *barkana2001WDMConstructorInternal%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2     &
         &                                            /0.15d0                                                                                           &
         &                                           )                                                                                                  &
         &                                      /(gX/1.5d0)                                                                                             &
         &                                      /(mX/1.0d0)                                                                                     **4    
    ! Read in the tabulated critical overdensity scaling.
    doc => parseFile(char(Galacticus_Input_Path())//"data/darkMatter/criticalOverdensityWarmDarkMatterBarkana.xml",iostat=ioStatus)
    if (ioStatus /= 0) call Galacticus_Error_Report('barkana2001WDMConstructorInternal','unable to find or parse the tabulated data')
    ! Extract the datum lists.
    element    => XML_Get_First_Element_By_Tag_Name(doc,"mass" )
    call XML_Array_Read(element,"datum",barkana2001WDMConstructorInternal%deltaTableMass )
    element    => XML_Get_First_Element_By_Tag_Name(doc,"delta")
    call XML_Array_Read(element,"datum",barkana2001WDMConstructorInternal%deltaTableDelta)
    barkana2001WDMConstructorInternal%deltaTableCount=size(barkana2001WDMConstructorInternal%deltaTableMass)
    ! Destroy the document.
    call destroy(doc)
    ! Convert tabulations to logarithmic versions.
    barkana2001WDMConstructorInternal%deltaTableMass =log(barkana2001WDMConstructorInternal%deltaTableMass )
    barkana2001WDMConstructorInternal%deltaTableDelta=log(barkana2001WDMConstructorInternal%deltaTableDelta)
    return
  end function barkana2001WDMConstructorInternal

  subroutine barkana2001WDMDestructor(self)
    !% Destructor for the barkana2001WDM critical overdensity for collapse class.
    implicit none
    type(criticalOverdensityBarkana2001WDM), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"  />
    !# <objectDestructor name="self%criticalOverdensityCDM"/>
    return
  end subroutine barkana2001WDMDestructor

  double precision function barkana2001WDMValue(self,time,expansionFactor,collapsing,mass)
    !% Returns a mass scaling for critical overdensities based on the results of \cite{barkana_constraints_2001}. This method
    !% assumes that their results for the original collapse barrier (i.e. the critical overdensity, and which they call $B_0$)
    !% scale with the effective Jeans mass of the warm dark matter particle as computed using their eqn.~(10).
    use Numerical_Interpolation
    use Table_Labels
    use FGSL
    implicit none
    class           (criticalOverdensityBarkana2001WDM), intent(inout)           :: self
    double precision                                   , intent(in   ), optional :: time                       , expansionFactor
    logical                                            , intent(in   ), optional :: collapsing
    double precision                                   , intent(in   ), optional :: mass
    double precision                                   , parameter               :: massScaleFreeMinimum=-10.d0
    double precision                                                             :: exponentialFit             , massScaleFree   , &
         &                                                                          powerLawFit                , smoothTransition

    ! Determine the scale-free mass.
    massScaleFree=log(mass/self%jeansMass)
    ! Compute the mass scaling via a fitting function or interpolation in tabulated results.
    if (self%useFittingFunction) then
       ! Impose a minimum to avoid divergence in fit.
       massScaleFree   =max(massScaleFree,massScaleFreeMinimum)
       smoothTransition=+1.0d0                                &
            &           /(                                    &
            &             +1.0d0                              &
            &             +exp(                               &
            &                  +(                             &
            &                    +massScaleFree               &
            &                    +barkana2001WDMFitParameterA &
            &                   )                             &
            &                  /  barkana2001WDMFitParameterB &
            &                 )                               &
            &            )
       powerLawFit     =+     barkana2001WDMFitParameterC & 
            &           /exp(                             &
            &                +barkana2001WDMFitParameterD &
            &                *massScaleFree               &
            &               )
       if (smoothTransition < 1.0d0) then
          exponentialFit=+exp(                                  &
               &              +     barkana2001WDMFitParameterE &
               &              /exp(                             &
               &                   +barkana2001WDMFitParameterF &
               &                   *massScaleFree               &
               &                  )                             &
               &             )
       else
          exponentialFit=+0.0d0
       end if
       barkana2001WDMValue=+       smoothTransition *powerLawFit    &
            &              +(1.0d0-smoothTransition)*exponentialFit
    else
       if      (massScaleFree > self%deltaTableMass(self%deltaTableCount)) then
          barkana2001WDMValue=1.0d0
       else if (massScaleFree < self%deltaTableMass(                   1)) then
          barkana2001WDMValue=exp(                                         &
               &                  +  self%deltaTableDelta(1)               &
               &                  +(                                       &
               &                    +massScaleFree                         &
               &                    -self%deltaTableMass (1)               &
               &                   )                                       &
               &                  *barkana2001WDMSmallMassLogarithmicSlope &
               &                 )
       else
          barkana2001WDMValue=exp(                                                                     &
               &                  Interpolate(                                                         &
               &                              self%deltaTableMass                                    , &
               &                              self%deltaTableDelta                                   , &
               &                              self%interpolationObject                               , &
               &                              self%interpolationAccelerator                          , &
               &                              massScaleFree                                          , &
               &                              extrapolationType            =extrapolationTypeFix     , &
               &                              reset                        =self%interpolationReset  , &
               &                              interpolationType            =FGSL_Interp_CSpline        &
               &                  )                                                                    &
               &                 )
       end if
    end if
    barkana2001WDMValue=barkana2001WDMValue*self%criticalOverdensityCDM%value(time,expansionFactor,collapsing,mass)
    return
  end function barkana2001WDMValue

  double precision function barkana2001WDMGradientTime(self,time,expansionFactor,collapsing,mass)
    !% Return the gradient with respect to time of critical overdensity at the given time and mass.
    implicit none
    class           (criticalOverdensityBarkana2001WDM), intent(inout)           :: self
    double precision                                   , intent(in   ), optional :: time      , expansionFactor
    logical                                            , intent(in   ), optional :: collapsing
    double precision                                   , intent(in   ), optional :: mass
    
    barkana2001WDMGradientTime=+self                       %value       (time,expansionFactor,collapsing,mass) &
         &                     *self%criticalOverdensityCDM%gradientTime(time,expansionFactor,collapsing,mass) &
         &                     /self%criticalOverdensityCDM%value       (time,expansionFactor,collapsing,mass)
    return
  end function barkana2001WDMGradientTime
  
  double precision function barkana2001WDMGradientMass(self,time,expansionFactor,collapsing,mass)
    !% Return the gradient with respect to mass of critical overdensity at the given time and mass.
    use Numerical_Interpolation
    use Table_Labels
    use FGSL
    implicit none
    class           (criticalOverdensityBarkana2001WDM), intent(inout)           :: self
    double precision                                   , intent(in   ), optional :: time                    , expansionFactor
    logical                                            , intent(in   ), optional :: collapsing
    double precision                                   , intent(in   ), optional :: mass
    double precision                                                             :: exponentialFit          , exponentialFitGradient, &
         &                                                                          massScaleFree           , powerLawFit           , &
         &                                                                          powerLawFitGradient     , smoothTransition      , &
         &                                                                          smoothTransitionGradient

    ! Determine the scale-free mass.
    massScaleFree=log(mass/self%jeansMass)
    ! Compute the mass scaling via a fitting function or interpolation in tabulated results.
    if (self%useFittingFunction) then
       smoothTransition          =+1.0d0                                   &
            &                     /(                                       &
            &                       +1.0d0                                 &
            &                       +exp(                                  &
            &                            +(                                &
            &                              +massScaleFree                  &
            &                              +barkana2001WDMFitParameterA    &
            &                             )                                &
            &                            /  barkana2001WDMFitParameterB    &
            &                           )                                  &
            &                      )
       powerLawFit               =+     barkana2001WDMFitParameterC        & 
            &                     /exp(                                    &
            &                          +barkana2001WDMFitParameterD        &
            &                          *massScaleFree                      &
            &                         )
       exponentialFit            =+exp(                                    &
            &                          +     barkana2001WDMFitParameterE   &
            &                          /exp(                               &
            &                               +barkana2001WDMFitParameterF   &
            &                               *massScaleFree                 &
            &                              )                               &
            &                         )
       powerLawFitGradient       =-barkana2001WDMFitParameterD             &
            &                     *powerLawFit                             &
            &                     /exp(                                    &
            &                          +massScaleFree                      &
            &                         )
       exponentialFitGradient    =-exponentialFit                          &
            &                     *barkana2001WDMFitParameterF             &
            &                     *barkana2001WDMFitParameterE             &
            &                     /exp(                                    &
            &                          +(                                  &
            &                            +1.0d0                            &
            &                            +barkana2001WDMFitParameterF      &
            &                           )                                  &
            &                          *massScaleFree                      &
            &                         )
       smoothTransitionGradient  =-exp(                                    &
            &                          +(                                  &
            &                            +     massScaleFree               &
            &                            +     barkana2001WDMFitParameterA &
            &                           )                                  &
            &                          /       barkana2001WDMFitParameterB &
            &                         )                                    &
            &                        /         barkana2001WDMFitParameterB &
            &                        /(                                    &
            &                          +1.0d0                              &
            &                          +exp(                               &
            &                               +(                             &
            &                                 +massScaleFree               &
            &                                 +barkana2001WDMFitParameterA &
            &                                )                             &
            &                               /  barkana2001WDMFitParameterB &
            &                              )                               &
            &                         )**2                                 &
            &                        /exp(                                 &
            &                             +massScaleFree                   &
            &                            )
       barkana2001WDMGradientMass=+(                                                                                         &
            &                       +       smoothTransition *powerLawFitGradient   +smoothTransitionGradient*powerLawFit    &
            &                       +(1.0d0-smoothTransition)*exponentialFitGradient-smoothTransitionGradient*exponentialFit &
            &                      )                                                                                         &
            &                     /self%jeansMass
    else
       if      (massScaleFree > self%deltaTableMass(self%deltaTableCount)) then
          barkana2001WDMGradientMass=+0.0d0
       else if (massScaleFree < self%deltaTableMass(                   1)) then
          barkana2001WDMGradientMass=+barkana2001WDMSmallMassLogarithmicSlope                                       &
               &                     *self%value(time,expansionFactor,collapsing,mass)                              &
               &                     /mass
       else
          barkana2001WDMGradientMass=+Interpolate_Derivative(                                                       &
               &                                             self%deltaTableMass                                  , &
               &                                             self%deltaTableDelta                                 , &
               &                                             self%interpolationObject                             , &
               &                                             self%interpolationAccelerator                        , &
               &                                             massScaleFree                                        , &
               &                                             extrapolationType            =extrapolationTypeFix   , &
               &                                             reset                        =self%interpolationReset, &
               &                                             interpolationType            =FGSL_Interp_CSpline      &
               &                                            )                                                       &
               &                     *self%value(time,expansionFactor,collapsing,mass)                              &
               &                     /mass
       end if
    end if
    ! Include gradient from CDM critical overdensity.
    barkana2001WDMGradientMass=+barkana2001WDMGradientMass                                                     &
         &                     +self                       %value       (time,expansionFactor,collapsing,mass) &
         &                     /self%criticalOverdensityCDM%value       (time,expansionFactor,collapsing,mass) &
         &                     *self%criticalOverdensityCDM%gradientMass(time,expansionFactor,collapsing,mass)
    return
  end function barkana2001WDMGradientMass

  subroutine barkana2001WDMDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class    (criticalOverdensityBarkana2001WDM), intent(inout) :: self
    type     (inputParameters                  ), intent(inout) :: descriptor
    type     (inputParameters                  )                :: subParameters
    character(len=10                           )                :: parameterLabel

    call descriptor%addParameter("criticalOverdensityMethod","barkana2001WDM")
    subParameters=descriptor%subparameters("criticalOverdensityMethod")
    write (parameterLabel,'(f10.6)') self%gX
    call subParameters%addParameter("gX"                ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(f10.6)') self%mX
    call subParameters%addParameter("mX"                ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(l1)   ') self%useFittingFunction
    call subParameters%addParameter("useFittingFunction",trim(adjustl(parameterLabel)))
    call self%criticalOverdensityCDM% descriptor(subParameters)
    return
  end subroutine barkana2001WDMDescriptor
