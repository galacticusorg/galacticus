!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a warm dark matter scaling of critical overdensities for collapse based on the work of \cite{barkana_constraints_2001}.

module Critical_Overdensity_Mass_Scalings_WDM
  !% Implements a warm dark matter scaling of critical overdensities for collapse based on the work of \cite{barkana_constraints_2001}.
  use Galacticus_Nodes
  use FGSL
  implicit none
  private
  public :: Critical_Overdensity_Mass_Scaling_WDM_Initialize

  ! The Jeans mass used in the warm dark matter fitting formula.
  double precision                                               :: jeansMass

  ! Tabulation of the critical overdensity scaling.
  integer                                                        :: deltaTableCount
  double precision                   , allocatable, dimension(:) :: deltaTableDelta                                              , deltaTableMass

  ! Interpolation objects.
  logical                                                        :: interpolationReset                                 =.true.
  type            (fgsl_interp_accel)                            :: interpolationAccelerator
  type            (fgsl_interp      )                            :: interpolationObject
  double precision                   , parameter                 :: smallMassLogarithmicSlope                          =-1.934d0

  ! Option controlling use of tabulated data or fitting function.
  logical                                                        :: warmDarkMatterCriticalOverdensityUseFittingFunction

  ! Parameters of the fitting function.
  double precision                   , parameter                 :: fitParameterA                                      =2.40000d0
  double precision                   , parameter                 :: fitParameterB                                      =0.10000d0
  double precision                   , parameter                 :: fitParameterC                                      =0.04000d0
  double precision                   , parameter                 :: fitParameterD                                      =2.30000d0
  double precision                   , parameter                 :: fitParameterE                                      =0.31687d0
  double precision                   , parameter                 :: fitParameterF                                      =0.80900d0

contains

  !# <criticalOverdensityMassScalingMethod>
  !#  <unitName>Critical_Overdensity_Mass_Scaling_WDM_Initialize</unitName>
  !# </criticalOverdensityMassScalingMethod>
  subroutine Critical_Overdensity_Mass_Scaling_WDM_Initialize(criticalOverdensityMassScalingMethod&
       &,Critical_Overdensity_Mass_Scaling_Get,Critical_Overdensity_Mass_Scaling_Gradient_Get)
    !% Initializes the ``warmDarkMatter'' critical overdensity mass scaling method.
    use ISO_Varying_String
    use Input_Parameters
    use Cosmological_Parameters
    use FoX_dom
    use Galacticus_Input_Paths
    use Galacticus_Error
    use Memory_Management
    implicit none
    type            (varying_string  ), intent(in   )          :: criticalOverdensityMassScalingMethod
    procedure       (double precision), intent(inout), pointer :: Critical_Overdensity_Mass_Scaling_Get, Critical_Overdensity_Mass_Scaling_Gradient_Get
    type            (Node            )               , pointer :: doc                                  , thisNode
    type            (NodeList        )               , pointer :: deltaDatumList                       , massDatumList                                 , &
         &                                                        thisList
    integer                                                    :: iDatum                               , ioErr
    double precision                                           :: matterRadiationEqualityRedshift      , warmDarkMatterCriticalOverdensityGX           , &
         &                                                        warmDarkMatterCriticalOverdensityMX

    if (criticalOverdensityMassScalingMethod == 'warmDarkMatter') then
       ! Return a pointer to our implementation of the mass scaling function.
       Critical_Overdensity_Mass_Scaling_Get          => Critical_Overdensity_Mass_Scaling_WDM
       Critical_Overdensity_Mass_Scaling_Gradient_Get => Critical_Overdensity_Mass_Scaling_Gradient_WDM

       ! Get warm dark matter particle properties.
       !@ <inputParameter>
       !@   <name>warmDarkMatterCriticalOverdensityGX</name>
       !@   <defaultValue>1.5</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The effective number of degrees of freedom for the warm dark matter particle.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('warmDarkMatterCriticalOverdensityGX',warmDarkMatterCriticalOverdensityGX,defaultValue=1.5d0)
       !@ <inputParameter>
       !@   <name>warmDarkMatterCriticalOverdensityMX</name>
       !@   <defaultValue>1.0 keV</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass (in keV) of the warm dark matter particle.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('warmDarkMatterCriticalOverdensityMX',warmDarkMatterCriticalOverdensityMX,defaultValue=1.0d0)
       !@ <inputParameter>
       !@   <name>warmDarkMatterCriticalOverdensityUseFittingFunction</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether the warm dark matter critical overdensity mass scaling should be computed from a fitting function or from tabulated data.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('warmDarkMatterCriticalOverdensityUseFittingFunction',warmDarkMatterCriticalOverdensityUseFittingFunction,defaultValue=.true.)

       ! Compute corresponding Jeans mass.
       matterRadiationEqualityRedshift=3600.0d0*(Omega_Matter()*Little_H_0()**2/0.15d0)-1.0d0
       jeansMass=3.06d8*((1.0d0+matterRadiationEqualityRedshift)/3000.0d0)**1.5d0*sqrt(Omega_Matter()*Little_H_0()**2/0.15d0)&
            &/(warmDarkMatterCriticalOverdensityGX/1.5d0)/(warmDarkMatterCriticalOverdensityMX/1.0d0)**4

       ! Read in the tabulated critical overdensity scaling.
       doc => parseFile(char(Galacticus_Input_Path())//"data/darkMatter/criticalOverdensityWarmDarkMatterBarkana.xml",iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Critical_Overdensity_Mass_Scaling_WDM_Initialize','unable to find or parse the tabulated data')
       ! Extract the datum lists.
       thisList       => getElementsByTagname(doc     ,"mass" )
       thisNode       => item(thisList,0)
       massDatumList  => getElementsByTagname(thisNode,"datum")
       thisList       => getElementsByTagname(doc     ,"delta")
       thisNode       => item(thisList,0)
       deltaDatumList => getElementsByTagname(thisNode,"datum")
       deltaTableCount=getLength(massDatumList)
       call Alloc_Array(deltaTableMass ,[deltaTableCount])
       call Alloc_Array(deltaTableDelta,[deltaTableCount])
       do iDatum=0,deltaTableCount-1
          thisNode => item(massDatumList ,iDatum)
          call extractDataContent(thisNode,deltaTableMass (iDatum+1))
          thisNode => item(deltaDatumList,iDatum)
          call extractDataContent(thisNode,deltaTableDelta(iDatum+1))
       end do
       ! Destroy the document.
       call destroy(doc)
       ! Convert tabulations to logarithmic versions.
       deltaTableMass =log(deltaTableMass )
       deltaTableDelta=log(deltaTableDelta)
    end if
    return
  end subroutine Critical_Overdensity_Mass_Scaling_WDM_Initialize

  double precision function Critical_Overdensity_Mass_Scaling_WDM(mass)
    !% Returns a mass scaling for critical overdensities based on the results of \cite{barkana_constraints_2001}. This method
    !% assumes that their results for the original collapse barrier (i.e. the critical overdensity, and which they call $B_0$)
    !% scale with the effective Jeans mass of the warm dark matter particle as computed using their eqn.~(10).
    use Numerical_Interpolation
    implicit none
    double precision, intent(in   ) :: mass
    double precision, parameter     :: massScaleFreeMinimum=-10.d0
    double precision                :: exponentialFit             , massScaleFree   , &
         &                             powerLawFit                , smoothTransition

    ! Determine the scale-free mass.
    massScaleFree=log(mass/jeansMass)

    ! Compute the mass scaling via a fitting function or interpolation in tabulated results.
    if (warmDarkMatterCriticalOverdensityUseFittingFunction) then
       massScaleFree   =max(massScaleFree,massScaleFreeMinimum) ! Impose a minimum to avoid divergence in fit.
       smoothTransition=1.0d0/(1.0d0+exp((massScaleFree+fitParameterA)/fitParameterB))
       powerLawFit     =fitParameterC/exp(fitParameterD*massScaleFree)
       if (smoothTransition < 1.0d0) then
          exponentialFit=exp(fitParameterE/exp(fitParameterF*massScaleFree))
       else
          exponentialFit=0.0d0
       end if
       Critical_Overdensity_Mass_Scaling_WDM=smoothTransition*powerLawFit+(1.0d0-smoothTransition)*exponentialFit
    else
       if (massScaleFree > deltaTableMass(deltaTableCount)) then
          Critical_Overdensity_Mass_Scaling_WDM=1.0d0
       else if (massScaleFree < deltaTableMass(1)) then
          Critical_Overdensity_Mass_Scaling_WDM=exp(deltaTableDelta(1)+(massScaleFree-deltaTableMass(1))*smallMassLogarithmicSlope)
       else
          Critical_Overdensity_Mass_Scaling_WDM=exp(Interpolate(deltaTableCount,deltaTableMass,deltaTableDelta,interpolationObject&
               &,interpolationAccelerator,massScaleFree,extrapolationType=extrapolationTypeFixed,reset=interpolationReset,interpolationType=fgsl_interp_cspline))
       end if
    end if
    return
  end function Critical_Overdensity_Mass_Scaling_WDM

  double precision function Critical_Overdensity_Mass_Scaling_Gradient_WDM(mass)
    !% Returns a mass scaling for critical overdensities based on the results of \cite{barkana_constraints_2001}. This method
    !% assumes that their results for the original collapse barrier (i.e. the critical overdensity, and which they call $B_0$)
    !% scale with the effective Jeans mass of the warm dark matter particle as computed using their eqn.~(10).
    use Numerical_Interpolation
    implicit none
    double precision, intent(in   ) :: mass
    double precision                :: exponentialFit          , exponentialFitGradient, &
         &                             massScaleFree           , powerLawFit           , &
         &                             powerLawFitGradient     , smoothTransition      , &
         &                             smoothTransitionGradient

    ! Determine the scale-free mass.
    massScaleFree=log(mass/jeansMass)

    ! Compute the mass scaling via a fitting function or interpolation in tabulated results.
    if (warmDarkMatterCriticalOverdensityUseFittingFunction) then
       smoothTransition        =1.0d0/(1.0d0+exp((massScaleFree+fitParameterA)/fitParameterB))
       powerLawFit             =fitParameterC/exp(fitParameterD*massScaleFree)
       exponentialFit          =exp(fitParameterE/exp(fitParameterF*massScaleFree))
       powerLawFitGradient     =-fitParameterD*powerLawFit/exp(massScaleFree)
       exponentialFitGradient  =-exponentialFit*fitParameterF*fitParameterE/exp((1.0d0+fitParameterF)*massScaleFree)
       smoothTransitionGradient=-exp((massScaleFree+fitParameterA)/fitParameterB)/fitParameterB/(1.0d0+exp((massScaleFree&
            &+fitParameterA)/fitParameterB))**2/exp(massScaleFree)
       Critical_Overdensity_Mass_Scaling_Gradient_WDM=(smoothTransition*powerLawFitGradient+smoothTransitionGradient*powerLawFit&
            &+(1.0d0-smoothTransition)*exponentialFitGradient-smoothTransitionGradient*exponentialFit)/jeansMass
    else
       if (massScaleFree > deltaTableMass(deltaTableCount)) then
          Critical_Overdensity_Mass_Scaling_Gradient_WDM=0.0d0
       else if (massScaleFree < deltaTableMass(1)) then
          Critical_Overdensity_Mass_Scaling_Gradient_WDM=smallMassLogarithmicSlope*Critical_Overdensity_Mass_Scaling_WDM(mass)/mass
       else
          Critical_Overdensity_Mass_Scaling_Gradient_WDM=Interpolate_Derivative(deltaTableCount,deltaTableMass,deltaTableDelta&
               &,interpolationObject,interpolationAccelerator,massScaleFree,extrapolationType=extrapolationTypeFixed,reset&
               &=interpolationReset,interpolationType=fgsl_interp_cspline)*Critical_Overdensity_Mass_Scaling_WDM(mass)/mass
       end if
    end if
    return
  end function Critical_Overdensity_Mass_Scaling_Gradient_WDM

end module Critical_Overdensity_Mass_Scalings_WDM
