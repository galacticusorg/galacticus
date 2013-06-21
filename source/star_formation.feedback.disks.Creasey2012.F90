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

!% Contains a module which implements the \cite{creasey_how_2012} model for star formation feedback in galactic disks.

module Star_Formation_Feedback_Disks_Creasey2012
  !% Implements the \cite{creasey_how_2012} model for star formation feedback in galactic disks.
  use Galacticus_Nodes
  implicit none
  private
  public :: Star_Formation_Feedback_Disks_Creasey2012_Initialize

  ! Parameters of the feedback model.
  double precision                    :: starFormationFeedbackDisksCreasy2012Beta0, starFormationFeedbackDisksCreasy2012Mu, & 
       &                                 starFormationFeedbackDisksCreasy2012Nu                                               
  
  ! Pointer to active node used in integral functions, plus variables needed by integral function.
  type            (treeNode), pointer :: activeNode                                                                           
  !$omp threadprivate(activeNode)
contains

  !# <starFormationFeedbackDisksMethod>
  !#  <unitName>Star_Formation_Feedback_Disks_Creasey2012_Initialize</unitName>
  !# </starFormationFeedbackDisksMethod>
  subroutine Star_Formation_Feedback_Disks_Creasey2012_Initialize(starFormationFeedbackDisksMethod,Star_Formation_Feedback_Disk_Outflow_Rate_Get)
    !% Initializes the ``Creasey et al. (2012)'' disk star formation feedback module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                                       ), intent(in   )          :: starFormationFeedbackDisksMethod              
    procedure(Star_Formation_Feedback_Disk_Outflow_Rate_Creasey2012), intent(inout), pointer :: Star_Formation_Feedback_Disk_Outflow_Rate_Get 
    
    if (starFormationFeedbackDisksMethod == 'Creasey2012') then
       Star_Formation_Feedback_Disk_Outflow_Rate_Get => Star_Formation_Feedback_Disk_Outflow_Rate_Creasey2012
       ! Get parameters of for the feedback calculation.
       !@ <inputParameter>
       !@   <name>starFormationFeedbackDisksCreasy2012Beta0</name>
       !@   <defaultValue>13</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The factor $\beta_0$ appearing in the \cite{creasey_how_2012} model for supernovae feedback.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationFeedbackDisksCreasy2012Beta0',starFormationFeedbackDisksCreasy2012Beta0,defaultValue=13.0d0)
       !@ <inputParameter>
       !@   <name>starFormationFeedbackDisksCreasy2012Mu</name>
       !@   <defaultValue>$1.15$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The factor $\mu$ appearing in the \cite{creasey_how_2012} model for supernovae feedback.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationFeedbackDisksCreasy2012Mu',starFormationFeedbackDisksCreasy2012Mu,defaultValue=1.15d0)
       !@ <inputParameter>
       !@   <name>starFormationFeedbackDisksCreasy2012Nu</name>
       !@   <defaultValue>$0.16$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The factor $\nu$ appearing in the \cite{creasey_how_2012} model for supernovae feedback.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationFeedbackDisksCreasy2012Nu',starFormationFeedbackDisksCreasy2012Nu,defaultValue=0.16d0)
    end if
    return
  end subroutine Star_Formation_Feedback_Disks_Creasey2012_Initialize

  double precision function Star_Formation_Feedback_Disk_Outflow_Rate_Creasey2012(thisNode,starFormationRate,energyInputRate)
    !% Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic disk of {\tt thisNode} using
    !% the model of \cite{creasey_how_2012}. The outflow rate is given by
    !% \begin{equation}
    !% \dot{M}_{\rm outflow} = \int_0^\infty \beta_0 \Sigma_{g,1}^{-\mu}(r) f_{\rm g}^\nu(r) \dot{\Sigma}_\star(r) 2 \pi r {\rm d}r,
    !% \end{equation}
    !% where $\Sigma_{g,1}(r)$ is the surface density of gas in units of $M_\odot$ pc$^{-2}$, $f_{\rm g}(r)$ is the gas fraction,
    !% $\dot{\Sigma}_\star(r)$ is the surface density of star formation rate, $\beta_0=${\tt
    !% [starFormationFeedbackDisksCreasy2012Beta0]}, $\mu=${\tt [starFormationFeedbackDisksCreasy2012Mu]}, and $\nu=${\tt
    !% [starFormationFeedbackDisksCreasy2012Nu]}.
    use               Galacticus_Nodes
    use               Numerical_Constants_Math
    use               Stellar_Feedback
    use               FGSL
    use               Numerical_Integration
    use, intrinsic :: ISO_C_Binding
    implicit none
    type            (treeNode                  ), intent(inout), pointer :: thisNode                                                           
    class           (nodeComponentDisk         )               , pointer :: thisDiskComponent                                                  
    double precision                            , intent(in   )          :: energyInputRate               , starFormationRate                  
    double precision                            , parameter              :: radiusInnerDimensionless=0.0d0, radiusOuterDimensionless=10.0d0    
    double precision                                                     :: diskScaleRadius               , gasMass                        , & 
         &                                                                  radiusInner                   , radiusOuter                    , & 
         &                                                                  stellarMass                                                        
    type            (c_ptr                     )                         :: parameterPointer                                                   
    type            (fgsl_function             )                         :: integrandFunction                                                  
    type            (fgsl_integration_workspace)                         :: integrationWorkspace                                               
    
    ! Get the disk properties.
    thisDiskComponent => thisNode         %disk       ()
    gasMass           =  thisDiskComponent%massGas    ()
    stellarMass       =  thisDiskComponent%massStellar()
    diskScaleRadius   =  thisDiskComponent%radius     ()
    ! Return immediately for a null disk.
    if (gasMass <= 0.0d0 .or. stellarMass <= 0.0d0 .or. diskScaleRadius <= 0.0d0) then
       Star_Formation_Feedback_Disk_Outflow_Rate_Creasey2012=0.0d0
       return
    end if
    ! Set a pointer to the node that is accessible by integral function.
    activeNode => thisNode
    ! Compute suitable limits for the integration.
    radiusInner=diskScaleRadius*radiusInnerDimensionless
    radiusOuter=diskScaleRadius*radiusOuterDimensionless
    ! Compute the outflow rate.
    Star_Formation_Feedback_Disk_Outflow_Rate_Creasey2012=                             &
         &  2.0d0*Pi                                                                   &
         & *starFormationFeedbackDisksCreasy2012Beta0                                  &
         & *Integrate(                                                                 &
         &            radiusInner                                                    , &
         &            radiusOuter                                                    , &
         &            Star_Formation_Feedback_Disk_Outflow_Rate_Creasey2012_Integrand, &
         &            parameterPointer                                               , &
         &            integrandFunction                                              , &
         &            integrationWorkspace                                           , &
         &            toleranceAbsolute=0.0d0                                        , &
         &            toleranceRelative=1.0d-3                                         &
         &           )                                                                 &
         & /starFormationRate                                                          &
         & *energyInputRate                                                            &
         & /feedbackEnergyInputAtInfinityCanonical
       call Integrate_Done(integrandFunction,integrationWorkspace)
    return
  end function Star_Formation_Feedback_Disk_Outflow_Rate_Creasey2012

  function Star_Formation_Feedback_Disk_Outflow_Rate_Creasey2012_Integrand(radius,parameterPointer) bind(c)
    !% Integrand function for the ``Creasey et al. (2012)'' supernovae feedback calculation.
    use               Galactic_Structure_Surface_Densities
    use               Galactic_Structure_Options
    use               Star_Formation_Rate_Surface_Density_Disks
    use               Numerical_Constants_Prefixes
    use, intrinsic :: ISO_C_Binding
    implicit none
    real            (kind=c_double)        :: Star_Formation_Feedback_Disk_Outflow_Rate_Creasey2012_Integrand                                     
    real            (kind=c_double), value :: radius                                                                                              
    type            (c_ptr        ), value :: parameterPointer                                                                                    
    double precision                       :: gasFraction                                                    , starFormationRateSurfaceDensity, & 
         &                                    surfaceDensityGas                                              , surfaceDensityStar                 
    
    ! Get gas surface density.
    surfaceDensityGas=Galactic_Structure_Surface_Density(activeNode,[radius,0.0d0,0.0d0],coordinateSystem&
         &=coordinateSystemCylindrical,componentType=componentTypeDisk,massType=massTypeGaseous)
    ! Get stellar surface density.
    surfaceDensityStar=Galactic_Structure_Surface_Density(activeNode,[radius,0.0d0,0.0d0],coordinateSystem&
         &=coordinateSystemCylindrical,componentType=componentTypeDisk,massType=massTypeStellar)
    ! Compute the gas fraction.
    gasFraction=surfaceDensityGas/(surfaceDensityGas+surfaceDensityStar)
    ! Convert gas surface density to correct units.
    surfaceDensityGas=surfaceDensityGas/mega**2
    ! Get the surface density of star formation rate.
    starFormationRateSurfaceDensity=Star_Formation_Rate_Surface_Density_Disk(activeNode,radius)
    ! Compute the outflow rate.
    Star_Formation_Feedback_Disk_Outflow_Rate_Creasey2012_Integrand=     &
         &  surfaceDensityGas**(-starFormationFeedbackDisksCreasy2012Mu) &
         & *gasFraction      **  starFormationFeedbackDisksCreasy2012Nu  &
         & *starFormationRateSurfaceDensity                              &
         & *radius
    return
  end function Star_Formation_Feedback_Disk_Outflow_Rate_Creasey2012_Integrand
  
end module Star_Formation_Feedback_Disks_Creasey2012
