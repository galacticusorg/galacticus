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

!% Contains a module which implements a simple cooling rate calculation in which the cooling rate equals the mass of hot gas
!% divided by a timescale which is a function of halo mass and redshift.

module Cooling_Rates_Simple_Scaling
  !% Implements a simple cooling rate calculation in which the cooling rate equals the mass of hot gas
  !% divided by a timescale which is a function of halo mass and redshift.
  use Galacticus_Nodes
  implicit none
  private
  public :: Cooling_Rate_Simple_Scaling_Initialize

  ! The fixed timescale for cooling.
  double precision :: coolingRateSimpleScalingTimescale     , coolingRateSimpleScalingTimescaleExponent, & 
       &              coolingRateSimpleScalingTransitionMass                                               
  
contains

  !# <coolingRateMethod>
  !#  <unitName>Cooling_Rate_Simple_Scaling_Initialize</unitName>
  !# </coolingRateMethod>
  subroutine Cooling_Rate_Simple_Scaling_Initialize(coolingRateMethod,Cooling_Rate_Get)
    !% Initializes the ``simple scaling'' cooling rate module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type     (varying_string             ), intent(in   )          :: coolingRateMethod 
    procedure(Cooling_Rate_Simple_Scaling), intent(inout), pointer :: Cooling_Rate_Get  
    
    if (coolingRateMethod == 'simpleScaling') then
       Cooling_Rate_Get => Cooling_Rate_Simple_Scaling

       ! Get cooling rate parameters.
       !@ <inputParameter>
       !@   <name>coolingRateSimpleScalingTimescale</name>
       !@   <defaultValue>1 Gyr</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The timescale (in Gyr) for cooling in low mass halos at $z=0$ in the simple scaling cooling rate model.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingRateSimpleScalingTimescale',coolingRateSimpleScalingTimescale,defaultValue=1.0d0)
       !@ <inputParameter>
       !@   <name>coolingRateSimpleScalingTimescaleExponent</name>
       !@   <defaultValue>$-1.5$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of $(1+z)$ in the cooling timescale for low mass halos in the simple scaling cooling rate model.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingRateSimpleScalingTimescaleExponent',coolingRateSimpleScalingTimescaleExponent,defaultValue=-1.5d0)
       !@ <inputParameter>
       !@   <name>coolingRateSimpleScalingTransitionMass</name>
       !@   <defaultValue>$10^{12}M_\odot$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The halo mass scale appearing in the exponential term for cooling timescale in the simple scaling cooling rate model.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingRateSimpleScalingTransitionMass',coolingRateSimpleScalingTransitionMass,defaultValue=1.0d12)

       ! Check that the properties we need are gettable.
       if (.not.defaultHotHaloComponent%massIsGettable())                             &
            & call Galacticus_Error_Report(                                           &
            &                              'Cooling_Rate_Simple_Scaling_Initialize' , &
            &                              'hotHalo component mass must be gettable'  &
            &                             )
       if (.not.defaultBasicComponent%massIsGettable  ())                             &
            & call Galacticus_Error_Report(                                           &
            &                              'Cooling_Rate_Simple_Scaling_Initialize' , &
            &                              'basic component mass must be gettable'    &
            &                             )
       if (.not.defaultBasicComponent%timeIsGettable  ())                             &
            & call Galacticus_Error_Report(                                           &
            &                              'Cooling_Rate_Simple_Scaling_Initialize' , &
            &                              'basic component time must be gettable'    &
            &                             )
    end if
    return
  end subroutine Cooling_Rate_Simple_Scaling_Initialize

  double precision function Cooling_Rate_Simple_Scaling(thisNode)
    !% Computes the mass cooling rate in a hot gas halo assuming a fixed timescale for cooling.
    use Cosmology_Functions
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode                                      
    double precision                      , parameter              :: massRatioMaximum    =100.0d0                  
    class           (nodeComponentBasic  )               , pointer :: thisBasicComponent                            
    class           (nodeComponentHotHalo)               , pointer :: thisHotHaloComponent                          
    double precision                                               :: coolingRate                 , expansionFactor 
    
    thisBasicComponent   => thisNode%basic  ()
    thisHotHaloComponent => thisNode%hotHalo()
    expansionFactor=Expansion_Factor(thisBasicComponent%time())
    coolingRate    = exp(-thisBasicComponent%mass()/coolingRateSimpleScalingTransitionMass) &
         &          *expansionFactor**coolingRateSimpleScalingTimescaleExponent             &
         &          /coolingRateSimpleScalingTimescale
    Cooling_Rate_Simple_Scaling=thisHotHaloComponent%mass()*coolingRate
    return
  end function Cooling_Rate_Simple_Scaling

end module Cooling_Rates_Simple_Scaling
