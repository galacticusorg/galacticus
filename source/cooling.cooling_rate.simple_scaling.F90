!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements a simple cooling rate calculation in which the cooling rate equals the mass of hot gas
!% divided by a timescale which is a function of halo mass and redshift.

module Cooling_Rates_Simple_Scaling
  !% Implements a simple cooling rate calculation in which the cooling rate equals the mass of hot gas
  !% divided by a timescale which is a function of halo mass and redshift.
  use Tree_Nodes
  implicit none
  private
  public :: Cooling_Rate_Simple_Scaling_Initialize

  ! The fixed timescale for cooling.
  double precision :: coolingRateSimpleScalingTimescale,coolingRateSimpleScalingTimescaleExponent,coolingRateSimpleScalingTransitionMass

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
    type(varying_string),                 intent(in   ) :: coolingRateMethod
    procedure(double precision), pointer, intent(inout) :: Cooling_Rate_Get
    
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
       if (.not.associated(Tree_Node_Hot_Halo_Mass        ))                                   &
            & call Galacticus_Error_Report(                                                    &
            &                              'Cooling_Rate_Simple_Scaling_Initialize'          , &
            &                              'Tree_Node_Hot_Halo_Mass must be gettable'          &
            &                             )
       if (.not.associated(Tree_Node_Mass                 ))                                   &
            & call Galacticus_Error_Report(                                                    &
            &                              'Cooling_Rate_Simple_Scaling_Initialize'          , &
            &                              'Tree_Node_Mass must be gettable'                   &
            &                             )
       if (.not.associated(Tree_Node_Time                 ))                                   &
            & call Galacticus_Error_Report(                                                    &
            &                              'Cooling_Rate_Simple_Scaling_Initialize'          , &
            &                              'Tree_Node_Time must be gettable'                   &
            &                             )
    end if
    return
  end subroutine Cooling_Rate_Simple_Scaling_Initialize

  double precision function Cooling_Rate_Simple_Scaling(thisNode)
    !% Computes the mass cooling rate in a hot gas halo assuming a fixed timescale for cooling.
    use Tree_Nodes
    use Cosmology_Functions
    implicit none
    type(treeNode)  , intent(inout), pointer :: thisNode
    double precision, parameter              :: massRatioMaximum=100.0d0
    double precision                         :: coolingRate,expansionFactor
    
    expansionFactor=Expansion_Factor(Tree_Node_Time(thisNode))
    coolingRate    = exp(-Tree_Node_Mass(thisNode)/coolingRateSimpleScalingTransitionMass) &
         &          *expansionFactor**coolingRateSimpleScalingTimescaleExponent            &
         &          /coolingRateSimpleScalingTimescale
    Cooling_Rate_Simple_Scaling=Tree_Node_Hot_Halo_Mass(thisNode)*coolingRate
    return
  end function Cooling_Rate_Simple_Scaling

end module Cooling_Rates_Simple_Scaling
