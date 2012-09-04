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


!% Contains a module which implements the \cite{parkinson_generating_2008} modifier of merger tree branching rates.

module Merger_Tree_Branching_Modifiers_Parkinson
  !% Implements the \cite{parkinson_generating_2008} modifier of merger tree branching rates.
  implicit none
  private
  public :: Merger_Tree_Branching_Modifiers_Parkinson_Initialize

  ! Parameters of the algorithm.
  double precision :: modifiedPressSchechterG0,modifiedPressSchechterGamma1,modifiedPressSchechterGamma2

contains
  
  !# <treeBranchingModifierMethod>
  !#  <unitName>Merger_Tree_Branching_Modifiers_Parkinson_Initialize</unitName>
  !# </treeBranchingModifierMethod>
  subroutine Merger_Tree_Branching_Modifiers_Parkinson_Initialize(treeBranchingModifierMethod,Merger_Tree_Branching_Modifier_Get)
    !% Initialize the null modifier method for merger tree branching rates.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: treeBranchingModifierMethod
    procedure(double precision), pointer, intent(inout) :: Merger_Tree_Branching_Modifier_Get
    
    if (treeBranchingModifierMethod == 'Parkinson-Cole-Helly2008') then
       Merger_Tree_Branching_Modifier_Get => Merger_Tree_Branching_Modifier_Parkinson
       !@ <inputParameter>
       !@   <name>modifiedPressSchechterG0</name>
       !@   <defaultValue>0.57</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The parameter $G_0$ appearing in the modified merger rate expression of \cite{parkinson_generating_2008}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('modifiedPressSchechterG0'                ,modifiedPressSchechterG0                ,defaultValue=&
            & 0.57d0)
       !@ <inputParameter>
       !@   <name>modifiedPressSchechterGamma1</name>
       !@   <defaultValue>0.38</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The parameter $\gamma_1$ appearing in the modified merger rate expression of \cite{parkinson_generating_2008}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('modifiedPressSchechterGamma1'            ,modifiedPressSchechterGamma1            ,defaultValue=&
            & 0.38d0)
       !@ <inputParameter>
       !@   <name>modifiedPressSchechterGamma2</name>
       !@   <defaultValue>-0.01</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The parameter $\gamma_2$ appearing in the modified merger rate expression of \cite{parkinson_generating_2008}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('modifiedPressSchechterGamma2'            ,modifiedPressSchechterGamma2            ,defaultValue=&
            &-0.01d0)
    end if
    return
  end subroutine Merger_Tree_Branching_Modifiers_Parkinson_Initialize

  double precision function Merger_Tree_Branching_Modifier_Parkinson(parentDelta,childSigma,parentSigma)
    !% Returns a modifier for merger tree branching rates using the \cite{parkinson_generating_2008} algorithm.
    implicit none
    double precision, intent(in) :: parentDelta,childSigma,parentSigma
    double precision, save       :: parentDeltaPrevious=-1.0d0,parentSigmaPrevious=-1.0d0,parentTerm
    !$omp threadprivate(parentDeltaPrevious,parentSigmaPrevious,parentTerm)

    ! Check if we need to update the "parent" term.
    if (parentDelta /= parentDeltaPrevious .or. parentSigma /= parentSigmaPrevious) then
       ! "Parent" term must be updated. Compute and store it for future re-use.
       parentDeltaPrevious=parentDelta
       parentSigmaPrevious=parentSigma
       parentTerm= modifiedPressSchechterG0                                  &
            &     *((parentDelta/parentSigma)**modifiedPressSchechterGamma2) &
            &     /(             parentSigma **modifiedPressSchechterGamma1)
    end if

    ! Compute the modifier.
    Merger_Tree_Branching_Modifier_Parkinson=parentTerm*(childSigma**modifiedPressSchechterGamma1)
    return
  end function Merger_Tree_Branching_Modifier_Parkinson

end module Merger_Tree_Branching_Modifiers_Parkinson
