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


!% Contains a module which implements calculations of merger tree branching probabilities.

module Merger_Tree_Branching
  !% Implements calculations of merger tree branching probabilities.
  use ISO_Varying_String
  implicit none
  private
  public :: Tree_Branching_Probability, Tree_Subresolution_Fraction, Tree_Branch_Mass, Tree_Maximum_Step

  ! Flag to indicate if this module has been initialized.  
  logical                                        :: treeBranchingInitialized=.false.

  ! Name of branching method used.
  type(varying_string)                           :: treeBranchingMethod

  ! Pointer to the functions that return branching probabilities and templates for the functions.
  procedure(Tree_Branching_Probability_Template ), pointer :: Tree_Branching_Probability_Function  => null()
  procedure(Tree_Subresolution_Fraction_Template), pointer :: Tree_Subresolution_Fraction_Function => null()
  procedure(Tree_Branch_Mass_Template           ), pointer :: Tree_Branch_Mass_Function            => null()
  procedure(Tree_Maximum_Step_Template          ), pointer :: Tree_Maximum_Step_Function           => null()
  interface Tree_Branching_Probability_Template
     double precision function Tree_Branching_Probability_Template(haloMass,deltaCritical,massResolution)
       double precision, intent(in) :: haloMass,deltaCritical,massResolution
     end function Tree_Branching_Probability_Template
  end interface
  interface Tree_Subresolution_Fraction_Template
     double precision function Tree_Subresolution_Fraction_Template(haloMass,deltaCritical,massResolution)
       double precision, intent(in) :: haloMass,deltaCritical,massResolution
     end function Tree_Subresolution_Fraction_Template
  end interface
  interface Tree_Branch_Mass_Template
     double precision function Tree_Branch_Mass_Template(haloMass,deltaCritical,massResolution,probabilityFraction)
       double precision, intent(in) :: haloMass,deltaCritical,massResolution,probabilityFraction
     end function Tree_Branch_Mass_Template
  end interface
  interface Tree_Maximum_Step_Template
     double precision function Tree_Maximum_Step_Template(haloMass,deltaCritical,massResolution)
       double precision, intent(in) :: haloMass,deltaCritical,massResolution
     end function Tree_Maximum_Step_Template
  end interface
 
contains

  double precision function Tree_Maximum_Step(haloMass,deltaCritical,massResolution)
    !% Return the maximum step in $\delta_{\rm crit}$ allowed for a halo in a merger tree.
    implicit none
    double precision, intent(in) :: haloMass,deltaCritical,massResolution

    ! Initialize if necessary.
    call Tree_Branching_Initialize

    ! Interpolate in the tabulated function and return a value.
    Tree_Maximum_Step=Tree_Maximum_Step_Function(haloMass,deltaCritical,massResolution)
    return
  end function Tree_Maximum_Step

  double precision function Tree_Branching_Probability(haloMass,deltaCritical,massResolution)
    !% Return the branching probability per unit $\delta_{\rm crit}$ for a halo in a merger tree.
    implicit none
    double precision, intent(in) :: haloMass,deltaCritical,massResolution

    ! Initialize if necessary.
    call Tree_Branching_Initialize

    ! Interpolate in the tabulated function and return a value.
    Tree_Branching_Probability=Tree_Branching_Probability_Function(haloMass,deltaCritical,massResolution)
    return
  end function Tree_Branching_Probability

  double precision function Tree_Subresolution_Fraction(haloMass,deltaCritical,massResolution)
    !% Return the fraction of mass accreted below the resolution limit per $\delta_{\rm crit}$ in a halo in a merger tree.
    implicit none
    double precision, intent(in) :: haloMass,deltaCritical,massResolution

    ! Initialize if necessary.
    call Tree_Branching_Initialize

    ! Interpolate in the tabulated function and return a value.
    Tree_Subresolution_Fraction=Tree_Subresolution_Fraction_Function(haloMass,deltaCritical,massResolution)
    return
  end function Tree_Subresolution_Fraction

  double precision function Tree_Branch_Mass(haloMass,deltaCritical,massResolution,probabilityFraction)
    !% Return the mass of a progenitor halo in a branch split.
    implicit none
    double precision, intent(in) :: haloMass,deltaCritical,massResolution,probabilityFraction

    ! Initialize if necessary.
    call Tree_Branching_Initialize

    ! Interpolate in the tabulated function and return a value.
    Tree_Branch_Mass=Tree_Branch_Mass_Function(haloMass,deltaCritical,massResolution,probabilityFraction)
    return
  end function Tree_Branch_Mass

  subroutine Tree_Branching_Initialize
    !% Initializes the tree branching module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="treeBranchingMethod" type="moduleUse">
    include 'merger_trees.branching_probability.modules.inc'
    !# </include>
    implicit none
 
    ! Initialize if necessary.
    !$omp critical(Tree_Branching_Initialization) 
    if (.not.treeBranchingInitialized) then
       ! Get the tree branching method parameter.
       !@ <inputParameter>
       !@   <name>treeBranchingMethod</name>
       !@   <defaultValue>modifiedPress-Schechter</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for computing merger tree branching probabilities when building merger trees.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('treeBranchingMethod',treeBranchingMethod,defaultValue='modifiedPress-Schechter')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="treeBranchingMethod" type="code" action="subroutine">
       !#  <subroutineArgs>treeBranchingMethod,Tree_Branching_Probability_Function,Tree_Subresolution_Fraction_Function,Tree_Branch_Mass_Function,Tree_Maximum_Step_Function</subroutineArgs>
       include 'merger_trees.branching_probability.inc'
       !# </include>
       if (       .not.associated(Tree_Branching_Probability_Function )  &
            & .or..not.associated(Tree_Subresolution_Fraction_Function)  &
            & .or..not.associated(Tree_Branch_Mass_Function           )  &
            & .or..not.associated(Tree_Maximum_Step_Function          )) &
            & call Galacticus_Error_Report('Tree_Branching_Initialize','method '//char(treeBranchingMethod)//' is unrecognized')
       
       treeBranchingInitialized=.true.
    end if
    !$omp end critical(Tree_Branching_Initialization)
    return
  end subroutine Tree_Branching_Initialize
  
end module Merger_Tree_Branching
