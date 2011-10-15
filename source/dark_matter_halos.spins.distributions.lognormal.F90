!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements a lognormal halo spin distribution.

module Halo_Spin_Distributions_Lognormal
  !% Implements a lognormal halo spin distribution.
  use FGSL
  use Tree_Nodes
  implicit none
  private
  public :: Halo_Spin_Distribution_Lognormal_Initialize, Halo_Spin_Distribution_Lognormal_Snapshot,&
       & Halo_Spin_Distribution_Lognormal_State_Store, Halo_Spin_Distribution_Lognormal_State_Retrieve

  ! Parameters of the spin distribution.
  double precision :: lognormalSpinDistributionMedian,lognormalSpinDistributionSigma

  ! Random number objects
  type(fgsl_rng) :: randomSequenceObject,clonedPseudoSequenceObject
  logical        :: resetRandomSequence=.true.,resetRandomSequenceSnapshot
  !$omp threadprivate(resetRandomSequence,randomSequenceObject)

contains

  !# <haloSpinDistributionMethod>
  !#  <unitName>Halo_Spin_Distribution_Lognormal_Initialize</unitName>
  !# </haloSpinDistributionMethod>
  subroutine Halo_Spin_Distribution_Lognormal_Initialize(haloSpinDistributionMethod,Halo_Spin_Sample_Get)
    !% Initializes the ``Lognormal'' halo spin distribution module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: haloSpinDistributionMethod
    procedure(double precision), pointer, intent(inout) :: Halo_Spin_Sample_Get
    
    if (haloSpinDistributionMethod == 'lognormal') then
       Halo_Spin_Sample_Get => Halo_Spin_Distribution_Lognormal
       !@ <inputParameter>
       !@   <name>lognormalSpinDistributionMedian</name>
       !@   <defaultValue>0.03687 \citep{bett_spin_2007}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The median in a lognormal halo spin distribution.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('lognormalSpinDistributionMedian',lognormalSpinDistributionMedian,defaultValue=0.03687d0)
       !@ <inputParameter>
       !@   <name>lognormalSpinDistributionSigma</name>
       !@   <defaultValue>0.2216 \citep{bett_spin_2007}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The dispersion in a lognormal halo spin distribution.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('lognormalSpinDistributionSigma' ,lognormalSpinDistributionSigma ,defaultValue=0.51025d0)
       lognormalSpinDistributionMedian=dlog(lognormalSpinDistributionMedian)
    end if
    return
  end subroutine Halo_Spin_Distribution_Lognormal_Initialize

  double precision function Halo_Spin_Distribution_Lognormal(thisNode)
    !% Return a halo spin from a lognormal distribution.
    use Tree_Nodes
    use Gaussian_Random
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: randomDeviate,logLambda

    randomDeviate=Gaussian_Random_Get(randomSequenceObject,lognormalSpinDistributionSigma,resetRandomSequence)
    logLambda=lognormalSpinDistributionMedian+randomDeviate
    Halo_Spin_Distribution_Lognormal=dexp(logLambda)
    return
  end function Halo_Spin_Distribution_Lognormal

  !# <galacticusStateSnapshotTask>
  !#  <unitName>Halo_Spin_Distribution_Lognormal_Snapshot</unitName>
  !# </galacticusStateSnapshotTask>
  subroutine Halo_Spin_Distribution_Lognormal_Snapshot
    !% Store a snapshot of the random number generator internal state.
    implicit none

    if (.not.resetRandomSequence) clonedPseudoSequenceObject=FGSL_Rng_Clone(randomSequenceObject)
    resetRandomSequenceSnapshot=resetRandomSequence
    return
  end subroutine Halo_Spin_Distribution_Lognormal_Snapshot
  
  !# <galacticusStateStoreTask>
  !#  <unitName>Halo_Spin_Distribution_Lognormal_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Halo_Spin_Distribution_Lognormal_State_Store(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    use Pseudo_Random
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) resetRandomSequenceSnapshot
    if (.not.resetRandomSequenceSnapshot) call Pseudo_Random_Store(clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine Halo_Spin_Distribution_Lognormal_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Halo_Spin_Distribution_Lognormal_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Halo_Spin_Distribution_Lognormal_State_Retrieve(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    use Pseudo_Random
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    read (stateFile) resetRandomSequence
    if (.not.resetRandomSequence) call Pseudo_Random_Retrieve(randomSequenceObject,fgslStateFile)
    return
  end subroutine Halo_Spin_Distribution_Lognormal_State_Retrieve
    
end module Halo_Spin_Distributions_Lognormal
