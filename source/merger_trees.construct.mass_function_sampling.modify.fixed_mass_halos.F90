!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Implements a modifier of halo mass samples which inserts a set of fixed mass halos.
 
  !# <massFunctionSamplingModifier name="massFunctionSamplingModifierFixedMassHalos">
  !#  <description>Modifies a halo mass sample by inserting a set of fixed mass halos.</description>
  !# </massFunctionSamplingModifier>

  type, extends(massFunctionSamplingModifierClass) :: massFunctionSamplingModifierFixedMassHalos
     !% A class implementing halo mass sample modification which inserts a set of fixed mass halos.
     private
     double precision, allocatable, dimension(:) :: haloMass
     integer         , allocatable, dimension(:) :: haloCount
   contains
     procedure :: modify => fixedMassHalosModify
  end type massFunctionSamplingModifierFixedMassHalos

  interface massFunctionSamplingModifierFixedMassHalos
     !% Constructors for the fixedMassHalos halo mass sample modifier class.
     module procedure fixedMassHalosDefaultConstructor
  end interface massFunctionSamplingModifierFixedMassHalos

  ! Initialization state.
  logical          :: fixedMassHalosInitialized=.false.

  ! Default settings.
  double precision, allocatable, dimension(:) :: haloMassSampleModifierFixedMassHalosMass
  integer         , allocatable, dimension(:) :: haloMassSampleModifierFixedMassHalosCount

contains

  function fixedMassHalosDefaultConstructor()
    !% Default constructor for the {\tt fixedMassHalos} halo mass sample modifier class.
    use Galacticus_Display
    use Input_Parameters
    use Galacticus_Error
    use Memory_Management
    implicit none
    type   (massFunctionSamplingModifierFixedMassHalos) :: fixedMassHalosDefaultConstructor
    integer                                             :: fixedHalosCount

    if (.not.fixedMassHalosInitialized) then
       !$omp critical (satelliteMergingTimescalesFixedMassHalosInitialize)
       if (.not.fixedMassHalosInitialized) then
          ! Find number of masses to insert.
          fixedHalosCount=-1
          if (Input_Parameter_Is_Present('haloMassSampleModifierFixedMassHalosMass' )) then
             if     (                                                                                                &
                  &   fixedHalosCount /= -1                                                                          &
                  &  .and.                                                                                           &
                  &   fixedHalosCount /= Get_Input_Parameter_Array_Size('haloMassSampleModifierFixedMassHalosMass' ) &
                  & )                                                                                                &
                  & call Galacticus_Error_Report(                                                                    &
                  &                              'fixedMassHalosDefaultConstructor',                                 &
                  &                              'parameter cardinality mismatch'                                    &
                  &                             )
             if (fixedHalosCount == -1) fixedHalosCount=Get_Input_Parameter_Array_Size('haloMassSampleModifierFixedMassHalosMass' )
          end if
          if (Input_Parameter_Is_Present('haloMassSampleModifierFixedMassHalosCount')) then
             if     (                                                                                                &
                  &   fixedHalosCount /= -1                                                                          &
                  &  .and.                                                                                           &
                  &   fixedHalosCount /= Get_Input_Parameter_Array_Size('haloMassSampleModifierFixedMassHalosCount') &
                  & )                                                                                                &
                  & call Galacticus_Error_Report(                                                                    &
                  &                              'fixedMassHalosDefaultConstructor',                                 &
                  &                              'parameter cardinality mismatch'                                    &
                  &                             )
             if (fixedHalosCount == -1) fixedHalosCount=Get_Input_Parameter_Array_Size('haloMassSampleModifierFixedMassHalosCount')
          end if
          if (fixedHalosCount == -1) fixedHalosCount=1
          call Alloc_Array(haloMassSampleModifierFixedMassHalosMass ,[fixedHalosCount])
          call Alloc_Array(haloMassSampleModifierFixedMassHalosCount,[fixedHalosCount])
          !@ <inputParameter>
          !@   <name>haloMassSampleModifierFixedMassHalosMass</name>
          !@   <defaultValue>$10^{12}M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies the masses of halos to insert into the halo mass sample when building halos.
          !@   </description>
          !@   <type>float</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('haloMassSampleModifierFixedMassHalosMass',haloMassSampleModifierFixedMassHalosMass,defaultValue=[1.0d12])
          !@ <inputParameter>
          !@   <name>haloMassSampleModifierFixedMassHalosCount</name>
          !@   <defaultValue>1</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies the number of halos to insert into the halo mass sample when building halos.
          !@   </description>
          !@   <type>float</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('haloMassSampleModifierFixedMassHalosCount',haloMassSampleModifierFixedMassHalosCount,defaultValue=[1])
          ! Record that we are now initialized.
          fixedMassHalosInitialized=.true.
       end if
       !$omp end critical (satelliteMergingTimescalesFixedMassHalosInitialize)
    end if
    fixedMassHalosDefaultConstructor=fixedMassHalosConstructor(haloMassSampleModifierFixedMassHalosMass,haloMassSampleModifierFixedMassHalosCount)
    return
  end function fixedMassHalosDefaultConstructor

  function fixedMassHalosConstructor(haloMass,haloCount)
    !% Generic constructor for the {\tt fixedMassHalos} halo mass sample modifier class.
    use Galacticus_Display
    use Input_Parameters
    implicit none
    type            (massFunctionSamplingModifierFixedMassHalos)                :: fixedMassHalosConstructor
    double precision                                            , intent(in   ), dimension(:) :: haloMass
    integer                                                     , intent(in   ), dimension(:) :: haloCount

    fixedMassHalosConstructor%haloMass =haloMass
    fixedMassHalosConstructor%haloCount=haloCount
    return
  end function fixedMassHalosConstructor

  subroutine fixedMassHalosModify(self,treeHaloMass)
    !% Modify a halo mass sample by inserting additional halos of fixed mass.
    use Memory_Management
    implicit none
    class           (massFunctionSamplingModifierFixedMassHalos)                           , intent(inout) :: self
    double precision                                            , allocatable, dimension(:), intent(inout) :: treeHaloMass
    double precision                                            , allocatable, dimension(:)                :: treeHaloMassTmp
    integer                                                                                                :: indexStart     , i

    call Move_Alloc (treeHaloMass,      treeHaloMassTmp                     )
    call Alloc_Array(treeHaloMass,shape(treeHaloMassTmp)+sum(self%haloCount))
    treeHaloMass(                      1:size(treeHaloMassTmp))=treeHaloMassTmp
    indexStart=size(treeHaloMassTmp)+1
    do i=1,size(self%haloCount)
       treeHaloMass(indexStart:indexStart+self%haloCount(i)-1)=self%haloMass(i)
       indexStart=indexStart+self%haloCount(i)
    end do
    call Dealloc_Array(treeHaloMassTmp)
    return
  end subroutine fixedMassHalosModify
  
