!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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
     double precision, allocatable, dimension(:) :: haloMass , haloRadius
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
  double precision, allocatable, dimension(:) :: haloMassSampleModifierFixedMassHalosMass     , haloMassSampleModifierFixedMassHalosRadius
  integer         , allocatable, dimension(:) :: haloMassSampleModifierFixedMassHalosCount
  logical                                     :: haloMassSampleModifierFixedMassHalosOverwrite

contains

  function fixedMassHalosDefaultConstructor()
    !% Default constructor for the {\normalfont \ttfamily fixedMassHalos} halo mass sample modifier class.
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
          if (Input_Parameter_Is_Present('haloMassSampleModifierFixedMassHalosRadius')) then
             if     (                                                                                                 &
                  &   fixedHalosCount /= -1                                                                           &
                  &  .and.                                                                                            &
                  &   fixedHalosCount /= Get_Input_Parameter_Array_Size('haloMassSampleModifierFixedMassHalosRadius') &
                  & )                                                                                                 &
                  & call Galacticus_Error_Report(                                                                     &
                  &                              'fixedMassHalosDefaultConstructor',                                  &
                  &                              'parameter cardinality mismatch'                                     &
                  &                             )
             if (fixedHalosCount == -1) fixedHalosCount=Get_Input_Parameter_Array_Size('haloMassSampleModifierFixedMassHalosRadius')
          end if
          if (fixedHalosCount == -1) fixedHalosCount=1
          call Alloc_Array(haloMassSampleModifierFixedMassHalosMass  ,[fixedHalosCount])
          call Alloc_Array(haloMassSampleModifierFixedMassHalosCount ,[fixedHalosCount])
          call Alloc_Array(haloMassSampleModifierFixedMassHalosRadius,[fixedHalosCount])
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
          call Get_Input_Parameter('haloMassSampleModifierFixedMassHalosMass',haloMassSampleModifierFixedMassHalosMass,defaultValue=spread(1.0d12,1,fixedHalosCount))
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
          call Get_Input_Parameter('haloMassSampleModifierFixedMassHalosCount',haloMassSampleModifierFixedMassHalosCount,defaultValue=spread(1,1,fixedHalosCount))
          !@ <inputParameter>
          !@   <name>haloMassSampleModifierFixedMassHalosRadius</name>
          !@   <defaultValue>1</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies the radii within which halo masses are specified when inserting fixed mass halos into the mass sample when building halos.
          !@   </description>
          !@   <type>float</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('haloMassSampleModifierFixedMassHalosRadius',haloMassSampleModifierFixedMassHalosRadius,defaultValue=spread(-1.0d0,1,fixedHalosCount))
          !@ <inputParameter>
          !@   <name>haloMassSampleModifierFixedMassHalosOverwrite</name>
          !@   <defaultValue>true</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     If {\normalfont \ttfamily true} the sample of halo masses will be overwritten instead of being added to.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('haloMassSampleModifierFixedMassHalosOverwrite',haloMassSampleModifierFixedMassHalosOverwrite,defaultValue=.false.)
          ! Record that we are now initialized.
          fixedMassHalosInitialized=.true.
       end if
       !$omp end critical (satelliteMergingTimescalesFixedMassHalosInitialize)
    end if
    fixedMassHalosDefaultConstructor=fixedMassHalosConstructor(haloMassSampleModifierFixedMassHalosMass,haloMassSampleModifierFixedMassHalosRadius,haloMassSampleModifierFixedMassHalosCount)
    return
  end function fixedMassHalosDefaultConstructor

  function fixedMassHalosConstructor(haloMass,haloRadius,haloCount)
    !% Generic constructor for the {\normalfont \ttfamily fixedMassHalos} halo mass sample modifier class.
    use Galacticus_Display
    use Input_Parameters
    implicit none
    type            (massFunctionSamplingModifierFixedMassHalos)                              :: fixedMassHalosConstructor
    double precision                                            , intent(in   ), dimension(:) :: haloMass                 , haloRadius
    integer                                                     , intent(in   ), dimension(:) :: haloCount

    fixedMassHalosConstructor%haloMass  =haloMass
    fixedMassHalosConstructor%haloRadius=haloRadius
    fixedMassHalosConstructor%haloCount =haloCount
    return
  end function fixedMassHalosConstructor

  subroutine fixedMassHalosModify(self,treeHaloMass,treeBaseTime)
    !% Modify a halo mass sample by inserting additional halos of fixed mass.
    use Memory_Management
    use Root_Finder
    use Galacticus_Calculations_Resets
    use Dark_Matter_Halo_Scales
    use Cosmology_Parameters
    use Galacticus_Nodes
    implicit none
    class           (massFunctionSamplingModifierFixedMassHalos)                           , intent(inout) :: self
    double precision                                            , allocatable, dimension(:), intent(inout) :: treeHaloMass
    double precision                                            , allocatable, dimension(:)                :: treeHaloMassTmp
    double precision                                                                       , intent(in   ) :: treeBaseTime
    type            (treeNode                                  ), pointer                                  :: node
    class           (nodeComponentBasic                        ), pointer                                  :: basic
    class           (cosmologyParametersClass                  ), pointer                                  :: cosmologyParameters_
    class           (darkMatterHaloScaleClass                  ), pointer                                  :: darkMatterHaloScale_
    integer                                                                                                :: indexStart     , i
    type            (rootFinder                                )                                           :: finder
    
    ! Get required objects
    cosmologyParameters_ => cosmologyParameters()
    darkMatterHaloScale_ => darkMatterHaloScale()
    ! Create temporary node.
    node                 => treeNode           ()
    basic                => node%basic(autoCreate=.true.)
    call basic%timeSet(treeBaseTime)
    do i=1,size(self%haloCount)
       ! Set the halo mass.
       call basic%massSet(self%haloMass(i))
       call Galacticus_Calculations_Reset(node)
       ! Set unspecified radii to virial radii.
       if (self%haloRadius(i) <= 0.0d0) self%haloRadius(i)=darkMatterHaloScale_%virialRadius(node)
       ! Convert masses to virial masses.
       call finder%tolerance(1.0d-6,1.0d-6)
       call finder%rangeExpand(rangeExpandUpward=2.0d0,rangeExpandDownward=0.5d0,rangeExpandType=rangeExpandMultiplicative)
       call finder%rootFunction(fixedEnclosedMass)
       self%haloMass(i)=finder%find(rootGuess=self%haloMass(i))
    end do
    call node%destroy()
    ! Insert halos into sample.
    if (haloMassSampleModifierFixedMassHalosOverwrite) then
       call Dealloc_Array(treeHaloMass                      )
       call Alloc_Array  (treeHaloMass,[sum(self%haloCount)])
       indexStart=1
       do i=1,size(self%haloCount)
          treeHaloMass(indexStart:indexStart+self%haloCount(i)-1)=self%haloMass(i)
          indexStart=indexStart+self%haloCount(i)
       end do       
    else
       call Move_Alloc (treeHaloMass,      treeHaloMassTmp                     )
       call Alloc_Array(treeHaloMass,shape(treeHaloMassTmp)+sum(self%haloCount))
       treeHaloMass(                      1:size(treeHaloMassTmp))=treeHaloMassTmp
       indexStart=size(treeHaloMassTmp)+1
       do i=1,size(self%haloCount)
          treeHaloMass(indexStart:indexStart+self%haloCount(i)-1)=self%haloMass(i)
          indexStart=indexStart+self%haloCount(i)
       end do
       call Dealloc_Array(treeHaloMassTmp)
    end if
    return

  contains
    
  double precision function fixedEnclosedMass(haloMass)
    !% Root finding function used to set the host halo mass for Local Group satellite mass function calculations.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    double precision, intent(in   ) :: haloMass

    call basic%massSet(haloMass)
    call Galacticus_Calculations_Reset(node)
    fixedEnclosedMass=                                                                      &
         & +Galactic_Structure_Enclosed_Mass(node,self%haloRadius(i),massType=massTypeDark) &
         & *  cosmologyParameters_%OmegaMatter()                                            &
         & /(                                                                               &
         &   +cosmologyParameters_%OmegaMatter()                                            &
         &   -cosmologyParameters_%OmegaBaryon()                                            &
         &  )                                                                               &
         & -self%haloMass(i)
    return
  end function fixedEnclosedMass

end subroutine fixedMassHalosModify
  
