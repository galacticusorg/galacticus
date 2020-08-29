!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
!!    Andrew Benson <abenson@carnegiescience.edu>
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

  !% Implementation of a merger tree masses class which uses a fixed mass for trees.
  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !# <mergerTreeBuildMasses name="mergerTreeBuildMassesFixedMass">
  !#  <description>A merger tree masses class which uses a fixed mass for trees.</description>
  !# </mergerTreeBuildMasses>
  type, extends(mergerTreeBuildMassesClass) :: mergerTreeBuildMassesFixedMass
     !% Implementation of a merger tree masses class which samples masses from a distribution.
     private
     class           (cosmologyParametersClass), pointer                   :: cosmologyParameters_   => null()
     class           (darkMatterHaloScaleClass), pointer                   :: darkMatterHaloScale_   => null()
     double precision                          , allocatable, dimension(:) :: massTree                        , radiusTree
     integer                                   , allocatable, dimension(:) :: treeCount
     double precision                                                      :: massIntervalFractional
   contains
     final     ::              fixedMassDestructor
     procedure :: construct => fixedMassConstruct
  end type mergerTreeBuildMassesFixedMass

  interface mergerTreeBuildMassesFixedMass
     module procedure fixedMassConstructorParameters
     module procedure fixedMassConstructorInternal
  end interface mergerTreeBuildMassesFixedMass

contains

  function fixedMassConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily fixedMass} merger tree masses class which takes a parameter set as
    !% input.
    use :: Galacticus_Error , only : Galacticus_Error_Report
    use :: Input_Parameters , only : inputParameter         , inputParameters
    use :: Memory_Management, only : allocateArray
    implicit none
    type            (mergerTreeBuildMassesFixedMass)                            :: self
    type            (inputParameters               ), intent(inout)             :: parameters
    double precision                                , allocatable, dimension(:) :: massTree              , radiusTree
    integer                                         , allocatable, dimension(:) :: treeCount
    class           (cosmologyParametersClass      ), pointer                   :: cosmologyParameters_
    class           (darkMatterHaloScaleClass      ), pointer                   :: darkMatterHaloScale_
    double precision                                                            :: massIntervalFractional
    integer                                                                     :: fixedHalosCount

    fixedHalosCount=-1
    if (parameters%isPresent('massTree' )) then
       if     (                                                                  &
            &   fixedHalosCount /= -1                                            &
            &  .and.                                                             &
            &   fixedHalosCount /= parameters%count('massTree'  )                &
            & )                                                                  &
            & call Galacticus_Error_Report(                                      &
            &                              'parameter cardinality mismatch'//    &
            &                              {introspection:location}              &
            &                             )
       if (fixedHalosCount == -1) fixedHalosCount=parameters%count('massTree' )
    end if
    if (parameters%isPresent('treeCount')) then
       if     (                                                                  &
            &   fixedHalosCount /= -1                                            &
            &  .and.                                                             &
            &   fixedHalosCount /= parameters%count('treeCount' )                &
            & )                                                                  &
            & call Galacticus_Error_Report(                                      &
            &                              'parameter cardinality mismatch'//    &
            &                              {introspection:location}              &
            &                             )
       if (fixedHalosCount == -1) fixedHalosCount=parameters%count('treeCount')
    end if
    if (parameters%isPresent('radiusTree')) then
       if     (                                                                  &
            &   fixedHalosCount /= -1                                            &
            &  .and.                                                             &
            &   fixedHalosCount /= parameters%count('radiusTree')                &
            & )                                                                  &
            & call Galacticus_Error_Report(                                      &
            &                              'parameter cardinality mismatch'//    &
            &                              {introspection:location}              &
            &                             )
       if (fixedHalosCount == -1) fixedHalosCount=parameters%count('radiusTree')
    end if
    if (fixedHalosCount == -1) fixedHalosCount=1
    call allocateArray(massTree  ,[fixedHalosCount])
    call allocateArray(treeCount ,[fixedHalosCount])
    call allocateArray(radiusTree,[fixedHalosCount])
    !# <inputParameter>
    !#   <name>massTree</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>spread(1.0d12,1,fixedHalosCount)</defaultValue>
    !#   <description>Specifies the masses of halos to use when building halos.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>treeCount</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>spread(1,1,fixedHalosCount)</defaultValue>
    !#   <description>Specifies the number of halos to use when building halos.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>radiusTree</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>spread(-1.0d0,1,fixedHalosCount)</defaultValue>
    !#   <description>Specifies the radii within which halo masses are specified when building halos.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massIntervalFractional</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <description>The fractional mass interval occupied by the trees. Where the intervals of trees of different mass would overlap this interval will be truncated.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=mergerTreeBuildMassesFixedMass(massTree,radiusTree,treeCount,massIntervalFractional,cosmologyParameters_,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="darkMatterHaloScale_"/>
    return
  end function fixedMassConstructorParameters

  function fixedMassConstructorInternal(massTree,radiusTree,treeCount,massIntervalFractional,cosmologyParameters_,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily fixedMass} merger tree masses class.
    implicit none
    type            (mergerTreeBuildMassesFixedMass)                              :: self
    double precision                                , intent(in   ), dimension(:) :: massTree              , radiusTree
    integer                                         , intent(in   ), dimension(:) :: treeCount
    double precision                                , intent(in   )               :: massIntervalFractional
    class           (cosmologyParametersClass      ), intent(in   ), target       :: cosmologyParameters_
    class           (darkMatterHaloScaleClass      ), intent(in   ), target       :: darkMatterHaloScale_
    !# <constructorAssign variables="massTree, radiusTree, treeCount, massIntervalFractional, *cosmologyParameters_, *darkMatterHaloScale_"/>

    return
  end function fixedMassConstructorInternal

  subroutine fixedMassDestructor(self)
    !% Destructor for the {\normalfont \ttfamily fixedMass} merger tree masses class.
    implicit none
    type(mergerTreeBuildMassesFixedMass), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    return
  end subroutine fixedMassDestructor

  subroutine fixedMassConstruct(self,time,mass,massMinimum,massMaximum,weight)
    !% Construct a set of merger tree masses by sampling from a distribution.
    use :: Galacticus_Calculations_Resets, only : Galacticus_Calculations_Reset
    use :: Galacticus_Nodes              , only : nodeComponentBasic           , treeNode
    use :: Memory_Management             , only : allocateArray
    use :: Root_Finder                   , only : rangeExpandMultiplicative    , rootFinder
    use :: Sorting                       , only : sort
    implicit none
    class           (mergerTreeBuildMassesFixedMass), intent(inout)                            :: self
    double precision                                , intent(in   )                            :: time
    double precision                                , intent(  out), allocatable, dimension(:) :: mass       , weight     , &
         &                                                                                        massMinimum, massMaximum
    type            (rootFinder                    )                                           :: finder
    type            (treeNode                      ), pointer                                  :: node
    class           (nodeComponentBasic            ), pointer                                  :: basic
    integer                                                                                    :: indexStart , i
    !$GLC attributes unused :: weight

    node  => treeNode      (                 )
    basic => node    %basic(autoCreate=.true.)
    call basic%timeSet            (time)
    call basic%timeLastIsolatedSet(time)
    do i=1,size(self%treeCount)
       if (self%radiusTree(i) > 0.0d0) then
          ! Set the halo mass.
          call Galacticus_Calculations_Reset(node)
          ! Convert masses to virial masses.
          call finder%tolerance   (                                               &
               &                   toleranceAbsolute  =1.0d-6                   , &
               &                   toleranceRelative  =1.0d-6                     &
               &                  )
          call finder%rangeExpand (                                               &
               &                   rangeExpandUpward  =2.0d+0                   , &
               &                   rangeExpandDownward=0.5d+0                   , &
               &                   rangeExpandType    =rangeExpandMultiplicative  &
               &                  )
          call finder%rootFunction(                                               &
               &                                       massEnclosed               &
               &                  )
          self%massTree(i)=finder%find(rootGuess=self%massTree(i))
       end if
    end do
    call node%destroy()
    deallocate(node)
    call allocateArray(mass       ,[sum(self%treeCount)])
    call allocateArray(massMinimum,[sum(self%treeCount)])
    call allocateArray(massMaximum,[sum(self%treeCount)])
    indexStart=1
    do i=1,size(self%treeCount)
       mass      (indexStart:indexStart+self%treeCount(i)-1)=+self%massTree (i)
       indexStart                                           =+self%treeCount(i) &
            &                                                +indexStart
    end do
    call sort(mass)
    indexStart=1
    do i=1,size(self%treeCount)
       massMinimum(indexStart:indexStart+self%treeCount(i)-1)=+self%massTree (i)/sqrt(1.0d0+self%massIntervalFractional)
       massMaximum(indexStart:indexStart+self%treeCount(i)-1)=+self%massTree (i)*sqrt(1.0d0+self%massIntervalFractional)
       if (i > 1 .and. massMinimum(indexStart) < massMaximum(indexStart-1)) then
          massMinimum(indexStart                    :indexStart+self%treeCount(i)-1)=sqrt(                    &
               &                                                                          +self%massTree(i-1) &
               &                                                                          *self%massTree(i  ) &
               &                                                                         )
          massMaximum(indexStart-self%treeCount(i-1):indexStart                  -1)=sqrt(                    &
               &                                                                          +self%massTree(i-1) &
               &                                                                          *self%massTree(i  ) &
               &                                                                         )
       end if
       indexStart=+self%treeCount(i)+indexStart
    end do
    return

  contains

    double precision function massEnclosed(massTree)
      !% Root finding function used to set the halo mass given the halo radius.
      use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
      use :: Galactic_Structure_Options        , only : massTypeDark
      implicit none
      double precision, intent(in   ) :: massTree

      call basic%massSet(massTree)
      call Galacticus_Calculations_Reset(node)
      massEnclosed=+Galactic_Structure_Enclosed_Mass(                              &
           &                                                   node              , &
           &                                                   self%radiusTree(i), &
           &                                          massType=massTypeDark        &
           &                                         )                             &
           &       *  self%cosmologyParameters_%OmegaMatter()                      &
           &       /(                                                              &
           &         +self%cosmologyParameters_%OmegaMatter()                      &
           &         -self%cosmologyParameters_%OmegaBaryon()                      &
           &        )                                                              &
           &       -self%massTree(i)
      return
    end function massEnclosed

  end subroutine fixedMassConstruct
