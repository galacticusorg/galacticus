!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

!!{
Implements an orbital velocity output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorVelocityOrbital">
   <description>
    An orbital velocity output analysis property extractor class. Specifically, the orbital velocity is defined relative to the
    top-level halo in any sub-halo hierarchy. That is, relative to the host halo which is itself not a sub-halo of any other
    halo. If the velocity of a (sub)$^i$-halo with respect to the center of its (sub)$^{i-1}$-halo host is $\mathbf{v}_i$ then
    the orbital velocity computed by this class is
    \begin{equation}
     \mathbf{v} = sum_{i=1}^N \mathbf{v}_i,
    \end{equation}
    where $N$ is the depth of the node in the sub-halo hierarchy.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorVelocityOrbital
     !!{
     An orbital velocity property extractor output analysis class.
     !!}
     private
   contains
     procedure :: elementCount => velocityOrbitalElementCount
     procedure :: extract      => velocityOrbitalExtract
     procedure :: names        => velocityOrbitalNames
     procedure :: descriptions => velocityOrbitalDescriptions
     procedure :: unitsInSI    => velocityOrbitalUnitsInSI
  end type nodePropertyExtractorVelocityOrbital

  interface nodePropertyExtractorVelocityOrbital
     !!{
     Constructors for the \refClass{nodePropertyExtractorVelocityOrbital} output analysis class.
     !!}
     module procedure velocityOrbitalConstructorParameters
  end interface nodePropertyExtractorVelocityOrbital

contains

  function velocityOrbitalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorVelocityOrbital} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorVelocityOrbital)                :: self
    type(inputParameters                     ), intent(inout) :: parameters

    self=nodePropertyExtractorVelocityOrbital()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function velocityOrbitalConstructorParameters

  integer function velocityOrbitalElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily velocityOrbital} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorVelocityOrbital), intent(inout) :: self
    double precision                                      , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    velocityOrbitalElementCount=3
    return
  end function velocityOrbitalElementCount

  function velocityOrbitalExtract(self,node,time,instance) result(velocity)
    !!{
    Implement a velocityOrbital output analysis.
    !!}
    use :: Galacticus_Nodes    , only : nodeComponentSatellite, nodeComponentBasic
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    double precision                                      , dimension(:) , allocatable :: velocity
    class           (nodePropertyExtractorVelocityOrbital), intent(inout), target      :: self
    type            (treeNode                            ), intent(inout), target      :: node
    double precision                                      , intent(in   )              :: time
    type            (multiCounter                        ), intent(inout), optional    :: instance
    type            (treeNode                            ), pointer                    :: nodeWork
    class           (nodeComponentBasic                  ), pointer                    :: basic
    class           (nodeComponentSatellite              ), pointer                    :: satellite
    !$GLC attributes unused :: self, instance, time

    allocate(velocity(3))
    velocity =  0.0d0
    nodeWork => node
    ! Walk up through all host halos of this node, accumulating velocity offsets from the host node center.
    do while (associated(nodeWork))
       basic     =>  nodeWork %basic    ()
       satellite =>  nodeWork %satellite()
       velocity  =  +          velocity   &
            &       +satellite%velocity ()
       if (nodeWork%isSatellite()) then
          ! Current node is a satellite, simply move to its parent.
          nodeWork => nodeWork%parent
       else
          ! Node is a host halo.          
          if (nodeWork%isOnMainBranch()) then
             ! We are on the main branch - so we are done.
             nodeWork => null()
          else
             ! We are not on the main branch - find the host halo at this time time on the branch with which we will merge.
             do while (nodeWork%isPrimaryProgenitor())
                nodeWork => nodeWork%parent
             end do
             nodeWork => nodeWork%parent
             basic    => nodeWork%basic ()
             do while (associated(nodeWork%firstChild).and..not.Values_Agree(basic%time(),time,relTol=1.0d-6))
                nodeWork => nodeWork%firstChild
                basic    => nodeWork%basic     ()
             end do
          end if
       end if
    end do
    return
  end function velocityOrbitalExtract
  
  subroutine velocityOrbitalNames(self,time,names)
    !!{
    Return the name of the velocityOrbital property.
    !!}
    implicit none
    class(nodePropertyExtractorVelocityOrbital), intent(inout)                             :: self
    double precision                           , intent(in   )                             :: time
    type (varying_string                      ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(3))
    names(1)=var_str('velocityOrbitalX')
    names(2)=var_str('velocityOrbitalY')
    names(3)=var_str('velocityOrbitalZ')
    return
  end subroutine velocityOrbitalNames

  subroutine velocityOrbitalDescriptions(self,time,descriptions)
    !!{
    Return a description of the velocityOrbital property.
    !!}
    implicit none
    class(nodePropertyExtractorVelocityOrbital), intent(inout)                            :: self
    double precision                           , intent(in   )                            :: time
    type (varying_string                      ), intent(inout), dimension(:), allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(3))
    descriptions(1)=var_str('The orbital x-velocity of the halo relative to the top-level host halo (i.e. the host which is not a sub-halo of any other halo).')
    descriptions(2)=var_str('The orbital y-velocity of the halo relative to the top-level host halo (i.e. the host which is not a sub-halo of any other halo).')
    descriptions(3)=var_str('The orbital z-velocity of the halo relative to the top-level host halo (i.e. the host which is not a sub-halo of any other halo).')
    return
  end subroutine velocityOrbitalDescriptions

  function velocityOrbitalUnitsInSI(self,time)
    !!{
    Return the units of the velocityOrbital property in the SI system.
    !!}
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    double precision                                      , dimension(:) , allocatable :: velocityOrbitalUnitsInSI
    class           (nodePropertyExtractorVelocityOrbital), intent(inout)              :: self
    double precision                                      , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(velocityOrbitalUnitsInSI(3))
    velocityOrbitalUnitsInSI=kilo
    return
  end function velocityOrbitalUnitsInSI
