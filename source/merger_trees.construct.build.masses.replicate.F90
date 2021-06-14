!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

  !% Implementation of a merger tree masses class which replicates masses from another class a specified number of times.

  !# <mergerTreeBuildMasses name="mergerTreeBuildMassesReplicate">
  !#  <description>A merger tree masses class which replicates masses from another class a specified number of times.</description>
  !# </mergerTreeBuildMasses>
  type, extends(mergerTreeBuildMassesClass) :: mergerTreeBuildMassesReplicate
     !% Implementation of a merger tree masses class which replicates masses from another class a specified number of times.
     private
     class  (mergerTreeBuildMassesClass), pointer :: mergerTreeBuildMasses_ => null()
     integer(c_size_t                  )          :: replicationCount
   contains
     final     ::              replicateDestructor
     procedure :: construct => replicateConstruct
  end type mergerTreeBuildMassesReplicate

  interface mergerTreeBuildMassesReplicate
     module procedure replicateConstructorParameters
     module procedure replicateConstructorInternal
  end interface mergerTreeBuildMassesReplicate

contains

  function replicateConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily replicate} merger tree masses class which takes a parameter set
    !% as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (mergerTreeBuildMassesReplicate)                :: self
    type   (inputParameters               ), intent(inout) :: parameters
    class  (mergerTreeBuildMassesClass    ), pointer       :: mergerTreeBuildMasses_
    integer(c_size_t                      )                :: replicationCount

    !# <inputParameter>
    !#   <name>replicationCount</name>
    !#   <source>parameters</source>
    !#   <description>The number of times to replicate each halo mass.</description>
    !# </inputParameter>
    !# <objectBuilder class="mergerTreeBuildMasses" name="mergerTreeBuildMasses_" source="parameters"/>
    self=mergerTreeBuildMassesReplicate(replicationCount,mergerTreeBuildMasses_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="mergerTreeBuildMasses_"/>
    return
  end function replicateConstructorParameters

  function replicateConstructorInternal(replicationCount,mergerTreeBuildMasses_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily replicate} merger tree masses class.
    implicit none
    type   (mergerTreeBuildMassesReplicate)                        :: self
    class  (mergerTreeBuildMassesClass    ), target, intent(in   ) :: mergerTreeBuildMasses_
    integer(c_size_t                      )        , intent(in   ) :: replicationCount
    !# <constructorAssign variables="replicationCount, *mergerTreeBuildMasses_"/>

    return
  end function replicateConstructorInternal

  subroutine replicateDestructor(self)
    !% Destructor for the merger tree mergerTreeBuildMasses function class.
    implicit none
    type(mergerTreeBuildMassesReplicate), intent(inout) :: self

    !# <objectDestructor name="self%mergerTreeBuildMasses_"/>
    return
  end subroutine replicateDestructor

  subroutine replicateConstruct(self,time,mass,massMinimum,massMaximum,weight)
    !% Construct a set of merger tree masses by sampling from a distribution.
    use            :: Galacticus_Error , only : Galacticus_Error_Report
    use, intrinsic :: ISO_C_Binding    , only : c_size_t
    use            :: Memory_Management, only : allocateArray          , deallocateArray
    implicit none
    class           (mergerTreeBuildMassesReplicate), intent(inout)                            :: self
    double precision                                , intent(in   )                            :: time
    double precision                                , intent(  out), allocatable, dimension(:) :: mass          , weight        , &
         &                                                                                        massMinimum   , massMaximum
    double precision                                               , allocatable, dimension(:) :: massTmp       , massMinimumTmp, &
         &                                                                                        massMaximumTmp, weightTmp
    integer         (c_size_t                       )                                          :: i

    call self%mergerTreeBuildMasses_%construct(time,massTmp,massMinimumTmp,massMaximumTmp,weightTmp)
    if (allocated(massTmp).and.allocated(weightTmp)) then
       call allocateArray(mass  ,shape(massTmp  )*self%replicationCount)
       call allocateArray(weight,shape(weightTmp)*self%replicationCount)
       do i=1,size(massTmp)
          mass  ((i-1)*self%replicationCount+1:i*self%replicationCount)=massTmp  (i)
          weight((i-1)*self%replicationCount+1:i*self%replicationCount)=weightTmp(i)/dble(self%replicationCount)
       end do
       call deallocateArray(massTmp  )
       call deallocateArray(weightTmp)
    else
       call Galacticus_Error_Report('masses and weights are required'//{introspection:location})
    end if
    return
  end subroutine replicateConstruct
