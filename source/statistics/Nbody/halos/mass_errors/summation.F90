!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Implements an N-body dark matter halo mass error class which sums over other error models.
!!}

  !![
  <nbodyHaloMassError name="nbodyHaloMassErrorSummation" docformat="rst">
   <description>
   An N-body dark matter halo mass error class which combines any number of other error models, assumed to be independent. The fractional mass errors of the individual models are summed in quadrature, :math:`\sigma^2 = \sum_i \sigma_i^2`, and the correlation between the mass errors of a pair of halos is :math:`C_{12} = \sum_i C_{i,12} \sigma_{i,1} \sigma_{i,2} / \sigma_1 \sigma_2`. When evaluating the correlation for a pair of halos at different times, the later halo is used as the reference halo for the earlier halo (see :galacticus-class:`nbodyHaloMassErrorTreeConstruction`).
   </description>
   <linkedList type="nbodyHaloMassErrorList" variable="nbodyHaloMassErrors" next="next" object="nbodyHaloMassError_" objectType="nbodyHaloMassErrorClass"/>
  </nbodyHaloMassError>
  !!]

  type, public :: nbodyHaloMassErrorList
     class(nbodyHaloMassErrorClass), pointer :: nbodyHaloMassError_ => null()
     type (nbodyHaloMassErrorList ), pointer :: next                => null()
  end type nbodyHaloMassErrorList

  type, extends(nbodyHaloMassErrorClass) :: nbodyHaloMassErrorSummation
     !!{RST
     An N-body halo mass error class which sums over other error models.
     !!}
     private
     type(nbodyHaloMassErrorList), pointer :: nbodyHaloMassErrors => null()
   contains
     final     ::                    summationDestructor
     procedure :: errorFractional => summationErrorFractional
     procedure :: correlation     => summationCorrelation
     procedure :: errorZeroAlways => summationErrorZeroAlways
  end type nbodyHaloMassErrorSummation

  interface nbodyHaloMassErrorSummation
     !!{RST
     Constructors for the :galacticus-class:`nbodyHaloMassErrorSummation` N-body halo mass error class.
     !!}
     module procedure summationConstructorParameters
     module procedure summationConstructorInternal
  end interface nbodyHaloMassErrorSummation

contains

  function summationConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nbodyHaloMassErrorSummation` N-body halo mass error class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nbodyHaloMassErrorSummation), target        :: self
    type   (inputParameters            ), intent(inout) :: parameters
    type   (nbodyHaloMassErrorList     ), pointer       :: member
    integer                                             :: i

    member => null()
    do i=1,parameters%copiesCount('nbodyHaloMassError',zeroIfNotPresent=.true.)
       if (associated(member)) then
          allocate(member%next)
          member => member%next
       else
          allocate(self%nbodyHaloMassErrors)
          member => self%nbodyHaloMassErrors
       end if
       !![
       <objectBuilder class="nbodyHaloMassError" name="member%nbodyHaloMassError_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="nbodyHaloMassError"/>
    !!]
    return
  end function summationConstructorParameters

  function summationConstructorInternal(nbodyHaloMassErrors) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nbodyHaloMassErrorSummation` N-body halo mass error class.
    !!}
    implicit none
    type(nbodyHaloMassErrorSummation)                        :: self
    type(nbodyHaloMassErrorList     ), target, intent(in   ) :: nbodyHaloMassErrors
    type(nbodyHaloMassErrorList     ), pointer               :: member

    self%nbodyHaloMassErrors => nbodyHaloMassErrors
    member                   => nbodyHaloMassErrors
    do while (associated(member))
       !![
       <referenceCountIncrement owner="member" object="nbodyHaloMassError_"/>
       !!]
       member => member%next
    end do
    return
  end function summationConstructorInternal

  subroutine summationDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nbodyHaloMassErrorSummation` N-body halo mass error class.
    !!}
    implicit none
    type(nbodyHaloMassErrorSummation), intent(inout) :: self
    type(nbodyHaloMassErrorList     ), pointer       :: member, memberNext

    member => self%nbodyHaloMassErrors
    do while (associated(member))
       memberNext => member%next
       !![
       <objectDestructor name="member%nbodyHaloMassError_"/>
       !!]
       deallocate(member)
       member => memberNext
    end do
    return
  end subroutine summationDestructor

  double precision function summationErrorFractional(self,node,nodeReference)
    !!{RST
    Return the fractional error on the mass of an N-body halo, summing the errors of all member models in quadrature.
    !!}
    implicit none
    class(nbodyHaloMassErrorSummation), intent(inout)           :: self
    type (treeNode                   ), intent(inout)           :: node
    type (treeNode                   ), intent(inout), optional :: nodeReference
    type (nbodyHaloMassErrorList     ), pointer                 :: member

    summationErrorFractional=0.0d0
    member => self%nbodyHaloMassErrors
    do while (associated(member))
       summationErrorFractional=+summationErrorFractional                                          &
            &                   +member%nbodyHaloMassError_%errorFractional(node,nodeReference)**2
       member => member%next
    end do
    summationErrorFractional=sqrt(summationErrorFractional)
    return
  end function summationErrorFractional

  double precision function summationCorrelation(self,node1,node2)
    !!{RST
    Return the correlation of the masses of a pair of N-body halos, combining the covariances of all member models. For halos at
    different times the later halo is used as the reference halo for the earlier halo, consistent with the errors evaluated by
    consumers which trace progenitor halos back from their parent halos.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (nbodyHaloMassErrorSummation), intent(inout) :: self
    type            (treeNode                   ), intent(inout) :: node1     , node2
    type            (nbodyHaloMassErrorList     ), pointer       :: member
    class           (nodeComponentBasic         ), pointer       :: basic1    , basic2
    double precision                                             :: error1    , error2   , &
         &                                                          variance1 , variance2, &
         &                                                          covariance

    basic1     => node1%basic()
    basic2     => node2%basic()
    variance1  =  0.0d0
    variance2  =  0.0d0
    covariance =  0.0d0
    member     => self%nbodyHaloMassErrors
    do while (associated(member))
       if      (basic1%time() < basic2%time()) then
          error1=member%nbodyHaloMassError_%errorFractional(node1,nodeReference=node2)
          error2=member%nbodyHaloMassError_%errorFractional(node2                    )
       else if (basic2%time() < basic1%time()) then
          error1=member%nbodyHaloMassError_%errorFractional(node1                    )
          error2=member%nbodyHaloMassError_%errorFractional(node2,nodeReference=node1)
       else
          error1=member%nbodyHaloMassError_%errorFractional(node1                    )
          error2=member%nbodyHaloMassError_%errorFractional(node2                    )
       end if
       variance1=+variance1+error1**2
       variance2=+variance2+error2**2
       if (error1 > 0.0d0 .and. error2 > 0.0d0)                               &
            & covariance=+covariance                                          &
            &            +member%nbodyHaloMassError_%correlation(node1,node2) &
            &            *error1                                              &
            &            *error2
       member => member%next
    end do
    if (variance1 > 0.0d0 .and. variance2 > 0.0d0) then
       ! Limit to the range [-1,1] to guard against rounding errors.
       summationCorrelation=max(-1.0d0,min(+1.0d0,covariance/sqrt(variance1*variance2)))
    else
       summationCorrelation=0.0d0
    end if
    return
  end function summationCorrelation

  logical function summationErrorZeroAlways(self)
    !!{RST
    Return true if the errors of all member models are always zero.
    !!}
    implicit none
    class(nbodyHaloMassErrorSummation), intent(inout) :: self
    type (nbodyHaloMassErrorList     ), pointer       :: member

    summationErrorZeroAlways=.true.
    member => self%nbodyHaloMassErrors
    do while (associated(member))
       summationErrorZeroAlways=member%nbodyHaloMassError_%errorZeroAlways()
       if (.not.summationErrorZeroAlways) return
       member => member%next
    end do
    return
  end function summationErrorZeroAlways
