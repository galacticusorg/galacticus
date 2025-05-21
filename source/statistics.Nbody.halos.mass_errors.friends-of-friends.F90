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
  Implements an N-body dark matter halo mass error class using a fit appropriate for friends-of-friends
  group finders.
  !!}

  !![
  <nbodyHaloMassError name="nbodyHaloMassErrorFriendsOfFriends">
   <description>An N-body dark matter halo mass error class which uses a fit appropriate for friends-of-friends group finders.</description>
  </nbodyHaloMassError>
  !!]
  type, extends(nbodyHaloMassErrorPowerLaw) :: nbodyHaloMassErrorFriendsOfFriends
     !!{
     An N-body halo mass error class which uses a fit appropriate for friends-of-friends group finders.
     !!}
     private
     double precision :: massParticle
   contains
  end type nbodyHaloMassErrorFriendsOfFriends

  interface nbodyHaloMassErrorFriendsOfFriends
     !!{
     Constructors for the \refClass{nbodyHaloMassErrorFriendsOfFriends} N-body halo mass error class.
     !!}
     module procedure nbodyHaloMassErrorFriendsOfFriendsParameters
     module procedure nbodyHaloMassErrorFriendsOfFriendsInternal
  end interface nbodyHaloMassErrorFriendsOfFriends

contains

  function nbodyHaloMassErrorFriendsOfFriendsParameters(parameters)
    !!{
    Constructor for the \refClass{nbodyHaloMassErrorFriendsOfFriends} N-body halo mass error class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyHaloMassErrorFriendsOfFriends)                :: nbodyHaloMassErrorFriendsOfFriendsParameters
    type            (inputParameters                   ), intent(inout) :: parameters
    double precision                                                    :: massParticle

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>massParticle</name>
      <source>parameters</source>
      <variable>massParticle</variable>
      <description>The mass of the particle in the N-body simulation in which friends-of-friends groups were found.</description>
    </inputParameter>
    !!]
    nbodyHaloMassErrorFriendsOfFriendsParameters=nbodyHaloMassErrorFriendsOfFriends(massParticle)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nbodyHaloMassErrorFriendsOfFriendsParameters

  function nbodyHaloMassErrorFriendsOfFriendsInternal(massParticle)
    !!{
    Internal constructor for the \refClass{nbodyHaloMassErrorFriendsOfFriends} N-body halo mass error class.
    !!}
    implicit none
    type            (nbodyHaloMassErrorFriendsOfFriends)                :: nbodyHaloMassErrorFriendsOfFriendsInternal
    double precision                                    , intent(in   ) :: massParticle
    double precision                                    , parameter     :: exponent         =-0.5d00
    double precision                                    , parameter     :: normalization    =+1.25d00
    double precision                                    , parameter     :: errorHighMass    =+0.022d00
    !![
    <constructorAssign variables="massParticle"/>
    !!]
    
    ! Convert from a model defined in terms of particle number (in which the fractional error is 1.2/sqrt(N)) to one defined in
    ! terms of halo mass as used in our parent class.
    nbodyHaloMassErrorFriendsOfFriendsInternal%normalizationSquared          =(normalization*(massReference/massParticle)**exponent)**2
    nbodyHaloMassErrorFriendsOfFriendsInternal%exponent                      =exponent
    nbodyHaloMassErrorFriendsOfFriendsInternal%fractionalErrorHighMassSquared=errorHighMass                                         **2
    return
  end function nbodyHaloMassErrorFriendsOfFriendsInternal
