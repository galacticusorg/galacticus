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

  !% An implementation of dark matter halo virial density contrasts based on a friends-of-friends linking length.

  !# <virialDensityContrast name="virialDensityContrastFriendsOfFriends">
  !#  <description>Dark matter halo virial density contrasts based on the friends-of-friends algorithm linking length.</description>
  !# </virialDensityContrast>

  type, extends(virialDensityContrastClass) :: virialDensityContrastFriendsOfFriends
     !% A dark matter halo virial density contrast class based on the friends-of-friends algorithm linking length.
     private
     double precision :: linkingLength, densityRatio
   contains
     procedure :: densityContrast             => friendsOfFriendsDensityContrast
     procedure :: densityContrastRateOfChange => friendsOfFriendsDensityContrastRateOfChange
  end type virialDensityContrastFriendsOfFriends

  interface virialDensityContrastFriendsOfFriends
     !% Constructors for the {\tt friendsOfFriends} dark matter halo virial density contrast class.
     module procedure friendsOfFriendsDefaultConstructor
     module procedure friendsOfFriendsConstructor
  end interface virialDensityContrastFriendsOfFriends

  ! Initialization state.
  logical          :: fofInitialized=.false.

  ! Default value of the linking length and density ratio parameters.
  double precision :: virialDensityContrastFoFLinkingLength, virialDensityContrastFoFDensityRatio

contains

  function friendsOfFriendsDefaultConstructor()
    !% Default constructor for the {\tt friendsOfFriends} dark matter halo virial density contrast class.
    use Input_Parameters
    implicit none
    type (virialDensityContrastFriendsOfFriends), target  :: friendsOfFriendsDefaultConstructor
    
    if (.not.fofInitialized) then
       !$omp critical(virialDensityContrastFriendsOfFriendsInitialize)
       if (.not.fofInitialized) then
          ! Get the linking length to use.
          !@ <inputParameter>
          !@   <name>virialDensityContrastFoFLinkingLength</name>
          !@   <defaultValue>0.2</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The friends-of-friends linking length algorithm to use in computing virial density contrast.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("virialDensityContrastFoFLinkingLength",virialDensityContrastFoFLinkingLength,defaultValue=0.2d0)
          !@ <inputParameter>
          !@   <name>virialDensityContrastFoFDensityRatio</name>
          !@   <defaultValue>4.688 (value appropriate for an \gls{nfw} profile with concentration $c=6.88$ which is the concentration found by \cite{prada_halo_2011} for halos with $\sigma=1.686$ which is the approximate critical overdensity for collapse).</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The ratio of mean virial density to density at the virial radius to assume when setting virial density contrasts in the friends-of-friends model.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("virialDensityContrastFoFDensityRatio",virialDensityContrastFoFDensityRatio,defaultValue=4.688d0)
          ! Record initialization.
          fofInitialized=.true.
       end if
       !$omp end critical(virialDensityContrastFriendsOfFriendsInitialize)
    end if
    friendsOfFriendsDefaultConstructor=friendsOfFriendsConstructor(virialDensityContrastFoFLinkingLength,virialDensityContrastFoFDensityRatio)
    return
  end function friendsOfFriendsDefaultConstructor

  function friendsOfFriendsConstructor(linkingLength,densityRatio)
    !% Generic constructor for the {\tt friendsOfFriends} dark matter halo virial density contrast class.
    use Input_Parameters
    implicit none
    type            (virialDensityContrastFriendsOfFriends), target        :: friendsOfFriendsConstructor
    double precision                                       , intent(in   ) :: linkingLength              , densityRatio

    friendsOfFriendsConstructor%linkingLength=linkingLength
    friendsOfFriendsConstructor%densityRatio =densityRatio
    return
  end function friendsOfFriendsConstructor

  double precision function friendsOfFriendsDensityContrast(self,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, based on the friends-of-friends algorithm linking length.
    use Numerical_Constants_Math
    implicit none
    class           (virialDensityContrastFriendsOfFriends), intent(inout)           :: self
    double precision                                       , intent(in   ), optional :: time                          , expansionFactor
    logical                                                , intent(in   ), optional :: collapsing
    double precision                                                                 :: boundingSurfaceDensityContrast

    boundingSurfaceDensityContrast =3.0d0/2.0d0/Pi/self%linkingLength**3
    friendsOfFriendsDensityContrast=self%densityRatio*boundingSurfaceDensityContrast
    return
  end function friendsOfFriendsDensityContrast

  double precision function friendsOfFriendsDensityContrastRateOfChange(self,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, based on the friends-of-friends algorithm linking length.
    use Numerical_Constants_Math
    implicit none
    class           (virialDensityContrastFriendsOfFriends), intent(inout)           :: self
    double precision                                       , intent(in   ), optional :: time      , expansionFactor
    logical                                                , intent(in   ), optional :: collapsing

    friendsOfFriendsDensityContrastRateOfChange=0.0d0
    return
  end function friendsOfFriendsDensityContrastRateOfChange
