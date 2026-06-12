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

  !!{RST
  An implementation of dark matter halo virial density contrasts based on a friends-of-friends linking length.
  !!}

  !![
  <virialDensityContrast name="virialDensityContrastFriendsOfFriends" docformat="rst">
   <description>
   Dark matter halo virial density contrasts based on the friends-of-friends algorithm linking length. The virial density contrast is computed consistent with a given friends-of-friends algorithm linking length. According to :cite:t:`lacey_merger_1994`, the friends-of-friends algorithm selects objects bounded by an isodensity contour of density contrast

   .. math::

      \Delta_\mathrm{iso} = {3 \over 2 \pi b^3},

   where :math:`b=`\ ``[virialDensityContrastFoFLinkingLength]`` is the dimensionless linking length of the algorithm (i.e. the linking length in units of the mean interparticle spacing). The virial density contrast is then given by:

   .. math::

      \Delta_\mathrm{vir} = {\bar{\rho}_\mathrm{vir} \over \rho(r_\mathrm{vir})} \Delta_\mathrm{iso},

   where :math:`\bar{\rho}_\mathrm{vir}` is the mean density inside the virial radius and :math:`\rho(r_\mathrm{vir})` is the density at the virial radius. The ratio :math:`\bar{\rho}_\mathrm{vir} / \rho(r_\mathrm{vir})` is specified via the parameter ``[virialDensityContrastFoFDensityRatio]``. Its default value of :math:`4.688` is appropriate for an :term:`NFW` halo of concentration :math:`c=6.88` which is the concentration found by :cite:t:`prada_halo_2011` for halos with :math:`\sigma=1.686` which is the approximate critical overdensity for collapse).
   </description>
  </virialDensityContrast>
  !!]
  type, extends(virialDensityContrastClass) :: virialDensityContrastFriendsOfFriends
     !!{RST
     A dark matter halo virial density contrast class based on the friends-of-friends algorithm linking length.
     !!}
     private
     double precision :: linkingLength, densityRatio
   contains
     procedure :: densityContrast             => friendsOfFriendsDensityContrast
     procedure :: densityContrastRateOfChange => friendsOfFriendsDensityContrastRateOfChange
  end type virialDensityContrastFriendsOfFriends

  interface virialDensityContrastFriendsOfFriends
     !!{RST
     Constructors for the :galacticus-class:`virialDensityContrastFriendsOfFriends` dark matter halo virial density contrast class.
     !!}
     module procedure friendsOfFriendsConstructorParameters
     module procedure friendsOfFriendsConstructorInternal
  end interface virialDensityContrastFriendsOfFriends

contains

  function friendsOfFriendsConstructorParameters(parameters) result(self)
    !!{RST
    Default constructor for the ``friendsOfFriends`` dark matter halo virial density contrast class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (virialDensityContrastFriendsOfFriends)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    double precision                                                       :: linkingLength, densityRatio

    !![
    <inputParameter docformat="rst">
     <name>linkingLength</name>
     <source>parameters</source>
     <defaultValue>0.2d0</defaultValue>
     <description>
     The friends-of-friends linking length algorithm to use in computing virial density contrast.
     </description>
    </inputParameter>
    <inputParameter docformat="rst">
     <name>densityRatio</name>
     <source>parameters</source>
     <defaultValue>4.688d0</defaultValue>
     <defaultSource>
     Value appropriate for an :term:`NFW` profile with concentration :math:`c=6.88` which is the concentration found by :cite:t:`prada_halo_2011` for halos with :math:`\sigma=1.686` which is the approximate critical overdensity for collapse.
     </defaultSource>
     <description>
     The ratio of mean virial density to density at the virial radius to assume when setting virial density contrasts in the friends-of-friends model.
     </description>
    </inputParameter>
    !!]
    self=virialDensityContrastFriendsOfFriends(linkingLength,densityRatio)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function friendsOfFriendsConstructorParameters

  function friendsOfFriendsConstructorInternal(linkingLength,densityRatio) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`virialDensityContrastFriendsOfFriends` dark matter halo virial density contrast class.
    !!}
    implicit none
    type            (virialDensityContrastFriendsOfFriends), target        :: self
    double precision                                       , intent(in   ) :: linkingLength, densityRatio
    !![
    <constructorAssign variables="linkingLength, densityRatio"/>
    !!]

    return
  end function friendsOfFriendsConstructorInternal

  double precision function friendsOfFriendsDensityContrast(self,mass,time,expansionFactor,collapsing)
    !!{RST
    Return the virial density contrast at the given epoch, based on the friends-of-friends algorithm linking length.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (virialDensityContrastFriendsOfFriends), intent(inout)           :: self
    double precision                                       , intent(in   )           :: mass
    double precision                                       , intent(in   ), optional :: time                          , expansionFactor
    logical                                                , intent(in   ), optional :: collapsing
    double precision                                                                 :: boundingSurfaceDensityContrast
    !$GLC attributes unused :: mass, time, expansionFactor, collapsing

    boundingSurfaceDensityContrast =3.0d0/2.0d0/Pi/self%linkingLength**3
    friendsOfFriendsDensityContrast=self%densityRatio*boundingSurfaceDensityContrast
    return
  end function friendsOfFriendsDensityContrast

  double precision function friendsOfFriendsDensityContrastRateOfChange(self,mass,time,expansionFactor,collapsing)
    !!{RST
    Return the virial density contrast at the given epoch, based on the friends-of-friends algorithm linking length.
    !!}
    implicit none
    class           (virialDensityContrastFriendsOfFriends), intent(inout)           :: self
    double precision                                       , intent(in   )           :: mass
    double precision                                       , intent(in   ), optional :: time      , expansionFactor
    logical                                                , intent(in   ), optional :: collapsing
    !$GLC attributes unused :: self, mass, time, expansionFactor, collapsing

    friendsOfFriendsDensityContrastRateOfChange=0.0d0
    return
  end function friendsOfFriendsDensityContrastRateOfChange
