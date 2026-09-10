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

!+    Contributions to this file made by:  Yu Zhao.

  !!{RST
    Implementation of a satellite tidal radius class which limits the tidal radius to be no smaller than the radius of the solitonic core.
  !!}

  !![
  <satelliteTidalStrippingRadius name="satelliteTidalStrippingRadiusSolitonLimited" docformat="rst">
   <description>
   A satellite tidal radius class which limits the tidal radius returned by a wrapped
   :galacticus-class:`satelliteTidalStrippingRadiusClass` to be no smaller than the radius of the solitonic core, so that
   stripping of the outer halo can not remove material from within the core itself. Intended for use with the
   :galacticus-class:`darkMatterProfileDMOSolitonNFW` and :galacticus-class:`darkMatterProfileDMOSolitonNFWHeated` dark matter
   profiles, which create the ``solitonRadiusSoliton`` meta-property that this class reads. In any model where that meta-property
   is not created---i.e. one containing no fuzzy dark matter---and for halos in which no soliton formed, or which have already
   been stripped down to their core, the wrapped tidal radius is returned unchanged.
   </description>
  </satelliteTidalStrippingRadius>
  !!]
  type, extends(satelliteTidalStrippingRadiusClass) :: satelliteTidalStrippingRadiusSolitonLimited
     !!{RST
     Implementation of a satellite tidal radius class which limits the tidal radius to be no smaller than the solitonic core radius.
     !!}
     private
     class  (satelliteTidalStrippingRadiusClass), pointer :: satelliteTidalStrippingRadius_ => null()
     integer                                              :: radiusSolitonID
   contains
     final     ::           solitonLimitedDestructor
     procedure :: radius => solitonLimitedRadius
  end type satelliteTidalStrippingRadiusSolitonLimited

  interface satelliteTidalStrippingRadiusSolitonLimited
     !!{RST
     Constructors for the :galacticus-class:`satelliteTidalStrippingRadiusSolitonLimited` satellite tidal stripping class.
     !!}
     module procedure solitonLimitedConstructorParameters
     module procedure solitonLimitedConstructorInternal
  end interface satelliteTidalStrippingRadiusSolitonLimited

contains

  function solitonLimitedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`satelliteTidalStrippingRadiusSolitonLimited` satellite tidal stripping class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteTidalStrippingRadiusSolitonLimited)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (satelliteTidalStrippingRadiusClass  ), pointer       :: satelliteTidalStrippingRadius_

    !![
    <objectBuilder class="satelliteTidalStrippingRadius" name="satelliteTidalStrippingRadius_" source="parameters"/>
    !!]
    self=satelliteTidalStrippingRadiusSolitonLimited(satelliteTidalStrippingRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteTidalStrippingRadius_"/>
    !!]
    return
  end function solitonLimitedConstructorParameters

  function solitonLimitedConstructorInternal(satelliteTidalStrippingRadius_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`satelliteTidalStrippingRadiusSolitonLimited` satellite tidal stripping class.
    !!}
    implicit none
    type            (satelliteTidalStrippingRadiusSolitonLimited)                     :: self
    class           (satelliteTidalStrippingRadiusClass  ), intent(in), target :: satelliteTidalStrippingRadius_
    !![
    <constructorAssign variables="*satelliteTidalStrippingRadius_"/>
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusSoliton"  id="self%radiusSolitonID"  isEvolvable="no"  isCreator="no"/>
    !!]
    return
  end function solitonLimitedConstructorInternal

  subroutine solitonLimitedDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`satelliteTidalStrippingRadiusSolitonLimited` satellite tidal stripping class.
    !!}
    implicit none
    type            (satelliteTidalStrippingRadiusSolitonLimited), intent(inout)      :: self

    !![
    <objectDestructor name="self%satelliteTidalStrippingRadius_"/>
    !!]
    return
  end subroutine solitonLimitedDestructor

  double precision function solitonLimitedRadius(self,node)
    !!{RST
    Return the tidal radius, limited to be no smaller than the radius of the solitonic core.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (satelliteTidalStrippingRadiusSolitonLimited), intent(inout), target :: self
    type            (treeNode                            ), intent(inout), target :: node
    class           (nodeComponentDarkMatterProfile      ), pointer               :: darkMatterProfile
    double precision                                                              :: radiusTidal      , radiusSoliton

    ! Compute the tidal radius using the wrapped tidal-radius model.
    radiusTidal=self%satelliteTidalStrippingRadius_%radius(node)
    ! Find the radius of the solitonic core. The meta-property always has a valid ID, but has storage only if some class creates
    ! it - which is not the case in a model containing no fuzzy dark matter. A non-positive radius indicates a halo in which no
    ! soliton formed, or one which has already been stripped down to its core so that no NFW envelope remains to be stripped
    ! (the soliton profiles store -1 in both cases). In all of these cases the tidal radius is returned unchanged.
    radiusSoliton     =  0.0d0
    darkMatterProfile => node%darkMatterProfile()
    if (darkMatterProfile%floatRank0MetaPropertyIsCreated(self%radiusSolitonID)) &
         & radiusSoliton=darkMatterProfile%floatRank0MetaPropertyGet(self%radiusSolitonID)
    ! Prevent the outer-halo stripping radius from entering the solitonic core.
    if (radiusSoliton > 0.0d0) then
       solitonLimitedRadius=max(radiusTidal,radiusSoliton)
    else
       solitonLimitedRadius=    radiusTidal
    end if
    return
  end function solitonLimitedRadius
