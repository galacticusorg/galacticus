!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% An abstract implementation of the orphan satellite distribution which assumes an isotropic distribution with randomly
  !% assigned positions.

  use Dark_Matter_Halo_Scales
  
  !# <satelliteOrphanDistribution name="satelliteOrphanDistributionTraceDarkMatter">
  !#  <description>An orphan satellite distribution which assumes an isotropic, random distribution with orphans tracing the radial distribution of dark matter.</description>
  !# </satelliteOrphanDistribution>
  type, extends(satelliteOrphanDistributionRandomIsotropic) :: satelliteOrphanDistributionTraceDarkMatter
     !% An orphan satellite distribution which assumes an isotropic, random distribution with orphans tracing the radial distribution of dark matter.
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_
   contains
     final     ::                                        traceDarkMatterDestructor
     procedure :: extent                              => traceDarkMatterExtent
     procedure :: inverseCumulativeMassFunctionRadial => traceDarkMatterInverseCMFRadial
  end type satelliteOrphanDistributionTraceDarkMatter

  interface satelliteOrphanDistributionTraceDarkMatter
     !% Constructors for the {\normalfont \ttfamily traceDarkMatter} orphan satellite distribution class.
     module procedure traceDarkMatterConstructorParameters
     module procedure traceDarkMatterConstructorInternal
  end interface satelliteOrphanDistributionTraceDarkMatter
  
contains

  function traceDarkMatterConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily traceDarkMatter} orphan satellite distribution class which takes a parameter
    !% list as input.
    use Input_Parameters2
    implicit none
    type(satelliteOrphanDistributionTraceDarkMatter)                :: self
    type(inputParameters                           ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <objectBuilder class="darkMatterHaloScale" name="self%darkMatterHaloScale_" source="parameters"/>
    return
  end function traceDarkMatterConstructorParameters

  function traceDarkMatterConstructorInternal(darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily traceDarkMatter} orphan satellite distribution class.
    implicit none
    type (satelliteOrphanDistributionTraceDarkMatter)                        :: self
    class(darkMatterHaloScaleClass                  ), intent(in   ), target :: darkMatterHaloScale_

    self%darkMatterHaloScale_ => darkMatterHaloScale_
    return
  end function traceDarkMatterConstructorInternal

  subroutine traceDarkMatterDestructor(self)
    !% Destructor for the {\normalfont \ttfamily traceDarkMatter} orphan satellite distribution class.
    implicit none
    type(satelliteOrphanDistributionTraceDarkMatter), intent(inout) :: self
    
    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    return
  end subroutine traceDarkMatterDestructor

  double precision function traceDarkMatterExtent(self,node)
    !% Return the maximum extent of the orphan satellite population for the {\normalfont \ttfamily traceDarkMatter} orphan
    !% satellite distribution class.
    implicit none
    class(satelliteOrphanDistributionTraceDarkMatter), intent(inout) :: self
    type (treeNode                                  ), intent(inout) :: node

    ! For this distribution the maximum extent is the virial radius.
    traceDarkMatterExtent=self%darkMatterHaloScale_%virialRadius(node)
    return
  end function traceDarkMatterExtent

  double precision function traceDarkMatterInverseCMFRadial(self,node,fraction)
    !% Return the radial coordinate within which the given {\normalfont \ttfamily fraction} of orphan satellites are found.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    class           (satelliteOrphanDistributionTraceDarkMatter), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    double precision                                            , intent(in   ) :: fraction
    type            (treeNode                                  ), pointer       :: nodeHost
    !GCC$ attributes unused :: self

    nodeHost                        => node%parent
    traceDarkMatterInverseCMFRadial =  Galactic_Structure_Radius_Enclosing_Mass(                                 &
         &                                                                      nodeHost                       , &
         &                                                                      fractionalMass=fraction        , &
         &                                                                      componentType =componentTypeAll, &
         &                                                                      massType      =massTypeDark    , &
         &                                                                      haloLoaded    =.true.            &
         &                                                                     )
    return
  end function traceDarkMatterInverseCMFRadial
