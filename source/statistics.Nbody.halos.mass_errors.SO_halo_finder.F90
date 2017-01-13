!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements an N-body dark matter halo mass error class which
!% implements a model for errors in spherical overdensity halo finders.

  use Dark_Matter_Halo_Scales
  use Dark_Matter_Profiles
  
  !# <nbodyHaloMassError name="nbodyHaloMassErrorSOHaloFinder">
  !#  <description>An N-body dark matter halo mass error class which implements a model for errors in spherical overdensity halo finders.</description>
  !# </nbodyHaloMassError>
  type, extends(nbodyHaloMassErrorClass) :: nbodyHaloMassErrorSOHaloFinder
     !% An N-body halo mass error class which implements a model for errors in spherical
     !% overdensity halo finders.
     private
     double precision                                    :: massParticle
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_
     class           (darkMatterProfileClass  ), pointer :: darkMatterProfile_
   contains
     procedure :: errorFractional => soHaloFinderErrorFractional
  end type nbodyHaloMassErrorSOHaloFinder

  interface nbodyHaloMassErrorSOHaloFinder
     !% Constructors for the {\normalfont \ttfamily soHaloFinder} N-body halo mass error class.
     module procedure nbodyHaloMassErrorSOHaloFinderParameters
     module procedure nbodyHaloMassErrorSOHaloFinderInternal
  end interface nbodyHaloMassErrorSOHaloFinder

contains

  function nbodyHaloMassErrorSOHaloFinderParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily soHaloFinder} N-body halo mass error class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(nbodyHaloMassErrorSOHaloFinder)                :: nbodyHaloMassErrorSOHaloFinderParameters
    type(inputParameters               ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)
    !# <inputParameter>
    !#   <name>massParticle</name>
    !#   <source>parameters</source>
    !#   <variable>nbodyHaloMassErrorSOHaloFinderParameters%massParticle</variable>
    !#   <description>Mass of particle in the simulation to which the spherical overdensity algorithm was applied.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    nbodyHaloMassErrorSOHaloFinderParameters%darkMatterHaloScale_ => darkMatterHaloScale()
    nbodyHaloMassErrorSOHaloFinderParameters%darkMatterProfile_   => darkMatterProfile  ()
    return
  end function nbodyHaloMassErrorSOHaloFinderParameters

  function nbodyHaloMassErrorSOHaloFinderInternal(darkMatterHaloScale_,darkMatterProfile_,massParticle)
    !% Internal constructor for the {\normalfont \ttfamily soHaloFinder} N-body halo mass error class.
    implicit none
    type            (nbodyHaloMassErrorSOHaloFinder)                        :: nbodyHaloMassErrorSOHaloFinderInternal
    class           (darkMatterHaloScaleClass      ), target, intent(in   ) :: darkMatterHaloScale_
    class           (darkMatterProfileClass        ), target, intent(in   ) :: darkMatterProfile_
    double precision                                        , intent(in   ) :: massParticle
 
    nbodyHaloMassErrorSOHaloFinderInternal%darkMatterHaloScale_   => darkMatterHaloScale_
    nbodyHaloMassErrorSOHaloFinderInternal%darkMatterProfile_     => darkMatterProfile_
    nbodyHaloMassErrorSOHaloFinderInternal%massParticle           =  massParticle
    return
  end function nbodyHaloMassErrorSOHaloFinderInternal

  double precision function soHaloFinderErrorFractional(self,node)
    !% Return the fractional error on the mass of an N-body halo in the power-law error model.
    use Numerical_Constants_Math
    implicit none
    class           (nbodyHaloMassErrorSOHaloFinder), intent(inout)            :: self
    type            (treeNode                      ), intent(inout), pointer   :: node
    class           (nodeComponentBasic            )               , pointer   :: basic
    double precision                                                           :: radiusHalo                   , densityOuterRadius, &
         &                                                                        densityRatioInternalToSurface, particleCount     , &
         &                                                                        errorFractionalFixedSphere

    ! Get the basic component of the node.
    basic                         =>  node                     %basic       (               )
    ! Determine number of particles in the halo.
    particleCount                 =  +basic%mass        () &
         &                           /self %massParticle
    ! Fractional error in mass wihtin fixed sphere, assuming Poisson statistics (which should be valid for a halo which contains a
    ! fraction of all particles in the simulation that is much less than unity).
    errorFractionalFixedSphere    =  +1.0d0               &
         &                           /sqrt(particleCount)
    ! Get the outer radius of the halo.
    radiusHalo                    =  +self%darkMatterHaloScale_%virialRadius(node           )
    ! Get the density at the edge of the halo.
    densityOuterRadius            =  +self%darkMatterProfile_  %density     (node,radiusHalo)
    ! Find the ratio of the mean interior density in the halo to the density at the halo outer radius.
    densityRatioInternalToSurface =  +3.0d0                 &
         &                           *basic%mass()          &
         &                           /4.0d0                 &
         &                           /Pi                    &
         &                           /radiusHalo        **3 &
         &                           /densityOuterRadius
    ! Compute the total error accounting for the change in size of the spherical region defining the halo.
    soHaloFinderErrorFractional   =  +errorFractionalFixedSphere        &
         &                           *(                                 &
         &                             +1.0d0                           &
         &                             +1.0d0                           &
         &                             /(                               &
         &                               +densityRatioInternalToSurface &
         &                               -1.0d0                         &
         &                              )                               &
         &                            )
    return
  end function soHaloFinderErrorFractional
  
