!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version. is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!% Contains a module which implements conversions of total halo mass between different definitions.

module Dark_Matter_Profile_Mass_Definitions
  !% Implements calculations of dark matter profile scale radii from concentrations.
  private
  public :: Dark_Matter_Profile_Mass_Definition

contains

  double precision function Dark_Matter_Profile_Mass_Definition(node,densityContrast,radius,velocity)
    !% Compute the mass of {\tt node} under the given density contrast definition.
    use Galacticus_Nodes
    use Cosmology_Parameters
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    implicit none    
    type            (treeNode                ), pointer , intent(inout) :: node
    double precision                                    , intent(in   ) :: densityContrast
    double precision                          , optional, intent(  out) :: radius              , velocity
    class           (nodeComponentBasic      ), pointer                 :: basic
    class           (cosmologyParametersClass), pointer                 :: cosmologyParameters_
    double precision                                                    :: darkMatterFraction  , radiusHalo

    ! Get required objects.
    cosmologyParameters_ => cosmologyParameters()
    ! Compute the dark matter fraction.
    darkMatterFraction=(1.0d0-cosmologyParameters_%omegaBaryon()/cosmologyParameters_%omegaMatter())
    ! Compute masses as necessary.
    basic                              => node%basic()
    radiusHalo                         =  Galactic_Structure_Radius_Enclosing_Density(node,densityContrast=densityContrast*darkMatterFraction,massType=massTypeDark,haloLoaded=.false.)
    Dark_Matter_Profile_Mass_Definition=  Galactic_Structure_Enclosed_Mass           (node,radiusHalo                                        ,massType=massTypeDark,haloLoaded=.false.) &
         &                               /darkMatterFraction
    if (present(radius  )) radius  =radiusHalo
    if (present(velocity)) velocity=sqrt(gravitationalConstantGalacticus*Dark_Matter_Profile_Mass_Definition/radiusHalo)
    return
  end function Dark_Matter_Profile_Mass_Definition

end module Dark_Matter_Profile_Mass_Definitions
