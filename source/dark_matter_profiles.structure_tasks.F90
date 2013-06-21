!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements structure tasks related to the dark matter halo density profile.

module Dark_Matter_Profile_Structure_Tasks
  !% Implements structure tasks related to the dark matter halo density profile.
  use Galacticus_Nodes
  use Dark_Matter_Profiles
  private
  public :: Dark_Matter_Profile_Enclosed_Mass_Task,Dark_Matter_Profile_Density_Task, Dark_Matter_Profile_Density&
       &,Dark_Matter_Profile_Rotation_Curve_Task ,Dark_Matter_Profile_Potential_Task&
       &,Dark_Matter_Profile_Rotation_Curve_Gradient_Task

contains

  !# <enclosedMassTask>
  !#  <unitName>Dark_Matter_Profile_Enclosed_Mass_Task</unitName>
  !# </enclosedMassTask>
  double precision function Dark_Matter_Profile_Enclosed_Mass_Task(thisNode,radius,componentType,massType,weightBy,weightIndex,haloLoaded)
    !% Computes the mass within a given radius for a dark matter profile.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Cosmological_Parameters
    use Galactic_Structure_Initial_Radii
    implicit none
    type            (treeNode          ), intent(inout), pointer  :: thisNode
    integer                             , intent(in   )           :: componentType     , massType     , weightBy, &
         &                                                           weightIndex
    double precision                    , intent(in   )           :: radius
    logical                             , intent(in   ), optional :: haloLoaded
    class           (nodeComponentBasic)               , pointer  :: thisBasicComponent
    double precision                                              :: darkMatterFraction, radiusInitial
    logical                                                       :: haloLoadedActual

    Dark_Matter_Profile_Enclosed_Mass_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeDark         )) return
    if (.not.(weightBy      == weightByMass                                                )) return

    ! Determine the dark matter fraction.
    darkMatterFraction=(Omega_Matter()-Omega_b())/Omega_Matter()
    ! Test radius.
    if (radius >= radiusLarge) then
       ! Return the total mass of the halo in this case.
       thisBasicComponent => thisNode%basic()
       Dark_Matter_Profile_Enclosed_Mass_Task=thisBasicComponent%mass()
    else if (radius <= 0.0d0) then
       ! Zero radius. Return zero mass.
       Dark_Matter_Profile_Enclosed_Mass_Task=0.0d0
    else
       ! Determine if we need to account for halo loading (a.k.a. adiabatic contraction, a.k.a. baryonic pinching).
       haloLoadedActual=.true.
       if (present(haloLoaded)) haloLoadedActual=haloLoaded
       if (haloLoadedActual) then
          ! Halo loading is to be accounted for - get the initial radius in the dark matter halo.
          radiusInitial=Galactic_Structure_Radius_Initial(thisNode,radius)
       else
          ! Halo loading is not to be accounted for. The radius to use is simply the radius we were given.
          radiusInitial=radius
       end if
       ! Return the mass within the initial radius.
       Dark_Matter_Profile_Enclosed_Mass_Task=Dark_Matter_Profile_Enclosed_Mass(thisNode,radiusInitial)
    end if
    ! Scale to account for just the dark component.
    Dark_Matter_Profile_Enclosed_Mass_Task=Dark_Matter_Profile_Enclosed_Mass_Task*darkMatterFraction
    return
  end function Dark_Matter_Profile_Enclosed_Mass_Task

  !# <rotationCurveTask>
  !#  <unitName>Dark_Matter_Profile_Rotation_Curve_Task</unitName>
  !# </rotationCurveTask>
  double precision function Dark_Matter_Profile_Rotation_Curve_Task(thisNode,radius,componentType,massType,haloLoaded)
    !% Computes the rotation curve at a given radius for a dark matter profile.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    implicit none
    type            (treeNode), intent(inout), pointer  :: thisNode
    integer                   , intent(in   )           :: componentType, massType
    double precision          , intent(in   )           :: radius
    logical                   , intent(in   ), optional :: haloLoaded
    double precision                                    :: componentMass

    ! Set to zero by default.
    Dark_Matter_Profile_Rotation_Curve_Task=0.0d0

    ! Compute if a spheroid is present.
    if (radius > 0.0d0) then
       componentMass=Dark_Matter_Profile_Enclosed_Mass_Task(thisNode,radius,componentType,massType,weightByMass,weightIndexNull&
            &,haloLoaded)
       if (componentMass > 0.0d0) Dark_Matter_Profile_Rotation_Curve_Task=sqrt(gravitationalConstantGalacticus*componentMass)&
            &/sqrt(radius)
    end if
    return
  end function Dark_Matter_Profile_Rotation_Curve_Task

  !# <densityTask>
  !#  <unitName>Dark_Matter_Profile_Density_Task</unitName>
  !# </densityTask>
  double precision function Dark_Matter_Profile_Density_Task(thisNode,positionSpherical,componentType,massType,weightBy,weightIndex,haloLoaded)
    !% Computes the density at a given position for a dark matter profile.
    use Galactic_Structure_Options
    use Numerical_Constants_Math
    use Galacticus_Error
    use Galactic_Structure_Initial_Radii
    use Cosmological_Parameters
    implicit none
    type            (treeNode), intent(inout), pointer  :: thisNode
    integer                   , intent(in   )           :: componentType        , massType     , weightBy      , &
         &                                                 weightIndex
    double precision          , intent(in   )           :: positionSpherical (3)
    logical                   , intent(in   ), optional :: haloLoaded
    logical                                             :: haloLoadedActual
    double precision                                    :: darkMatterFraction   , radiusInitial, radiusJacobian

    ! Return zero if the component and mass type is not matched.
    Dark_Matter_Profile_Density_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeDark         )) return
    if (.not.(weightBy      == weightByMass                                                )) return
    ! Determine the dark matter fraction.
    darkMatterFraction=(Omega_Matter()-Omega_b())/Omega_Matter()
    ! Determine if we need to account for halo loading (a.k.a. adiabatic contraction, a.k.a. baryonic pinching).
    haloLoadedActual=.true.
    if (present(haloLoaded)) haloLoadedActual=haloLoaded
    if (haloLoadedActual) then
       ! Halo loading is to be accounted for - get the initial radius in the dark matter halo, and the Jacobian.
       radiusInitial = Galactic_Structure_Radius_Initial           (thisNode,positionSpherical(1))
       radiusJacobian= Galactic_Structure_Radius_Initial_Derivative(thisNode,positionSpherical(1)) &
            &         *(radiusInitial/positionSpherical(1))**2
    else
       ! Halo loading is not to be accounted for. The radius to use is simply the radius we were given.
       radiusInitial =positionSpherical(1)
       radiusJacobian=1.0d0
    end if
    ! Return the mass within the initial radius.
    Dark_Matter_Profile_Density_Task=Dark_Matter_Profile_Density(thisNode,radiusInitial)
    ! Account for the Jacobian.
    Dark_Matter_Profile_Density_Task=Dark_Matter_Profile_Density_Task*radiusJacobian
    ! Scale to account for just the dark component.
    Dark_Matter_Profile_Density_Task=Dark_Matter_Profile_Density_Task*darkMatterFraction
    return
  end function Dark_Matter_Profile_Density_Task

  !# <rotationCurveGradientTask>
  !#  <unitName>Dark_Matter_Profile_Rotation_Curve_Gradient_Task</unitName>
  !# </rotationCurveGradientTask>
  double precision function Dark_Matter_Profile_Rotation_Curve_Gradient_Task(thisNode,radius,componentType,massType,haloLoaded)
    !% Computes the rotation curve gradient for the dark matter.
    use Galacticus_Nodes
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Math
    use Galacticus_Error
    implicit none
    type            (treeNode), intent(inout), pointer  :: thisNode
    integer                   , intent(in   )           :: componentType   , massType
    double precision          , intent(in   )           :: radius
    logical                   , intent(in   ), optional :: haloLoaded
    double precision                                    :: componentDensity, componentMass, positionSpherical(3)

    ! Set to zero by default.
    Dark_Matter_Profile_Rotation_Curve_Gradient_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(     massType == massTypeAll      .or.      massType == massTypeDark         )) return
    if (radius <= 0.0d0) return

    positionSpherical=[radius,0.0d0,0.0d0]
    componentMass   =Dark_Matter_Profile_Enclosed_Mass_Task(thisNode,radius           ,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
    componentDensity=Dark_Matter_Profile_Density_Task      (thisNode,positionSpherical,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
    if (componentMass ==0.0d0 .or. componentDensity == 0.0d0) return
    Dark_Matter_Profile_Rotation_Curve_Gradient_Task=           &
         &                   gravitationalConstantGalacticus    &
         &                  *(-componentMass/radius**2          &
         &                    +4.0d0*Pi*radius*componentDensity &
         &                   )
    return
  end function Dark_Matter_Profile_Rotation_Curve_Gradient_Task

  !# <potentialTask>
  !#  <unitName>Dark_Matter_Profile_Potential_Task</unitName>
  !# </potentialTask>
  double precision function Dark_Matter_Profile_Potential_Task(thisNode,radius,componentType,massType,haloLoaded)
    !% Return the potential due to dark matter.
    use Galactic_Structure_Options
    use Galacticus_Error
    implicit none
    type            (treeNode), intent(inout), pointer  :: thisNode
    integer                   , intent(in   )           :: componentType, massType
    double precision          , intent(in   )           :: radius
    logical                   , intent(in   ), optional :: haloLoaded

    Dark_Matter_Profile_Potential_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeDark         )) return

    if (present(haloLoaded)) then
       if (haloLoaded) call Galacticus_Error_Report('Dark_Matter_Profile_Potential_Task','dark matter potential not available for baryon loaded halos')
    end if

    Dark_Matter_Profile_Potential_Task=Dark_Matter_Profile_Potential(thisNode,radius)
    return
  end function Dark_Matter_Profile_Potential_Task

end module Dark_Matter_Profile_Structure_Tasks
