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

!!{
Contains a module which implements conversions of total halo mass between different definitions.
!!}

module Dark_Matter_Profile_Mass_Definitions
  !!{
  Implements calculations of dark matter profile scale radii from concentrations.
  !!}
  private
  public :: Dark_Matter_Profile_Mass_Definition

contains

  function Dark_Matter_Profile_Mass_Definition(node,densityContrast,radius,velocity,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterProfileDMO_,useLastIsolatedTime) result(massHalo)
    !!{
    Compute the mass of {\normalfont \ttfamily node} under the given density contrast definition.
    !!}
    use :: Cosmology_Functions             , only : cosmologyFunctionsClass
    use :: Cosmology_Parameters            , only : cosmologyParametersClass
    use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMOClass
    use :: Galacticus_Nodes                , only : nodeComponentBasic            , treeNode
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Math_Exponentiation             , only : cubeRoot
    use :: Numerical_Comparison            , only : Values_Agree
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Virial_Density_Contrast         , only : virialDensityContrastClass
    use :: Galactic_Structure_Options      , only : componentTypeDarkMatterOnly   , massTypeDark
    implicit none
    double precision                                                      :: massHalo
    type            (treeNode                  )          , intent(inout) :: node
    double precision                                      , intent(in   ) :: densityContrast
    double precision                            , optional, intent(  out) :: radius                , velocity
    logical                                     , optional, intent(in   ) :: useLastIsolatedTime
    class           (cosmologyParametersClass  )          , intent(inout) :: cosmologyParameters_
    class           (cosmologyFunctionsClass   )          , intent(inout) :: cosmologyFunctions_
    class           (virialDensityContrastClass)          , intent(inout) :: virialDensityContrast_
    class           (darkMatterProfileDMOClass ), optional, intent(inout) :: darkMatterProfileDMO_
    class           (massDistributionClass     ), pointer                 :: massDistribution_
    class           (nodeComponentBasic        ), pointer                 :: basic
    double precision                                                      :: radiusHalo            , density , &
         &                                                                   time
    !![
    <optionalArgument name="useLastIsolatedTime" defaultsTo=".false." />
    !!]

    ! Compute the density from the density contrast.
    basic   =>  node                %basic           (    )
    if (useLastIsolatedTime_) then
       time =  basic                %timeLastIsolated(    )
    else
       time =  basic                %time            (    )
    end if
    density =  +                     densityContrast           &
         &     *cosmologyParameters_%omegaMatter     (    )    &
         &     *cosmologyParameters_%densityCritical (    )    &
         &     /cosmologyFunctions_ %expansionFactor (time)**3
    ! Check if requested density contrast matches the virial density contrast.
    if     (                                                                                                    &
         &   Values_Agree(                                                                                      &
         &                       virialDensityContrast_%densityContrast(basic%mass(),basic%timeLastIsolated()), &
         &                                              densityContrast                                       , &
         &                relTol=1.0d-4                                                                         &
         &               )                                                                                      &
         &  .and.                                                                                               &
         &   time == basic%timeLastIsolated()                                                                   &
         & ) then
       ! Requested density contrast is just the virial density contrast.
       massHalo=basic%mass()
       if (present(radius).or.present(velocity)) &
            & radiusHalo=+cubeRoot(              &
            &                      +3.0d0        &
            &                      /4.0d0        &
            &                      /Pi           &
            &                      *massHalo     &
            &                      /density      &
            &                     )
    else
       ! Mismatched density contrast definitions - compute the mass directly.
       ! Get the radius in the halo enclosing this density.
       if (present(darkMatterProfileDMO_)) then
          massDistribution_ => darkMatterProfileDMO_%get                   (node                                    )
       else
          massDistribution_ => node                 %massDistribution      (componentTypeDarkMatterOnly,massTypeDark)
       end if
       radiusHalo           =  massDistribution_    %radiusEnclosingDensity(density                                 )
       !![
       <objectDestructor name="massDistribution_"/>
       !!]       
       ! Find the mass within that radius - this is computable directly from the mean density and the radius enclosing that mean
       ! density.
       massHalo  =+4.0d0         &
            &     *Pi            &
            &     *density       &
            &     *radiusHalo**3 &
            &     /3.0d0
    end if
    ! If necessary, return the radius and circular velocity also.
    if (present(radius  )) radius=radiusHalo
    if (present(velocity)) then
       if (radiusHalo > 0.0d0) then
          velocity=sqrt(gravitationalConstant_internal*massHalo/radiusHalo)
       else
          velocity=0.0d0
       end if
    end if
    return
  end function Dark_Matter_Profile_Mass_Definition

end module Dark_Matter_Profile_Mass_Definitions
