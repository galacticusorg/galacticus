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
  An implementation of dark matter halo profile concentrations using a fixed value for concentration.
  !!}

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationFixed">
   <description>Dark matter halo concentrations are computed using a fixed value for concentration.</description>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationFixed
     !!{
     A dark matter halo profile concentration class implementing the algorithm of \cite{gao_redshift_2008}.
     !!}
     private
     class           (virialDensityContrastClass), pointer :: virialDensityContrast_ => null()
     class           (darkMatterProfileDMOClass ), pointer :: darkMatterProfileDMO_  => null()
     double precision                                      :: concentration_
   contains
     final     ::                                   fixedDestructor
     procedure :: concentration                  => fixedConcentration
     procedure :: densityContrastDefinition      => fixedDensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => fixedDarkMatterProfileDefinition
  end type darkMatterProfileConcentrationFixed

  interface darkMatterProfileConcentrationFixed
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationFixed} dark matter halo profile concentration class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface darkMatterProfileConcentrationFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileConcentrationFixed} dark matter halo profile concentration class which takes a parameter
    list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileConcentrationFixed)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (virialDensityContrastClass         ), pointer       :: virialDensityContrast_
    class           (darkMatterProfileDMOClass          ), pointer       :: darkMatterProfileDMO_
    double precision                                                     :: concentration

    !![
    <inputParameter>
      <name>concentration</name>
      <source>parameters</source>
      <description>The concentration.</description>
    </inputParameter>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"  source="parameters"/>
    !!]
    self=darkMatterProfileConcentrationFixed(concentration,virialDensityContrast_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="virialDensityContrast_"/>
    <objectDestructor name="darkMatterProfileDMO_" />
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(concentration_,virialDensityContrast_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileConcentrationFixed} dark matter halo profile concentration class.
    !!}
    implicit none
    type            (darkMatterProfileConcentrationFixed)                        :: self
    double precision                                     , intent(in   )         :: concentration_
    class           (virialDensityContrastClass         ), intent(in   ), target :: virialDensityContrast_
    class           (darkMatterProfileDMOClass          ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="concentration_, *virialDensityContrast_, *darkMatterProfileDMO_"/>
    !!]

    return
  end function fixedConstructorInternal

  subroutine fixedDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationFixed} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationFixed), intent(inout) :: self

    !![
    <objectDestructor name="self%virialDensityContrast_"/>
    <objectDestructor name="self%darkMatterProfileDMO_" />
    !!]
    return
  end subroutine fixedDestructor

  double precision function fixedConcentration(self,node)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(darkMatterProfileConcentrationFixed), intent(inout), target  :: self
    type (treeNode                           ), intent(inout), target  :: node
    !$GLC attributes unused :: node

    fixedConcentration=self%concentration_
    return
  end function fixedConcentration

  function fixedDensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of concentration.
    !!}
    implicit none
    class(virialDensityContrastClass         ), pointer       :: fixedDensityContrastDefinition
    class(darkMatterProfileConcentrationFixed), intent(inout) :: self

    fixedDensityContrastDefinition => self%virialDensityContrast_
    return
  end function fixedDensityContrastDefinition

  function fixedDarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of concentration.
    !!}
    implicit none
    class(darkMatterProfileDMOClass          ), pointer       :: fixedDarkMatterProfileDefinition
    class(darkMatterProfileConcentrationFixed), intent(inout) :: self

    fixedDarkMatterProfileDefinition => self%darkMatterProfileDMO_
    return
  end function fixedDarkMatterProfileDefinition

