!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module which defines the abundances structure used for describing elemental abundances in \glc. (Currently a dumb
!% implementation.)

module Abundances_Structure
  !% Defines the abundances structure used for describing elemental abundances in \glc. (Currently a dumb implementation.)
  use Numerical_Constants_Astronomical
  private
  public :: abundancesStructure, Abundances_Names, Abundances_Property_Count, Abundances_Get_Metallicity
  

  type abundancesStructure
     !% The abundances structure used for describing elemental abundances in \glc. (Currently a dumb implementation.)
     double precision          :: metallicityValue
   contains
     ! Pack/unpack methods.
     !@ <objectMethods>
     !@   <object>abundancesStructure</object>
     !@   <objectMethod>
     !@     <method>pack</method>
     !@     <description>Packs abundance data into a 1-D array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>unpack(array)</method>
     !@     <description>Unpacks abundance data from a 1-D array.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                 :: pack                   => Abundances_Pack
     procedure                 :: unpack                 => Abundances_Unpack
     ! Metallicity methods.
     !@ <objectMethods>
     !@   <object>abundancesStructure</object>
     !@   <objectMethod>
     !@     <method>metallicity</method>
     !@     <description>Returns the metallicity.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>metallicitySet(metallicity)</method>
     !@     <description>Sets the metallicity to {\tt metallicity}.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                 :: metallicity            => Abundances_Get_Metallicity
     procedure                 :: metallicitySet         => Abundances_Set_Metallicity
     ! Hydrogen/helium methods.
     !@ <objectMethods>
     !@   <object>abundancesStructure</object>
     !@   <objectMethod>
     !@     <method>hydrogenNumberFraction</method>
     !@     <description>Returns the hydrogen fraction by number.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>hydrogenMassFraction</method>
     !@     <description>Returns the hydrogen fraction by mass.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>heliumMassFraction</method>
     !@     <description>Returns the helium fraction by mass.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                 :: hydrogenNumberFraction => Abundances_Hydrogen_Number_Fraction
     procedure                 :: hydrogenMassFraction   => Abundances_Hydrogen_Mass_Fraction
     procedure                 :: heliumMassFraction     => Abundances_Helium_Mass_Fraction
  end type abundancesStructure

  ! Count of the number of elements being tracked.
  integer,          parameter         :: elementsCount=0
  integer,          parameter         :: propertyCount=elementsCount+1

  ! Type of metallicity/abundance measure required.
  integer,          parameter, public :: linearByMass            =0
  integer,          parameter, public :: linearByNumber          =1
  integer,          parameter, public :: logarithmicByMassSolar  =2
  integer,          parameter, public :: logarithmicByNumberSolar=3
  
  ! Value used to indicate a zero metallicity on logarithmic scales.
  double precision, parameter, public :: logMetallicityZero      =-99.0d0

contains


  integer function Abundances_Property_Count()
    !% Return the number of properties required to track abundances. This is equal to the number of elements tracked, {\tt
    !% elementsCount}, plus one since we always track a total metallicity.
    implicit none

    Abundances_Property_Count=propertyCount
    return
  end function Abundances_Property_Count


  function Abundances_Names(index)
    !% Return a name for the specified entry in the abundances structure.
    use ISO_Varying_String
    implicit none
    type(varying_string)             :: Abundances_Names
    integer,              intent(in) :: index

    select case (index)
    case (1)
       Abundances_Names="Metals"
    end select
    return
  end function Abundances_Names


  subroutine Abundances_Pack(abundances,abundancesArray)
    !% Pack abundances from an array into an abundances structure.
    implicit none
    type(abundancesStructure), intent(out)              :: abundances
    double precision,          intent(in), dimension(:) :: abundancesArray

    abundances%metallicityValue=abundancesArray(1)
    return
  end subroutine Abundances_Pack

  subroutine Abundances_Unpack(abundances,abundancesArray)
    !% Pack abundances from an array into an abundances structure.
    implicit none
    double precision,          intent(out), dimension(:) :: abundancesArray(:)
    type(abundancesStructure), intent(in)                :: abundances

    abundancesArray(1)=abundances%metallicityValue
    return
  end subroutine Abundances_Unpack


  double precision function Abundances_Get_Metallicity(abundances,metallicityType)
    !% Return the metallicity of the {\tt abundances} structure.
    use Numerical_Constants_Astronomical
    use Galacticus_Error
    implicit none
    type(abundancesStructure), intent(in)           :: abundances
    integer,                   intent(in), optional :: metallicityType

    Abundances_Get_Metallicity=abundances%metallicityValue
    if (present(metallicityType)) then
       select case (metallicityType)
       case (linearByMass)
          ! Do nothing, this is what we compute by default.
       case (logarithmicByMassSolar)
          ! Convert to a logarithmic metallicity by mass relative to Solar.
          if (Abundances_Get_Metallicity > 0.0d0) then
             Abundances_Get_Metallicity=dlog10(Abundances_Get_Metallicity/metallicitySolar)
          else
             Abundances_Get_Metallicity=logMetallicityZero
          end if
       case default
          call Galacticus_Error_Report('Abundances_Get_Metallicity','metallicity type not supported')
       end select
    end if
    return
  end function Abundances_Get_Metallicity

  subroutine Abundances_Set_Metallicity(abundances,metallicity,metallicityType)
    !% Set the metallicity of the {\tt abundances} structure to {\tt metallicity}.
    use Galacticus_Error
    implicit none
    type(abundancesStructure), intent(inout)        :: abundances
    double precision,          intent(in)           :: metallicity
    integer,                   intent(in), optional :: metallicityType

    abundances%metallicityValue=metallicity
    if (present(metallicityType)) then
       select case (metallicityType)
          case (linearByMass)
             ! Do nothing, this is how we store metallicity.
          case (logarithmicByMassSolar)
             abundances%metallicityValue=(10.0d0**abundances%metallicityValue)*metallicitySolar
          case default
             call Galacticus_Error_Report('Abundances_Set_Metallicity','type not supported')
          end select
       end if
    return
  end subroutine Abundances_Set_Metallicity


  double precision function Abundances_Hydrogen_Mass_Fraction(abundances)
    !% Returns the mass fraction of hydrogen.
    implicit none
    type(abundancesStructure), intent(in) :: abundances

    Abundances_Hydrogen_Mass_Fraction=(abundances%metallicityValue/metallicitySolar)*(hydrogenByMassSolar-hydrogenByMassPrimordial)&
         &+hydrogenByMassPrimordial
    return
  end function Abundances_Hydrogen_Mass_Fraction

  double precision function Abundances_Helium_Mass_Fraction(abundances)
    !% Returns the mass fraction of helium.
    implicit none
    type(abundancesStructure), intent(in) :: abundances

    Abundances_Helium_Mass_Fraction=(abundances%metallicityValue/metallicitySolar)*(heliumByMassSolar-heliumByMassPrimordial)&
         &+heliumByMassPrimordial
    return
  end function Abundances_Helium_Mass_Fraction

  double precision function Abundances_Hydrogen_Number_Fraction(abundances)
    !% Returns the number fraction of hydrogen.
    implicit none
    type(abundancesStructure), intent(in) :: abundances
    double precision                      :: numberHydrogen,numberHelium

    numberHydrogen=Abundances_Hydrogen_Mass_Fraction(abundances)/atomicMassHydrogen
    numberHelium  =Abundances_Helium_Mass_Fraction(abundances)/atomicMassHelium
    Abundances_Hydrogen_Number_Fraction=numberHydrogen/(numberHydrogen+numberHelium)
    return
  end function Abundances_Hydrogen_Number_Fraction

end module Abundances_Structure
