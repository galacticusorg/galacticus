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

!% Contains a module which implements useful cosmological functions.

module Cosmology_Functions
  !% Implements useful cosmological functions.
  use, intrinsic :: ISO_C_Binding
  use ISO_Varying_String
  implicit none
  private
  public :: Cosmology_Age, Expansion_Factor, Hubble_Parameter, Early_Time_Density_Scaling, Expansion_Factor_Is_Valid,&
       & Cosmic_Time_Is_Valid, Omega_Matter_Total, Omega_Dark_Energy, Expansion_Rate, Epoch_of_Matter_Dark_Energy_Equality,&
       & Epoch_of_Matter_Domination, Expansion_Factor_from_Redshift, Redshift_from_Expansion_Factor, CMB_Temperature,&
       & Comoving_Distance, Time_From_Comoving_Distance, Comoving_Distance_Conversion, Comoving_Volume_Element_Redshift,&
       & Comoving_Volume_Element_Time, Cosmology_Dark_Energy_Equation_Of_State, Cosmology_Dark_Energy_Exponent

  ! Flag to indicate if this module has been initialized.
  logical                                                                 :: cosmologyInitialized                       =.false.

  ! Name of cosmology functions method used.
  type     (varying_string                                     )          :: cosmologyMethod

  ! Pointer to the functions that actually do the calculations.
  procedure(Cosmology_Logical_Function_Double_Template         ), pointer :: Expansion_Factor_Is_Valid_Get              =>null()
  procedure(Cosmology_Logical_Function_Double_Template         ), pointer :: Cosmic_Time_Is_Valid_Get                   =>null()
  procedure(Cosmology_Double_Function_dCollapse_Template       ), pointer :: Cosmology_Age_Get                          =>null()
  procedure(Cosmology_Double_Function_Double_Template          ), pointer :: Expansion_Factor_Get                       =>null()
  procedure(Cosmology_Double_Function_ddCollapse_Template      ), pointer :: Hubble_Parameter_Get                       =>null()
  procedure(Cosmology_Density_Scaling_Template                 ), pointer :: Early_Time_Density_Scaling_Get             =>null()
  procedure(Cosmology_Double_Function_ddCollapse_Template      ), pointer :: Omega_Matter_Total_Get                     =>null()
  procedure(Cosmology_Double_Function_ddCollapse_Template      ), pointer :: Omega_Dark_Energy_Get                      =>null()
  procedure(Cosmology_Double_Function_ddCollapse_Template      ), pointer :: CMB_Temperature_Get                        =>null()
  procedure(Cosmology_Double_Function_Double_Template          ), pointer :: Expansion_Rate_Get                         =>null()
  procedure(Cosmology_Double_Function_Optional_Integer_Template), pointer :: Epoch_of_Matter_Dark_Energy_Equality_Get   =>null()
  procedure(Cosmology_Double_Function_Optional_Integer_Template), pointer :: Epoch_of_Matter_Curvature_Equality_Get     =>null()
  procedure(Cosmology_Double_Function_Double_Template          ), pointer :: Epoch_of_Matter_Domination_Get             =>null()
  procedure(Cosmology_Double_Function_Double_Template          ), pointer :: Comoving_Distance_Get                      =>null()
  procedure(Cosmology_Double_Function_Double_Template          ), pointer :: Time_From_Comoving_Distance_Get            =>null()
  procedure(Cosmology_Comoving_Distance_Conversion_Template    ), pointer :: Comoving_Distance_Conversion_Get           =>null()
  procedure(Cosmology_Dark_Energy_Equation_Of_State            ), pointer :: Cosmology_Dark_Energy_Equation_Of_State_Get=>null()
  procedure(Cosmology_Dark_Energy_Exponent                     ), pointer :: Cosmology_Dark_Energy_Exponent_Get         =>null()

  abstract interface
     double precision function Cosmology_Comoving_Distance_Conversion_Template(output,distanceModulus,redshift)
       integer         , intent(in   )           :: output
       double precision, intent(in   ), optional :: distanceModulus, redshift
     end function Cosmology_Comoving_Distance_Conversion_Template
  end interface

  abstract interface
     double precision function Cosmology_Double_Function_Optional_Integer_Template(inputParameter)
       integer, intent(in   ), optional :: inputParameter
     end function Cosmology_Double_Function_Optional_Integer_Template
  end interface

  abstract interface
     double precision function Cosmology_Double_Function_Double_Template(inputParameter)
       double precision, intent(in   ) :: inputParameter
     end function Cosmology_Double_Function_Double_Template
  end interface

  abstract interface
     double precision function Cosmology_Double_Function_dCollapse_Template(inputParameter,collapsingPhase)
       double precision, intent(in   )           :: inputParameter
       logical         , intent(in   ), optional :: collapsingPhase
     end function Cosmology_Double_Function_dCollapse_Template
  end interface

  abstract interface
     double precision function Cosmology_Double_Function_ddCollapse_Template(inputParameter1,inputParameter2,collapsingPhase)
       double precision, intent(in   ), optional :: inputParameter1, inputParameter2
       logical         , intent(in   ), optional :: collapsingPhase
     end function Cosmology_Double_Function_ddCollapse_Template
  end interface

  abstract interface
     logical function Cosmology_Logical_Function_Double_Template(inputParameter)
       double precision, intent(in   ) :: inputParameter
     end function Cosmology_Logical_Function_Double_Template
  end interface

  abstract interface
     subroutine Cosmology_Density_Scaling_Template(dominateFactor,densityPower,aDominant,Omega_Dominant)
       double precision, intent(in   )           :: dominateFactor
       double precision, intent(  out)           :: aDominant     , densityPower
       double precision, intent(  out), optional :: Omega_Dominant
     end subroutine Cosmology_Density_Scaling_Template
  end interface

contains

  subroutine Cosmology_Functions_Initialize
    !% Initialize the cosmology functions module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="cosmologyMethod" type="moduleUse">
    include 'cosmology_functions.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.cosmologyInitialized) then
       !$omp critical(Cosmology_Functions_Initialization)
       if (.not.cosmologyInitialized) then
          ! Get the disk star formation feedback method parameter.
          !@ <inputParameter>
          !@   <name>cosmologyMethod</name>
          !@   <defaultValue>matter-lambda</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for cosmology calculations.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('cosmologyMethod',cosmologyMethod,defaultValue='matter-lambda')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="cosmologyMethod" type="functionCall" functionType="void">
          !#  <functionArgs>cosmologyMethod,Expansion_Factor_Is_Valid_Get,Cosmic_Time_Is_Valid_Get,Cosmology_Age_Get,Expansion_Factor_Get,Hubble_Parameter_Get,Early_Time_Density_Scaling_Get,Omega_Matter_Total_Get,Omega_Dark_Energy_Get,Expansion_Rate_Get,Epoch_of_Matter_Dark_Energy_Equality_Get,Epoch_of_Matter_Domination_Get,Epoch_of_Matter_Curvature_Equality_Get,CMB_Temperature_Get,Comoving_Distance_Get,Time_From_Comoving_Distance_Get,Comoving_Distance_Conversion_Get,Cosmology_Dark_Energy_Equation_Of_State_Get,Cosmology_Dark_Energy_Exponent_Get</functionArgs>
          include 'cosmology_functions.inc'
          !# </include>
          if     (                                                                   &
               &  .not.(                                                             &
               &             associated(Expansion_Factor_Is_Valid_Get              ) &
               &        .and.associated(Cosmic_Time_Is_Valid_Get                   ) &
               &        .and.associated(Cosmology_Age_Get                          ) &
               &        .and.associated(Expansion_Factor_Get                       ) &
               &        .and.associated(Hubble_Parameter_Get                       ) &
               &        .and.associated(Early_Time_Density_Scaling_Get             ) &
               &        .and.associated(Omega_Matter_Total_Get                     ) &
               &        .and.associated(Omega_Dark_Energy_Get                      ) &
               &        .and.associated(Expansion_Rate_Get                         ) &
               &        .and.associated(Epoch_of_Matter_Dark_Energy_Equality_Get   ) &
               &        .and.associated(Epoch_of_Matter_Domination_Get             ) &
               &        .and.associated(Epoch_of_Matter_Curvature_Equality_Get     ) &
               &        .and.associated(CMB_Temperature_Get                        ) &
               &        .and.associated(Comoving_Distance_Get                      ) &
               &        .and.associated(Time_From_Comoving_Distance_Get            ) &
               &        .and.associated(Comoving_Distance_Conversion_Get           ) &
               &        .and.associated(Cosmology_Dark_Energy_Equation_Of_State_Get) &
               &        .and.associated(Cosmology_Dark_Energy_Exponent_Get         ) &
               &       )                                                             &
               & )                                                                   &
               & call Galacticus_Error_Report('Cosmology_Functions','method '//char(cosmologyMethod)//' is unrecognized')
          cosmologyInitialized=.true.
       end if
       !$omp end critical(Cosmology_Functions_Initialization)
    end if
    return
  end subroutine Cosmology_Functions_Initialize

  logical function Expansion_Factor_Is_Valid(aExpansion)
    !% Returns true if the given expansion factor is valid one for this cosmology.
    implicit none
    double precision, intent(in   ) :: aExpansion

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    Expansion_Factor_Is_Valid=Expansion_Factor_Is_Valid_Get(aExpansion)

    return
  end function Expansion_Factor_Is_Valid

  logical function Cosmic_Time_Is_Valid(time)
    !% Returns true if the given cosmic time is valid one for this cosmology.
    implicit none
    double precision, intent(in   ) :: time

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    Cosmic_Time_Is_Valid=Cosmic_Time_Is_Valid_Get(time)

    return
  end function Cosmic_Time_Is_Valid

  double precision function Cosmology_Age(aExpansion,collapsingPhase)
    !% Return the age of the universe (in Gyr) at given expansion factor.
    implicit none
    double precision, intent(in   )           :: aExpansion
    logical         , intent(in   ), optional :: collapsingPhase

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    Cosmology_Age=Cosmology_Age_Get(aExpansion,collapsingPhase)

    return
  end function Cosmology_Age

  double precision function Expansion_Factor(tCosmological)
    !% Returns the expansion factor at cosmological time {\tt tCosmological}.
    implicit none
    double precision, intent(in   ) :: tCosmological

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    Expansion_Factor=Expansion_Factor_Get(tCosmological)

    return
  end function Expansion_Factor

  function Expansion_Factor_C(tCosmological) bind(c,name="Expansion_Factor")
    !% A C-bound wrapper function for \hyperlink{cosmology.functions.F90:cosmology_functions:expansion_factor}{\tt Expansion\_Factor()}.
    implicit none
    real(kind=c_double)        :: Expansion_Factor_C
    real(kind=c_double), value :: tCosmological

    Expansion_Factor_C=Expansion_Factor(tCosmological)
    return
  end function Expansion_Factor_C

  double precision function Expansion_Rate(aExpansion)
    !% Returns the cosmological expansion rate, $\dot{a}/a$ at expansion factor {\tt aExpansion}.
    implicit none
    double precision, intent(in   ) :: aExpansion

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    Expansion_Rate=Expansion_Rate_Get(aExpansion)

    return
  end function Expansion_Rate

  double precision function Hubble_Parameter(tCosmological,aExpansion,collapsingPhase)
    !% Returns the Hubble parameter at the request cosmological time, {\tt tCosmological}, or expansion factor, {\tt aExpansion}.
    implicit none
    double precision, intent(in   ), optional :: aExpansion     , tCosmological
    logical         , intent(in   ), optional :: collapsingPhase

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    Hubble_Parameter=Hubble_Parameter_Get(tCosmological,aExpansion,collapsingPhase)

    return
  end function Hubble_Parameter

  subroutine Early_Time_Density_Scaling(dominateFactor,densityPower,aDominant,Omega_Dominant)
    !% Compute the scaling of density with expansion factor at early times in the universe.
    implicit none
    double precision, intent(in   )           :: dominateFactor
    double precision, intent(  out)           :: aDominant     , densityPower
    double precision, intent(  out), optional :: Omega_Dominant

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    call Early_Time_Density_Scaling_Get(dominateFactor,densityPower,aDominant,Omega_Dominant)

    return
  end subroutine Early_Time_Density_Scaling

  double precision function Epoch_of_Matter_Domination(dominateFactor)
    !% Compute the epoch at which matter dominates over other forms of energy by a given factor.
    implicit none
    double precision, intent(in   ) :: dominateFactor

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    Epoch_of_Matter_Domination=Epoch_of_Matter_Domination_Get(dominateFactor)

    return
  end function Epoch_of_Matter_Domination

  double precision function Omega_Matter_Total(tCosmological,aExpansion,collapsingPhase)
    !% Return the matter density parameter at expansion factor {\tt aExpansion}.
    implicit none
    double precision, intent(in   ), optional :: aExpansion     , tCosmological
    logical         , intent(in   ), optional :: collapsingPhase

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    Omega_Matter_Total=Omega_Matter_Total_Get(tCosmological,aExpansion,collapsingPhase)

    return
  end function Omega_Matter_Total

  double precision function Omega_Dark_Energy(tCosmological,aExpansion,collapsingPhase)
    !% Return the dark energy density parameter at expansion factor {\tt aExpansion}.
    implicit none
    double precision, intent(in   ), optional :: aExpansion     , tCosmological
    logical         , intent(in   ), optional :: collapsingPhase

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    Omega_Dark_Energy=Omega_Dark_Energy_Get(tCosmological,aExpansion,collapsingPhase)

    return
  end function Omega_Dark_Energy

  double precision function CMB_Temperature(tCosmological,aExpansion,collapsingPhase)
    !% Return the temperature of the cosmic microwave background at {\tt aExpansion}.
    implicit none
    double precision, intent(in   ), optional :: aExpansion     , tCosmological
    logical         , intent(in   ), optional :: collapsingPhase

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    CMB_Temperature=CMB_Temperature_Get(tCosmological,aExpansion,collapsingPhase)

    return
  end function CMB_Temperature

  double precision function Epoch_of_Matter_Dark_Energy_Equality(requestType)
    !% Return the epoch of matter-dark energy magnitude equality (either expansion factor or cosmic time).
    implicit none
    integer, intent(in   ), optional :: requestType

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    Epoch_of_Matter_Dark_Energy_Equality=Epoch_of_Matter_Dark_Energy_Equality_Get(requestType)

    return
  end function Epoch_of_Matter_Dark_Energy_Equality

  double precision function Epoch_of_Matter_Curvature_Equality(requestType)
    !% Return the epoch of matter-curvature magnitude equality (either expansion factor or cosmic time).
    implicit none
    integer, intent(in   ) :: requestType

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    Epoch_of_Matter_Curvature_Equality=Epoch_of_Matter_Curvature_Equality_Get(requestType)

    return
  end function Epoch_of_Matter_Curvature_Equality

  elemental double precision function Redshift_from_Expansion_Factor(expansionFactor)
    !% Returns redshift for a given expansion factor.
    implicit none
    double precision, intent(in   ) :: expansionFactor

    Redshift_from_Expansion_Factor=1.0d0/expansionFactor-1.0d0
    return
  end function Redshift_from_Expansion_Factor

  elemental double precision function Expansion_Factor_from_Redshift(redshift)
    !% Returns expansion factor given a redshift.
    implicit none
    double precision, intent(in   ) :: redshift

    Expansion_Factor_from_Redshift=1.0d0/(1.0d0+redshift)
    return
  end function Expansion_Factor_from_Redshift

  double precision function Comoving_Volume_Element_Redshift(time)
    !% Returns the differential comoving volume element ${\rm d}V/{\rm d}z = r_{\rm c}^2(t) {\rm c} H^{-1}(t)$ (where $r_{\rm c}$
    !% is the comoving distance to time $t$ and $H(t)$ is the Hubble parameter at that time) for unit solid angle at the
    !% specified {\tt time}.
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Physical
    implicit none
    double precision, intent(in   ) :: time

    Comoving_Volume_Element_Redshift=Comoving_Distance(time)**2*(speedLight/kilo)/Hubble_Parameter(tCosmological=time)
    return
  end function Comoving_Volume_Element_Redshift

  double precision function Comoving_Volume_Element_Time(time)
    !% Returns the differential comoving volume element ${\rm d}V/{\rm d}t = r_{\rm c}^2(t) {\rm c} a(t)$ (where $r_{\rm c}$
    !% is the comoving distance to time $t$ and $a(t)$ is the expansion at that time) for unit solid angle at the
    !% specified {\tt time}.
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    implicit none
    double precision, intent(in   ) :: time

    Comoving_Volume_Element_Time=Comoving_Distance(time)**2*(gigaYear*speedLight/megaParsec)/Expansion_Factor(time)
    return
  end function Comoving_Volume_Element_Time

  double precision function Comoving_Distance(time)
    !% Return the comoving distance to the given cosmic {\tt time}.
    implicit none
    double precision, intent(in   ) :: time

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    Comoving_Distance=Comoving_Distance_Get(time)

    return
  end function Comoving_Distance

  double precision function Time_From_Comoving_Distance(comovingDistance)
    !% Return the cosmic time corresponding to the given {\tt comovingDistance}.
    implicit none
    double precision, intent(in   ) :: comovingDistance

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    Time_From_Comoving_Distance=Time_From_Comoving_Distance_Get(comovingDistance)

    return
  end function Time_From_Comoving_Distance

  double precision function Comoving_Distance_Conversion(output,distanceModulus,redshift)
    !% Convert between different measures of comoving distance.
    implicit none
    integer         , intent(in   )           :: output
    double precision, intent(in   ), optional :: distanceModulus, redshift

    ! Initialize the module.
    call Cosmology_Functions_Initialize

    ! Get the answer using the selected method.
    Comoving_Distance_Conversion=Comoving_Distance_Conversion_Get(output,distanceModulus,redshift)

    return
  end function Comoving_Distance_Conversion

  double precision function Cosmology_Dark_Energy_Equation_Of_State(time,expansionFactor)
    !% Return the dark energy equation of state.
    implicit none
    double precision, intent(in   ), optional :: expansionFactor, time

    ! Initialize the module.
    call Cosmology_Functions_Initialize
    ! Get the answer using the selected method.
    Cosmology_Dark_Energy_Equation_Of_State=Cosmology_Dark_Energy_Equation_Of_State_Get(time,expansionFactor)
    return
  end function Cosmology_Dark_Energy_Equation_Of_State

  double precision function Cosmology_Dark_Energy_Exponent(time,expansionFactor)
    !% Return the dark energy equation of state.
    implicit none
    double precision, intent(in   ), optional :: expansionFactor, time

    ! Initialize the module.
    call Cosmology_Functions_Initialize
    ! Get the answer using the selected method.
    Cosmology_Dark_Energy_Exponent=Cosmology_Dark_Energy_Exponent_Get(time,expansionFactor)
    return
  end function Cosmology_Dark_Energy_Exponent

end module Cosmology_Functions
