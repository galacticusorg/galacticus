!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which handles cosmological parameters.

module Cosmological_Parameters
  !% Implements cosmological parameters and related derived quantities. Default parameter values are taken from
  !% \cite{komatsu_seven-year_2010}.
  use Input_Parameters
  implicit none
  private
  public :: Omega_b, Omega_Matter, Omega_DE, Omega_Radiation, Omega_K, T_CMB, H_0, H_0_invGyr, Little_H_0, Critical_Density

  ! Stored values of cosmological parameters.
  logical          :: Omega_b_Is_Set=.false., Omega_Matter_Is_Set=.false., Omega_DE_Is_Set=.false., Omega_Radiation_Is_Set=.false.,&
       & Omega_K_Is_Set=.false., T_CMB_Is_Set=.false., H_0_Is_Set =.false., H_0_invGyr_Is_Set=.false., Critical_Density_Is_Set&
       &=.false.
  double precision :: Omega_b_Value,Omega_Matter_Value,Omega_DE_Value,Omega_Radiation_Value,Omega_K_Value,T_CMB_Value,H_0_Value&
       &,H_0_invGyr_Value,Critical_Density_Value

contains

  double precision function Omega_b()
    !% Returns the value of $\Omega_{\rm b}$, reading it in first if necessary.
    implicit none

    if (.not.Omega_b_Is_Set) then
       !$omp critical (Omega_b_Initialization)
       if (.not.Omega_b_Is_Set) then
          !@ <inputParameter>
          !@   <name>Omega_b</name>
          !@   <defaultValue>0.04568 (\citealt{story_measurement_2012}; CMB$+H_0+$BAO)</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The density of baryons in the Universe in units of the critical density.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>cosmology</group>
          !@ </inputParameter>
          call Get_Input_Parameter('Omega_b',Omega_b_Value,defaultValue=0.04568d0)
          Omega_b_Is_Set=.true.
       end if
       !$omp end critical (Omega_b_Initialization)
    end if

    Omega_b=Omega_b_Value
    return
  end function Omega_b

  double precision function Omega_Matter()
    !% Returns the value of $\Omega_{\rm b}$, reading it in first if necessary.
    use Galacticus_Error
    implicit none

    if (.not.Omega_Matter_Is_Set) then
       !$omp critical (Omega_Matter_Initialization)
       if (.not.Omega_Matter_Is_Set) then
          
          ! Check for deprecated parameter name.
          if (Input_Parameter_Is_Present('Omega_0')) call Galacticus_Error_Report('Omega_Matter','use of "Omega_0" in parameter file is deprecated - use "Omega_Matter" instead')
          
          !@ <inputParameter>
          !@   <name>Omega_Matter</name>
          !@   <defaultValue>0.2833 (\citealt{story_measurement_2012}; CMB$+H_0+$BAO)</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The density of matter in the Universe in units of the critical density.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>cosmology</group>
          !@ </inputParameter>
          call Get_Input_Parameter('Omega_Matter',Omega_Matter_Value,defaultValue=0.2833d0)
          Omega_Matter_Is_Set=.true.
       end if
       !$omp end critical (Omega_Matter_Initialization)
    end if

    Omega_Matter=Omega_Matter_Value
    return
  end function Omega_Matter

  double precision function Omega_DE()
    !% Returns the value of $\Omega_{\rm b}$, reading it in first if necessary.
    implicit none

    if (.not.Omega_DE_Is_Set) then
       !$omp critical (Omega_DE_Initialization)
       if (.not.Omega_DE_Is_Set) then
          !@ <inputParameter>
          !@   <name>Omega_DE</name>
          !@   <defaultValue>0.7167 (\citealt{story_measurement_2012}; CMB$+H_0+$BAO)</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The density of dark energy in the Universe in units of the critical density.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>cosmology</group>
          !@ </inputParameter>
          call Get_Input_Parameter('Omega_DE',Omega_DE_Value,defaultValue=0.7167d0)
          Omega_DE_Is_Set=.true.
       end if
       !$omp end critical (Omega_DE_Initialization)
    end if

    Omega_DE=Omega_DE_Value
    return
  end function Omega_DE

  double precision function T_CMB()
    !% Returns the value of $T_{\rm CMB}$, reading it in first if necessary.
    implicit none

    if (.not.T_CMB_Is_Set) then
       !$omp critical (T_CMB_Initialization)
       if (.not.T_CMB_Is_Set) then
          !@ <inputParameter>
          !@   <name>T_CMB</name>
          !@   <defaultValue>2.72548 \citep{fixsen_temperature_2009}</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The present day temperature of the \gls{cmb} in units of Kelvin.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>cosmology</group>
          !@ </inputParameter>
          call Get_Input_Parameter('T_CMB',T_CMB_Value,defaultValue=2.72548d0)
          T_CMB_Is_Set=.true.
       end if
       !$omp end critical (T_CMB_Initialization)
    end if

    T_CMB=T_CMB_Value
    return
  end function T_CMB

  double precision function Omega_Radiation()
    !% Returns the value of $\Omega_{\rm r}$, computing it first if necessary.
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    implicit none

    if (.not.Omega_Radiation_Is_Set) then
       !$omp critical (Omega_Radiation_Initialization)
       if (.not.Omega_Radiation_Is_Set) then
          Omega_Radiation_Value=radiationConstant*(T_CMB()**4)*megaParsec**3/massSolar/speedLight**2/Critical_Density()
          Omega_Radiation_Is_Set=.true.
       end if
       !$omp end critical (Omega_Radiation_Initialization)
    end if
    
    Omega_Radiation=Omega_Radiation_Value
    return
  end function Omega_Radiation

  double precision function Omega_K()
    !% Returns the value of $\Omega_{\rm K}$, computing it first if necessary.
    implicit none

    if (.not.Omega_K_Is_Set) then
       !$omp critical (Omega_K_Initialization)
       if (.not.Omega_K_Is_Set) then
          Omega_K_Value=1.0d0-Omega_Matter()-Omega_DE()
          Omega_K_Is_Set=.true.
       end if
       !$omp end critical (Omega_K_Initialization)
    end if

    Omega_K=Omega_K_Value
    return
  end function Omega_K

  double precision function Little_H_0()
    !% Returns $h_0=H_0/100$km/s/Mpc.
    implicit none
    double precision, parameter :: Big_H_0=100.0 ! km/s/Mpc.

    Little_H_0=H_0()/Big_H_0
    return
  end function Little_H_0

  double precision function H_0()
    !% Returns the value of $H_0$, reading it in first if necessary.
    use Galacticus_Display
    implicit none

    if (.not.H_0_Is_Set) then
       !$omp critical (H_0_Initialization)
       if (.not.H_0_Is_Set) then
          !@ <inputParameter>
          !@   <name>H_0</name>
          !@   <defaultValue>69.62 (\citealt{story_measurement_2012}; CMB$+H_0+$BAO)</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The present day value of the Hubble parameter in units of km/s/Mpc.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>cosmology</group>
          !@ </inputParameter>
          call Get_Input_Parameter('H_0',H_0_Value,defaultValue=69.62d0)
          ! Validate the input value.
          if (H_0_Value <= 0.0d0) call Galacticus_Display_Message("WARNING: H_0<=0 - are you sure this is what you wanted?",verbosityWarn)
          ! Record that H_0 is now set.
          H_0_Is_Set=.true.
       end if
       !$omp end critical (H_0_Initialization)
    end if

    H_0=H_0_Value
    return
  end function H_0

  double precision function H_0_invGyr()
    !% Returns the value of $H_0$, in units of inverse Gyr.
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    implicit none

    if (.not.H_0_invGyr_Is_Set) then
       !$omp critical (H_0_invGyr_Initialization)
       if (.not.H_0_invGyr_Is_Set) then
          H_0_invGyr_Value=H_0()*gigaYear*kilo/megaParsec
          H_0_invGyr_Is_Set=.true.
       end if
       !$omp end critical (H_0_invGyr_Initialization)
    end if

    H_0_invGyr=H_0_invGyr_Value
    return
  end function H_0_invGyr

  double precision function Critical_Density()
    !% Returns the critical density in units of $M_\odot/$Mpc$^3$.
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    implicit none

    if (.not.Critical_Density_Is_Set) then
       !$omp critical (Critical_Density_Initialization)
       if (.not.Critical_Density_Is_Set) then
          Critical_Density_Value=3.0d0*(H_0()**2)/8.0d0/Pi/gravitationalConstantGalacticus
          Critical_Density_Is_Set=.true.
       end if
       !$omp end critical (Critical_Density_Initialization)
    end if

    Critical_Density=Critical_Density_Value
    return
  end function Critical_Density

end module Cosmological_Parameters
