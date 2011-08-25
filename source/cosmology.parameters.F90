!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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

    !$omp critical (Omega_b_Initialization)
    if (.not.Omega_b_Is_Set) then
       !@ <inputParameter>
       !@   <name>Omega_b</name>
       !@   <defaultValue>0.0455 \citep{komatsu_seven-year_2010}</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The density of baryons in the Universe in units of the critical density.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('Omega_b',Omega_b_Value,defaultValue=0.0455d0)
       Omega_b_Is_Set=.true.
    end if
    !$omp end critical (Omega_b_Initialization)

    Omega_b=Omega_b_Value
    return
  end function Omega_b

  double precision function Omega_Matter()
    !% Returns the value of $\Omega_{\rm b}$, reading it in first if necessary.
    use Galacticus_Error
    implicit none

    !$omp critical (Omega_Matter_Initialization)
    if (.not.Omega_Matter_Is_Set) then

       ! Check for deprecated parameter name.
       if (Input_Parameter_Is_Present('Omega_0')) call Galacticus_Error_Report('Omega_Matter','use of "Omega_0" in parameter file is deprecated - use "Omega_Matter" instead')

       !@ <inputParameter>
       !@   <name>Omega_Matter</name>
       !@   <defaultValue>0.2725 \citep{komatsu_seven-year_2010}</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The density of matter in the Universe in units of the critical density.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('Omega_Matter',Omega_Matter_Value,defaultValue=0.2725d0)
       Omega_Matter_Is_Set=.true.
    end if
    !$omp end critical (Omega_Matter_Initialization)

    Omega_Matter=Omega_Matter_Value
    return
  end function Omega_Matter

  double precision function Omega_DE()
    !% Returns the value of $\Omega_{\rm b}$, reading it in first if necessary.
    implicit none

    !$omp critical (Omega_DE_Initialization)
    if (.not.Omega_DE_Is_Set) then
       !@ <inputParameter>
       !@   <name>Omega_DE</name>
       !@   <defaultValue>0.7275 \citep{komatsu_seven-year_2010}</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The density of dark energy in the Universe in units of the critical density.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('Omega_DE',Omega_DE_Value,defaultValue=0.7275d0)
       Omega_DE_Is_Set=.true.
    end if
    !$omp end critical (Omega_DE_Initialization)

    Omega_DE=Omega_DE_Value
    return
  end function Omega_DE

  double precision function T_CMB()
    !% Returns the value of $T_{\rm CMB}$, reading it in first if necessary.
    implicit none

    !$omp critical (T_CMB_Initialization)
    if (.not.T_CMB_Is_Set) then
       !@ <inputParameter>
       !@   <name>T_CMB</name>
       !@   <defaultValue>2.72548 \citep{fixsen_temperature_2009}</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The present day temperature of the \CMB\ in units of Kelvin.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('T_CMB',T_CMB_Value,defaultValue=2.72548d0)
       T_CMB_Is_Set=.true.
    end if
    !$omp end critical (T_CMB_Initialization)

    T_CMB=T_CMB_Value
    return
  end function T_CMB

  double precision function Omega_Radiation()
    !% Returns the value of $\Omega_{\rm r}$, computing it first if necessary.
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    implicit none

    !$omp critical (Omega_Radiation_Initialization)
    if (.not.Omega_Radiation_Is_Set) then
       Omega_Radiation_Value=radiationConstant*(T_CMB()**4)*megaParsec**3/massSolar/speedLight**2/Critical_Density()
       Omega_Radiation_Is_Set=.true.
    end if
    !$omp end critical (Omega_Radiation_Initialization)

    Omega_Radiation=Omega_Radiation_Value
    return
  end function Omega_Radiation

  double precision function Omega_K()
    !% Returns the value of $\Omega_{\rm K}$, computing it first if necessary.
    implicit none

    !$omp critical (Omega_K_Initialization)
    if (.not.Omega_K_Is_Set) then
       Omega_K_Value=1.0d0-Omega_Matter()-Omega_DE()
       Omega_K_Is_Set=.true.
    end if
    !$omp end critical (Omega_K_Initialization)

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
    implicit none

    !$omp critical (H_0_Initialization)
    if (.not.H_0_Is_Set) then
       !@ <inputParameter>
       !@   <name>H_0</name>
       !@   <defaultValue>70.2 \citep{komatsu_seven-year_2010}</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The present day value of the Hubble parameter in units of km/s/Mpc.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('H_0',H_0_Value,defaultValue=70.2d0)
       H_0_Is_Set=.true.
    end if
    !$omp end critical (H_0_Initialization)

    H_0=H_0_Value
    return
  end function H_0

  double precision function H_0_invGyr()
    !% Returns the value of $H_0$, in units of inverse Gyr.
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    implicit none

    !$omp critical (H_0_invGyr_Initialization)
    if (.not.H_0_invGyr_Is_Set) then
       H_0_invGyr_Value=H_0()*gigaYear*kilo/megaParsec
       H_0_invGyr_Is_Set=.true.
    end if
    !$omp end critical (H_0_invGyr_Initialization)

    H_0_invGyr=H_0_invGyr_Value
    return
  end function H_0_invGyr

  double precision function Critical_Density()
    !% Returns the critical density in units of $M_\odot/$Mpc$^3$.
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    implicit none

    !$omp critical (Critical_Density_Initialization)
    if (.not.Critical_Density_Is_Set) then
       Critical_Density_Value=3.0d0*(H_0()**2)/8.0d0/Pi/gravitationalConstantGalacticus
    end if
    !$omp end critical (Critical_Density_Initialization)

    Critical_Density=Critical_Density_Value
    return
  end function Critical_Density

end module Cosmological_Parameters
