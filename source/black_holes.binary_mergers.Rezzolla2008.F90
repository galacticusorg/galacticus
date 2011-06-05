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


!% Contains a module which implements calculations of black hole mass and spin resulting from binary mergers utilizing the
!% approximations of \cite{rezzolla_final_2008}.

module Black_Hole_Binary_Mergers_Rezzolla
  !% Implements calculations of black hole mass and spin resulting from binary mergers utilizing the approximations of
  !% \cite{rezzolla_final_2008}.
  private
  public :: Black_Hole_Binary_Merger_Rezzolla_Initialize

contains

  !# <blackHoleBinaryMergersMethod>
  !#  <unitName>Black_Hole_Binary_Merger_Rezzolla_Initialize</unitName>
  !# </blackHoleBinaryMergersMethod>
  subroutine Black_Hole_Binary_Merger_Rezzolla_Initialize(blackHoleBinaryMergersMethod,Black_Hole_Binary_Merger_Do)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: blackHoleBinaryMergersMethod
    procedure(),          pointer, intent(inout) :: Black_Hole_Binary_Merger_Do
    
    if (blackHoleBinaryMergersMethod == 'Rezzolla2008') Black_Hole_Binary_Merger_Do => Black_Hole_Binary_Merger_Rezzolla
    return
  end subroutine Black_Hole_Binary_Merger_Rezzolla_Initialize

  subroutine Black_Hole_Binary_Merger_Rezzolla(blackHoleMassA,blackHoleMassB,blackHoleSpinA,blackHoleSpinB,blackHoleMassFinal&
       &,blackHoleSpinFinal)
    !% Computes the mass and spin of a black hole resulting from a binary merger utilizing the approximations of
    !% \cite{rezzolla_final_2008}.
    implicit none
    double precision, intent(in)  :: blackHoleMassA,blackHoleMassB,blackHoleSpinA,blackHoleSpinB
    double precision, intent(out) :: blackHoleMassFinal,blackHoleSpinFinal
    ! Parameters of the fitting functions used by Rezzolla et al.
    double precision, parameter   :: s4=-0.129d0, s5=-0.384d0, t0=-2.686d0, t2=-3.454d0, t3=2.353d0
    ! Assumed fixed values for the various alignment angles.
    double precision, parameter   :: cosinePhi=1.0d0, cosineTheta=1.0d0, cosineXi=1.0d0
    double precision              :: blackHoleMass1,blackHoleMass2,blackHoleSpin1,blackHoleSpin2,massRatio,symmetricMassRatio&
         &,orbitalAngularMomentum

    ! Find which is the more massive of the two black holes.
    if (blackHoleMassA < blackHoleMassB) then
       blackHoleMass1=blackHoleMassB
       blackHoleMass2=blackHoleMassA
       blackHoleSpin1=blackHoleSpinB
       blackHoleSpin2=blackHoleSpinA
    else
       blackHoleMass1=blackHoleMassA
       blackHoleMass2=blackHoleMassB
       blackHoleSpin1=blackHoleSpinA
       blackHoleSpin2=blackHoleSpinB
    end if

    ! Assume negligible mass loss (Rezzolla et al. report 5-7% mass loss via gravitational radiation).
    blackHoleMassFinal=blackHoleMass1+blackHoleMass2

    ! Compute the black hole mass ratio.
    massRatio=blackHoleMass2/blackHoleMass1
    symmetricMassRatio=massRatio/(1.0d0+massRatio)**2

    ! Evaluate the fitting formula for the orbital angular momentum (which accounts for losses to gravitational waves).
    orbitalAngularMomentum=s4*(blackHoleSpin1**2+blackHoleSpin2**2*massRatio**4+2.0d0*blackHoleSpin1*blackHoleSpin2*massRatio**2&
         &*cosinePhi)/(1.0d0+massRatio**2)**2+((s5*symmetricMassRatio+t0+2.0d0)/(1.0d0+massRatio**2))*(blackHoleSpin1*cosineTheta&
         &+blackHoleSpin2*massRatio**2*cosinePhi)+2.0d0*dsqrt(3.0d0)+t2*symmetricMassRatio+t3*symmetricMassRatio**2

    ! Compute the spin of the final black hole by evaluating the fitting formula.
    blackHoleSpinFinal=dsqrt(blackHoleSpin1**2+blackHoleSpin2**2*massRatio**4+2.0d0*blackHoleSpin2*blackHoleSpin1*massRatio**2&
         &*cosinePhi+2.0d0*(blackHoleSpin1*cosineTheta+blackHoleSpin2*massRatio**2*cosineXi)*orbitalAngularMomentum*massRatio&
         &+(orbitalAngularMomentum*massRatio)**2)/(1.0d0+massRatio)**2
    
    return
  end subroutine Black_Hole_Binary_Merger_Rezzolla

end module Black_Hole_Binary_Mergers_Rezzolla
