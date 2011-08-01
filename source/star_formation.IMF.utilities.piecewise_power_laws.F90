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


!% Contains a module which implements calculations of piecewise power-law initial mass functions.

module Star_Formation_IMF_PPL
  !% Implements calculations of piecewise power-law initial mass functions.
  implicit none
  private
  public :: Piecewise_Power_Law_IMF_Normalize, Piecewise_Power_Law_IMF_Phi
  
  ! Generic interface to Phi functions.
  interface Piecewise_Power_Law_IMF_Phi
     module procedure Piecewise_Power_Law_IMF_Phi_Scalar
     module procedure Piecewise_Power_Law_IMF_Phi_Array
  end interface

contains

  subroutine Piecewise_Power_Law_IMF_Normalize(massLower,massUpper,massExponent,imfNormalization)
    !% Computes normalizations for the pieces of a piecewise power-law IMF such that the IMF is continuous and normalized to unit
    !% stellar mass.
    implicit none
    double precision, intent(in)  :: massLower(:),massUpper(:),massExponent(:)
    double precision, intent(out) :: imfNormalization(:)
    integer                       :: iPiece
    double precision              :: totalMass

    ! Loop over each piece of the IMF.
    totalMass=0.0d0
    do iPiece=1,size(massLower)
       ! Ensure continuity in the IMF.
       if (iPiece == 1) then
          ! Arbitrary normalization for first piece.
          imfNormalization(iPiece)=1.0d0
       else
          ! Normalize to get continuity at the lower boundary.
          imfNormalization(iPiece)=imfNormalization(iPiece-1)*(massUpper(iPiece-1)**massExponent(iPiece-1))/(massLower(iPiece)&
               &**massExponent(iPiece))
       end if
       ! Sum the mass contributed by this range.
       if (massExponent(iPiece) == -2.0d0) then
          ! Integral is logarithmic in this case.
          totalMass=totalMass+imfNormalization(iPiece)*dlog(massUpper(iPiece)/massLower(iPiece))
       else
          ! Integral is a power-law in all other cases.
          totalMass=totalMass+imfNormalization(iPiece)*(massUpper(iPiece)**(2.0d0+massExponent(iPiece))-massLower(iPiece) &
               &**(2.0d0+massExponent(iPiece)))/(2.0d0+massExponent(iPiece))
       end if
    end do
    ! Divide through by the total mass to get the correct normalization.
    imfNormalization=imfNormalization/totalMass
    return
  end subroutine Piecewise_Power_Law_IMF_Normalize
  
  pure double precision function Piecewise_Power_Law_IMF_Phi_Scalar(massLower,massUpper,massExponent,imfNormalization,mass)
    !% Returns the IMF at given {\tt mass}.
    implicit none
    double precision, intent(in) :: massLower(:),massUpper(:),massExponent(:),imfNormalization(:),mass
    integer                      :: iPiece

    ! Assume zero value by default.
    Piecewise_Power_Law_IMF_Phi_Scalar=0.0d0
    ! Loop over all pieces of the IMF.
    do iPiece=1,size(massLower)
       ! Check if mass is within the range for this piece.
       if (mass >= massLower(iPiece) .and. mass < massUpper(iPiece)) then
          ! It is, so compute IMF and exit the loop.
          Piecewise_Power_Law_IMF_Phi_Scalar=imfNormalization(iPiece)*(mass**massExponent(iPiece))
          exit
       end if
    end do
    return
  end function Piecewise_Power_Law_IMF_Phi_Scalar

  function Piecewise_Power_Law_IMF_Phi_Array(massLower,massUpper,massExponent,imfNormalization,mass)
    !% Returns the IMF at given {\tt mass()}.
    implicit none
    double precision, intent(in) :: massLower(:),massUpper(:),massExponent(:),imfNormalization(:),mass(:)
    double precision             :: Piecewise_Power_Law_IMF_Phi_Array(size(mass))
    integer                      :: iMass

    ! Compute IMF for all elements of array by calling the scalar version of this routine.
    forall(iMass=1:size(mass))
       Piecewise_Power_Law_IMF_Phi_Array(iMass)=Piecewise_Power_Law_IMF_Phi_Scalar(massLower,massUpper,massExponent&
            &,imfNormalization,mass(iMass))
    end forall
    return
  end function Piecewise_Power_Law_IMF_Phi_Array

end module Star_Formation_IMF_PPL
