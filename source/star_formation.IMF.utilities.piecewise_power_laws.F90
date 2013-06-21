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
    double precision, intent(in   ) :: massExponent    (:), massLower(:), massUpper(:)
    double precision, intent(  out) :: imfNormalization(:)
    integer                         :: iPiece
    double precision                :: totalMass

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
          totalMass=totalMass+imfNormalization(iPiece)*log(massUpper(iPiece)/massLower(iPiece))
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
    double precision, intent(in   ) :: imfNormalization(:), mass        , massExponent(:), &
         &                             massLower       (:), massUpper(:)
    integer                         :: iPiece

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
    double precision, intent(in   ) :: imfNormalization                 (:         ), mass     (:), &
         &                             massExponent                     (:         ), massLower(:), &
         &                             massUpper                        (:         )
    double precision                :: Piecewise_Power_Law_IMF_Phi_Array(size(mass))
    integer                         :: iMass

    ! Compute IMF for all elements of array by calling the scalar version of this routine.
    forall(iMass=1:size(mass))
       Piecewise_Power_Law_IMF_Phi_Array(iMass)=Piecewise_Power_Law_IMF_Phi_Scalar(massLower,massUpper,massExponent&
            &,imfNormalization,mass(iMass))
    end forall
    return
  end function Piecewise_Power_Law_IMF_Phi_Array

end module Star_Formation_IMF_PPL
