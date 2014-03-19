!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of virial overdensity based on friends-of-friends linking length.

module Virial_Densities_Friends_of_Friends
  !% Implements calculations of virial overdensity based on friends-of-friends linking length.
  implicit none
  private
  public :: Virial_Density_Friends_of_Friends_Initialize

  ! Value of the linking length and density ratio parameters.
  double precision            :: virialDensityContrastFoFLinkingLength       , virialDensityContrastFoFDensityRatio

  ! Variables to hold the tabulated critical overdensity data.
  double precision            :: deltaTableTimeMaximum                =20.0d0, deltaTableTimeMinimum               =1.0d0
  integer         , parameter :: deltaTableNPointsPerDecade           =100

contains

  !# <virialDensityContrastMethod>
  !#  <unitName>Virial_Density_Friends_of_Friends_Initialize</unitName>
  !# </virialDensityContrastMethod>
  subroutine Virial_Density_Friends_of_Friends_Initialize(virialDensityContrastMethod,Virial_Density_Contrast_Tabulate)
    !% Initializes the $\Delta_{\rm vir}$ calculation for the fixed value implementation.
    use Input_Parameters
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    type     (varying_string                   ), intent(in   )          :: virialDensityContrastMethod
    procedure(Virial_Density_Friends_of_Friends), intent(inout), pointer :: Virial_Density_Contrast_Tabulate

    if (virialDensityContrastMethod == 'friendsOfFriends') then
       ! Return a pointer to our tabulation function.
       Virial_Density_Contrast_Tabulate => Virial_Density_Friends_of_Friends
       ! Get the linking length to use.
       !@ <inputParameter>
       !@   <name>virialDensityContrastFoFLinkingLength</name>
       !@   <defaultValue>0.2</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The friends-of-friends linking length algorithm to use in computing virial density contrast.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("virialDensityContrastFoFLinkingLength",virialDensityContrastFoFLinkingLength,defaultValue=0.2d0)
       !@ <inputParameter>
       !@   <name>virialDensityContrastFoFDensityRatio</name>
       !@   <defaultValue>4.688 (value appropriate for an \gls{nfw} profile with concentration $c=6.88$ which is the concentration found by \cite{prada_halo_2011} for halos with $\sigma=1.686$ which is the approximate critical overdensity for collapse).</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The ratio of mean virial density to density at the virial radius to assume when setting virial density contrasts in the friends-of-friends model.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("virialDensityContrastFoFDensityRatio",virialDensityContrastFoFDensityRatio,defaultValue=4.688d0)
    end if
    return
  end subroutine Virial_Density_Friends_of_Friends_Initialize

  subroutine Virial_Density_Friends_of_Friends(time,deltaVirialTable)
    !% Tabulate the virial density contrast assuming a fixed value.
    use Cosmology_Functions
    use Tables
    use Numerical_Constants_Math
    implicit none
    double precision                                      , intent(in   ) :: time
    class           (table1D                ), allocatable, intent(inout) :: deltaVirialTable
    class           (cosmologyFunctionsClass), pointer                    :: cosmologyFunctionsDefault
    integer                                                               :: deltaTableNumberPoints   , iTime
    double precision                                                      :: boundingSurfaceDensityContrast, virialDensityContrast

    ! Find minimum and maximum times to tabulate.
    deltaTableTimeMinimum=min(deltaTableTimeMinimum,time/2.0d0)
    deltaTableTimeMaximum=max(deltaTableTimeMaximum,time*2.0d0)

    ! Determine number of points to tabulate.
    deltaTableNumberPoints=int(log10(deltaTableTimeMaximum/deltaTableTimeMinimum)&
         &*dble(deltaTableNPointsPerDecade))

    ! Deallocate table if currently allocated.
    if (allocated(deltaVirialTable)) then
       call deltaVirialTable%destroy()
       deallocate(deltaVirialTable)
    end if
    allocate(table1DLogarithmicLinear :: deltaVirialTable)
    select type (deltaVirialTable)
    type is (table1DLogarithmicLinear)
       ! Get the default cosmology functions object.
       cosmologyFunctionsDefault => cosmologyFunctions()
       ! Create the table.
       call deltaVirialTable%create(deltaTableTimeMinimum,deltaTableTimeMaximum,deltaTableNumberPoints)
       ! Populate the table.
       do iTime=1,deltaTableNumberPoints
          boundingSurfaceDensityContrast=3.0d0/2.0d0/Pi/virialDensityContrastFoFLinkingLength**3
          virialDensityContrast=virialDensityContrastFoFDensityRatio*boundingSurfaceDensityContrast
          call deltaVirialTable%populate(virialDensityContrast,iTime)
       end do
    end select
    return
  end subroutine Virial_Density_Friends_of_Friends

end module Virial_Densities_Friends_of_Friends
