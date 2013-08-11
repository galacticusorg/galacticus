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

!% Contains a module which generates a tabulated Tinker2008 halo mass function.

module Halo_Mass_Function_Tinker2008
  !% Implements generation of a tabulated power-law primordial power spectrum.
  implicit none
  private
  public :: Halo_Mass_Function_Tinker2008_Initialize

  ! Variables to hold the table of parameters vs. overdensity.
  integer                                     :: deltaTableNumberPoints
  double precision, allocatable, dimension(:) :: deltaTableA            , deltaTableB    , &
       &                                         deltaTableC            , deltaTableDelta, &
       &                                         deltaTableNormalization

contains

  !# <haloMassFunctionMethod>
  !#  <unitName>Halo_Mass_Function_Tinker2008_Initialize</unitName>
  !# </haloMassFunctionMethod>
  subroutine Halo_Mass_Function_Tinker2008_Initialize(haloMassFunctionMethod,Halo_Mass_Function_Differential_Get)
    !% Initializes the ``Tinker2008 mass functon'' module.
    use ISO_Varying_String
    use FoX_dom
    use IO_XML
    use Galacticus_Error
    use Galacticus_Input_Paths
    implicit none
    type            (varying_string  ), intent(in   )          :: haloMassFunctionMethod
    procedure       (double precision), intent(inout), pointer :: Halo_Mass_Function_Differential_Get
    type            (Node            )               , pointer :: columnElement                      , columnsElement, &
         &                                                        doc
    integer                                                    :: ioErr

    if (haloMassFunctionMethod == 'Tinker2008') then
       Halo_Mass_Function_Differential_Get => Halo_Mass_Function_Differential_Tinker2008
       ! Read the data file which gives fitting parameters as a function of halo overdensity.
       !$omp critical (FoX_DOM_Access)
       doc => parseFile(char(Galacticus_Input_Path())//"data/darkMatter/Halo_Mass_Function_Parameters_Tinker_2008.xml",iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Halo_Mass_Function_Tinker2008_Initialize','Unable to find data file')
       columnsElement => XML_Get_First_Element_By_Tag_Name(doc           ,"columns"      )
       columnElement  => XML_Get_First_Element_By_Tag_Name(columnsElement,"overdensity"  )
       call XML_Array_Read(columnElement,"data",deltaTableDelta        )
       columnElement  => XML_Get_First_Element_By_Tag_Name(columnsElement,"normalization")
       call XML_Array_Read(columnElement,"data",deltaTableNormalization)
       columnElement  => XML_Get_First_Element_By_Tag_Name(columnsElement,"parameterA"   )
       call XML_Array_Read(columnElement,"data",deltaTableA            )
       columnElement  => XML_Get_First_Element_By_Tag_Name(columnsElement,"parameterB"   )
       call XML_Array_Read(columnElement,"data",deltaTableB            )
       columnElement  => XML_Get_First_Element_By_Tag_Name(columnsElement,"parameterC"   )
       call XML_Array_Read(columnElement,"data",deltaTableC            )
       deltaTableNumberPoints=size(deltaTableDelta)
       ! Destroy the document.
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)
    end if
    return
  end subroutine Halo_Mass_Function_Tinker2008_Initialize

  double precision function Halo_Mass_Function_Differential_Tinker2008(time,mass)
    !% Compute the \cite{tinker_towardhalo_2008} halo mass function.
    use Power_Spectra
    use Virial_Density_Contrast
    use Cosmological_Parameters
    use Cosmology_Functions
    use FGSL
    use Numerical_Interpolation
    use Linear_Growth
    implicit none
    double precision                   , intent(in   ) :: mass                           , time
    double precision                                   :: alpha                          , sigma
    type            (fgsl_interp      ), save          :: interpolationObject
    type            (fgsl_interp_accel), save          :: interpolationAccelerator
    logical                            , save          :: resetInterpolation      =.true.
    !$omp threadprivate(interpolationObject,interpolationAccelerator,resetInterpolation)
    double precision                   , save          :: timePrevious            =-1.0d0
    !$omp threadprivate(timePrevious)
    double precision                   , save          :: Delta                          , a             , &
         &                                                a0                             , alphaDelta    , &
         &                                                b                              , b0            , &
         &                                                c                              , c0            , &
         &                                                expansionFactor                , growthFactor  , &
         &                                                normalization                  , normalization0
    !$omp threadprivate(expansionFactor,Delta,growthFactor,normalization0,a0,b0,c0,normalization,a,alphaDelta,b,c)
    ! Update fitting function parameters if the time differs from that on the previous call.
    if (time /= timePrevious) then
       ! Get halo virial density contrast, expansion factor and growth factor.
       expansionFactor=Expansion_Factor            (time)
       Delta          =Halo_Virial_Density_Contrast(time)
       growthFactor   =Linear_Growth_Factor        (time)

       ! Compute coefficients of fitting function.
       normalization0=Interpolate(deltaTableNumberPoints,deltaTableDelta,deltaTableNormalization &
            &,interpolationObject,interpolationAccelerator,Delta,reset=resetInterpolation,extrapolationType=extrapolationTypeLinear)
       a0            =Interpolate(deltaTableNumberPoints,deltaTableDelta,deltaTableA &
            &,interpolationObject,interpolationAccelerator,Delta,reset=resetInterpolation,extrapolationType=extrapolationTypeLinear)
       b0            =Interpolate(deltaTableNumberPoints,deltaTableDelta,deltaTableB &
            &,interpolationObject,interpolationAccelerator,Delta,reset=resetInterpolation,extrapolationType=extrapolationTypeLinear)
       c0            =Interpolate(deltaTableNumberPoints,deltaTableDelta,deltaTableC &
            &,interpolationObject,interpolationAccelerator,Delta,reset=resetInterpolation,extrapolationType=extrapolationTypeLinear)

       ! Extrapolate to higher redshift using redshift scalings given by Tinker et al. (2008; eqns. 5-8).
       normalization=normalization0*expansionFactor**0.14d0
       a=a0*expansionFactor**0.06d0
       alphaDelta=10.0d0**(-(0.75d0/log10(Delta/75.0d0))**1.2d0)
       b=b0*expansionFactor**alphaDelta
       c=c0

       ! Store the time.
       timePrevious=time
    end if

    ! Compute the mass function.
    sigma=Cosmological_Mass_Root_Variance(mass)*growthFactor
    alpha=abs(Cosmological_Mass_Root_Variance_Logarithmic_Derivative(mass))
    Halo_Mass_Function_Differential_Tinker2008=(Omega_Matter()*Critical_Density()/mass**2)*alpha*normalization*exp(-c/sigma**2)&
         &*(1.0d0+(b/sigma)**a)
    return
  end function Halo_Mass_Function_Differential_Tinker2008

end module Halo_Mass_Function_Tinker2008
