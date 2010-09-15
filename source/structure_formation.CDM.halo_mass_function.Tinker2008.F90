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


!% Contains a module which generates a tabulated Tinker2008 halo mass function.

module Halo_Mass_Function_Tinker2008
  !% Implements generation of a tabulated power-law primordial power spectrum.
  private
  public :: Halo_Mass_Function_Tinker2008_Initialize

  ! Parameters controlling the gridding of the power spectrum and default wavenumber range.
  integer,          parameter :: nPointsPerDecade=1000
  double precision            :: logMassMinimum=dlog(1.0d9), logMassMaximum=dlog(1.0d15)

  ! Variables to hold the table of parameters vs. overdensity.
  integer                                     :: deltaTableNumberPoints
  double precision, allocatable, dimension(:) :: deltaTableDelta,deltaTableNormalization,deltaTableA,deltaTableB,deltaTableC

contains
  
  !# <haloMassFunctionMethod>
  !#  <unitName>Halo_Mass_Function_Tinker2008_Initialize</unitName>
  !# </haloMassFunctionMethod>
  subroutine Halo_Mass_Function_Tinker2008_Initialize(haloMassFunctionMethod,Halo_Mass_Function_Tabulate)
    !% Initializes the ``Tinker2008 mass functon'' module.
    use ISO_Varying_String
    use FoX_dom
    use Galacticus_Error
    use Memory_Management
    implicit none
    type(varying_string),          intent(in)    :: haloMassFunctionMethod
    procedure(),          pointer, intent(inout) :: Halo_Mass_Function_Tabulate
    type(Node),           pointer                :: doc,columnsElement,columnElement,datum
    type(NodeList),       pointer                :: deltaList,normalizationList,parameterAList,parameterBList,parameterCList
    integer                                      :: iDatum,ioErr
    double precision                             :: datumValue

    if (haloMassFunctionMethod.eq.'Tinker2008') then
       Halo_Mass_Function_Tabulate => Halo_Mass_Function_Tinker2008_Tabulate
       
       ! Read the data file which gives fitting parameters as a function of halo overdensity.
       !$omp critical (FoX_DOM_Access)
       doc => parseFile("./data/Halo_Mass_Function_Parameters_Tinker_2008.xml",iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Halo_Mass_Function_Tinker2008_Initialize','Unable to find data file')
       columnsElement    => item(getElementsByTagname(doc           ,"columns"      ),0)
       columnElement     => item(getElementsByTagname(columnsElement,"overdensity"  ),0)
       deltaList         =>      getElementsByTagname(columnElement ,"data"         )
       columnElement     => item(getElementsByTagname(columnsElement,"normalization"),0)
       normalizationList =>      getElementsByTagname(columnElement ,"data"         )
       columnElement     => item(getElementsByTagname(columnsElement,"parameterA"   ),0)
       parameterAList    =>      getElementsByTagname(columnElement ,"data"         )
       columnElement     => item(getElementsByTagname(columnsElement,"parameterB"   ),0)
       parameterBList    =>      getElementsByTagname(columnElement ,"data"         )
       columnElement     => item(getElementsByTagname(columnsElement,"parameterC"   ),0)
       parameterCList    =>      getElementsByTagname(columnElement ,"data"         )
       deltaTableNumberPoints=getLength(deltaList)
       call Alloc_Array(deltaTableDelta        ,[deltaTableNumberPoints])
       call Alloc_Array(deltaTableNormalization,[deltaTableNumberPoints])
       call Alloc_Array(deltaTableA            ,[deltaTableNumberPoints])
       call Alloc_Array(deltaTableB            ,[deltaTableNumberPoints])
       call Alloc_Array(deltaTableC            ,[deltaTableNumberPoints])
       do iDatum=0,getLength(deltaList)-1
          datum => item(deltaList,iDatum)
          call extractDataContent(datum,datumValue)
          deltaTableDelta(iDatum+1)=datumValue
          datum => item(normalizationList,iDatum)
          call extractDataContent(datum,datumValue)
          deltaTableNormalization(iDatum+1)=datumValue
          datum => item(parameterAList,iDatum)
          call extractDataContent(datum,datumValue)
          deltaTableA(iDatum+1)=datumValue
          datum => item(parameterBList,iDatum)
          call extractDataContent(datum,datumValue)
          deltaTableB(iDatum+1)=datumValue
          datum => item(parameterCList,iDatum)
          call extractDataContent(datum,datumValue)
          deltaTableC(iDatum+1)=datumValue
       end do
       ! Destroy the document.
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)
    end if
    return
  end subroutine Halo_Mass_Function_Tinker2008_Initialize

  subroutine Halo_Mass_Function_Tinker2008_Tabulate(time,logMass,haloMassFunctionNumberPoints,haloMassFunctionLogMass &
       &,haloMassFunctionLogAbundance)
    !% Tabulate a \cite{tinker_towardhalo_2008} halo mass function.
    use Memory_Management
    use Numerical_Ranges
    use Numerical_Constants_Math
    use CDM_Power_Spectrum
    use Virial_Density_Contrast
    use Cosmological_Parameters
    use Cosmology_Functions
    use FGSL
    use Numerical_Interpolation
    use Linear_Growth
    implicit none
    double precision,                            intent(in)    :: time,logMass
    double precision, allocatable, dimension(:), intent(inout) :: haloMassFunctionLogMass,haloMassFunctionLogAbundance
    integer,                                     intent(out)   :: haloMassFunctionNumberPoints
    integer                                                    :: iMass
    double precision                                           :: mass,sigma,alpha,normalization,a,b,c,expansionFactor,Delta,a0&
         &,b0,c0,alphaDelta,normalization0,growthFactor
    type(fgsl_interp),       save                              :: interpolationObject
    type(fgsl_interp_accel), save                              :: interpolationAccelerator
    logical,                 save                              :: resetInterpolation=.true.
    !$omp threadprivate(interpolationObject,interpolationAccelerator,resetInterpolation)

    ! Determine range of masss required.
    logMassMinimum=min(logMassMinimum,logMass-ln10)
    logMassMaximum=max(logMassMaximum,logMass+ln10)
    
    ! Determine number of points to tabulate.
    haloMassFunctionNumberPoints=int((logMassMaximum-logMassMinimum)*dble(nPointsPerDecade)/ln10)

    ! Deallocate arrays if currently allocated.
    if (allocated(haloMassFunctionLogMass))      call Dealloc_Array(haloMassFunctionLogMass     )
    if (allocated(haloMassFunctionLogAbundance)) call Dealloc_Array(haloMassFunctionLogAbundance)
    ! Allocate the arrays to current required size.
    call Alloc_Array(haloMassFunctionLogMass     ,[haloMassFunctionNumberPoints])
    call Alloc_Array(haloMassFunctionLogAbundance,[haloMassFunctionNumberPoints])

    expansionFactor=Expansion_Factor            (time)
    Delta          =Halo_Virial_Density_Contrast(time)
    growthFactor   =Linear_Growth_Factor        (time)

    normalization0=Interpolate(deltaTableNumberPoints,deltaTableDelta,deltaTableNormalization &
         &,interpolationObject,interpolationAccelerator,Delta,reset=resetInterpolation)   
    a0            =Interpolate(deltaTableNumberPoints,deltaTableDelta,deltaTableA &
         &,interpolationObject,interpolationAccelerator,Delta,reset=resetInterpolation) 
    b0            =Interpolate(deltaTableNumberPoints,deltaTableDelta,deltaTableB &
         &,interpolationObject,interpolationAccelerator,Delta,reset=resetInterpolation) 
    c0            =Interpolate(deltaTableNumberPoints,deltaTableDelta,deltaTableC &
         &,interpolationObject,interpolationAccelerator,Delta,reset=resetInterpolation) 
    
    ! Extrapolate to higher redshift using redshift scalings given by Tinker et al. (2008; eqns. 5-8).
    normalization=normalization0*expansionFactor**0.14d0
    a=a0*expansionFactor**0.06
    alphaDelta=10.0d0**(-(0.75d0/dlog10(Delta/75.0d0))**1.2d0)
    b=b0*expansionFactor**alphaDelta
    c=c0
    
    ! Tabulate the function.
    haloMassFunctionLogMass=Make_Range(logMassMinimum,logMassMaximum,haloMassFunctionNumberPoints,rangeTypeLinear)
    do iMass=1,haloMassFunctionNumberPoints
       mass=dexp(haloMassFunctionLogMass(iMass))
       sigma=sigma_CDM(mass)*growthFactor
       alpha=dabs(sigma_CDM_Logarithmic_Derivative(mass))
       haloMassFunctionLogAbundance(iMass)=(Omega_0()*Critical_Density()/mass**2)*alpha*normalization*dexp(-c/sigma**2)*(1.0d0+(b&
            &/sigma)**a)
    end do
    where (haloMassFunctionLogAbundance > 0.0d0)
       haloMassFunctionLogAbundance=dlog(haloMassFunctionLogAbundance)
    elsewhere
       haloMassFunctionLogAbundance=-1000.0d0
    end where
    return
  end subroutine Halo_Mass_Function_Tinker2008_Tabulate
  
end module Halo_Mass_Function_Tinker2008
