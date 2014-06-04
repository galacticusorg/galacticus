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

!% Contains a module which provides an object that implements stellar spectra postprocessors.

module Stellar_Population_Spectra_Postprocess
  !% Provides an object that implements stellar spectra postprocessors.
  use, intrinsic :: ISO_C_Binding
  use               ISO_Varying_String
  !# <include directive="spectraPostprocessor" type="functionModules" >
  include 'spectraPostprocessor.functionModules.inc'
  !# </include>
  private
  public :: Stellar_Population_Spectrum_Postprocess, Stellar_Population_Spectrum_Postprocess_Index, Stellar_Population_Spectrum_Postprocess_Chain_Methods

  ! A chain of postprocessing algorithms.
  type postprocessor
     class(spectraPostprocessorClass), pointer :: method
  end type postprocessor
  type postprocessors
     type(varying_string)                            :: methodsLabel
     type(postprocessor ), allocatable, dimension(:) :: postprocess
  end type postprocessors
  ! Initialization state.
  logical                                            :: stellarPopulationSpectraPostprocessInitialized=.false.
  ! Array of postprocessing chains.
  type   (postprocessors), allocatable, dimension(:) :: postprocessingChains
  type   (varying_string), allocatable, dimension(:) :: postprocessingChainNames

  !# <include directive="spectraPostprocessor" type="function" >
  !#  <descriptiveName>Spectra Postprocessor</descriptiveName>
  !#  <description>Object providing postprocessors for spectra.</description>
  !#  <default>null</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <stateful>no</stateful>
  !#  <method name="apply" >
  !#   <description>Apply postprocessing to a spectrum.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: wavelength, age, redshift</argument>
  !#   <argument>double precision, intent(inout) :: modifier</argument>
  !#  </method>
  include 'spectraPostprocessor.type.inc'
  !# </include>

  function Stellar_Population_Spectrum_Postprocess_Chain_Methods(postprocessingChainIndex)
    !% Return a label describing the postprocessing methods used in the indicated chain.
    implicit none
    type   (varying_string)                :: Stellar_Population_Spectrum_Postprocess_Chain_Methods
    integer                , intent(in   ) :: postprocessingChainIndex

    Stellar_Population_Spectrum_Postprocess_Chain_Methods=postprocessingChains(postprocessingChainIndex)%methodsLabel
    return
  end function Stellar_Population_Spectrum_Postprocess_Chain_Methods

  integer function Stellar_Population_Spectrum_Postprocess_Index(postprocessingChain)
    !% Return the index to the specified postprocessing chain.
    use String_Handling
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="stellarPopulationSpectraPostprocessInitialize" type="moduleUse">
    include 'stellar_populations.spectra.postprocess.initialize.modules.inc'
    !# </include>
    implicit none
    type   (varying_string), intent(in   )               :: postprocessingChain
    type   (varying_string), allocatable  , dimension(:) :: postprocessingChainNamesTemporary
    type   (postprocessors), allocatable  , dimension(:) :: postprocessingChainsTemporary
    integer                                              :: i                                , methodCount
    type   (varying_string)                              :: parameterName

    ! Check whether we already have this chain loaded.
    if (allocated(postprocessingChainNames)) then
       if (.not.any(postprocessingChainNames == postprocessingChain)) then
          call Move_Alloc(postprocessingChainNames,postprocessingChainNamesTemporary)
          call Move_Alloc(postprocessingChains    ,postprocessingChainsTemporary    )
          allocate(postprocessingChainNames(size(postprocessingChainNamesTemporary)+1))
          allocate(postprocessingChains    (size(postprocessingChainNamesTemporary)+1))
          postprocessingChainNames(1:size(postprocessingChainNamesTemporary))=postprocessingChainNamesTemporary
          postprocessingChains    (1:size(postprocessingChainNamesTemporary))=postprocessingChainsTemporary
          postprocessingChainNames(size(postprocessingChainNames))=postprocessingChain
          Stellar_Population_Spectrum_Postprocess_Index=size(postprocessingChainNames)
          deallocate(postprocessingChainNamesTemporary)
          deallocate(postprocessingChainsTemporary    )
       else
          do Stellar_Population_Spectrum_Postprocess_Index=1,size(postprocessingChainNames)
             if (postprocessingChainNames(Stellar_Population_Spectrum_Postprocess_Index) == postprocessingChain) return
          end do
       end if
    else
       allocate(postprocessingChainNames(1))
       allocate(postprocessingChains    (1))
       postprocessingChainNames(1)=postprocessingChain
       Stellar_Population_Spectrum_Postprocess_Index=1
    end if
    ! We do not have this chain loaded. Load it now.
    ! Get the stellar population postprocessing methods parameter.
    !@ <inputParameter>
    !@   <regEx>stellarPopulationSpectraPostprocess[a-zA-Z0-9_]+Methods</regEx>
    !@   <defaultValue>inoue2014 (for ``Default'' chain)</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The name of methods to be used for post-processing of stellar population spectra.
    !@   </description>
    !@   <type>string</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    ! Construct the name of the parameter to read.
    if (postprocessingChain == "default") then
       parameterName='stellarPopulationSpectraPostprocess'//String_Upper_Case_First(char(postprocessingChain))//'Methods'
    else
       parameterName='stellarPopulationSpectraPostprocess'//String_Upper_Case_First(char(postprocessingChain))//'Methods'
    end if
    ! Determine how many methods are to be applied.
    methodCount=Get_Input_Parameter_Array_Size(char(parameterName))
    ! Allocate methods array and read method names.
    if (methodCount > 0) then
       allocate(postprocessingChainNamesTemporary(methodCount))
       call Get_Input_Parameter(char(parameterName),postprocessingChainNamesTemporary)
    else
       if (postprocessingChain == "default") then
          methodCount=1
          allocate(postprocessingChainNamesTemporary(methodCount))
          call Get_Input_Parameter(char(parameterName),postprocessingChainNamesTemporary,defaultValue=['inoue2014'])
       else
          call Galacticus_Error_Report('Stellar_Population_Spectrum_Postprocess_Index','parameter ['//parameterName//'] is not present in parameter file')
       end if
    end if
    ! Store the list of method names.
    postprocessingChains(Stellar_Population_Spectrum_Postprocess_Index)%methodsLabel=String_Join(postprocessingChainNamesTemporary,":")
    ! Allocate postprocessors.
    allocate(postprocessingChains(Stellar_Population_Spectrum_Postprocess_Index)%postprocess(methodCount))
    do i=1,methodCount
       nullify(postprocessingChains(Stellar_Population_Spectrum_Postprocess_Index)%postprocess(i)%method)
       ! Create the appropriate postprocessor.
       postprocessingChains(Stellar_Population_Spectrum_Postprocess_Index)%postprocess(i)%method => &
            & spectraPostprocessor(char(postprocessingChainNamesTemporary(i)))
    end do
    deallocate(postprocessingChainNamesTemporary)
    return
  end function Stellar_Population_Spectrum_Postprocess_Index

  subroutine Stellar_Population_Spectrum_Postprocess_Initialize
    !% Initialize the stellar population spectra postprocessing module
    implicit none
    integer                 :: defaultIndex
    type   (varying_string) :: postprocessingChain

    ! Initialize if necessary.
    if (.not.stellarPopulationSpectraPostprocessInitialized) then
       !$omp critical(Stellar_Population_Spectrum_Postprocess_Initialization)
       if (.not.stellarPopulationSpectraPostprocessInitialized) then
          ! Initialize by ensuring that the default postprocessing chain is loaded.
          postprocessingChain="default"
          defaultIndex=Stellar_Population_Spectrum_Postprocess_Index(postprocessingChain)
          stellarPopulationSpectraPostprocessInitialized=.true.
       end if
       !$omp end critical(Stellar_Population_Spectrum_Postprocess_Initialization)
    end if
    return
  end subroutine Stellar_Population_Spectrum_Postprocess_Initialize

  double precision function Stellar_Population_Spectrum_Postprocess(postprocessingChainIndex,wavelength,age,redshift)
    !% Return a multiplicative factor by which a stellar population spectrum should be modified by any postprocessing.
    implicit none
    integer         , intent(in   ) :: postprocessingChainIndex
    double precision, intent(in   ) :: age                     , redshift, &
         &                             wavelength
    integer                         :: i

    ! Initialize the module.
    call Stellar_Population_Spectrum_Postprocess_Initialize()

    ! Compute the postprocessing factor.
    Stellar_Population_Spectrum_Postprocess=1.0d0
    do i=1,size(postprocessingChains(postprocessingChainIndex)%postprocess)
       call postprocessingChains(postprocessingChainIndex)%postprocess(i)%method%apply(wavelength,age,redshift,Stellar_Population_Spectrum_Postprocess)
    end do
    return
  end function Stellar_Population_Spectrum_Postprocess

end module Stellar_Population_Spectra_Postprocess
