!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!!{
Contains a module of globally-accessible functions supporting the \refClass{accretionDiskSpectraClass} class.
!!}

module Accretion_Disk_Spectra_Utilities
  !!{
  Provides globally-accessible functions supporting the \refClass{accretionDiskSpectraClass} class.
  !!}
  private
  public :: accretionDiskSpectraConstruct       , accretionDiskSpectraDestruct    , accretionDiskSpectraDeepCopy  , accretionDiskSpectraDeepCopyReset, &
       &    accretionDiskSpectraDeepCopyFinalize, accretionDiskSpectraStateRestore, accretionDiskSpectraStateStore, accretionDiskSpectraSpectrumNode

  ! Module-scope pointer to our task object. This is used for reference counting so that debugging information is consistent
  ! between the increments and decrements.
  class(*), pointer :: accretionDiskSpectra__
  !$omp threadprivate(accretionDiskSpectra__)

contains

  !![
  <functionGlobal>
   <unitName>accretionDiskSpectraConstruct</unitName>
   <type>void</type>
   <module>Input_Parameters, only : inputParameters</module>
   <arguments>type (inputParameters), intent(inout), target  :: parameters           </arguments>
   <arguments>class(*              ), intent(  out), pointer :: accretionDiskSpectra_</arguments>
  </functionGlobal>
  !!]
  subroutine accretionDiskSpectraConstruct(parameters,accretionDiskSpectra_)
    !!{
    Build a {\normalfont \ttfamily accretionDiskSpectra} object from a given parameter set. This is a globally-callable function
    to allow us to subvert the class/module hierarchy.
    !!}
    use :: Error                  , only : Error_Report
    use :: Input_Parameters       , only : inputParameter           , inputParameters
    use :: Accretion_Disk_Spectra , only : accretionDiskSpectraClass, accretionDiskSpectra
    implicit none
    type (inputParameters), intent(inout), target  :: parameters
    class(*              ), intent(  out), pointer :: accretionDiskSpectra_
    type (inputParameters)               , pointer :: parametersCurrent

    parametersCurrent => parameters
    do while (.not.parametersCurrent%isPresent('accretionDiskSpectra').and.associated(parametersCurrent%parent))
       parametersCurrent => parametersCurrent%parent
    end do
    if (.not.parametersCurrent%isPresent('accretionDiskSpectra')) parametersCurrent => parameters
    accretionDiskSpectra__ => accretionDiskSpectra(parametersCurrent)
    select type (accretionDiskSpectra__)
    class is (accretionDiskSpectraClass) 
       !![
       <referenceCountIncrement object="accretionDiskSpectra__"/>
       !!]
       call accretionDiskSpectra__%autoHook()
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select
    accretionDiskSpectra_ => accretionDiskSpectra__
    return
  end subroutine accretionDiskSpectraConstruct

  !![
  <functionGlobal>
   <unitName>accretionDiskSpectraSpectrumNode</unitName>
   <type>double precision</type>
   <module>Galacticus_Nodes, only : treeNode</module>
   <arguments>class           (*       ), intent(inout) :: accretionDiskSpectra_</arguments>
   <arguments>type            (treeNode), intent(inout) :: node                 </arguments>
   <arguments>double precision          , intent(in   ) :: wavelength           </arguments>
  </functionGlobal>
  !!]
  double precision function accretionDiskSpectraSpectrumNode(accretionDiskSpectra_,node,wavelength) result(spectrum)
    !!{
    Evaluate AGN spectra using a {\normalfont \ttfamily accretionDiskSpectra} object passed to us as an unlimited polymorphic object.
    !!}
    use :: Error                 , only : Error_Report
    use :: Accretion_Disk_Spectra, only : accretionDiskSpectraClass
    use :: Galacticus_Nodes      , only : treeNode
    implicit none
    class           (*       ), intent(inout) :: accretionDiskSpectra_
    type            (treeNode), intent(inout) :: node
    double precision          , intent(in   ) :: wavelength
    
    select type (accretionDiskSpectra_)
    class is (accretionDiskSpectraClass)
       spectrum=accretionDiskSpectra_%spectrum(node,wavelength)
    class default
       spectrum=0.0d0
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end function accretionDiskSpectraSpectrumNode
  
  !![
  <functionGlobal>
   <unitName>accretionDiskSpectraDestruct</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout), pointer :: accretionDiskSpectra_</arguments>
  </functionGlobal>
  !!]
  subroutine accretionDiskSpectraDestruct(accretionDiskSpectra_)
    !!{
    Destruct a {\normalfont \ttfamily accretionDiskSpectra} object passed to us as an unlimited polymorphic object.
    !!}
    use :: Error                 , only : Error_Report
    use :: Accretion_Disk_Spectra, only : accretionDiskSpectraClass
    implicit none
    class(*), intent(inout), pointer :: accretionDiskSpectra_

    accretionDiskSpectra__ => accretionDiskSpectra_
    select type (accretionDiskSpectra__)
    class is (accretionDiskSpectraClass)
       !![
       <objectDestructor name="accretionDiskSpectra__"/>
       !!]
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end subroutine accretionDiskSpectraDestruct

  !![
  <functionGlobal>
   <unitName>accretionDiskSpectraStateRestore</unitName>
   <type>void</type>
   <module>ISO_C_Binding, only : c_ptr, c_size_t</module>
   <arguments>class  (*       ), intent(inout) :: self            </arguments>
   <arguments>integer          , intent(in   ) :: stateFile       </arguments>
   <arguments>type   (c_ptr   ), intent(in   ) :: gslStateFile    </arguments>
   <arguments>integer(c_size_t), intent(in   ) :: stateOperationID</arguments>
   </functionGlobal>
  !!]
  subroutine accretionDiskSpectraStateRestore(self,stateFile,gslStateFile,stateOperationID)
    !!{
    Perform a deep copy of galactic structure objects.
    !!}
    use, intrinsic :: ISO_C_Binding         , only : c_ptr                    , c_size_t
    use            :: Error                 , only : Error_Report
    use            :: Accretion_Disk_Spectra, only : accretionDiskSpectraClass
    implicit none
    class  (*       ), intent(inout) :: self
    integer          , intent(in   ) :: stateFile
    type   (c_ptr   ), intent(in   ) :: gslStateFile
    integer(c_size_t), intent(in   ) :: stateOperationID

    select type (self)
    class is (accretionDiskSpectraClass)
       call self%stateRestore(stateFile,gslStateFile,stateOperationID)
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine accretionDiskSpectraStateRestore

  !![
  <functionGlobal>
   <unitName>accretionDiskSpectraStateStore</unitName>
   <type>void</type>
   <module>ISO_C_Binding, only : c_ptr, c_size_t</module>
   <arguments>class  (*       ), intent(inout) :: self            </arguments>
   <arguments>integer          , intent(in   ) :: stateFile       </arguments>
   <arguments>type   (c_ptr   ), intent(in   ) :: gslStateFile    </arguments>
   <arguments>integer(c_size_t), intent(in   ) :: stateOperationID</arguments>
  </functionGlobal>
  !!]
  subroutine accretionDiskSpectraStateStore(self,stateFile,gslStateFile,stateOperationID)
    !!{
    Perform a deep copy of galactic structure objects.
    !!}
    use, intrinsic :: ISO_C_Binding         , only : c_ptr                    , c_size_t
    use            :: Error                 , only : Error_Report
    use            :: Accretion_Disk_Spectra, only : accretionDiskSpectraClass
    implicit none
    class  (*       ), intent(inout) :: self
    integer          , intent(in   ) :: stateFile
    type   (c_ptr   ), intent(in   ) :: gslStateFile
    integer(c_size_t), intent(in   ) :: stateOperationID

    select type (self)
    class is (accretionDiskSpectraClass)
       call self%stateStore(stateFile,gslStateFile,stateOperationID)
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine accretionDiskSpectraStateStore

  !![
  <functionGlobal>
   <unitName>accretionDiskSpectraDeepCopyReset</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self</arguments>
  </functionGlobal>
  !!]
  subroutine accretionDiskSpectraDeepCopyReset(self)
    !!{
    Perform a deep copy of galactic structure objects.
    !!}
    use :: Error                 , only : Error_Report
    use :: Accretion_Disk_Spectra, only : accretionDiskSpectraClass
    implicit none
    class(*), intent(inout) :: self
    
    select type (self)
    class is (accretionDiskSpectraClass)
       call self%deepCopyReset()
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine accretionDiskSpectraDeepCopyReset
  
  !![
  <functionGlobal>
   <unitName>accretionDiskSpectraDeepCopyFinalize</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self</arguments>
  </functionGlobal>
  !!]
  subroutine accretionDiskSpectraDeepCopyFinalize(self)
    !!{
    Finalize a deep copy of galactic structure objects.
    !!}
    use :: Error                 , only : Error_Report
    use :: Accretion_Disk_Spectra, only : accretionDiskSpectraClass
    implicit none
    class(*), intent(inout) :: self
    
    select type (self)
    class is (accretionDiskSpectraClass)
       call self%deepCopyFinalize()
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine accretionDiskSpectraDeepCopyFinalize
  
  !![
  <functionGlobal>
   <unitName>accretionDiskSpectraDeepCopy</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self, destination</arguments>
  </functionGlobal>
  !!]
  subroutine accretionDiskSpectraDeepCopy(self,destination)
    !!{
    Perform a deep copy of galactic structure objects.
    !!}
    use :: Error                 , only : Error_Report
    use :: Accretion_Disk_Spectra, only : accretionDiskSpectraClass
    implicit none
    class(*), intent(inout) :: self, destination

    select type (self)
    class is (accretionDiskSpectraClass)
       select type (destination)
       class is (accretionDiskSpectraClass)
          call self%deepCopy(destination)
       class default
          call Error_Report("unexpected class"//{introspection:location})
       end select
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine accretionDiskSpectraDeepCopy

end module Accretion_Disk_Spectra_Utilities
