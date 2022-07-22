!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module of globally-accessible functions supporting the galactic structure class.
!!}

module Galactic_Structure_Utilities
  !!{
  Provides globally-accessible functions supporting the galactic structure class.
  !!}
  private
  public :: galacticStructureConstruct       , galacticStructureMassEnclosed            , galacticStructureDestruct  , galacticStructureDeepCopy    , &
       &    galacticStructureDeepCopyReset   , galacticStructureDeepCopyFinalize        , galacticStructureStateStore, galacticStructureStateRestore, &
       &    galacticStructureVelocityRotation, galacticStructureVelocityRotationGradient, galacticStructureDensity

  ! Module-scope pointer to our task object. This is used for reference counting so that debugging information is consistent
  ! between the increments and decrements.
  class(*), pointer :: galacticStructure__
  !$omp threadprivate(galacticStructure__)

contains

  !![
  <functionGlobal>
   <unitName>galacticStructureConstruct</unitName>
   <type>void</type>
   <module>Input_Parameters, only : inputParameters</module>
   <arguments>type (inputParameters), intent(inout)          :: parameters</arguments>
   <arguments>class(*              ), intent(  out), pointer :: galacticStructure_</arguments>
  </functionGlobal>
  !!]
  subroutine galacticStructureConstruct(parameters,galacticStructure_)
    !!{
    Build a {\normalfont \ttfamily galacticStructure} object from a given parameter set. This is a globally-callable function
    to allow us to subvert the class/module hierarchy.
    !!}
    use :: Error             , only : Error_Report
    use :: Input_Parameters  , only : inputParameter        , inputParameters
    use :: Galactic_Structure, only : galacticStructureClass, galacticStructure
    implicit none
    type (inputParameters), intent(inout)          :: parameters
    class(*              ), intent(  out), pointer :: galacticStructure_

    galacticStructure__ => galacticStructure(parameters)
    select type (galacticStructure__)
    class is (galacticStructureClass)
       !![
       <referenceCountIncrement object="galacticStructure__"/>
       !!]
       call galacticStructure__%autoHook()
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select
    galacticStructure_ => galacticStructure__
    return
  end subroutine galacticStructureConstruct

  !![
  <functionGlobal>
   <unitName>galacticStructureMassEnclosed</unitName>
   <type>double precision</type>
   <module>Galacticus_Nodes          , only : treeNode</module>
   <module>Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType</module>
   <arguments>class           (*                           ), intent(inout)           :: galacticStructure_</arguments>
   <arguments>type            (treeNode                    ), intent(inout)           :: node</arguments>
   <arguments>double precision                              , intent(in   ), optional :: radius</arguments>
   <arguments>type            (enumerationComponentTypeType), intent(in   ), optional :: componentType</arguments>
   <arguments>type            (enumerationMassTypeType     ), intent(in   ), optional :: massType</arguments>
   <arguments>type            (enumerationWeightByType     ), intent(in   ), optional :: weightBy</arguments>
   <arguments>integer                                       , intent(in   ), optional :: weightIndex</arguments>
  </functionGlobal>
  !!]
  double precision function galacticStructureMassEnclosed(galacticStructure_,node,radius,componentType,massType,weightBy,weightIndex)
    !!{
    Compute the mass enclosed for a {\normalfont \ttfamily galacticStructure} object passed to us as an unlimited polymorphic object.
    !!}
    use :: Error                     , only : Error_Report
    use :: Galactic_Structure        , only : galacticStructureClass
    use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
    use :: Galacticus_Nodes          , only : treeNode
    implicit none
    class           (*                           ), intent(inout)           :: galacticStructure_
    type            (treeNode                    ), intent(inout)           :: node
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (enumerationWeightByType     ), intent(in   ), optional :: weightBy
    integer                                       , intent(in   ), optional :: weightIndex
    double precision                              , intent(in   ), optional :: radius
 
    select type (galacticStructure_)
    class is (galacticStructureClass)
       galacticStructureMassEnclosed=galacticStructure_%massEnclosed(node,radius,componentType,massType,weightBy,weightIndex)
    class default
       galacticStructureMassEnclosed=0.0d0
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end function galacticStructureMassEnclosed

  !![
  <functionGlobal>
    <unitName>galacticStructureDensity</unitName>
    <type>double precision</type>
    <module>Galacticus_Nodes          , only : treeNode</module>
    <module>Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType, enumerationCoordinateSystemType</module>
    <arguments>class           (*                              ), intent(inout)               :: galacticStructure_</arguments>
    <arguments>type            (treeNode                       ), intent(inout)               :: node</arguments>
    <arguments>double precision                                 , intent(in   ), dimension(3) :: position</arguments>
    <arguments>type            (enumerationCoordinateSystemType), intent(in   ), optional     :: coordinateSystem</arguments>
    <arguments>type            (enumerationComponentTypeType   ), intent(in   ), optional     :: componentType</arguments>
    <arguments>type            (enumerationMassTypeType        ), intent(in   ), optional     :: massType</arguments>
    <arguments>type            (enumerationWeightByType        ), intent(in   ), optional     :: weightBy</arguments>
    <arguments>integer                                          , intent(in   ), optional     :: weightIndex</arguments>
  </functionGlobal>
  !!]
  double precision function galacticStructureDensity(galacticStructure_,node,position,coordinateSystem,componentType,massType,weightBy,weightIndex)
    !!{
    Compute the density for a {\normalfont \ttfamily galacticStructure} object passed to us as an unlimited polymorphic object.
    !!}
    use :: Error                     , only : Error_Report
    use :: Galactic_Structure        , only : galacticStructureClass
    use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType, enumerationCoordinateSystemType
    use :: Galacticus_Nodes          , only : treeNode
    implicit none
    class           (*                              ), intent(inout)               :: galacticStructure_
    type            (treeNode                       ), intent(inout)               :: node
    type            (enumerationComponentTypeType   ), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional     :: massType
    type            (enumerationWeightByType        ), intent(in   ), optional     :: weightBy
    type            (enumerationCoordinateSystemType), intent(in   ), optional     :: coordinateSystem
    integer                                          , intent(in   ), optional     :: weightIndex
    double precision                                 , intent(in   ), dimension(3) :: position
 
    select type (galacticStructure_)
    class is (galacticStructureClass)
       galacticStructureDensity=galacticStructure_%density(node,position,coordinateSystem,componentType,massType,weightBy,weightIndex)
    class default
       galacticStructureDensity=0.0d0
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end function galacticStructureDensity
  
  !![
  <functionGlobal>
   <unitName>galacticStructureVelocityRotation</unitName>
   <type>double precision</type>
   <module>Galacticus_Nodes          , only : treeNode</module>
   <module>Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType</module>
   <arguments>class           (*                           ), intent(inout)           :: galacticStructure_</arguments>
   <arguments>type            (treeNode                    ), intent(inout)           :: node</arguments>
   <arguments>double precision                              , intent(in   ), optional :: radius</arguments>
   <arguments>type            (enumerationComponentTypeType), intent(in   ), optional :: componentType</arguments>
   <arguments>type            (enumerationMassTypeType     ), intent(in   ), optional :: massType</arguments>
  </functionGlobal>
  !!]
  double precision function galacticStructureVelocityRotation(galacticStructure_,node,radius,componentType,massType)
    !!{
    Compute the rotation velocity for a {\normalfont \ttfamily galacticStructure} object passed to us as an unlimited polymorphic object.
    !!}
    use :: Error                     , only : Error_Report
    use :: Galactic_Structure        , only : galacticStructureClass
    use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType
    use :: Galacticus_Nodes          , only : treeNode
    implicit none
    class           (*                           ), intent(inout)           :: galacticStructure_
    type            (treeNode                    ), intent(inout)           :: node
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                              , intent(in   ), optional :: radius
 
    select type (galacticStructure_)
    class is (galacticStructureClass)
       galacticStructureVelocityRotation=galacticStructure_%velocityRotation(node,radius,componentType,massType)
    class default
       galacticStructureVelocityRotation=0.0d0
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end function galacticStructureVelocityRotation

  !![
  <functionGlobal>
   <unitName>galacticStructureVelocityRotationGradient</unitName>
   <type>double precision</type>
   <module>Galacticus_Nodes          , only : treeNode</module>
   <module>Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType</module>
   <arguments>class           (*                           ), intent(inout)           :: galacticStructure_</arguments>
   <arguments>type            (treeNode                    ), intent(inout)           :: node</arguments>
   <arguments>double precision                              , intent(in   ), optional :: radius</arguments>
   <arguments>type            (enumerationComponentTypeType), intent(in   ), optional :: componentType</arguments>
   <arguments>type            (enumerationMassTypeType     ), intent(in   ), optional :: massType</arguments>
  </functionGlobal>
  !!]
  double precision function galacticStructureVelocityRotationGradient(galacticStructure_,node,radius,componentType,massType)
    !!{
    Compute the rotation velocity gradient for a {\normalfont \ttfamily galacticStructure} object passed to us as an unlimited polymorphic object.
    !!}
    use :: Error                     , only : Error_Report
    use :: Galactic_Structure        , only : galacticStructureClass
    use :: Galacticus_Nodes          , only : treeNode
    use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType
    implicit none
    class           (*                           ), intent(inout)           :: galacticStructure_
    type            (treeNode                    ), intent(inout)           :: node
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                              , intent(in   ), optional :: radius

    select type (galacticStructure_)
    class is (galacticStructureClass)
       galacticStructureVelocityRotationGradient=galacticStructure_%velocityRotationGradient(node,radius,componentType,massType)
    class default
       galacticStructureVelocityRotationGradient=0.0d0
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end function galacticStructureVelocityRotationGradient
  
  !![
  <functionGlobal>
   <unitName>galacticStructureDestruct</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout), pointer :: galacticStructure_</arguments>
  </functionGlobal>
  !!]
  subroutine galacticStructureDestruct(galacticStructure_)
    !!{
    Destruct a {\normalfont \ttfamily taskEvolveForests} object passed to us as an unlimited polymorphic object.
    !!}
    use :: Error             , only : Error_Report
    use :: Galactic_Structure, only : galacticStructureClass
    implicit none
    class(*), intent(inout), pointer :: galacticStructure_

    galacticStructure__ => galacticStructure_
    select type (galacticStructure__)
    class is (galacticStructureClass)
       !![
       <objectDestructor name="galacticStructure__"/>
       !!]
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end subroutine galacticStructureDestruct

  !![
  <functionGlobal>
   <unitName>galacticStructureStateRestore</unitName>
   <type>void</type>
   <module>ISO_C_Binding, only : c_ptr, c_size_t</module>
   <arguments>class  (*       ), intent(inout) :: self</arguments>
   <arguments>integer          , intent(in   ) :: stateFile</arguments>
   <arguments>type   (c_ptr   ), intent(in   ) :: gslStateFile</arguments>
   <arguments>integer(c_size_t), intent(in   ) :: stateOperationID</arguments>
   </functionGlobal>
  !!]
  subroutine galacticStructureStateRestore(self,stateFile,gslStateFile,stateOperationID)
    !!{
    Perform a deep copy of galactic structure objects.
    !!}
    use, intrinsic :: ISO_C_Binding     , only : c_ptr                 , c_size_t
    use            :: Error             , only : Error_Report
    use            :: Galactic_Structure, only : galacticStructureClass
    implicit none
    class  (*       ), intent(inout) :: self
    integer          , intent(in   ) :: stateFile
    type   (c_ptr   ), intent(in   ) :: gslStateFile
    integer(c_size_t), intent(in   ) :: stateOperationID

    select type (self)
    class is (galacticStructureClass)
       call self%stateRestore(stateFile,gslStateFile,stateOperationID)
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine galacticStructureStateRestore

  !![
  <functionGlobal>
   <unitName>galacticStructureStateStore</unitName>
   <type>void</type>
   <module>ISO_C_Binding, only : c_ptr, c_size_t</module>
   <arguments>class  (*       ), intent(inout) :: self</arguments>
   <arguments>integer          , intent(in   ) :: stateFile</arguments>
   <arguments>type   (c_ptr   ), intent(in   ) :: gslStateFile</arguments>
   <arguments>integer(c_size_t), intent(in   ) :: stateOperationID</arguments>
  </functionGlobal>
  !!]
  subroutine galacticStructureStateStore(self,stateFile,gslStateFile,stateOperationID)
    !!{
    Perform a deep copy of galactic structure objects.
    !!}
    use, intrinsic :: ISO_C_Binding     , only : c_ptr                  , c_size_t
    use            :: Error             , only : Error_Report
    use            :: Galactic_Structure, only : galacticStructureClass
    implicit none
    class  (*       ), intent(inout) :: self
    integer          , intent(in   ) :: stateFile
    type   (c_ptr   ), intent(in   ) :: gslStateFile
    integer(c_size_t), intent(in   ) :: stateOperationID

    select type (self)
    class is (galacticStructureClass)
       call self%stateStore(stateFile,gslStateFile,stateOperationID)
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine galacticStructureStateStore

  !![
  <functionGlobal>
   <unitName>galacticStructureDeepCopyReset</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self</arguments>
  </functionGlobal>
  !!]
  subroutine galacticStructureDeepCopyReset(self)
    !!{
    Perform a deep copy of galactic structure objects.
    !!}
    use :: Error             , only : Error_Report
    use :: Galactic_Structure, only : galacticStructureClass
    implicit none
    class(*), intent(inout) :: self
    
    select type (self)
    class is (galacticStructureClass)
       call self%deepCopyReset()
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine galacticStructureDeepCopyReset
  
  !![
  <functionGlobal>
   <unitName>galacticStructureDeepCopyFinalize</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self</arguments>
  </functionGlobal>
  !!]
  subroutine galacticStructureDeepCopyFinalize(self)
    !!{
    Finalize a deep copy of galactic structure objects.
    !!}
    use :: Error             , only : Error_Report
    use :: Galactic_Structure, only : galacticStructureClass
    implicit none
    class(*), intent(inout) :: self
    
    select type (self)
    class is (galacticStructureClass)
       call self%deepCopyFinalize()
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine galacticStructureDeepCopyFinalize
  
  !![
  <functionGlobal>
   <unitName>galacticStructureDeepCopy</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self, destination</arguments>
  </functionGlobal>
  !!]
  subroutine galacticStructureDeepCopy(self,destination)
    !!{
    Perform a deep copy of galactic structure objects.
    !!}
    use :: Error             , only : Error_Report
    use :: Galactic_Structure, only : galacticStructureClass
    implicit none
    class(*), intent(inout) :: self, destination

    select type (self)
    class is (galacticStructureClass)
       select type (destination)
       class is (galacticStructureClass)
          call self%deepCopy(destination)
       class default
          call Error_Report("unexpected class"//{introspection:location})
       end select
    class default
       call Error_Report("unexpected class"//{introspection:location})
    end select
    return
  end subroutine galacticStructureDeepCopy

end module Galactic_Structure_Utilities
