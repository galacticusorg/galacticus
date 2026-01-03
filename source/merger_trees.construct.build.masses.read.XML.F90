!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  Implementation of a merger tree masses class which reads masses from an XML file.
  !!}

  !![
  <mergerTreeBuildMasses name="mergerTreeBuildMassesReadXML">
   <description>
    A merger tree build masses class which reads masses from an XML file. The XML file should have the following form:
    \begin{verbatim}
     &lt;mergerTrees>
      &lt;treeRootMass>13522377303.5998&lt;/treeRootMass>
      &lt;treeRootMass>19579530191.8709&lt;/treeRootMass>
      &lt;treeRootMass>21061025282.9613&lt;/treeRootMass>
      .
      .
      .
      &lt;treeWeight>13522377303.5998&lt;/treeWeight>
      &lt;treeWeight>19579530191.8709&lt;/treeWeight>
      &lt;treeWeight>21061025282.9613&lt;/treeWeight>
      .
      .
      .
     &lt;/mergerTrees>
    \end{verbatim}
    where each {\normalfont \ttfamily treeRootMass} element gives the mass (in Solar masses) of the root halo of a tree to
    generate, and the (optional) {\normalfont \ttfamily treeWeight} elements give the corresponding weight (in units of
    Mpc$^{-3}$) to assign to each tree.
   </description>
  </mergerTreeBuildMasses>
  !!]
  type, extends(mergerTreeBuildMassesRead) :: mergerTreeBuildMassesReadXML
     !!{
     Implementation of a merger tree masses class which reads masses from an XML file.
     !!}
     private
   contains
     procedure :: read => readXMLRead
  end type mergerTreeBuildMassesReadXML

  interface mergerTreeBuildMassesReadXML
     module procedure readXMLConstructorParameters
     module procedure readXMLConstructorInternal
  end interface mergerTreeBuildMassesReadXML

contains

  function readXMLConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeBuildMassesReadXML} merger tree masses class which takes a parameter set
    as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeBuildMassesReadXML)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    type            (varying_string              )                :: fileName
    double precision                                              :: massIntervalFractional

    !![
    <inputParameter>
      <name>fileName</name>
      <description>The name of the file from which to read merger tree masses.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massIntervalFractional</name>
      <defaultValue>0.1d0</defaultValue>
      <description>The fractional mass interval occupied by the trees. Where the intervals of trees of different mass would overlap this interval will be truncated.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=mergerTreeBuildMassesReadXML(fileName,massIntervalFractional)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function readXMLConstructorParameters

  function readXMLConstructorInternal(fileName,massIntervalFractional) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeBuildMassesReadXML} merger tree masses class.
    !!}
    implicit none
    type            (mergerTreeBuildMassesReadXML)                :: self
    type            (varying_string              ), intent(in   ) :: fileName
    double precision                              , intent(in   ) :: massIntervalFractional
    !![
    <constructorAssign variables="fileName, massIntervalFractional"/>
    !!]

    return
  end function readXMLConstructorInternal

  subroutine readXMLRead(self,mass,weight)
    !!{
    Read merger tree masses from file.
    !!}
    use :: FoX_DOM, only : destroy       , getDocumentElement   , node
    use :: Error  , only : Error_Report
    use :: IO_XML , only : XML_Array_Read, XML_Array_Read_Static, XML_Path_Exists, XML_Parse
    implicit none
    class           (mergerTreeBuildMassesReadXML), intent(inout)                            :: self
    double precision                              , intent(  out), allocatable, dimension(:) :: mass , weight
    type            (node                        ), pointer                                  :: doc  , rootNode
    integer                                                                                  :: ioErr

    !$omp critical (FoX_DOM_Access)
    doc => XML_Parse(self%fileName,iostat=ioErr)
    if (ioErr /= 0) call Error_Report('unable to read or parse merger tree root mass file'//{introspection:location})
    rootNode => getDocumentElement(doc)
    ! Read all tree masses.
    call XML_Array_Read(doc,"treeRootMass",mass)
    ! Extract tree weights if available.
    if (XML_Path_Exists(rootNode,"treeWeight")) then
       call XML_Array_Read(doc,"treeWeight",weight)
       if (size(weight) /= size(mass)) call Error_Report('`mass` and `weight` arrays must have the same size'//{introspection:location})
    end if
    ! Finished - destroy the XML document.
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    return
  end subroutine readXMLRead

