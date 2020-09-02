!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implements a node operator class that outputs data on mergers between galaxies.

  use :: IO_HDF5         , only : hdf5Object
  use :: Galactic_Filters, only : galacticFilterClass
  
  !# <nodeOperator name="nodeOperatorGalaxyMergersOutput">
  !#  <description>A node operator class that outputs data on mergers between galaxies.</description>
  !# </nodeOperator>
  type, extends(nodeOperatorClass) :: nodeOperatorGalaxyMergersOutput
     !% A node operator class that shifts node indices at node promotion.
     private
     type (varying_string     )          :: mergersGroupName
     type (hdf5Object         )          :: mergersGroup
     class(galacticFilterClass), pointer :: galacticFilterSatellite_ => null(), galacticFilterCentral_ => null()
   contains
     final     ::                  galaxyMergersDestructor
     procedure :: galaxiesMerge => galaxyMergersOutputGalaxiesMerge
  end type nodeOperatorGalaxyMergersOutput
  
  interface nodeOperatorGalaxyMergersOutput
     !% Constructors for the {\normalfont \ttfamily galaxyMergersOutput} node operator class.
     module procedure galaxyMergersOutputConstructorParameters
     module procedure galaxyMergersOutputConstructorInternal
  end interface nodeOperatorGalaxyMergersOutput
  
contains
  
  function galaxyMergersOutputConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily galaxyMergersOutput} node operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorGalaxyMergersOutput)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(galacticFilterClass            ), pointer       :: galacticFilterSatellite_, galacticFilterCentral_
    type (varying_string                 )                :: mergersGroupName

    !# <inputParameter>
    !#   <name>mergersGroupName</name>
    !#   <defaultValue>var_str('galaxyMergers')</defaultValue>
    !#   <description>The name of the HDF5 group to which galaxy merger data should be output.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="galacticFilter" parameterName="galacticFilterSatellite" name="galacticFilterSatellite_" source="parameters">
    !#  <default>
    !#   <galacticFilterSatellite value="always"/>
    !#  </default>
    !# </objectBuilder>
    !# <objectBuilder class="galacticFilter" parameterName="galacticFilterCentral"   name="galacticFilterCentral_"   source="parameters">
    !#  <default>
    !#   <galacticFilterCentral value="always"/>
    !#  </default>
    !# </objectBuilder>
    self=nodeOperatorGalaxyMergersOutput(mergersGroupName,galacticFilterSatellite_,galacticFilterCentral_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="galacticFilterSatellite_"/>
    !# <objectDestructor name="galacticFilterCentral_"  />
    return
  end function galaxyMergersOutputConstructorParameters

  function galaxyMergersOutputConstructorInternal(mergersGroupName,galacticFilterSatellite_,galacticFilterCentral_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily galaxyMergersOutput} node operator class.
    implicit none
    type (nodeOperatorGalaxyMergersOutput)                        :: self
    class(galacticFilterClass            ), intent(in   ), target :: galacticFilterSatellite_, galacticFilterCentral_
    type (varying_string                 ), intent(in   )         :: mergersGroupName
    !# <constructorAssign variables="mergersGroupName, *galacticFilterSatellite_, *galacticFilterCentral_"/>

    return
  end function galaxyMergersOutputConstructorInternal
  
  subroutine galaxyMergersDestructor(self)
    !% Destructor for the {\normalfont \ttfamily galaxyMergersOutput} node operator class.
    !$ use :: IO_HDF5         , only : hdf5Access
    implicit none
    type(nodeOperatorGalaxyMergersOutput), intent(inout) :: self

    !# <objectDestructor name="self%galacticFilterSatellite_"/>
    !# <objectDestructor name="self%galacticFilterCentral_" />
    !$ call hdf5Access%set  ()
    if (self%mergersGroup%isOpen()) call self%mergersGroup%close()
    !$ call hdf5Access%unset()
    return
  end subroutine galaxyMergersDestructor
  
  subroutine galaxyMergersOutputGalaxiesMerge(self,node)
    !% Act on a merger between galaxies.
    use    :: Galacticus_HDF5 , only : galacticusOutputFile
    use    :: Galacticus_Nodes, only : nodeComponentBasic
    !$ use :: IO_HDF5         , only : hdf5Access
    implicit none
    class(nodeOperatorGalaxyMergersOutput), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node
    type (treeNode                       ), pointer       :: nodeHost
    class(nodeComponentBasic             ), pointer       :: basic   , basicHost

    nodeHost  => node    %mergesWith()
    basic     => node    %basic     ()
    basicHost => nodeHost%basic     ()
    if (.not.(self%galacticFilterSatellite_%passes(node).and.self%galacticFilterCentral_%passes(nodeHost))) return
    !$ call hdf5Access%set  ()
    if (.not.self%mergersGroup%isOpen()) self%mergersGroup=galacticusOutputFile%openGroup(char(self%mergersGroupName),"Satellite mergers data.")
    call self%mergersGroup%writeDataset([node              %index()],"indexSatellite","Index of the satellite."               ,appendTo=.true.)
    call self%mergersGroup%writeDataset([nodeHost          %index()],"indexCentral"  ,"Index of the central."                 ,appendTo=.true.)
    call self%mergersGroup%writeDataset([node     %hostTree%index  ],"indexTree"     ,"Index of the tree."                    ,appendTo=.true.)
    call self%mergersGroup%writeDataset([basic             %time ()],"time"          ,"Time at which the merger occurs [Gyr].",appendTo=.true.)
    call self%mergersGroup%writeDataset([basic             %mass ()],"massSatellite" ,"Mass of the satellite halo [M☉]."      ,appendTo=.true.)
    call self%mergersGroup%writeDataset([basicHost         %mass ()],"massCentral"   ,"Mass of the central halo [M☉]."        ,appendTo=.true.)
    !$ call hdf5Access%unset()
    return
  end subroutine galaxyMergersOutputGalaxiesMerge
