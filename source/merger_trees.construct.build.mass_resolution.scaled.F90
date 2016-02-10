!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of a merger tree builder mass resolution which assumes a resolution which scales with tree mass.

  !# <mergerTreeMassResolution name="mergerTreeMassResolutionScaled">
  !#  <description>Provides a mass resolution for merger tree building which scales with tree mass.</description>
  !# </mergerTreeMassResolution>
  type, extends(mergerTreeMassResolutionClass) :: mergerTreeMassResolutionScaled
     !% A merger tree mass resolution class which assumes a mass resolution that scales with tree mass.
     private
     double precision :: massResolutionMinimum, massResolutionFractional
   contains
     procedure :: resolution => scaledResolution
  end type mergerTreeMassResolutionScaled

  interface mergerTreeMassResolutionScaled
     !% Constructors for the {\normalfont \ttfamily scaled} merger tree resolution class.
     module procedure scaledConstructorParameters
     module procedure scaledConstructorInternal
  end interface mergerTreeMassResolutionScaled
  
contains

  function scaledConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily scaled} merger tree building mass resolution class which reads parameters from a
    !% provided parameter list.
    implicit none
    type(mergerTreeMassResolutionScaled)                :: scaledConstructorParameters    
    type(inputParameters               ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>massResolutionMinimum</name>
    !#   <source>parameters</source>
    !#   <variable>scaledConstructorParameters%massResolutionMinimum</variable>
    !#   <defaultValue>5.0d9</defaultValue>
    !#   <description>The minimum mass resolution to use when building merger trees.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massResolutionFractional</name>
    !#   <source>parameters</source>
    !#   <variable>scaledConstructorParameters%massResolutionFractional</variable>
    !#   <defaultValue>1.0d-3</defaultValue>
    !#   <description>The fraction of the tree's root node mass to be used for the mass resolution when building merger trees.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    return
  end function scaledConstructorParameters

  function scaledConstructorInternal(massResolutionMinimum,massResolutionFractional)
    !% Internal constructor for the {\normalfont \ttfamily scaled} merger tree building mass resolution class.
    implicit none
    type            (mergerTreeMassResolutionScaled)                :: scaledConstructorInternal
    double precision                                , intent(in   ) :: massResolutionMinimum    , massResolutionFractional

    scaledConstructorInternal%massResolutionMinimum   =massResolutionMinimum
    scaledConstructorInternal%massResolutionFractional=massResolutionFractional
    return
  end function scaledConstructorInternal

  double precision function scaledResolution(self,tree)
    !% Returns a scaled mass resolution to use when building merger trees.
    implicit none
    class(mergerTreeMassResolutionScaled), intent(inout) :: self
    type (mergerTree                    ), intent(in   ) :: tree
    class(nodeComponentBasic            ), pointer       :: baseBasic

    baseBasic        => tree%baseNode%basic()
    scaledResolution =  max(                                                &
         &                  self%massResolutionMinimum                    , &
         &                  self%massResolutionFractional*baseBasic%mass()  &
         &                 )
    return
  end function scaledResolution
