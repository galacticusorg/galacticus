!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

  !# <nodePropertyExtractor name="nodePropertyExtractorScalar" abstract="yes">
  !#  <description>An abstract output analysis property extractor class which provieds a scalar floating point property.</description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorClass), abstract :: nodePropertyExtractorScalar
     !% A scalar output analysis class.
     private
   contains
     !# <methods>
     !#   <method description="Extract the property from the given {\normalfont \ttfamily node}." method="extract" pass="yes" />
     !#   <method description="Return the name of the property extracted." method="name" pass="yes" />
     !#   <method description="Return a description of the property extracted." method="description" pass="yes" />
     !#   <method description="Return the units of the property extracted in the SI system." method="unitsInSI" pass="yes" />
     !# </methods>
     procedure(scalarExtract  ), deferred :: extract
     procedure(scalarName     ), deferred :: name
     procedure(scalarName     ), deferred :: description
     procedure(scalarUnitsInSI), deferred :: unitsInSI
  end type nodePropertyExtractorScalar

  abstract interface
     double precision function scalarExtract(self,node,instance)
       !% Interface for scalar property extraction.
       import nodePropertyExtractorScalar, treeNode, multiCounter
       class(nodePropertyExtractorScalar), intent(inout)           :: self
       type (treeNode                   ), intent(inout), target   :: node
       type (multiCounter               ), intent(inout), optional :: instance
     end function scalarExtract
  end interface

  abstract interface
     function scalarName(self)
       !% Interface for scalar property name.
       import varying_string, nodePropertyExtractorScalar
       type (varying_string             )                :: scalarName
       class(nodePropertyExtractorScalar), intent(inout) :: self
     end function scalarName
  end interface

  abstract interface
     double precision function scalarUnitsInSI(self)
       !% Interface for scalar property units.
       import nodePropertyExtractorScalar
       class(nodePropertyExtractorScalar), intent(inout) :: self
     end function scalarUnitsInSI
  end interface
