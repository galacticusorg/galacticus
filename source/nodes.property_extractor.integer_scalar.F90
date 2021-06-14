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

  use :: Kind_Numbers, only : kind_int8

  !# <nodePropertyExtractor name="nodePropertyExtractorIntegerScalar" abstract="yes">
  !#  <description>An abstract output analysis property extractor class which provieds a scalar integer property.</description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorClass), abstract :: nodePropertyExtractorIntegerScalar
     !% A scalar integer node property extractor class.
     private
   contains
     !# <methods>
     !#   <method description="Extract the property from the given {\normalfont \ttfamily node}." method="extract" pass="yes" />
     !#   <method description="Return the name of the property extracted." method="name" pass="yes" />
     !#   <method description="Return a description of the property extracted." method="description" pass="yes" />
     !#   <method description="Return the units of the property extracted in the SI system." method="unitsInSI" pass="yes" />
     !# </methods>
     procedure(integerScalarExtract), deferred :: extract
     procedure(integerScalarName   ), deferred :: name
     procedure(integerScalarName   ), deferred :: description
     procedure                                 :: unitsInSI   => integerScalarUnitsInSI
  end type nodePropertyExtractorIntegerScalar

  abstract interface
     function integerScalarExtract(self,node,time,instance)
       !% Interface for integerScalar property extraction.
       import nodePropertyExtractorIntegerScalar, treeNode, multiCounter, kind_int8
       integer         (kind_int8                         )                          :: integerScalarExtract
       class           (nodePropertyExtractorIntegerScalar), intent(inout)           :: self
       type            (treeNode                          ), intent(inout), target   :: node
       double precision                                    , intent(in   )           :: time
       type            (multiCounter                      ), intent(inout), optional :: instance
     end function integerScalarExtract
  end interface

  abstract interface
     function integerScalarName(self)
       !% Interface for integerScalar property name.
       import varying_string, nodePropertyExtractorIntegerScalar
       type (varying_string                    )                :: integerScalarName
       class(nodePropertyExtractorIntegerScalar), intent(inout) :: self
     end function integerScalarName
  end interface

contains

  double precision function integerScalarUnitsInSI(self)
    !% Interface for integerScalar property units.
    implicit none
    class(nodePropertyExtractorIntegerScalar), intent(inout) :: self
    !$GLC attributes unused :: self

    integerScalarUnitsInSI=0.0d0
    return
  end function integerScalarUnitsInSI
