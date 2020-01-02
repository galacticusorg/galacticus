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

  use :: Kind_Numbers, only : kind_int8

  !# <nodePropertyExtractor name="nodePropertyExtractorIntegerTuple" abstract="yes">
  !#  <description>An abstract output analysis property extractor class which provides a tuple of integer properties.</description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorClass), abstract :: nodePropertyExtractorIntegerTuple
     !% A integerTuple property extractor.
     private
   contains
     !@ <objectMethods>
     !@  <object>nodePropertyExtractorIntegerTuple</object>
     !@  <objectMethod>
     !@   <method>elementCount</method>
     !@   <description>Return the number of properties in the tuple.</description>
     !@   <type>\intzero</type>
     !@   <pass>yes</pass>
     !@   <arguments>\doublezero\ time\argin</arguments>
     !@  </objectMethod>
     !@  <objectMethod>
     !@   <method>extract</method>
     !@   <description>Extract the properties from the given {\normalfont \ttfamily node}.</description>
     !@   <type>\intone</type>
     !@   <pass>yes</pass>
     !@   <arguments>\textcolor{red}{\textless type(treeNode)\textgreater} node\argin, \doublezero\ time\argin, \textcolor{red}{\textless type(multiCounter)\textgreater} [instance]\argin</arguments>
     !@  </objectMethod>
     !@  <objectMethod>
     !@   <method>names</method>
     !@   <description>Return the names of the properties extracted.</description>
     !@   <type>\textcolor{red}{\textless type(varying\_string)\textgreater}(:)</type>
     !@   <pass>yes</pass>
     !@   <arguments>\doublezero\ time\argin</arguments>
     !@  </objectMethod>
     !@  <objectMethod>
     !@   <method>descriptions</method>
     !@   <description>Return descriptions of the properties extracted.</description>
     !@   <type>\textcolor{red}{\textless type(varying\_string)\textgreater}(:)</type>
     !@   <pass>yes</pass>
     !@   <arguments>\doublezero\ time\argin</arguments>
     !@  </objectMethod>
     !@  <objectMethod>
     !@   <method>unitsInSI</method>
     !@   <description>Return the units of the properties extracted in the SI system.</description>
     !@   <type>\doubleone</type>
     !@   <pass>yes</pass>
     !@   <arguments>\doublezero\ time\argin</arguments>
     !@  </objectMethod>
     !@ </objectMethods>
     procedure(integerTupleElementCount), deferred :: elementCount
     procedure(integerTupleExtract     ), deferred :: extract
     procedure(integerTupleNames       ), deferred :: names
     procedure(integerTupleNames       ), deferred :: descriptions
     procedure(integerTupleUnitsInSI   ), deferred :: unitsInSI
  end type nodePropertyExtractorIntegerTuple

  abstract interface
     function integerTupleExtract(self,node,time,instance)
       !% Interface for {\normalfont \ttfamily integerTuple} property extraction.
       import nodePropertyExtractorIntegerTuple, treeNode, multiCounter, kind_int8
       integer         (kind_int8                        ), dimension(:) , allocatable :: integerTupleExtract
       class           (nodePropertyExtractorIntegerTuple), intent(inout)              :: self
       type            (treeNode                         ), intent(inout)              :: node
       double precision                                   , intent(in   )              :: time
       type            (multiCounter                     ), intent(inout), optional    :: instance
     end function integerTupleExtract
  end interface

  abstract interface
     function integerTupleNames(self,time)
       !% Interface for {\normalfont \ttfamily integerTuple} property names.
       import varying_string, nodePropertyExtractorIntegerTuple
       type            (varying_string                   ), dimension(:) , allocatable :: integerTupleNames
       class           (nodePropertyExtractorIntegerTuple), intent(inout)              :: self
       double precision                                   , intent(in   )              :: time
     end function integerTupleNames
  end interface

  abstract interface
     function integerTupleUnitsInSI(self,time)
       !% Interface for {\normalfont \ttfamily integerTuple property units.
       import nodePropertyExtractorIntegerTuple
       double precision                                    , dimension(:) , allocatable :: integerTupleUnitsInSI
       class           (nodePropertyExtractorIntegerTuple), intent(inout)              :: self
       double precision                                   , intent(in   )              :: time
     end function integerTupleUnitsInSI
  end interface

  abstract interface
     integer function integerTupleElementCount(self,time)
       !% Interface for {\normalfont \ttfamily integerTuple} element count.
       import nodePropertyExtractorIntegerTuple
       class           (nodePropertyExtractorIntegerTuple), intent(inout) :: self
       double precision                                   , intent(in   ) :: time
     end function integerTupleElementCount
  end interface
