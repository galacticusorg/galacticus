!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module that implements a class of parameter mapping functions.

module Constraints_Mappings
  !% Implements a class of parameter mapping functions.
  use Constraints_Priors
  private
  public :: mappingNew

  ! Define the basic mapping type.
  type, abstract, public :: mapping
   contains
     !@ <objectMethods>
     !@   <object>mapping</object>
     !@   <objectMethod>
     !@     <method>mapPrior</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless class(prior)\textgreater} thisPrior\arginout</arguments>
     !@     <description>Attempt to apply the mapping to a prior.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>unmap</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero x\argin</arguments>
     !@     <description>Unmap the given {\tt x}.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure(mappingMapPrior), deferred :: mapPrior
     procedure(mappingUnmap   ), deferred :: unmap
  end type mapping
  
  abstract interface
     subroutine mappingMapPrior(self,thisPrior)
       import :: mapping, prior
       class(mapping), intent(in   ) :: self
       class(prior  ), intent(inout) :: thisPrior
     end subroutine mappingMapPrior
  end interface
  
  abstract interface
     double precision function mappingUnmap(self,x)
       import :: mapping
       class           (mapping), intent(in   ) :: self
       double precision         , intent(in   ) :: x
     end function mappingUnmap
  end interface
  
  ! Define a list of mappings.
  type, public :: mappingList
     class(mapping), pointer :: thisMapping
  end type mappingList

  ! Include all distribution types.
  include 'constraints.mappings.linear.type.inc'
  include 'constraints.mappings.logarithmic.type.inc'

contains

  function mappingNew(definition) result (newMapping)
    !% Create a new mapping from an XML definition.
    use FoX_DOM
    use IO_XML
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    class           (mapping), pointer                :: newMapping
    type            (node   ), pointer, intent(in   ) :: definition
  
    select case (char(XML_Extract_Text(XML_Get_First_Element_By_Tag_Name(definition,"type"))))
    case ("linear")
       allocate(mappingLinear :: newMapping)
       select type (newMapping)
       type is (mappingLinear)
          newMapping=mappingLinear()
       end select
    case ("logarithmic")
       allocate(mappingLogarithmic :: newMapping)
       select type (newMapping)
       type is (mappingLogarithmic)
          newMapping=mappingLogarithmic()
       end select
    case default
       call Galacticus_Error_Report('mappingNew','mapping type is unrecognized')
    end select
    return
  end function mappingNew

  ! Include all mapping methods.
  include 'constraints.mappings.linear.methods.inc'
  include 'constraints.mappings.logarithmic.methods.inc'
  
end module Constraints_Mappings
