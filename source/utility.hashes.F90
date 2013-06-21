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

!% Contains a module which implements ``hashes'' (i.e. associative arrays).

module Hashes
  !% Implements ``hashes'' (i.e. associative arrays).
  use ISO_Varying_String
  implicit none
  private
  public :: integerScalarHash

  type integerScalarHash
     !% Derived type for integer hashes.
     private
     integer                              :: allocatedSize   , elementCount 
     integer                , allocatable :: hashValues   (:)               
     type   (varying_string), allocatable :: hashKeys     (:)               
   contains
     !@ <objectMethods>
     !@   <object>integerScalarHash</object>
     !@   <objectMethod>
     !@     <method>initialize</method>
     !@     <description>Initialize the hash.</description>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>set</method>
     !@     <description>Set the value of a key in the hash.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless (character(len=*)|varying\_string)\textgreater} key\argin, \intzero\ value\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>delete</method>
     !@     <description>Delete a key from the hash.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless (character(len=*)|varying\_string)\textgreater} key\argin, \intzero\ value\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>value</method>
     !@     <description>Return the value for the given key.</description>
     !@     <type>\intzero</type>
     !@     <arguments>\textcolor{red}{\textless (character(len=*)|varying\_string|\intzero)\textgreater} key\argin, \intzero\ value\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>key</method>
     !@     <description>Return the key of the {\tt indexValue}$^{\rm th}$ entry in the hash.</description>
     !@     <type>\textcolor{red}{\textless type(varying\_string)\textgreater}</type>
     !@     <arguments>\intzero\ indexValue\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>keys</method>
     !@     <description>Return an array of all keys in the hash.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(varying\_string)[:]\textgreater} keys\arginout</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>values</method>
     !@     <description>Return an array of all values in the hash.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless integer[:]\textgreater} values\arginout</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>exists</method>
     !@     <description>Return true if the specified key exists in the hash.</description>
     !@     <type>\logicalzero</type>
     !@     <arguments>\textcolor{red}{\textless (character(len=*)|varying\_string|\intzero)\textgreater} key\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>size</method>
     !@     <description>Return the number of keys in the hash.</description>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: initialize           =>Initialize_Integer_Scalar 
     procedure :: Set_Integer_Scalar_VS                            
     procedure :: Set_Integer_Scalar_CH                            
     generic   :: set        => Set_Integer_Scalar_VS,Set_Integer_Scalar_CH
     procedure :: Delete_Integer_Scalar_VS 
     procedure :: Delete_Integer_Scalar_CH 
     generic   :: delete     => Delete_Integer_Scalar_VS,Delete_Integer_Scalar_CH
     procedure :: Value_Integer_Scalar_VS 
     procedure :: Value_Integer_Scalar_CH 
     procedure :: Value_Integer_Scalar_I  
     generic   :: value      => Value_Integer_Scalar_VS,Value_Integer_Scalar_CH,Value_Integer_Scalar_I
     procedure :: key                     =>Key_Integer_Scalar_I  
     procedure :: Exists_Integer_Scalar_VS                        
     procedure :: Exists_Integer_Scalar_CH                        
     procedure :: keys                    =>Keys_Integer_Scalar   
     procedure :: values                  =>Values_Integer_Scalar 
     generic   :: exists     => Exists_Integer_Scalar_VS,Exists_Integer_Scalar_CH
     procedure :: size=>Size_Integer_Scalar 
  end type integerScalarHash

  ! The number of new elements by which to extend hashes that need to grow.
  integer, parameter :: hashSizeIncrement=128 
  
contains

  subroutine Initialize_Integer_Scalar(thisHash)
    !% Routine to initialize (or re-initialize) an integer hash.
    implicit none
    class(integerScalarHash), intent(  out) :: thisHash 
    
    select type (thisHash)
    type is (integerScalarHash)
       thisHash%elementCount =0
       thisHash%allocatedSize=0
       if (allocated(thisHash%hashValues)) deallocate(thisHash%hashValues)
       if (allocated(thisHash%hashKeys  )) deallocate(thisHash%hashKeys  )
    end select
    return
  end subroutine Initialize_Integer_Scalar

  integer function Size_Integer_Scalar(thisHash)
    !% Returns the number of elements in the specified {\tt Hash}.
    implicit none
    class(integerScalarHash), intent(in   ) :: thisHash 
    
    select type (thisHash)
    type is (integerScalarHash)
       Size_Integer_Scalar=thisHash%elementCount
    end select
    return
  end function Size_Integer_Scalar

  logical function Exists_Integer_Scalar_CH(thisHash,keyCH)
    !% Returns true if the specified {\tt key} exists in the specified {\tt thisHash}, false otherwise.
    implicit none
    class    (integerScalarHash), intent(in   ) :: thisHash 
    character(len=*            ), intent(in   ) :: keyCH    
    type     (varying_string   ), save          :: key      
    !$omp threadprivate(key)
    key=trim(keyCH)
    Exists_Integer_Scalar_CH=Exists_Integer_Scalar_VS(thisHash,key)
    return
  end function Exists_Integer_Scalar_CH

  logical function Exists_Integer_Scalar_VS(thisHash,key)
    !% Returns true if the specified {\tt key} exists in the specified {\tt thisHash}, false otherwise.
    implicit none
    class(integerScalarHash), intent(in   ) :: thisHash 
    type (varying_string   ), intent(in   ) :: key      
    
    select type (thisHash)
    type is (integerScalarHash)
       if (thisHash%elementCount > 0) then
          Exists_Integer_Scalar_VS=any(thisHash%hashKeys(1:thisHash%elementCount) == key)
       else
          Exists_Integer_Scalar_VS=.false.
       end if
    end select
    return
  end function Exists_Integer_Scalar_VS
  
  subroutine Delete_Integer_Scalar_CH(thisHash,keyCH)
    !% Deletes entry {\tt key} from {\tt thisHash}.
    implicit none
    character(len=*            ), intent(in   ) :: keyCH    
    class    (integerScalarHash), intent(inout) :: thisHash 
    type     (varying_string   ), save          :: key      
    !$omp threadprivate(key)
    key=trim(keyCH)
    call Delete_Integer_Scalar_VS(thisHash,key)
    return
  end subroutine Delete_Integer_Scalar_CH

  subroutine Delete_Integer_Scalar_VS(thisHash,key)
    !% Deletes entry {\tt key} from {\tt Hash}.
    use Arrays_Search
    use Galacticus_Error
    implicit none
    type   (varying_string   ), intent(in   ) :: key      
    class  (integerScalarHash), intent(inout) :: thisHash 
    integer                   , save          :: iKey     
    
    select type (thisHash)
    type is (integerScalarHash)
       if (Exists_Integer_Scalar_VS(thisHash,key)) then
          iKey=Search_Array(thisHash%hashKeys(1:thisHash%elementCount),key)
          thisHash%hashKeys  (ikey:thisHash%elementCount-1)=thisHash%hashKeys  (ikey+1:thisHash%elementCount)
          thisHash%hashValues(ikey:thisHash%elementCount-1)=thisHash%hashValues(ikey+1:thisHash%elementCount)
          thisHash%elementCount                        =thisHash%elementCount-1
       else
          call Galacticus_Error_Report('Delete_Integer_Scalar_VS','key '''//char(key)//''' does not exist in hash')
       end if
    end select
    return
  end subroutine Delete_Integer_Scalar_VS

  function Key_Integer_Scalar_I(thisHash,indexValue) result (key)
    !% Returns the key of entry number {\tt index} in {\tt thisHash}.
    use ISO_Varying_String
    implicit none
    type   (varying_string   )                :: key        
    integer                   , intent(in   ) :: indexValue 
    class  (integerScalarHash), intent(in   ) :: thisHash   
    
    select type (thisHash)
    type is (integerScalarHash)
       key=thisHash%hashKeys(indexValue)
    end select
    return
  end function Key_Integer_Scalar_I

  subroutine Keys_Integer_Scalar(thisHash,keys)
    !% Returns an array of all keys in {\tt thisHash}.
    implicit none
    type (varying_string   ), allocatable, dimension(:), intent(inout) :: keys     
    class(integerScalarHash)                           , intent(in   ) :: thisHash 
    
    select type (thisHash)
    type is (integerScalarHash)
       if (allocated(keys)) deallocate(keys)
       allocate(keys(thisHash%elementCount))
       keys=thisHash%hashKeys(1:thisHash%elementCount)
    end select
    return
  end subroutine Keys_Integer_Scalar

  subroutine Values_Integer_Scalar(thisHash,values)
    !% Returns an array of all values in {\tt thisHash}.
    implicit none
    integer                   , allocatable, dimension(:), intent(inout) :: values   
    class  (integerScalarHash)                           , intent(in   ) :: thisHash 
    
    select type (thisHash)
    type is (integerScalarHash)
       if (allocated(values)) deallocate(values)
       allocate(values(thisHash%elementCount))
       values=thisHash%hashValues(1:thisHash%elementCount)
    end select
    return
  end subroutine Values_Integer_Scalar

  integer function Value_Integer_Scalar_I(thisHash,indexValue)
    !% Returns the value of entry number {\tt index} in {\tt Hash}.
    implicit none
    class  (integerScalarHash), intent(in   ) :: thisHash   
    integer                   , intent(in   ) :: indexValue 
    
    select type (thisHash)
    type is (integerScalarHash)
       Value_Integer_Scalar_I=thisHash%hashValues(indexValue)
    end select
    return
  end function Value_Integer_Scalar_I

  integer function Value_Integer_Scalar_CH(thisHash,keyCH)
    !% Returns the value of {\tt Key} in {\tt Hash}.
    implicit none
    character(len=*            ), intent(in   ) :: keyCH    
    class    (integerScalarHash), intent(in   ) :: thisHash 
    type     (varying_string   ), save          :: key      
    !$omp threadprivate(key)
    key=trim(keyCH)
    Value_Integer_Scalar_CH=Value_Integer_Scalar_VS(thisHash,key)
    return
  end function Value_Integer_Scalar_CH

  integer function Value_Integer_Scalar_VS(thisHash,key)
    !% Returns the value of {\tt key} in {\tt thisHash}.
    use Arrays_Search
    use Galacticus_Error
    implicit none
    class  (integerScalarHash), intent(in   ) :: thisHash 
    type   (varying_string   ), intent(in   ) :: key      
    integer                                   :: iKey     
    
    select type (thisHash)
    type is (integerScalarHash)
       if (Exists_Integer_Scalar_VS(thisHash,key)) then
          iKey=Search_Array(thisHash%hashKeys(1:thisHash%elementCount),key)
          Value_Integer_Scalar_VS=thisHash%hashValues(iKey)
       else
          call Galacticus_Error_Report('Value_Integer_Scalar','key '''//char(key)//''' does not exist in hash')
       end if
    end select
    return
  end function Value_Integer_Scalar_VS

  subroutine Set_Integer_Scalar_CH(thisHash,keyCH,value)
    !% Sets the value of {\tt key} in {\tt thisHash} to {\tt value}.
    use Arrays_Search
    implicit none
    integer                     , intent(in   ) :: value    
    character(len=*            ), intent(in   ) :: keyCH    
    class    (integerScalarHash), intent(inout) :: thisHash 
    type     (varying_string   ), save          :: key      
    
    key=trim(keyCH)
    call Set_Integer_Scalar_VS(thisHash,key,value)
    return
  end subroutine Set_Integer_Scalar_CH

  subroutine Set_Integer_Scalar_VS(thisHash,key,value)
    !% Sets the value of {\tt key} in {\tt thisHash} to {\tt value}.
    use Arrays_Search
    implicit none
    integer                                              , intent(in   ) :: Value           
    type   (varying_string   )                           , intent(in   ) :: Key             
    class  (integerScalarHash)                           , intent(inout) :: thisHash        
    integer                                                              :: iKey            
    logical                                                              :: keyExists       
    integer                   , allocatable, dimension(:)                :: valuesTemporary 
    type   (varying_string   ), allocatable, dimension(:)                :: keysTemporary   
    
    select type (thisHash)
    type is (integerScalarHash)
       ! Check if key already exists.
       if (thisHash%elementCount > 0) then
          keyExists=any(thisHash%hashKeys(1:thisHash%elementCount) == key)
       else
          keyExists=.false.
       end if
       if (keyExists) then
          iKey=Search_Array(thisHash%hashKeys(1:thisHash%elementCount),key)
          thisHash%hashValues(iKey)=value
       else
          ! Increase hash size if necessary.
          if (thisHash%elementCount == thisHash%allocatedSize) then
             if (thisHash%allocatedSize > 0) then
                allocate(valuesTemporary(thisHash%allocatedSize))
                allocate(keysTemporary  (thisHash%allocatedSize))
                valuesTemporary=thisHash%hashValues
                keysTemporary  =thisHash%hashKeys
                deallocate(thisHash%hashValues)
                deallocate(thisHash%hashKeys  )
                thisHash%allocatedSize=thisHash%allocatedSize+hashSizeIncrement
                allocate(thisHash%hashValues(thisHash%allocatedSize))
                allocate(thisHash%hashKeys  (thisHash%allocatedSize))
                thisHash%hashValues(1:size(valuesTemporary))=valuesTemporary
                thisHash%hashKeys  (1:size(valuesTemporary))=keysTemporary
                deallocate(valuesTemporary)
                deallocate(keysTemporary  )
             else
                thisHash%allocatedSize=hashSizeIncrement
                allocate(thisHash%hashValues(thisHash%allocatedSize))
                allocate(thisHash%hashKeys  (thisHash%allocatedSize))
             end if
          end if
          if (thisHash%elementCount > 0) then
             iKey=Search_Array(thisHash%hashKeys(1:thisHash%elementCount),key)
          else
             iKey=1
          end if
          if (iKey > thisHash%elementCount) then
             ! Insert at end.
             thisHash%elementCount                 =thisHash%elementCount+1
             thisHash%hashKeys  (thisHash%elementCount)=key
             thisHash%hashValues(thisHash%elementCount)=value
          else
             ! Shift array then insert.
             thisHash%hashKeys        (iKey+2:thisHash%elementCount+1)=thisHash%hashKeys  (iKey+1:thisHash%elementCount)
             thisHash%hashValues      (iKey+2:thisHash%elementCount+1)=thisHash%hashValues(iKey+1:thisHash%elementCount)
             thisHash%hashKeys        (iKey+1                        )=key
             thisHash%hashValues      (iKey+1                        )=value
             thisHash%elementCount                                    =thisHash%elementCount+1
          end if
       end if
    end select
    return
  end subroutine Set_Integer_Scalar_VS

end module Hashes
