!% Contains a module which implements useful operations on arrays.

module Array_Utilities
  !% Contains routines which implement useful operations on arrays.
  private
  public :: Array_Reverse

  interface Array_Reverse
     !% Interface to generic routines which reverse the direction of an array.
     module procedure Array_Reverse_Real
     module procedure Array_Reverse_Double
  end interface

contains

  function Array_Reverse_Real(array) result (reversedArray)
    !% Reverses the direction of a real array.
    implicit none
    real, intent(in)             :: array(:)
    real, dimension(size(array)) :: reversedArray
    integer                      :: i
    
    forall (i=1:size(array))
       reversedArray(i)=array(size(array)+1-i)
    end forall
    return
  end function Array_Reverse_Real

  function Array_Reverse_Double(array) result (reversedArray)
    !% Reverses the direction of a double precision array.
    implicit none
    double precision, intent(in)             :: array(:)
    double precision, dimension(size(array)) :: reversedArray
    integer                                  :: i
    
    forall (i=1:size(array))
       reversedArray(i)=array(size(array)+1-i)
    end forall
    return
  end function Array_Reverse_Double

end module Array_Utilities
