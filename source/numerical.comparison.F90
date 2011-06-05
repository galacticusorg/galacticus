!% Contains a module which implements comparisons of values.

module Numerical_Comparison
  !% Implements comparisons of values.
  private
  public :: Values_Differ

contains

  logical function Values_Differ(value1,value2,absTol,relTol)
    !% Returns true if {\tt value1} and {\tt value2} differ by more than {\tt absTol} in absolute terms, or {\tt relTol} in
    !% relative terms.
    implicit none
    double precision, intent(in)           :: value1,value2
    double precision, intent(in), optional :: absTol,relTol
    
    Values_Differ=.false.
    if (present(absTol)) Values_Differ=(abs(value1-value2) > absTol)
    if (present(relTol)) Values_Differ=Values_Differ.or.(abs(value1-value2) > 0.5d0*abs(value1+value2)*relTol)
    if (.not.(present(absTol).or.present(relTol))) Values_Differ=(value1 /= value2)
    return
  end function Values_Differ
  
end module Numerical_Comparison
