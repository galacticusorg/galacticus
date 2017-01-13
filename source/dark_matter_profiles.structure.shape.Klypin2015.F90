!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

  !% An implementation of dark matter halo profile shapes  using the \cite{klypin_multidark_2014} algorithm.

  !# <darkMatterProfileShape name="darkMatterProfileShapeKlypin2015">
  !#  <description>Dark matter halo shape parameters are computed using the algorithm of \cite{klypin_multidark_2014}.</description>
  !# </darkMatterProfileShape>

  type, extends(darkMatterProfileShapeClass) :: darkMatterProfileShapeKlypin2015
     !% A dark matter halo profile shape parameter class implementing the algorithm of \cite{klypin_multidark_2014}.
     private
     integer   :: sample
   contains
     procedure :: shape => klypin2015Shape
  end type darkMatterProfileShapeKlypin2015

  interface darkMatterProfileShapeKlypin2015
     !% Constructors for the {\normalfont \ttfamily klypin2015} dark matter halo profile shape parameter class.
     module procedure klypin2015DefaultConstructor
     module procedure klypin2015Constructor
  end interface darkMatterProfileShapeKlypin2015

  ! Labels for sample selection.
  integer, parameter :: klypin2015SampleAll          =0
  integer, parameter :: klypin2015SampleRelaxed      =1
  
  ! Default sample.
  integer            :: klypin2015ConcentrationSample

  ! Initialization status.
  logical            :: klypin2015Initialized        =.false.
  
contains

  function klypin2015DefaultConstructor()
    !% Default constructor for the {\normalfont \ttfamily klypin2015} dark matter halo profile shape parameter class.
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type(darkMatterProfileShapeKlypin2015), target :: klypin2015DefaultConstructor
    type(varying_string                  )         :: klypin2015ShapeSampleText
    
    if (.not.klypin2015Initialized) then
       !$omp critical(klypin2015DefaultInitialize)
       if (.not.klypin2015Initialized) then
          ! Get parameters of the model.
          !@ <inputParameter>
          !@   <name>klypin2015ShapeSample</name>
          !@   <defaultValue>all</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The sample to use for the halo shape parameter algorithm of \cite{klypin_multidark_2014}.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("klypin2015ShapeSample",klypin2015ShapeSampleText,defaultValue='all')
          select case (char(klypin2015ShapeSampleText))
          case ('all'     )
             klypin2015ConcentrationSample=klypin2015SampleAll
          case ('relaxed' )
             klypin2015ConcentrationSample=klypin2015SampleRelaxed
          case default
             call Galacticus_Error_Report('klypin2015DefaultConstructor','unrecognized sample')
          end select
          ! Record that method is now initialized.
          klypin2015Initialized=.true.
       end if
       !$omp end critical(klypin2015DefaultInitialize)
    end if
    ! Construct the object.
    klypin2015DefaultConstructor=klypin2015Constructor(klypin2015ConcentrationSample)
    return
  end function klypin2015DefaultConstructor
  
  function klypin2015Constructor(sample)
    !% Constructor for the {\normalfont \ttfamily klypin2015} dark matter halo profile shape parameter class.
    use Galacticus_Error
    implicit none
    type   (darkMatterProfileShapeKlypin2015)                :: klypin2015Constructor
    integer                                  , intent(in   ) :: sample
    
    select case (sample)
    case (klypin2015SampleAll    )
       klypin2015Constructor%sample=klypin2015SampleAll
    case (klypin2015SampleRelaxed)
       klypin2015Constructor%sample=klypin2015SampleRelaxed
    case default
       call Galacticus_Error_Report('klypin2015DefaultConstructor','unrecognized sample')
    end select
    return
  end function klypin2015Constructor
  
  double precision function klypin2015Shape(self,node)
    !% Return the Einasto profile shape parameter of the dark matter halo profile of {\normalfont \ttfamily node} using the
    !% \cite{klypin_multidark_2014} algorithm.
    use Galacticus_Nodes
    use Galacticus_Error
    use Cosmological_Mass_Variance
    use Critical_Overdensities
    implicit none
    class           (darkMatterProfileShapeKlypin2015), intent(inout)          :: self
    type            (treeNode                        ), intent(inout), pointer :: node
    class           (criticalOverdensityClass        )               , pointer :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass   )               , pointer :: cosmologicalMassVariance_
    class           (nodeComponentBasic              )               , pointer :: basic
    double precision                                                           :: nu
    
    ! Get default objects.
    criticalOverdensity_      => criticalOverdensity     ()
    cosmologicalMassVariance_ => cosmologicalMassVariance()
    ! Get the basic component.
    basic => node%basic()
    ! Compute the shape parameter.
    nu     =+criticalOverdensity_     %value       (time=basic%time(),mass=basic%mass()) &
         &  /cosmologicalMassVariance_%rootVariance(                       basic%mass())
    select case (self%sample)
    case (klypin2015SampleAll    )
       klypin2015Shape=0.115d0+0.0165d0*nu**2
    case (klypin2015SampleRelaxed)
       klypin2015Shape=0.115d0+0.0140d0*nu**2
    case default
       klypin2015Shape=0.0d0
       call Galacticus_Error_Report('klypin2015Shape','unknown sample')
    end select
    return
  end function klypin2015Shape
