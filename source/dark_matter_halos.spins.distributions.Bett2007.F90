!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% An implementation of the dark matter halo spin distribution which uses the fitting function proposed by
  !% \cite{bett_spin_2007}.

  use Tables
  use Table_Labels

  !# <haloSpinDistribution name="haloSpinDistributionBett2007">
  !#  <description>A halo spin distribution using the fitting formula of \cite{bett_spin_2007}.</description>
  !# </haloSpinDistribution>
  type, extends(haloSpinDistributionClass) :: haloSpinDistributionBett2007
     !% A dark matter halo spin distribution class which assumes a \cite{bett_spin_2007} distribution.
     private
     double precision                                        :: alpha               , lambda0
     type            (table1DLogarithmicLinear)              :: distribution
     class           (table1D                 ), allocatable :: distributionInverse
     type            (fgsl_rng                )              :: clonedPseudoSequence, randomSequence
     logical                                                 :: resetRandomSequence , resetRandomSequenceSnapshot
   contains
     final     ::                  bett2007Destructor
     procedure :: sample        => bett2007Sample
     procedure :: stateSnapshot => bett2007StateSnapshot
     procedure :: stateStore    => bett2007StateStore
     procedure :: stateRestore  => bett2007StateRestore
  end type haloSpinDistributionBett2007
  
  interface haloSpinDistributionBett2007
     !% Constructors for the {\normalfont \ttfamily bett2007} dark matter halo spin
     !% distribution class.
     module procedure bett2007ConstructorParameters
     module procedure bett2007ConstructorInternal
  end interface haloSpinDistributionBett2007

  ! Tabulation parameters.
  integer         , parameter :: bett2007TabulationPointsCount=1000
  double precision, parameter :: bett2007SpinMaximum          =   0.2d+0
  double precision, parameter :: bett2007SpinMinimum          =   1.0d-6

contains

  function bett2007ConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily bett2007} dark matter halo spin
    !% distribution class which takes a parameter list as input.
    use Input_Parameters2
    implicit none
    type            (haloSpinDistributionBett2007)                :: bett2007ConstructorParameters
    type            (inputParameters             ), intent(in   ) :: parameters
    double precision                                              :: lambda0                      , alpha
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>lambda0</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.04326d0</defaultValue>
    !#   <defaultSource>\citep{bett_spin_2007}</defaultSource>
    !#   <description>The parameter $\lambda_0$ in the halo spin distribution of \cite{bett_spin_2007}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>alpha</name>
    !#   <source>parameters</source>
    !#   <defaultValue>2.509d0</defaultValue>
    !#   <defaultSource>\citep{bett_spin_2007}</defaultSource>
    !#   <description>The parameter $\alpha$ in the halo spin distribution of \cite{bett_spin_2007}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    bett2007ConstructorParameters=bett2007ConstructorInternal(lambda0,alpha)
    return
  end function bett2007ConstructorParameters

  function bett2007ConstructorInternal(lambda0,alpha)
    !% Internal constructor for the {\normalfont \ttfamily bett2007} dark matter halo spin
    !% distribution class.
    use, intrinsic :: ISO_C_Binding
    use Gamma_Functions
    implicit none
    type            (haloSpinDistributionBett2007)                :: bett2007ConstructorInternal
    double precision                              , intent(in   ) :: lambda0                    , alpha
    double precision                                              :: spinDimensionless          , tableMaximum
    integer                                                       :: iSpin

    bett2007ConstructorInternal%lambda0                    =lambda0
    bett2007ConstructorInternal%alpha                      =alpha
    bett2007ConstructorInternal%resetRandomSequence        =.true.
    bett2007ConstructorInternal%resetRandomSequenceSnapshot=.true.
    if (FGSL_Well_Defined(bett2007ConstructorInternal%      randomSequence))                &
         & call FGSL_Obj_C_Ptr(bett2007ConstructorInternal%      randomSequence,C_Null_Ptr)
    if (FGSL_Well_Defined(bett2007ConstructorInternal%clonedPseudoSequence))                &
         & call FGSL_Obj_C_Ptr(bett2007ConstructorInternal%clonedPseudoSequence,C_Null_Ptr)
    ! Tabulate the cumulative distribution.
    tableMaximum=(bett2007SpinMaximum/lambda0)**(3.0d0/alpha)
    call bett2007ConstructorInternal%distribution%destroy()
    call bett2007ConstructorInternal%distribution%create (                                                                                        &
         &                                                lambda0*bett2007SpinMinimum**(alpha/3.0d0),                                             &
         &                                                lambda0*tableMaximum       **(alpha/3.0d0)                                            , &
         &                                                bett2007TabulationPointsCount                                                         , &
         &                                                extrapolationType                         =[extrapolationTypeFix,extrapolationTypeFix]  &
         &                                               )
    ! Compute the cumulative probability distribution.
    do iSpin=1,bett2007TabulationPointsCount
       spinDimensionless=(                                                   &
            &             +bett2007ConstructorInternal%distribution%x(iSpin) &
            &             /lambda0                                           &
            &            )**(3.0d0/alpha)
       call bett2007ConstructorInternal%distribution%populate(                                                            &
            &                                                 Gamma_Function_Incomplete_Complementary(                    &
            &                                                                                         +alpha            , &
            &                                                                                         +alpha              &
            &                                                                                         *spinDimensionless  &
            &                                                                                        )                  , &
            &                                                 iSpin                                                       &
            &                                                )
    end do
    call bett2007ConstructorInternal%distribution%reverse(bett2007ConstructorInternal%distributionInverse)
    return
  end function bett2007ConstructorInternal

  subroutine bett2007Destructor(self)
    !% Destructor for the {\normalfont \ttfamily bett2007} dark matter halo spin
    !% distribution class.
    use Pseudo_Random
    implicit none
    type(haloSpinDistributionBett2007), intent(inout) :: self

    if (.not.self%resetRandomSequence        ) call Pseudo_Random_Free(self%randomSequence      )
    if (.not.self%resetRandomSequenceSnapshot) call Pseudo_Random_Free(self%clonedPseudoSequence)
    return
  end subroutine bett2007Destructor

  double precision function bett2007Sample(self,node)
    !% Sample from a \cite{bett_spin_2007} spin parameter distribution for the given {\normalfont
    !% \ttfamily node}.
    use Pseudo_Random
    implicit none
    class(haloSpinDistributionBett2007), intent(inout)          :: self
    type (treeNode                    ), intent(inout), pointer :: node

    bett2007Sample=self%distributionInverse%interpolate(Pseudo_Random_Get(self%randomSequence,self%resetRandomSequence))
    return
  end function bett2007Sample

  subroutine bett2007StateSnapshot(self)
    !% Store a snapshot of the random number generator internal state.
    use Pseudo_Random
    implicit none
    class(haloSpinDistributionBett2007), intent(inout) :: self
    
    if (.not.self%resetRandomSequence) then
       if (FGSL_Well_Defined(self%clonedPseudoSequence)) call Pseudo_Random_Free(self%clonedPseudoSequence)
       self%clonedPseudoSequence=FGSL_Rng_Clone(self%randomSequence)
    end if
    self%resetRandomSequenceSnapshot=self%resetRandomSequence
    return
  end subroutine bett2007StateSnapshot

  subroutine bett2007StateStore(self,stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use Pseudo_Random
    implicit none
    class  (haloSpinDistributionBett2007), intent(inout) :: self
    integer                               , intent(in   ) :: stateFile
    type   (fgsl_file                    ), intent(in   ) :: fgslStateFile

    write (stateFile) self%resetRandomSequenceSnapshot
    if (.not.self%resetRandomSequenceSnapshot) call Pseudo_Random_Store(self%clonedPseudoSequence,fgslStateFile)
    return
  end subroutine bett2007StateStore

  subroutine bett2007StateRestore(self,stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use Pseudo_Random
    implicit none
    class  (haloSpinDistributionBett2007), intent(inout) :: self
    integer                               , intent(in   ) :: stateFile
    type   (fgsl_file                    ), intent(in   ) :: fgslStateFile

    read (stateFile) self%resetRandomSequence
    if (.not.self%resetRandomSequence) call Pseudo_Random_Retrieve(self%randomSequence,fgslStateFile)
    return
  end subroutine bett2007StateRestore
