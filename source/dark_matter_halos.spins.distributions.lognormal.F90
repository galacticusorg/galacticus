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

  !% An implementation of the dark matter halo spin distribution which assumes a
  !% log-normal distribution.

  !# <haloSpinDistribution name="haloSpinDistributionLogNormal">
  !#  <description>A log-normal halo spin distribution.</description>
  !# </haloSpinDistribution>
  type, extends(haloSpinDistributionClass) :: haloSpinDistributionLogNormal
     !% A dark matter halo spin distribution concentration class which assumes a
     !% log-normal distribution.
     private
     double precision           :: median              , sigma
     type            (fgsl_rng) :: clonedPseudoSequence, randomSequence
     logical                    :: resetRandomSequence , resetRandomSequenceSnapshot
   contains
     final     ::                  logNormalDestructor
     procedure :: sample        => logNormalSample
     procedure :: stateSnapshot => logNormalStateSnapshot
     procedure :: stateStore    => logNormalStateStore
     procedure :: stateRestore  => logNormalStateRestore
  end type haloSpinDistributionLogNormal
  
  interface haloSpinDistributionLogNormal
     !% Constructors for the {\normalfont \ttfamily logNormal} dark matter halo spin
     !% distribution class.
     module procedure logNormalConstructorParameters
     module procedure logNormalConstructorInternal
  end interface haloSpinDistributionLogNormal

contains

  function logNormalConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily logNormal} dark matter halo spin
    !% distribution class which takes a parameter list as input.
    use Input_Parameters2
    implicit none
    type(haloSpinDistributionLogNormal)                :: logNormalConstructorParameters
    type(inputParameters              ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>median</name>
    !#   <source>parameters</source>
    !#   <variable>logNormalConstructorParameters%median</variable>
    !#   <defaultValue>0.03687d0</defaultValue>
    !#   <defaultSource>\citep{bett_spin_2007}</defaultSource>
    !#   <description>The median spin in a log-normal spin distribution.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sigma</name>
    !#   <source>parameters</source>
    !#   <variable>logNormalConstructorParameters%sigma</variable>
    !#   <defaultValue>0.2216d0</defaultValue>
    !#   <defaultSource>\citep{bett_spin_2007}</defaultSource>
    !#   <description>The width of a log-normal spin distribution.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    logNormalConstructorParameters%median=log(logNormalConstructorParameters%median)
    return
  end function logNormalConstructorParameters

  function logNormalConstructorInternal(median,sigma)
    !% Internal constructor for the {\normalfont \ttfamily logNormal} dark matter halo spin
    !% distribution class.
    implicit none
    type            (haloSpinDistributionLogNormal)                :: logNormalConstructorInternal
    double precision                               , intent(in   ) :: median                      , sigma

    logNormalConstructorInternal%median=log(median)
    logNormalConstructorInternal%sigma =    sigma
    return
  end function logNormalConstructorInternal

  subroutine logNormalDestructor(self)
    !% Destructor for the {\normalfont \ttfamily logNormal} dark matter halo spin
    !% distribution class.
    use Pseudo_Random
    implicit none
    type(haloSpinDistributionLogNormal), intent(inout) :: self

    if (.not.self%resetRandomSequence        ) call Pseudo_Random_Free(self%randomSequence      )
    if (.not.self%resetRandomSequenceSnapshot) call Pseudo_Random_Free(self%clonedPseudoSequence)
    return
  end subroutine logNormalDestructor

  double precision function logNormalSample(self,node)
    !% Sample from a log-normal spin parameter distribution for the given {\normalfont
    !% \ttfamily node}.
    use Gaussian_Random
    implicit none
    class(haloSpinDistributionLogNormal), intent(inout)          :: self
    type (treeNode                     ), intent(inout), pointer :: node
    !GCC$ attributes unused :: node

    logNormalSample=exp(                                               &
         &              +self%median                                   &
         &              +Gaussian_Random_Get(                          &
         &                                   self%randomSequence     , &
         &                                   self%sigma              , &
         &                                   self%resetRandomSequence  &
         &                                   )                         &
         &             )
    return
  end function logNormalSample

  subroutine logNormalStateSnapshot(self)
    !% Store a snapshot of the random number generator internal state.
    implicit none
    class(haloSpinDistributionLogNormal), intent(inout) :: self

    if (.not.self%resetRandomSequence) self%clonedPseudoSequence=FGSL_Rng_Clone(self%randomSequence)
    self%resetRandomSequenceSnapshot=self%resetRandomSequence
    return
  end subroutine logNormalStateSnapshot

  subroutine logNormalStateStore(self,stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use Galacticus_Display
    use Pseudo_Random
    implicit none
    class  (haloSpinDistributionLogNormal), intent(inout) :: self
    integer                               , intent(in   ) :: stateFile
    type   (fgsl_file                    ), intent(in   ) :: fgslStateFile

    call Galacticus_Display_Message('Storing state for: haloSpinDistribution -> logNormal',verbosity=verbosityInfo)
    write (stateFile) self%resetRandomSequenceSnapshot
    if (.not.self%resetRandomSequenceSnapshot) call Pseudo_Random_Store(self%clonedPseudoSequence,fgslStateFile)
    return
  end subroutine logNormalStateStore

  subroutine logNormalStateRestore(self,stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use Galacticus_Display
    use Pseudo_Random
    implicit none
    class  (haloSpinDistributionLogNormal), intent(inout) :: self
    integer                               , intent(in   ) :: stateFile
    type   (fgsl_file                    ), intent(in   ) :: fgslStateFile

    call Galacticus_Display_Message('Retrieving state for: haloSpinDistribution -> logNormal',verbosity=verbosityInfo)
    read (stateFile) self%resetRandomSequence
    if (.not.self%resetRandomSequence) call Pseudo_Random_Retrieve(self%randomSequence,fgslStateFile)
    return
  end subroutine logNormalStateRestore
