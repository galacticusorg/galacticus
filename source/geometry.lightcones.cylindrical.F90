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

  !% An implementation of the lightcone geometry class which assumes a cylindrical ``cone''.

  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Kind_Numbers            , only : kind_int8
  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  use :: Output_Times            , only : outputTimesClass

  !# <geometryLightcone name="geometryLightconeCylindrical">
  !#  <description>
  !#   A lightcone geometry class which assumes a cylindrical ``cone'', i.e. defined such that a point $(x,y,z)$ is in the survey
  !#   if $\sqrt{x^2+y^2} &lt; r$, where $r$ is the radius of the ``cone''.
  !#  </description>
  !# </geometryLightcone>
  type, extends(geometryLightconeClass) :: geometryLightconeCylindrical
     !% A lightcone geometry class which assumes a cylindrical ``cone''.
     private
     class           (cosmologyFunctionsClass   ), pointer                     :: cosmologyFunctions_    => null()
     class           (outputTimesClass          ), pointer                     :: outputTimes_           => null()
     class           (randomNumberGeneratorClass), pointer                     :: randomNumberGenerator_ => null()
     double precision                            , dimension(:  ), allocatable :: distanceMinimum                 , distanceMaximum     , &
          &                                                                       volume
     double precision                                                          :: radiusCylinderComoving          , radiusBufferComoving
     integer         (kind_int8                 )                              :: activeCenterUniqueID            , activeUniqueID
     integer                                                                   :: activeCenterCount               , activeCount
     logical                                     , dimension(:  ), allocatable :: activeInCylinder
     integer                                     , dimension(:  ), allocatable :: activeInstances
     double precision                            , dimension(:,:), allocatable :: activeCenterPosition            , activePosition
     integer         (c_size_t                  )                              :: activeOutput
   contains
     !@ <objectMethods>
     !@   <object>geometryLightconeCylindrical</object>
     !@   <objectMethod>
     !@     <method>sampleNode</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(treeNode)\textgreater} node\arginout</arguments>
     !@     <description>Determine if, and how many times, the given node appears in the ``lightcone'', and choose positions for each instance.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                     cylindricalDestructor
     procedure :: isInLightcone    => cylindricalIsInLightcone
     procedure :: replicationCount => cylindricalReplicationCount
     procedure :: solidAngle       => cylindricalSolidAngle
     procedure :: position         => cylindricalPosition
     procedure :: velocity         => cylindricalVelocity
     procedure :: sampleNode       => cylindricalSampleNode
  end type geometryLightconeCylindrical

  interface geometryLightconeCylindrical
     !% Constructors for the {\normalfont \ttfamily cylindrical} dark matter halo spin distribution class.
     module procedure cylindricalConstructorParameters
     module procedure cylindricalConstructorInternal
  end interface geometryLightconeCylindrical

  ! Tolerance for matching to output times.
  double precision, parameter :: timeTolerance =1.0d-3

contains

  function cylindricalConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily cylindrical} lightcone geometry distribution class which takes a parameter list as
    !% input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (geometryLightconeCylindrical)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (outputTimesClass            ), pointer       :: outputTimes_
    class           (cosmologyFunctionsClass     ), pointer       :: cosmologyFunctions_
    class           (randomNumberGeneratorClass  ), pointer       :: randomNumberGenerator_
    double precision                                              :: radiusCylinderComoving, radiusBufferComoving

    !# <inputParameter>
    !#   <name>radiusCylinderComoving</name>
    !#   <source>parameters</source>
    !#   <description>The comoving radius of the cylinder to populate.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>radiusBufferComoving</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The comoving buffer radius to add around the cylinder. This is used to ensure that the sample within the cylinder is complete.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    !# <objectBuilder class="outputTimes"           name="outputTimes_"           source="parameters"/>
    !# <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    self=geometryLightconeCylindrical(radiusCylinderComoving,radiusBufferComoving,cosmologyFunctions_,outputTimes_,randomNumberGenerator_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="outputTimes_"          />
    !# <objectDestructor name="cosmologyFunctions_"   />
    !# <objectDestructor name="randomNumberGenerator_"/>
    return
  end function cylindricalConstructorParameters

  function cylindricalConstructorInternal(radiusCylinderComoving,radiusBufferComoving,cosmologyFunctions_,outputTimes_,randomNumberGenerator_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily cylindrical} lightcone geometry distribution class.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (geometryLightconeCylindrical)                        :: self
    class           (cosmologyFunctionsClass     ), target, intent(in   ) :: cosmologyFunctions_
    class           (outputTimesClass            ), target, intent(in   ) :: outputTimes_
    class           (randomNumberGeneratorClass  ), target, intent(in   ) :: randomNumberGenerator_
    double precision                                      , intent(in   ) :: radiusCylinderComoving, radiusBufferComoving
    integer         (c_size_t                    )                        :: output
    double precision                                                      :: timeMinimum           , timeMaximum
    !# <constructorAssign variables="radiusCylinderComoving, radiusBufferComoving, *cosmologyFunctions_, *outputTimes_, *randomNumberGenerator_"/>
    
    ! Find the minimum and maximum distance associated with each output time.
    allocate(self%distanceMinimum(self%outputTimes_%count()))
    allocate(self%distanceMaximum(self%outputTimes_%count()))
    allocate(self%volume         (self%outputTimes_%count()))
    do output=1,self%outputTimes_%count()
       if (output == 1                        ) then
          timeMinimum=                                      self%outputTimes_%time(output)
       else
          timeMinimum=sqrt(self%outputTimes_%time(output-1)*self%outputTimes_%time(output))
       end if
       if (output == self%outputTimes_%count()) then
          timeMaximum=                                      self%outputTimes_%time(output)
       else
          timeMaximum=sqrt(self%outputTimes_%time(output+1)*self%outputTimes_%time(output))
       end if
       self%distanceMinimum(output)=self%cosmologyFunctions_%distanceComoving(timeMaximum)
       self%distanceMaximum(output)=self%cosmologyFunctions_%distanceComoving(timeMinimum)
    end do
    ! Compute the comoving volume associated with each output time, including the buffer region around our cylinder.
    self%volume=+Pi                            &
         &      *(                             &
         &        +self%radiusCylinderComoving &
         &        +self%radiusBufferComoving   &
         &       )**2                          &
         &      *(                             &
         &        +self%distanceMaximum        &
         &        -self%distanceMinimum        &
         &      )
    ! Initialize the unique ID of the active node to an impossible value.
    self%activeCenterUniqueID=-1_kind_int8
    self%activeUniqueID      =-1_kind_int8
    self%activeOutput        =-1_c_size_t
    return
  end function cylindricalConstructorInternal

  subroutine cylindricalDestructor(self)
    !% Destructor for the {\normalfont \ttfamily cylindrical} lightcone geometry distribution class.
    implicit none
    type(geometryLightconeCylindrical), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"   />
    !# <objectDestructor name="self%outputTimes_"          />
    !# <objectDestructor name="self%randomNumberGenerator_"/>
    return
  end subroutine cylindricalDestructor

  function cylindricalReplicationCount(self,node)
    !% Determine the number of times {\normalfont \ttfamily node} appears in the lightcone.
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    integer(c_size_t                    )                :: cylindricalReplicationCount
    class  (geometryLightconeCylindrical), intent(inout) :: self
    type   (treeNode                    ), intent(inout) :: node

    call self%sampleNode(node)
    cylindricalReplicationCount=self%activeCount
    return
  end function cylindricalReplicationCount
  
  logical function cylindricalIsInLightcone(self,node,atPresentEpoch,radiusBuffer)
    !% Determine if the given {\normalfont \ttfamily node} lies within the lightcone.
    use :: Galacticus_Nodes    , only : nodeComponentBasic 
    use :: Numerical_Comparison, only : Values_Agree
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: String_Handling     , only : operator(//)
    implicit none
    class           (geometryLightconeCylindrical), intent(inout)            :: self
    type            (treeNode                    ), intent(inout)            :: node
    logical                                       , intent(in   ) , optional :: atPresentEpoch
    double precision                              , intent(in   ) , optional :: radiusBuffer
    class           (nodeComponentBasic          )                , pointer  :: basic 
    character       (len=10                      )                           :: label
    type            (varying_string              )                           :: message
    integer         (c_size_t                    )                           :: output
    !# <optionalArgument name="atPresentEpoch" defaultsTo=".true." />

    ! Get the basic component.
    basic => node%basic()
    ! Assume not in the lightcone by default.
    cylindricalIsInLightcone=.false.
    ! Check if this node exists prior to any lightcone time. If it does it will not be output.
    if (basic%time() < self%outputTimes_%time(1_c_size_t)*(1.0d0-timeTolerance)) return
    ! Determine to which output this galaxy corresponds.
    output=self%outputTimes_%index(basic%time(),findClosest=.true.)
    if (atPresentEpoch_) then
       ! We want to check only the current time for this node. Check that the node exists precisely at a lightcone snapshot time,
       ! and then set the maximum output to check to equal to minimum, such that we test only the current time.
       if (.not.Values_Agree(self%outputTimes_%time(output),basic%time(),relTol=timeTolerance)) then
          message=         'failed to find matching time in lightcone'                      //char(10)
          write (label,'(f6.3)') basic%time()
          message=message//'               node time = '//trim(label)//' Gyr'               //char(10)
          write (label,'(f6.3)') self%outputTimes_%time(output)
          message=message//'  closest lightcone time = '//trim(label)//' Gyr ['//output//']'
          call Galacticus_Error_Report(message//{introspection:location})
       end if
       call self%sampleNode(node)
       cylindricalIsInLightcone=self%activeCount > 0_c_size_t
    else
       cylindricalIsInLightcone=.false.
       call Galacticus_Error_Report('not well-defined'//{introspection:location})  
    end if
    return
  end function cylindricalIsInLightcone

  double precision function cylindricalSolidAngle(self)
    !% Return the solid angle (in steradians) of a cylindrical lightcone.
    implicit none
    class(geometryLightconeCylindrical), intent(inout) :: self
    !$GLC attributes unused :: self

    ! Solid angle is not well-defined for this "lightcone" class. Simply return zero.
    cylindricalSolidAngle=0.0d0
    return
  end function cylindricalSolidAngle

  function cylindricalPosition(self,node,instance)
    !% Return the position of the node in lightcone coordinates.
    use            :: Galacticus_Error    , only : Galacticus_Error_Report
    use            :: Galacticus_Nodes    , only : nodeComponentBasic
    use, intrinsic :: ISO_C_Binding       , only : c_size_t
    use            :: ISO_Varying_String  , only : varying_string
    use            :: Numerical_Comparison, only : Values_Agree
    use            :: String_Handling     , only : operator(//)
    implicit none
    double precision                              , dimension(3)           :: cylindricalPosition
    class           (geometryLightconeCylindrical), intent(inout)          :: self
    type            (treeNode                    ), intent(inout), target  :: node
    integer         (c_size_t                    ), intent(in   )          :: instance
    class           (nodeComponentBasic          )               , pointer :: basic
    integer         (c_size_t                    )                         :: output
    character       (len=10                      )                         :: label
    type            (varying_string              )                         :: message

    ! Get the basic component.
    basic => node%basic()
    ! Determine to which output this node corresponds.
    output=self%outputTimes_%index(basic%time(),findClosest=.true.)
    if (.not.Values_Agree(self%outputTimes_%time(output),basic%time(),relTol=timeTolerance)) then
       message=         'failed to find matching time in lightcone'                      //char(10)
       write (label,'(f6.3)') basic%time()
       message=message//'               node time = '//trim(label)//' Gyr'               //char(10)
       write (label,'(f6.3)') self%outputTimes_%time(output)
       message=message//'  closest lightcone time = '//trim(label)//' Gyr ['//output//']'
       call Galacticus_Error_Report(message//{introspection:location})
    end if
    call self%sampleNode(node)
    cylindricalPosition=self%activePosition(:,self%activeInstances(instance))
    return
  end function cylindricalPosition

  function cylindricalVelocity(self,node,instance)
    !% Return the velocity of the node in lightcone coordinates.
    use :: Galacticus_Nodes, only : nodeComponentPosition, treeNode
    implicit none
    double precision                              , dimension(3)  :: cylindricalvelocity
    class           (geometryLightconeCylindrical), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    integer         (c_size_t                    ), intent(in   ) :: instance
    !$GLC attributes unused :: self, node, instance

    ! Currently this class provides no model for velocities.
    cylindricalVelocity=0.0d0
    return
  end function cylindricalVelocity

  subroutine cylindricalSampleNode(self,node)
    !% Determine how many times the given node appears in the ``lightcone''.
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Comparison    , only : Values_Agree
    use :: Galacticus_Error        , only : Galacticus_Error_Report
    use :: Galacticus_Nodes        , only : nodeComponentBasic     , nodeComponentSatellite
    use :: String_Handling         , only : operator(//)
    implicit none
    class           (geometryLightconeCylindrical), intent(inout)          :: self
    type            (treeNode                    ), intent(inout), target  :: node
    class           (nodeComponentBasic          )               , pointer :: basicHost
    class           (nodeComponentSatellite      )               , pointer :: satellite
    type            (treeNode                    )               , pointer :: nodeHost
    double precision                              , dimension(3)           :: positionOffset
    double precision                                                       :: distanceComoving, radiusComoving, &
         &                                                                    theta           , numberMean
    integer                                                                :: i               , j
    character       (len=10                      )                         :: label
    type            (varying_string              )                         :: message
    integer         (c_size_t                    )                         :: output
    
    ! Find the isolated node.
    nodeHost => node
    do while (nodeHost%isSatellite())
       nodeHost => nodeHost%parent
    end do
    ! Get the basic component.
    basicHost => nodeHost%basic()
    ! Determine to which output this isolated node corresponds.
    output=self%outputTimes_%index(basicHost%time(),findClosest=.true.)
    if (.not.Values_Agree(self%outputTimes_%time(output),basicHost%time(),relTol=timeTolerance)) then
       message=         'failed to find matching time in lightcone'                      //char(10)
       write (label,'(f6.3)') basicHost%time()
       message=message//'               node time = '//trim(label)//' Gyr'               //char(10)
       write (label,'(f6.3)') self%outputTimes_%time(output)
       message=message//'  closest lightcone time = '//trim(label)//' Gyr ['//output//']'
       call Galacticus_Error_Report(message//{introspection:location})
    end if
    ! Determine if the output has changed.
    if (output /= self%activeOutput) then
       self%activeOutput        =output
       self%activeCenterUniqueID=-1_c_size_t
       self%activeUniqueID      =-1_c_size_t
    end if
    ! Check if we already have results for this isolated node stored.
    if (nodeHost%uniqueID() /= self%activeCenterUniqueID) then
       self%activeCenterUniqueID=nodeHost%uniqueID()
       self%activeUniqueID      =-1_c_size_t
       ! Compute the number of instances of this halo which appear in the cylinder.
       numberMean            =+self                           %volume       (output    ) &
            &                 *nodeHost%hostTree              %volumeWeight
       self%activeCenterCount=+self    %randomNumberGenerator_%poissonSample(numberMean)
       ! Generate random positions for each instance.
       if (allocated(self%activeCenterPosition)) deallocate(self%activeCenterPosition)
       if (allocated(self%activePosition      )) deallocate(self%activePosition      )
       if (allocated(self%activeInCylinder    )) deallocate(self%activeInCylinder    )
       if (allocated(self%activeInstances     )) deallocate(self%activeInstances     )
       allocate(self%activeCenterPosition(3,self%activeCenterCount))
       allocate(self%activePosition      (3,self%activeCenterCount))
       allocate(self%activeInCylinder    (  self%activeCenterCount))
       allocate(self%activeInstances     (  self%activeCenterCount))
       do i=1,self%activeCenterCount          
          radiusComoving  =sqrt(                                                                                  &
               &                self%randomNumberGenerator_%uniformSample()*(                                     &
               &                                                             +self%radiusCylinderComoving         &
               &                                                             +self%radiusBufferComoving           &
               &                                                            )**2                                  &
               &               )
          theta           =     self%randomNumberGenerator_%uniformSample()*  2.0d0                               &
               &                                                           *  Pi
          distanceComoving=     self%randomNumberGenerator_%uniformSample()*(                                     &
               &                                                             +self%distanceMaximum       (output) &
               &                                                             -self%distanceMinimum       (output) &
               &                                                            )                                     &
               &                                                             +self%distanceMinimum       (output)
          self%activeCenterPosition(1,i)=radiusComoving  *cos(theta)
          self%activeCenterPosition(2,i)=radiusComoving  *sin(theta)
          self%activeCenterPosition(3,i)=distanceComoving
       end do
    end if
    ! Check if we already have results for this specific node stored.
    if (node%uniqueID() /= self%activeUniqueID) then
       self%activeUniqueID=node%uniqueID()
       ! Calculate offset from the isolated node center.
       positionOffset=0.0d0
       nodeHost => node
       do while (nodeHost%isSatellite())
          satellite      =>  nodeHost      %satellite                          (                )
          positionOffset =  +positionOffset                                                       &
               &            +satellite                         %position       (                ) &
               &            /self          %cosmologyFunctions_%expansionFactor(basicHost%time())
          nodeHost       =>  nodeHost      %parent
       end do
       do i=1,self%activeCenterCount          
          self%activePosition(:,i)=self%activeCenterPosition(:,i)+positionOffset
       end do
       self%activeInCylinder=sum  (self%activePosition  (1:2,:)**2,dim=1) <= self%radiusCylinderComoving**2
       self%activeCount     =count(self%activeInCylinder                )
       j                    =0
       do i=1,self%activeCenterCount
          if (self%activeInCylinder(i)) then
             j                      =j+1
             self%activeInstances(j)=i
          end if
       end do
    end if
    return
  end subroutine cylindricalSampleNode
