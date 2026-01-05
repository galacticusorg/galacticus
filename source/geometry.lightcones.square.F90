!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{
  An implementation of the lightcone geometry class which assumes a square field of view.
  !!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass
  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Output_Times        , only : outputTimesClass

  !![
  <geometryLightcone name="geometryLightconeSquare">
   <description>
    A lightcone geometry class which assumes a square field of view., i.e. defined such that a point $(x,y,z)$ is in the survey
    angular mask if $|\hbox{atan2}(y,x)| &lt; \psi/2$ and $|\hbox{atan2}(z,x)| &lt; \psi/2$ where $\hbox{atan2}()$ is the
    quadrant-aware inverse tangent function, and $\psi$ is the angular size of the field, we compute the solid angle of the
    lightcone as follows. Define a spherical coordinate system $(\theta,\phi)$ with the pole ($\theta=0$) aligned with the
    $x$-axis. The solid angle of the field is then
    \begin{equation}
     \Omega = 2 \pi \int_0^{\psi/2} \sin\theta \mathrm{d}\theta + 8 \int_{\psi/2}^{\tan^{-1}(\sqrt{2}\tan(\psi/2))} \mathrm{d}\theta \sin\theta \int_{\cos^{-1}(\tan(\psi/2)/\tan\theta)}^{\pi/4} \mathrm{d}\phi,
    \end{equation}
    which is
    \begin{equation}
     \Omega = 2 \pi [1-\cos(\psi/2)] + 8 \int_{\psi/2}^{\tan^{-1}(\sqrt{2}\tan(\psi/2))} \mathrm{d}\theta \sin\theta \left[ {\pi\over 4} - \cos^{-1}\left({\tan(\psi/2)\over \tan\theta}\right)\right],
    \end{equation}
    or
    \begin{equation}
     \Omega = 2 \pi [1 - \cos(\tan^{-1}(\sqrt{2}\tan(\psi/2)))] - 8 \int_{\psi/2}^{\tan^{-1}(\sqrt{2}\tan(\psi/2))} \mathrm{d}\theta \sin\theta \cos^{-1}\left({\tan(\psi/2)\over \tan\theta}\right),
    \end{equation}
    The final integral can be evaluated (using Mathematica for example) to give
    \begin{eqnarray}
     \Omega &amp;=&amp; 2 \pi [3 - \cos(\tan^{-1}(\sqrt{2}\tan(\psi/2)))] - 8 \sin(x) \left( \sqrt{(a^2+1)\cos(2x)+a^2-1}(\log(a(\sqrt{2}\sqrt{2a^2\cos^2(x)+\cos(2x)-1} \right. \nonumber \\
    &amp; &amp;  +2a))-\log(\sqrt{\cos(2x)-1}))\sqrt{\csc^2(x)(-((a^2+1)\cos(2x)+a^2-1))}-\cot(x)((a^2+1)\cos(2x)+a^2-1) \nonumber \\
    &amp; &amp; \left. \cos^{-1}(a \cot(x)) \right) / [(a^2+1)\cos(2x)+a^2-1],
    \end{eqnarray}
    where $a=\tan(\psi/2)$ and $x=\tan^{-1}[\sqrt{2}\tan (\psi/2)]$.
  
    Various sub-parameters specify the details of the lightcone geometry. The {\normalfont \ttfamily lengthReplication} parameter
    should give the length of the simulation box (the box will be replicated to span the volume covered by the lightcone),
    with the {\normalfont \ttfamily lengthUnitsInSI} parameter giving the length unit in SI units and {\normalfont \ttfamily
    lengthHubbleExponent} giving the exponent of $h$ that appears in the length unit. The {\normalfont \ttfamily angularSize}
    parameter of {\normalfont \ttfamily fieldOfView} should gives the length of the side of the square field of view in
    degrees. The {\normalfont \ttfamily origin} element must contain the $x$, $y$, $z$ coordinates of the origin of the
    lightcone within the simulation box, while the {\normalfont \ttfamily unitVectorX} parameters must give unit vectors which
    point along the lightcone (for {\normalfont \ttfamily X}$=1$), and in the two directions perpendicular to the lightcone
    (for {\normalfont \ttfamily X}$=2$ and 3). The {\normalfont \ttfamily redshift} parameters must list the redshifts of
    available outputs.
   </description>
   <stateStore>
     <stateStore variables="nodeOperator_" store="nodeOperatorStateStore_" restore="nodeOperatorStateRestore_" module="Functions_Global"/>
   </stateStore>
   <deepCopy>
     <ignore variables="nodeOperator_"/>
   </deepCopy>
  </geometryLightcone>
  !!]
  type, extends(geometryLightconeClass) :: geometryLightconeSquare
     !!{
     A lightcone geometry class which assumes a square field of view.
     !!}
     private
     class           (cosmologyParametersClass), pointer                   :: cosmologyParameters_      => null()
     class           (cosmologyFunctionsClass ), pointer                   :: cosmologyFunctions_       => null()
     class           (outputTimesClass        ), pointer                   :: outputTimes_              => null()
     class           (*                       ), pointer                   :: nodeOperator_             => null()
     double precision                          , dimension(3,3)            :: unitVector
     double precision                          , dimension(3  )            :: origin                             , nodePositionCrossing   , &
          &                                                                   nodeVelocityCrossing               , unitVector1            , &
          &                                                                   unitVector2                        , unitVector3
     integer         (kind_int8               ), dimension(:), allocatable :: nodeIndicesReport
     double precision                          , dimension(:), allocatable :: distanceMinimum                    , distanceMaximum        , &
          &                                                                   timeMinimum_                       , timeMaximum_           , &
          &                                                                   outputTimes
     double precision                                                      :: lengthReplication                  , angularSize            , &
          &                                                                   solidAngleSubtended                , lengthUnitsInSI
     integer                                                               :: lengthHubbleExponent
     logical                                                               :: timeEvolvesAlongLightcone          , positionGettableChecked, &
          &                                                                   timeOfMergingGettableChecked
     integer         (kind_int8               )                            :: nodeUniqueIDCrossing
   contains
     !![
     <methods>
      <method method="positionAtOutput">
       <description>Returns the position of a point, {\normalfont \ttfamily nodePosition} (given in physical coordinates within the primary replicant volume), in comoving coordinates in the replicant volume in which it appears in the lightcone. If the point is \emph{not} in the lightcone the returned position is set to the largest possible negative number in each coordinate. If the optional {\normalfont \ttfamily positionFound} argument is given it will be set to true or false to indicate whether or not the point was found in the lightcone volume.</description>
      </method>
      <method method="replicants">
       <description>
        Performs various actions related to replicants of nodes appearing in lightcone output, depending on the value of the {\normalfont \ttfamily action} argument:
        \begin{description}
        \item[{\normalfont \ttfamily replicantActionCount}] returns in {\normalfont \ttfamily count} the number of replicants in which the node appears in the lightcone;
        \item[{\normalfont \ttfamily replicantActionAny}] returns true in {\normalfont \ttfamily isInLightcone} if the given position appears in \emph{any} replicant in the lightcone;
        \item[{\normalfont \ttfamily replicantActionInstance}] returns in {\normalfont \ttfamily position} the position in the {\normalfont \ttfamily instance}$^\mathrm{th}$ replicant in which this position appears in the lightcone.
        \end{description}
       </description>
      </method>
      <method method="periodicRange">
       <description>Computes the range of periodic replicants which could contribute to the lightcone in the given interval.</description>
      </method>
      <method method="nodePositionReplicant">
       <description>Computes the comoving position of a node in the specified replicant.</description>
      </method>
      <method method="nodeVelocityReplicant">
       <description>Computes the physical velocity of a node in the specified replicant.</description>
      </method>
      <method method="replicantLightConeCrossing">
	<description>Return the indices of the replicant for which the given node is crossing the lightcone.</description>
      </method>
      <method method="isInFieldOfView">
	<description>Return true if the given position is in the field of view.</description>
      </method>
     </methods>
     !!]
     final     ::                               squareDestructor
     procedure :: timeMinimum                => squareTimeMinimum
     procedure :: timeMaximum                => squareTimeMaximum
     procedure :: isInLightcone              => squareIsInLightcone
     procedure :: replicationCount           => squareReplicationCount
     procedure :: solidAngle                 => squareSolidAngle
     procedure :: position                   => squarePosition
     procedure :: velocity                   => squareVelocity
     procedure :: timeLightconeCrossing      => squareTimeLightconeCrossing
     procedure :: positionLightconeCrossing  => squarePositionLightconeCrossing
     procedure :: velocityLightconeCrossing  => squareVelocityLightconeCrossing
     procedure :: positionAtOutput           => squarePositionAtOutput
     procedure :: replicants                 => squareReplicants
     procedure :: periodicRange              => squarePeriodicRange
     procedure :: nodePositionReplicant      => squareNodePositionReplicant
     procedure :: nodeVelocityReplicant      => squareNodeVelocityReplicant
     procedure :: replicantLightConeCrossing => squareReplicantLightConeCrossing
     procedure :: isInFieldOfView            => squareIsInFieldOfView
     procedure :: deepCopy                   => squareDeepCopy
     procedure :: deepCopyReset              => squareDeepCopyReset
     procedure :: deepCopyFinalize           => squareDeepCopyFinalize
  end type geometryLightconeSquare

  interface geometryLightconeSquare
     !!{
     Constructors for the \refClass{geometryLightconeSquare} dark matter halo spin distribution class.
     !!}
     module procedure squareConstructorParameters
     module procedure squareConstructorInternal
  end interface geometryLightconeSquare

  ! Enumeration describing actions for the replicant method.
  !![
  <enumeration>
   <name>replicantAction</name>
   <description>Used to specify type of action required from the replicant method.</description>
   <visibility>private</visibility>
   <entry label="count"   />
   <entry label="any"     />
   <entry label="instance"/>
  </enumeration>
  !!]

  ! Sub-module scope variables used in root finding.
  class           (geometryLightconeSquare), pointer                   :: self_
  type            (treeNode               ), pointer                   :: node_
  integer                                                              :: i_       , j_          , &
       &                                                                  k_       , countSeek
  logical                                                              :: reporting
  double precision                         , allocatable, dimension(:) :: timeSeek , distanceSeek
  !$omp threadprivate(self_,node_,i_,j_,k_,timeSeek,distanceSeek,countSeek,reporting)
  
contains

  function squareConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{geometryLightconeSquare} lightcone geometry distribution class which takes a parameter list as
    input.
    !!}
    use :: Cosmology_Parameters            , only : cosmologyParameters   , cosmologyParametersClass, hubbleUnitsLittleH
    use :: Error                           , only : Error_Report
    use :: Input_Parameters                , only : inputParameter        , inputParameters
    use :: Numerical_Constants_Astronomical, only : degreesToRadians      , megaParsec
    use :: Functions_Global                , only : nodeOperatorConstruct_, nodeOperatorDestruct_
    implicit none
    type            (geometryLightconeSquare )                                :: self
    type            (inputParameters         ), intent(inout)                 :: parameters
    double precision                                         , dimension(3,3) :: unitVector
    double precision                                         , dimension(3  ) :: origin                   , unitVector1         , &
         &                                                                       unitVector2              , unitVector3
    integer         (kind_int8               ), allocatable  , dimension(:  ) :: nodeIndicesReport
    class           (cosmologyParametersClass), pointer                       :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer                       :: cosmologyFunctions_
    class           (outputTimesClass        ), pointer                       :: outputTimes_
    class           (*                       ), pointer                       :: nodeOperator_
    double precision                                                          :: lengthReplication        , angularSize         , &
         &                                                                       lengthUnitsInSI          , unitConversionLength
    integer                                                                   :: lengthHubbleExponent
    logical                                                                   :: timeEvolvesAlongLightcone

    allocate(nodeIndicesReport(parameters%count('nodeIndicesReport',zeroIfNotPresent=.true.)))
    if (parameters%isPresent('nodeIndicesReport')) then
       !![
       <inputParameter>
	 <name>nodeIndicesReport</name>
	 <description>A list of node indices for which reporting should be performed.</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
    end if
    !![
    <inputParameter>
      <name>origin</name>
      <source>parameters</source>
      <variable>origin</variable>
      <description>The origin for the lightcone.</description>
    </inputParameter>
    <inputParameter>
      <name>unitVector1</name>
      <source>parameters</source>
      <description>The first (radial) unit vector defining the lightcone geometry.</description>
    </inputParameter>
    !!]
    unitVector(:,1)=unitVector1
    !![
    <inputParameter>
      <name>unitVector2</name>
      <source>parameters</source>
      <description>The second (angular) unit vector defining the lightcone geometry.</description>
    </inputParameter>
    !!]
    unitVector(:,2)=unitVector2
    !![
    <inputParameter>
      <name>unitVector3</name>
      <source>parameters</source>
      <description>The third (angular) unit vector defining the lightcone geometry.</description>
    </inputParameter>
    !!]
    unitVector(:,3)=unitVector3
    !![
    <inputParameter>
      <name>lengthReplication</name>
      <source>parameters</source>
      <variable>lengthReplication</variable>
      <description>The length of the simulation box being used to construct the lightcone.</description>
    </inputParameter>
    <inputParameter>
      <name>lengthUnitsInSI</name>
      <source>parameters</source>
      <variable>lengthUnitsInSI</variable>
      <description>The units of the box length in the SI system.</description>
    </inputParameter>
    <inputParameter>
      <name>lengthHubbleExponent</name>
      <source>parameters</source>
      <variable>lengthHubbleExponent</variable>
      <description>The exponent of the ``little-$h$'' parameter used in the definition of the box length.</description>
    </inputParameter>
    <inputParameter>
      <name>angularSize</name>
      <source>parameters</source>
      <variable>angularSize</variable>
      <description>The angular size (i.e. side length) of the square field of view of the lightcone (in units of degrees).</description>
    </inputParameter>
    <inputParameter>
      <name>timeEvolvesAlongLightcone</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <variable>timeEvolvesAlongLightcone</variable>
      <description>If {\normalfont \ttfamily true}, cosmic time evolves along the lightcone as expected. Otherwise, time is fixed at the present epoch throughout the lightcone. This allows construction of lightcones with no evolution.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="outputTimes"         name="outputTimes_"         source="parameters"/>
    !!]
    call nodeOperatorConstruct_(parameters,nodeOperator_)
    ! Convert angle to radians.
    angularSize=angularSize*degreesToRadians
    ! Convert lengths units internal units.
    unitConversionLength=+lengthUnitsInSI                                                               &
         &               *cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**lengthHubbleExponent &
         &               /megaParsec
    origin              =origin           *unitConversionLength
    lengthReplication   =lengthReplication*unitConversionLength
    ! Construct the object.
    self                     =geometryLightconeSquare(origin,unitVector,angularSize,lengthReplication,timeEvolvesAlongLightcone,nodeIndicesReport,cosmologyParameters_,cosmologyFunctions_,outputTimes_,nodeOperator_)
    self%lengthUnitsInSI     =lengthUnitsInSI
    self%lengthHubbleExponent=lengthHubbleExponent
    !![
    <inputParametersValidate source="parameters" extraAllowedNames="nodeOperator"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="outputTimes_"        />
    !!]
    if (associated(nodeOperator_)) call nodeOperatorDestruct_(nodeOperator_)
    return
  end function squareConstructorParameters

  function squareConstructorInternal(origin,unitVector,angularSize,lengthReplication,timeEvolvesAlongLightcone,nodeIndicesReport,cosmologyParameters_,cosmologyFunctions_,outputTimes_,nodeOperator_) result(self)
    !!{
    Internal constructor for the \refClass{geometryLightconeSquare} lightcone geometry distribution class.
    !!}
    use :: Error                           , only : Error_Report
    use :: ISO_Varying_String              , only : var_str         , varying_string
    use :: Numerical_Constants_Math        , only : Pi              , e
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: String_Handling                 , only : operator(//)
    use :: Vectors                         , only : Vector_Magnitude
    implicit none
    type            (geometryLightconeSquare )                                :: self
    double precision                          , dimension(3,3), intent(in   ) :: unitVector
    double precision                          , dimension(3)  , intent(in   ) :: origin
    class           (cosmologyParametersClass), target        , intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), target        , intent(in   ) :: cosmologyFunctions_
    class           (outputTimesClass        ), target        , intent(in   ) :: outputTimes_
    class           (*                       ), target        , intent(in   ) :: nodeOperator_
    double precision                                          , intent(in   ) :: lengthReplication                 , angularSize
    logical                                                   , intent(in   ) :: timeEvolvesAlongLightcone
    integer         (kind_int8               ), dimension(:)  , intent(in   ) :: nodeIndicesReport
    double precision                          , parameter                     :: orthogonalityTolerance   =1.000d-6
    double precision                          , parameter                     :: angularSizeSmall         =0.017d+0
    integer                                                                   :: i                                 , j
    integer         (c_size_t                )                                :: iOutput
    double precision                                                          :: timeMinimum                       , timeMaximum    , &
         &                                                                       tanHalfAngle                      , outputTime     , &
         &                                                                       distanceMinimum                   , distanceMaximum, &
         &                                                                       magnitude
    character       (len=12                  )                                :: label
    type            (varying_string          )                                :: message
    !![
    <constructorAssign variables="origin, unitVector, angularSize, lengthReplication, timeEvolvesAlongLightcone, nodeIndicesReport, *cosmologyParameters_, *cosmologyFunctions_, *outputTimes_, *nodeOperator_"/>
    !!]

    ! Store unit vectors.
    self%unitVector1=unitVector(:,1)
    self%unitVector2=unitVector(:,2)
    self%unitVector3=unitVector(:,3)
    ! Extract times from the outputTimes object.
    allocate(self%outputTimes(self%outputTimes_%count()))
    do iOutput=1_c_size_t,self%outputTimes_%count()
       self%outputTimes(iOutput)=self%outputTimes_%time(iOutput)
    end do
    ! Find the minimum and maximum distance associated with each output time.
    allocate(self%distanceMinimum(size(self%outputTimes)))
    allocate(self%distanceMaximum(size(self%outputTimes)))
    allocate(self%timeMinimum_   (size(self%outputTimes)))
    allocate(self%timeMaximum_   (size(self%outputTimes)))
    do iOutput=1,size(self%outputTimes)
       if (iOutput == 1                     ) then
          timeMinimum=                                 self%outputTimes(iOutput)
       else
          timeMinimum=sqrt(self%outputTimes(iOutput-1)*self%outputTimes(iOutput))
       end if
       if (iOutput == size(self%outputTimes)) then
          timeMaximum=                                 self%outputTimes(iOutput)
       else
          timeMaximum=sqrt(self%outputTimes(iOutput+1)*self%outputTimes(iOutput))
       end if
       self%timeMinimum_   (iOutput)=                                          timeMinimum
       self%timeMaximum_   (iOutput)=                                          timeMaximum
       self%distanceMinimum(iOutput)=self%cosmologyFunctions_%distanceComoving(timeMaximum)
       self%distanceMaximum(iOutput)=self%cosmologyFunctions_%distanceComoving(timeMinimum)
    end do
    ! Normalize unit vectors.
    do i=1,3
       if (Vector_Magnitude(self%unitVector(:,i)) == 0.0d0) &
            & call Error_Report('null unit vector is not permitted'//{introspection:location})
       self%unitVector(:,i)=+                 self%unitVector(:,i)  &
            &               /Vector_Magnitude(self%unitVector(:,i))
    end do
    ! Test for orthogonality.
    do i=1,2
       do j=i+1,3
          magnitude=abs(Dot_Product(self%unitVector(:,i),self%unitVector(:,j)))
          if (magnitude > orthogonalityTolerance) then
             write (label,'(e12.6)') magnitude
             message=var_str('unit vectors ')//i//' and '//j//' are not orthogonal (|êᵢ⨯êⱼ|='//label//')'
             call Error_Report(message//{introspection:location})
          end if
       end do
    end do
    ! If requested, reduce the list of lightcone output times to the final time, and have it incorporate the full radial length of
    ! the lightcone. This allows the construction of lightcones with no evolution along the cone.
    if (.not.timeEvolvesAlongLightcone) then
       outputTime     =self%outputTimes    (size(self%outputTimes))
       distanceMinimum=self%distanceMinimum(size(self%outputTimes))
       distanceMaximum=self%distanceMaximum(                    1 )
       timeMaximum    =self%timeMinimum_   (size(self%outputTimes))
       timeMinimum    =self%timeMaximum_   (                    1 )
       deallocate(self%outputTimes        )
       deallocate(self%distanceMinimum    )
       deallocate(self%distanceMaximum    )
       allocate(self%outputTimes    (1))
       allocate(self%distanceMinimum(1))
       allocate(self%distanceMaximum(1))
       self%outputTimes    =outputTime
       self%timeMinimum_   =timeMinimum
       self%timeMaximum_   =timeMaximum
       self%distanceMinimum=distanceMinimum
       self%distanceMaximum=distanceMaximum
    end if
    ! Compute the solid angle of the lightcone.
    if (angularSize < angularSizeSmall) then
       self%solidAngleSubtended=angularSize**2
    else
       tanHalfAngle=tan(angularSize/2.0d0)
       self%solidAngleSubtended=+2.0d0                                     &
            &                   *Pi                                        &
            &                   *(                                         &
            &                     +3.0d0                                   &
            &                     -cos(                                    &
            &                          atan(                               &
            &                                sqrt(2.0d0)                   &
            &                               *tanHalfAngle                  &
            &                              )                               &
            &                         )                                    &
            &                    )                                         &
            &                   -8.0d0                                     &
            &                   *inverseCosineIntegral(                    &
            &                                          tanHalfAngle      , &
            &                                          atan(               &
            &                                                sqrt(2.0d0)   &
            &                                               *tanHalfAngle  &
            &                                              )               &
            &                                         )
    end if
    ! Set property attribute check status.
    self%positionGettableChecked     =.false.
    self%timeOfMergingGettableChecked=.false.
    ! Initialize lightcone crossing state.
    self%nodeUniqueIDCrossing=-1_kind_int8
    self%nodePositionCrossing=0.0d0
    self%nodeVelocityCrossing=0.0d0
    ! Set unit conversion rates (to defaults).
    self%lengthUnitsInSI     =megaParsec
    self%lengthHubbleExponent=0
    return

  contains

    double precision function inverseCosineIntegral(a,x)
      !!{
      Integral of $\sin(x)*\cos^{-1}[a/\tan(x)]$ evaluated using Wolfram Alpha.
      !!}
      use :: Trigonometric_Functions, only : cosec, cot
      implicit none
      double precision, intent(in) :: a , x
      double complex               :: aa, xx

      aa                   =a
      xx                   =x
      inverseCosineIntegral=real(                                                                                     &
           &                     +sin(xx)                                                                             &
           &                     *(                                                                                   &
           &                       +sqrt(                (aa**2+1.0d0)*cos(2.0d0*xx)+aa**2-1.0d0)                     &
           &                       *sqrt(cosec(xx)**2*(-((aa**2+1.0d0)*cos(2.0d0*xx)+aa**2-1.0d0)))                   &
           &                       *(                                                                                 &
           &                         +log(aa*(sqrt(2.0d0)*sqrt(2.0d0*aa**2*cos(xx)**2+cos(2.0d0*xx)-1.0d0)+2.0d0*aa)) &
           &                         -log(                sqrt(                       cos(2.0d0*xx)-1.0d0)          ) &
           &                        )                                                                                 &
           &                       -     cot  (xx)   *(+ (aa**2+1.0d0)*cos(2.0d0*xx)+aa**2-1.0d0)*acos(aa*cot(xx))    &
           &                      )                                                                                   &
           &                     /                    (+ (aa**2+1.0d0)*cos(2.0d0*xx)+aa**2-1.0d0)                     &
           &                    )
      return
    end function inverseCosineIntegral

  end function squareConstructorInternal

  subroutine squareDestructor(self)
    !!{
    Destructor for the \refClass{geometryLightconeSquare} lightcone geometry distribution class.
    !!}
    use :: Functions_Global, only : nodeOperatorDestruct_
    implicit none
    type(geometryLightconeSquare), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%outputTimes_"        />
    !!]
    if (associated(self%nodeOperator_)) call nodeOperatorDestruct_(self%nodeOperator_)
    return
  end subroutine squareDestructor

  function squareReplicationCount(self,node)
    !!{
    Determine the number of times {\normalfont \ttfamily node} appears in the lightcone.
    !!}
    use            :: Arrays_Search   , only : searchArrayClosest
    use            :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentPosition, treeNode
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    implicit none
    integer(c_size_t               )                :: squareReplicationCount
    class  (geometryLightconeSquare), intent(inout) :: self
    type   (treeNode               ), intent(inout) :: node
    class  (nodeComponentBasic     ), pointer       :: basic
    class  (nodeComponentPosition  ), pointer       :: position
    integer(c_size_t               )                :: output

    basic    => node%basic   ()
    position => node%position()
    output   =  searchArrayClosest(self%outputTimes,basic%time())
    call self%replicants(output,position%position(),replicantActionCount,count=squareReplicationCount)
    return
  end function squareReplicationCount

  double precision function squareTimeMinimum(self)
    !!{
    Return the minimum time in the lightcone.
    !!}
    implicit none
    class(geometryLightconeSquare), intent(inout) :: self

    squareTimeMinimum=self%outputTimes(1)
    return
  end function squareTimeMinimum

  double precision function squareTimeMaximum(self)
    !!{
    Return the minimum time in the lightcone.
    !!}
    implicit none
    class(geometryLightconeSquare), intent(inout) :: self

    squareTimeMaximum=self%outputTimes(size(self%outputTimes))
    return
  end function squareTimeMaximum

  logical function squareIsInLightcone(self,node,atPresentEpoch,radiusBuffer)
    !!{
    Determine if the given {\normalfont \ttfamily node} lies within the lightcone. Note that, when called with {\normalfont
    \ttfamily atPresentEpoch=false} this function returns true if the node is in the lightcone at any point during its
    existence. However, this check is made assuming that each node remains at fixed comoving coordinates between each output
    time---there is no consideration of movement between output times. It is therefore recommended that some buffer is added to
    catch any nodes which may briefly enter the lightcone between output times.
    !!}
    use            :: Arrays_Search       , only : searchArray
    use            :: Error               , only : Component_List          , Error_Report
    use            :: Galacticus_Nodes    , only : defaultPositionComponent, defaultSatelliteComponent, nodeComponentBasic , nodeComponentPosition, &
          &                                        nodeComponentSatellite  , treeNode
    use, intrinsic :: ISO_C_Binding       , only : c_size_t
    use            :: ISO_Varying_String  , only : varying_string
    use            :: Numerical_Comparison, only : Values_Agree
    use            :: String_Handling     , only : operator(//)
    implicit none
    class           (geometryLightconeSquare), intent(inout)               :: self
    type            (treeNode               ), intent(inout)               :: node
    logical                                  , intent(in   ) , optional    :: atPresentEpoch
    double precision                         , intent(in   ) , optional    :: radiusBuffer
    class           (nodeComponentBasic     )                , pointer     :: basic               , basicParent
    class           (nodeComponentPosition  )                , pointer     :: position
    class           (nodeComponentSatellite )                , pointer     :: satellite
    double precision                         , parameter                   :: timeTolerance=1.0d-3
    double precision                         , dimension(3  )              :: nodePosition
    double precision                         , dimension(:,:), allocatable :: nodePositionHistory
    integer         (c_size_t               )                                 outputMinimum       , outputMaximum, &
         &                                                                    output
    double precision                                                       :: timeCurrent         , timeMerge    , &
         &                                                                    timeFinal
    character       (len=10                 )                              :: label
    type            (varying_string         )                              :: message
    !![
    <optionalArgument name="atPresentEpoch" defaultsTo=".true." />
    !!]

    ! Get the basic component.
    basic => node%basic()
    ! Assume not in the lightcone by default.
    squareIsInLightcone=.false.
    ! If no output time exists after our node, so it can not be in the lightcone.
    if (self%timeMaximum_(size(self%timeMaximum_)) < basic%time()) return
    ! Check that we can get the position of a node.
    if (.not.self%positionGettableChecked) then
       self%positionGettableChecked=.true.
       if (.not.defaultPositionComponent%positionIsGettable())                                                                          &
            & call Error_Report                                                                                                         &
            &      (                                                                                                                    &
            &       'testing if a node is in the lightcone requires that the position property of the position component be gettable'// &
            &       Component_List(                                                                                                     &
            &                      'position'                                                                                        ,  &
            &                       defaultPositionComponent%positionAttributeMatch(requireGettable=.true.)                             &
            &                     )                                                                                                  // &
            &       {introspection:location}                                                                                            &
            &      )
    end if
    ! Determine to which output this galaxy corresponds.
    if (basic%time() > self%timeMinimum_(size(self%timeMinimum_))) then
       outputMinimum=size(self%timeMinimum_)
    else
       outputMinimum=searchArray(self%timeMinimum_,basic%time())
    end if
    if (atPresentEpoch_) then
       ! We want to check only the current time for this node. Check that the node exists precisely at a lightcone snapshot time,
       ! and then set the maximum output to check to equal to minimum, such that we test only the current time.
       if (.not.Values_Agree(self%outputTimes(outputMinimum),basic%time(),relTol=timeTolerance)) then
          message=         'failed to find matching time in lightcone'                       //char(10)
          write (label,'(f6.3)') basic%time()
          message=message//'               node time = '//trim(label)//' Gyr'                //char(10)
          write (label,'(f6.3)') self%outputTimes(outputMinimum)
          message=message//'  closest lightcone time = '//trim(label)//' Gyr ['//outputMinimum//']'
          call Error_Report(message//{introspection:location})
       end if
       outputMaximum=outputMinimum
    else
       ! We need to test if the node is in the lightcone at any point during its existence. Determine for how long this node will
       ! persist. This requires that the merging time be defined. To perform the test also requires that a position exists.  Check
       ! that we can get the time of merging and the position of a node.
       if (.not.self%timeOfMergingGettableChecked) then
          self%timeOfMergingGettableChecked=.true.
          if (.not.defaultSatelliteComponent%timeOfMergingIsGettable())                                                                               &
               & call Error_Report                                                                                                                    &
               &      (                                                                                                                               &
               &       'testing if a node is ever in the lightcone requires that the timeOfMerging property of the satellite component be gettable'// &
               &       Component_List(                                                                                                                &
               &                      'satellite'                                                                                                   , &
               &                       defaultSatelliteComponent%timeOfMergingAttributeMatch(requireGettable=.true.)                                  &
               &                     )                                                                                                             // &
               &       {introspection:location}                                                                                                       &
               &      )
          if (.not.defaultPositionComponent%positionIsGettable())                                                                                     &
               & call Error_Report                                                                                                                    &
               &      (                                                                                                                               &
               &       'testing if a node is ever in the lightcone requires that the position property of the position component be gettable'      // &
               &       Component_List(                                                                                                                &
               &                      'position'                                                                                                    , &
               &                       defaultPositionComponent %positionAttributeMatch     (requireGettable=.true.)                                  &
               &                     )                                                                                                             // &
               &       {introspection:location}                                                                                                       &
               &      )
       end if
       ! Get node position component.
       position => node%position()
       ! Test for primary progenitor status.
       allocate(nodePositionHistory(0,0))
       if (node%isPrimaryProgenitor().or..not.associated(node%parent)) then
          ! Node is the primary progenitor (or is the root node in its tree), so exists only until the time of its parent and has
          ! fixed position.
          if (associated(node%parent)) then
             basicParent => node%parent%basic()
             timeFinal   =  basicParent%time ()
          else
             timeFinal   =  basic      %time ()
          end if
          if (timeFinal > self%timeMinimum_(size(self%timeMinimum_))) then
             outputMaximum=size(self%timeMinimum_)
          else
             outputMaximum=searchArray(self%timeMinimum_,timeFinal)
          end if
          ! If the earliest output time is after the parent time, so the node can not be in the lightcone.
          if (self%timeMinimum_(1) > timeFinal) return
          ! Construct array of positions at output times.
          deallocate(nodePositionHistory)
          allocate(nodePositionHistory(3_c_size_t,outputMaximum-outputMinimum+1))
          do output=1,outputMaximum-outputMinimum+1
             nodePositionHistory(:,output)=position%position()
          end do
       else
          ! Node is not the primary progenitor, so will become a satellite. We therefore need to know for how long it will persist
          ! as a satellite and its position during that time.
          satellite => node     %satellite    ()
          timeMerge =  satellite%timeOfMerging()
          if (timeMerge < basic%time()) call Error_Report('can not determine if satellite is in lightcone without knowledge of its time of merging'//{introspection:location})
          outputMaximum=searchArray(self%timeMinimum_,timeMerge)
          ! Construct array of positions at output times.
          deallocate(nodePositionHistory)
          allocate(nodePositionHistory(3_c_size_t,outputMaximum-outputMinimum+1))
          timeCurrent=basic%time()
          do output=1,outputMaximum-outputMinimum+1
             call basic%timeSet(self%outputTimes(output+outputMinimum-1))
             nodePositionHistory(:,output)=position%position()
          end do
          call basic%timeSet(timeCurrent)
       end if
    end if
    ! Iterate over possible output times.
    do output=outputMinimum,outputMaximum
       ! Get position of galaxy in physical coordinates.
       if (atPresentEpoch_) then
          position     => node    %position()
          nodePosition =  position%position()
       else
          nodePosition =  nodePositionHistory(:,output-outputMinimum+1)
       end if
       ! Find the position (if possible) of the node in the lightcone for this output.
       call self%replicants(output,nodePosition,replicantActionAny,isInLightcone=squareIsInLightcone,radiusBuffer=radiusBuffer)
       ! Return as soon as we find that the node is in the lightcone.
       if (squareIsInLightcone) return
    end do
    return
  end function squareIsInLightcone

  double precision function squareSolidAngle(self)
    !!{
    Return the solid angle (in steradians) of a square lightcone.
    !!}
    implicit none
    class(geometryLightconeSquare), intent(inout) :: self

    squareSolidAngle=self%solidAngleSubtended
    return
  end function squareSolidAngle

  function squarePosition(self,node,instance)
    !!{
    Return the position of the node in lightcone coordinates.
    !!}
    use            :: Arrays_Search       , only : searchArrayClosest
    use            :: Error               , only : Error_Report
    use            :: Galacticus_Nodes    , only : nodeComponentBasic, nodeComponentPosition, treeNode
    use, intrinsic :: ISO_C_Binding       , only : c_size_t
    use            :: ISO_Varying_String  , only : varying_string
    use            :: Numerical_Comparison, only : Values_Agree
    use            :: String_Handling     , only : operator(//)
    implicit none
    double precision                         , dimension(3)          :: squarePosition
    class           (geometryLightconeSquare), intent(inout)         :: self
    type            (treeNode               ), intent(inout), target :: node
    integer         (c_size_t               ), intent(in   )         :: instance
    class           (nodeComponentBasic     ), pointer               :: basic
    class           (nodeComponentPosition  ), pointer               :: position
    double precision                         , parameter             :: timeTolerance =1.0d-3
    integer         (c_size_t               )                        :: output
    character       (len=10                 )                        :: label
    type            (varying_string         )                        :: message

    ! Get the basic component.
    basic    => node%basic   ()
    position => node%position()
    ! Determine to which output this node corresponds.
    output=searchArrayClosest(self%outputTimes,basic%time())
    if (.not.Values_Agree(self%outputTimes(output),basic%time(),relTol=timeTolerance)) then
       message=         'failed to find matching time in lightcone'                       //char(10)
       write (label,'(f6.3)') basic%time()
       message=message//'               node time = '//trim(label)//' Gyr'                //char(10)
       write (label,'(f6.3)') self%outputTimes(output)
       message=message//'  closest lightcone time = '//trim(label)//' Gyr ['//output//']'
       call Error_Report(message//{introspection:location})
    end if
    squarePosition=self%positionAtOutput(output,position%position(),instance)
    return
  end function squarePosition

  function squareVelocity(self,node,instance)
    !!{
    Return the velocity of the node in lightcone coordinates.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentPosition, treeNode
    implicit none
    double precision                         , dimension(3)  :: squarevelocity
    class           (geometryLightconeSquare), intent(inout) :: self
    type            (treeNode               ), intent(inout) :: node
    integer         (c_size_t               ), intent(in   ) :: instance
    class           (nodeComponentPosition  ), pointer       :: position
    integer                                                  :: i
    !$GLC attributes unused :: instance

    ! Get the position component.
    position => node%position()
    ! Compute velocity of galaxy in lightcone coordinate system.
    do i=1,3
       squareVelocity(i)=Dot_Product(position%velocity(),self%unitVector(:,i))
    end do
    return
  end function squareVelocity

  subroutine squareReplicants(self,output,nodePosition,action,count,isInLightcone,radiusBuffer,instance,position)
    !!{
    Compute quantities related to the number of replicants in which a node appears.
    !!}
    use            :: Error        , only : Error_Report
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: Vectors      , only : Vector_Magnitude
    implicit none
    class           (geometryLightconeSquare       ), intent(inout)                           :: self
    integer         (c_size_t                      ), intent(in   )                           :: output
    double precision                                , intent(in   ), dimension(3  )           :: nodePosition
    type            (enumerationReplicantActionType), intent(in   )                           :: action
    integer         (c_size_t                      ), intent(  out)                , optional :: count
    integer         (c_size_t                      ), intent(in   )                , optional :: instance
    logical                                         , intent(  out)                , optional :: isInLightcone
    double precision                                , intent(in   )                , optional :: radiusBuffer
    double precision                                , intent(  out), dimension(3  ), optional :: position
    double precision                                               , dimension(3  )           :: nodePositionComoving, nodePositionReplicant, &
         &                                                                                       origin
    integer                                                        , dimension(3,2)           :: periodicRange
    integer                                                                                   :: i                   , j                    , &
         &                                                                                       k                   , iAxis                , &
         &                                                                                       replicantCount
    logical                                                                                   :: isInFieldOfView     , found
    double precision                                                                          :: distanceMinimum     , distanceMaximum      , &
         &                                                                                       distanceRadial

    ! Validate input.
    select case (action%ID)
    case (replicantActionCount   %ID)
       if (.not.present(count        )) call Error_Report('count argument not provided'        //{introspection:location})
    case (replicantActionAny     %ID)
       if (.not.present(isInLightcone)) call Error_Report('isInLightcone argument not provided'//{introspection:location})
    case (replicantActionInstance%ID)
       if (.not.present(instance     )) call Error_Report('instance argument not provided'     //{introspection:location})
       if (.not.present(position     )) call Error_Report('position argument not provided'     //{introspection:location})
    case default
       call Error_Report('unknown action'//{introspection:location})
    end select
    ! Determine range of possible replicants of this galaxy which could be in the lightcone.
    periodicRange=self%periodicRange(self%distanceMinimum(output),self%distanceMaximum(output),radiusBuffer,origin,distanceMinimum,distanceMaximum)
    ! Get comoving position.
    nodePositionComoving=nodePosition/self%cosmologyFunctions_%expansionFactor(self%outputTimes(output))
    ! Iterate over all replicants.
    found         =.false.
    replicantCount=0
    do i=periodicRange(1,1),periodicRange(1,2)
       do j=periodicRange(2,1),periodicRange(2,2)
          do k=periodicRange(3,1),periodicRange(3,2)
             ! Compute position of node in lightcone coordinate system.
             forall(iAxis=1:3)
                nodePositionReplicant(iAxis)=Dot_Product(nodePositionComoving-origin+self%lengthReplication*dble([i,j,k]),self%unitVector(:,iAxis))
             end forall
             ! Test if node lies within the correct angular window.
             isInFieldOfView=self%isInFieldOfView(nodePositionReplicant)
             ! Test if node lies within appropriate radial range.
             distanceRadial=Vector_Magnitude(nodePositionReplicant)
             found         =      isInFieldOfView                    &
                  &         .and. distanceRadial  >  distanceMinimum &
                  &         .and. distanceRadial  <= distanceMaximum
             if (found) then
                if (action == replicantActionAny) then
                   ! If we are just testing for presence in the lightcone, then we can return right now.
                   isInLightcone=.true.
                   return
                end if
                replicantCount=replicantCount+1
                if (action == replicantActionInstance .and. replicantCount == instance) then
                   position=nodePositionReplicant
                   return
                end if
             end if
          end do
       end do
    end do
    ! Set output.
    select case (action%ID)
    case (replicantActionCount   %ID)
       count        =replicantCount
    case (replicantActionAny     %ID)
       isInLightcone=.false.
    case (replicantActionInstance%ID)
       call Error_Report('instance not found'//{introspection:location})
    case default
       call Error_Report('unknown action'    //{introspection:location})
    end select
    return
  end subroutine squareReplicants

  function squarePositionAtOutput(self,output,nodePosition,instance)
    !!{
    Return the position of the node in lightcone coordinates.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    double precision                                        , dimension(3) :: squarePositionAtOutput
    class           (geometryLightconeSquare), intent(inout)               :: self
    integer         (c_size_t               ), intent(in   )               :: output
    double precision                         , intent(in   ), dimension(3) :: nodePosition
    integer         (c_size_t               ), intent(in   )               :: instance

    call self%replicants(output,nodePosition,replicantActionInstance,instance=instance,position=squarePositionAtOutput)
    return
  end function squarePositionAtOutput

  double precision function squareTimeLightconeCrossing(self,node,timeStart,timeEnd,timesCrossing)
    !!{
    Return the time of the next lightcone crossing for this node.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic                      , nodeComponentPosition
    use :: Display                         , only : displayIndent                           , displayUnindent      , displayMessage
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr                       , gigaYear             , megaParsec
    use :: Numerical_Constants_Physical    , only : speedLight
    use :: Vectors                         , only : Vector_Magnitude
    use :: Root_Finder                     , only : rootFinder
    use :: Sorting                         , only : sort
    use :: Functions_Global                , only : nodeOperatorPredeterminedSolveAnalytics_
    implicit none
    class           (geometryLightconeSquare), intent(inout)                , target                :: self
    type            (treeNode               ), intent(inout)                , target                :: node
    double precision                         , intent(in   )                                        :: timeStart                   , timeEnd
    double precision                         , intent(inout), dimension(:  ), allocatable, optional :: timesCrossing
    double precision                                        , dimension(:  ), allocatable           :: timesCrossingTmp
    integer                                                 , dimension(3,2)                        :: periodicRange
    class           (nodeComponentBasic     )                               , pointer               :: basic
    class           (nodeComponentPosition  )                               , pointer               :: position
    double precision                                        , dimension(3  )                        :: positionReference           , nodePositionStart , &
         &                                                                                             nodePositionEnd
    double precision                         , parameter                                            :: timeSeekStep         =1.0d-3 ! Timestep used in reporting of crossing times.
    double precision                         , parameter                                            :: speedMaximum         =2000.0 ! Maximum plausible physical speed for any node.
    double precision                                                                                :: distanceMinimum             , distanceMaximum   , &
         &                                                                                             distanceNodeStart           , distanceNodeEnd   , &
         &                                                                                             radiusBuffer                , timeCrossing      , &
         &                                                                                             timeOriginal                , timeStart_        , &
         &                                                                                             distanceCrossingSeek
    logical                                                                                         :: isInFieldOfView
    integer                                                                                         :: i                           , j                 , &
         &                                                                                             k                           , countSteps        , &
         &                                                                                             iSeek
    type            (rootFinder             ), save                                                 :: finder
    logical                                  , save                                                 :: finderConstructed   =.false.
    !$omp threadprivate(finder,finderConstructed)
    character       (len= 64                )                                                       :: label

    reporting=any(node%index() == self%nodeIndicesReport)
    if (reporting) then
       write (label,'(i12)') node%index()
       call displayIndent('Lightcone crossing report for node: '//trim(adjustl(label)))
    end if
    basic                       => node    %basic                               (          )
    position                    => node    %position                            (          ) 
    positionReference           =  position%position                            (          )
    timeOriginal                =  basic   %time                                (          )
    timeStart_                  =  max(timeStart,timeOriginal)
    distanceMinimum             =  self    %cosmologyFunctions_%distanceComoving(timeEnd   )
    distanceMaximum             =  self    %cosmologyFunctions_%distanceComoving(timeStart_)
    radiusBuffer                =  +(                       &
         &                           +timeEnd               &
         &                           -timeStart_            &
         &                          )                       &
         &                         *speedMaximum            &
         &                         /MpcPerKmPerSToGyr
    ! Find the range of replicants in which this node might cross the lightcone.
    periodicRange=self%periodicRange(distanceMinimum,distanceMaximum,radiusBuffer)
    if (reporting) then
       write (label,'(3(e12.6,1x))') positionReference
       call displayMessage('Position reference: '//trim(adjustl(label)))
       write (label,'(1(e12.6,1x))') timeOriginal
       call displayMessage('Time original:      '//trim(adjustl(label)))
       write (label,'(1(e12.6,1x))') timeStart
       call displayMessage('Time start:         '//trim(adjustl(label)))
       write (label,'(2(e12.6,1x))') distanceMinimum,distanceMaximum
       call displayMessage('Distance min/max:   '//trim(adjustl(label)))
       write (label,'(1(e12.6,1x))') radiusBuffer
       call displayMessage('Buffer radius:      '//trim(adjustl(label)))
       write (label,'(6(i4,1x))') periodicRange
       call displayMessage('Periodic range:     '//trim(adjustl(label)))
    end if
    ! Iterate over replicants of interest.
    if (.not.finderConstructed) then 
       finder           =rootFinder(rootFunction=timeCrossingRoot,toleranceRelative=1.0d-9)
       finderConstructed=.true.
    end if
    self_ => self
    node_ => node
    squareTimeLightconeCrossing=huge(0.0d0)
    do i=periodicRange(1,1),periodicRange(1,2)
       do j=periodicRange(2,1),periodicRange(2,2)
          do k=periodicRange(3,1),periodicRange(3,2)
             if (reporting) then
                write (label,'(3(i4,1x))') i,j,k
                call displayIndent("Periodic replicant: "//trim(adjustl(label)))
             end if
             ! Compute position of node in lightcone coordinate system.
             nodePositionStart=self%nodePositionReplicant(node             ,timeStart_,self%origin,[i,j,k],setTime=.true.,report=reporting)
             nodePositionEnd  =self%nodePositionReplicant(node             ,timeEnd   ,self%origin,[i,j,k],setTime=.true.,report=reporting)
             distanceNodeStart=Vector_Magnitude          (nodePositionStart                                              )
             distanceNodeEnd  =Vector_Magnitude          (nodePositionEnd                                                )
             if (reporting) then
                write (label,'(3(e12.6,1x),"| ",e12.6)') nodePositionStart,distanceNodeStart
                call displayMessage("Postion start: "//trim(adjustl(label)))
                write (label,'(3(e12.6,1x),"| ",e12.6)') nodePositionEnd  ,distanceNodeEnd
                call displayMessage("Postion end:   "//trim(adjustl(label)))
             end if
             if     (                                      &
                  &   distanceNodeStart <  distanceMaximum &
                  &  .and.                                 &
                  &   distanceNodeEnd   >= distanceMinimum &
                  & ) then
                ! Find the precise time of lightcone crossing.
                i_          =i
                j_          =j
                k_          =k
                countSeek   =0
                timeCrossing=finder%find(rootRange=[timeStart_,timeEnd])
                if (reporting) then
                   call displayIndent("Lightcone crossing seek:")
                   call displayMessage("time         distance root")
                   call displayMessage("__________________________")
                   call sort(timeSeek(1:countSeek),distanceSeek(1:countSeek))
                   do iSeek=1,countSeek
                      write (label,'(e12.6,1x,e12.6)') timeSeek(iSeek),distanceSeek(iSeek)
                      call displayMessage(trim(adjustl(label)))
                   end do
                   call displayUnindent("")
                   deallocate(timeSeek    )
                   deallocate(distanceSeek)
                end if
                ! Check that the node is in the field of view at this time, that this is the earliest crossing, and that the
                ! crossing occurs at least some small time after the current time of the node. (This last condition is to ensure
                ! that a node which was stopped at precisely the time of lightcone crossing is not marked to be crossing the
                ! lightcone again at that same time.)
                nodePositionStart=self%nodePositionReplicant(node             ,timeCrossing,self%origin,[i,j,k],setTime=.true.)
                isInFieldOfView  =self%isInFieldOfView      (nodePositionStart                                                )
                if (reporting) then
                   write (label,'(e12.6)') timeCrossing
                   call displayMessage("Time crossing: "//trim(adjustl(label)))
                   write (label,'(l1)'   ) isInFieldOfView
                   call displayMessage("In FOV:        "//trim(adjustl(label)))
                end if
                if (isInFieldOfView) then
                   if (reporting) then
                      countSeek=0
                      countSteps=int((timeEnd-timeStart_)/timeSeekStep)+1
                      do iSeek=1,countSteps
                         distanceCrossingSeek=timeCrossingRoot(min(timeStart_+dble(iSeek-1)*timeSeekStep,timeEnd))
                      end do
                      call displayIndent("Lightcone crossing seek fine:")
                      call displayMessage("time         distance root")
                      call displayMessage("__________________________")
                     do iSeek=1,countSteps
                         write (label,'(e12.6,1x,e12.6)') timeSeek(iSeek),distanceSeek(iSeek)
                         call displayMessage(trim(adjustl(label)))
                      end do
                      call displayUnindent("")
                      deallocate(timeSeek    )
                      deallocate(distanceSeek)
                   end if
                   ! Only set this crossing as the result if it is the earliest crossing time found so far.
                   if (timeCrossing < squareTimeLightconeCrossing) then
                      squareTimeLightconeCrossing                     =timeCrossing
                      self                       %nodeUniqueIDCrossing=node%uniqueID()
                      self                       %nodePositionCrossing=self%nodePositionReplicant(node,timeCrossing,self%origin,[i,j,k],setTime=.true.)
                      self                       %nodeVelocityCrossing=self%nodeVelocityReplicant(node,timeCrossing            ,[i,j,k],setTime=.true.)
                   end if
                   if (present(timesCrossing)) then
                      ! Append this crossing time to the list of all crossing times.
                      if (allocated(timesCrossing)) then
                         call move_alloc(timesCrossing,timesCrossingTmp)
                         allocate(timesCrossing(size(timesCrossingTmp)+1))
                         timesCrossing(1:size(timesCrossingTmp))=timesCrossingTmp
                         deallocate(timesCrossingTmp)
                      else
                         allocate(timesCrossing(1))
                      end if
                      timesCrossing(size(timesCrossing))=timeCrossing
                   end if
                end if
             end if
             if (reporting) call displayUnindent("")
          end do
       end do
    end do
    if (reporting) call displayUnindent('')
    ! Must reset position and velocity as these can be used in the pre-determined solution.
    call nodeOperatorPredeterminedSolveAnalytics_(self%nodeOperator_,node,timeOriginal)
    ! Sort crossing times if necessary.
    if (present(timesCrossing)) then
       if (allocated(timesCrossing)) then
          if (size(timesCrossing) > 1) call sort(timesCrossing)
       else
          allocate(timesCrossing(0))
       end if
    end if    
    return
  end function squareTimeLightconeCrossing
  
  double precision function timeCrossingRoot(time)
    !!{
    Function used to find the time at which a node crosses the lightcone.
    !!}
    use :: Vectors, only : Vector_Magnitude
    implicit none
    double precision, intent(in   )             :: time
    double precision, dimension(3)              :: positionNode
    double precision, dimension(:), allocatable :: timeSeek_   , distanceSeek_
    double precision                            :: distanceNode
    
    positionNode    =self_%nodePositionReplicant(node_,time,self_%origin,[i_,j_,k_],setTime=.true.)
    distanceNode    =Vector_Magnitude(positionNode)
    timeCrossingRoot=+                          distanceNode           &
         &           -self_%cosmologyFunctions_%distanceComoving(time)
      if (reporting) then
         if (countSeek == 0) then
            allocate(timeSeek    (32))
            allocate(distanceSeek(32))
         else if (countSeek == size(timeSeek)) then
            call move_alloc(timeSeek    ,timeSeek_    )
            call move_alloc(distanceSeek,distanceSeek_)
            allocate(timeSeek    (countSeek*2))
            allocate(distanceSeek(countSeek*2))
            timeSeek    (1:countSeek)=timeSeek_
            distanceSeek(1:countSeek)=distanceSeek_
            deallocate(timeSeek_    )
            deallocate(distanceSeek_)
         end if
         countSeek=countSeek+1
         timeSeek    (countSeek)=time
         distanceSeek(countSeek)=timeCrossingRoot
      end if
    return
  end function timeCrossingRoot
    
  function squareIsInFieldOfView(self,position) result(isInFieldOfView)
    !!{
    Return true if the given position is in the field of view of the lightcone.
    !!}
    implicit none
    logical                                                                :: isInFieldOfView
    class           (geometryLightconeSquare), intent(inout)               :: self
    double precision                         , intent(in   ), dimension(3) :: position

    isInFieldOfView= abs(atan2(position(2),position(1))) < 0.5d0*self%angularsize &
         &          .and.                                                         &
         &           abs(atan2(position(3),position(1))) < 0.5d0*self%angularsize
    return
  end function squareIsInFieldOfView

  function squareReplicantLightConeCrossing(self,node) result(replicant)
    !!{
    Return the indices of the replicant for which the given {\normalfont \ttfamily node} is crossing the lightcone.
    !!}
    use :: Error               , only : Error_Report
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    use :: Vectors             , only : Vector_Magnitude
    use :: Numerical_Comparison, only : Values_Agree
    use :: ISO_Varying_String  , only : var_str
    use :: String_Handling     , only : operator(//)
    implicit none
    integer                                  , dimension(3  ) :: replicant
    class           (geometryLightconeSquare), intent(inout)  :: self
    type            (treeNode               ), intent(inout)  :: node
    class           (nodeComponentBasic     ), pointer        :: basic
    integer                                  , dimension(3,2) :: periodicRange
    double precision                         , dimension(3  ) :: positionNode
    integer                                                   :: i                , j           , &
         &                                                       k
    double precision                                          :: distanceLightcone, distanceNode, &
         &                                                       time
    
    basic             => node                     %basic           (                                                      )
    time              =  basic                    %time            (                                                      )
    distanceLightcone =  self %cosmologyFunctions_%distanceComoving(time                                                  )
    periodicRange     =  self                     %periodicRange   (distanceLightcone,distanceLightcone,radiusBuffer=0.0d0)
    do       i=periodicRange(1,1),periodicRange(1,2)
       do    j=periodicRange(2,1),periodicRange(2,2)
          do k=periodicRange(3,1),periodicRange(3,2)
             positionNode   =self%nodePositionReplicant(node,time,self%origin,[i,j,k],setTime=.true.)
             distanceNode   =Vector_Magnitude(positionNode)
             if (self%isInFieldOfView(positionNode) .and. Values_Agree(distanceNode,distanceLightcone,relTol=1.0d-6)) then
                replicant=[i,j,k]
                return
             end if
          end do
       end do
    end do
    replicant=0
    call Error_Report(var_str('lightcone crossing not found for node ')//node%index()//{introspection:location})    
    return
  end function squareReplicantLightConeCrossing
  
  function squarePositionLightconeCrossing(self,node) result(position)
    !!{
    Return the position at lightcone crossing for this node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    double precision                         , dimension(3)  :: position
    class           (geometryLightconeSquare), intent(inout) :: self
    type            (treeNode               ), intent(inout) :: node
    integer                                  , dimension(3)  :: replicant
    class           (nodeComponentBasic     ), pointer       :: basic

    replicant =  self%replicantLightconeCrossing(node                                                  )
    basic     => node%basic                     (                                                      )
    position  =  self%nodePositionReplicant     (node,basic%time(),self%origin,replicant,setTime=.true.)
    return
  end function squarePositionLightconeCrossing
  
  function squareVelocityLightconeCrossing(self,node) result(velocity)
    !!{
    Return the velocity at lightcone crossing for this node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    double precision                         , dimension(3)  :: velocity
    class           (geometryLightconeSquare), intent(inout) :: self
    type            (treeNode               ), intent(inout) :: node
    integer                                  , dimension(3)  :: replicant
    class           (nodeComponentBasic     ), pointer       :: basic
    
    replicant =  self%replicantLightconeCrossing(node                                      )
    basic     => node%basic                     (                                          )
    velocity  =  self%nodeVelocityReplicant     (node,basic%time(),replicant,setTime=.true.)    
    return
  end function squareVelocityLightconeCrossing
  
  function squarePeriodicRange(self,distanceMinimum,distanceMaximum,radiusBuffer,originBuffered,distanceMinimumBuffered,distanceMaximumBuffered)
    !!{
    Compute the range of possible replicants of the volume which could be in the lightcone in the given interval.
    !!}
    implicit none
    integer                                                 , dimension(3,2)           :: squarePeriodicRange
    class           (geometryLightconeSquare), intent(inout)                           :: self
    double precision                         , intent(in   )                           :: distanceMinimum        , distanceMaximum
    double precision                         , intent(in   )                , optional :: radiusBuffer
    double precision                         , intent(  out)                , optional :: distanceMinimumBuffered, distanceMaximumBuffered
    double precision                         , intent(  out), dimension(3  ), optional :: originBuffered
    double precision                                        , dimension(3  )           :: origin
    double precision                                                                   :: distanceMinimum_       , distanceMaximum_       , &
         &                                                                                lengthOffsetOrigin
    !![
    <optionalArgument name="radiusBuffer" defaultsTo="0.0d0" />
    !!]

    ! If a buffer radius (i.e. a point is to be considered to be inside the lightcone if it is within that distance from an edge
    ! of the cone) is being specified then adjust origin. We shift the origin backward along the lightcone principal axis such
    ! that the cone is effectively broadened by the desired radius. In the diagram below, the inner cone is the true lightcone
    ! geometry, the outer cone is shifted back by a distance L, such that the edge of the outer cone is offset from the inner cone
    ! by the desired buffer radius, R. From the geometry it is clear that L=R/sin(θ/2):
    !
    !  \ \         / /
    !   \ \       /R/
    !    \ \     /⤡/
    !     \ \ ∡θ/ /      L=R/sin(θ/2)
    !      \ \ / /
    !       \ ↑L/
    !        \↓/
    !
    ! The minimum and maximum distances for this section of the cone must then be offset by +L to counteract the shift in
    ! origin. Finally, the minimum/maximum distance is decreased/increased by R such that the caps of the cone section are
    ! buffered by a distance R also.
    if (radiusBuffer_ > 0.0d0) then
       lengthOffsetOrigin=+radiusBuffer_           &
            &             /sin(                    &
            &                  +0.5d0              &
            &                  *self%angularSize   &
            &                 )
       origin            =+self%origin             &
            &             -self%unitVector   (:,1) &
            &             *lengthOffsetOrigin
       distanceMinimum_  =+distanceMinimum+lengthOffsetOrigin-radiusBuffer_
       distanceMaximum_  =+distanceMaximum+lengthOffsetOrigin+radiusBuffer_
    else
       origin             =+self%origin
       distanceMinimum_   =+distanceMinimum
       distanceMaximum_   =+distanceMaximum
    end if
    if (present(originBuffered         )) originBuffered         =origin
    if (present(distanceMinimumBuffered)) distanceMinimumBuffered=distanceMinimum_
    if (present(distanceMaximumBuffered)) distanceMaximumBuffered=distanceMaximum_
    ! Determine range of possible replicants of this galaxy which could be in the lightcone.
    squarePeriodicRange(:,1)=+floor  (                                                                          &
         &                            +(                                                                        &
         &                              +origin                                                                 &
         &                              +cos(0.5d0*self%angularSize)*    distanceMinimum_*self%unitVector(:,1)  &
         &                              -sin(0.5d0*self%angularSize)                                            &
         &                              *(                                                                      &
         &                                +                          abs(distanceMaximum_*self%unitVector(:,2)) &
         &                                +                          abs(distanceMaximum_*self%unitVector(:,3)) &
         &                               )                                                                      &
         &                             )                                                                        &
         &                            /self%lengthReplication                                                   &
         &                          )                                                                           &
         &                   -1
    squarePeriodicRange(:,2)=+ceiling(                                                                          &
         &                            +(                                                                        &
         &                              +origin                                                                 &
         &                              +                                distanceMaximum_*self%unitVector(:,1)  &
         &                              +sin(0.5d0*self%angularSize)                                            &
         &                              *(                                                                      &
         &                                +                          abs(distanceMaximum_*self%unitVector(:,2)) &
         &                                +                          abs(distanceMaximum_*self%unitVector(:,3)) &
         &                               )                                                                      &
         &                             )                                                                        &
         &                            /self%lengthReplication                                                   &
         &                          )                                                                           &
         &                   +0
    return
  end function squarePeriodicRange

  function squareNodePositionReplicant(self,node,time,origin,replicant,setTime,report)
    !!{
    Compute the comoving position of the given node in the given replicant.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic                      , nodeComponentPosition
    use :: Functions_Global, only : nodeOperatorPredeterminedSolveAnalytics_
    use :: Display         , only : displayMessage                          , displayIndent        , displayUnindent
    implicit none
    double precision                                        , dimension(3)           :: squareNodePositionReplicant
    class           (geometryLightconeSquare), intent(inout)                         :: self
    type            (treeNode               ), intent(inout)                         :: node
    double precision                         , intent(in   )                         :: time
    double precision                         , intent(in   ), dimension(3)           :: origin
    integer                                  , intent(in   ), dimension(3)           :: replicant
    logical                                  , intent(in   )              , optional :: setTime                    , report
    class           (nodeComponentBasic     ), pointer                               :: basic
    class           (nodeComponentPosition  ), pointer                               :: position
    double precision                                        , dimension(3)           :: positionComovingNode
    integer                                                                          :: i
    double precision                                                                 :: timeOriginal               , expansionFactor
    character       (len=64                 )                                        :: label
    !![
    <optionalArgument name="setTime" defaultsTo=".false." />
    <optionalArgument name="report"  defaultsTo=".false." />
    !!]

    if (setTime_) then
       basic        => node %basic()
       timeOriginal =  basic%time ()
       call nodeOperatorPredeterminedSolveAnalytics_(self%nodeOperator_,node,time)
    end if
    position             =>  node    %position                           (    )
    expansionFactor      =  +self    %cosmologyFunctions_%expansionFactor(time)
    positionComovingNode =  +position%position                           (    ) &
         &                  /                             expansionFactor
    if (report_) then
       call displayIndent("Replicant position")
       write (label,'(3(e12.6,1x))') time
       call displayMessage("Time:              "//trim(adjustl(label)))
       write (label,'(3(e12.6,1x))') expansionFactor
       call displayMessage("Expansion factor:  "//trim(adjustl(label)))
       write (label,'(3(e12.6,1x))') positionComovingNode
       call displayMessage("Position comoving: "//trim(adjustl(label)))
       call displayUnindent("")
    end if
    do i=1,3
       squareNodePositionReplicant(i)=Dot_Product(positionComovingNode-origin+self%lengthReplication*dble(replicant),self%unitVector(:,i))
    end do
    if (setTime_) call nodeOperatorPredeterminedSolveAnalytics_(self%nodeOperator_,node,timeOriginal)
    return
  end function squareNodePositionReplicant

  function squareNodeVelocityReplicant(self,node,time,replicant,setTime)
    !!{
    Compute the physical velocity of the given node in the given replicant.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentPosition
    implicit none
    double precision                                        , dimension(3)           :: squareNodeVelocityReplicant
    class           (geometryLightconeSquare), intent(inout)                         :: self
    type            (treeNode               ), intent(inout)                         :: node
    double precision                         , intent(in   )                         :: time
    integer                                  , intent(in   ), dimension(3)           :: replicant
    logical                                  , intent(in   )              , optional :: setTime
    class           (nodeComponentBasic     ), pointer                               :: basic
    class           (nodeComponentPosition  ), pointer                               :: position
    double precision                                        , dimension(3)           :: velocity
    integer                                                                          :: i
    double precision                                                                 :: timeOriginal
    !![
    <optionalArgument name="setTime" defaultsTo=".false." />
    !!]

    if (setTime_) then
       basic        => node %basic()
       timeOriginal =  basic%time ()
       call basic%timeSet(time)
    end if
    position => node    %position()
    velocity =  position%velocity()
    do i=1,3
       squareNodeVelocityReplicant(i)=Dot_Product(velocity,self%unitVector(:,i))
    end do     
    if (setTime_) call basic%timeSet(timeOriginal)
    return
  end function squareNodeVelocityReplicant

  subroutine squareDeepCopyReset(self)
    !!{
    Perform a deep copy reset of the object.
    !!}
    use :: Functions_Global, only : nodeOperatorDeepCopyReset_
    implicit none
    class(geometryLightconeSquare), intent(inout) :: self
    
    self%copiedSelf => null()
    if (associated(self%cosmologyparameters_)) call self%cosmologyparameters_%deepCopyReset()
    if (associated(self%cosmologyfunctions_ )) call self%cosmologyfunctions_ %deepCopyReset()
    if (associated(self%outputtimes_        )) call self%outputtimes_        %deepCopyReset()
    if (associated(self%nodeOperator_       )) call nodeOperatorDeepCopyReset_(self%nodeOperator_)
    return
  end subroutine squareDeepCopyReset
  
  subroutine squareDeepCopyFinalize(self)
    !!{
    Finalize a deep reset of the object.
    !!}
    use :: Functions_Global, only : nodeOperatorDeepCopyFinalize_
    implicit none
    class(geometryLightconeSquare), intent(inout) :: self

    if (associated(self%cosmologyparameters_)) call self%cosmologyparameters_%deepCopyFinalize()
    if (associated(self%cosmologyfunctions_ )) call self%cosmologyfunctions_ %deepCopyFinalize()
    if (associated(self%outputtimes_        )) call self%outputtimes_        %deepCopyFinalize()
    if (associated(self%nodeOperator_       )) call nodeOperatorDeepCopyFinalize_(self%nodeOperator_)
    return
  end subroutine squareDeepCopyFinalize
  
  subroutine squareDeepCopy(self,destination)
    !!{
    Perform a deep copy of the object.
    !!}
    use :: Functions_Global, only : nodeOperatorDeepCopy_
    implicit none
    class(geometryLightconeSquare), intent(inout), target :: self
    class(geometryLightconeClass ), intent(inout)         :: destination

    call self%deepCopy_(destination)
    select type (destination)
    type is (geometryLightconeSquare)
       nullify(destination%nodeOperator_)
       if (associated(self%nodeOperator_)) then
          allocate(destination%nodeOperator_,mold=self%nodeOperator_)
          call nodeOperatorDeepCopy_(self%nodeOperator_,destination%nodeOperator_)
#ifdef OBJECTDEBUG
          if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): galacticstructure : [destination] : ')//loc(destination)//' : '//loc(destination%nodeOperator_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
       end if
    class default
       call Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine squareDeepCopy

