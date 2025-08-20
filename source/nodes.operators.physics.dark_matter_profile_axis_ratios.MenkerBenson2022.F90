!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implements a node operator class that initializes halo axis ratios using the model of \cite{menker_random_2022}.
  !!}
  
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Virial_Orbits           , only : virialOrbitClass

  !![
  <nodeOperator name="nodeOperatorHaloAxisRatiosMenkerBenson2022">
   <description>
    A node operator class that initializes halo axis ratios using the model of \cite{menker_random_2022}.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorHaloAxisRatiosMenkerBenson2022
     !!{
     A node operator class that initializes halo axis ratios using the model of \cite{menker_random_2022}.
     !!}
     private
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_               => null()
     class           (virialOrbitClass         ), pointer :: virialOrbit_                        => null()
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_                => null()
     double precision                                     :: timescaleSphericalizationFractional          , exponentMass, &
          &                                                  energyBoost
     integer                                              :: ellipsoidEigenvaluesID
   contains
     !![
     <methods>
       <method method="energyTensorEigenvalues" description="Compute the energy tensor eigenvalues of a given node."/>
       <method method="energyTensor"            description="Compute the energy tensor of a given node."            />
       <method method="energyTensorOrbital"     description="Compute the orbital energy tensor of a given node."    />
     </methods>
     !!]
     final     ::                            haloAxisRatiosMenkerBenson2022Destructor
     procedure :: nodeTreeInitialize      => haloAxisRatiosMenkerBenson2022NodeTreeInitialize
     procedure :: energyTensorEigenvalues => haloAxisRatiosMenkerBenson2022EnergyTensorEigenvalues
     procedure :: energyTensor            => haloAxisRatiosMenkerBenson2022EnergyTensor
     procedure :: energyTensorOrbital     => haloAxisRatiosMenkerBenson2022EnergyTensorOrbital
  end type nodeOperatorHaloAxisRatiosMenkerBenson2022
  
  interface nodeOperatorHaloAxisRatiosMenkerBenson2022
     !!{
     Constructors for the \refClass{nodeOperatorHaloAxisRatiosMenkerBenson2022} node operator class.
     !!}
     module procedure haloAxisRatiosMenkerBenson2022ConstructorParameters
     module procedure haloAxisRatiosMenkerBenson2022ConstructorInternal
  end interface nodeOperatorHaloAxisRatiosMenkerBenson2022
  
contains
  
  function haloAxisRatiosMenkerBenson2022ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorHaloAxisRatiosMenkerBenson2022} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorHaloAxisRatiosMenkerBenson2022)                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass                 ), pointer       :: darkMatterProfileDMO_
    class           (virialOrbitClass                          ), pointer       :: virialOrbit_
    class           (darkMatterHaloScaleClass                  ), pointer       :: darkMatterHaloScale_
    double precision                                                            :: timescaleSphericalizationFractional, exponentMass, &
         &                                                                         energyBoost

    !![
    <inputParameter>
      <name>energyBoost</name>
      <defaultValue>0.673d0</defaultValue>
      <defaultSource>\citep{johnson_random_2021}</defaultSource>
      <source>parameters</source>
      <description>A boost to the energy.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentMass</name>
      <defaultValue>1.518d0</defaultValue>
      <source>parameters</source>
      <description>The exponent of mass ratio appearing in the orbital angular momentum term in the axis ratio model.</description>
    </inputParameter>
    <inputParameter>
      <name>timescaleSphericalizationFractional</name>
      <defaultValue>1.75d-3</defaultValue>
      <defaultSource>\citep{menker_random_2022}</defaultSource>
      <source>parameters</source>
      <description>The timescale (in units of the halo dynamical time) for sphericalization.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="virialOrbit"          name="virialOrbit_"          source="parameters"/>
    !!]
    self=nodeOperatorHaloAxisRatiosMenkerBenson2022(timescaleSphericalizationFractional,energyBoost,exponentMass,darkMatterProfileDMO_,darkMatterHaloScale_,virialOrbit_)
    !![
    <objectDestructor name="darkMatterProfileDMO_"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="virialOrbit_"         />
    !!]
    return
  end function haloAxisRatiosMenkerBenson2022ConstructorParameters

  function haloAxisRatiosMenkerBenson2022ConstructorInternal(timescaleSphericalizationFractional,energyBoost,exponentMass,darkMatterProfileDMO_,darkMatterHaloScale_,virialOrbit_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorHaloAxisRatiosMenkerBenson2022} node operator class.
    !!}
    use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMONFW
    use :: Error                   , only : Error_Report
    implicit none
    type            (nodeOperatorHaloAxisRatiosMenkerBenson2022)                        :: self
    class           (darkMatterProfileDMOClass                 ), intent(in   ), target :: darkMatterProfileDMO_
    class           (virialOrbitClass                          ), intent(in   ), target :: virialOrbit_
    class           (darkMatterHaloScaleClass                  ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                            , intent(in   )         :: timescaleSphericalizationFractional, energyBoost, &
         &                                                                                 exponentMass
    !![
    <constructorAssign variables="timescaleSphericalizationFractional, energyBoost, exponentMass, *darkMatterProfileDMO_, *darkMatterHaloScale_, *virialOrbit_"/>
    !!]
    
    !![
    <addMetaProperty component="darkMatterProfile" name="ellipsoidEigenvaluesID" id="self%ellipsoidEigenvaluesID" rank="1" isEvolvable="no" isCreator="yes"/>
    !!]
    ! Validate the dark matter profile.
    select type (darkMatterProfileDMO_)
    type is (darkMatterProfileDMONFW)
       ! This is as expected. Nothing to do.
    class default
       call Error_Report('this model is applicable only to NFW dark matter profiles'//{introspection:location})
    end select
    return
  end function haloAxisRatiosMenkerBenson2022ConstructorInternal

  subroutine haloAxisRatiosMenkerBenson2022Destructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorHaloAxisRatiosMenkerBenson2022} node operator class.
    !!}
    implicit none
    type(nodeOperatorHaloAxisRatiosMenkerBenson2022), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%virialOrbit_"         />
     !!]
    return
  end subroutine haloAxisRatiosMenkerBenson2022Destructor

  subroutine haloAxisRatiosMenkerBenson2022NodeTreeInitialize(self,node)
    !!{
    Initialize the spin of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentDarkMatterProfile, nodeComponentBasic
    use :: Linear_Algebra            , only : matrix                        , vector            , matrixRotationRandom, assignment(=), &
         &                                    operator(*)
    use :: Multidimensional_Minimizer, only : multiDMinimizer
    use :: Sorting                   , only : sort
    use :: Array_Utilities           , only : Array_Reverse
    implicit none
    class           (nodeOperatorHaloAxisRatiosMenkerBenson2022), intent(inout), target      :: self
    type            (treeNode                                  ), intent(inout), target      :: node
    type            (treeNode                                  )               , pointer     :: nodeSibling        , nodeChild
    type            (multiDMinimizer                           )               , allocatable :: minimizer_
    class           (nodeComponentDarkMatterProfile            )               , pointer     :: darkMatterProfile  , darkMatterProfileChild
    class           (nodeComponentBasic                        )               , pointer     :: basic              , basicChild             , &
         &                                                                                      basicSibling
    double precision                                            , dimension(3)               :: axisRatios         , eigenvalues_
    type            (matrix                                    )                             :: energyTensorPrimary, energyTensorSecondary  , &
         &                                                                                      energyTensorOrbital, matrixRotationSecondary, &
         &                                                                                      eigenvectors       , matrixIntermediate0    , &
         &                                                                                      matrixIntermediate1, matrixIntermediate2    , &
         &                                                                                      matrixIntermediate3, matrixIntermediate4    , &
         &                                                                                      matrixIntermediate5
    type            (vector                                    )                             :: eigenvalues
    double precision                                                                         :: frequencyPrecession
    integer                                                                                  :: iteration
    logical                                                                                  :: converged
    
    ! Create a dark matter profile in the node.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    ! If this node has no children, assume it to be spherical.
    if (.not.associated(node%firstChild)) then
       axisRatios  =[1.0d0,1.0d0,1.0d0]
       call darkMatterProfile%axisRatiosSet            (                            axisRatios  )
       eigenvalues_=self%energyTensorEigenvalues(node)
       call darkMatterProfile%floatRank1MetaPropertySet(self%ellipsoidEigenvaluesID,eigenvalues_)
    else
       ! Node has children, so accumulate their energy tensors, and apply sphericalization.
       !! Find energy tensor of the primary progenitor halo.
       !! NOTE: No correction for unresolved accretion is made here. A more detailed model of this process should be added.
       nodeChild           => node %firstChild
       basicChild          => nodeChild%basic       (         )
       energyTensorPrimary =  self     %energyTensor(nodeChild)
       ! Iterate over any secondary progenitor halos.
       nodeSibling         => nodeChild
       do while(associated(nodeSibling%sibling))
          ! Move to the next sibling (secondary progenitor).
          nodeSibling  => nodeSibling%sibling
          basicSibling => nodeSibling%basic  () 
          ! Compute the energy tensor of the secondary progenitor.
          energyTensorSecondary=self%energyTensor(nodeSibling)
          ! Generate a random rotation matrix for the secondary.
          matrixRotationSecondary=matrixRotationRandom(nodeSibling%hostTree%randomNumberGenerator_)
          ! Compute the orbital energy tensor of the secondary progenitor.
          energyTensorOrbital=self%energyTensorOrbital(nodeSibling)
          ! Sum energy tensors (Menker & Benson 2022, equation 19).                
          matrixIntermediate0= matrixRotationSecondary%transpose()
          matrixIntermediate1= energyTensorSecondary  *matrixIntermediate0
          matrixIntermediate2= matrixRotationSecondary*matrixIntermediate1
          matrixIntermediate3= matrixIntermediate2    +energyTensorOrbital
          matrixIntermediate4= matrixIntermediate3                                                &
               &              *(                                                                  &
               &                +(1.0d0+self%energyBoost                     )                    &
               &                /(1.0d0+basicSibling%mass()/basicChild%mass())**self%exponentMass &
               &               )                                                                  &
               &              *nodeSibling%subsamplingWeight()
          matrixIntermediate5= energyTensorPrimary    +matrixIntermediate4
          energyTensorPrimary= matrix(matrixIntermediate5)
       end do
       ! Perform an eigen-decomposition of the energy tensor.
       call energyTensorPrimary%eigensystem(eigenvectors,eigenvalues)
       eigenvalues_=eigenvalues
       ! Apply sphericalization.
       basic                  =>  node                  %basic            ()
       darkMatterProfileChild =>  nodeChild             %darkMatterProfile()
       axisRatios             =   darkMatterProfileChild%axisRatios       ()
       !! Menker & Benson (2022; equation 20).
       frequencyPrecession    =  +(                                                       &
            &                      +1.0d0                                                 &
            &                      -(                                                     &
            &                        +axisRatios(2)                                       &
            &                        *axisRatios(3)                                       &
            &                       )**0.25d0                                             &
            &                     )                                                       &
            &                    /self%darkMatterHaloScale_%timescaleDynamical(nodeChild)
       !! Menker & Benson (2022; equation 23).
       eigenvalues_=+              sum(eigenvalues_)/3.0d0                                &
            &       +(eigenvalues_-sum(eigenvalues_)/3.0d0)                               &
            &       *exp(                                                                 &
            &            -self%timescaleSphericalizationFractional*frequencyPrecession    &
            &            *(                                                               &
            &              +basic     %time()                                             &
            &              -basicChild%time()                                             &
            &             )                                                               &
            &           )
       call sort(eigenvalues_)
       ! Test boundedness of the halo.
       if (all(eigenvalues_ < 0.0d0)) then
          ! Energy tensor eigenvalues are all negative (i.e. halo is "bound" along all three principle axes). Compute the
          ! axis ratios. To do this we use a multi-dimensional minimizer to find the axis ratios that correspond to the
          ! current eigenvalue ratios.
          !! Build a minimizer object that we will use to find axis ratios given the eigenvalues of the energy tensor.
          allocate(minimizer_)
          minimizer_=multiDMinimizer(2_c_size_t,axisRatioCost)
          call minimizer_%set(x=[0.5d0,0.5d0],stepSize=[0.01d0,0.01d0])
          iteration=0
          converged=.false.
          do while (.not.converged .and. iteration < 100)
             call minimizer_%iterate()
             iteration=iteration+1
             converged=minimizer_%testSize(toleranceAbsolute=1.0d-12)
          end do
          axisRatios(2:3)=minimizer_%x()
          axisRatios(1  )=    1.0d0
          axisRatios(2  )=min(1.0d0,axisRatios(2))
          axisRatios(3  )=min(1.0d0,axisRatios(3))
          call sort(axisRatios)
          axisRatios=Array_Reverse(axisRatios)
          ! If any axis ratios are negative, simply do not update the axis ratios.
          if (any(axisRatios <= 0.0d0)) axisRatios=darkMatterProfileChild%axisRatios()
       else
          ! Some eigenvalues are non-negative - indicating the halo is "unbound" along one or more axes. Do not update the
          ! axis ratios.
          axisRatios=darkMatterProfileChild%axisRatios()
       end if
       call darkMatterProfile%axisRatiosSet            (                            axisRatios  )
       call darkMatterProfile%floatRank1MetaPropertySet(self%ellipsoidEigenvaluesID,eigenvalues_)
    end if
    return
    
  contains
    
    double precision function axisRatioCost(axisRatios_)
      !!{
      Evaluate the cost function used in solving for the ellipsoidal axis ratios given the eigenvalues of the energy tensor.
      !!}
      implicit none
      double precision, intent(in   ), dimension(:) :: axisRatios_
      double precision               , dimension(3) :: coefficientsA_
      
      coefficientsA_=coefficientsA([1.0d0,max(1.0d-6,min(1.0d0,axisRatios_(1))),max(1.0d-6,min(1.0d0,axisRatios_(1),axisRatios_(2)))])
      axisRatioCost =+(axisRatios_(1)**2*coefficientsA_(2)/coefficientsA_(1)-eigenvalues_(2)/eigenvalues_(1))**2 &
           &         +(axisRatios_(2)**2*coefficientsA_(3)/coefficientsA_(1)-eigenvalues_(3)/eigenvalues_(1))**2
      return
    end function axisRatioCost

  end subroutine haloAxisRatiosMenkerBenson2022NodeTreeInitialize

  function coefficientsA(axisRatios)
    !!{
    Evaluate the coefficients, $A_i$ appearing in the expression for the potential energy tensor. We follow the convention of
    \cite[][their equation~12]{menker_random_2022}.
    !!}
    use :: Elliptic_Integrals, only : Incomplete_Elliptic_Integral_E , Incomplete_Elliptic_Integral_F
    implicit none
    double precision               , dimension(3) :: coefficientsA
    double precision, intent(in   ), dimension(3) :: axisRatios
    double precision, parameter                   :: tolerance    =1.0d-2
    double precision                              :: theta                , phi         , &
         &                                           sinTheta             , eccentricity

    ! Handle special cases where some (or all) axes are equal.
    if (axisRatios(2) > (1.0d0-tolerance)*axisRatios(1) .and. axisRatios(3) > (1.0d0-tolerance)*axisRatios(1)) then
       ! Spherical. Since ∑ᵢ₌₁³ Aᵢ = 2 (S. Chandrasekhar, "Ellipsoidal Figures of Equilibrium", Yale University Press, 1969,
       ! Chapter 2, equation 24), then by symmetry Aᵢ=⅔ in the spherical case.
       coefficientsA=2.0d0/3.0d0
    else if (axisRatios(2) > (1.0d0-tolerance)*axisRatios(1)) then
       ! Oblate spheroid (S. Chandrasekhar, "Ellipsoidal Figures of Equilibrium", Yale University Press, 1969, Chapter 2, equation
       ! 36-37), but without the a₂a₃/a₁² pre-factor.
       eccentricity =sqrt(1.0d0-(axisRatios(3)/axisRatios(1))**2)
       coefficientsA=[                                                                                                               &
            &               +sqrt(1.0d0-eccentricity**2)*asin(eccentricity)/eccentricity**3-(1.0d0-eccentricity**2)/eccentricity**2, &
            &               +sqrt(1.0d0-eccentricity**2)*asin(eccentricity)/eccentricity**3-(1.0d0-eccentricity**2)/eccentricity**2, &
            &         -2.0d0*sqrt(1.0d0-eccentricity**2)*asin(eccentricity)/eccentricity**3+ 2.0d0                 /eccentricity**2  &
            &        ]/axisRatios(2)/axisRatios(3)
    else if (axisRatios(3) > (1.0d0-tolerance)*axisRatios(2)) then
       ! Prolate spheroid (S. Chandrasekhar, "Ellipsoidal Figures of Equilibrium", Yale University Press, 1969, Chapter 2,
       ! equation 38-39), but without the a₂a₃/a₁² pre-factor.
       eccentricity =sqrt(1.0d0-(axisRatios(3)/axisRatios(1))**2)
       coefficientsA=[                                                                                                                                             &
            &         +(1.0d0-eccentricity**2)*log((1.0d0+eccentricity)/(1.0d0-eccentricity))      /eccentricity**3-2.0d0*(1.0d0-eccentricity**2)/eccentricity**2, &
            &         -(1.0d0-eccentricity**2)*log((1.0d0+eccentricity)/(1.0d0-eccentricity))/2.0d0/eccentricity**3+1.0d0                        /eccentricity**2, &
            &         -(1.0d0-eccentricity**2)*log((1.0d0+eccentricity)/(1.0d0-eccentricity))/2.0d0/eccentricity**3+1.0d0                        /eccentricity**2  &
            &        ]/axisRatios(2)/axisRatios(3)
    else
       ! General ellipsoid (S. Chandrasekhar, "Ellipsoidal Figures of Equilibrium", Yale University Press, 1969, Chapter 2,
       ! equation 31-35), but without the a₂a₃/a₁² pre-factor.
       sinTheta     =sqrt((1.0d0-(axisRatios(2)/axisRatios(1))**2)/(1.0d0-(axisRatios(3)/axisRatios(1))**2))
       theta        =asin(sinTheta)
       phi          =acos(axisRatios(3)/axisRatios(1))
       coefficientsA=[                                                                                                                                                                                                                    &
            &         +2.0d0*(+Incomplete_Elliptic_Integral_F(phi,sinTheta**2)              -Incomplete_Elliptic_Integral_E(phi,sinTheta**2)                                                   )/sin(phi)**3/sin(theta)**2              , &
            &         +2.0d0*(-Incomplete_Elliptic_Integral_F(phi,sinTheta**2)*cos(theta)**2+Incomplete_Elliptic_Integral_E(phi,sinTheta**2)-axisRatios(3)/axisRatios(2)*sin(theta)**2*sin(phi))/sin(phi)**3/sin(theta)**2/cos(theta)**2, &
            &         +2.0d0*(                                                              -Incomplete_Elliptic_Integral_E(phi,sinTheta**2)+axisRatios(2)/axisRatios(3)              *sin(phi))/sin(phi)**3/cos(theta)**2                &
            &        ]
    end if
    return
  end function coefficientsA
  
  function haloAxisRatiosMenkerBenson2022EnergyTensorEigenvalues(self,node) result(eigenvalues)
    !!{
    Compute the energy tensor eigenvalues of the given node.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Galacticus_Nodes                , only : nodeComponentBasic            , nodeComponentDarkMatterProfile
    implicit none
    double precision                                            , dimension(3)  :: eigenvalues
    class           (nodeOperatorHaloAxisRatiosMenkerBenson2022), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    class           (nodeComponentBasic                        ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile            ), pointer       :: darkMatterProfile
    double precision                                            , dimension(3)  :: axisRatios       , coefficientA
    integer                                                                     :: i
    double precision                                                            :: radiusScale      , radiusVirial, &
         &                                                                         densityScale     , massScale   , &
         &                                                                         concentration    , gamma

    ! Extract axis ratios.
    basic             => node             %basic            ()
    darkMatterProfile => node             %darkMatterProfile()
    axisRatios        =  darkMatterProfile%axisRatios       ()
    ! Compute the potential energy tensor.
    !! NOTE: This is specific to an NFW profile as in Menker & Benson (2022). Ideally this should be generalized to an arbitrary
    !! profile.
    coefficientA=coefficientsA(axisRatios)
    do i=1,3
       eigenvalues(i)=+coefficientA(i)    &
            &        *axisRatios  (i)**2
    end do
    radiusScale  =+     darkMatterProfile   %scale       (    )
    radiusVirial =+self%darkMatterHaloScale_%radiusVirial(node)
    concentration=+radiusVirial                                              &
         &        /radiusScale
    gamma        =+9.0d0                                                     &
         &        /4.0d0                                                     &
         &        *(                                                         &
         &          +1.0d0                                                   &
         &          -1.0d0                         /(1.0d0+concentration**2) &
         &          -2.0d0*log(1.0d0+concentration)/(1.0d0+concentration   ) &
         &         )
    densityScale =+3.0d0                                                     &
         &        *basic%mass()                                              &
         &        /4.0d0                                                     &
         &        /Pi                                                        &
         &        /(                                                         &
         &          +              log(1.0d0+concentration)                  &
         &          -concentration/   (1.0d0+concentration)                  &
         &         )                                                         &
         &        /radiusScale**3
    massScale    =+4.0d0                                                     &
         &        *Pi                                                        &
         &        /3.0d0                                                     &
         &        *product(axisRatios)                                       &
         &        *densityScale                                              &
         &        *radiusScale**3
    eigenvalues  =-gamma                                                     &
         &        *gravitationalConstant_internal                            &
         &        *massScale**2                                              &
         &        /radiusScale                                               &
         &        *eigenvalues                                               &
         &        /2.0d0
    return
  end function haloAxisRatiosMenkerBenson2022EnergyTensorEigenvalues

  function haloAxisRatiosMenkerBenson2022EnergyTensor(self,node) result(energyTensor)
    !!{
    Compute the energy tensor of the given node.
    !!}
    use :: Linear_Algebra  , only : matrix
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    type            (matrix                                    )                 :: energyTensor
    class           (nodeOperatorHaloAxisRatiosMenkerBenson2022), intent(inout)  :: self
    type            (treeNode                                  ), intent(inout)  :: node
    class           (nodeComponentDarkMatterProfile            ), pointer        :: darkMatterProfile
    double precision                                            , dimension(3  ) :: eigenvalues
    double precision                                            , dimension(3,3) :: energyTensor_
    integer                                                                      :: i

    ! The energy tensor is constructed directly from the (stored) eigenvalues.
    darkMatterProfile => node             %darkMatterProfile        (                           )
    eigenvalues       =  darkMatterProfile%floatRank1MetaPropertyGet(self%ellipsoidEigenvaluesID)
    energyTensor_     =  0.0d0
    do i=1,3
       energyTensor_(i,i)=eigenvalues(i)
    end do
    energyTensor=matrix(energyTensor_)
    return
  end function haloAxisRatiosMenkerBenson2022EnergyTensor

  function haloAxisRatiosMenkerBenson2022EnergyTensorOrbital(self,node) result(energyTensorOrbital)
    !!{
    Compute the orbital energy tensor of the given node.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Linear_Algebra                  , only : matrix                        , assignment(=)
    use :: Galacticus_Nodes                , only : nodeComponentBasic            , nodeComponentSatellite
    use :: Kepler_Orbits                   , only : keplerOrbit
    use :: Vectors                         , only : Vector_Product                , Vector_Magnitude
    implicit none
    type            (matrix                                    )                 :: energyTensorOrbital
    class           (nodeOperatorHaloAxisRatiosMenkerBenson2022), intent(inout)  :: self
    type            (treeNode                                  ), intent(inout)  :: node
    type            (treeNode                                  ), pointer        :: nodePrimary
    class           (nodeComponentBasic                        ), pointer        :: basicPrimary        , basicSecondary
    class           (nodeComponentSatellite                    ), pointer        :: satellite
    double precision                                            , dimension(3,3) :: energyTensorOrbital_, matrixRotation_             , &
         &                                                                          basisUnrotated      , basisRotated
    integer                                                                      :: i                   , j
    type            (keplerOrbit                               )                 :: orbit
    double precision                                                             :: theta               , phi                         , &
         &                                                                          psi                 , massReduced
    type            (matrix                                    )                 :: matrixRotation      , energyTensorOrbitalUnrotated, &
         &                                                                          matrixIntermediate1 , matrixIntermediate2

    ! Get the orbit.
    nodePrimary    => node       %parent%firstChild
    basicPrimary   => nodePrimary       %basic      (                 )
    basicSecondary => node              %basic      (                 )
    satellite      => node              %satellite  (autoCreate=.true.)
    orbit          =  satellite         %virialOrbit(                 )
    ! Construct the coordinate system for the merger, and the rotation matrix to rotate it into the coordinate system of the
    ! primary halo.
    theta              =+acos(+2.0d0   *node%hostTree%randomNumberGenerator_%uniformSample()-1.0d0)
    phi                =      +2.0d0*Pi*node%hostTree%randomNumberGenerator_%uniformSample()
    psi                =      +2.0d0*Pi*node%hostTree%randomNumberGenerator_%uniformSample()
    basisUnrotated     =reshape(                    &
         &                      [                   &
         &                       1.0d0,0.0d0,0.0d0, &
         &                       0.0d0,1.0d0,0.0d0, &
         &                       0.0d0,0.0d0,1.0d0  &
         &                      ]                 , &
         &                      [3,3]               &
         &                     )
    basisRotated(:,1)=[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(phi)]
    basisRotated(:,2)=Vector_Product(basisRotated(:,1),[1.0d0,0.0d0,0.0d0])
    basisRotated(:,3)=Vector_Product(basisRotated(:,1),basisRotated(:,2)  )
    do i=1,3
       basisRotated(:,i)=+                 basisRotated(:,i)  &
            &            /Vector_Magnitude(basisRotated(:,i))
    end do
    do i=1,3
       do j=1,3
          matrixRotation_(i,j)=Dot_Product(basisUnrotated(:,i),basisRotated(:,j))
       end do
    end do
    matrixRotation=matrix(matrixRotation_)
    ! Construct the orbital energy tensor in the merger coordinate system.
    massReduced                 =  +basicSecondary%mass()*basicPrimary%mass()  &
         &                       /(+basicSecondary%mass()+basicPrimary%mass())
    energyTensorOrbital_        =+0.0d0
    energyTensorOrbital_(1,1)   =+energyTensorOrbital_(1,1)                                     &
         &                       -gravitationalConstant_internal                                &
         &                       *basicPrimary                       %mass        (           ) &
         &                       *basicSecondary                     %mass        (           ) &
         &                       /self          %darkMatterHaloScale_%radiusVirial(nodePrimary)
    energyTensorOrbital_(1,1)   =+energyTensorOrbital_(1,1)                                &
         &                       +0.5d0*massReduced*orbit%velocityRadial    ()**2
    energyTensorOrbital_(2,2)   =+energyTensorOrbital_(2,2)                                &
         &                       +0.5d0*massReduced*orbit%velocityTangential()**2*cos(psi)
    energyTensorOrbital_(3,3   )=+energyTensorOrbital_(3,3)                                &
         &                       +0.5d0*massReduced*orbit%velocityTangential()**2*sin(psi)
    energyTensorOrbitalUnrotated= matrix(energyTensorOrbital_)
    ! Rotate to the coordinate system of the primary halo.
    matrixIntermediate1=matrixRotation%transpose()
    matrixIntermediate2=energyTensorOrbitalUnrotated*matrixIntermediate1
    energyTensorOrbital=matrixRotation              *matrixIntermediate2
    return
  end function haloAxisRatiosMenkerBenson2022EnergyTensorOrbital
