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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which implements a class of emulators of model predictions.
!!}

module Statistics_Emulators
  !!{RST
  Implements a class of emulators of model predictions.
  !!}
  use :: ISO_Varying_String, only : varying_string
  private

  !![
  <functionClass docformat="rst">
   <name>emulator</name>
   <descriptiveName>Emulators</descriptiveName>
   <description>
   Class providing emulators of the predictions of a model---fast approximations, trained on a set of runs of the model, to
   the prediction of some observable (for example, the binned values of a stellar mass function) as a function of the
   model's parameters (see :ref:`manual-sec-Emulation`). An emulator's inputs are the prior quantiles of the parameters on
   which it depends, in the order given by ``inputNames``. For each bin of its observable it predicts a mean and a
   variance, the latter quantifying the uncertainty of the emulation itself.
   </description>
   <default>gaussianProcess</default>
   <method name="countInputs" >
    <description>
    Return the number of inputs to the emulator.
    </description>
    <type>integer</type>
    <pass>yes</pass>
   </method>
   <method name="inputNames" >
    <description>
    Return the names of the parameters on which the emulator depends, in the order of its inputs.
    </description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(varying_string), intent(  out), allocatable, dimension(:) :: names</argument>
   </method>
   <method name="countOutputs" >
    <description>
    Return the number of outputs (bins) of the emulated observable.
    </description>
    <type>integer</type>
    <pass>yes</pass>
   </method>
   <method name="outputs" >
    <description>
    Return the abscissae (e.g. the bin centers) of the outputs of the emulated observable.
    </description>
    <type>double precision, allocatable, dimension(:)</type>
    <pass>yes</pass>
   </method>
   <method name="target" >
    <description>
    Return the data to which the emulated observable is compared, and its covariance, in the (transformed) space in which the
    observable is emulated.
    </description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(  out), allocatable, dimension(:  ) :: yTarget</argument>
    <argument>double precision, intent(  out), allocatable, dimension(:,:) :: covarianceTarget</argument>
   </method>
   <method name="trainingPoints" >
    <description>
    Return the prior ``quantiles`` of the emulator's inputs (in the order of ``inputNames``) at each point on which the emulator
    was trained, and the corresponding parameter ``values``. Both arrays have shape ``(countInputs, countPoints)``. These allow
    a user of the emulator to check that its priors are those under which the emulator was trained.
    </description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(  out), allocatable, dimension(:,:) :: quantiles, values</argument>
   </method>
   <method name="predict" >
    <description>
    Predict the mean and variance of the emulated observable, in each of its bins, given the prior ``quantiles`` of the
    emulator's inputs (in the order of ``inputNames``). The arrays ``mean`` and ``variance`` must have size ``countOutputs``.
    </description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), dimension(:) :: quantiles</argument>
    <argument>double precision, intent(  out), dimension(:) :: mean     , variance</argument>
   </method>
  </functionClass>
  !!]

end module Statistics_Emulators
