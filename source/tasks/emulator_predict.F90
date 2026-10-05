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
  Implements a task which evaluates an emulator at a set of points.
  !!}

  use :: Statistics_Emulators, only : emulatorClass

  !![
  <task name="taskEmulatorPredict" docformat="rst">
   <description>
   A task which evaluates an ``[emulator]`` at a set of points, given as the prior quantiles of its inputs in the dataset
   ``quantiles`` (of shape :math:`M \times d` for :math:`M` points and :math:`d` inputs, in the order of the emulator's
   inputs) of the HDF5 file ``[quantilesFileName]``. The predicted mean and variance of each output at each point are
   written to datasets ``mean`` and ``variance`` (each of shape :math:`M \times B` for :math:`B` outputs), together with the
   abscissae of the outputs (``x``), to the HDF5 file ``[outputFileName]``.
   </description>
  </task>
  !!]
  type, extends(taskClass) :: taskEmulatorPredict
     !!{RST
     Implementation of a task which evaluates an emulator at a set of points.
     !!}
     private
     class(emulatorClass ), pointer :: emulator_         => null()
     type (varying_string)          :: quantilesFileName          , outputFileName
   contains
     final     ::                       emulatorPredictDestructor
     procedure :: perform            => emulatorPredictPerform
     procedure :: requiresOutputFile => emulatorPredictRequiresOutputFile
  end type taskEmulatorPredict

  interface taskEmulatorPredict
     !!{RST
     Constructors for the :galacticus-class:`taskEmulatorPredict` task.
     !!}
     module procedure emulatorPredictConstructorParameters
     module procedure emulatorPredictConstructorInternal
  end interface taskEmulatorPredict

contains

  function emulatorPredictConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`taskEmulatorPredict` task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (taskEmulatorPredict)                :: self
    type (inputParameters    ), intent(inout) :: parameters
    class(emulatorClass      ), pointer       :: emulator_
    type (varying_string     )                :: quantilesFileName, outputFileName

    !![
    <inputParameter docformat="rst">
      <name>quantilesFileName</name>
      <description>
      The name of the HDF5 file containing the points (as prior quantiles) at which to evaluate the emulator.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>outputFileName</name>
      <description>
      The name of the HDF5 file to which the predictions are written.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="emulator" name="emulator_" source="parameters"/>
    !!]
    self=taskEmulatorPredict(quantilesFileName,outputFileName,emulator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="emulator_"/>
    !!]
    return
  end function emulatorPredictConstructorParameters

  function emulatorPredictConstructorInternal(quantilesFileName,outputFileName,emulator_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`taskEmulatorPredict` task class.
    !!}
    implicit none
    type (taskEmulatorPredict)                        :: self
    type (varying_string     ), intent(in   )         :: quantilesFileName, outputFileName
    class(emulatorClass      ), intent(in   ), target :: emulator_
    !![
    <constructorAssign variables="quantilesFileName, outputFileName, *emulator_"/>
    !!]

    return
  end function emulatorPredictConstructorInternal

  subroutine emulatorPredictDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`taskEmulatorPredict` task class.
    !!}
    implicit none
    type(taskEmulatorPredict), intent(inout) :: self

    !![
    <objectDestructor name="self%emulator_"/>
    !!]
    return
  end subroutine emulatorPredictDestructor

  logical function emulatorPredictRequiresOutputFile(self)
    !!{RST
    Specifies that this task does not require the main output file.
    !!}
    implicit none
    class(taskEmulatorPredict), intent(inout) :: self
    !$GLC attributes unused :: self

    emulatorPredictRequiresOutputFile=.false.
    return
  end function emulatorPredictRequiresOutputFile

  subroutine emulatorPredictPerform(self,status)
    !!{RST
    Evaluate the emulator at each point, and write the predictions.
    !!}
    use :: Display           , only : displayIndent, displayUnindent
    use :: Error             , only : Error_Report , errorStatusSuccess
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File
    use :: ISO_Varying_String, only : char
    implicit none
    class           (taskEmulatorPredict), intent(inout), target                 :: self
    integer                              , intent(  out), optional               :: status
    double precision                     , allocatable  , dimension(:,:)         :: quantiles, mean, &
         &                                                                          variance
    integer                                                                      :: i

    call displayIndent('Begin task: emulator predictions')
    !$ call hdf5Access%set()
    block
      type(hdf5File) :: file
      file=hdf5File(char(self%quantilesFileName),readOnly=.true.)
      call file%readDataset('quantiles',quantiles)
    end block
    !$ call hdf5Access%unset()
    if (size(quantiles,dim=1) /= self%emulator_%countInputs()) call Error_Report('the number of quantiles per point does not match the number of inputs to the emulator'//{introspection:location})
    allocate(mean    (self%emulator_%countOutputs(),size(quantiles,dim=2)))
    allocate(variance(self%emulator_%countOutputs(),size(quantiles,dim=2)))
    do i=1,size(quantiles,dim=2)
       call self%emulator_%predict(quantiles(:,i),mean(:,i),variance(:,i))
    end do
    !$ call hdf5Access%set()
    block
      type(hdf5File) :: file
      file=hdf5File(char(self%outputFileName),overWrite=.true.,readOnly=.false.)
      call file%writeDataset(self%emulator_%outputs(),'x'       ,'The abscissae of the outputs.'                  )
      call file%writeDataset(mean                    ,'mean'    ,'The predicted mean of each output at each point.'    )
      call file%writeDataset(variance                ,'variance','The predicted variance of each output at each point.')
    end block
    !$ call hdf5Access%unset()
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: emulator predictions')
    return
  end subroutine emulatorPredictPerform
