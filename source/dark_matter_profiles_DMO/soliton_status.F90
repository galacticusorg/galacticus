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

  !!{RST
  Contains a module which provides an enumeration describing the structural state of a :term:`FDM` halo in the soliton+NFW dark
  matter profile classes.
  !!}

  module Dark_Matter_Profiles_Soliton_Status
    !!{RST
    Provides an enumeration describing the structural state of a :term:`FDM` halo in the soliton+NFW dark matter profile classes.
    !!}
    implicit none
    private

    !![
    <enumeration docformat="rst">
     <name>solitonStatus</name>
     <description>
     The structural state of a :term:`FDM` halo in the :galacticus-class:`darkMatterProfileDMOSolitonNFW` and
     :galacticus-class:`darkMatterProfileDMOSolitonNFWHeated` dark matter profile classes. The state is cached in the
     ``solitonStatus`` meta-property of the ``darkMatterProfile`` component so that a halo is treated consistently throughout its
     evolution.

     A halo begins ``uninitialized``. On the first attempt to construct its mass distribution it becomes either ``solitonNFW`` (a
     solitonic core embedded in an NFW envelope) if a solution for the soliton transition radius is found, or ``nfwOnly`` if no
     such solution exists. In the heated class a ``solitonNFW`` halo may subsequently become ``solitonOnly``, once tidal stripping
     has removed the NFW envelope and only the solitonic core survives. The ``solitonOnly`` and ``nfwOnly`` states are terminal.

     Note that ``uninitialized`` must remain the first entry of this enumeration, so that it takes the value zero: meta-properties
     are zero-valued until first set, and the code relies on that to detect a halo whose state has not yet been determined.
     </description>
     <entry label="uninitialized"/>
     <entry label="solitonNFW"   />
     <entry label="solitonOnly"  />
     <entry label="nfwOnly"      />
    </enumeration>
    !!]

  end module Dark_Matter_Profiles_Soliton_Status
