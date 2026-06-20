Reionization Calculations
=========================

Galacticus can self-consistently solve for the evolution of the as it becomes
photoionized by light emitted by stars and AGN. To activate this
calculation, include the following in your parameters file:

.. code-block:: xml

   <!-- IGM evolver -->
   <intergalacticMediumState value="internal"/>
   <universeOperator value="intergalacticMediumStateEvolve">
     <timeCountPerDecade value=" 10"/>
     <redshiftMaximum    value="150"/>
   </universeOperator>

   <!-- Background radiation -->
   <radiationFieldIntergalacticBackground value="intergalacticBackgroundInternal">
     <wavelengthCountPerDecade value="50"/>
     <timeCountPerDecade       value="10"/>
   </radiationFieldIntergalacticBackground>

   <!-- Halo accretion options -->
   <accretionHalo value="naozBarkana2007"/>

The first block of parameters switches Galacticus to using an internal calculation
for the state of the IGM, instructs it to solve for properties as a
function of time, and specifies that properties should be updated 10
times per decade of cosmic time. Specifically, at each of these time
intervals, solving of galaxy evolution is halted and the IGM evolved up to
this time using the currently computed photoionizing background
spectrum.

The second block of parameters activates an internal calculation of
cosmic background radiation, in which the background is computed from
the emissivities of model galaxies and AGN. The number of points at
which to tabulate the background per decade of wavelength and cosmic
time are specified.

Finally, the third block of parameters tells Galacticus to use the
`Naoz & Barkana (2007 <https://ui.adsabs.harvard.edu/abs/2007MNRAS.377..667N>`_; `naozBarkana2007 <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_) prescription for computing gas accretion into halos
from the IGM. This prescription uses the filtering mass to determine
accretion rates, and will take the filtering mass from the internal
evolution calculation.

With these three sets of configurations, Galacticus will perform a self-consistent
evolution of the IGM—in the sense that the is ionized by photons emitted by
model galaxies and AGN, while galaxy evolution is affected by the
computed state of the model IGM. Note that, when run in this way, Galacticus needs to
keep all merger trees in memory simultaneously (as they are run
synchronously to allow the properties to evolved alongside galaxy
properties).

Once completed, data on the and background radiation are written to the
output file in the ``igmProperties`` and
``backgroundRadiation`` groups respectively.
