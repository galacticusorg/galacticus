Spheroid sizes during bar instabilities
=======================================

Galacticus transfers mass from an unstable disk to a spheroid continuously.
The optional ``barInstabilitySpheroidAngularMomentum`` model specifies the
associated rate of change of the spheroid's pseudo-angular momentum. It does
not change the instability criterion or make the transfer a discrete event.
The default model, ``angularMomentumTransfer``, preserves the historical
prescription.

Motivation
----------

A pressure-supported spheroid need not obtain its support from net rotation.
Assigning its pseudo-angular momentum as a fraction of the physical angular
momentum transferred from a disk can produce very compact spheroids at small
retention fractions. In a self-gravitating limit, a structural relation of
the form :math:`j_\mathrm{s}^2\propto G M_\mathrm{s}R_\mathrm{s}` amplifies this
sensitivity. Here :math:`j_\mathrm{s}` is specific pseudo-angular momentum;
the evolved component property :math:`J_\mathrm{s}=M_\mathrm{s}j_\mathrm{s}`
is its total value.

The ``bindingEnergy`` model instead chooses the bar-driven rate of
:math:`J_\mathrm{s}` from an energy constraint, retaining pseudo-angular
momentum as the variable used by the equilibrium structure solver.

Energy constraint
-----------------

Motivated by the energy-based remnant-size approach of :cite:t:`cole_hierarchical_2000`,
define a positive binding-energy proxy

.. math::

   \mathcal B = c_\mathrm d {G M_\mathrm d^2\over R_\mathrm d}
              +c_\mathrm s {G M_\mathrm s^2\over R_\mathrm s}
              +f_\mathrm{int}{G M_\mathrm d M_\mathrm s\over R_\mathrm d+R_\mathrm s}.

Masses include gas and stars; radii are baryonic half-mass radii, not disk
scale lengths or projected effective radii. The defaults are
:math:`c_\mathrm d=c_\mathrm s=0.5` and :math:`f_\mathrm{int}=1`.
This is an adopted structural proxy, not an exact calculation of the full
galaxy-plus-halo energy. In particular there is no explicit halo interaction
term or radiative dissipation term, although the equilibrium radii respond
to the halo potential.

Let :math:`\dot{\boldsymbol q}_\mathrm{bar}` denote the known rates of disk
and spheroid mass transfer and disk angular-momentum change. The energy
branch chooses

.. math::

   \dot J_\mathrm{s,E} = -{D_{\dot{\boldsymbol q}_\mathrm{bar}}\mathcal B
                               \over \partial\mathcal B/\partial J_\mathrm s},

so that the bar contribution to :math:`\dot{\mathcal B}` vanishes to the
accuracy of the derivative estimates and ODE integration. Other processes
remain free to change :math:`\mathcal B`; it is not conserved over the entire
galaxy history. The derivatives include the equilibrium radii's response.
They are evaluated by finite perturbations of the component state and
repeated structure solves, after which the original state is restored.
``stepFiniteDifferenceRelative`` defaults to 0.01 and controls two independent
forward-difference probes. For the known bar-driven changes, the trial interval
is chosen so that the largest fractional change among disk gas mass, disk
stellar mass and disk angular momentum is 1%. The separate spheroid probe
increases total spheroid pseudo-angular momentum by 1%, with all other
properties fixed. These probes do not advance the galaxy in time. Transferring
1% of the disk can change a small spheroid's mass by much more than 1%, so
convergence should be checked together with the evolution tolerances,
particularly near the transition described below.

The derivatives can be no more accurate than the equilibrium radii, and the
probe step must be small compared with the scales over which the derivatives
change, so both ``stepFiniteDifferenceRelative`` and the structure solver's
``solutionTolerance`` matter. On the isolated tree of
``testSuite/test-bar-instability-binding-energy.py``, the bar changes
:math:`\mathcal B` between the time the spheroid passes
:math:`\mu_\mathrm s=0.01` (defined below) and the end of the run by:

* -42% at the default settings (step 0.01, and the ``equilibrium`` solver's
  default ``solutionTolerance`` of 0.01);
* -11% with ``solutionTolerance`` reduced to :math:`10^{-4}` or
  :math:`10^{-7}`;
* -4% with the step at :math:`10^{-2.5}` and ``solutionTolerance`` at
  :math:`10^{-6}` or :math:`10^{-7}`;
* -41% and -13% with the step at :math:`10^{-3.5}` and ``solutionTolerance``
  at 0.01 and :math:`10^{-4}` respectively;
* -0.8% to -0.9% with the step at :math:`10^{-3.5}` and ``solutionTolerance``
  between :math:`10^{-5}` and :math:`10^{-7}`.

For comparison, the angular-momentum-transfer prescription changes it by -73%.
The solver tolerance must be well below the step (here by a factor of about
30 or more) for the probes to resolve the response. Set
``stepFiniteDifferenceRelative`` to :math:`10^{-3.5}` and ``solutionTolerance``
to :math:`10^{-7}` (as in the example below), and check convergence in both
parameters for the model at hand.

Nascent spheroids
-----------------

When the spheroid is very small, the derivative with respect to its
pseudo-angular momentum becomes poorly conditioned. Define

.. math::

   \mu_\mathrm s={M_\mathrm s\over M_\mathrm d+M_\mathrm s}.

Below :math:`\mu_\mathrm s=0.003`, use the historical angular-momentum-transfer rate.
Above :math:`\mu_\mathrm s=0.010`, use the energy rate. Between these limits,
with :math:`x=(\mu_\mathrm s-0.003)/(0.010-0.003)`, interpolate as

.. math::

   \dot J_\mathrm s=(1-w)\dot J_\mathrm{s,transfer}+w\dot J_\mathrm{s,E},
   \qquad w=3x^2-2x^3.

This regularization is a modeling choice, not a consequence of energy
conservation. The blended rate does not in general conserve the proxy.
Zero-mass or zero-angular-momentum spheroids, and invalid or numerically
degenerate energy derivatives, fall back to the angular-momentum-transfer prescription.
Consequently the spheroid retention factor still controls seeding and
fallbacks even when the established-spheroid energy rate is insensitive to it.

Using the model
---------------

The following fragment belongs inside a complete parameter file. The
``nodeOperator`` belongs within the existing multi-operator list; do not add
a second bar-instability operator. ``fractionAngularMomentumRetainedDisk``
should retain the intended value from the surrounding model independently
of the spheroid factor; the value below is an example, not a calibration.

.. code-block:: xml

   <galacticStructureSolver value="equilibrium">
     <solutionTolerance value="1.0e-7"/>
   </galacticStructureSolver>
   <nodeOperator value="barInstability">
     <galacticDynamicsBarInstability value="efstathiou1982">
       <fractionAngularMomentumRetainedSpheroid value="1.0"/>
       <fractionAngularMomentumRetainedDisk value="0.8"/>
     </galacticDynamicsBarInstability>
     <barInstabilitySpheroidAngularMomentum value="bindingEnergy">
       <stepFiniteDifferenceRelative value="3.1622776601683794e-4"/>
       <fractionBaryonicSpheroidTransitionMinimum value="0.003"/>
       <fractionBaryonicSpheroidTransitionMaximum value="0.010"/>
     </barInstabilitySpheroidAngularMomentum>
   </nodeOperator>

The transition parameters shown explicitly are the defaults. Omitting
``barInstabilitySpheroidAngularMomentum`` altogether selects
``angularMomentumTransfer`` and preserves the previous prescription. Existing
parameter files need no modification to retain that behavior, and the
existing disk and spheroid retention parameters keep their names and meanings.
Changing this prescription changes the physical model and may require
recalibration; agreement with observations is not guaranteed at parameters
calibrated with another size model.

Multiplying all three energy coefficients by the same positive constant leaves
the ideal energy-conserving angular-momentum rate unchanged. Only their ratios
matter, so at least one coefficient should remain fixed during calibration.
For example, fixing :math:`c_\mathrm d` leaves :math:`c_\mathrm s/c_\mathrm d`
and :math:`f_\mathrm{int}/c_\mathrm d` as the two independent parameters.
These change the relative contributions of spheroid self-binding and
disk-spheroid interaction energy; they do not simply multiply the resulting
spheroid radius by a constant. At extreme coefficient normalizations, numerical
fallbacks can affect this scale invariance.
