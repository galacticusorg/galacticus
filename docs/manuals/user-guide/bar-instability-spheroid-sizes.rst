Spheroid sizes during bar instabilities
======================================

Galacticus transfers mass from an unstable disk to a spheroid continuously.
The optional ``barInstabilitySpheroidAngularMomentum`` model specifies the
associated rate of change of the spheroid's pseudo-angular momentum. It does
not change the instability criterion or make the transfer a discrete event.
The default model, ``retained``, preserves the historical prescription.

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

Motivated by the energy-based remnant-size approach of Cole et al. (2000),
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
``finiteDifferenceStepRelative`` defaults to 0.01. Its convergence should be
tested together with the evolution tolerances for the application of interest.

Nascent spheroids
-----------------

When the spheroid is very small, the derivative with respect to its
pseudo-angular momentum becomes poorly conditioned. Define

.. math::

   \mu_\mathrm s={M_\mathrm s\over M_\mathrm d+M_\mathrm s}.

Below :math:`\mu_\mathrm s=0.003`, use the historical retained rate.
Above :math:`\mu_\mathrm s=0.010`, use the energy rate. Between these limits,
with :math:`x=(\mu_\mathrm s-0.003)/(0.010-0.003)`, interpolate as

.. math::

   \dot J_\mathrm s=(1-w)\dot J_\mathrm{s,retained}+w\dot J_\mathrm{s,E},
   \qquad w=3x^2-2x^3.

This regularization is a modelling choice, not a consequence of energy
conservation. The blended rate does not in general conserve the proxy.
Zero-mass or zero-angular-momentum spheroids, and invalid or numerically
degenerate energy derivatives, fall back to the retained prescription.
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

   <galacticStructureSolver value="equilibrium"/>
   <nodeOperator value="barInstability">
     <galacticDynamicsBarInstability value="efstathiou1982">
       <fractionAngularMomentumRetainedSpheroid value="1.0"/>
       <fractionAngularMomentumRetainedDisk value="0.8"/>
     </galacticDynamicsBarInstability>
     <barInstabilitySpheroidAngularMomentum value="bindingEnergy">
       <spheroidBaryonFractionTransitionMinimum value="0.003"/>
       <spheroidBaryonFractionTransitionMaximum value="0.010"/>
     </barInstabilitySpheroidAngularMomentum>
   </nodeOperator>

The transition parameters shown explicitly are the defaults. Omitting
``barInstabilitySpheroidAngularMomentum`` altogether still selects ``retained``.
Changing this prescription changes the physical model and may require
recalibration; agreement with observations is not guaranteed at parameters
calibrated with another size model.

Only two ratios of the three energy coefficients affect the ideal energy
rate: a common positive multiplier cancels between numerator and denominator.
For calibration, one coefficient should therefore remain fixed. Numerical
safeguards can break this invariance at extreme normalizations. The direction
and magnitude of the size response to either ratio should be measured rather
than interpreted as a direct multiplicative radius adjustment.
