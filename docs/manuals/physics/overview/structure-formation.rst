Structure Formation
===================

Below is a flowchart indicating the ingredients of Galacticus structure formation model. (Galacticus is highly modular - many different ingredients can be included and excluded - this is intended just as a typical example.)

.. mermaid::

   flowchart LR
       Cosmology[<a href='https://galacticus.readthedocs.io/en/latest/physics/cosmologyParameters.html' style='text-decoration: none'>Cosmology</a>]
       Power[<a href='https://galacticus.readthedocs.io/en/latest/physics/powerSpectrumPrimordial.html' style='text-decoration: none'>Power spectrum</a>]
       Transfer[<a href='https://galacticus.readthedocs.io/en/latest/physics/transferFunction.html' style='text-decoration: none'>Transfer function</a>]
       Window[<a href='https://galacticus.readthedocs.io/en/latest/physics/powerSpectrumWindowFunction.html' style='text-decoration: none'>Window function</a>]
       Linear[<a href='https://galacticus.readthedocs.io/en/latest/physics/linearGrowth.html' style='text-decoration: none'>Linear growth</a>]
       Variance["<a href='https://galacticus.readthedocs.io/en/latest/physics/cosmologicalMassVariance.html' style='text-decoration: none'>σ(M)</a>"]
       Critical["<a href='https://galacticus.readthedocs.io/en/latest/physics/criticalOverdensity.html' style='text-decoration: none'>δ<sub>c</sub>(t)</a>"]
       Environment["<a href='https://galacticus.readthedocs.io/en/latest/physics/haloEnvironment.html' style='text-decoration: none'>Environment</a>"]
       HMF[<a href='https://galacticus.readthedocs.io/en/latest/physics/haloMassFunction.html' style='text-decoration: none'>Halo mass function</a>]
       Cosmology --> Linear
       Cosmology --> Transfer
       Power --> Transfer
       Environment --> Variance
       Transfer --> Variance
       Window --> Variance
       Linear --> Critical
       Critical --> HMF
       Cosmology --> HMF
       Variance --> HMF
