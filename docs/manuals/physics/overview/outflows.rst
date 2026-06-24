Outflows
========

Below is a flowchart indicating the ingredients of Galacticus outflows model. (Galacticus is highly modular - many different ingredients can be included and excluded - this is intended just as a typical example.)

.. mermaid::

   flowchart LR
      ISM
      CGM
      IGM
      ISM -->|"<a href='https://galacticus.readthedocs.io/en/latest/physics/stellarFeedbackOutflows.html' style='text-decoration: none'>outflow (ejective)</a>"| CGM
      ISM -->|"<a href='https://galacticus.readthedocs.io/en/latest/physics/stellarFeedbackOutflows.html' style='text-decoration: none'>outflow (expulsive)</a>"| IGM
