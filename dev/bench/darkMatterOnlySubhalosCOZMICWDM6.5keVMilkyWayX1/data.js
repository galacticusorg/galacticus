window.BENCHMARK_DATA = {
  "lastUpdate": 1761478809337,
  "repoUrl": "https://github.com/galacticusorg/galacticus",
  "entries": {
    "Dark matter-only subhalos benchmarks (COZMIC Milky Way WDM 6.5keV resolutionX1)": [
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ce82e1607cbd130bdfe09df824f2ebce1b1a8170",
          "message": "Merge pull request #951 from galacticusorg/featCOZMICWDMValidation\n\nAdd validation runs for all COZMIC WDM models",
          "timestamp": "2025-10-20T14:52:17Z",
          "tree_id": "f5f08b5585e485c2fbc67efa92041b9b512314b3",
          "url": "https://github.com/galacticusorg/galacticus/commit/ce82e1607cbd130bdfe09df824f2ebce1b1a8170"
        },
        "date": 1760996788190,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4675536705526684,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.4672380251014334,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 49.267571961179186,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "Andrew Benson",
            "username": "abensonca",
            "email": "abenson@carnegiescience.edu"
          },
          "committer": {
            "name": "Andrew Benson",
            "username": "abensonca",
            "email": "abenson@carnegiescience.edu"
          },
          "id": "cb9801c036f9c2442eca6d08df4d91f0bd733185",
          "message": "feat: Add an option to the emission line table creation script to allow stopping at the outer radius of the cloud\n\nIf the `--stopOuterRadius` option is set (allowed only if the `'--normalization massStellar` option is also used), then the outer radius of the cloud is computed (based on its mass and density), and Cloudy is instructed to stop at this radius if reached. This accounts for the finite extent of clouds and, in principle, will allow for estimation of the emergent spectrum.",
          "timestamp": "2025-10-20T23:27:10Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/cb9801c036f9c2442eca6d08df4d91f0bd733185"
        },
        "date": 1761037648190,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4618752753440183,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.403510468162678,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 43.9578074887433,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a8684d46308dcd337834a248a072b5f4653b7d70",
          "message": "Merge pull request #953 from galacticusorg/fixConstrainedTreeBuild\n\nAvoid overwriting earlier stages during constrained tree build",
          "timestamp": "2025-10-21T15:26:52Z",
          "tree_id": "382e32c7c1b05994d804edcba303d5ab4dcec49f",
          "url": "https://github.com/galacticusorg/galacticus/commit/a8684d46308dcd337834a248a072b5f4653b7d70"
        },
        "date": 1761097416159,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4306827222561054,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.680998180677329,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 34.160006084438166,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "520abf8874e55ae5faca38acfd76be4e737977e2",
          "message": "Merge pull request #954 from galacticusorg/fixUnrecognizedParameters\n\nFor named parameters in `objectBuilder` directives, add them to the list of recognized parameters for each executable",
          "timestamp": "2025-10-22T14:29:40Z",
          "tree_id": "359aae7be7afe41a9251ff58385863bd86751613",
          "url": "https://github.com/galacticusorg/galacticus/commit/520abf8874e55ae5faca38acfd76be4e737977e2"
        },
        "date": 1761168053299,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4248591671761626,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 3.293865970839539,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 44.260796411930286,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "19548c248a266d2592dfb8068af563e88e717f6f",
          "message": "Merge pull request #956 from galacticusorg/fixIgnoreWellOrdering\n\nUnlink node from tree before destroying its branch",
          "timestamp": "2025-10-23T17:09:22Z",
          "tree_id": "803eaef9fa56b3abf90cfe57d007dcd75038ffe4",
          "url": "https://github.com/galacticusorg/galacticus/commit/19548c248a266d2592dfb8068af563e88e717f6f"
        },
        "date": 1761277078834,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4691897726538756,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.4926843440653328,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 44.71366290588569,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ee1a915ee5a43a63d05ea315805440c3e100da9c",
          "message": "Merge pull request #957 from galacticusorg/fixSubhaloValidationStats\n\nIncrease the number of realizations for some subhalo validation models",
          "timestamp": "2025-10-26T04:30:42Z",
          "tree_id": "e23aaa5814520f99aebe825fab83e68f5250201f",
          "url": "https://github.com/galacticusorg/galacticus/commit/ee1a915ee5a43a63d05ea315805440c3e100da9c"
        },
        "date": 1761478807983,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4382788710366126,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 3.3251186976343705,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 43.58429148674663,
            "unit": "-logℒ"
          }
        ]
      }
    ]
  }
}