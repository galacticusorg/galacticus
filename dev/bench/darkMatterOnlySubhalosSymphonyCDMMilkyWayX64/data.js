window.BENCHMARK_DATA = {
  "lastUpdate": 1761478783985,
  "repoUrl": "https://github.com/galacticusorg/galacticus",
  "entries": {
    "Dark matter-only subhalos benchmarks (Symphony Milky Way resolutionX64)": [
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
        "date": 1760996766316,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.7976931348623154e+302,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 13.479846199676029,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 1.7976931348623156e+292,
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
        "date": 1761037626352,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.7976931348623154e+302,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 18.349774340790727,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 1.7976931348623156e+292,
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
        "date": 1761097394316,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.7976931348623154e+302,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 13.479846199676029,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 1.7976931348623156e+292,
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
        "date": 1761168027355,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.7976931348623154e+302,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 16.13982235146408,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 1.7976931348623156e+292,
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
        "date": 1761277059312,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.7976931348623154e+302,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 18.349774340790727,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 1.7976931348623156e+292,
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
        "date": 1761478783136,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 4.743324844465189,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.176948316773921,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 400.6022114728399,
            "unit": "-logℒ"
          }
        ]
      }
    ]
  }
}