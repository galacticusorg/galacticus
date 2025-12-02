window.BENCHMARK_DATA = {
  "lastUpdate": 1764714658201,
  "repoUrl": "https://github.com/galacticusorg/galacticus",
  "entries": {
    "Dark matter-only subhalos benchmarks (COZMIC Milky Way WDM 5keV resolutionX8)": [
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
        "date": 1760996805478,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 21.015019987678674,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 83.07044476812875,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761037665693,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 20.115951757342383,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 81.97619861280694,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761097433595,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 24.08128785729561,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 94.31755975688495,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761168073381,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 22.477485309948978,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 90.45904329692159,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761277094487,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 19.749453303775137,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 85.80466944994843,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761478828149,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 21.976092663460552,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 84.63025308540752,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 122.72717221947354,
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
          "id": "dfb8b9b2de370a6ef3a65666d4830d920f4a5713",
          "message": "Merge pull request #962 from galacticusorg/featSEDErrorMessage\n\nAdd improve error handling for SED integration",
          "timestamp": "2025-11-06T15:40:33Z",
          "tree_id": "44e9b6b02b5390b9fa323156b189f7dd7d7b46af",
          "url": "https://github.com/galacticusorg/galacticus/commit/dfb8b9b2de370a6ef3a65666d4830d920f4a5713"
        },
        "date": 1762486867364,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 22.037844715828335,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 84.82742026836235,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 122.29778821365795,
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
          "id": "54595bd309e540a815fada1b722970ca9a7b48d6",
          "message": "Merge pull request #964 from galacticusorg/fixSubhaloMassFunctionEmptyBins\n\nIgnore empty bins in subhalo mass function for Symphony X64 validation model",
          "timestamp": "2025-11-10T15:55:14Z",
          "tree_id": "abcd863f7d8cf137d7221fb00bae449fd0d43cca",
          "url": "https://github.com/galacticusorg/galacticus/commit/54595bd309e540a815fada1b722970ca9a7b48d6"
        },
        "date": 1762815746789,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 22.335207660976295,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 84.95680938681161,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 125.62916674204952,
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
          "id": "13441a27d008de4ecc507342cf4a81f59bd52a6b",
          "message": "Merge pull request #969 from galacticusorg/fixObsoletedParameters\n\nRemove obsoleted parameters from test suite files",
          "timestamp": "2025-11-20T15:48:46Z",
          "tree_id": "1f17d57d08c3096f32cd7f8b5ca3c3ae02e72875",
          "url": "https://github.com/galacticusorg/galacticus/commit/13441a27d008de4ecc507342cf4a81f59bd52a6b"
        },
        "date": 1763688367057,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 23.716684436620383,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 89.27660473833812,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 127.034457575858,
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
          "id": "e32181bf27fc720d78222dc337fbc680730fb48d",
          "message": "Merge pull request #970 from galacticusorg/featPostprocessing\n\nAdd postprocessing functionality",
          "timestamp": "2025-11-21T17:03:14Z",
          "tree_id": "21302a250d4dd2e10fb9bbb3e04fe3385f70093f",
          "url": "https://github.com/galacticusorg/galacticus/commit/e32181bf27fc720d78222dc337fbc680730fb48d"
        },
        "date": 1763767824418,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 23.72935887981499,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 88.66761737891122,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 131.35116053053707,
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
          "id": "2473ecce8a8e6d5be099483957f8dfe86943f8d5",
          "message": "Merge pull request #976 from galacticusorg/fixLightconeInterpolation\n\nUse the correct expansion factor when converting from comoving to physical coordinates",
          "timestamp": "2025-12-02T15:36:04Z",
          "tree_id": "07087afd048ec27cc7798e60526d2e712bcc368b",
          "url": "https://github.com/galacticusorg/galacticus/commit/2473ecce8a8e6d5be099483957f8dfe86943f8d5"
        },
        "date": 1764714656912,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 23.54333554452427,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 87.11029019214978,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 127.00395758425589,
            "unit": "-logℒ"
          }
        ]
      }
    ]
  }
}