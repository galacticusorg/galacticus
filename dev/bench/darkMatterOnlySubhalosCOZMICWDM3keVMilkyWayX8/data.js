window.BENCHMARK_DATA = {
  "lastUpdate": 1782334671097,
  "repoUrl": "https://github.com/galacticusorg/galacticus",
  "entries": {
    "Dark matter-only subhalos benchmarks (COZMIC Milky Way WDM 3keV resolutionX8)": [
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
        "date": 1760996796896,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 47.45148097930923,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 69.62094569995877,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761037656792,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 45.05640531637072,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 67.40129376469717,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761097424807,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 45.542259987668864,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 69.81873427638199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761168063342,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 45.57842475755169,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 70.13740946952583,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761277086192,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 46.7407434564065,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 70.81617112150276,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761478817780,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 47.975410840605846,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 71.4028957130627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 702.5955105157317,
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
        "date": 1762486859223,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 45.30113915643369,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 70.07354178075005,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 685.1932622530119,
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
        "date": 1762815738474,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 45.56796597512482,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 69.83727906761335,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 716.782220124175,
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
        "date": 1763688358345,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 47.67870356959134,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 70.55897112880288,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 714.5974612893018,
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
        "date": 1763767816368,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.328157270649825,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 70.62015841167697,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 718.5239489845014,
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
        "date": 1764714647187,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 49.889583350642205,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 74.36536479370525,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 713.7877064790192,
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
          "id": "49cf1a71a19accbd4cf914d402e3a33f069f7bb5",
          "message": "Merge pull request #977 from galacticusorg/featStrictParameters\n\nAdd an option for strict parameter processing",
          "timestamp": "2025-12-07T04:48:47Z",
          "tree_id": "66e3165b875a6d9daaf7233dd46852a90bb1efc0",
          "url": "https://github.com/galacticusorg/galacticus/commit/49cf1a71a19accbd4cf914d402e3a33f069f7bb5"
        },
        "date": 1765130161766,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.59911205351246,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 74.05200115546667,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 712.5947040132327,
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
          "id": "2b1a2183c2aa3152b395d6eb85f0e3965064b0bc",
          "message": "Merge pull request #979 from yzhaoastro/heating\n\nUpdate FDM core–halo relation and add validation test",
          "timestamp": "2025-12-12T16:59:27Z",
          "tree_id": "e347447241a24f9f8f18e6cca88c4b55e91b910f",
          "url": "https://github.com/galacticusorg/galacticus/commit/2b1a2183c2aa3152b395d6eb85f0e3965064b0bc"
        },
        "date": 1765585418192,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.518309604111536,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 73.37118145517896,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 695.3300548874421,
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
          "id": "0ce30cf062010a1b1f0690f6f59939321afc6af7",
          "message": "Merge pull request #985 from galacticusorg/featMergerTreeExport\n\nImprove merger tree export",
          "timestamp": "2025-12-19T14:29:29Z",
          "tree_id": "2c51e6711564a3da558ab0596c7209f2f0db9d67",
          "url": "https://github.com/galacticusorg/galacticus/commit/0ce30cf062010a1b1f0690f6f59939321afc6af7"
        },
        "date": 1766179992201,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 46.75956475696381,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 69.93234189790707,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 678.1054479124455,
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
          "id": "372971d60deb4602455fde9f3cfe24ec3f41ac2c",
          "message": "Merge pull request #986 from galacticusorg/fixRefactorClassInterface\n\nRefactor the CLASS interface code",
          "timestamp": "2025-12-21T00:42:29Z",
          "tree_id": "9eca5b304fe04dc6b1435b0da975e351147d57fd",
          "url": "https://github.com/galacticusorg/galacticus/commit/372971d60deb4602455fde9f3cfe24ec3f41ac2c"
        },
        "date": 1766302816346,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 45.7791186436787,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 69.27517333204989,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 716.9506799600475,
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
          "id": "b8361eea84b8e7b35d57b0c9a3e6262055dc458c",
          "message": "Merge pull request #980 from galacticusorg/featOutputAnalysisGalaxySizes\n\nAdd an output analysis for galaxy sizes as a function of stellar mass",
          "timestamp": "2025-12-23T03:56:10Z",
          "tree_id": "b10ce957b37aa8a05c89ea0e227a7e58e50c76e7",
          "url": "https://github.com/galacticusorg/galacticus/commit/b8361eea84b8e7b35d57b0c9a3e6262055dc458c"
        },
        "date": 1766485362100,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 49.37650567254301,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 75.62023585789667,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 711.5255930232892,
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
          "id": "5992c02a80404d08b56dd616df76f7b88c5bde9b",
          "message": "Merge pull request #988 from cgannonucm/master\n\nfeat: Implement tidal \"nodePropertyExtractorTidalField\"",
          "timestamp": "2025-12-30T20:22:03Z",
          "tree_id": "0025540a37d8a03003eb2841a3cab52bb681f28f",
          "url": "https://github.com/galacticusorg/galacticus/commit/5992c02a80404d08b56dd616df76f7b88c5bde9b"
        },
        "date": 1767151246941,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 44.2611381516301,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 68.48613233909846,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 722.433370152796,
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
          "id": "3663cfdcec92fc7d14cf19eba4940f9570bf82bb",
          "message": "Merge pull request #989 from galacticusorg/featEmptyConstructorDestructor\n\nAdd static analysis to detect empty constructors/destructors",
          "timestamp": "2026-01-02T02:11:44Z",
          "tree_id": "6e0992befdab4277362e358f3753faf2b7b2a789",
          "url": "https://github.com/galacticusorg/galacticus/commit/3663cfdcec92fc7d14cf19eba4940f9570bf82bb"
        },
        "date": 1767344602390,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.50532508210253,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 71.40416893795332,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 720.4831166287656,
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
          "id": "681a3d8f9320bd6b4a27dd44b993aaf4252790da",
          "message": "Merge pull request #991 from galacticusorg/featBump2026\n\nBump copyright year to 2026",
          "timestamp": "2026-01-03T19:48:29Z",
          "tree_id": "d7813ec3cfebc382df1caa719ef7363127a53ad2",
          "url": "https://github.com/galacticusorg/galacticus/commit/681a3d8f9320bd6b4a27dd44b993aaf4252790da"
        },
        "date": 1767494953421,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 46.85937366107072,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 69.76883456462427,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 690.0821363949452,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "distinct": true,
          "id": "e314c1a6751cf8d2221a38e051cfd5ec6bd25d49",
          "message": "fix(style): Formatting only",
          "timestamp": "2026-01-12T07:30:58-08:00",
          "tree_id": "64ecdb654c47b8637646b211eaeb8102d9894c24",
          "url": "https://github.com/galacticusorg/galacticus/commit/e314c1a6751cf8d2221a38e051cfd5ec6bd25d49"
        },
        "date": 1768277203982,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 47.40739025305113,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 70.25147735278297,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 688.7442228499705,
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
          "id": "5ec48cc7307039e6c5a31cc169fb32e62402f8d9",
          "message": "Merge pull request #992 from yzhaoastro/tidal-merge\n\nAdd solitonNFW heating model",
          "timestamp": "2026-01-16T14:12:17Z",
          "tree_id": "cfb0e600067218dbc60dd48069e2e49290666bb3",
          "url": "https://github.com/galacticusorg/galacticus/commit/5ec48cc7307039e6c5a31cc169fb32e62402f8d9"
        },
        "date": 1768597535803,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 49.510368382183316,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 75.72230816735751,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 697.4411524454847,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "distinct": true,
          "id": "e7835ec31f6c58d0209161ff1e6d43195e90448b",
          "message": "fix: Replace broken link with Internet Archive copy",
          "timestamp": "2026-01-20T07:23:27-08:00",
          "tree_id": "77917c8e549676617f812aa4a44b9cb743b1f285",
          "url": "https://github.com/galacticusorg/galacticus/commit/e7835ec31f6c58d0209161ff1e6d43195e90448b"
        },
        "date": 1768946575670,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 44.43382825300437,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 68.50357061641131,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 719.1561911421618,
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
          "id": "de968065160e1fe4744ab2b62b718636e4f5247a",
          "message": "Merge pull request #995 from galacticusorg/fixLightconeCrossingBugs\n\nFix bugs in determination of lightcone crossing",
          "timestamp": "2026-01-22T03:47:11Z",
          "tree_id": "077bc7dc6e0b7440af306afeb070ca94ffd730b0",
          "url": "https://github.com/galacticusorg/galacticus/commit/de968065160e1fe4744ab2b62b718636e4f5247a"
        },
        "date": 1769076972784,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 47.66314672355689,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 70.74452420811636,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 685.3232243709866,
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
          "id": "2b31d251742f809d3c94017b9fad2b7901ab7132",
          "message": "Merge pull request #996 from galacticusorg/fixMergerTreeWeightScale\n\nScale weights in all trees in a forest",
          "timestamp": "2026-01-31T05:39:36Z",
          "tree_id": "a23781d6b2fd048c06bd0f77a30a74b0810e61c3",
          "url": "https://github.com/galacticusorg/galacticus/commit/2b31d251742f809d3c94017b9fad2b7901ab7132"
        },
        "date": 1769863223642,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 49.95321892728718,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 75.77880014048357,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 723.5580229979253,
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
          "id": "f3934061f79a8904e4b9d18349853645e1df57dc",
          "message": "Merge pull request #997 from galacticusorg/fixAGNEmissionLines\n\nCorrect normalization of AGN emission lines",
          "timestamp": "2026-02-04T02:49:03Z",
          "tree_id": "d6d7a42150ad5ae480f94381748d839a1b64fb33",
          "url": "https://github.com/galacticusorg/galacticus/commit/f3934061f79a8904e4b9d18349853645e1df57dc"
        },
        "date": 1770223916163,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.91867476524538,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 71.05824542960944,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 722.8534514373971,
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
          "id": "d62cf7bbc6d30c7bdbb6a224356cbe61aff67938",
          "message": "Merge pull request #998 from galacticusorg/fixCloudyCIELocks\n\nUse `dynamic` datasets path for Cloudy table lock files",
          "timestamp": "2026-02-09T15:36:52Z",
          "tree_id": "f99fe40018df6727ab2f8ce343f4401a78b77dbf",
          "url": "https://github.com/galacticusorg/galacticus/commit/d62cf7bbc6d30c7bdbb6a224356cbe61aff67938"
        },
        "date": 1770698074289,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 45.483632671429525,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 67.62575691830747,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 797.4567742359702,
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
          "id": "a366adc5a89e34f4faa06751d26c0f7fa6535bec",
          "message": "Merge pull request #999 from galacticusorg/fixCloudyTableExtrapolation\n\nAllow extrapolation to smaller ages when integrating over Cloudy emission line tables",
          "timestamp": "2026-02-13T16:27:41Z",
          "tree_id": "96ac99849fc179140b63ec5bcaab6294a8b79eca",
          "url": "https://github.com/galacticusorg/galacticus/commit/a366adc5a89e34f4faa06751d26c0f7fa6535bec"
        },
        "date": 1771024645432,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 50.44654434191357,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 75.89133405974047,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 721.5227222258975,
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
          "id": "62bdb858689c23aece580508e1399efc7b67a23a",
          "message": "Merge pull request #1000 from galacticusorg/fixHopins2007LockLocation\n\nUse the `dynamic` datasets path for a lock file",
          "timestamp": "2026-02-17T15:15:20Z",
          "tree_id": "704d7ea1e4022aa578979412e960d11685d80f9c",
          "url": "https://github.com/galacticusorg/galacticus/commit/62bdb858689c23aece580508e1399efc7b67a23a"
        },
        "date": 1771385621639,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.113819394607795,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 70.93275774893947,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 710.566228260795,
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
          "id": "2f1bed7f74d0a1eacd8595a1282cebe480a31b8f",
          "message": "Merge pull request #1003 from galacticusorg/fixJohnson2021WithRegridding\n\nHandle parent-child halo pairs with identical masses in the Johnson et al. (2021) concentration model",
          "timestamp": "2026-02-21T06:52:46Z",
          "tree_id": "e4922a19dc04e710348d508c5243ba5d6afbfaa8",
          "url": "https://github.com/galacticusorg/galacticus/commit/2f1bed7f74d0a1eacd8595a1282cebe480a31b8f"
        },
        "date": 1771682521896,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 44.46515399854183,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 67.77176817626984,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 794.2649862566267,
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
          "id": "c6ab017920ddac6869393e6c6ae981c1952e7f73",
          "message": "Merge pull request #1004 from galacticusorg/boundMassInitialization\n\nAllow initializing satellite’s bound mass in a more flexible way",
          "timestamp": "2026-02-27T17:21:54Z",
          "tree_id": "f1d773cf1212928165af0df2c5904ce7d160aed4",
          "url": "https://github.com/galacticusorg/galacticus/commit/c6ab017920ddac6869393e6c6ae981c1952e7f73"
        },
        "date": 1772248589711,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.88488208490951,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 70.82277286264164,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 706.345725429209,
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
          "id": "f8fb4af0a7b8d6300c636159e4db5c41cd35bf28",
          "message": "Merge pull request #1005 from galacticusorg/fixParticleGeneratorBug\n\nFix a bug related to drawing the coordinates of the N-body particles",
          "timestamp": "2026-02-28T17:20:10Z",
          "tree_id": "4587e9fd209a195018ecc2e7b39d349513a814a2",
          "url": "https://github.com/galacticusorg/galacticus/commit/f8fb4af0a7b8d6300c636159e4db5c41cd35bf28"
        },
        "date": 1772321945992,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.90918402184421,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 70.79251704663777,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 703.6390551628251,
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
          "id": "18d6c3f638b0e215939598b4579dcd2c4863fd97",
          "message": "Merge pull request #1009 from galacticusorg/fixEvolveStatusPrivate\n\nUse a private `status` variable when evolving merger trees",
          "timestamp": "2026-03-03T15:17:54Z",
          "tree_id": "310c3daf949327a3fc535e0704cd3cdd71962516",
          "url": "https://github.com/galacticusorg/galacticus/commit/18d6c3f638b0e215939598b4579dcd2c4863fd97"
        },
        "date": 1772575840325,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 43.90918974626865,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 68.58833062908143,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 710.7180549358065,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "Andrew Benson",
            "username": "abensonca",
            "email": "abensonca@gmail.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "18d6c3f638b0e215939598b4579dcd2c4863fd97",
          "message": "Merge pull request #1009 from galacticusorg/fixEvolveStatusPrivate\n\nUse a private `status` variable when evolving merger trees",
          "timestamp": "2026-03-03T15:17:54Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/18d6c3f638b0e215939598b4579dcd2c4863fd97"
        },
        "date": 1772689690142,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.92482483625811,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 76.42611952112604,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 706.589164767238,
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
          "id": "18d733d92fa17af9767ad8d445f75613faae7bec",
          "message": "Merge pull request #1011 from galacticusorg/claude/standardize-monospace-formatting-wtb7l\n\nReplace {\\normalfont \\ttfamily} with \\mono{} macro in documentation",
          "timestamp": "2026-03-08T21:58:01Z",
          "tree_id": "86d3f7dbee77cf3de35d051c4b1d06ba7f60a647",
          "url": "https://github.com/galacticusorg/galacticus/commit/18d733d92fa17af9767ad8d445f75613faae7bec"
        },
        "date": 1773066849076,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 49.06132605911826,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 75.27855930713966,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 693.8791293686503,
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
          "id": "55c3949b9b3b3d77ced790d7150e1b15640557d7",
          "message": "Merge pull request #1019 from galacticusorg/dependabot/github_actions/docker/build-push-action-7\n\nbuild(deps): bump docker/build-push-action from 6 to 7",
          "timestamp": "2026-03-12T15:53:30Z",
          "tree_id": "7d465a628881c1f165ed62680c02a2d41d9ce9c1",
          "url": "https://github.com/galacticusorg/galacticus/commit/55c3949b9b3b3d77ced790d7150e1b15640557d7"
        },
        "date": 1773392405169,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 44.255423010420024,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 68.8923620582117,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 715.3012486813005,
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
          "id": "086b0bf51ce385e409c2fabac3cc9b8f271b0405",
          "message": "Merge pull request #1027 from galacticusorg/depup/cmbant/CAMB\n\nfix(deps): update cmbant/CAMB to 1.6.6",
          "timestamp": "2026-03-14T05:32:10Z",
          "tree_id": "370e497f5d3069dfc6b5bf0d39b6a63db7674bec",
          "url": "https://github.com/galacticusorg/galacticus/commit/086b0bf51ce385e409c2fabac3cc9b8f271b0405"
        },
        "date": 1773499777639,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 50.17476528129611,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 74.8216516898927,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 723.9659517689574,
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
          "id": "5cb357d0be452a1de0a9833b1734a7094017fdd8",
          "message": "Merge pull request #1028 from galacticusorg/claude/perl-to-python-port-3rZpu\n\nMigrate auxiliary scripts from Perl to Python",
          "timestamp": "2026-03-15T15:41:01Z",
          "tree_id": "df818e76c2fd19f939a55c634a1c078c4209f788",
          "url": "https://github.com/galacticusorg/galacticus/commit/5cb357d0be452a1de0a9833b1734a7094017fdd8"
        },
        "date": 1773647610884,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.997134334368894,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 74.84379989012947,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 704.8787756063059,
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
          "id": "72da9a5daf2dd07d5c5afaf4faac671c5b90658a",
          "message": "Merge pull request #1031 from galacticusorg/featAdjustableTimeOfCollapseTabulation\n\nAllow the resolution of the time-of-collapse table to be adjustable via a parameter",
          "timestamp": "2026-03-17T01:19:01Z",
          "tree_id": "75635a1b0f25c95df0875f61729c6446a580ea46",
          "url": "https://github.com/galacticusorg/galacticus/commit/72da9a5daf2dd07d5c5afaf4faac671c5b90658a"
        },
        "date": 1773734957858,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 50.31487350029792,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 76.05247018361857,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 710.976835828633,
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
          "id": "ad099f281dea4e336fb66ff7b7586f2aa1813b29",
          "message": "Merge pull request #1029 from galacticusorg/copilot/add-schema-for-directives\n\nAdd XSD schemas for 21 previously unvalidated directives",
          "timestamp": "2026-03-20T14:54:01Z",
          "tree_id": "f3b311233a9a5085aca647a84294924423c26bb4",
          "url": "https://github.com/galacticusorg/galacticus/commit/ad099f281dea4e336fb66ff7b7586f2aa1813b29"
        },
        "date": 1774071015163,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 44.596227498778276,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 66.875503467912,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 773.4800858270746,
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
          "id": "00dcebaab02563e2d87be5078ecf40aa6cf15d29",
          "message": "Merge pull request #1037 from galacticusorg/featReplicateFullySpecified\n\nAllow multiple realizations of fully-specified merger trees to be generated",
          "timestamp": "2026-03-22T19:36:33Z",
          "tree_id": "379bba0723fb30eb4a71ffc615490789f1c94d23",
          "url": "https://github.com/galacticusorg/galacticus/commit/00dcebaab02563e2d87be5078ecf40aa6cf15d29"
        },
        "date": 1774232510309,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 44.76204539463324,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 69.06863429512627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 731.4757742430502,
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
          "id": "4c7b97a65ab8bf017036bd8062ee61eb3217a70f",
          "message": "Merge pull request #1040 from galacticusorg/featCollapseHDF5Diffs\n\nMake HDF5 file diff PR comments collapsible",
          "timestamp": "2026-03-26T01:41:19Z",
          "tree_id": "2d48a657b6c57353b2dd7a823d91db9d0e5e53b8",
          "url": "https://github.com/galacticusorg/galacticus/commit/4c7b97a65ab8bf017036bd8062ee61eb3217a70f"
        },
        "date": 1774516384422,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 47.88639770782956,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 70.38490790255166,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 700.0011118100068,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abenson@carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "abenson@carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "distinct": true,
          "id": "a80676f14060689c7c3647752a25a494713a27f4",
          "message": "fix(perf): Optimize the Python `queuemanager` class\n\nFor Slurm managers, prioritize submitting jobs over postprocessing them (to ensure that the queue is always as full as possible). When retrieving status information on finished jobs, do it in a single batch for all finished jobs, instead of one at a time.",
          "timestamp": "2026-03-26T13:46:33-07:00",
          "tree_id": "71ae2e8fe961bcb42d8bfced4ec4f9a258538fe0",
          "url": "https://github.com/galacticusorg/galacticus/commit/a80676f14060689c7c3647752a25a494713a27f4"
        },
        "date": 1774795864352,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 44.56850637509743,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 68.82049205310385,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 741.0572778215193,
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
          "id": "704502e230b94e4f5f20af7f48e096d7742a853b",
          "message": "Merge pull request #1051 from galacticusorg/featCountCopiesErrorReporting\n\nAdd some error reporting when the `countCopies()` method fails",
          "timestamp": "2026-04-01T14:25:27Z",
          "tree_id": "7b64c0165b34be86cebbc1f140c53dc2f40d6517",
          "url": "https://github.com/galacticusorg/galacticus/commit/704502e230b94e4f5f20af7f48e096d7742a853b"
        },
        "date": 1775079324746,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 43.95695198022223,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 66.03535946496714,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 778.6762945729755,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "distinct": true,
          "id": "439f3eceebc2a3db1c10363c57f3a910934f5660",
          "message": "fix: Merge branch 'claude/port-linkchecker-perl-python-9CtmA'",
          "timestamp": "2026-04-02T08:54:42-07:00",
          "tree_id": "332a8c9838a34624a4b599ac58524ca35e205fa4",
          "url": "https://github.com/galacticusorg/galacticus/commit/439f3eceebc2a3db1c10363c57f3a910934f5660"
        },
        "date": 1775391593684,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.0412682957806,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 70.63197956972753,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 724.633078931407,
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
          "id": "bbd2fdeb3488409651e9b6250f75a3ce920703b8",
          "message": "Merge pull request #1058 from galacticusorg/fixObsoletedScripts\n\nRemove obsoleted scripts",
          "timestamp": "2026-04-09T20:37:56Z",
          "tree_id": "f28e76f9412a59ad7dc8a43c78bd45fde1ba1ac8",
          "url": "https://github.com/galacticusorg/galacticus/commit/bbd2fdeb3488409651e9b6250f75a3ce920703b8"
        },
        "date": 1775793508507,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 49.78042662870587,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 73.29411559629602,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 744.1492264545302,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "Andrew Benson",
            "username": "abensonca",
            "email": "abenson@obs.carnegiescience.edu"
          },
          "committer": {
            "name": "Andrew Benson",
            "username": "abensonca",
            "email": "abenson@obs.carnegiescience.edu"
          },
          "id": "57a796d08f966b2354f92b36e553fdc0fddb3517",
          "message": "fix: Reduce resolution of a test model\n\nThe adapative star formation history dataset length test is slow to run. It runs much faster with lower resolution and still allows for dataset lengths to be tested.",
          "timestamp": "2026-04-10T22:04:25Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/57a796d08f966b2354f92b36e553fdc0fddb3517"
        },
        "date": 1775932970639,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 44.43217316386177,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 69.75066845004773,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 721.1100586417577,
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
          "id": "5c441d11c9474971575f552d6974cf2b49d0eae6",
          "message": "Merge pull request #1065 from galacticusorg/fixTidalTruncationRadii\n\nSet an infinite truncation radius for profiles that are precisely NFW (i.e. are not truncated)",
          "timestamp": "2026-04-13T14:44:44Z",
          "tree_id": "33b5f185688e8fbd43a8f7b3f2c34e84128c18ae",
          "url": "https://github.com/galacticusorg/galacticus/commit/5c441d11c9474971575f552d6974cf2b49d0eae6"
        },
        "date": 1776116267919,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 47.77063440670083,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 71.25761690516617,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 728.2272566751257,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "distinct": true,
          "id": "b7a9e3228e0641d6ac649d3c86680de36f730793",
          "message": "fix: Remove obsoleted Perl modules",
          "timestamp": "2026-04-21T09:22:24-07:00",
          "tree_id": "eda9a0713cb42cf0fad1f528d07bd585def03410",
          "url": "https://github.com/galacticusorg/galacticus/commit/b7a9e3228e0641d6ac649d3c86680de36f730793"
        },
        "date": 1776890919272,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 44.733297507363034,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 68.52186797282828,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 731.5940793072261,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "Andrew Benson",
            "username": "abensonca",
            "email": "abenson@obs.carnegiescience.edu"
          },
          "committer": {
            "name": "Andrew Benson",
            "username": "abensonca",
            "email": "abenson@obs.carnegiescience.edu"
          },
          "id": "b61d7087fd808659956cee81e17422f4f8991b3a",
          "message": "fix: Merge branch 'master' of github.com:galacticusorg/galacticus",
          "timestamp": "2026-04-24T16:09:58Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/b61d7087fd808659956cee81e17422f4f8991b3a"
        },
        "date": 1777090022762,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.23927165868746,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 72.17163498999987,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 707.8795272729797,
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
          "id": "c1a5bddef352f69f112bacd1b3d5707f451c115c",
          "message": "fix: Use correct secret names in notarizatiom workflow",
          "timestamp": "2026-04-25T21:35:20Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/c1a5bddef352f69f112bacd1b3d5707f451c115c"
        },
        "date": 1777177919598,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 45.19120478297126,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 68.46733388100078,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 722.6533583793206,
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
          "id": "e3727aa117f9afecb9008ea052a90e9d26255ced",
          "message": "Merge pull request #1088 from galacticusorg/claude/mpi-unique-node-ids-mUbcE\n\nAdd MPI-aware unique ID assignment for tree nodes",
          "timestamp": "2026-04-27T14:19:53Z",
          "tree_id": "0212181d39ece2f3617244d4e61fd9e6ba103fa5",
          "url": "https://github.com/galacticusorg/galacticus/commit/e3727aa117f9afecb9008ea052a90e9d26255ced"
        },
        "date": 1777338591849,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 44.78462270117357,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 69.60493487078031,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 723.2333082053694,
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
          "id": "309719ffd333e5238520098542e0d28f09beddd2",
          "message": "Merge pull request #1094 from galacticusorg/fixLMCLabel\n\nCorrect constrained MW+LMC tree parameter file to label the correct branch",
          "timestamp": "2026-04-29T15:52:01Z",
          "tree_id": "22578b4386cc64939d95ec9969fa17351c00657a",
          "url": "https://github.com/galacticusorg/galacticus/commit/309719ffd333e5238520098542e0d28f09beddd2"
        },
        "date": 1777533617019,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 44.912574219376765,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 67.67013338550814,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 698.5454936325876,
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
          "id": "0a2758ccb239b8fe8e3125f29d2c8afe416d742e",
          "message": "Merge pull request #1104 from galacticusorg/claude/refactor-disk-node-operator-X5uR8\n\nrefactor: migrate verySimple disk analytic solver to a nodeOperator",
          "timestamp": "2026-05-06T14:23:54Z",
          "tree_id": "bb26e4daecaff6740e945719b4b84a599da79a5d",
          "url": "https://github.com/galacticusorg/galacticus/commit/0a2758ccb239b8fe8e3125f29d2c8afe416d742e"
        },
        "date": 1778139660570,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 43.888861706790856,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 67.2120343960273,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 702.7320245727111,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "distinct": true,
          "id": "211b550a51594169472fd70565771a1f8fe6ced1",
          "message": "fix: Update comments to reflect names of Python-ported scripts",
          "timestamp": "2026-05-11T07:47:53-07:00",
          "tree_id": "422cdfb052de02fd13920d55835c42c3f564d181",
          "url": "https://github.com/galacticusorg/galacticus/commit/211b550a51594169472fd70565771a1f8fe6ced1"
        },
        "date": 1778543463197,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 49.7582720698572,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 75.49601192834989,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 730.5737905430428,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abenson@carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "abenson@carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "distinct": true,
          "id": "40dae4d6b28b3f643ae5fbde551c362eaeb94608",
          "message": "feat: Expand the parameter description for disk angular momentum factor\n\nNow includes a note on what this would be for a fully self-gravitating, razor-thin exponential disk.",
          "timestamp": "2026-05-14T10:17:06-07:00",
          "tree_id": "680e78977f28e87593aa4d48d3968c7acaa3a508",
          "url": "https://github.com/galacticusorg/galacticus/commit/40dae4d6b28b3f643ae5fbde551c362eaeb94608"
        },
        "date": 1778825571807,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.315366580141976,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 71.39252741592648,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 719.2898704682389,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "82462a1fbf030125f0b1538612bc433191ccf531",
          "message": "Merge pull request #1117 from galacticusorg/featValidatePkShape\n\nAdds validation of the `power` dataset",
          "timestamp": "2026-05-15T08:37:36-07:00",
          "tree_id": "5b06bd69f57b6bd8e7b15b42307be4c4db586ec1",
          "url": "https://github.com/galacticusorg/galacticus/commit/82462a1fbf030125f0b1538612bc433191ccf531"
        },
        "date": 1778940779958,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 46.07629898833202,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 68.95700118905756,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 797.7050652737028,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "Claude",
            "username": "claude",
            "email": "noreply@anthropic.com"
          },
          "committer": {
            "name": "Claude",
            "username": "claude",
            "email": "noreply@anthropic.com"
          },
          "id": "970a6b1c6c7d3a381cb92ee3da3ebe408343e098",
          "message": "fix: Use the deploy app token for the bleeding-edge tag push\n\nThe Update-Release-Tag job generated an app token and tried to use it by\nsetting GITHUB_TOKEN in the step env, but git itself doesn't read\nGITHUB_TOKEN — actions/checkout had already persisted the default\nworkflow token into .git/config as an http.extraheader, and that is what\ngit push was actually authenticating with.  Non-workflow tag pushes\nsucceeded because the default token carries contents:write; pushes whose\ncommit touched .github/workflows/* failed because that token never\ncarries workflows:write.\n\nGenerate the token before checkout and pass it via the checkout step's\n`token:` input so it is the one persisted into .git/config.  Also\ndeclare permission-contents: write explicitly — the v3 action narrows\nthe issued token to the permission-* inputs supplied, so workflows alone\nwould have dropped contents.",
          "timestamp": "2026-05-16T15:26:08Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/970a6b1c6c7d3a381cb92ee3da3ebe408343e098"
        },
        "date": 1778974498087,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 49.32412196432414,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 72.18522319988473,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 719.4198290919978,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "distinct": true,
          "id": "d6e53e0a189efe9404011acc080503a4f481eff7",
          "message": "chore: remove stale gfortran workaround note in lookup builder\n\nKeep the actual workaround since it actually fits with code style.",
          "timestamp": "2026-05-26T18:43:35-07:00",
          "tree_id": "decafb77a2b7dee2a2821e35cfeed7d1b31b7320",
          "url": "https://github.com/galacticusorg/galacticus/commit/d6e53e0a189efe9404011acc080503a4f481eff7"
        },
        "date": 1779870480956,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.459315708186864,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 72.19753915792775,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 720.7731057621187,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "distinct": true,
          "id": "ea894079180f7f0814daecb7380ed55a48532645",
          "message": "chore: remove stale gfortran workaround note in stellar populations\n\nKeep the actual workaround as it complies with codeing standards.",
          "timestamp": "2026-05-27T07:26:55-07:00",
          "tree_id": "0c29a2802f7cf482dabc796672eac3fab3ab4057",
          "url": "https://github.com/galacticusorg/galacticus/commit/ea894079180f7f0814daecb7380ed55a48532645"
        },
        "date": 1779948048203,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 46.813094663000065,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 69.40912055929225,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 785.5861525586361,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d5877a89b6b73b008e860f86df9ea90a5f5d3ba8",
          "message": "Merge pull request #1136 from galacticusorg/fixBHMSigmaReport\n\nSet default weights when decoding radius specifiers",
          "timestamp": "2026-05-29T16:52:07-07:00",
          "tree_id": "6fbea6c2521222f5881da4c52a04369922605075",
          "url": "https://github.com/galacticusorg/galacticus/commit/d5877a89b6b73b008e860f86df9ea90a5f5d3ba8"
        },
        "date": 1780134058218,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.14893697731221,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 70.30563026124146,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 709.6347419533101,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abenson@obs.carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b9b0b4917d1fcb2b05e0bb45212bc73419651e5e",
          "message": "Merge pull request #1010 from galacticusorg/dmConstraintPipeline\n\nAdd infrastructure for calibration of halo mass functions",
          "timestamp": "2026-06-03T16:22:28-07:00",
          "tree_id": "675d4dd1593eb741c844e671ec04811581cc203b",
          "url": "https://github.com/galacticusorg/galacticus/commit/b9b0b4917d1fcb2b05e0bb45212bc73419651e5e"
        },
        "date": 1780589178905,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 49.09476394940827,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 70.47435968047937,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 754.373220228267,
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
          "id": "8c2f8ff09d7d756ee144289d9ed7dada25eb30f5",
          "message": "Merge pull request #1157 from galacticusorg/claude/optimize-simulations-analysis-6C8sT\n\nAdd halo formation/crossing-time operators and hierarchical simulation-analysis layout",
          "timestamp": "2026-06-12T02:40:09Z",
          "tree_id": "403582a518c0c064c361ee28af6e8d0c88e62c67",
          "url": "https://github.com/galacticusorg/galacticus/commit/8c2f8ff09d7d756ee144289d9ed7dada25eb30f5"
        },
        "date": 1781276878879,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 43.922795497127844,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 63.74590425884215,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 748.621536518211,
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
          "id": "a6aa58c04a5287b97bf8bf96385be9984c201014",
          "message": "Merge pull request #1162 from galacticusorg/fix/fdm-soliton-heated-memory-leak\n\nfix(fdm): free prior massDistributionHeated_ before reallocating in solitonNFWHeated",
          "timestamp": "2026-06-15T14:55:04Z",
          "tree_id": "1b3c4687493bb018f79003ae174dbce4b4b770ff",
          "url": "https://github.com/galacticusorg/galacticus/commit/a6aa58c04a5287b97bf8bf96385be9984c201014"
        },
        "date": 1781574640061,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.596331991434326,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 69.40018002076702,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 697.2426779767901,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abenson@carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "abenson@carnegiescience.edu",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "distinct": true,
          "id": "931f580aca5a858da4c57ed11cf5a4bd7ae2b24c",
          "message": "docs(pool): refer to objectPool without broken cross-reference\n\nResolve undefined `\\refClass{objectPool}` to a plain `\\mono` reference\nsince the class lacks a corresponding cross-reference target.",
          "timestamp": "2026-06-16T15:12:57-07:00",
          "tree_id": "c350a980a9b204bb05b7646ff63744432822807d",
          "url": "https://github.com/galacticusorg/galacticus/commit/931f580aca5a858da4c57ed11cf5a4bd7ae2b24c"
        },
        "date": 1781679868041,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 44.42736420289374,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 64.23196739961803,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 743.40054986834,
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
          "id": "0f9ace2318ca44c96c4069b5833e5f2bbf94dff3",
          "message": "Merge pull request #1178 from galacticusorg/ci-macos-homebrew-gcc16\n\nci(buildMacOS): switch back to Homebrew GCC 16 install",
          "timestamp": "2026-06-20T04:56:26Z",
          "tree_id": "bfb3e63143b59a60d866c50f5e8c96e4121b5aea",
          "url": "https://github.com/galacticusorg/galacticus/commit/0f9ace2318ca44c96c4069b5833e5f2bbf94dff3"
        },
        "date": 1781963014492,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 43.55559003481266,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 64.95303778985003,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 713.3897020266232,
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
          "id": "e50f9e06c34e0511034666c205ab1cc4d1dfd554",
          "message": "Merge pull request #1179 from galacticusorg/docs/review-phase0\n\ndocs: documentation review — Phase 0/1 fixes (links, compiler version, docstring RST, extractor parameter-drop)",
          "timestamp": "2026-06-21T02:55:37Z",
          "tree_id": "50bda24deea511920857b26d5e63630bd9d84fe7",
          "url": "https://github.com/galacticusorg/galacticus/commit/e50f9e06c34e0511034666c205ab1cc4d1dfd554"
        },
        "date": 1782112794439,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.15252513803347,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 69.6465957504745,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 699.000536458748,
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
          "id": "55c235cce48356077131ec411afdea96a9cf4eac",
          "message": "Merge pull request #1183 from galacticusorg/dependabot/pip/sphinxcontrib-bibtex-gte-2.6.5\n\nchore(deps): update sphinxcontrib-bibtex requirement from >=2.6 to >=2.7.0",
          "timestamp": "2026-06-23T14:28:21Z",
          "tree_id": "06dcf18a2d65e8feffef5bbf9c041d9687527d12",
          "url": "https://github.com/galacticusorg/galacticus/commit/55c235cce48356077131ec411afdea96a9cf4eac"
        },
        "date": 1782258513176,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 48.26337438480359,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 69.00282509521493,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 707.8071100657734,
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
          "id": "72702669f0a79e7ee835e7e194f51b18e94cbc44",
          "message": "Merge pull request #1188 from galacticusorg/feature/metaproperty-getreference-pointer-getter\n\nperf: pointer-returning meta-property getter to avoid per-call array copy",
          "timestamp": "2026-06-24T14:24:56Z",
          "tree_id": "85a96f6c41c239f94657f088e34259d3daaed848",
          "url": "https://github.com/galacticusorg/galacticus/commit/72702669f0a79e7ee835e7e194f51b18e94cbc44"
        },
        "date": 1782334670255,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 47.98421188224765,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 68.93505113822809,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:3keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 721.0251421020442,
            "unit": "-logℒ"
          }
        ]
      }
    ]
  }
}