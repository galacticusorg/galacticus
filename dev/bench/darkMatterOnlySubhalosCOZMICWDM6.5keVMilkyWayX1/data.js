window.BENCHMARK_DATA = {
  "lastUpdate": 1767494946389,
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
        "date": 1762486851236,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.506712590055584,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.344180800961353,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 45.126816155122526,
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
        "date": 1762815729821,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4700806271905011,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.3852722860307785,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 43.697897636010055,
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
        "date": 1763688349578,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4449990223122477,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.3505418616359623,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 44.050078176505835,
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
        "date": 1763767808404,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4380005281683226,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 3.282293687516268,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 43.404799562346305,
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
        "date": 1764714637032,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4267828736472683,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.664941710603113,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 33.78518705788059,
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
        "date": 1765130151817,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4669064913948726,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.3371325904806772,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 42.88505234397701,
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
        "date": 1765585408350,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4504123403390428,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.3904483782719903,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 43.49997770569844,
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
        "date": 1766179982784,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4975411329711381,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.4716015771540671,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 43.33235278030962,
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
        "date": 1766302808576,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4311471558447824,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.6930413848921344,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 33.94282995046092,
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
        "date": 1766485352403,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.5089116086542882,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.324184550047825,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 43.17461786512376,
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
        "date": 1767151238622,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4628820482560774,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 3.216784885820383,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 43.38950098382144,
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
        "date": 1767344594320,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4665006647073944,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.3532317125251105,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 44.165194140624074,
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
        "date": 1767494945480,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.5019047280979982,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 3.2676952776629467,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:6.5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 44.65537117881257,
            "unit": "-logℒ"
          }
        ]
      }
    ]
  }
}