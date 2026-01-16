window.BENCHMARK_DATA = {
  "lastUpdate": 1768597519208,
  "repoUrl": "https://github.com/galacticusorg/galacticus",
  "entries": {
    "Dark matter-only subhalos benchmarks (COZMIC Milky Way WDM 5keV resolutionX1)": [
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
        "date": 1760996779642,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4801859109506597,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.7861351436762738,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 106.3946125531284,
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
        "date": 1761037639360,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.469341391573151,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.7943045694509476,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 107.83147960871173,
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
        "date": 1761097407428,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4586007974418838,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.706218298579079,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 105.08574753971837,
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
        "date": 1761168043041,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4198946906138203,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.8096540694346461,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 101.77612585356441,
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
        "date": 1761277071074,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.495000012039342,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.7819658472703037,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 101.46535858428597,
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
        "date": 1761478797720,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4780019188499882,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.7745198403106746,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 14.468308391488288,
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
        "date": 1762486843867,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.4257270113078588,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.12002814923446081,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 13.567428600942037,
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
        "date": 1762815720532,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.6333871030750182,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.0570231993148238,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 27.60018129074501,
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
        "date": 1763688341190,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.6312775347439707,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.9443541076633362,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 27.16172816424804,
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
        "date": 1763767800371,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.602332718046354,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.0819450365931487,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 26.191231878248512,
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
        "date": 1764714627056,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.6085548693794585,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.036745278614593,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 26.777325306687068,
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
        "date": 1765130141930,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.6519982389758363,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.065490910691278,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 26.62662768256888,
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
        "date": 1765585398033,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.6134618605301938,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.0715123388118717,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 26.026486152491422,
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
        "date": 1766179973375,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.6497691621453834,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.0246073388160224,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 26.715692519568158,
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
        "date": 1766302800994,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.6459666405387763,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.958138405227015,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 26.16422195812043,
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
        "date": 1766485342016,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.6421272653866343,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.0187558362688063,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 27.043097156556403,
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
        "date": 1767151230914,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.6319666425378867,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.958280586703522,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 25.600237753661,
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
        "date": 1767344585936,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.6441043346693278,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.9563794105201464,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 26.60990094629248,
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
        "date": 1767494937462,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.5915432124529427,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.0892126688624777,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 25.92283177614482,
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
        "date": 1768277182346,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.6223450406964823,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.1835896207240117,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 25.081971181748244,
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
        "date": 1768597518348,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.6049000807995395,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.0040519593136454,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:5keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 25.83744874202849,
            "unit": "-logℒ"
          }
        ]
      }
    ]
  }
}