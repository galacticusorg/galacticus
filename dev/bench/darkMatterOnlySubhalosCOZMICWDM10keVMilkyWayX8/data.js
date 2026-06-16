window.BENCHMARK_DATA = {
  "lastUpdate": 1781574664255,
  "repoUrl": "https://github.com/galacticusorg/galacticus",
  "entries": {
    "Dark matter-only subhalos benchmarks (COZMIC Milky Way WDM 10keV resolutionX8)": [
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
        "date": 1760996818389,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 3.2283299096767992,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.041578905063373,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761037678622,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 3.3618634052729544,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.7976931348623154e+302,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761097447370,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 3.1908872975052143,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.08168404317912081,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761168088247,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 3.041486296952316,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.9888054498215357,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761277106324,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.8771127755416206,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 2.1318401246114025,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
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
        "date": 1761478842555,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 3.659376067681287,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.8322000505588489,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 85.59914427958104,
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
        "date": 1762486878929,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 3.676858286678377,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.010480618113493989,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 91.46793049585978,
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
        "date": 1762815759779,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.428291266676867,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.21651781577333473,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 108.62907997609582,
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
        "date": 1763688379753,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.386992830706112,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.6047197848087507,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 107.72302430643006,
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
        "date": 1763767836374,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.4613210847328535,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.24808167435039863,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 116.83640953799072,
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
        "date": 1764714672339,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.2737472519142123,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.6314703511291295,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 108.7344805348136,
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
        "date": 1765130186697,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.0193917893621167,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.22894703262141824,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 118.51361514177009,
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
        "date": 1765585443807,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.7321447760683488,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.2057827063211155,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 110.12505331544016,
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
        "date": 1766180015629,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.131128506593111,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.01091657792979106,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 110.87967010492432,
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
        "date": 1766302836286,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.1527814872916524,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.9624807613954358,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 108.20412844932888,
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
        "date": 1766485387148,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.8098021106086202,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.1545205555312208,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 114.80703217938152,
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
        "date": 1767151267925,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.3106735545999255,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.08310624072392425,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 118.49032243189848,
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
        "date": 1767344621836,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.121927927448173,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.009007542931576573,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 111.06434468330082,
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
        "date": 1767494973528,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.9874485206782715,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.2872824899000692,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 109.12771059099421,
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
        "date": 1768277229503,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.408577777535764,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.0065659097432968,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 104.77235730652825,
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
        "date": 1768597559331,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.2807003375566355,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.028656606292711828,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 107.71915855069534,
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
        "date": 1768946601438,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.3298292060830734,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.012873969143097241,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 108.67205647504575,
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
        "date": 1769076992637,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.3384730672842435,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.17765874397693882,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 105.83795194186659,
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
        "date": 1769863243488,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.445954852099816,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.12578210685973978,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 108.72262689242667,
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
        "date": 1770223937746,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.2532653529490902,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.620063105831381,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 113.75840261518002,
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
        "date": 1770698100968,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.8081684194331655,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.04454173528883043,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 111.38991425146578,
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
        "date": 1771024665247,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.388195928948666,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.7629871321899193,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 104.71112298916728,
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
        "date": 1771385642958,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.0871836080745423,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.08858757519630522,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 114.7310169936543,
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
        "date": 1771682543626,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.2898160852448903,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.8919614955684274,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 107.27956848149026,
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
        "date": 1772248616661,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.108454847521762,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.926895593654617,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 117.57592360846954,
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
        "date": 1772321969435,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.489968931645798,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.11741562933328353,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 114.47646406755625,
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
        "date": 1772575865208,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.1491839636727628,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.7903703411589449,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 109.1998033318965,
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
        "date": 1772689712797,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.7447959572658334,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.855111752254091,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 107.9543728870454,
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
        "date": 1773066873413,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.321411063374362,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.9483326777489367,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 117.32820091936924,
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
        "date": 1773392429977,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.03614548208889,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.260842544572224,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 109.19651505258287,
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
        "date": 1773499800067,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.4300349557094654,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.3153432825887199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 117.94882409036843,
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
        "date": 1773647638374,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.4798665004880958,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.0879774509999827,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 117.22258048922001,
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
        "date": 1773734982568,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.5086146620696734,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.0414201070027964,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 114.55591861829504,
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
        "date": 1774071043451,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.4484418462548785,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.04413437530213393,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 117.46946234252424,
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
        "date": 1774232530520,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.061499012131212,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.7582827308683919,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 106.30634957681534,
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
        "date": 1774516411408,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.440345802364927,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.6553703182826496,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 107.54767617318795,
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
        "date": 1774795885952,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.317574799849085,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.0660676045759967,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 106.57504797281682,
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
        "date": 1775079352153,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.0877827752391367,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.11099627324238792,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 117.09529594039074,
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
        "date": 1775391627089,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.115317541239518,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.0475228567594623,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 111.06880452331599,
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
        "date": 1775793531952,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.0173317419351067,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.15074308785820367,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 108.89433610878133,
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
        "date": 1775932992960,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.7741165650432327,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.06003753109366983,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 111.41197845698265,
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
        "date": 1776116290087,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.0627584976364357,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.16270807954705635,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 113.868677356272,
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
        "date": 1776890944074,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.4212266725127307,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.09221329710941983,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 106.5611221204002,
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
        "date": 1777090050326,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.028683807135079,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.03822320793437406,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 116.12502615661597,
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
        "date": 1777177948991,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.5198540432161316,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.7800835259657459,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 116.032930668355,
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
        "date": 1777338617327,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.0752320741412884,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.2626377937087101,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 117.97438620107314,
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
        "date": 1777533641808,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.0449656326558356,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.9814322049854489,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 117.2361092399704,
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
        "date": 1778139688308,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.033344380723473,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.051957357029309925,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 111.6676541835717,
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
        "date": 1778543493041,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.1270032528379685,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.1396205215630273,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 116.04934571956157,
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
        "date": 1778825594767,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.827989290819939,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.08618422723391372,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 112.5152640666631,
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
        "date": 1778940803295,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.139821241564805,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.08646317873727494,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 113.82639263613507,
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
        "date": 1778974522086,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.0717634678737538,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.05229317766628583,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 114.1100462311903,
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
        "date": 1779870508824,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.2748393296762153,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.23832285794884456,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 106.7660677452767,
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
        "date": 1779948073160,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.264999523175394,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.8522036740827534,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 107.39137636834694,
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
        "date": 1780134081495,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.3528201728692615,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 0.01305737738759305,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 107.76191754385644,
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
        "date": 1780589205268,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.8910963061722799,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.024587203162036,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 50.584918770686436,
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
        "date": 1781276912606,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 2.0176558136909026,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.2268905119508648,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 50.61819397767006,
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
        "date": 1781574663402,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.8948432603147742,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 1.0802062729499822,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:10keV resolution X8 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 47.4286724528893,
            "unit": "-logℒ"
          }
        ]
      }
    ]
  }
}