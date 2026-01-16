window.BENCHMARK_DATA = {
  "lastUpdate": 1768597505663,
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
        "date": 1762486830660,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 1.7976931348623154e+302,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.03479387721847,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 317.0274207936602,
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
        "date": 1762815707743,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 7.754934371852208,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.763163991916394,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 248.20962751029316,
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
        "date": 1763688328176,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 7.577074360276393,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.593055953429549,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 257.38247159469574,
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
        "date": 1763767788413,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 7.328734762678806,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.7697311077807045,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 254.93504278722526,
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
        "date": 1764714611779,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 7.675290983223745,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.823379779723212,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 251.6735614526815,
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
        "date": 1765130127127,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 7.336917140736256,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.980095205745266,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 257.6531614170191,
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
        "date": 1765585373283,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 7.336917140736256,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.980095205745209,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 257.653161417019,
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
        "date": 1766179958859,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 7.404469833683873,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.2207552835193995,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 255.75697675229503,
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
        "date": 1766302789544,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 7.144614532445136,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.0522180174468305,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 255.79749847865375,
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
        "date": 1766485326919,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 7.289190804166218,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.670477200830054,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 255.1605276554006,
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
        "date": 1767151217697,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 8.779034653822466,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.031263718389938,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 255.5358477021642,
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
        "date": 1767344573607,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 7.7996783683031765,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.978003945693573,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 251.4338282429502,
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
        "date": 1767494925281,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 7.144614532445136,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.0522180174468305,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 255.79749847865375,
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
        "date": 1768277166124,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 7.4603397261965725,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.873990160691645,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 255.73020295396185,
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
        "date": 1768597505055,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 7.459705447273258,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.869930003609682,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (Symphony CDM resolution X64 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 255.6946188568181,
            "unit": "-logℒ"
          }
        ]
      }
    ]
  }
}