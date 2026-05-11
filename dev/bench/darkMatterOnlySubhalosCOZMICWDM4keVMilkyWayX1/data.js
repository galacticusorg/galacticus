window.BENCHMARK_DATA = {
  "lastUpdate": 1778543436216,
  "repoUrl": "https://github.com/galacticusorg/galacticus",
  "entries": {
    "Dark matter-only subhalos benchmarks (COZMIC Milky Way WDM 4keV resolutionX1)": [
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
        "date": 1760996775476,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.07630530094742882,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.128148512671542,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 89.2402703602242,
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
        "date": 1761037635227,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.20044598633400312,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.027486641207112,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 89.05407763841555,
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
        "date": 1761097402924,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.016430930071591465,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.003009299466298,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 65.84922839688284,
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
        "date": 1761168037240,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.1765671343309645,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 6.132133327221135,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 90.74919901336182,
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
        "date": 1761277067348,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.21707483247959125,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 7.011113956633588,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 89.61961222023209,
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
        "date": 1761478793088,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.0036917323875629293,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.071520689497477,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 52.407063092171484,
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
        "date": 1762486839301,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.06476149567582856,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.941437405544189,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 157.37872897161347,
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
        "date": 1762815716440,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.03256638673725376,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.865587351139611,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 121.33586889515043,
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
        "date": 1763688336773,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.13666211020207197,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.045606783399787,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 146.22568208047693,
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
        "date": 1763767796456,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.134266203195323,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.8388430233734265,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 141.3960528273522,
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
        "date": 1764714622112,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.024676064620933746,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.955699347229061,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 133.71245819786043,
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
        "date": 1765130137146,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.11731670917670556,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.016992401997484,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 135.22452864698647,
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
        "date": 1765585392626,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.12768098074293976,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.0064886406850455,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 148.4041111843438,
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
        "date": 1766179968268,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.018984127732819434,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.084493052594777,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 115.45035152082984,
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
        "date": 1766302797074,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.1354529792923891,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.012382882846042,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 140.61547570785527,
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
        "date": 1766485337079,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.1490474568279445,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.9705298994452445,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 145.78488885389214,
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
        "date": 1767151226539,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.13222558371473858,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.981792884401603,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 139.47309539706592,
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
        "date": 1767344581591,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.1412628543371477,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.989630630882198,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 144.11920948521475,
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
        "date": 1767494933872,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.11292564649304748,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.037923041322813,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 142.27424200365382,
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
        "date": 1768277176850,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.13325612363756245,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.98734093624099,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 141.7369426062628,
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
        "date": 1768597513683,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.07575722200876167,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.751008868203135,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 137.3374045286988,
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
        "date": 1768946550369,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.14657762333801005,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.903605110860477,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 148.99720993792434,
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
        "date": 1769076951646,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.13422300175564428,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.86144831763165,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 149.73029608347085,
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
        "date": 1769863203223,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.10099436284886787,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.188389691884989,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 150.63672174985197,
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
        "date": 1770223895635,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.024015689806104024,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.289183150715457,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 129.3467068420693,
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
        "date": 1770698047009,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.11032540266239577,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.89970368904422,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 154.69303730090792,
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
        "date": 1771024625513,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.1259562788533446,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.989177256245258,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 153.20726376623276,
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
        "date": 1771385599989,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.017653603864404066,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.936032760038726,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 117.13856598272913,
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
        "date": 1771682499594,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.12203972568221721,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.956824290747868,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 138.09828283003014,
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
        "date": 1772248563785,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.14269621410518796,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.003862937705749,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 141.99207802161988,
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
        "date": 1772321922461,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.0042529622688124435,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.098598759759767,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 123.20985701784195,
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
        "date": 1772575815003,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.03011985159434516,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.971483179317005,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 118.15685433245844,
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
        "date": 1772689667038,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.1257693073054036,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.910001470096573,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 150.24507420840362,
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
        "date": 1773066824278,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.13618748473048092,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.980064191318263,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 152.22066301529986,
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
        "date": 1773392380106,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.08173745720818626,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.825833939009321,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 137.8169534619577,
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
        "date": 1773499755665,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.13929704179596603,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.984388167235958,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 142.06628947461084,
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
        "date": 1773647583279,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.07316683454017048,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.855014376834069,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 125.7324308187961,
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
        "date": 1773734934925,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.12556063392790207,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.079876408652222,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 131.7816711789902,
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
        "date": 1774070987434,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.1197271090309735,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.999968202724082,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 142.12256728266814,
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
        "date": 1774232489829,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.07826054493878354,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.87420663427203,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 155.483299912865,
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
        "date": 1774516354804,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.10861629830922703,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.062833746265833,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 150.63278597696097,
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
        "date": 1774795842991,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.010496747353131108,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.15986147634294,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 112.05732556854542,
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
        "date": 1775079297421,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.020408413803347902,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 6.052696516192842,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 121.54715867739054,
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
        "date": 1775391566260,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.051303489084987386,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.165559673657903,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 127.46094536323645,
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
        "date": 1775793483790,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.13120881956576746,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.955326732301957,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 138.8866546674024,
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
        "date": 1775932948178,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.036217561858874814,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.9317990771830775,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 118.48161176979524,
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
        "date": 1776116245933,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.04279473596829042,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.064281937880528,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 130.56784366633673,
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
        "date": 1776890893966,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.12218347054551904,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.000246642389557,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 149.9250303756386,
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
        "date": 1777089995072,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.023299067432348664,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.048223829766535,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 129.69516316001352,
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
        "date": 1777177890043,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.1279738458719375,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.993788217445617,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 152.86297494301073,
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
        "date": 1777338567205,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.14645591535598024,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.911304995443242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 155.47858399214022,
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
        "date": 1777533592529,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.08890643172795554,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.979006674290915,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 143.53972413211704,
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
        "date": 1778139631763,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.10218350980875779,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 4.972144453551779,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 145.64851636513515,
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
        "date": 1778543435015,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloMassFunction",
            "value": 0.1296819034272696,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloRadialDistribution",
            "value": 5.837820222854896,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos (COZMIC WDM:4keV resolution X1 Milky Way) - Likelihood - subhaloVelocityMaximumMean",
            "value": 147.9232339629192,
            "unit": "-logℒ"
          }
        ]
      }
    ]
  }
}