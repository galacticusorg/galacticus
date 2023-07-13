window.BENCHMARK_DATA = {
  "lastUpdate": 1689226029383,
  "repoUrl": "https://github.com/galacticusorg/galacticus",
  "entries": {
    "Dark matter-only subhalos benchmarks": [
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
          "id": "d961f9f4b89ee922ef58bb7ef8b97b5b4b01cf22",
          "message": "Merge pull request #328 from galacticusorg/chandraIntegrals\n\nPass a half-mass radius to Chandrasekhar integral functions",
          "timestamp": "2022-11-14T22:19:06-08:00",
          "tree_id": "9e04dfbaec9029b954effe6798c60e4646535466",
          "url": "https://github.com/galacticusorg/galacticus/commit/d961f9f4b89ee922ef58bb7ef8b97b5b4b01cf22"
        },
        "date": 1668514022306,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.1223603836366,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.7948185912783,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26209.4879648531,
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
          "id": "6eab8997cd73cb0a474228ade542d133890ad138",
          "message": "fix: Avoid out-of-range error when integrating Vega spectrum\n\nThe interpolation of the Vega spectrum used for computing AB-Vega offsets previously had no `extrapolationType` defined, causing it to fail if integration extended outside of the tabulated range (e.g. for some ionizing luminosity filters). The `extrapolationType` is now set to `zero` so that zero flux is assumed outside of the tabulated range.",
          "timestamp": "2022-11-15T08:28:30-08:00",
          "tree_id": "f719a6b9bcd1a3c271446d6abae1156039d42e13",
          "url": "https://github.com/galacticusorg/galacticus/commit/6eab8997cd73cb0a474228ade542d133890ad138"
        },
        "date": 1668539615152,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 47.907,
            "unit": "seconds",
            "range": 0.666351333757287
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
          "id": "6eab8997cd73cb0a474228ade542d133890ad138",
          "message": "fix: Avoid out-of-range error when integrating Vega spectrum\n\nThe interpolation of the Vega spectrum used for computing AB-Vega offsets previously had no `extrapolationType` defined, causing it to fail if integration extended outside of the tabulated range (e.g. for some ionizing luminosity filters). The `extrapolationType` is now set to `zero` so that zero flux is assumed outside of the tabulated range.",
          "timestamp": "2022-11-15T08:28:30-08:00",
          "tree_id": "f719a6b9bcd1a3c271446d6abae1156039d42e13",
          "url": "https://github.com/galacticusorg/galacticus/commit/6eab8997cd73cb0a474228ade542d133890ad138"
        },
        "date": 1668539625629,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.4771891211537,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.0053889680488,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26391.8664916722,
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
          "id": "5707e0484802f86d2aa519dc4bbec27faa438fdf",
          "message": "fix: Correct parameters in reference models",
          "timestamp": "2022-11-15T21:34:00Z",
          "tree_id": "61be8a4fed82fa482177de47efa5eb0e8fed9488",
          "url": "https://github.com/galacticusorg/galacticus/commit/5707e0484802f86d2aa519dc4bbec27faa438fdf"
        },
        "date": 1668557708288,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 53.887,
            "unit": "seconds",
            "range": 0.424068508616203
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
          "id": "5707e0484802f86d2aa519dc4bbec27faa438fdf",
          "message": "fix: Correct parameters in reference models",
          "timestamp": "2022-11-15T21:34:00Z",
          "tree_id": "61be8a4fed82fa482177de47efa5eb0e8fed9488",
          "url": "https://github.com/galacticusorg/galacticus/commit/5707e0484802f86d2aa519dc4bbec27faa438fdf"
        },
        "date": 1668557718306,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.4008817695147,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9680154730719,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26281.0360034663,
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
          "id": "9bf9a9f6f3c9ba7625f3b4b4c4e7e97e6b871184",
          "message": "fix(style): Formatting only",
          "timestamp": "2022-11-17T00:07:42Z",
          "tree_id": "a8e3ff3a93a3e5f134df6ed443a4d4f296804d24",
          "url": "https://github.com/galacticusorg/galacticus/commit/9bf9a9f6f3c9ba7625f3b4b4c4e7e97e6b871184"
        },
        "date": 1668676844671,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.754,
            "unit": "seconds",
            "range": 0.098622512642903
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
          "id": "9bf9a9f6f3c9ba7625f3b4b4c4e7e97e6b871184",
          "message": "fix(style): Formatting only",
          "timestamp": "2022-11-17T00:07:42Z",
          "tree_id": "a8e3ff3a93a3e5f134df6ed443a4d4f296804d24",
          "url": "https://github.com/galacticusorg/galacticus/commit/9bf9a9f6f3c9ba7625f3b4b4c4e7e97e6b871184"
        },
        "date": 1668676851702,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.1223603836366,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.7948185912783,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26209.4879648531,
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
          "id": "84c486c856f3e700304334525bb0db6fd07d2634",
          "message": "Merge pull request #331 from galacticusorg/missingDestructsFix\n\nAdd missing destructors",
          "timestamp": "2022-11-17T07:12:34-08:00",
          "tree_id": "664050cf32daa91fd589a01d922132f179a807f7",
          "url": "https://github.com/galacticusorg/galacticus/commit/84c486c856f3e700304334525bb0db6fd07d2634"
        },
        "date": 1668709384625,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 50.991,
            "unit": "seconds",
            "range": 0.0517387668964034
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
          "id": "84c486c856f3e700304334525bb0db6fd07d2634",
          "message": "Merge pull request #331 from galacticusorg/missingDestructsFix\n\nAdd missing destructors",
          "timestamp": "2022-11-17T07:12:34-08:00",
          "tree_id": "664050cf32daa91fd589a01d922132f179a807f7",
          "url": "https://github.com/galacticusorg/galacticus/commit/84c486c856f3e700304334525bb0db6fd07d2634"
        },
        "date": 1668709394645,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.1526876081635,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8349286542639,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26574.7073451904,
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
          "id": "ce515f9e4525c3464161e5183cbd7bfdc56c4e0d",
          "message": "Merge pull request #332 from galacticusorg/componentSubParameters\n\nMake components read their options from their own subparameters",
          "timestamp": "2022-11-18T08:33:45-08:00",
          "tree_id": "204010dd628ad22cc403b4085d8678a797ab5248",
          "url": "https://github.com/galacticusorg/galacticus/commit/ce515f9e4525c3464161e5183cbd7bfdc56c4e0d"
        },
        "date": 1668806872681,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 58.754,
            "unit": "seconds",
            "range": 0.0249479458058025
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
          "id": "ce515f9e4525c3464161e5183cbd7bfdc56c4e0d",
          "message": "Merge pull request #332 from galacticusorg/componentSubParameters\n\nMake components read their options from their own subparameters",
          "timestamp": "2022-11-18T08:33:45-08:00",
          "tree_id": "204010dd628ad22cc403b4085d8678a797ab5248",
          "url": "https://github.com/galacticusorg/galacticus/commit/ce515f9e4525c3464161e5183cbd7bfdc56c4e0d"
        },
        "date": 1668806886354,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.1223603836366,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.7948185912783,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26209.4879648531,
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
          "id": "380311c81f358602f5bfc449c3a59549db810f16",
          "message": "Merge pull request #335 from galacticusorg/satelliteDistanceMinimum\n\nAdd a `nodeOperator` and `nodePropertyExtractor` to track the minimum distance of approach of a satellite to the center of its host halo",
          "timestamp": "2022-11-28T21:33:33-08:00",
          "tree_id": "7cb9bf93ba71ef9efadac4301d6ea1f0d00b8053",
          "url": "https://github.com/galacticusorg/galacticus/commit/380311c81f358602f5bfc449c3a59549db810f16"
        },
        "date": 1669721034608,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 64.678,
            "unit": "seconds",
            "range": 0.17965411211561
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
          "id": "380311c81f358602f5bfc449c3a59549db810f16",
          "message": "Merge pull request #335 from galacticusorg/satelliteDistanceMinimum\n\nAdd a `nodeOperator` and `nodePropertyExtractor` to track the minimum distance of approach of a satellite to the center of its host halo",
          "timestamp": "2022-11-28T21:33:33-08:00",
          "tree_id": "7cb9bf93ba71ef9efadac4301d6ea1f0d00b8053",
          "url": "https://github.com/galacticusorg/galacticus/commit/380311c81f358602f5bfc449c3a59549db810f16"
        },
        "date": 1669721043431,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.4771891211537,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.0053889680488,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26391.8664916722,
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
          "id": "c4e8e9541952d9b1e4b05d9c3fb2155e0bf4aa99",
          "message": "Merge pull request #336 from galacticusorg/krumholz2009Fix\n\nCatch disks with tiny gas densities",
          "timestamp": "2022-11-30T19:59:12-08:00",
          "tree_id": "68a6fc00f9fd673fd464341b7230db190d006249",
          "url": "https://github.com/galacticusorg/galacticus/commit/c4e8e9541952d9b1e4b05d9c3fb2155e0bf4aa99"
        },
        "date": 1669876299679,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 50.39,
            "unit": "seconds",
            "range": 0.100399203184397
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
          "id": "c4e8e9541952d9b1e4b05d9c3fb2155e0bf4aa99",
          "message": "Merge pull request #336 from galacticusorg/krumholz2009Fix\n\nCatch disks with tiny gas densities",
          "timestamp": "2022-11-30T19:59:12-08:00",
          "tree_id": "68a6fc00f9fd673fd464341b7230db190d006249",
          "url": "https://github.com/galacticusorg/galacticus/commit/c4e8e9541952d9b1e4b05d9c3fb2155e0bf4aa99"
        },
        "date": 1669876308195,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.1223603836366,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.7948185912783,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26209.4879648531,
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
          "id": "7be5401636bedde25e4d0352d00fa7acf5a86c5e",
          "message": "Merge pull request #337 from galacticusorg/constructorResultSelf\n\nUse `self` as the result object in all constructors",
          "timestamp": "2022-12-01T21:37:10-08:00",
          "tree_id": "81b8960e176f65948b522f2ba8650d258663116e",
          "url": "https://github.com/galacticusorg/galacticus/commit/7be5401636bedde25e4d0352d00fa7acf5a86c5e"
        },
        "date": 1669969130298,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 52.35,
            "unit": "seconds",
            "range": 0.0422610932176075
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
          "id": "7be5401636bedde25e4d0352d00fa7acf5a86c5e",
          "message": "Merge pull request #337 from galacticusorg/constructorResultSelf\n\nUse `self` as the result object in all constructors",
          "timestamp": "2022-12-01T21:37:10-08:00",
          "tree_id": "81b8960e176f65948b522f2ba8650d258663116e",
          "url": "https://github.com/galacticusorg/galacticus/commit/7be5401636bedde25e4d0352d00fa7acf5a86c5e"
        },
        "date": 1669969137795,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.1223603836366,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.7948185912783,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26209.4879648531,
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
          "id": "6f3a157d52f9d699060947b92874ecf5cfc7dae8",
          "message": "fix: Use unicode characters in comments where possible",
          "timestamp": "2022-12-02T15:05:58-08:00",
          "tree_id": "29b6fcd195d92414cc44636cf69d97127f95b7ec",
          "url": "https://github.com/galacticusorg/galacticus/commit/6f3a157d52f9d699060947b92874ecf5cfc7dae8"
        },
        "date": 1670032543216,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.535,
            "unit": "seconds",
            "range": 0.028469281689976
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
          "id": "6f3a157d52f9d699060947b92874ecf5cfc7dae8",
          "message": "fix: Use unicode characters in comments where possible",
          "timestamp": "2022-12-02T15:05:58-08:00",
          "tree_id": "29b6fcd195d92414cc44636cf69d97127f95b7ec",
          "url": "https://github.com/galacticusorg/galacticus/commit/6f3a157d52f9d699060947b92874ecf5cfc7dae8"
        },
        "date": 1670032552139,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.1223603836366,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.7948185912783,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26209.4879648531,
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
          "id": "5ea7873f6a3a1ef08f67d114710b77cca7d9529b",
          "message": "Merge pull request #338 from galacticusorg/constrainedTrees\n\nAdd functionality for building constrained merger trees",
          "timestamp": "2022-12-05T20:54:55-08:00",
          "tree_id": "5c8098f5f07ea2fb1ff6f841dfb6033304849fc3",
          "url": "https://github.com/galacticusorg/galacticus/commit/5ea7873f6a3a1ef08f67d114710b77cca7d9529b"
        },
        "date": 1670322839671,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 63.607,
            "unit": "seconds",
            "range": 0.123134479330149
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
          "id": "5ea7873f6a3a1ef08f67d114710b77cca7d9529b",
          "message": "Merge pull request #338 from galacticusorg/constrainedTrees\n\nAdd functionality for building constrained merger trees",
          "timestamp": "2022-12-05T20:54:55-08:00",
          "tree_id": "5c8098f5f07ea2fb1ff6f841dfb6033304849fc3",
          "url": "https://github.com/galacticusorg/galacticus/commit/5ea7873f6a3a1ef08f67d114710b77cca7d9529b"
        },
        "date": 1670322848989,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3752255203618,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.1474686204082,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27887.2181311663,
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
          "id": "a09b0fa26436b205f555778ad89f1ff941544530",
          "message": "Merge pull request #340 from galacticusorg/massDistributions\n\nAdd new mass distributions",
          "timestamp": "2022-12-07T07:37:12-08:00",
          "tree_id": "3781fa2371af076f5713d0183c00ddd4b0f4f83d",
          "url": "https://github.com/galacticusorg/galacticus/commit/a09b0fa26436b205f555778ad89f1ff941544530"
        },
        "date": 1670444812907,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 59.932,
            "unit": "seconds",
            "range": 0.179971108792713
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
          "id": "a09b0fa26436b205f555778ad89f1ff941544530",
          "message": "Merge pull request #340 from galacticusorg/massDistributions\n\nAdd new mass distributions",
          "timestamp": "2022-12-07T07:37:12-08:00",
          "tree_id": "3781fa2371af076f5713d0183c00ddd4b0f4f83d",
          "url": "https://github.com/galacticusorg/galacticus/commit/a09b0fa26436b205f555778ad89f1ff941544530"
        },
        "date": 1670444822652,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3752255203618,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.1474686204082,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27887.2181311663,
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
          "id": "93588f314d5a73a1ae63a5bf908ffe0d7d1acefe",
          "message": "Merge pull request #341 from galacticusorg/h2Clumping\n\nAccount for hot halo density profile in molecular hydrogen calculations",
          "timestamp": "2022-12-08T05:45:42-08:00",
          "tree_id": "9ba658a1a08a5d54d49e1eb5284880b78fb206c3",
          "url": "https://github.com/galacticusorg/galacticus/commit/93588f314d5a73a1ae63a5bf908ffe0d7d1acefe"
        },
        "date": 1670521157234,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 50.527,
            "unit": "seconds",
            "range": 0.0570622467131572
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
          "id": "93588f314d5a73a1ae63a5bf908ffe0d7d1acefe",
          "message": "Merge pull request #341 from galacticusorg/h2Clumping\n\nAccount for hot halo density profile in molecular hydrogen calculations",
          "timestamp": "2022-12-08T05:45:42-08:00",
          "tree_id": "9ba658a1a08a5d54d49e1eb5284880b78fb206c3",
          "url": "https://github.com/galacticusorg/galacticus/commit/93588f314d5a73a1ae63a5bf908ffe0d7d1acefe"
        },
        "date": 1670521165691,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3752255203618,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.1474686204082,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27887.2181311663,
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
          "id": "c0ceae113abac02e4eab71eb137360c49308d11e",
          "message": "Merge pull request #343 from galacticusorg/massDistributions\n\nAdd support for component and mass types in `massDistribution` classes",
          "timestamp": "2022-12-12T17:02:30-08:00",
          "tree_id": "719b329ab361427c0147b1960d842a23ff63d174",
          "url": "https://github.com/galacticusorg/galacticus/commit/c0ceae113abac02e4eab71eb137360c49308d11e"
        },
        "date": 1670903236616,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 53.885,
            "unit": "seconds",
            "range": 0.371906574289683
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
          "id": "c0ceae113abac02e4eab71eb137360c49308d11e",
          "message": "Merge pull request #343 from galacticusorg/massDistributions\n\nAdd support for component and mass types in `massDistribution` classes",
          "timestamp": "2022-12-12T17:02:30-08:00",
          "tree_id": "719b329ab361427c0147b1960d842a23ff63d174",
          "url": "https://github.com/galacticusorg/galacticus/commit/c0ceae113abac02e4eab71eb137360c49308d11e"
        },
        "date": 1670903244989,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.4588651544499,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.3170227176869,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27856.5941320004,
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
          "id": "fb2e85397c804ff72669d750aedeb5db934a848f",
          "message": "Merge pull request #345 from galacticusorg/baryonicModelDebug\n\nDebug various problems related to the constrained baryonic physics model",
          "timestamp": "2022-12-13T17:24:58-08:00",
          "tree_id": "6c946e0bfbb3622509dd0b60570e208bb148ecce",
          "url": "https://github.com/galacticusorg/galacticus/commit/fb2e85397c804ff72669d750aedeb5db934a848f"
        },
        "date": 1671029388482,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 59.97,
            "unit": "seconds",
            "range": 0.102625532885841
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
          "id": "fb2e85397c804ff72669d750aedeb5db934a848f",
          "message": "Merge pull request #345 from galacticusorg/baryonicModelDebug\n\nDebug various problems related to the constrained baryonic physics model",
          "timestamp": "2022-12-13T17:24:58-08:00",
          "tree_id": "6c946e0bfbb3622509dd0b60570e208bb148ecce",
          "url": "https://github.com/galacticusorg/galacticus/commit/fb2e85397c804ff72669d750aedeb5db934a848f"
        },
        "date": 1671029397435,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3752255203618,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.1474686204082,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27887.2181311663,
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
          "id": "a9f8f5b8fc960e2b52a6075fa7adfe9fbd5aa5bb",
          "message": "Merge pull request #346 from galacticusorg/igbrHDF5\n\nSupport HDF5 format for intergalactic background radiation data files",
          "timestamp": "2022-12-17T05:44:03-08:00",
          "tree_id": "04872b9274be1b8317f3a4ea85c66bb50c64500b",
          "url": "https://github.com/galacticusorg/galacticus/commit/a9f8f5b8fc960e2b52a6075fa7adfe9fbd5aa5bb"
        },
        "date": 1671294881056,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.898,
            "unit": "seconds",
            "range": 0.0714114836699856
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
          "id": "a9f8f5b8fc960e2b52a6075fa7adfe9fbd5aa5bb",
          "message": "Merge pull request #346 from galacticusorg/igbrHDF5\n\nSupport HDF5 format for intergalactic background radiation data files",
          "timestamp": "2022-12-17T05:44:03-08:00",
          "tree_id": "04872b9274be1b8317f3a4ea85c66bb50c64500b",
          "url": "https://github.com/galacticusorg/galacticus/commit/a9f8f5b8fc960e2b52a6075fa7adfe9fbd5aa5bb"
        },
        "date": 1671294888270,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3752255203618,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.1474686204082,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27887.2181311663,
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
          "id": "429de06c773f660e4d551bfdf44adcff5407204e",
          "message": "fix: Deallocate `elementsToTrack` if necessary before reallocating",
          "timestamp": "2022-12-20T17:18:46Z",
          "tree_id": "3f54c2f22764414cab1ef835b54f1e62704cf345",
          "url": "https://github.com/galacticusorg/galacticus/commit/429de06c773f660e4d551bfdf44adcff5407204e"
        },
        "date": 1671565944333,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 52.151,
            "unit": "seconds",
            "range": 0.0413388437191904
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
          "id": "429de06c773f660e4d551bfdf44adcff5407204e",
          "message": "fix: Deallocate `elementsToTrack` if necessary before reallocating",
          "timestamp": "2022-12-20T17:18:46Z",
          "tree_id": "3f54c2f22764414cab1ef835b54f1e62704cf345",
          "url": "https://github.com/galacticusorg/galacticus/commit/429de06c773f660e4d551bfdf44adcff5407204e"
        },
        "date": 1671565955723,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3752255203618,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.1474686204082,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27887.2181311663,
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
          "id": "41073b970bbcceef060e5d689f7fc81502c552d3",
          "message": "fix: `mpiCounter` destructor must reset initialization state\n\nOtherwise such an object can be finalizef twice causing a segfault.",
          "timestamp": "2022-12-22T17:54:28Z",
          "tree_id": "dd72f9c5bdeda938043f148e7570c946bc8b6811",
          "url": "https://github.com/galacticusorg/galacticus/commit/41073b970bbcceef060e5d689f7fc81502c552d3"
        },
        "date": 1671742146629,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 52.371,
            "unit": "seconds",
            "range": 0.0567000881838577
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
          "id": "41073b970bbcceef060e5d689f7fc81502c552d3",
          "message": "fix: `mpiCounter` destructor must reset initialization state\n\nOtherwise such an object can be finalizef twice causing a segfault.",
          "timestamp": "2022-12-22T17:54:28Z",
          "tree_id": "dd72f9c5bdeda938043f148e7570c946bc8b6811",
          "url": "https://github.com/galacticusorg/galacticus/commit/41073b970bbcceef060e5d689f7fc81502c552d3"
        },
        "date": 1671742155213,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3635405347359,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.0294603368507,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26684.9112442229,
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
          "id": "eb8a10b2903ae1458484d5b955a0a583d6e005b5",
          "message": "fix: Add missing optional arguments",
          "timestamp": "2022-12-25T18:27:06Z",
          "tree_id": "a277d937a12b62a526e262fe81e510550914a55f",
          "url": "https://github.com/galacticusorg/galacticus/commit/eb8a10b2903ae1458484d5b955a0a583d6e005b5"
        },
        "date": 1672009607865,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 52.296,
            "unit": "seconds",
            "range": 0.0257371327058633
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
          "id": "eb8a10b2903ae1458484d5b955a0a583d6e005b5",
          "message": "fix: Add missing optional arguments",
          "timestamp": "2022-12-25T18:27:06Z",
          "tree_id": "a277d937a12b62a526e262fe81e510550914a55f",
          "url": "https://github.com/galacticusorg/galacticus/commit/eb8a10b2903ae1458484d5b955a0a583d6e005b5"
        },
        "date": 1672009617300,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3752255203618,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.1474686204082,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27887.2181311663,
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
          "id": "a655520045c0c9c418fbc040a4d395937907b203",
          "message": "Merge pull request #347 from galacticusorg/mpiCounterReset\n\nAdd a reset method to MPI counters",
          "timestamp": "2022-12-26T08:53:03-07:00",
          "tree_id": "82465daf4378f4c8dc58ffc31cf0ec8d161b6f57",
          "url": "https://github.com/galacticusorg/galacticus/commit/a655520045c0c9c418fbc040a4d395937907b203"
        },
        "date": 1672080008501,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.033,
            "unit": "seconds",
            "range": 0.0292933439539586
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
          "id": "a655520045c0c9c418fbc040a4d395937907b203",
          "message": "Merge pull request #347 from galacticusorg/mpiCounterReset\n\nAdd a reset method to MPI counters",
          "timestamp": "2022-12-26T08:53:03-07:00",
          "tree_id": "82465daf4378f4c8dc58ffc31cf0ec8d161b6f57",
          "url": "https://github.com/galacticusorg/galacticus/commit/a655520045c0c9c418fbc040a4d395937907b203"
        },
        "date": 1672080017354,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3635405347359,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.0294603368507,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26684.9112442229,
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
          "id": "5ff1b5f390e2744ba1a00c869bde318f1ba4973c",
          "message": "fix: Update copyright year",
          "timestamp": "2023-01-09T10:43:32-08:00",
          "tree_id": "1d1886ad70a8dc36a51dd065f274185290aece54",
          "url": "https://github.com/galacticusorg/galacticus/commit/5ff1b5f390e2744ba1a00c869bde318f1ba4973c"
        },
        "date": 1673314909979,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 60.785,
            "unit": "seconds",
            "range": 0.0386069941844015
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
          "id": "5ff1b5f390e2744ba1a00c869bde318f1ba4973c",
          "message": "fix: Update copyright year",
          "timestamp": "2023-01-09T10:43:32-08:00",
          "tree_id": "1d1886ad70a8dc36a51dd065f274185290aece54",
          "url": "https://github.com/galacticusorg/galacticus/commit/5ff1b5f390e2744ba1a00c869bde318f1ba4973c"
        },
        "date": 1673314918306,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3752255203618,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.1474686204082,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27887.2181311663,
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
          "id": "618eb5b594676b31af0db361aef52aa77085632c",
          "message": "Merge pull request #348 from galacticusorg/nullInitializePointers\n\nEnsure all pointers in derived-type objects all null initialized",
          "timestamp": "2023-01-10T15:55:55-08:00",
          "tree_id": "0fd46a2948334ecdb1ece10b2c13616a3b7df721",
          "url": "https://github.com/galacticusorg/galacticus/commit/618eb5b594676b31af0db361aef52aa77085632c"
        },
        "date": 1673412835822,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.348,
            "unit": "seconds",
            "range": 0.0414921679358762
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
          "id": "618eb5b594676b31af0db361aef52aa77085632c",
          "message": "Merge pull request #348 from galacticusorg/nullInitializePointers\n\nEnsure all pointers in derived-type objects all null initialized",
          "timestamp": "2023-01-10T15:55:55-08:00",
          "tree_id": "0fd46a2948334ecdb1ece10b2c13616a3b7df721",
          "url": "https://github.com/galacticusorg/galacticus/commit/618eb5b594676b31af0db361aef52aa77085632c"
        },
        "date": 1673412844400,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3752255203618,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.1474686204082,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27887.2181311663,
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
          "id": "1fd62974f45f3a114efd9dba474113e7e151c5a0",
          "message": "Merge pull request #349 from galacticusorg/orbitingNBodyTrees\n\nReset satellite timescales when a subhalo is promoted to an isolated halo",
          "timestamp": "2023-01-11T06:10:38-08:00",
          "tree_id": "db4a3c4fc16467627d563c488dd475656a3cc7ab",
          "url": "https://github.com/galacticusorg/galacticus/commit/1fd62974f45f3a114efd9dba474113e7e151c5a0"
        },
        "date": 1673461423403,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 60.165,
            "unit": "seconds",
            "range": 0.171185571821762
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
          "id": "1fd62974f45f3a114efd9dba474113e7e151c5a0",
          "message": "Merge pull request #349 from galacticusorg/orbitingNBodyTrees\n\nReset satellite timescales when a subhalo is promoted to an isolated halo",
          "timestamp": "2023-01-11T06:10:38-08:00",
          "tree_id": "db4a3c4fc16467627d563c488dd475656a3cc7ab",
          "url": "https://github.com/galacticusorg/galacticus/commit/1fd62974f45f3a114efd9dba474113e7e151c5a0"
        },
        "date": 1673461432514,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3752255203618,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.1474686204082,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27887.2181311663,
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
          "id": "ebb98121bba118e49fac2cfec2d68dea4fbedba2",
          "message": "fix: Use copy constructors when walking the parameter tree\n\nThis avoids segmentation faults due to finalization on assignment.",
          "timestamp": "2023-01-12T11:22:51-08:00",
          "tree_id": "9fa7a16930cefaac01d0a539162b84ddc6123b1f",
          "url": "https://github.com/galacticusorg/galacticus/commit/ebb98121bba118e49fac2cfec2d68dea4fbedba2"
        },
        "date": 1673575181189,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 59.602,
            "unit": "seconds",
            "range": 0.0799099493172472
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
          "id": "ebb98121bba118e49fac2cfec2d68dea4fbedba2",
          "message": "fix: Use copy constructors when walking the parameter tree\n\nThis avoids segmentation faults due to finalization on assignment.",
          "timestamp": "2023-01-12T11:22:51-08:00",
          "tree_id": "9fa7a16930cefaac01d0a539162b84ddc6123b1f",
          "url": "https://github.com/galacticusorg/galacticus/commit/ebb98121bba118e49fac2cfec2d68dea4fbedba2"
        },
        "date": 1673575191058,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3752255203618,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.1474686204082,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27887.2181311663,
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
          "id": "edf56bd05e127325a2415787df342c4e72f15671",
          "message": "feat: Allow for smoothed transitions in high-pass filters\n\nMake use of these in the COSMOS SHMR likelihood class to avoid having sharp transitions in the likelihood function.",
          "timestamp": "2023-01-13T17:13:09Z",
          "tree_id": "b1c263cc0161db07919d832a57f94b36c905bbc8",
          "url": "https://github.com/galacticusorg/galacticus/commit/edf56bd05e127325a2415787df342c4e72f15671"
        },
        "date": 1673651805661,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 63.681,
            "unit": "seconds",
            "range": 0.315732323337358
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
          "id": "edf56bd05e127325a2415787df342c4e72f15671",
          "message": "feat: Allow for smoothed transitions in high-pass filters\n\nMake use of these in the COSMOS SHMR likelihood class to avoid having sharp transitions in the likelihood function.",
          "timestamp": "2023-01-13T17:13:09Z",
          "tree_id": "b1c263cc0161db07919d832a57f94b36c905bbc8",
          "url": "https://github.com/galacticusorg/galacticus/commit/edf56bd05e127325a2415787df342c4e72f15671"
        },
        "date": 1673651814888,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3752255203618,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.1474686204082,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27887.2181311663,
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
          "id": "baab79ac9a0bba2e7f9c4e4c4cb864ca535bbd41",
          "message": "Merge pull request #351 from galacticusorg/sigmaMIntegration\n\nMake integration of σ(M) more robust",
          "timestamp": "2023-01-13T17:02:28-08:00",
          "tree_id": "359f77eaa4f9f8f97befa8164933b1db07d74300",
          "url": "https://github.com/galacticusorg/galacticus/commit/baab79ac9a0bba2e7f9c4e4c4cb864ca535bbd41"
        },
        "date": 1673668324675,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.322,
            "unit": "seconds",
            "range": 0.0219909072109672
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
          "id": "baab79ac9a0bba2e7f9c4e4c4cb864ca535bbd41",
          "message": "Merge pull request #351 from galacticusorg/sigmaMIntegration\n\nMake integration of σ(M) more robust",
          "timestamp": "2023-01-13T17:02:28-08:00",
          "tree_id": "359f77eaa4f9f8f97befa8164933b1db07d74300",
          "url": "https://github.com/galacticusorg/galacticus/commit/baab79ac9a0bba2e7f9c4e4c4cb864ca535bbd41"
        },
        "date": 1673668333121,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3844543134178,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9186267741694,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26366.2430061311,
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
          "id": "9adeb9a47c16ec662d76136db78461fd1526e16c",
          "message": "Merge pull request #352 from galacticusorg/trapZeroOuterRadius\n\nTrap zero hot halo outer radius",
          "timestamp": "2023-01-15T14:02:27-08:00",
          "tree_id": "5092d2fc43c7741617c69dce444b7af551a39a77",
          "url": "https://github.com/galacticusorg/galacticus/commit/9adeb9a47c16ec662d76136db78461fd1526e16c"
        },
        "date": 1673829653141,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 50.184,
            "unit": "seconds",
            "range": 0.0399799949970758
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
          "id": "9adeb9a47c16ec662d76136db78461fd1526e16c",
          "message": "Merge pull request #352 from galacticusorg/trapZeroOuterRadius\n\nTrap zero hot halo outer radius",
          "timestamp": "2023-01-15T14:02:27-08:00",
          "tree_id": "5092d2fc43c7741617c69dce444b7af551a39a77",
          "url": "https://github.com/galacticusorg/galacticus/commit/9adeb9a47c16ec662d76136db78461fd1526e16c"
        },
        "date": 1673829663631,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.2586127359834,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.0610118696272,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26810.9028287877,
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
          "id": "7c8da1a05472814ae8eae84412d2d7a147abf5f7",
          "message": "feat: Add a function to return GSL error details\n\nUseful for when GSL errors are trapped and we want to provide context in the error message.",
          "timestamp": "2023-01-17T23:40:10Z",
          "tree_id": "4ec9f1edae72fe248d7247152788eccaab7519f5",
          "url": "https://github.com/galacticusorg/galacticus/commit/7c8da1a05472814ae8eae84412d2d7a147abf5f7"
        },
        "date": 1674044456797,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 60.085,
            "unit": "seconds",
            "range": 0.123622408972253
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
          "id": "7c8da1a05472814ae8eae84412d2d7a147abf5f7",
          "message": "feat: Add a function to return GSL error details\n\nUseful for when GSL errors are trapped and we want to provide context in the error message.",
          "timestamp": "2023-01-17T23:40:10Z",
          "tree_id": "4ec9f1edae72fe248d7247152788eccaab7519f5",
          "url": "https://github.com/galacticusorg/galacticus/commit/7c8da1a05472814ae8eae84412d2d7a147abf5f7"
        },
        "date": 1674044465898,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3844543134178,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9186267741694,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26366.2430061311,
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
          "id": "e2cb6f3570438aa22bbdd4d6bf4ed18668d99434",
          "message": "Merge pull request #356 from galacticusorg/font2008RoundingErrors\n\nfix: Avoid root search range failure in `font2008` hot halo ram pressure stripping class",
          "timestamp": "2023-01-18T07:42:56-08:00",
          "tree_id": "5c2c250b769545c02df451188ebbe2436c09018c",
          "url": "https://github.com/galacticusorg/galacticus/commit/e2cb6f3570438aa22bbdd4d6bf4ed18668d99434"
        },
        "date": 1674070322350,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 61.492,
            "unit": "seconds",
            "range": 0.29433586257872
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
          "id": "e2cb6f3570438aa22bbdd4d6bf4ed18668d99434",
          "message": "Merge pull request #356 from galacticusorg/font2008RoundingErrors\n\nfix: Avoid root search range failure in `font2008` hot halo ram pressure stripping class",
          "timestamp": "2023-01-18T07:42:56-08:00",
          "tree_id": "5c2c250b769545c02df451188ebbe2436c09018c",
          "url": "https://github.com/galacticusorg/galacticus/commit/e2cb6f3570438aa22bbdd4d6bf4ed18668d99434"
        },
        "date": 1674070332616,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3844543134178,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9186267741694,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26366.2430061311,
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
          "id": "6b1d1a78f5991eaf1c8b3906fea05bee8179c2e1",
          "message": "Merge pull request #357 from galacticusorg/excursionSetsProgressBarFix\n\nFix an issue that the progress bar is inaccurate when solving the excursion sets problem with MPI",
          "timestamp": "2023-01-18T21:09:36-08:00",
          "tree_id": "5c02a2059732b8f299be3b725a8f536a1df6c8a8",
          "url": "https://github.com/galacticusorg/galacticus/commit/6b1d1a78f5991eaf1c8b3906fea05bee8179c2e1"
        },
        "date": 1674114130063,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.939,
            "unit": "seconds",
            "range": 0.089481282958862
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
          "id": "6b1d1a78f5991eaf1c8b3906fea05bee8179c2e1",
          "message": "Merge pull request #357 from galacticusorg/excursionSetsProgressBarFix\n\nFix an issue that the progress bar is inaccurate when solving the excursion sets problem with MPI",
          "timestamp": "2023-01-18T21:09:36-08:00",
          "tree_id": "5c02a2059732b8f299be3b725a8f536a1df6c8a8",
          "url": "https://github.com/galacticusorg/galacticus/commit/6b1d1a78f5991eaf1c8b3906fea05bee8179c2e1"
        },
        "date": 1674114139509,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3844543134178,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9186267741694,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26366.2430061311,
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
          "id": "c17eea79b5ad78e21a97d2ab909b4a5b26a0af3a",
          "message": "feat: Refactor the Gaussian emulator model likelihood class\n\n* Makes this class more efficient;\n\n* Handles cases where the simulator has intrinsic variance (i.e. exploits the fact that the emulator can not possibly have smaller variance than the simulator);\n\n* Adds a variogram class to provide more flexibility in variogram models.",
          "timestamp": "2023-01-20T19:41:51Z",
          "tree_id": "12e20e3efa8f3c824d5c0462c7bd52d7cd17e12b",
          "url": "https://github.com/galacticusorg/galacticus/commit/c17eea79b5ad78e21a97d2ab909b4a5b26a0af3a"
        },
        "date": 1674268920224,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 59.694,
            "unit": "seconds",
            "range": 0.065881712181075
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
          "id": "c17eea79b5ad78e21a97d2ab909b4a5b26a0af3a",
          "message": "feat: Refactor the Gaussian emulator model likelihood class\n\n* Makes this class more efficient;\n\n* Handles cases where the simulator has intrinsic variance (i.e. exploits the fact that the emulator can not possibly have smaller variance than the simulator);\n\n* Adds a variogram class to provide more flexibility in variogram models.",
          "timestamp": "2023-01-20T19:41:51Z",
          "tree_id": "12e20e3efa8f3c824d5c0462c7bd52d7cd17e12b",
          "url": "https://github.com/galacticusorg/galacticus/commit/c17eea79b5ad78e21a97d2ab909b4a5b26a0af3a"
        },
        "date": 1674268928087,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3844543134178,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9186267741694,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26366.2430061311,
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
          "id": "2a133f61cef37420abaf49b2a4e9b275ffbbccad",
          "message": "Merge pull request #359 from galacticusorg/inputParameterFixes\n\nEnsure FoX accesses are inside appropriate critical sections",
          "timestamp": "2023-01-20T22:44:28-08:00",
          "tree_id": "55f6c619a26ee490ca918d5efa53b585f42254ad",
          "url": "https://github.com/galacticusorg/galacticus/commit/2a133f61cef37420abaf49b2a4e9b275ffbbccad"
        },
        "date": 1674294080545,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.998,
            "unit": "seconds",
            "range": 0.0233152310711439
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
          "id": "2a133f61cef37420abaf49b2a4e9b275ffbbccad",
          "message": "Merge pull request #359 from galacticusorg/inputParameterFixes\n\nEnsure FoX accesses are inside appropriate critical sections",
          "timestamp": "2023-01-20T22:44:28-08:00",
          "tree_id": "55f6c619a26ee490ca918d5efa53b585f42254ad",
          "url": "https://github.com/galacticusorg/galacticus/commit/2a133f61cef37420abaf49b2a4e9b275ffbbccad"
        },
        "date": 1674294087777,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3844543134178,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9186267741694,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26366.2430061311,
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
          "id": "b9e1110209090aebbab11164a64a25d45b51945e",
          "message": "Merge pull request #358 from galacticusorg/namespaceVariables\n\nRemove prefix-based namespacing of submodule-scope variables",
          "timestamp": "2023-01-21T08:32:54-08:00",
          "tree_id": "83d720b7565a9f25568c5e060562f0c9546ad532",
          "url": "https://github.com/galacticusorg/galacticus/commit/b9e1110209090aebbab11164a64a25d45b51945e"
        },
        "date": 1674328006820,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.97,
            "unit": "seconds",
            "range": 0.0379736750932567
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
          "id": "b9e1110209090aebbab11164a64a25d45b51945e",
          "message": "Merge pull request #358 from galacticusorg/namespaceVariables\n\nRemove prefix-based namespacing of submodule-scope variables",
          "timestamp": "2023-01-21T08:32:54-08:00",
          "tree_id": "83d720b7565a9f25568c5e060562f0c9546ad532",
          "url": "https://github.com/galacticusorg/galacticus/commit/b9e1110209090aebbab11164a64a25d45b51945e"
        },
        "date": 1674328014992,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3149518526459,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.848186606352,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27713.5248462996,
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
          "id": "5664c7f70e5ae715ead9b390257b8dcb5e68ab4a",
          "message": "fix: Trap out of range values and floating points errors in the `bhattacharya2011` halo mass function class",
          "timestamp": "2023-01-23T09:06:27-08:00",
          "tree_id": "50b7f975db924da05b2ae0ea8583d2a312c178f8",
          "url": "https://github.com/galacticusorg/galacticus/commit/5664c7f70e5ae715ead9b390257b8dcb5e68ab4a"
        },
        "date": 1674506634443,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 64.296,
            "unit": "seconds",
            "range": 0.184278050781548
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
          "id": "5664c7f70e5ae715ead9b390257b8dcb5e68ab4a",
          "message": "fix: Trap out of range values and floating points errors in the `bhattacharya2011` halo mass function class",
          "timestamp": "2023-01-23T09:06:27-08:00",
          "tree_id": "50b7f975db924da05b2ae0ea8583d2a312c178f8",
          "url": "https://github.com/galacticusorg/galacticus/commit/5664c7f70e5ae715ead9b390257b8dcb5e68ab4a"
        },
        "date": 1674506643479,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3844543134178,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9186267741694,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26366.2430061311,
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
          "id": "ae48d60d531a85eacf32f66f2a7930be70f1960f",
          "message": "Merge pull request #360 from galacticusorg/interpolator2D\n\nAdd a simple 2D interpolator class",
          "timestamp": "2023-01-26T06:19:50-08:00",
          "tree_id": "e8501ca80058ef1fa74c8dc9a5704cbc3d6e0157",
          "url": "https://github.com/galacticusorg/galacticus/commit/ae48d60d531a85eacf32f66f2a7930be70f1960f"
        },
        "date": 1674764753106,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 59.895,
            "unit": "seconds",
            "range": 0.219873827455829
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
          "id": "ae48d60d531a85eacf32f66f2a7930be70f1960f",
          "message": "Merge pull request #360 from galacticusorg/interpolator2D\n\nAdd a simple 2D interpolator class",
          "timestamp": "2023-01-26T06:19:50-08:00",
          "tree_id": "e8501ca80058ef1fa74c8dc9a5704cbc3d6e0157",
          "url": "https://github.com/galacticusorg/galacticus/commit/ae48d60d531a85eacf32f66f2a7930be70f1960f"
        },
        "date": 1674764760886,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3844543134178,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9186267741694,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26366.2430061311,
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
          "id": "f6d60ffbc3145e9bb417842402d4b1262efbbdb2",
          "message": "fix: Fix the `mergerTreeBuildMassesUnion` class\n\n* Correctly handle deep copy and state store of member classes;\n\n* Mark that multiple `mergerTreeBuildMasses` members are allowed;\n\n* Fix typo which caused double allocation of one of the masses arrays.",
          "timestamp": "2023-01-27T00:50:53Z",
          "tree_id": "5eab63d7834b71dea95ddf5b49e823f0dbca1c07",
          "url": "https://github.com/galacticusorg/galacticus/commit/f6d60ffbc3145e9bb417842402d4b1262efbbdb2"
        },
        "date": 1674799574195,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.901,
            "unit": "seconds",
            "range": 0.0214685816921121
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
          "id": "f6d60ffbc3145e9bb417842402d4b1262efbbdb2",
          "message": "fix: Fix the `mergerTreeBuildMassesUnion` class\n\n* Correctly handle deep copy and state store of member classes;\n\n* Mark that multiple `mergerTreeBuildMasses` members are allowed;\n\n* Fix typo which caused double allocation of one of the masses arrays.",
          "timestamp": "2023-01-27T00:50:53Z",
          "tree_id": "5eab63d7834b71dea95ddf5b49e823f0dbca1c07",
          "url": "https://github.com/galacticusorg/galacticus/commit/f6d60ffbc3145e9bb417842402d4b1262efbbdb2"
        },
        "date": 1674799583357,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.2577034798009,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.1512854931292,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26350.9375297304,
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
          "id": "7e0956d088d0d34abfcaadec5ca77daa33d99768",
          "message": "Merge pull request #362 from galacticusorg/satelliteDestructThreshold\n\nAllow satellite destruction mass threshold to additionally scale with the tree mass",
          "timestamp": "2023-01-27T06:20:50-08:00",
          "tree_id": "f30ad2ecc7bd2ec643af54568478558189cc348c",
          "url": "https://github.com/galacticusorg/galacticus/commit/7e0956d088d0d34abfcaadec5ca77daa33d99768"
        },
        "date": 1674847878455,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.561,
            "unit": "seconds",
            "range": 0.058265770396311
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
          "id": "7e0956d088d0d34abfcaadec5ca77daa33d99768",
          "message": "Merge pull request #362 from galacticusorg/satelliteDestructThreshold\n\nAllow satellite destruction mass threshold to additionally scale with the tree mass",
          "timestamp": "2023-01-27T06:20:50-08:00",
          "tree_id": "f30ad2ecc7bd2ec643af54568478558189cc348c",
          "url": "https://github.com/galacticusorg/galacticus/commit/7e0956d088d0d34abfcaadec5ca77daa33d99768"
        },
        "date": 1674847886187,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.5967294919036,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.4513198173429,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27682.9945390393,
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
          "id": "79b7325d1563d2760811f8624ec6fbb1a4fb6d91",
          "message": "Merge pull request #363 from galacticusorg/chandraDensity\n\nAvoid computing Chandrasekhar integrals if density is zero",
          "timestamp": "2023-01-27T14:24:33-08:00",
          "tree_id": "f8ba0b4169808a697f7051550c1ffa7c00304fad",
          "url": "https://github.com/galacticusorg/galacticus/commit/79b7325d1563d2760811f8624ec6fbb1a4fb6d91"
        },
        "date": 1674882104785,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 55.861,
            "unit": "seconds",
            "range": 0.307936519432103
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
          "id": "79b7325d1563d2760811f8624ec6fbb1a4fb6d91",
          "message": "Merge pull request #363 from galacticusorg/chandraDensity\n\nAvoid computing Chandrasekhar integrals if density is zero",
          "timestamp": "2023-01-27T14:24:33-08:00",
          "tree_id": "f8ba0b4169808a697f7051550c1ffa7c00304fad",
          "url": "https://github.com/galacticusorg/galacticus/commit/79b7325d1563d2760811f8624ec6fbb1a4fb6d91"
        },
        "date": 1674882113163,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.5967294919036,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.4513198173429,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27682.9945390393,
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
          "id": "a8a286ec1628d80447bafa24fa50974574b3fbf7",
          "message": "Merge pull request #365 from galacticusorg/takahashi2011FileLock\n\nAdd file locking to the `gravitationalLensingTakahashi2011` class to avoid conflicts between processes",
          "timestamp": "2023-01-27T22:20:25-08:00",
          "tree_id": "add031ae1efa8288047b8bb200decc3f61dca0fe",
          "url": "https://github.com/galacticusorg/galacticus/commit/a8a286ec1628d80447bafa24fa50974574b3fbf7"
        },
        "date": 1674897001675,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.658,
            "unit": "seconds",
            "range": 0.0354908438890592
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
          "id": "a8a286ec1628d80447bafa24fa50974574b3fbf7",
          "message": "Merge pull request #365 from galacticusorg/takahashi2011FileLock\n\nAdd file locking to the `gravitationalLensingTakahashi2011` class to avoid conflicts between processes",
          "timestamp": "2023-01-27T22:20:25-08:00",
          "tree_id": "add031ae1efa8288047b8bb200decc3f61dca0fe",
          "url": "https://github.com/galacticusorg/galacticus/commit/a8a286ec1628d80447bafa24fa50974574b3fbf7"
        },
        "date": 1674897010008,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.4721000776937,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9954078427452,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 29402.9078379971,
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
          "id": "7a094cbaef51d934685b780dbe5a1f2db830e1c6",
          "message": "fix: Remove an unnecessary module `use`",
          "timestamp": "2023-01-31T19:34:04Z",
          "tree_id": "435df98dda03f985cd4f627db7a12018a7edf773",
          "url": "https://github.com/galacticusorg/galacticus/commit/7a094cbaef51d934685b780dbe5a1f2db830e1c6"
        },
        "date": 1675210334462,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 50.33,
            "unit": "seconds",
            "range": 0.0749533188051942
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
          "id": "7a094cbaef51d934685b780dbe5a1f2db830e1c6",
          "message": "fix: Remove an unnecessary module `use`",
          "timestamp": "2023-01-31T19:34:04Z",
          "tree_id": "435df98dda03f985cd4f627db7a12018a7edf773",
          "url": "https://github.com/galacticusorg/galacticus/commit/7a094cbaef51d934685b780dbe5a1f2db830e1c6"
        },
        "date": 1675210342175,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.5967294919036,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.4513198173429,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 27682.9945390393,
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
          "id": "558eac414f363eda8870ebf6922b8b21386aace6",
          "message": "fix: Correct parameter names for subhalo destruction mass thresholds",
          "timestamp": "2023-02-02T00:02:39Z",
          "tree_id": "245e18a0dbce1212e9525f1b9913c0e8d5c53dc1",
          "url": "https://github.com/galacticusorg/galacticus/commit/558eac414f363eda8870ebf6922b8b21386aace6"
        },
        "date": 1675306352766,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.337,
            "unit": "seconds",
            "range": 0.0221833270707651
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
          "id": "558eac414f363eda8870ebf6922b8b21386aace6",
          "message": "fix: Correct parameter names for subhalo destruction mass thresholds",
          "timestamp": "2023-02-02T00:02:39Z",
          "tree_id": "245e18a0dbce1212e9525f1b9913c0e8d5c53dc1",
          "url": "https://github.com/galacticusorg/galacticus/commit/558eac414f363eda8870ebf6922b8b21386aace6"
        },
        "date": 1675306361482,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3844543134178,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9186267741694,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26366.2430061311,
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
          "id": "a6b432a7e58272dd1d5a781ce8ff9fa440057d6f",
          "message": "Merge pull request #371 from galacticusorg/oscillatingTk\n\nAllow σ(M) integral to be split at local minima of the transfer function",
          "timestamp": "2023-02-07T18:53:32-08:00",
          "tree_id": "baf78ba1ea22233d882d7d78982086b8466b10eb",
          "url": "https://github.com/galacticusorg/galacticus/commit/a6b432a7e58272dd1d5a781ce8ff9fa440057d6f"
        },
        "date": 1675834176953,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.701,
            "unit": "seconds",
            "range": 0.0333301665162854
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
          "id": "a6b432a7e58272dd1d5a781ce8ff9fa440057d6f",
          "message": "Merge pull request #371 from galacticusorg/oscillatingTk\n\nAllow σ(M) integral to be split at local minima of the transfer function",
          "timestamp": "2023-02-07T18:53:32-08:00",
          "tree_id": "baf78ba1ea22233d882d7d78982086b8466b10eb",
          "url": "https://github.com/galacticusorg/galacticus/commit/a6b432a7e58272dd1d5a781ce8ff9fa440057d6f"
        },
        "date": 1675834184886,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3844543134178,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9186267741694,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 26366.2430061311,
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
          "id": "c24b46c497e12677d6cd492d586cf92757badaaa",
          "message": "Merge pull request #373 from galacticusorg/massDefinitions\n\nAdd functionality to convert between halo mass definitions based on infall time (instead of just current time)",
          "timestamp": "2023-02-09T05:18:20Z",
          "tree_id": "ee1a4954a131b6887721024ae064af7b34b8d0b8",
          "url": "https://github.com/galacticusorg/galacticus/commit/c24b46c497e12677d6cd492d586cf92757badaaa"
        },
        "date": 1675929124063,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 59.683,
            "unit": "seconds",
            "range": 0.137673890044293
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
          "id": "c24b46c497e12677d6cd492d586cf92757badaaa",
          "message": "Merge pull request #373 from galacticusorg/massDefinitions\n\nAdd functionality to convert between halo mass definitions based on infall time (instead of just current time)",
          "timestamp": "2023-02-09T05:18:20Z",
          "tree_id": "ee1a4954a131b6887721024ae064af7b34b8d0b8",
          "url": "https://github.com/galacticusorg/galacticus/commit/c24b46c497e12677d6cd492d586cf92757badaaa"
        },
        "date": 1675929133222,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3717638717242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.689755316789,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25371.9758848307,
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
          "id": "9fb8d320f10a42e07881040c4fc65536c19ed54f",
          "message": "fix: Update CAMB version in Dockerfile",
          "timestamp": "2023-02-09T15:20:32Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/9fb8d320f10a42e07881040c4fc65536c19ed54f"
        },
        "date": 1675965851212,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 50.603,
            "unit": "seconds",
            "range": 0.0215429802962522
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
          "id": "9fb8d320f10a42e07881040c4fc65536c19ed54f",
          "message": "fix: Update CAMB version in Dockerfile",
          "timestamp": "2023-02-09T15:20:32Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/9fb8d320f10a42e07881040c4fc65536c19ed54f"
        },
        "date": 1675965861962,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "f3540a4c84d6aac85141e8b18ece4c201d2bd1ba",
          "message": "Merge pull request #374 from galacticusorg/axionTk\n\nAdd the Passaglia & Hu (2022) transfer function for axionic dark matter",
          "timestamp": "2023-02-10T04:21:22Z",
          "tree_id": "35d77f35302e67f37bda2a2b429e033047f17c09",
          "url": "https://github.com/galacticusorg/galacticus/commit/f3540a4c84d6aac85141e8b18ece4c201d2bd1ba"
        },
        "date": 1676016276861,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 52.117,
            "unit": "seconds",
            "range": 0.0347577329516728
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
          "id": "f3540a4c84d6aac85141e8b18ece4c201d2bd1ba",
          "message": "Merge pull request #374 from galacticusorg/axionTk\n\nAdd the Passaglia & Hu (2022) transfer function for axionic dark matter",
          "timestamp": "2023-02-10T04:21:22Z",
          "tree_id": "35d77f35302e67f37bda2a2b429e033047f17c09",
          "url": "https://github.com/galacticusorg/galacticus/commit/f3540a4c84d6aac85141e8b18ece4c201d2bd1ba"
        },
        "date": 1676016287568,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "c6538528d7165a7765add222e2a830bbcef77827",
          "message": "Merge pull request #375 from galacticusorg/font2008PhysPlaus\n\nIn `hotHaloRamPressureStrippingFont2008` only solve ram pressure stripping for physically-plausible systems",
          "timestamp": "2023-02-10T21:42:23Z",
          "tree_id": "a82af3a6861ece24eca14dcce054e2e2006d2532",
          "url": "https://github.com/galacticusorg/galacticus/commit/c6538528d7165a7765add222e2a830bbcef77827"
        },
        "date": 1676076003605,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.283,
            "unit": "seconds",
            "range": 0.0962917441943925
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
          "id": "c6538528d7165a7765add222e2a830bbcef77827",
          "message": "Merge pull request #375 from galacticusorg/font2008PhysPlaus\n\nIn `hotHaloRamPressureStrippingFont2008` only solve ram pressure stripping for physically-plausible systems",
          "timestamp": "2023-02-10T21:42:23Z",
          "tree_id": "a82af3a6861ece24eca14dcce054e2e2006d2532",
          "url": "https://github.com/galacticusorg/galacticus/commit/c6538528d7165a7765add222e2a830bbcef77827"
        },
        "date": 1676076013614,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "75a0ac0981ddda3703102dca6150278517e269ba",
          "message": "Merge pull request #376 from galacticusorg/buildProfile\n\nAdd build profiling to the CI/CD workflow",
          "timestamp": "2023-02-11T06:22:22Z",
          "tree_id": "811523dcfe2f10d750b54e40b65b8efd270b642a",
          "url": "https://github.com/galacticusorg/galacticus/commit/75a0ac0981ddda3703102dca6150278517e269ba"
        },
        "date": 1676105769134,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 64.088,
            "unit": "seconds",
            "range": 0.333987424913239
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
          "id": "75a0ac0981ddda3703102dca6150278517e269ba",
          "message": "Merge pull request #376 from galacticusorg/buildProfile\n\nAdd build profiling to the CI/CD workflow",
          "timestamp": "2023-02-11T06:22:22Z",
          "tree_id": "811523dcfe2f10d750b54e40b65b8efd270b642a",
          "url": "https://github.com/galacticusorg/galacticus/commit/75a0ac0981ddda3703102dca6150278517e269ba"
        },
        "date": 1676105776901,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "5dfbff91ce1e09e118e4811239f99467212a6881",
          "message": "Merge pull request #378 from galacticusorg/coolingFunctionSummationFix\n\nAvoid unnecessary divide/multiply in summation cooling function",
          "timestamp": "2023-02-16T00:50:23Z",
          "tree_id": "e05944b249e790eeb048954c55cc42b3fc1211aa",
          "url": "https://github.com/galacticusorg/galacticus/commit/5dfbff91ce1e09e118e4811239f99467212a6881"
        },
        "date": 1676524750888,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 50.132,
            "unit": "seconds",
            "range": 0.0795210663915013
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
          "id": "5dfbff91ce1e09e118e4811239f99467212a6881",
          "message": "Merge pull request #378 from galacticusorg/coolingFunctionSummationFix\n\nAvoid unnecessary divide/multiply in summation cooling function",
          "timestamp": "2023-02-16T00:50:23Z",
          "tree_id": "e05944b249e790eeb048954c55cc42b3fc1211aa",
          "url": "https://github.com/galacticusorg/galacticus/commit/5dfbff91ce1e09e118e4811239f99467212a6881"
        },
        "date": 1676524758795,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3717638717242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.689755316789,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25371.9758848307,
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
          "id": "6c39c5914b27d870f1465fff2f117631f7e6a149",
          "message": "Merge pull request #379 from galacticusorg/linAlgMemoryLeak\n\nAvoid memory leaks with interpolators and matrices",
          "timestamp": "2023-02-17T04:03:11Z",
          "tree_id": "877820fc00d2fff92c170531a3689e9c6456d55d",
          "url": "https://github.com/galacticusorg/galacticus/commit/6c39c5914b27d870f1465fff2f117631f7e6a149"
        },
        "date": 1676616342200,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.731,
            "unit": "seconds",
            "range": 0.0328161545569748
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
          "id": "6c39c5914b27d870f1465fff2f117631f7e6a149",
          "message": "Merge pull request #379 from galacticusorg/linAlgMemoryLeak\n\nAvoid memory leaks with interpolators and matrices",
          "timestamp": "2023-02-17T04:03:11Z",
          "tree_id": "877820fc00d2fff92c170531a3689e9c6456d55d",
          "url": "https://github.com/galacticusorg/galacticus/commit/6c39c5914b27d870f1465fff2f117631f7e6a149"
        },
        "date": 1676616350599,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3717638717242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.689755316789,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25371.9758848307,
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
          "id": "33856823b9a43d7a1103f883a408daa1c50fd484",
          "message": "Merge pull request #380 from galacticusorg/threadIOSafety\n\nImprovement to internal I/O thread safety",
          "timestamp": "2023-02-17T23:06:54Z",
          "tree_id": "e70ef42d6d1115147630e5fe5d2c9ca06fe4ac7e",
          "url": "https://github.com/galacticusorg/galacticus/commit/33856823b9a43d7a1103f883a408daa1c50fd484"
        },
        "date": 1676685532854,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.776,
            "unit": "seconds",
            "range": 0.025345611060205
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
          "id": "33856823b9a43d7a1103f883a408daa1c50fd484",
          "message": "Merge pull request #380 from galacticusorg/threadIOSafety\n\nImprovement to internal I/O thread safety",
          "timestamp": "2023-02-17T23:06:54Z",
          "tree_id": "e70ef42d6d1115147630e5fe5d2c9ca06fe4ac7e",
          "url": "https://github.com/galacticusorg/galacticus/commit/33856823b9a43d7a1103f883a408daa1c50fd484"
        },
        "date": 1676685541144,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3717638717242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.689755316789,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25371.9758848307,
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
          "id": "ffd38134fc6be768eb7091621c4c84dfd0b273e4",
          "message": "Merge pull request #381 from galacticusorg/constrainedBranches\n\nFix constrained branch indicator propagation",
          "timestamp": "2023-02-19T17:09:52Z",
          "tree_id": "a7aacc54af6d0937f4279a0767c0c6da855cd283",
          "url": "https://github.com/galacticusorg/galacticus/commit/ffd38134fc6be768eb7091621c4c84dfd0b273e4"
        },
        "date": 1676836666666,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 50.196,
            "unit": "seconds",
            "range": 0.0442538133941168
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
          "id": "ffd38134fc6be768eb7091621c4c84dfd0b273e4",
          "message": "Merge pull request #381 from galacticusorg/constrainedBranches\n\nFix constrained branch indicator propagation",
          "timestamp": "2023-02-19T17:09:52Z",
          "tree_id": "a7aacc54af6d0937f4279a0767c0c6da855cd283",
          "url": "https://github.com/galacticusorg/galacticus/commit/ffd38134fc6be768eb7091621c4c84dfd0b273e4"
        },
        "date": 1676836674642,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "5df11a52eb6d0323b8ed5cf4a54d5c66442cec4b",
          "message": "fix: Avoid floating point exceptions for empty mass distributions",
          "timestamp": "2023-02-22T14:54:24Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/5df11a52eb6d0323b8ed5cf4a54d5c66442cec4b"
        },
        "date": 1677089037847,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.415,
            "unit": "seconds",
            "range": 0.0724327274097912
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
          "id": "5df11a52eb6d0323b8ed5cf4a54d5c66442cec4b",
          "message": "fix: Avoid floating point exceptions for empty mass distributions",
          "timestamp": "2023-02-22T14:54:24Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/5df11a52eb6d0323b8ed5cf4a54d5c66442cec4b"
        },
        "date": 1677089047935,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "be4aad9836b68a5a98948480fb8e76054d353acf",
          "message": "Merge pull request #382 from galacticusorg/fixSpelling\n\nFix spelling mistakes and typos in documentation",
          "timestamp": "2023-02-23T14:14:20Z",
          "tree_id": "ba78be22f9c570aa2dd292b6b9935a4fb70f663d",
          "url": "https://github.com/galacticusorg/galacticus/commit/be4aad9836b68a5a98948480fb8e76054d353acf"
        },
        "date": 1677171077297,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 52.464,
            "unit": "seconds",
            "range": 0.0570823965848302
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
          "id": "be4aad9836b68a5a98948480fb8e76054d353acf",
          "message": "Merge pull request #382 from galacticusorg/fixSpelling\n\nFix spelling mistakes and typos in documentation",
          "timestamp": "2023-02-23T14:14:20Z",
          "tree_id": "ba78be22f9c570aa2dd292b6b9935a4fb70f663d",
          "url": "https://github.com/galacticusorg/galacticus/commit/be4aad9836b68a5a98948480fb8e76054d353acf"
        },
        "date": 1677171088740,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3717638717242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.689755316789,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25371.9758848307,
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
          "id": "c92abfa868f500de89dd190dc8329b96df64219e",
          "message": "Merge pull request #386 from galacticusorg/dependencies\n\nUse a file to store dependency versions",
          "timestamp": "2023-02-25T16:43:43Z",
          "tree_id": "59751bd8805093a1dcdfdb639eb8391c74873a67",
          "url": "https://github.com/galacticusorg/galacticus/commit/c92abfa868f500de89dd190dc8329b96df64219e"
        },
        "date": 1677359080363,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 60.721,
            "unit": "seconds",
            "range": 0.107111624019647
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
          "id": "c92abfa868f500de89dd190dc8329b96df64219e",
          "message": "Merge pull request #386 from galacticusorg/dependencies\n\nUse a file to store dependency versions",
          "timestamp": "2023-02-25T16:43:43Z",
          "tree_id": "59751bd8805093a1dcdfdb639eb8391c74873a67",
          "url": "https://github.com/galacticusorg/galacticus/commit/c92abfa868f500de89dd190dc8329b96df64219e"
        },
        "date": 1677359087455,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3717638717242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.689755316789,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25371.9758848307,
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
          "id": "631a147700ebbf8bc4744a98bad28567a3dabf6c",
          "message": "fix: Merge branch 'master' of github.com:galacticusorg/galacticus",
          "timestamp": "2023-02-27T07:30:46-08:00",
          "tree_id": "a7057bc0913f2a42e0a192bed24e712be9fa715e",
          "url": "https://github.com/galacticusorg/galacticus/commit/631a147700ebbf8bc4744a98bad28567a3dabf6c"
        },
        "date": 1677524645621,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 54.088,
            "unit": "seconds",
            "range": 0.0894851943059495
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
          "id": "631a147700ebbf8bc4744a98bad28567a3dabf6c",
          "message": "fix: Merge branch 'master' of github.com:galacticusorg/galacticus",
          "timestamp": "2023-02-27T07:30:46-08:00",
          "tree_id": "a7057bc0913f2a42e0a192bed24e712be9fa715e",
          "url": "https://github.com/galacticusorg/galacticus/commit/631a147700ebbf8bc4744a98bad28567a3dabf6c"
        },
        "date": 1677524654101,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "fe442e0f2c47f5592a92406409eafed4a970141e",
          "message": "fix: Read dependencies versions from file when building containerized Galacticus",
          "timestamp": "2023-02-27T20:49:24Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/fe442e0f2c47f5592a92406409eafed4a970141e"
        },
        "date": 1677567482028,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 61.037,
            "unit": "seconds",
            "range": 0.125339937769529
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
          "id": "fe442e0f2c47f5592a92406409eafed4a970141e",
          "message": "fix: Read dependencies versions from file when building containerized Galacticus",
          "timestamp": "2023-02-27T20:49:24Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/fe442e0f2c47f5592a92406409eafed4a970141e"
        },
        "date": 1677567491399,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "5da960e711eba8a2b71d2d8e9318fe985a87d4bd",
          "message": "Merge pull request #392 from galacticusorg/allowedParametersScaling\n\nImprove run-time scaling of `allowedParameters()` methods for very large classes",
          "timestamp": "2023-02-28T21:49:58Z",
          "tree_id": "7193342c4ca93fe0a1dca18d1ceecd0943db1548",
          "url": "https://github.com/galacticusorg/galacticus/commit/5da960e711eba8a2b71d2d8e9318fe985a87d4bd"
        },
        "date": 1677656063500,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 52.221,
            "unit": "seconds",
            "range": 0.0257856549265766
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
          "id": "5da960e711eba8a2b71d2d8e9318fe985a87d4bd",
          "message": "Merge pull request #392 from galacticusorg/allowedParametersScaling\n\nImprove run-time scaling of `allowedParameters()` methods for very large classes",
          "timestamp": "2023-02-28T21:49:58Z",
          "tree_id": "7193342c4ca93fe0a1dca18d1ceecd0943db1548",
          "url": "https://github.com/galacticusorg/galacticus/commit/5da960e711eba8a2b71d2d8e9318fe985a87d4bd"
        },
        "date": 1677656071693,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "a46cb17c8d15102eb0497fc06938117c4091a061",
          "message": "fix: Provide a `massDistribution` method for the \"very simple size\" disk component",
          "timestamp": "2023-02-28T23:29:59Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/a46cb17c8d15102eb0497fc06938117c4091a061"
        },
        "date": 1677656073960,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 52.042,
            "unit": "seconds",
            "range": 0.0313942669927601
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
          "id": "a46cb17c8d15102eb0497fc06938117c4091a061",
          "message": "fix: Provide a `massDistribution` method for the \"very simple size\" disk component",
          "timestamp": "2023-02-28T23:29:59Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/a46cb17c8d15102eb0497fc06938117c4091a061"
        },
        "date": 1677656083160,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3684662387286,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.6611853229731,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25619.8050208255,
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
          "id": "5da960e711eba8a2b71d2d8e9318fe985a87d4bd",
          "message": "Merge pull request #392 from galacticusorg/allowedParametersScaling\n\nImprove run-time scaling of `allowedParameters()` methods for very large classes",
          "timestamp": "2023-02-28T21:49:58Z",
          "tree_id": "7193342c4ca93fe0a1dca18d1ceecd0943db1548",
          "url": "https://github.com/galacticusorg/galacticus/commit/5da960e711eba8a2b71d2d8e9318fe985a87d4bd"
        },
        "date": 1677684037544,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 52.221,
            "unit": "seconds",
            "range": 0.0257856549265766
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
          "id": "5da960e711eba8a2b71d2d8e9318fe985a87d4bd",
          "message": "Merge pull request #392 from galacticusorg/allowedParametersScaling\n\nImprove run-time scaling of `allowedParameters()` methods for very large classes",
          "timestamp": "2023-02-28T21:49:58Z",
          "tree_id": "7193342c4ca93fe0a1dca18d1ceecd0943db1548",
          "url": "https://github.com/galacticusorg/galacticus/commit/5da960e711eba8a2b71d2d8e9318fe985a87d4bd"
        },
        "date": 1677684046852,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "8eac31dae114401ae9a0116562e45dd151b2a12f",
          "message": "fix: Use correct name for default \"very simple size\" disk mass distribution parameter",
          "timestamp": "2023-03-01T15:33:56Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/8eac31dae114401ae9a0116562e45dd151b2a12f"
        },
        "date": 1677709723620,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 59.608,
            "unit": "seconds",
            "range": 0.168260512301223
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
          "id": "8eac31dae114401ae9a0116562e45dd151b2a12f",
          "message": "fix: Use correct name for default \"very simple size\" disk mass distribution parameter",
          "timestamp": "2023-03-01T15:33:56Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/8eac31dae114401ae9a0116562e45dd151b2a12f"
        },
        "date": 1677709732851,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3717638717242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.689755316789,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25371.9758848307,
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
          "id": "c2b0d3661e653ceacfae8a22db2abe26b15951a1",
          "message": "fix: Import required `uniq` function",
          "timestamp": "2023-03-01T16:50:23-08:00",
          "tree_id": "81a13b3a78dec542ca02e9598d0f1aa48278237b",
          "url": "https://github.com/galacticusorg/galacticus/commit/c2b0d3661e653ceacfae8a22db2abe26b15951a1"
        },
        "date": 1677740570234,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 60.385,
            "unit": "seconds",
            "range": 0.023590252225957
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
          "id": "c2b0d3661e653ceacfae8a22db2abe26b15951a1",
          "message": "fix: Import required `uniq` function",
          "timestamp": "2023-03-01T16:50:23-08:00",
          "tree_id": "81a13b3a78dec542ca02e9598d0f1aa48278237b",
          "url": "https://github.com/galacticusorg/galacticus/commit/c2b0d3661e653ceacfae8a22db2abe26b15951a1"
        },
        "date": 1677740578936,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3717638717242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.689755316789,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25371.9758848307,
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
          "id": "ca3da93e4020b654e59f77aef5fdd28621f61f22",
          "message": "Merge pull request #396 from galacticusorg/unknownParameters\n\nReport the path to a parameter when resporting an unknown parameter name",
          "timestamp": "2023-03-02T15:35:14Z",
          "tree_id": "5d3e4fa9a7502eefb30d2aa2d62ac6b8e3e23a87",
          "url": "https://github.com/galacticusorg/galacticus/commit/ca3da93e4020b654e59f77aef5fdd28621f61f22"
        },
        "date": 1677781094396,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 50.767,
            "unit": "seconds",
            "range": 0.0312745903261058
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
          "id": "ca3da93e4020b654e59f77aef5fdd28621f61f22",
          "message": "Merge pull request #396 from galacticusorg/unknownParameters\n\nReport the path to a parameter when resporting an unknown parameter name",
          "timestamp": "2023-03-02T15:35:14Z",
          "tree_id": "5d3e4fa9a7502eefb30d2aa2d62ac6b8e3e23a87",
          "url": "https://github.com/galacticusorg/galacticus/commit/ca3da93e4020b654e59f77aef5fdd28621f61f22"
        },
        "date": 1677781104813,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "4c055d6316922ca14b4a03eb253d08e43d6f33a5",
          "message": "Merge pull request #398 from galacticusorg/sidmTabulationFix\n\nAllow isothermal SIDM density profile to extend the range of tabulated solutions as needed",
          "timestamp": "2023-03-02T23:35:00Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/4c055d6316922ca14b4a03eb253d08e43d6f33a5"
        },
        "date": 1677825874729,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 50.406,
            "unit": "seconds",
            "range": 0.0461779168009253
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
          "id": "4c055d6316922ca14b4a03eb253d08e43d6f33a5",
          "message": "Merge pull request #398 from galacticusorg/sidmTabulationFix\n\nAllow isothermal SIDM density profile to extend the range of tabulated solutions as needed",
          "timestamp": "2023-03-02T23:35:00Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/4c055d6316922ca14b4a03eb253d08e43d6f33a5"
        },
        "date": 1677825885717,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "6249994d04f1b4e8ee292404787c4f2b95fadbe8",
          "message": "Merge pull request #399 from galacticusorg/memoryLeakFix\n\nFix memory leaks",
          "timestamp": "2023-03-07T00:23:44Z",
          "tree_id": "fb353d02fb8451de9400f8bd8b37aac9b385e309",
          "url": "https://github.com/galacticusorg/galacticus/commit/6249994d04f1b4e8ee292404787c4f2b95fadbe8"
        },
        "date": 1678168150780,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 60.909,
            "unit": "seconds",
            "range": 0.229283449031782
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
          "id": "6249994d04f1b4e8ee292404787c4f2b95fadbe8",
          "message": "Merge pull request #399 from galacticusorg/memoryLeakFix\n\nFix memory leaks",
          "timestamp": "2023-03-07T00:23:44Z",
          "tree_id": "fb353d02fb8451de9400f8bd8b37aac9b385e309",
          "url": "https://github.com/galacticusorg/galacticus/commit/6249994d04f1b4e8ee292404787c4f2b95fadbe8"
        },
        "date": 1678168158533,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.399234412822,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.7769036125367,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23989.8072123588,
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
          "id": "b97188142b9e3ec8e4e296e6592871911355231e",
          "message": "fix: Catch zero mass black holes",
          "timestamp": "2023-03-07T15:46:30Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/b97188142b9e3ec8e4e296e6592871911355231e"
        },
        "date": 1678218637000,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.636,
            "unit": "seconds",
            "range": 0.441405029423275
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
          "id": "b97188142b9e3ec8e4e296e6592871911355231e",
          "message": "fix: Catch zero mass black holes",
          "timestamp": "2023-03-07T15:46:30Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/b97188142b9e3ec8e4e296e6592871911355231e"
        },
        "date": 1678218647546,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "dd67e743f57f3ff053a94226584827feb92a7003",
          "message": "Merge pull request #400 from galacticusorg/galacticStructureStack\n\nUse a linked-list of `galacticStructureStandard` state",
          "timestamp": "2023-03-09T15:21:48Z",
          "tree_id": "e5d9c4f0c6cd551a7ca6b415d54f4884a46ac5f8",
          "url": "https://github.com/galacticusorg/galacticus/commit/dd67e743f57f3ff053a94226584827feb92a7003"
        },
        "date": 1678399903242,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 57.967,
            "unit": "seconds",
            "range": 0.109709160966726
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
          "id": "dd67e743f57f3ff053a94226584827feb92a7003",
          "message": "Merge pull request #400 from galacticusorg/galacticStructureStack\n\nUse a linked-list of `galacticStructureStandard` state",
          "timestamp": "2023-03-09T15:21:48Z",
          "tree_id": "e5d9c4f0c6cd551a7ca6b415d54f4884a46ac5f8",
          "url": "https://github.com/galacticusorg/galacticus/commit/dd67e743f57f3ff053a94226584827feb92a7003"
        },
        "date": 1678399911203,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3717638717242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.689755316789,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25371.9758848307,
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
          "id": "b844981ce953c7b895be99c7f94b98e28367e81c",
          "message": "feat: Add a polynomial systematic in halo mass for the COSMO SHMR output analysis",
          "timestamp": "2023-03-10T17:44:20Z",
          "tree_id": "8b2468ff0a9d8e851babf87d6303135bf488f42d",
          "url": "https://github.com/galacticusorg/galacticus/commit/b844981ce953c7b895be99c7f94b98e28367e81c"
        },
        "date": 1678483317094,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.542,
            "unit": "seconds",
            "range": 0.109332520322458
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
          "id": "b844981ce953c7b895be99c7f94b98e28367e81c",
          "message": "feat: Add a polynomial systematic in halo mass for the COSMO SHMR output analysis",
          "timestamp": "2023-03-10T17:44:20Z",
          "tree_id": "8b2468ff0a9d8e851babf87d6303135bf488f42d",
          "url": "https://github.com/galacticusorg/galacticus/commit/b844981ce953c7b895be99c7f94b98e28367e81c"
        },
        "date": 1678483325531,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.399234412822,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.7769036125367,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23989.8072123588,
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
          "id": "2b7e265e029d5d690719f32c67055d03deb4a843",
          "message": "feat: Use `massDistribution` directly in several `galacticStructure` methods",
          "timestamp": "2023-03-11T06:29:37Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/2b7e265e029d5d690719f32c67055d03deb4a843"
        },
        "date": 1678527946441,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 50.826,
            "unit": "seconds",
            "range": 0.0281496003524077
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
          "id": "2b7e265e029d5d690719f32c67055d03deb4a843",
          "message": "feat: Use `massDistribution` directly in several `galacticStructure` methods",
          "timestamp": "2023-03-11T06:29:37Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/2b7e265e029d5d690719f32c67055d03deb4a843"
        },
        "date": 1678527957062,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "7fd32e0c9ba487562e5f59c37469bbf00fcb2e90",
          "message": "fix: Correct iterator limit for allowed parameters",
          "timestamp": "2023-03-14T17:30:11Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/7fd32e0c9ba487562e5f59c37469bbf00fcb2e90"
        },
        "date": 1678827296595,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 55.895,
            "unit": "seconds",
            "range": 0.0558972271219394
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
          "id": "7fd32e0c9ba487562e5f59c37469bbf00fcb2e90",
          "message": "fix: Correct iterator limit for allowed parameters",
          "timestamp": "2023-03-14T17:30:11Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/7fd32e0c9ba487562e5f59c37469bbf00fcb2e90"
        },
        "date": 1678827307258,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "aff1fe9062b45f6bc502a6c2c42022abe2619640",
          "message": "fix: Set permissions for Docker push",
          "timestamp": "2023-03-17T23:17:20Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/aff1fe9062b45f6bc502a6c2c42022abe2619640"
        },
        "date": 1679112489072,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 62.622,
            "unit": "seconds",
            "range": 0.133265149232817
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
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "aff1fe9062b45f6bc502a6c2c42022abe2619640",
          "message": "fix: Set permissions for Docker push",
          "timestamp": "2023-03-17T23:17:20Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/aff1fe9062b45f6bc502a6c2c42022abe2619640"
        },
        "date": 1679112497153,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3684662387286,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.6611853229731,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25619.8050208255,
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
          "id": "37b460897d0ece5c27aa90ab4a35048eab3ff438",
          "message": "fix: Switch from DockerHub to GitHub Container Registry",
          "timestamp": "2023-03-18T20:13:18Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/37b460897d0ece5c27aa90ab4a35048eab3ff438"
        },
        "date": 1679181520021,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.524,
            "unit": "seconds",
            "range": 0.0352476949601105
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
          "id": "37b460897d0ece5c27aa90ab4a35048eab3ff438",
          "message": "fix: Switch from DockerHub to GitHub Container Registry",
          "timestamp": "2023-03-18T20:13:18Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/37b460897d0ece5c27aa90ab4a35048eab3ff438"
        },
        "date": 1679181529845,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "5612535ab93bdce2be244d3f0b7cdd56cb01776d",
          "message": "fix: Correct spelling error",
          "timestamp": "2023-03-21T07:22:27-07:00",
          "tree_id": "4b6e681d1a462c26f7757cc3cc4ffde1977e189c",
          "url": "https://github.com/galacticusorg/galacticus/commit/5612535ab93bdce2be244d3f0b7cdd56cb01776d"
        },
        "date": 1679418448437,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 62.524,
            "unit": "seconds",
            "range": 0.414275753574641
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
          "id": "5612535ab93bdce2be244d3f0b7cdd56cb01776d",
          "message": "fix: Correct spelling error",
          "timestamp": "2023-03-21T07:22:27-07:00",
          "tree_id": "4b6e681d1a462c26f7757cc3cc4ffde1977e189c",
          "url": "https://github.com/galacticusorg/galacticus/commit/5612535ab93bdce2be244d3f0b7cdd56cb01776d"
        },
        "date": 1679418456880,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "986b15bbc3b1271fb5f837083f6f5574861397d2",
          "message": "fix: Update Dockerfile to pull build environment from GitHub Container Registry",
          "timestamp": "2023-03-21T19:49:52Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/986b15bbc3b1271fb5f837083f6f5574861397d2"
        },
        "date": 1679438172408,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 48.971,
            "unit": "seconds",
            "range": 0.0357337375591751
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
          "id": "986b15bbc3b1271fb5f837083f6f5574861397d2",
          "message": "fix: Update Dockerfile to pull build environment from GitHub Container Registry",
          "timestamp": "2023-03-21T19:49:52Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/986b15bbc3b1271fb5f837083f6f5574861397d2"
        },
        "date": 1679438182466,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "2217d31cd95e04211fa073ad037d35a67ea0b552",
          "message": "fix: Correct units and description for SED output",
          "timestamp": "2023-03-21T15:37:44-07:00",
          "tree_id": "bdc526d57e50154546203f74d0c0476f4068b6bc",
          "url": "https://github.com/galacticusorg/galacticus/commit/2217d31cd95e04211fa073ad037d35a67ea0b552"
        },
        "date": 1679448281063,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 50.916,
            "unit": "seconds",
            "range": 0.0385019480018966
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
          "id": "2217d31cd95e04211fa073ad037d35a67ea0b552",
          "message": "fix: Correct units and description for SED output",
          "timestamp": "2023-03-21T15:37:44-07:00",
          "tree_id": "bdc526d57e50154546203f74d0c0476f4068b6bc",
          "url": "https://github.com/galacticusorg/galacticus/commit/2217d31cd95e04211fa073ad037d35a67ea0b552"
        },
        "date": 1679448289212,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "ad2c49917e449b5c27d1127d9c9830c82f446c00",
          "message": "feat: Update the empirical elliptical galaxy model to allow the radius to be specified directly",
          "timestamp": "2023-03-22T14:57:44-07:00",
          "tree_id": "cb5709c9cf0203aa8561c5253f18027e3442d7a1",
          "url": "https://github.com/galacticusorg/galacticus/commit/ad2c49917e449b5c27d1127d9c9830c82f446c00"
        },
        "date": 1679531801877,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.372,
            "unit": "seconds",
            "range": 0.0325207625988818
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
          "id": "ad2c49917e449b5c27d1127d9c9830c82f446c00",
          "message": "feat: Update the empirical elliptical galaxy model to allow the radius to be specified directly",
          "timestamp": "2023-03-22T14:57:44-07:00",
          "tree_id": "cb5709c9cf0203aa8561c5253f18027e3442d7a1",
          "url": "https://github.com/galacticusorg/galacticus/commit/ad2c49917e449b5c27d1127d9c9830c82f446c00"
        },
        "date": 1679531810153,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "9a481bc04207a47828a92ae6485f377a2f2794fa",
          "message": "feat: Refactor density, spherically-averaged-density, and surface density in `galacticStructure` class\n\nThese now make use of the `massDistribution` object.",
          "timestamp": "2023-03-23T21:52:10Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/9a481bc04207a47828a92ae6485f377a2f2794fa"
        },
        "date": 1679619742612,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 55.094,
            "unit": "seconds",
            "range": 0.0616636035278551
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
          "id": "9a481bc04207a47828a92ae6485f377a2f2794fa",
          "message": "feat: Refactor density, spherically-averaged-density, and surface density in `galacticStructure` class\n\nThese now make use of the `massDistribution` object.",
          "timestamp": "2023-03-23T21:52:10Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/9a481bc04207a47828a92ae6485f377a2f2794fa"
        },
        "date": 1679619750858,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.235741215839,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.7337765334936,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25102.793154175,
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
          "id": "d6f6c31b9dbba869c7d30ffcc0db1bec2a58640d",
          "message": "feat: Add an option to make ram pressure stripping radius solver failures non-fatal\n\nAlso removes debugging code.",
          "timestamp": "2023-03-24T00:02:11Z",
          "tree_id": "12af966ee52995c3ba0684884d2c1ea0e7905571",
          "url": "https://github.com/galacticusorg/galacticus/commit/d6f6c31b9dbba869c7d30ffcc0db1bec2a58640d"
        },
        "date": 1679628774319,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 48.875,
            "unit": "seconds",
            "range": 0.0451054320452721
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
          "id": "d6f6c31b9dbba869c7d30ffcc0db1bec2a58640d",
          "message": "feat: Add an option to make ram pressure stripping radius solver failures non-fatal\n\nAlso removes debugging code.",
          "timestamp": "2023-03-24T00:02:11Z",
          "tree_id": "12af966ee52995c3ba0684884d2c1ea0e7905571",
          "url": "https://github.com/galacticusorg/galacticus/commit/d6f6c31b9dbba869c7d30ffcc0db1bec2a58640d"
        },
        "date": 1679628783196,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3717638717242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.689755316789,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25371.9758848307,
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
          "id": "09d0c2ca3c666123da18e8fd807809e23f8ae7fa",
          "message": "fix: Workaround perturber extent calculation until dark matter profiles are implemented as `massDistributionClass` objects",
          "timestamp": "2023-03-30T19:38:09Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/09d0c2ca3c666123da18e8fd807809e23f8ae7fa"
        },
        "date": 1680216438252,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 54.091,
            "unit": "seconds",
            "range": 0.0473803756854681
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
          "id": "09d0c2ca3c666123da18e8fd807809e23f8ae7fa",
          "message": "fix: Workaround perturber extent calculation until dark matter profiles are implemented as `massDistributionClass` objects",
          "timestamp": "2023-03-30T19:38:09Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/09d0c2ca3c666123da18e8fd807809e23f8ae7fa"
        },
        "date": 1680216445992,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "0f60eb30a53a2afeac453171d801e6dd63621bae",
          "message": "fix: Remove FFTW3 library from linker if not available",
          "timestamp": "2023-03-30T19:23:34-07:00",
          "tree_id": "e009512541cfffe168d7818f2e8653403017b5b0",
          "url": "https://github.com/galacticusorg/galacticus/commit/0f60eb30a53a2afeac453171d801e6dd63621bae"
        },
        "date": 1680250295207,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.57,
            "unit": "seconds",
            "range": 0.025139610178064
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
          "id": "0f60eb30a53a2afeac453171d801e6dd63621bae",
          "message": "fix: Remove FFTW3 library from linker if not available",
          "timestamp": "2023-03-30T19:23:34-07:00",
          "tree_id": "e009512541cfffe168d7818f2e8653403017b5b0",
          "url": "https://github.com/galacticusorg/galacticus/commit/0f60eb30a53a2afeac453171d801e6dd63621bae"
        },
        "date": 1680250306195,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "bfef8ffd1beab55042e0d5a024dee0b00a14c60a",
          "message": "fix(style): Formatting only",
          "timestamp": "2023-03-31T02:40:29Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/bfef8ffd1beab55042e0d5a024dee0b00a14c60a"
        },
        "date": 1680268653756,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 52.983,
            "unit": "seconds",
            "range": 0.0290189593181578
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
          "id": "bfef8ffd1beab55042e0d5a024dee0b00a14c60a",
          "message": "fix(style): Formatting only",
          "timestamp": "2023-03-31T02:40:29Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/bfef8ffd1beab55042e0d5a024dee0b00a14c60a"
        },
        "date": 1680268664446,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3717638717242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.689755316789,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25371.9758848307,
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
          "id": "cc4eee10982e01c5d5bced0d79a894ac1462c74d",
          "message": "fix: Filter out unrealistically tiny galaxies when computing the Local Group mass-size relation",
          "timestamp": "2023-04-03T23:54:26Z",
          "tree_id": "645ad338920c918458b3e3e346279b8f9b1cddbc",
          "url": "https://github.com/galacticusorg/galacticus/commit/cc4eee10982e01c5d5bced0d79a894ac1462c74d"
        },
        "date": 1680584660955,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.956,
            "unit": "seconds",
            "range": 0.071835924160388
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
          "id": "cc4eee10982e01c5d5bced0d79a894ac1462c74d",
          "message": "fix: Filter out unrealistically tiny galaxies when computing the Local Group mass-size relation",
          "timestamp": "2023-04-03T23:54:26Z",
          "tree_id": "645ad338920c918458b3e3e346279b8f9b1cddbc",
          "url": "https://github.com/galacticusorg/galacticus/commit/cc4eee10982e01c5d5bced0d79a894ac1462c74d"
        },
        "date": 1680584668833,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "801cb4cb8e9fe6722d5ffb6599dd22d4a88ee456",
          "message": "fix: Correct iterator limit for allowed parameters",
          "timestamp": "2023-04-04T18:15:15Z",
          "tree_id": "a15d38aacfc4abde91fd0b62a8f74a893cadeced",
          "url": "https://github.com/galacticusorg/galacticus/commit/801cb4cb8e9fe6722d5ffb6599dd22d4a88ee456"
        },
        "date": 1680641646153,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 57.748,
            "unit": "seconds",
            "range": 0.237772159850382
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
          "id": "801cb4cb8e9fe6722d5ffb6599dd22d4a88ee456",
          "message": "fix: Correct iterator limit for allowed parameters",
          "timestamp": "2023-04-04T18:15:15Z",
          "tree_id": "a15d38aacfc4abde91fd0b62a8f74a893cadeced",
          "url": "https://github.com/galacticusorg/galacticus/commit/801cb4cb8e9fe6722d5ffb6599dd22d4a88ee456"
        },
        "date": 1680641654270,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.399234412822,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.7769036125367,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23989.8072123588,
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
          "id": "f24b6658837bbb9e810e1f71acf27311444e2220",
          "message": "fix: Filter out unphysical galaxies in the Local Group mass-velocity dispersion analysis",
          "timestamp": "2023-04-06T21:04:42Z",
          "tree_id": "90614c3c6c4000f4605125120c7ae64d82fe8387",
          "url": "https://github.com/galacticusorg/galacticus/commit/f24b6658837bbb9e810e1f71acf27311444e2220"
        },
        "date": 1680825349415,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.412,
            "unit": "seconds",
            "range": 0.0662389613449946
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
          "id": "f24b6658837bbb9e810e1f71acf27311444e2220",
          "message": "fix: Filter out unphysical galaxies in the Local Group mass-velocity dispersion analysis",
          "timestamp": "2023-04-06T21:04:42Z",
          "tree_id": "90614c3c6c4000f4605125120c7ae64d82fe8387",
          "url": "https://github.com/galacticusorg/galacticus/commit/f24b6658837bbb9e810e1f71acf27311444e2220"
        },
        "date": 1680825358516,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "defddd0c6c83876ed0c11959fbb15f94f9aeaa16",
          "message": "fix(style): Fix typo",
          "timestamp": "2023-04-13T14:34:28-07:00",
          "tree_id": "eb4177091b67b35b2b69f68a20c0a9dc252c6e43",
          "url": "https://github.com/galacticusorg/galacticus/commit/defddd0c6c83876ed0c11959fbb15f94f9aeaa16"
        },
        "date": 1681431714943,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.47,
            "unit": "seconds",
            "range": 0.125849116007961
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
          "id": "defddd0c6c83876ed0c11959fbb15f94f9aeaa16",
          "message": "fix(style): Fix typo",
          "timestamp": "2023-04-13T14:34:28-07:00",
          "tree_id": "eb4177091b67b35b2b69f68a20c0a9dc252c6e43",
          "url": "https://github.com/galacticusorg/galacticus/commit/defddd0c6c83876ed0c11959fbb15f94f9aeaa16"
        },
        "date": 1681431725268,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "00e8626cfc7b7c5273a28fef3e5316063c4d249c",
          "message": "fix: Add `--allow-run-as-root` for `valgrind`",
          "timestamp": "2023-04-16T23:48:34Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/00e8626cfc7b7c5273a28fef3e5316063c4d249c"
        },
        "date": 1681699207359,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 52.933,
            "unit": "seconds",
            "range": 0.160268836646381
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
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "00e8626cfc7b7c5273a28fef3e5316063c4d249c",
          "message": "fix: Add `--allow-run-as-root` for `valgrind`",
          "timestamp": "2023-04-16T23:48:34Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/00e8626cfc7b7c5273a28fef3e5316063c4d249c"
        },
        "date": 1681699219006,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "96daacc6391e2035fee11e4c21c405fac90746dd",
          "message": "fix: Handle malformed XML\n\nSometimes the XML produced by `valgrind` running under MPI can be malformed. Handle such cases by ignoring them.",
          "timestamp": "2023-04-17T05:07:44Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/96daacc6391e2035fee11e4c21c405fac90746dd"
        },
        "date": 1681718293029,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.487,
            "unit": "seconds",
            "range": 0.0289153938235245
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
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "96daacc6391e2035fee11e4c21c405fac90746dd",
          "message": "fix: Handle malformed XML\n\nSometimes the XML produced by `valgrind` running under MPI can be malformed. Handle such cases by ignoring them.",
          "timestamp": "2023-04-17T05:07:44Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/96daacc6391e2035fee11e4c21c405fac90746dd"
        },
        "date": 1681718302210,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3717638717242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.689755316789,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25371.9758848307,
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
          "id": "37b2de62363d5d06bec6c741ff4c79e3cadfe3bb",
          "message": "fix: Add a non-static build for use with `valgrind`",
          "timestamp": "2023-04-17T15:36:52Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/37b2de62363d5d06bec6c741ff4c79e3cadfe3bb"
        },
        "date": 1681755628910,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.506,
            "unit": "seconds",
            "range": 0.0605012396566364
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
          "id": "37b2de62363d5d06bec6c741ff4c79e3cadfe3bb",
          "message": "fix: Add a non-static build for use with `valgrind`",
          "timestamp": "2023-04-17T15:36:52Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/37b2de62363d5d06bec6c741ff4c79e3cadfe3bb"
        },
        "date": 1681755639592,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "6ebc45062686ef0809ad3ced3fb32692a6fa0cef",
          "message": "feat: Add option to ignore integration errors for line-of-sight velocity dispersions",
          "timestamp": "2023-04-18T02:44:45Z",
          "tree_id": "f4c3ac104865e3e7477666f39fdf15250169b01b",
          "url": "https://github.com/galacticusorg/galacticus/commit/6ebc45062686ef0809ad3ced3fb32692a6fa0cef"
        },
        "date": 1681796155707,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.081,
            "unit": "seconds",
            "range": 0.140871927650447
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
          "id": "6ebc45062686ef0809ad3ced3fb32692a6fa0cef",
          "message": "feat: Add option to ignore integration errors for line-of-sight velocity dispersions",
          "timestamp": "2023-04-18T02:44:45Z",
          "tree_id": "f4c3ac104865e3e7477666f39fdf15250169b01b",
          "url": "https://github.com/galacticusorg/galacticus/commit/6ebc45062686ef0809ad3ced3fb32692a6fa0cef"
        },
        "date": 1681796164389,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "b60e818869ea0bad7e2fcc2b9320cabbe02cf550",
          "message": "fix: Avoid outputting state if failures are non-fatal\n\nDebug state should only be output if we are about to abort because of a radius solver failure. If such failures are being ignored there is no need to output state.",
          "timestamp": "2023-04-19T23:49:35Z",
          "tree_id": "143810e6de1c30da807e3a8537a6544a3cd88959",
          "url": "https://github.com/galacticusorg/galacticus/commit/b60e818869ea0bad7e2fcc2b9320cabbe02cf550"
        },
        "date": 1681960460303,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.729,
            "unit": "seconds",
            "range": 0.0369986486244301
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
          "id": "b60e818869ea0bad7e2fcc2b9320cabbe02cf550",
          "message": "fix: Avoid outputting state if failures are non-fatal\n\nDebug state should only be output if we are about to abort because of a radius solver failure. If such failures are being ignored there is no need to output state.",
          "timestamp": "2023-04-19T23:49:35Z",
          "tree_id": "143810e6de1c30da807e3a8537a6544a3cd88959",
          "url": "https://github.com/galacticusorg/galacticus/commit/b60e818869ea0bad7e2fcc2b9320cabbe02cf550"
        },
        "date": 1681960468743,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3717638717242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.689755316789,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25371.9758848307,
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
          "id": "3828701f378f6883a1975d81949e9f28238302d9",
          "message": "Merge pull request #403 from galacticusorg/cloudyDependency\n\nUse the dependency file to specify the version of Cloudy to use",
          "timestamp": "2023-04-25T14:23:37Z",
          "tree_id": "5b111a285d94e7e0fd2f6d3bf2ed7d110ae8debf",
          "url": "https://github.com/galacticusorg/galacticus/commit/3828701f378f6883a1975d81949e9f28238302d9"
        },
        "date": 1682451995075,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.155,
            "unit": "seconds",
            "range": 0.0732700484502023
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
          "id": "3828701f378f6883a1975d81949e9f28238302d9",
          "message": "Merge pull request #403 from galacticusorg/cloudyDependency\n\nUse the dependency file to specify the version of Cloudy to use",
          "timestamp": "2023-04-25T14:23:37Z",
          "tree_id": "5b111a285d94e7e0fd2f6d3bf2ed7d110ae8debf",
          "url": "https://github.com/galacticusorg/galacticus/commit/3828701f378f6883a1975d81949e9f28238302d9"
        },
        "date": 1682452003198,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3717638717242,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.689755316789,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25371.9758848307,
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
          "id": "7f0f10f4c016f8cf0333daf81a4c784b78121df6",
          "message": "fix: Add `libdl` option in linker",
          "timestamp": "2023-04-26T05:24:18Z",
          "tree_id": "9b88298dff920d4e7fa11b7080540d83216f4f94",
          "url": "https://github.com/galacticusorg/galacticus/commit/7f0f10f4c016f8cf0333daf81a4c784b78121df6"
        },
        "date": 1682505946437,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.391,
            "unit": "seconds",
            "range": 0.0727523195503033
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
          "id": "7f0f10f4c016f8cf0333daf81a4c784b78121df6",
          "message": "fix: Add `libdl` option in linker",
          "timestamp": "2023-04-26T05:24:18Z",
          "tree_id": "9b88298dff920d4e7fa11b7080540d83216f4f94",
          "url": "https://github.com/galacticusorg/galacticus/commit/7f0f10f4c016f8cf0333daf81a4c784b78121df6"
        },
        "date": 1682505955226,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "915ddca2724ddc775b68646639e4100210c37245",
          "message": "fix: Correct `Dockerfile` syntax",
          "timestamp": "2023-04-26T17:36:55Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/915ddca2724ddc775b68646639e4100210c37245"
        },
        "date": 1682546169256,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 70.197,
            "unit": "seconds",
            "range": 3.62315278452341
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
          "id": "915ddca2724ddc775b68646639e4100210c37245",
          "message": "fix: Correct `Dockerfile` syntax",
          "timestamp": "2023-04-26T17:36:55Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/915ddca2724ddc775b68646639e4100210c37245"
        },
        "date": 1682546179665,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3684662387286,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.6611853229731,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25619.8050208255,
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
          "id": "6790efb69d319a54319aedee5e34836118291824",
          "message": "fix: Set a tree volume weight (of `1.0`) in the `mergerTreeConstructorFullySpecified` class",
          "timestamp": "2023-04-27T16:05:58Z",
          "tree_id": "61c50f26e43f6942e9e0bca0b35035c58e4a047b",
          "url": "https://github.com/galacticusorg/galacticus/commit/6790efb69d319a54319aedee5e34836118291824"
        },
        "date": 1682622196638,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 67.898,
            "unit": "seconds",
            "range": 0.411921837245655
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
          "id": "6790efb69d319a54319aedee5e34836118291824",
          "message": "fix: Set a tree volume weight (of `1.0`) in the `mergerTreeConstructorFullySpecified` class",
          "timestamp": "2023-04-27T16:05:58Z",
          "tree_id": "61c50f26e43f6942e9e0bca0b35035c58e4a047b",
          "url": "https://github.com/galacticusorg/galacticus/commit/6790efb69d319a54319aedee5e34836118291824"
        },
        "date": 1682622205618,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "b8451cd19add1a6cac8e7866ca0c4b06ce7c07eb",
          "message": "fix: Point broken links to the Internet Archive",
          "timestamp": "2023-05-04T06:46:36-07:00",
          "tree_id": "0741ca16be1c47084de04770bbb117423e804547",
          "url": "https://github.com/galacticusorg/galacticus/commit/b8451cd19add1a6cac8e7866ca0c4b06ce7c07eb"
        },
        "date": 1683218493273,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 58.787,
            "unit": "seconds",
            "range": 0.0977246130715213
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
          "id": "b8451cd19add1a6cac8e7866ca0c4b06ce7c07eb",
          "message": "fix: Point broken links to the Internet Archive",
          "timestamp": "2023-05-04T06:46:36-07:00",
          "tree_id": "0741ca16be1c47084de04770bbb117423e804547",
          "url": "https://github.com/galacticusorg/galacticus/commit/b8451cd19add1a6cac8e7866ca0c4b06ce7c07eb"
        },
        "date": 1683218500983,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.2466314592324,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.6533248976655,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24200.3940817449,
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
          "id": "f833d28846fc288fd519bef8d6977b5163d1af11",
          "message": "fix: Avoid attempts to get `git` revision information in Cloudy builds\n\nAs we build Cloudy within the Galacticus `datasets` repo seeking the `git` revision will erroneously get the revision of the `datasets` repo, which can break the build.",
          "timestamp": "2023-05-04T14:48:48-07:00",
          "tree_id": "79e85e410b340d691849451589862e74af5c83c1",
          "url": "https://github.com/galacticusorg/galacticus/commit/f833d28846fc288fd519bef8d6977b5163d1af11"
        },
        "date": 1683266466402,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 60.725,
            "unit": "seconds",
            "range": 0.194762676095309
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
          "id": "f833d28846fc288fd519bef8d6977b5163d1af11",
          "message": "fix: Avoid attempts to get `git` revision information in Cloudy builds\n\nAs we build Cloudy within the Galacticus `datasets` repo seeking the `git` revision will erroneously get the revision of the `datasets` repo, which can break the build.",
          "timestamp": "2023-05-04T14:48:48-07:00",
          "tree_id": "79e85e410b340d691849451589862e74af5c83c1",
          "url": "https://github.com/galacticusorg/galacticus/commit/f833d28846fc288fd519bef8d6977b5163d1af11"
        },
        "date": 1683266476301,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.397560243055,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8403056003349,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25153.545938522,
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
          "id": "b3e73beaca12b4a769d59ce6d66bc022013bae30",
          "message": "fix: Merge branch 'master' of github.com:galacticusorg/galacticus",
          "timestamp": "2023-05-05T14:11:29Z",
          "tree_id": "ba226965ee134069f86b3f7928eb5c74008630ad",
          "url": "https://github.com/galacticusorg/galacticus/commit/b3e73beaca12b4a769d59ce6d66bc022013bae30"
        },
        "date": 1683306929571,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.615,
            "unit": "seconds",
            "range": 0.0456344168359849
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
          "id": "b3e73beaca12b4a769d59ce6d66bc022013bae30",
          "message": "fix: Merge branch 'master' of github.com:galacticusorg/galacticus",
          "timestamp": "2023-05-05T14:11:29Z",
          "tree_id": "ba226965ee134069f86b3f7928eb5c74008630ad",
          "url": "https://github.com/galacticusorg/galacticus/commit/b3e73beaca12b4a769d59ce6d66bc022013bae30"
        },
        "date": 1683306937582,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.3275525628901,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8647158159627,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24145.1797388743,
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
          "id": "6e8872ed1bb5b58cb9a6098a21bb2d3513d0b906",
          "message": "Merge pull request #409 from galacticusorg/preInfallOrbit\n\nAdd support for tidal evolution of halos pre-infall",
          "timestamp": "2023-05-05T23:35:50Z",
          "tree_id": "2f7f7c98dfce6d9e7b3f427368ab6b04e53848aa",
          "url": "https://github.com/galacticusorg/galacticus/commit/6e8872ed1bb5b58cb9a6098a21bb2d3513d0b906"
        },
        "date": 1683340040500,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 64.024,
            "unit": "seconds",
            "range": 0.0976237675980581
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
          "id": "6e8872ed1bb5b58cb9a6098a21bb2d3513d0b906",
          "message": "Merge pull request #409 from galacticusorg/preInfallOrbit\n\nAdd support for tidal evolution of halos pre-infall",
          "timestamp": "2023-05-05T23:35:50Z",
          "tree_id": "2f7f7c98dfce6d9e7b3f427368ab6b04e53848aa",
          "url": "https://github.com/galacticusorg/galacticus/commit/6e8872ed1bb5b58cb9a6098a21bb2d3513d0b906"
        },
        "date": 1683340049208,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.4991777473382,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.0216488557544,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24211.0957851519,
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
          "id": "0c31f3f185c6c48d9bd8d2e02a400ff6a1ac6fd5",
          "message": "Merge pull request #410 from galacticusorg/postprocessSED\n\nImprove SED calculation",
          "timestamp": "2023-05-09T14:38:32Z",
          "tree_id": "d33777a848d3d194c139fa6bdcfa422a1e4863f1",
          "url": "https://github.com/galacticusorg/galacticus/commit/0c31f3f185c6c48d9bd8d2e02a400ff6a1ac6fd5"
        },
        "date": 1683657498204,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 59.052,
            "unit": "seconds",
            "range": 0.0669298139846011
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
          "id": "0c31f3f185c6c48d9bd8d2e02a400ff6a1ac6fd5",
          "message": "Merge pull request #410 from galacticusorg/postprocessSED\n\nImprove SED calculation",
          "timestamp": "2023-05-09T14:38:32Z",
          "tree_id": "d33777a848d3d194c139fa6bdcfa422a1e4863f1",
          "url": "https://github.com/galacticusorg/galacticus/commit/0c31f3f185c6c48d9bd8d2e02a400ff6a1ac6fd5"
        },
        "date": 1683657506737,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.4991777473382,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.0216488557544,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24211.0957851519,
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
          "id": "775ba693fe53c6b5354263997e74dde82922ef39",
          "message": "fix(style): Formatting only",
          "timestamp": "2023-05-09T15:06:04-07:00",
          "tree_id": "5f258ee3915762da41eaccbb433b3037fd47ff6e",
          "url": "https://github.com/galacticusorg/galacticus/commit/775ba693fe53c6b5354263997e74dde82922ef39"
        },
        "date": 1683683130787,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 50.782,
            "unit": "seconds",
            "range": 0.0270481052944074
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
          "id": "775ba693fe53c6b5354263997e74dde82922ef39",
          "message": "fix(style): Formatting only",
          "timestamp": "2023-05-09T15:06:04-07:00",
          "tree_id": "5f258ee3915762da41eaccbb433b3037fd47ff6e",
          "url": "https://github.com/galacticusorg/galacticus/commit/775ba693fe53c6b5354263997e74dde82922ef39"
        },
        "date": 1683683139104,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.5908987186761,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.90287496491,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25180.0342547967,
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
          "id": "775ba693fe53c6b5354263997e74dde82922ef39",
          "message": "fix(style): Formatting only",
          "timestamp": "2023-05-09T22:06:04Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/775ba693fe53c6b5354263997e74dde82922ef39"
        },
        "date": 1683906879540,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 59.496,
            "unit": "seconds",
            "range": 0.0641903419534563
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
          "id": "775ba693fe53c6b5354263997e74dde82922ef39",
          "message": "fix(style): Formatting only",
          "timestamp": "2023-05-09T22:06:04Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/775ba693fe53c6b5354263997e74dde82922ef39"
        },
        "date": 1683906888078,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.4968367295476,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.8220562018781,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24179.0559197622,
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
          "id": "b2b82e792403bd7f38c2a3097b550c1e5e6d3963",
          "message": "feat: Add more informative labels in stellar mass function analysis classes",
          "timestamp": "2023-05-12T21:09:16Z",
          "tree_id": "a23e58cfc53d109b497605c0e5efeb3d54738a99",
          "url": "https://github.com/galacticusorg/galacticus/commit/b2b82e792403bd7f38c2a3097b550c1e5e6d3963"
        },
        "date": 1683935954338,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 62.627,
            "unit": "seconds",
            "range": 0.374489118667049
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
          "id": "b2b82e792403bd7f38c2a3097b550c1e5e6d3963",
          "message": "feat: Add more informative labels in stellar mass function analysis classes",
          "timestamp": "2023-05-12T21:09:16Z",
          "tree_id": "a23e58cfc53d109b497605c0e5efeb3d54738a99",
          "url": "https://github.com/galacticusorg/galacticus/commit/b2b82e792403bd7f38c2a3097b550c1e5e6d3963"
        },
        "date": 1683935964087,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.4991777473382,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.0216488557544,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24211.0957851519,
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
          "id": "8d96a6fc15b56091328735ffbd2150ce2e65edd3",
          "message": "fix(style): Formatting only",
          "timestamp": "2023-05-17T10:14:41-07:00",
          "tree_id": "e091de5a36524d7187f56da6f98ba8b6f1e8fc7a",
          "url": "https://github.com/galacticusorg/galacticus/commit/8d96a6fc15b56091328735ffbd2150ce2e65edd3"
        },
        "date": 1684358558264,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.649,
            "unit": "seconds",
            "range": 0.0976980040737357
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
          "id": "8d96a6fc15b56091328735ffbd2150ce2e65edd3",
          "message": "fix(style): Formatting only",
          "timestamp": "2023-05-17T10:14:41-07:00",
          "tree_id": "e091de5a36524d7187f56da6f98ba8b6f1e8fc7a",
          "url": "https://github.com/galacticusorg/galacticus/commit/8d96a6fc15b56091328735ffbd2150ce2e65edd3"
        },
        "date": 1684358566003,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.4991777473382,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.0216488557544,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24211.0957851519,
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
          "id": "f5b5769ec7ef9569ff31488ea2ec20fee5bb27b5",
          "message": "fix: Initialize physically plausibility state\n\nAvoids warnings about use of uninitialized data.",
          "timestamp": "2023-05-22T16:22:29-07:00",
          "tree_id": "4d7b190aa8de61c3fba50eda1c934af25fcedad0",
          "url": "https://github.com/galacticusorg/galacticus/commit/f5b5769ec7ef9569ff31488ea2ec20fee5bb27b5"
        },
        "date": 1684807916686,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 51.91,
            "unit": "seconds",
            "range": 0.172435495185903
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
          "id": "f5b5769ec7ef9569ff31488ea2ec20fee5bb27b5",
          "message": "fix: Initialize physically plausibility state\n\nAvoids warnings about use of uninitialized data.",
          "timestamp": "2023-05-22T16:22:29-07:00",
          "tree_id": "4d7b190aa8de61c3fba50eda1c934af25fcedad0",
          "url": "https://github.com/galacticusorg/galacticus/commit/f5b5769ec7ef9569ff31488ea2ec20fee5bb27b5"
        },
        "date": 1684807925114,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.4991777473382,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.0216488557544,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24211.0957851519,
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
          "id": "4bc073b416fd78a2c5511eeeb53e37e9c28f4bc3",
          "message": "feat: Expand range of memory reporting\n\nPreviously the largest suffix supported was GB, which lead to failure to format reported memory usage on larger memory machines. Now supports up to YB.",
          "timestamp": "2023-05-23T07:56:50-07:00",
          "tree_id": "8094ff323ea9aba6616cc4d8150ae0bf20fbd868",
          "url": "https://github.com/galacticusorg/galacticus/commit/4bc073b416fd78a2c5511eeeb53e37e9c28f4bc3"
        },
        "date": 1684864489726,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 64.672,
            "unit": "seconds",
            "range": 0.388434807914992
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
          "id": "4bc073b416fd78a2c5511eeeb53e37e9c28f4bc3",
          "message": "feat: Expand range of memory reporting\n\nPreviously the largest suffix supported was GB, which lead to failure to format reported memory usage on larger memory machines. Now supports up to YB.",
          "timestamp": "2023-05-23T07:56:50-07:00",
          "tree_id": "8094ff323ea9aba6616cc4d8150ae0bf20fbd868",
          "url": "https://github.com/galacticusorg/galacticus/commit/4bc073b416fd78a2c5511eeeb53e37e9c28f4bc3"
        },
        "date": 1684864499884,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.4646519971992,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.77177881396,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25496.6894737413,
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
          "id": "9a5a5b484917ceaa6b3080ddebc0fd97d62efeb1",
          "message": "fix: Reduce code duplication in `mergerTreeEvolver` class",
          "timestamp": "2023-05-23T20:39:32Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/9a5a5b484917ceaa6b3080ddebc0fd97d62efeb1"
        },
        "date": 1684886374431,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.051,
            "unit": "seconds",
            "range": 0.0287210724025466
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
          "id": "9a5a5b484917ceaa6b3080ddebc0fd97d62efeb1",
          "message": "fix: Reduce code duplication in `mergerTreeEvolver` class",
          "timestamp": "2023-05-23T20:39:32Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/9a5a5b484917ceaa6b3080ddebc0fd97d62efeb1"
        },
        "date": 1684886385327,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.4991777473382,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.0216488557544,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 24211.0957851519,
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
          "id": "aef062b6a32bd9371c923188c298cf19e6bfd206",
          "message": "fix: Correct comment",
          "timestamp": "2023-05-24T16:42:17-07:00",
          "tree_id": "7201af723f187860e72f6dc20e96f5bfc85bb605",
          "url": "https://github.com/galacticusorg/galacticus/commit/aef062b6a32bd9371c923188c298cf19e6bfd206"
        },
        "date": 1684990705041,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 50.412,
            "unit": "seconds",
            "range": 0.0785977098900179
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
          "id": "aef062b6a32bd9371c923188c298cf19e6bfd206",
          "message": "fix: Correct comment",
          "timestamp": "2023-05-24T16:42:17-07:00",
          "tree_id": "7201af723f187860e72f6dc20e96f5bfc85bb605",
          "url": "https://github.com/galacticusorg/galacticus/commit/aef062b6a32bd9371c923188c298cf19e6bfd206"
        },
        "date": 1684990715336,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6748844868819,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9524525639199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23715.206261202,
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
          "id": "8f480a1febb414ac44d775d13442a70f918a0d59",
          "message": "feat: Allow instances of the `mergerTreeEvolveProfilerSimple` class to automatically reduce onto the object from which they were deep-copied\n\nIf they were deep-copied reduction is performed (atomically) back onto the originating object. Otherwise, results are written to file. This allows multiple deep-copied `mergerTreeEvolveProfilerSimple` objects to accumulate data independently, and to then combine those results for eventual writing to file by the ultimate originating object.",
          "timestamp": "2023-05-25T16:34:23Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/8f480a1febb414ac44d775d13442a70f918a0d59"
        },
        "date": 1685044123866,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 48.834,
            "unit": "seconds",
            "range": 0.0334424879477016
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
          "id": "8f480a1febb414ac44d775d13442a70f918a0d59",
          "message": "feat: Allow instances of the `mergerTreeEvolveProfilerSimple` class to automatically reduce onto the object from which they were deep-copied\n\nIf they were deep-copied reduction is performed (atomically) back onto the originating object. Otherwise, results are written to file. This allows multiple deep-copied `mergerTreeEvolveProfilerSimple` objects to accumulate data independently, and to then combine those results for eventual writing to file by the ultimate originating object.",
          "timestamp": "2023-05-25T16:34:23Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/8f480a1febb414ac44d775d13442a70f918a0d59"
        },
        "date": 1685044134356,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.7413507633629,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.1157669635928,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23783.3379050868,
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
          "id": "d54e213504cc20609971cf7c5b4dd81f94d085a5",
          "message": "Merge pull request #413 from galacticusorg/massPeak\n\nTrack peak bound masses of nodes",
          "timestamp": "2023-05-26T18:26:29Z",
          "tree_id": "485810bd2f7a69fc501186d65ae0ed3b13b9fd45",
          "url": "https://github.com/galacticusorg/galacticus/commit/d54e213504cc20609971cf7c5b4dd81f94d085a5"
        },
        "date": 1685137966273,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 65.432,
            "unit": "seconds",
            "range": 0.147876975894723
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
          "id": "d54e213504cc20609971cf7c5b4dd81f94d085a5",
          "message": "Merge pull request #413 from galacticusorg/massPeak\n\nTrack peak bound masses of nodes",
          "timestamp": "2023-05-26T18:26:29Z",
          "tree_id": "485810bd2f7a69fc501186d65ae0ed3b13b9fd45",
          "url": "https://github.com/galacticusorg/galacticus/commit/d54e213504cc20609971cf7c5b4dd81f94d085a5"
        },
        "date": 1685137973907,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6748844868819,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9524525639199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23715.206261202,
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
          "id": "a2d8a7e74a9d874470f8c44826897859bcfa719c",
          "message": "Merge pull request #415 from galacticusorg/threadedTrees\n\nImplement threaded evolution of individual merger trees",
          "timestamp": "2023-06-01T01:38:17Z",
          "tree_id": "504baa893d40161972eb6332900cd5bf6dfd32d9",
          "url": "https://github.com/galacticusorg/galacticus/commit/a2d8a7e74a9d874470f8c44826897859bcfa719c"
        },
        "date": 1685593562454,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.676,
            "unit": "seconds",
            "range": 0.0299065210282907
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
          "id": "a2d8a7e74a9d874470f8c44826897859bcfa719c",
          "message": "Merge pull request #415 from galacticusorg/threadedTrees\n\nImplement threaded evolution of individual merger trees",
          "timestamp": "2023-06-01T01:38:17Z",
          "tree_id": "504baa893d40161972eb6332900cd5bf6dfd32d9",
          "url": "https://github.com/galacticusorg/galacticus/commit/a2d8a7e74a9d874470f8c44826897859bcfa719c"
        },
        "date": 1685593571248,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.7589536159056,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9651779509004,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23794.3123326806,
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
          "id": "53c2c16a35fa92560d38fe549298c2d7ef0d7a97",
          "message": "fix: Remove debugging statements",
          "timestamp": "2023-06-01T07:43:33-07:00",
          "tree_id": "8c8a55f9bd6f315c8586e3e098714e642571dfc2",
          "url": "https://github.com/galacticusorg/galacticus/commit/53c2c16a35fa92560d38fe549298c2d7ef0d7a97"
        },
        "date": 1685642318557,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.64,
            "unit": "seconds",
            "range": 0.0399499687105292
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
          "id": "53c2c16a35fa92560d38fe549298c2d7ef0d7a97",
          "message": "fix: Remove debugging statements",
          "timestamp": "2023-06-01T07:43:33-07:00",
          "tree_id": "8c8a55f9bd6f315c8586e3e098714e642571dfc2",
          "url": "https://github.com/galacticusorg/galacticus/commit/53c2c16a35fa92560d38fe549298c2d7ef0d7a97"
        },
        "date": 1685642326675,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6748844868819,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9524525639199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23715.206261202,
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
          "id": "28202ef2a833d86889c91e9836a855c38571600c",
          "message": "fix: Add default values for parameters",
          "timestamp": "2023-06-02T16:15:15Z",
          "tree_id": "d18b87816399eaa6a605894bff3e7d3a3174869c",
          "url": "https://github.com/galacticusorg/galacticus/commit/28202ef2a833d86889c91e9836a855c38571600c"
        },
        "date": 1685741300467,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 58.197,
            "unit": "seconds",
            "range": 0.305689548398269
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
          "id": "28202ef2a833d86889c91e9836a855c38571600c",
          "message": "fix: Add default values for parameters",
          "timestamp": "2023-06-02T16:15:15Z",
          "tree_id": "d18b87816399eaa6a605894bff3e7d3a3174869c",
          "url": "https://github.com/galacticusorg/galacticus/commit/28202ef2a833d86889c91e9836a855c38571600c"
        },
        "date": 1685741310298,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6748844868819,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9524525639199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23715.206261202,
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
          "id": "e7a0d75118e788f65c3d12eaa620cd8e3c42398f",
          "message": "feat: Implement the Viel et al. (2005) WDM transfer function fit",
          "timestamp": "2023-06-02T22:52:05Z",
          "tree_id": "d6f0e86daf3ffcfaea1a48aabd7b8b367db8d38c",
          "url": "https://github.com/galacticusorg/galacticus/commit/e7a0d75118e788f65c3d12eaa620cd8e3c42398f"
        },
        "date": 1685763715642,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.957,
            "unit": "seconds",
            "range": 0.085088777168261
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
          "id": "e7a0d75118e788f65c3d12eaa620cd8e3c42398f",
          "message": "feat: Implement the Viel et al. (2005) WDM transfer function fit",
          "timestamp": "2023-06-02T22:52:05Z",
          "tree_id": "d6f0e86daf3ffcfaea1a48aabd7b8b367db8d38c",
          "url": "https://github.com/galacticusorg/galacticus/commit/e7a0d75118e788f65c3d12eaa620cd8e3c42398f"
        },
        "date": 1685763726616,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6748844868819,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9524525639199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23715.206261202,
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
          "id": "25acc8c52e00af2a5680e44e93c3f6c03134fd08",
          "message": "Merge pull request #418 from galacticusorg/inactiveEfficiency\n\nImprove efficiency of non-inactive evolution",
          "timestamp": "2023-06-03T04:08:31Z",
          "tree_id": "2a3f106d0d1552fe4deb83bca0cf84e64dbfa5a2",
          "url": "https://github.com/galacticusorg/galacticus/commit/25acc8c52e00af2a5680e44e93c3f6c03134fd08"
        },
        "date": 1685775872258,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 62.488,
            "unit": "seconds",
            "range": 0.22963797595338
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
          "id": "25acc8c52e00af2a5680e44e93c3f6c03134fd08",
          "message": "Merge pull request #418 from galacticusorg/inactiveEfficiency\n\nImprove efficiency of non-inactive evolution",
          "timestamp": "2023-06-03T04:08:31Z",
          "tree_id": "2a3f106d0d1552fe4deb83bca0cf84e64dbfa5a2",
          "url": "https://github.com/galacticusorg/galacticus/commit/25acc8c52e00af2a5680e44e93c3f6c03134fd08"
        },
        "date": 1685775881370,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6748844868819,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9524525639199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23715.206261202,
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
          "id": "9bc93f77147ac2c9e30a6debc1767e11b70214c3",
          "message": "fix: Correct equation, and improve model description",
          "timestamp": "2023-06-08T08:39:10-07:00",
          "tree_id": "4c3a18f7c504d6ed8c1af559994b7cf600955a92",
          "url": "https://github.com/galacticusorg/galacticus/commit/9bc93f77147ac2c9e30a6debc1767e11b70214c3"
        },
        "date": 1686250139281,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 54.617,
            "unit": "seconds",
            "range": 0.322973837949787
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
          "id": "9bc93f77147ac2c9e30a6debc1767e11b70214c3",
          "message": "fix: Correct equation, and improve model description",
          "timestamp": "2023-06-08T08:39:10-07:00",
          "tree_id": "4c3a18f7c504d6ed8c1af559994b7cf600955a92",
          "url": "https://github.com/galacticusorg/galacticus/commit/9bc93f77147ac2c9e30a6debc1767e11b70214c3"
        },
        "date": 1686250148296,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.7589536159056,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9651779509004,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23794.3123326806,
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
          "id": "d7e8d701de321cad9797a17d560893304a79bdfb",
          "message": "fix: Serialize deadlock reporting in threaded tree evolution",
          "timestamp": "2023-06-08T13:50:40-07:00",
          "tree_id": "86a1940c08b46a1edd34e0bc30a6b2c7f898f1cb",
          "url": "https://github.com/galacticusorg/galacticus/commit/d7e8d701de321cad9797a17d560893304a79bdfb"
        },
        "date": 1686268513549,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 48.587,
            "unit": "seconds",
            "range": 0.0224521713859421
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
          "id": "d7e8d701de321cad9797a17d560893304a79bdfb",
          "message": "fix: Serialize deadlock reporting in threaded tree evolution",
          "timestamp": "2023-06-08T13:50:40-07:00",
          "tree_id": "86a1940c08b46a1edd34e0bc30a6b2c7f898f1cb",
          "url": "https://github.com/galacticusorg/galacticus/commit/d7e8d701de321cad9797a17d560893304a79bdfb"
        },
        "date": 1686268523393,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6748844868819,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9524525639199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23715.206261202,
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
          "id": "2978943e5d2eaa7b5a1f5b2031b1e280f72fec51",
          "message": "fix: Correct spelling of word",
          "timestamp": "2023-06-14T07:28:02-07:00",
          "tree_id": "a7222a0d7dda634b72d942c04787ac0c2e7015a9",
          "url": "https://github.com/galacticusorg/galacticus/commit/2978943e5d2eaa7b5a1f5b2031b1e280f72fec51"
        },
        "date": 1686763464098,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 57.514,
            "unit": "seconds",
            "range": 0.107565793818906
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
          "id": "2978943e5d2eaa7b5a1f5b2031b1e280f72fec51",
          "message": "fix: Correct spelling of word",
          "timestamp": "2023-06-14T07:28:02-07:00",
          "tree_id": "a7222a0d7dda634b72d942c04787ac0c2e7015a9",
          "url": "https://github.com/galacticusorg/galacticus/commit/2978943e5d2eaa7b5a1f5b2031b1e280f72fec51"
        },
        "date": 1686763474032,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6748844868819,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9524525639199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23715.206261202,
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
          "id": "5f5475eb504a185be42ac2e99153f9cdfca22d9b",
          "message": "fix: Correct equation in parameter definition",
          "timestamp": "2023-06-20T08:07:45-07:00",
          "tree_id": "f77da09d33f572b107ddd9451a109be09c43d388",
          "url": "https://github.com/galacticusorg/galacticus/commit/5f5475eb504a185be42ac2e99153f9cdfca22d9b"
        },
        "date": 1687288633451,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.942,
            "unit": "seconds",
            "range": 0.0246495436056456
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
          "id": "5f5475eb504a185be42ac2e99153f9cdfca22d9b",
          "message": "fix: Correct equation in parameter definition",
          "timestamp": "2023-06-20T08:07:45-07:00",
          "tree_id": "f77da09d33f572b107ddd9451a109be09c43d388",
          "url": "https://github.com/galacticusorg/galacticus/commit/5f5475eb504a185be42ac2e99153f9cdfca22d9b"
        },
        "date": 1687288642254,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6748844868819,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9524525639199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23715.206261202,
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
          "id": "50c07f709063ddfe99b7fd3b569745b0e29c3010",
          "message": "feat: Add a `galacticFilter` that passes only nodes on the constrained branch",
          "timestamp": "2023-06-21T15:52:30Z",
          "tree_id": "9315f145782b85c53b3e476562ad4e173eac805b",
          "url": "https://github.com/galacticusorg/galacticus/commit/50c07f709063ddfe99b7fd3b569745b0e29c3010"
        },
        "date": 1687382843805,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 63.845,
            "unit": "seconds",
            "range": 0.149219636777396
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
          "id": "50c07f709063ddfe99b7fd3b569745b0e29c3010",
          "message": "feat: Add a `galacticFilter` that passes only nodes on the constrained branch",
          "timestamp": "2023-06-21T15:52:30Z",
          "tree_id": "9315f145782b85c53b3e476562ad4e173eac805b",
          "url": "https://github.com/galacticusorg/galacticus/commit/50c07f709063ddfe99b7fd3b569745b0e29c3010"
        },
        "date": 1687382851715,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6748844868819,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9524525639199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23715.206261202,
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
          "id": "47ab81ae3e9a4086476941ec9dced579a87fd5a1",
          "message": "Merge pull request #420 from galacticusorg/satelliteDestructionHelp\n\nDetect failure to destroy satellite halos",
          "timestamp": "2023-06-21T23:31:05Z",
          "tree_id": "0cd8f196cc3fdd8ab20782df8a7526eb4191599a",
          "url": "https://github.com/galacticusorg/galacticus/commit/47ab81ae3e9a4086476941ec9dced579a87fd5a1"
        },
        "date": 1687401343623,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.707,
            "unit": "seconds",
            "range": 0.0168552662393373
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
          "id": "47ab81ae3e9a4086476941ec9dced579a87fd5a1",
          "message": "Merge pull request #420 from galacticusorg/satelliteDestructionHelp\n\nDetect failure to destroy satellite halos",
          "timestamp": "2023-06-21T23:31:05Z",
          "tree_id": "0cd8f196cc3fdd8ab20782df8a7526eb4191599a",
          "url": "https://github.com/galacticusorg/galacticus/commit/47ab81ae3e9a4086476941ec9dced579a87fd5a1"
        },
        "date": 1687401351715,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6556202706915,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.0940486680264,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23691.6828998057,
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
          "distinct": false,
          "id": "7d940aeeeaf9329a8f1906a02e3dfcfa95caf915",
          "message": "Merge pull request #426 from galacticusorg/metaPropertyComponentFix\n\nMake meta-property set/get functions respect the active componens list",
          "timestamp": "2023-06-23T16:40:35Z",
          "tree_id": "98b2947cebf7e0015c188941db9c3fb899e192ec",
          "url": "https://github.com/galacticusorg/galacticus/commit/7d940aeeeaf9329a8f1906a02e3dfcfa95caf915"
        },
        "date": 1687556871482,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.793,
            "unit": "seconds",
            "range": 0.0415222831753291
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
          "distinct": false,
          "id": "7d940aeeeaf9329a8f1906a02e3dfcfa95caf915",
          "message": "Merge pull request #426 from galacticusorg/metaPropertyComponentFix\n\nMake meta-property set/get functions respect the active componens list",
          "timestamp": "2023-06-23T16:40:35Z",
          "tree_id": "98b2947cebf7e0015c188941db9c3fb899e192ec",
          "url": "https://github.com/galacticusorg/galacticus/commit/7d940aeeeaf9329a8f1906a02e3dfcfa95caf915"
        },
        "date": 1687556880094,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.7589536159056,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9651779509004,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23794.3123326806,
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
          "id": "dc1bd8c1bf667c62e3dd4d1718caf77f48a67b72",
          "message": "fix: Ensure that checkpointing does not interfere with running multiple trees",
          "timestamp": "2023-06-29T00:02:29Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/dc1bd8c1bf667c62e3dd4d1718caf77f48a67b72"
        },
        "date": 1688014417079,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.718,
            "unit": "seconds",
            "range": 0.0196367003333407
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
          "id": "dc1bd8c1bf667c62e3dd4d1718caf77f48a67b72",
          "message": "fix: Ensure that checkpointing does not interfere with running multiple trees",
          "timestamp": "2023-06-29T00:02:29Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/dc1bd8c1bf667c62e3dd4d1718caf77f48a67b72"
        },
        "date": 1688014426864,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6748844868819,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9524525639199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23715.206261202,
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
          "distinct": false,
          "id": "d726a1c8686d37595881f2e118acb01a85942103",
          "message": "Merge pull request #427 from galacticusorg/desrializationFix\n\nDeserialization fix",
          "timestamp": "2023-06-29T14:29:07Z",
          "tree_id": "cb156267cb8650852e0bc944ea2e0187182799b0",
          "url": "https://github.com/galacticusorg/galacticus/commit/d726a1c8686d37595881f2e118acb01a85942103"
        },
        "date": 1688608357045,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.85,
            "unit": "seconds",
            "range": 0.0224499443206618
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
          "distinct": false,
          "id": "d726a1c8686d37595881f2e118acb01a85942103",
          "message": "Merge pull request #427 from galacticusorg/desrializationFix\n\nDeserialization fix",
          "timestamp": "2023-06-29T14:29:07Z",
          "tree_id": "cb156267cb8650852e0bc944ea2e0187182799b0",
          "url": "https://github.com/galacticusorg/galacticus/commit/d726a1c8686d37595881f2e118acb01a85942103"
        },
        "date": 1688608366365,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6748844868819,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9524525639199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23715.206261202,
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
          "distinct": false,
          "id": "212e634c8a6c24f5dc09ae49ad4c9004fcc76053",
          "message": "Merge pull request #428 from galacticusorg/checkpointing\n\nImplement (limited) checkpointing",
          "timestamp": "2023-07-07T15:58:53Z",
          "tree_id": "fe32efe68225e92bc677a4091a8f1b2cdabdea51",
          "url": "https://github.com/galacticusorg/galacticus/commit/212e634c8a6c24f5dc09ae49ad4c9004fcc76053"
        },
        "date": 1688772299065,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 62.773,
            "unit": "seconds",
            "range": 0.305388441169516
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
          "distinct": false,
          "id": "212e634c8a6c24f5dc09ae49ad4c9004fcc76053",
          "message": "Merge pull request #428 from galacticusorg/checkpointing\n\nImplement (limited) checkpointing",
          "timestamp": "2023-07-07T15:58:53Z",
          "tree_id": "fe32efe68225e92bc677a4091a8f1b2cdabdea51",
          "url": "https://github.com/galacticusorg/galacticus/commit/212e634c8a6c24f5dc09ae49ad4c9004fcc76053"
        },
        "date": 1688772308587,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6556202706915,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.0940486680264,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23691.6828998057,
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
          "distinct": false,
          "id": "c7f57b8bb906b62bbb314f01962b1be69cd65ea1",
          "message": "Merge pull request #430 from galacticusorg/nodeLabels\n\nImplement adding arbitrary labels to nodes",
          "timestamp": "2023-07-08T05:07:56Z",
          "tree_id": "6b3a98ab1ece20368604d970019b7965ed4fa606",
          "url": "https://github.com/galacticusorg/galacticus/commit/c7f57b8bb906b62bbb314f01962b1be69cd65ea1"
        },
        "date": 1688803547106,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 48.616,
            "unit": "seconds",
            "range": 0.541516758743409
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
          "distinct": false,
          "id": "c7f57b8bb906b62bbb314f01962b1be69cd65ea1",
          "message": "Merge pull request #430 from galacticusorg/nodeLabels\n\nImplement adding arbitrary labels to nodes",
          "timestamp": "2023-07-08T05:07:56Z",
          "tree_id": "6b3a98ab1ece20368604d970019b7965ed4fa606",
          "url": "https://github.com/galacticusorg/galacticus/commit/c7f57b8bb906b62bbb314f01962b1be69cd65ea1"
        },
        "date": 1688803555097,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6748844868819,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9524525639199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23715.206261202,
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
          "id": "f3c02aa2df60000d3b8f0240deb27124e8d9323a",
          "message": "fix: Update the `Merge_Models.pl` script to handle multi-dimensional datasets",
          "timestamp": "2023-07-10T09:15:26-07:00",
          "tree_id": "f5412ce4a4a2d468bc4e6f6c7d46a581d9b13a24",
          "url": "https://github.com/galacticusorg/galacticus/commit/f3c02aa2df60000d3b8f0240deb27124e8d9323a"
        },
        "date": 1689017066137,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 48.845,
            "unit": "seconds",
            "range": 0.0285043856271913
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
          "id": "f3c02aa2df60000d3b8f0240deb27124e8d9323a",
          "message": "fix: Update the `Merge_Models.pl` script to handle multi-dimensional datasets",
          "timestamp": "2023-07-10T09:15:26-07:00",
          "tree_id": "f5412ce4a4a2d468bc4e6f6c7d46a581d9b13a24",
          "url": "https://github.com/galacticusorg/galacticus/commit/f3c02aa2df60000d3b8f0240deb27124e8d9323a"
        },
        "date": 1689017074532,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 57.6748844868819,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9524525639199,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 23715.206261202,
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
          "distinct": false,
          "id": "82b8b6cff4bd312afea560b08232aca6148a1177",
          "message": "Merge pull request #429 from galacticusorg/parallelTreeBuild\n\nImplement parallel tree builds",
          "timestamp": "2023-07-11T14:37:11Z",
          "tree_id": "53965748a6f7c2e2b05f7df5b9889971900049a1",
          "url": "https://github.com/galacticusorg/galacticus/commit/82b8b6cff4bd312afea560b08232aca6148a1177"
        },
        "date": 1689106719548,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 49.343,
            "unit": "seconds",
            "range": 0.0281797799846736
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
          "distinct": false,
          "id": "82b8b6cff4bd312afea560b08232aca6148a1177",
          "message": "Merge pull request #429 from galacticusorg/parallelTreeBuild\n\nImplement parallel tree builds",
          "timestamp": "2023-07-11T14:37:11Z",
          "tree_id": "53965748a6f7c2e2b05f7df5b9889971900049a1",
          "url": "https://github.com/galacticusorg/galacticus/commit/82b8b6cff4bd312afea560b08232aca6148a1177"
        },
        "date": 1689106729298,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 58.5242080709765,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.4439006068673,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25257.8908941139,
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
          "id": "286b169fcc84637043f8092810cafef367fbbf09",
          "message": "fix: Add missing `close()` of HDF5 file",
          "timestamp": "2023-07-11T16:28:24-07:00",
          "tree_id": "7c0d242e9743bf848fbada41e92e8866e57a34f6",
          "url": "https://github.com/galacticusorg/galacticus/commit/286b169fcc84637043f8092810cafef367fbbf09"
        },
        "date": 1689129314745,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 56.017,
            "unit": "seconds",
            "range": 0.0405475030043149
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
          "id": "286b169fcc84637043f8092810cafef367fbbf09",
          "message": "fix: Add missing `close()` of HDF5 file",
          "timestamp": "2023-07-11T16:28:24-07:00",
          "tree_id": "7c0d242e9743bf848fbada41e92e8866e57a34f6",
          "url": "https://github.com/galacticusorg/galacticus/commit/286b169fcc84637043f8092810cafef367fbbf09"
        },
        "date": 1689129322561,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 58.5486581678349,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.4253490900446,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25236.854696677,
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
          "id": "51fa2531801a56b2903d10a89173cd55c304f020",
          "message": "fix: Avoid attempting to access a non-existant star formation history",
          "timestamp": "2023-07-12T15:09:49Z",
          "tree_id": "f75659c7441a8a75cc4504c7afe5aa43bd38cb46",
          "url": "https://github.com/galacticusorg/galacticus/commit/51fa2531801a56b2903d10a89173cd55c304f020"
        },
        "date": 1689189091892,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 57.827,
            "unit": "seconds",
            "range": 0.0628816348392391
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
          "id": "51fa2531801a56b2903d10a89173cd55c304f020",
          "message": "fix: Avoid attempting to access a non-existant star formation history",
          "timestamp": "2023-07-12T15:09:49Z",
          "tree_id": "f75659c7441a8a75cc4504c7afe5aa43bd38cb46",
          "url": "https://github.com/galacticusorg/galacticus/commit/51fa2531801a56b2903d10a89173cd55c304f020"
        },
        "date": 1689189101011,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 58.5231439321362,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 26.3590095525618,
            "unit": "-logℒ"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 25172.6986688112,
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
          "distinct": false,
          "id": "abbd5fd84eeb1f6ab3f6e3a790125849d7c05548",
          "message": "Merge pull request #435 from galacticusorg/labelingImprovements\n\nAdd further node labeling functionality",
          "timestamp": "2023-07-13T02:14:50Z",
          "tree_id": "a0ffa0dbe4b437ee455f25e8232ef280a5262dcf",
          "url": "https://github.com/galacticusorg/galacticus/commit/abbd5fd84eeb1f6ab3f6e3a790125849d7c05548"
        },
        "date": 1689226027712,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 57.315,
            "unit": "seconds",
            "range": 0.0194036079140963
          }
        ]
      }
    ]
  }
}