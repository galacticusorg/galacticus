window.BENCHMARK_DATA = {
  "lastUpdate": 1689270190813,
  "repoUrl": "https://github.com/galacticusorg/galacticus",
  "entries": {
    "Milky Way model benchmarks": [
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
        "date": 1668539630689,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.31979942141629,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.63147644711799,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.759755093109701,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7020542554721,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 290.888298717031,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 11.1286085085407,
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
        "date": 1668557713617,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 173.963,
            "unit": "seconds",
            "range": 0.170528883181682
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
        "date": 1668557722384,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.67490424988166,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.8291148261784,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.740264346117683,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.8108941711284,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 289.797926519564,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.38782096420975,
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
        "date": 1668676847957,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 187.108,
            "unit": "seconds",
            "range": 0.834208367255623
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
        "date": 1668676855353,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.42314718457122,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.78763440613205,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.717111170458921,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6411532833414,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 288.997452249363,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.46636447671286,
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
        "date": 1668709389159,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 213.298,
            "unit": "seconds",
            "range": 0.294899304847486
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
        "date": 1668709399042,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.43556443436774,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.1133219315734,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.980139986280522,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.21673702163,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 68.1164520447607,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.41545803532447,
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
        "date": 1668806882102,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 213.459,
            "unit": "seconds",
            "range": 0.201779334917814
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
        "date": 1668806890084,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.97204773396112,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.23891487131586,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.967930654500715,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 22.9941561260547,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 66.4322825794555,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.30543912211942,
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
        "date": 1669721039194,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 185.571,
            "unit": "seconds",
            "range": 0.30274229965706
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
        "date": 1669721047849,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.45982964385489,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.64974015742461,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.794502076751787,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6686716744117,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 65.2960772707066,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 11.0722562418528,
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
        "date": 1669876303953,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 262.986,
            "unit": "seconds",
            "range": 0.204108794518558
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
        "date": 1669876311958,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.98905712606957,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.32233700724143,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.05709405260061,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.1123147651418,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 286.992446362324,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.71239664232992,
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
        "date": 1669969134268,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 216.312,
            "unit": "seconds",
            "range": 0.541636040160487
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
        "date": 1669969141605,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.07151757285394,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.9198625142718,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.979868451256468,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.612207440699,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 69.272179966433,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.82672374082168,
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
        "date": 1670032548275,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 299.663,
            "unit": "seconds",
            "range": 0.467103949892123
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
        "date": 1670032556911,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.25183710138655,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.73810414639679,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.778343031788948,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.4784672434038,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 290.409659817015,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 11.2705584376751,
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
        "date": 1670322844624,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 183.63,
            "unit": "seconds",
            "range": 0.189425447074571
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
        "date": 1670322852958,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.53131197535584,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.27879851439277,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.756982573778964,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6855158302598,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 307.144185081575,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 3.14177759939532,
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
        "date": 1670444818127,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 194.267,
            "unit": "seconds",
            "range": 0.166961372780854
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
        "date": 1670444827235,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.67432845241746,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.07607313793598,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.788381258293099,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6983232343028,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 304.955173726459,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 3.0218737042497,
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
        "date": 1670521161664,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 265.072,
            "unit": "seconds",
            "range": 0.310653504728486
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
        "date": 1670521169793,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.71646062086302,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.0372700792549,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.974906688906732,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.450409249967,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 65.2342861452598,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.26426554522604,
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
        "date": 1670903240890,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 267.318,
            "unit": "seconds",
            "range": 0.346940340695513
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
        "date": 1670903249110,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.74719607234151,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 12.9230317019417,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.0735962547226,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 22.9401142163395,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 66.7953818935949,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 11.402561612827,
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
        "date": 1671029393363,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 257.554,
            "unit": "seconds",
            "range": 1.1553442776949
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
        "date": 1671029401532,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.88122664661788,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.17394958701615,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.861909274185707,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6742307131901,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 37.9605757332934,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 2.31616085287071,
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
        "date": 1671294885168,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 217.037,
            "unit": "seconds",
            "range": 0.196041067129384
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
        "date": 1671294891619,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.73491060550965,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 13.1763608475997,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.05758713178247,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 22.9150597842577,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 70.7268487670879,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 11.714061735364,
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
        "date": 1671565950630,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 192.941,
            "unit": "seconds",
            "range": 0.234505650255011
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
        "date": 1671565960639,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.6708451777911,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 13.44060374184,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.939897471452238,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.0189648016844,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 70.0203154760673,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 11.484144632184,
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
        "date": 1671742151131,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 217.835,
            "unit": "seconds",
            "range": 0.185344274262362
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
        "date": 1671742158857,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.96879569185002,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.77599062006793,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.719694732850795,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.8697963230676,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 75.8500959995063,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 2.45160734605733,
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
        "date": 1672009612727,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 186.258,
            "unit": "seconds",
            "range": 0.0965380753870462
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
        "date": 1672009620991,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.65270505878469,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 12.9545014982176,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.09516852807782,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 22.7378529753071,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 70.1728557758844,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 11.5015497850488,
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
        "date": 1672080013272,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 297.831,
            "unit": "seconds",
            "range": 0.514644440365331
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
        "date": 1672080021117,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.07714205044172,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.63624824485506,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.697736604654721,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0358685185567,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 40.473740760251,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 2.49187741591571,
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
        "date": 1673314914362,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 202.06,
            "unit": "seconds",
            "range": 2.25866819165616
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
        "date": 1673314922920,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.81072221097032,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.45076725456392,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.734678542731142,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7987908290354,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 68.9324830990023,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 2.3714606398896,
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
        "date": 1673412840108,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 168.512,
            "unit": "seconds",
            "range": 0.22697929420991
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
        "date": 1673412853808,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.0185649676824,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.42705059826988,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.710829329351966,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.618161115823,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 70.3174544128157,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 2.30808952113768,
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
        "date": 1673461428335,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 271.572,
            "unit": "seconds",
            "range": 0.24263058339916
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
        "date": 1673461437276,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.05164746874204,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.17391935197404,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.976019385325774,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.1207057807193,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 396.484873053692,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.38944120305051,
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
        "date": 1673575186612,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 289.176,
            "unit": "seconds",
            "range": 0.773898184517248
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
        "date": 1673575195790,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.95939795523904,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.4534409106027,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.88541450748119,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.3893277268888,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 31.6309966007773,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.31448769550719,
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
        "date": 1673651810627,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 225.519,
            "unit": "seconds",
            "range": 0.257388616685902
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
        "date": 1673651819063,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.76590512125523,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.41322430000064,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.836541445760069,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6273527469705,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 397.987190575674,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 2.5353450699136,
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
        "date": 1673668329016,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 224.453,
            "unit": "seconds",
            "range": 0.430441749833079
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
        "date": 1673668337367,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.22816408219212,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.91057149080348,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.35814741301699,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3211646851389,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 127.436197207536,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.70563247666536,
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
        "date": 1673829658667,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 257.527,
            "unit": "seconds",
            "range": 0.777258065253703
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
        "date": 1673829668397,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.38502391741538,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.66379844371419,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.40237170633146,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3021192290636,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 125.673926080859,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.68252811180385,
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
        "date": 1674044461666,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 237.975,
            "unit": "seconds",
            "range": 0.345427995391691
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
        "date": 1674044470340,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.19047688207088,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.64981516689162,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.32655354050141,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3953838429523,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 128.672768871945,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.65090468033682,
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
        "date": 1674070328209,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 210.149,
            "unit": "seconds",
            "range": 0.59596048526604
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
        "date": 1674070337704,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.23815187816524,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.40904067034189,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.29172590208447,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1826493284013,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 127.463911655394,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.43541678924816,
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
        "date": 1674114134923,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 215.806,
            "unit": "seconds",
            "range": 0.199535460509189
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
        "date": 1674114143277,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.1169506865339,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.75832153735309,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.42042735399537,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3751965182674,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 128.370883852034,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.62063195740025,
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
        "date": 1674268924395,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 193.588,
            "unit": "seconds",
            "range": 1.68454195554745
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
        "date": 1674268932270,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.10506800245251,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.53819607330649,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.33404255369956,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3429843749987,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 127.383069346823,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.4595580177165,
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
        "date": 1674294083794,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 169.036,
            "unit": "seconds",
            "range": 0.203018225782213
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
        "date": 1674294091290,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.19825214314605,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.37455923595765,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.3415896322779,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2914008947745,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 128.378002000386,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.57542206221595,
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
        "date": 1674328011349,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 198.737,
            "unit": "seconds",
            "range": 0.0393967003617005
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
        "date": 1674328019012,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.24401521702511,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.38636092377912,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.2229992582984,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2738358273315,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 166.494354419332,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.37805817800634,
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
        "date": 1674506639348,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 215.888,
            "unit": "seconds",
            "range": 0.140675513151767
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
        "date": 1674506647236,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.47407785599028,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 3.70240379789367,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.673723765898714,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2081458375965,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 59.984565525506,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 7.9905564307121,
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
        "date": 1674764757240,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 169.465,
            "unit": "seconds",
            "range": 0.190511154529273
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
        "date": 1674764765090,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.94777572825317,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.15675210202248,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.36315484000958,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0777624600845,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 143.261797500755,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.8679542626897,
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
        "date": 1674799579401,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 218.09,
            "unit": "seconds",
            "range": 0.160449368959167
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
        "date": 1674799587199,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.19998766468113,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.1911017486927,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.33243382789192,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.5523708949112,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 143.059749716486,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.71184258662052,
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
        "date": 1674847882541,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 222.415,
            "unit": "seconds",
            "range": 0.42409256065045
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
        "date": 1674847889779,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.38581413657457,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.09707388600789,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.27948101200236,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3517746574536,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 141.562867040534,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.59796194932996,
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
        "date": 1674882109370,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 258.598,
            "unit": "seconds",
            "range": 0.192181164519245
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
        "date": 1674882116326,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.24129233581603,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.00661558569096,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.26477982147318,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.5991915543534,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 142.579555078319,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.62058430371848,
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
        "date": 1674897006415,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 218.419,
            "unit": "seconds",
            "range": 0.0669096405610489
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
        "date": 1674897013901,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.20126175766276,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.52724069668987,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.3476075689826,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.5621452342043,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 141.687205717416,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.70350721329553,
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
        "date": 1675210338495,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 236.384,
            "unit": "seconds",
            "range": 0.140756527377953
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
        "date": 1675210345292,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.24743796106887,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.52271281983246,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.24029118341503,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.6869107720592,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 127.583723687279,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.51293095561443,
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
        "date": 1675306357575,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 210.372,
            "unit": "seconds",
            "range": 0.174938846460426
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
        "date": 1675306365743,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.23069962892554,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.62348085635882,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.29905942249526,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2433559609089,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 128.118786327513,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.18323187399699,
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
        "date": 1675834180918,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 203.646,
            "unit": "seconds",
            "range": 0.141938014642032
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
        "date": 1675834188891,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.24129233581603,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.00661558569096,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.26477982147318,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.5991915543534,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 142.579555078319,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.62058430371848,
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
        "date": 1675929128811,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 193.609,
            "unit": "seconds",
            "range": 0.188156583726062
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1675929137285,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.20680353037085,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.91802523710661,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.45331293497713,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9174768346761,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 140.934767261758,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.74975397136834,
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
        "date": 1675965856857,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 205.476,
            "unit": "seconds",
            "range": 0.110355788249943
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
        "date": 1675965867388,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.99354009002036,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.06357475544158,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.13330400835217,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.057374866708,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 64.9092665000299,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.64630488066786,
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
        "date": 1676016282777,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 166.563,
            "unit": "seconds",
            "range": 0.121268709892022
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1676016292526,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.99354009002036,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.06357475544158,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.13330400835217,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.057374866708,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 64.9092665000299,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.64630488066786,
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
        "date": 1676076008966,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 210.788,
            "unit": "seconds",
            "range": 0.111801610005691
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1676076018354,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.439903738332,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.78521774416878,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.35573635014961,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7431041256945,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 142.535108592588,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.72434381474404,
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
        "date": 1676105773024,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 198.352,
            "unit": "seconds",
            "range": 0.211280855734581
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1676105780556,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.23943917652469,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.73389638931614,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.40507615274271,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6881082203394,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 141.339118751903,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.74047832975763,
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
        "date": 1676524755280,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 191.988,
            "unit": "seconds",
            "range": 0.165213800876532
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1676524762854,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.43238662126871,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.79146054714478,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.42857131484659,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7345886615783,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 164.233717568923,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.4778507490274,
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
        "date": 1676616347186,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 229.258,
            "unit": "seconds",
            "range": 1.94321270065846
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1676616354069,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.30718783246319,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.16304336879793,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.40778345557581,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7562426697074,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 81.3855201915745,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.72276313987552,
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
        "date": 1676685537151,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 255.021,
            "unit": "seconds",
            "range": 2.06548999029326
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1676685545238,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.43603072557468,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.5330987601578,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.49206863496866,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.5144410208183,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 164.950608366688,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 7.07072873003251,
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
        "date": 1676836670905,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 240.287,
            "unit": "seconds",
            "range": 0.158883919894827
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1676836678131,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.24356588017375,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.71419848134857,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.40433805044331,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7474627244234,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 64.2713383008829,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.44883822229689,
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
        "date": 1677089043433,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 265.328,
            "unit": "seconds",
            "range": 0.282629793195501
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
        "date": 1677089052635,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.36617009483029,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.40486938438531,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.50462094988436,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.562177435062,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 65.650894084706,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.49685254772603,
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
        "date": 1677171082890,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 177.423,
            "unit": "seconds",
            "range": 0.212899271955418
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1677171092921,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.26735690234578,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.94853331071863,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.37654105701467,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.573277136367,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 142.893552185239,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.67828658728243,
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
        "date": 1677359083831,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 216.917,
            "unit": "seconds",
            "range": 0.202524319527258
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1677359091043,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.34868638993938,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.54221189534925,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.51016345109794,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.8885373102773,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 142.127144267418,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.67537446498706,
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
        "date": 1677524650352,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 238.797,
            "unit": "seconds",
            "range": 0.144665476189165
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
        "date": 1677524658060,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.29748740788134,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.91451257220672,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.64900437061172,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3952237775667,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 188.072024728008,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.1990526901257,
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
        "date": 1677567486982,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 283.452,
            "unit": "seconds",
            "range": 2.49529549352372
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
        "date": 1677567495756,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.34289205344315,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.09183775402284,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.41317954554462,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7308317220014,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 141.604775784843,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.73650795739818,
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
        "date": 1677656067869,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 289.239,
            "unit": "seconds",
            "range": 1.23868837889082
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
        "date": 1677656078300,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 242.639,
            "unit": "seconds",
            "range": 0.0955034030874793
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
        "date": 1677656087264,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.2868741377072,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.18677898060952,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.4353720002123,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7690786865207,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 143.338305672129,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.4578763192726,
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
        "date": 1677684042331,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 289.239,
            "unit": "seconds",
            "range": 1.23868837889082
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1677684051458,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.00004825328433,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.9050419106813,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.25728230506136,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1566038899932,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 142.984488336172,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.47215884706706,
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
        "date": 1677709728282,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 231.255,
            "unit": "seconds",
            "range": 0.199184587753079
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
        "date": 1677709737378,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.14126877981415,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.741638017373,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.4279200033209,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.5245949673294,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 142.171869633297,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.44404786307016,
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
        "date": 1677740574875,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 219.831,
            "unit": "seconds",
            "range": 0.294507894631044
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
        "date": 1677740583197,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.47364025092844,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.53758268871479,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.53218219851468,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.789260525407,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 64.2621357464532,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.73819027400434,
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
        "date": 1677781098929,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 175.017,
            "unit": "seconds",
            "range": 0.215657367134727
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1677781108997,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.26505273631688,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.37018318943062,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.40413434301698,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.4516442234532,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 141.956766478276,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.44431744185951,
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
        "date": 1677825880416,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 217.843,
            "unit": "seconds",
            "range": 0.192733235324866
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
        "date": 1677825890024,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.39057856495757,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.99290907676918,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.702771518012212,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2363599160367,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 317.861696751683,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.2914709476475,
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
        "date": 1678168155166,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 251.982,
            "unit": "seconds",
            "range": 0.529716527964554
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1678168161923,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.27297114387412,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.64501953034879,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.39004551853678,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7786895873609,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 67.1071494168893,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.56338874010157,
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
        "date": 1678218642863,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 249.897,
            "unit": "seconds",
            "range": 0.6016644413621
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
        "date": 1678218652300,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.24891831268653,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.80243197696866,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.838163012761256,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.608573425344,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 188.705020874503,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 7.86472142815646,
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
        "date": 1678399907491,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 176.052,
            "unit": "seconds",
            "range": 0.226007964462746
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1678399914762,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.34868638993938,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.54221189534925,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.51016345109794,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.8885373102773,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 142.127144267418,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.67537446498706,
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
        "date": 1678483321713,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 213.056,
            "unit": "seconds",
            "range": 0.122304537930657
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
        "date": 1678483329692,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.33483125507875,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.61978407452969,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.32133941547381,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.8470219564497,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 142.657639774832,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.71336170415854,
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
        "date": 1678527952587,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 296.383,
            "unit": "seconds",
            "range": 0.213274705490401
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
        "date": 1678527961601,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.35332041192434,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.5436592434515,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.40028004503177,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7432096851195,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 141.286764233104,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.80107183733161,
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
        "date": 1678827302057,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 312.955,
            "unit": "seconds",
            "range": 0.247484342937196
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
        "date": 1678827312449,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.20123519434001,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.45442870928173,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.41613121166454,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6438603673934,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 143.058448144314,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.50329009495148,
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
        "date": 1679112493357,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 194.874,
            "unit": "seconds",
            "range": 0.0769051363781869
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
        "date": 1679112501165,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.28110570570714,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.94626779330773,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.41075286456742,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6289779685293,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 64.6494200268977,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.46178080362455,
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
        "date": 1679181525402,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 291.915,
            "unit": "seconds",
            "range": 0.0961587229654007
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
        "date": 1679181534227,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.44306329035049,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.73841147144231,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.37600041972057,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6305453963666,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 140.873929119564,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.46573297111352,
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
        "date": 1679418452932,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 290.423,
            "unit": "seconds",
            "range": 0.374916123946998
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
        "date": 1679418460697,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.26779796554176,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.87606235212848,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.42197662955673,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6866793752517,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 141.168654015089,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.46257742770265,
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
        "date": 1679438177848,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 267.23,
            "unit": "seconds",
            "range": 0.121762063046193
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
        "date": 1679438187002,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.44306329035049,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.73841147144231,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.37600041972057,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6305453963666,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 140.873929119564,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.46573297111352,
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
        "date": 1679448285504,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 164.637,
            "unit": "seconds",
            "range": 0.383497196859046
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
        "date": 1679448293121,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.40383504329283,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.9489584060363,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.36946233156206,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.5595817769696,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 64.6921484520708,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.4383754375548,
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
        "date": 1679531806149,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 163.27,
            "unit": "seconds",
            "range": 0.0340881211104489
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
        "date": 1679531814116,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.44306329035049,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.73841147144231,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.37600041972057,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6305453963666,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 140.873929119564,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.46573297111352,
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
        "date": 1679619746928,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 222.957,
            "unit": "seconds",
            "range": 0.212250088337337
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
        "date": 1679619754844,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.4065645380475,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 11.0008878438288,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.44038215417029,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.8697639479956,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 64.8876428775087,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.8083511673837,
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
        "date": 1679628779334,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 214.689,
            "unit": "seconds",
            "range": 0.235080624467007
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
        "date": 1679628787282,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.17075252022424,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.86441096319093,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.46615875509521,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6205876828087,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 143.167759718882,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.1246291825675,
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
        "date": 1680216441954,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 205.8,
            "unit": "seconds",
            "range": 0.10509995242446
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
        "date": 1680216450143,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.42226129623765,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.23794918619196,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.34938858597536,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.602285406403,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 64.3549855937196,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.50721401915105,
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
        "date": 1680250301033,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 226.093,
            "unit": "seconds",
            "range": 0.30597075023678
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1680250311083,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.1604698102024,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.18049603499992,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.39540881878449,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.5809491844137,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 142.335278692501,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.74043316890731,
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
        "date": 1680268659794,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 225.277,
            "unit": "seconds",
            "range": 0.220613009586261
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
        "date": 1680268671557,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.27649920589013,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.60362937640474,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.43339231678352,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7539764404336,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 65.7228956730472,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.54274301479992,
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
        "date": 1680584665446,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 204.771,
            "unit": "seconds",
            "range": 0.343946071356451
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
        "date": 1680584673001,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.51724399426135,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.20989929280312,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.42914734439032,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6790135034038,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 61.3407273049249,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.71313261005903,
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
        "date": 1680641650021,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 208.294,
            "unit": "seconds",
            "range": 0.354711713930895
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
        "date": 1680641657928,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.37258905697622,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.30645050912398,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.43019387446046,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6241023203276,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 65.682013747538,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.69120733767007,
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
        "date": 1680825354005,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 374.089,
            "unit": "seconds",
            "range": 0.162323442556393
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
        "date": 1680825362360,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.79157256013101,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.9702721992907,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.49676212097557,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6964829504736,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 144.294084063239,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.66446053023421,
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
        "date": 1681431720560,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 205.923,
            "unit": "seconds",
            "range": 0.155968265999435
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
        "date": 1681431729792,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.37258905697622,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.30645050912398,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.43019387446882,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6241023203276,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 65.682013747538,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.69120733767007,
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
        "date": 1681699213475,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 207.323,
            "unit": "seconds",
            "range": 0.208360504892514
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
        "date": 1681699223550,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.1879854460878,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.28827507651342,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.41812133815718,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7889630707262,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 143.690697978052,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.6820950873926,
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
        "date": 1681718297669,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 146.768,
            "unit": "seconds",
            "range": 0.181833990223515
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
        "date": 1681718306126,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.95588558726728,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.90016548757795,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.31035125268713,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2518695796989,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 141.477773290414,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.58859052153581,
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
        "date": 1681755634796,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 232.328,
            "unit": "seconds",
            "range": 0.203040882576428
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
        "date": 1681755644043,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.15260348845849,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.88884192954164,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.39137870343073,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.800670104096,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 142.807040177065,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.72040107378985,
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
        "date": 1681796160102,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 232.177,
            "unit": "seconds",
            "range": 0.233195411613946
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
        "date": 1681796168040,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.40848388147287,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.90776852188703,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.40081185287867,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6940050229357,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 164.224255022934,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.43847832308561,
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
        "date": 1681960464829,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 164.502,
            "unit": "seconds",
            "range": 0.183601742911915
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
        "date": 1681960472456,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.95588558726728,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.90016548757795,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.31035125268713,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2518695796989,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 141.477773290414,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.58859052153581,
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
        "date": 1682451999142,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 215.223,
            "unit": "seconds",
            "range": 0.816993329225781
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1682452006485,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.37318010105767,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.16737114123167,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.40060454861053,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7555252763289,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 141.886141942288,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.59228627570489,
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
        "date": 1682505951004,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 230.772,
            "unit": "seconds",
            "range": 0.15328926902835
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
        "date": 1682505958956,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.35251320425936,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.91915099143715,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.37299492422494,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9472009751759,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 81.6797573817577,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.46241880763895,
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
        "date": 1682546174841,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 281.627,
            "unit": "seconds",
            "range": 1.61383025749316
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
        "date": 1682546184633,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.95588558726728,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.90016548757795,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.31035125268713,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2518695796989,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 141.477773290414,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.58859052153581,
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
        "date": 1682622201526,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 176.57,
            "unit": "seconds",
            "range": 3.18538851633524
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
        "date": 1682622209621,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.49003840458162,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.07968276939436,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.40316640045569,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6783636630002,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 143.213465154251,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.65349739357498,
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
        "date": 1683218497017,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 146.611,
            "unit": "seconds",
            "range": 0.126463038075957
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
        "date": 1683218504901,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.34005003963681,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.21152816469245,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.44655802394222,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.8230528475217,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 165.687788789493,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.71390650619909,
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
        "date": 1683266471541,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 221.373,
            "unit": "seconds",
            "range": 0.185876571949277
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
        "date": 1683266480498,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.37258905697622,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.30645050912398,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.43019387446882,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6241023203276,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 65.682013747538,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.69120733767007,
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
        "date": 1683306934234,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 310.196,
            "unit": "seconds",
            "range": 0.591158523576856
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
        "date": 1683306941022,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.40435489865117,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 7.44572980615909,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.50251714920327,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.8757319170346,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 142.267166699983,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.455149210829,
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
        "date": 1683340045109,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 225.929,
            "unit": "seconds",
            "range": 0.533931549920959
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1683340052768,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.93880700013043,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.8438243414472,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.05304502286597,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6916426681282,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 78.6006755767918,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.17751562217427,
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
        "date": 1683657503222,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 284.095,
            "unit": "seconds",
            "range": 2.01174712625691
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1683657510697,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.90491572568576,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.16497865799282,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.18952076461701,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.5643115318306,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 79.5988521572255,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.00396050595154,
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
        "date": 1683683134994,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 299.417,
            "unit": "seconds",
            "range": 1.82795352785531
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
        "date": 1683683143183,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.93880700013043,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.8438243414472,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.05304502286597,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6916426681282,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 78.6006755767918,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.17751562217427,
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
        "date": 1683906884033,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 256.943,
            "unit": "seconds",
            "range": 0.822666457319122
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
        "date": 1683906892080,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.04863710074242,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.06824981733348,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.30870350374352,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7547550014842,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 191.965325949026,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 5.95261724146583,
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
        "date": 1683935959593,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 219.881,
            "unit": "seconds",
            "range": 0.216492263141475
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
        "date": 1683935968790,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.2307457691348,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.12568708939384,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.01424135719686,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.8146062547659,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 191.211269239738,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.05564049562558,
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
        "date": 1684358561978,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 209.285,
            "unit": "seconds",
            "range": 0.0629483915639545
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
        "date": 1684358569996,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.8742310654479,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.06094083141281,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.995451347241398,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.756989650466,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 193.847396048984,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 5.87884747329162,
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
        "date": 1684807921073,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 264.7,
            "unit": "seconds",
            "range": 0.372206931691557
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
        "date": 1684807929142,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.05632311027689,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.05016706125999,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.00680032371244,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7547550014842,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 191.965325949026,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 5.95261724146583,
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
        "date": 1684864495397,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 201.658,
            "unit": "seconds",
            "range": 0.12158782834273
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
        "date": 1684864504605,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.96100493213371,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.05873751865678,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.06641185982112,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.8860743427807,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 191.27612442949,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.03679087034369,
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
        "date": 1684886380770,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 202.226,
            "unit": "seconds",
            "range": 0.187409711596133
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
        "date": 1684886390174,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.9885462323776,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.15532531605401,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.15091374811842,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.8643067302427,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 192.55717236619,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.09061791726236,
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
        "date": 1684990710612,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 216.595,
            "unit": "seconds",
            "range": 0.375377809681642
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
        "date": 1684990719284,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.84030994559918,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.17935225487445,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.21849582539602,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.4869673669854,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 154.865069771587,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 5.96004795640386,
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
        "date": 1685044129420,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 173.431,
            "unit": "seconds",
            "range": 0.726925649568686
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
        "date": 1685044139210,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.00784990314007,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.16085731827904,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.239033051456,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.4563334681843,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 192.545155919715,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 5.97223540104491,
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
        "date": 1685137970137,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 246.78,
            "unit": "seconds",
            "range": 0.658889975640739
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1685137977579,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.56739232860638,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.18903947791977,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.26468036174333,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.4748029251001,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 192.956452060495,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 5.94140602502992,
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
        "date": 1685593566923,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 244,
            "unit": "seconds",
            "range": 0.233109416369715
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1685593574952,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.96235372065076,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.1804968025526,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.29545962139899,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6772787267708,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 192.193256122352,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 5.88977757093192,
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
        "date": 1685642323280,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 204.925,
            "unit": "seconds",
            "range": 0.245231523258018
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
        "date": 1685642330145,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.83058893955056,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.98039956663122,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.06059606539572,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.5675275525897,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 191.415674308323,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.12459068172068,
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
        "date": 1685741305977,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 273.225,
            "unit": "seconds",
            "range": 0.386911488585884
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
        "date": 1685741313988,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.60851727937052,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.13758541862713,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.2257324706347,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.462949209936,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 192.207263280988,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.05139084441407,
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
        "date": 1685763721615,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 248.36,
            "unit": "seconds",
            "range": 1.31814490857328
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
        "date": 1685763731611,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.83058893955056,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.98039956663122,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.06059606539572,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.5675275525897,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 191.415674308323,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.12459068172068,
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
        "date": 1685775877231,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 209.822,
            "unit": "seconds",
            "range": 0.101378498702854
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1685775885259,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.83002874168314,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.97438924285676,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.10074127188926,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.3570982265789,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 83.3494872322984,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.0851804415356,
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
        "date": 1686250144502,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 218.031,
            "unit": "seconds",
            "range": 0.175940046606163
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
        "date": 1686250152313,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.93266246974713,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.39545386054292,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.01358311157705,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.453618351017,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 69.6690965239668,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 5.93200964271884,
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
        "date": 1686268519294,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 244.325,
            "unit": "seconds",
            "range": 0.570580844400586
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
        "date": 1686268527552,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.08825951277246,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.99432492813471,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.15780130070103,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.2615342294598,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 186.11700056719,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.09629604787087,
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
        "date": 1686763469603,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 166.254,
            "unit": "seconds",
            "range": 0.243134530658155
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
        "date": 1686763478303,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.83002874168314,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.97438924285676,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.10074127188926,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.3570982265789,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 83.3494872322984,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.0851804415356,
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
        "date": 1687288638444,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 248.26,
            "unit": "seconds",
            "range": 0.839116201727462
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
        "date": 1687288646096,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.72074708948263,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.31954438430172,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.1119041179893,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.4996413764913,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 192.273984977573,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.06084840812869,
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
        "date": 1687382848221,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 242.237,
            "unit": "seconds",
            "range": 0.281517495016391
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
        "date": 1687382855461,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.57802337380474,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.11858503200993,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.04971583402274,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.4064182556131,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 157.154730450498,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.00787386508876,
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
        "date": 1687401347851,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 149.17,
            "unit": "seconds",
            "range": 0.179571712692209
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1687401355568,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.76948199735895,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 6.31035785788167,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.08971092734705,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.5585731490364,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 193.148716101583,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.11743884655478,
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
        "date": 1687556876324,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 121.512,
            "unit": "seconds",
            "range": 0.389830219454792
          }
        ]
      },
      {
        "commit": {
          "author": {
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
        "date": 1687556884210,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.95594381129803,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.97339955428138,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.07471215103281,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0520351452134,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 72.5282581458872,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.06542955244775,
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
        "date": 1688014422630,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 229.766,
            "unit": "seconds",
            "range": 1.49921592841006
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
        "date": 1688014431566,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.54202032277833,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.76302937052992,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.2081664473679,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7022799249084,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 191.445301173917,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.20197792067725,
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
        "date": 1688608361961,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 196.044,
            "unit": "seconds",
            "range": 0.223898191148928
          }
        ]
      },
      {
        "commit": {
          "author": {
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
        "date": 1688608369986,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.59452945357607,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.26756694781851,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.14434504832383,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.5157624963898,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 156.391873945702,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 5.90987743828344,
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
        "date": 1688772304355,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 235.051,
            "unit": "seconds",
            "range": 0.901518108525551
          }
        ]
      },
      {
        "commit": {
          "author": {
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
        "date": 1688772313281,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.05015729768744,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.8220792019277,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.08524353349698,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6705567212596,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 190.95681502939,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 5.95070663917398,
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
        "date": 1688803551570,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 198.02,
            "unit": "seconds",
            "range": 0.203612376829696
          }
        ]
      },
      {
        "commit": {
          "author": {
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
        "date": 1688803559105,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.66366645570138,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.3945807106886,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.34484984684749,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6420163611984,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 192.797481179007,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.04556142629178,
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
        "date": 1689017070665,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 267.609,
            "unit": "seconds",
            "range": 0.211160839167807
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
        "date": 1689017078120,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.62874201209761,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.68123260560176,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.29532804831089,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6322889682008,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 193.034302029219,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 6.30119565041707,
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
        "date": 1689106724832,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 185.673,
            "unit": "seconds",
            "range": 0.107889295119519
          }
        ]
      },
      {
        "commit": {
          "author": {
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
        "date": 1689106733693,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.7347365361624,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.59737490251823,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.826032170943457,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9391299921939,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 62.4722673081793,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.34492160456599,
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
        "date": 1689129319166,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 255.502,
            "unit": "seconds",
            "range": 1.2568904486862
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
        "date": 1689129326017,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.81013956968771,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.72175777570985,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.842397693470851,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.8484174753897,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 64.2141604748231,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.47962594839969,
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
        "date": 1689189097407,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 203.231,
            "unit": "seconds",
            "range": 0.662550299976914
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
        "date": 1689189104931,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.7347365361624,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.59737490251823,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.826032170943457,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9391299921939,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 62.4722673081793,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.34492160456599,
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
        "date": 1689226032522,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 208.792,
            "unit": "seconds",
            "range": 0.139869939582979
          }
        ]
      },
      {
        "commit": {
          "author": {
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
        "date": 1689226041572,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.35731384777102,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.2372091657514,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.829918672072195,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0522420518317,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 69.9114796424333,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.26406860997775,
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
          "id": "9874ca3b4ade45f41c7c185abe202378dd152829",
          "message": "fix: Correct class reference",
          "timestamp": "2023-07-13T07:18:43-07:00",
          "tree_id": "8ea083cf4817a7613663bff91c3156ba9c7a195e",
          "url": "https://github.com/galacticusorg/galacticus/commit/9874ca3b4ade45f41c7c185abe202378dd152829"
        },
        "date": 1689270181165,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 309.241,
            "unit": "seconds",
            "range": 0.701655827311403
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
          "id": "9874ca3b4ade45f41c7c185abe202378dd152829",
          "message": "fix: Correct class reference",
          "timestamp": "2023-07-13T07:18:43-07:00",
          "tree_id": "8ea083cf4817a7613663bff91c3156ba9c7a195e",
          "url": "https://github.com/galacticusorg/galacticus/commit/9874ca3b4ade45f41c7c185abe202378dd152829"
        },
        "date": 1689270189761,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.7347365361624,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.59737490251823,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.826032170943457,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9391299921939,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 62.4722673081793,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.34492160456599,
            "unit": "-logℒ"
          }
        ]
      }
    ]
  }
}