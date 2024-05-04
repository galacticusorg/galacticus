window.BENCHMARK_DATA = {
  "lastUpdate": 1714798062050,
  "repoUrl": "https://github.com/galacticusorg/galacticus",
  "entries": {
    "Milky Way model benchmarks": [
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
          "id": "7a5e7e98ef62bc89ffc1c4004efd57549b53b16b",
          "message": "fix: Patch `libmatheval-1.1.12` for MacOS builds",
          "timestamp": "2023-07-26T21:36:04Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/7a5e7e98ef62bc89ffc1c4004efd57549b53b16b"
        },
        "date": 1690418876550,
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
          "id": "f46fbef1ee3c056e0a9e74c0458c48ce17976b4c",
          "message": "Merge pull request #442 from galacticusorg/buildRefactor\n\nReduce code duplication in `functionClass` `descriptor()` method builder",
          "timestamp": "2023-07-28T05:17:04Z",
          "tree_id": "397fa2bc4c2f47d56f46760567aa2c89f83e68d4",
          "url": "https://github.com/galacticusorg/galacticus/commit/f46fbef1ee3c056e0a9e74c0458c48ce17976b4c"
        },
        "date": 1690532341325,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 206.667,
            "unit": "seconds",
            "range": 0.489981734353731
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "f46fbef1ee3c056e0a9e74c0458c48ce17976b4c",
          "message": "Merge pull request #442 from galacticusorg/buildRefactor\n\nReduce code duplication in `functionClass` `descriptor()` method builder",
          "timestamp": "2023-07-28T05:17:04Z",
          "tree_id": "397fa2bc4c2f47d56f46760567aa2c89f83e68d4",
          "url": "https://github.com/galacticusorg/galacticus/commit/f46fbef1ee3c056e0a9e74c0458c48ce17976b4c"
        },
        "date": 1690532352555,
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
            "name": "Andrew Benson",
            "username": "abensonca",
            "email": "abenson@obs.carnegiescience.edu"
          },
          "committer": {
            "name": "Andrew Benson",
            "username": "abensonca",
            "email": "abenson@obs.carnegiescience.edu"
          },
          "id": "f93400effb7083a742d6565a5b634cee23f8deab",
          "message": "fix: Test for build warnings only in the non-static build\n\nAvoids any possible problems with warnings from the static build regarding linking of static libraries.",
          "timestamp": "2023-07-28T22:37:53Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/f93400effb7083a742d6565a5b634cee23f8deab"
        },
        "date": 1690602355639,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 198.511,
            "unit": "seconds",
            "range": 0.142831719164083
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
          "id": "f93400effb7083a742d6565a5b634cee23f8deab",
          "message": "fix: Test for build warnings only in the non-static build\n\nAvoids any possible problems with warnings from the static build regarding linking of static libraries.",
          "timestamp": "2023-07-28T22:37:53Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/f93400effb7083a742d6565a5b634cee23f8deab"
        },
        "date": 1690602364596,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.83622231819236,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.02183701424305,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.858044315919163,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7932954184213,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 68.8080477368147,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.46954705432508,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "bd3358fc3d5e9193b16179658b145d8cda342de3",
          "message": "Merge pull request #443 from galacticusorg/buildRefactor\n\nReduce code duplication in state store/restore method builder",
          "timestamp": "2023-07-29T16:33:18Z",
          "tree_id": "4368532e9b45c02424efdff762c0c3fe55b81f79",
          "url": "https://github.com/galacticusorg/galacticus/commit/bd3358fc3d5e9193b16179658b145d8cda342de3"
        },
        "date": 1690675733094,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 207.49,
            "unit": "seconds",
            "range": 0.101173118961713
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "bd3358fc3d5e9193b16179658b145d8cda342de3",
          "message": "Merge pull request #443 from galacticusorg/buildRefactor\n\nReduce code duplication in state store/restore method builder",
          "timestamp": "2023-07-29T16:33:18Z",
          "tree_id": "4368532e9b45c02424efdff762c0c3fe55b81f79",
          "url": "https://github.com/galacticusorg/galacticus/commit/bd3358fc3d5e9193b16179658b145d8cda342de3"
        },
        "date": 1690675740470,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.83622231819236,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.02183701424305,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.858044315919163,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7932954184213,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 68.8080477368147,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.46954705432508,
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
          "id": "3fea25a0a779c2b879057cf7328e7c70e15ec77c",
          "message": "fix: Deallocate object to avoid memory leak",
          "timestamp": "2023-07-31T03:33:40Z",
          "tree_id": "22001050019317d78ef5424a28a1eb411d04349e",
          "url": "https://github.com/galacticusorg/galacticus/commit/3fea25a0a779c2b879057cf7328e7c70e15ec77c"
        },
        "date": 1690786449521,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 196.53,
            "unit": "seconds",
            "range": 0.463158720095747
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
          "id": "3fea25a0a779c2b879057cf7328e7c70e15ec77c",
          "message": "fix: Deallocate object to avoid memory leak",
          "timestamp": "2023-07-31T03:33:40Z",
          "tree_id": "22001050019317d78ef5424a28a1eb411d04349e",
          "url": "https://github.com/galacticusorg/galacticus/commit/3fea25a0a779c2b879057cf7328e7c70e15ec77c"
        },
        "date": 1690786456726,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.83622231819236,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.02183701424305,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.858044315919163,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7932954184213,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 68.8080477368147,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.46954705432508,
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
          "id": "f7f8371f7d5d628d316d6eb76a4c12e965dfa7c8",
          "message": "fix: Fix syntax",
          "timestamp": "2023-07-31T08:18:56-07:00",
          "tree_id": "621d7216d3e9eaf25349b9c6cac9ddccdf2b7acd",
          "url": "https://github.com/galacticusorg/galacticus/commit/f7f8371f7d5d628d316d6eb76a4c12e965dfa7c8"
        },
        "date": 1690830149213,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 225.964,
            "unit": "seconds",
            "range": 0.121739065210309
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
          "id": "f7f8371f7d5d628d316d6eb76a4c12e965dfa7c8",
          "message": "fix: Fix syntax",
          "timestamp": "2023-07-31T08:18:56-07:00",
          "tree_id": "621d7216d3e9eaf25349b9c6cac9ddccdf2b7acd",
          "url": "https://github.com/galacticusorg/galacticus/commit/f7f8371f7d5d628d316d6eb76a4c12e965dfa7c8"
        },
        "date": 1690830156844,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.27601298855862,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.50668579880912,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.856109814363159,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.969225876387,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 70.7473742298182,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.35131210320676,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2da1da445ed3f84caabd4ae61c2d73c9fe015a3b",
          "message": "Merge pull request #444 from galacticusorg/projectedDensity\n\nInclude mass outside of the virial radius in projected density calculations",
          "timestamp": "2023-08-01T01:37:00Z",
          "tree_id": "3dcca2c14dbe14011d83c7e060dc9e0dbc4e34f4",
          "url": "https://github.com/galacticusorg/galacticus/commit/2da1da445ed3f84caabd4ae61c2d73c9fe015a3b"
        },
        "date": 1690866904119,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 169.642,
            "unit": "seconds",
            "range": 0.178167337073197
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2da1da445ed3f84caabd4ae61c2d73c9fe015a3b",
          "message": "Merge pull request #444 from galacticusorg/projectedDensity\n\nInclude mass outside of the virial radius in projected density calculations",
          "timestamp": "2023-08-01T01:37:00Z",
          "tree_id": "3dcca2c14dbe14011d83c7e060dc9e0dbc4e34f4",
          "url": "https://github.com/galacticusorg/galacticus/commit/2da1da445ed3f84caabd4ae61c2d73c9fe015a3b"
        },
        "date": 1690866911928,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.35082248599075,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.1104275370683,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.725536547467159,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0132634550444,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 67.8364127988294,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.25211122078835,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "55cac102888b7b0324e90e11ddc2e673fa851131",
          "message": "Merge pull request #445 from galacticusorg/buildRefactor\n\nReplace `include` directives with events in some trivial cases",
          "timestamp": "2023-08-02T01:22:21Z",
          "tree_id": "8e731ac7269422465ea73bfbf9eb0741167bf769",
          "url": "https://github.com/galacticusorg/galacticus/commit/55cac102888b7b0324e90e11ddc2e673fa851131"
        },
        "date": 1690951072163,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 223.701,
            "unit": "seconds",
            "range": 0.0787965735289962
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "55cac102888b7b0324e90e11ddc2e673fa851131",
          "message": "Merge pull request #445 from galacticusorg/buildRefactor\n\nReplace `include` directives with events in some trivial cases",
          "timestamp": "2023-08-02T01:22:21Z",
          "tree_id": "8e731ac7269422465ea73bfbf9eb0741167bf769",
          "url": "https://github.com/galacticusorg/galacticus/commit/55cac102888b7b0324e90e11ddc2e673fa851131"
        },
        "date": 1690951080792,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.27601298855862,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.50668579880912,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.856109814363159,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.969225876387,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 70.7473742298182,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.35131210320676,
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
          "id": "d7c95ab107aca859ee58abcbda553d5d5df46b33",
          "message": "fix: Mark variable as `threadprivate`\n\nVariable `treeTimeLatest` in `taskEvolveForests` class was not `threadprivate`, but is used independently without locking by all threads. Marking it `threadprivate` avoids possible race conditions.",
          "timestamp": "2023-08-02T14:27:27Z",
          "tree_id": "c66d3c0a485581f3a76d1b2f49ce139a509b137e",
          "url": "https://github.com/galacticusorg/galacticus/commit/d7c95ab107aca859ee58abcbda553d5d5df46b33"
        },
        "date": 1690998488115,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 267.341,
            "unit": "seconds",
            "range": 0.508524237374967
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
          "id": "d7c95ab107aca859ee58abcbda553d5d5df46b33",
          "message": "fix: Mark variable as `threadprivate`\n\nVariable `treeTimeLatest` in `taskEvolveForests` class was not `threadprivate`, but is used independently without locking by all threads. Marking it `threadprivate` avoids possible race conditions.",
          "timestamp": "2023-08-02T14:27:27Z",
          "tree_id": "c66d3c0a485581f3a76d1b2f49ce139a509b137e",
          "url": "https://github.com/galacticusorg/galacticus/commit/d7c95ab107aca859ee58abcbda553d5d5df46b33"
        },
        "date": 1690998497947,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.34325456084193,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.3830440654232,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.723674755626745,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9950179067304,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 69.0038856999168,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.27697524047164,
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
          "id": "5b2eae9907a38e091f3228cf60551515d2280cfa",
          "message": "fix: Improve formatting of error message",
          "timestamp": "2023-08-03T16:59:23-07:00",
          "tree_id": "927595231056b72ecab922d996b75d34c175bb00",
          "url": "https://github.com/galacticusorg/galacticus/commit/5b2eae9907a38e091f3228cf60551515d2280cfa"
        },
        "date": 1691118282223,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 175.838,
            "unit": "seconds",
            "range": 0.261315135421407
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
          "id": "5b2eae9907a38e091f3228cf60551515d2280cfa",
          "message": "fix: Improve formatting of error message",
          "timestamp": "2023-08-03T16:59:23-07:00",
          "tree_id": "927595231056b72ecab922d996b75d34c175bb00",
          "url": "https://github.com/galacticusorg/galacticus/commit/5b2eae9907a38e091f3228cf60551515d2280cfa"
        },
        "date": 1691118290360,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.98928711762356,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.63688368400778,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.770036956075677,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7236663735307,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 214.003002527713,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.32369340966966,
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
          "id": "a1ec083424e1f623e08cd4969cdbe96e7779dbd4",
          "message": "fix: Split long YAML list across lines for easier `diff`ing",
          "timestamp": "2023-08-08T15:30:50Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/a1ec083424e1f623e08cd4969cdbe96e7779dbd4"
        },
        "date": 1691522419028,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 265.236,
            "unit": "seconds",
            "range": 0.575427145692556
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
          "id": "a1ec083424e1f623e08cd4969cdbe96e7779dbd4",
          "message": "fix: Split long YAML list across lines for easier `diff`ing",
          "timestamp": "2023-08-08T15:30:50Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/a1ec083424e1f623e08cd4969cdbe96e7779dbd4"
        },
        "date": 1691522429800,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.35082248599075,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.1104275370683,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.725536547467159,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0132634550444,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 67.8364127988294,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.25211122078835,
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
          "id": "0b289fc3dc9e5087db1763272d0ebca2f14b097d",
          "message": "fix: Add default for parameter",
          "timestamp": "2023-08-09T17:54:31Z",
          "tree_id": "3bbb529955d104986c83a309df39e9924f0097d2",
          "url": "https://github.com/galacticusorg/galacticus/commit/0b289fc3dc9e5087db1763272d0ebca2f14b097d"
        },
        "date": 1691615747939,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 149.551,
            "unit": "seconds",
            "range": 0.182375711102372
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
          "id": "0b289fc3dc9e5087db1763272d0ebca2f14b097d",
          "message": "fix: Add default for parameter",
          "timestamp": "2023-08-09T17:54:31Z",
          "tree_id": "3bbb529955d104986c83a309df39e9924f0097d2",
          "url": "https://github.com/galacticusorg/galacticus/commit/0b289fc3dc9e5087db1763272d0ebca2f14b097d"
        },
        "date": 1691615757908,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.74297669185497,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.64751495885572,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.12206627571949,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2296904135044,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 89.6220240964809,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.844703095162,
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
          "id": "a375dcde2702f204d63d2ed2c58fdf7cea779ec0",
          "message": "feat: Add functionality to the `mergerTreeOperatorRegrid` class\n\nAdds an option to prevent removal of nodes at  ungridded times (i.e. new nodes are inserted at the grid times, but no nodes are removed). Also refactors the code to allow a single grid time when ungridded nodes are not to be removed.",
          "timestamp": "2023-08-09T23:26:36Z",
          "tree_id": "3bc6c70419311c62580075b4a266277eae16c01c",
          "url": "https://github.com/galacticusorg/galacticus/commit/a375dcde2702f204d63d2ed2c58fdf7cea779ec0"
        },
        "date": 1691635995179,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 212.132,
            "unit": "seconds",
            "range": 0.605443308659505
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
          "id": "a375dcde2702f204d63d2ed2c58fdf7cea779ec0",
          "message": "feat: Add functionality to the `mergerTreeOperatorRegrid` class\n\nAdds an option to prevent removal of nodes at  ungridded times (i.e. new nodes are inserted at the grid times, but no nodes are removed). Also refactors the code to allow a single grid time when ungridded nodes are not to be removed.",
          "timestamp": "2023-08-09T23:26:36Z",
          "tree_id": "3bc6c70419311c62580075b4a266277eae16c01c",
          "url": "https://github.com/galacticusorg/galacticus/commit/a375dcde2702f204d63d2ed2c58fdf7cea779ec0"
        },
        "date": 1691636001889,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.27601298855862,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.50668579880912,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.856109814363159,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.969225876387,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 70.7473742298182,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.35131210320676,
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
          "id": "03f1be418ea2846f56cee7c299f4ae208dc6c216",
          "message": "fix(style): Correct a typo in a comment",
          "timestamp": "2023-08-10T23:33:15Z",
          "tree_id": "fda400cd223a48cfd05abfdb04f4c298cdbaf9e7",
          "url": "https://github.com/galacticusorg/galacticus/commit/03f1be418ea2846f56cee7c299f4ae208dc6c216"
        },
        "date": 1691722084532,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 192.464,
            "unit": "seconds",
            "range": 0.426851730697404
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
          "id": "03f1be418ea2846f56cee7c299f4ae208dc6c216",
          "message": "fix(style): Correct a typo in a comment",
          "timestamp": "2023-08-10T23:33:15Z",
          "tree_id": "fda400cd223a48cfd05abfdb04f4c298cdbaf9e7",
          "url": "https://github.com/galacticusorg/galacticus/commit/03f1be418ea2846f56cee7c299f4ae208dc6c216"
        },
        "date": 1691722093507,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.95668999492098,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.71149374820828,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.780308933436236,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6022197906029,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 88.9508519609884,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.34122822663242,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "7f9ef8648b3d04c59a53383d47dd1b7c07703fdf",
          "message": "fix: Ensure count of node labels is correct if none exist",
          "timestamp": "2023-08-10T22:16:51-07:00",
          "tree_id": "a43aadcef4d33ab8eeee927946b7965737eba2a9",
          "url": "https://github.com/galacticusorg/galacticus/commit/7f9ef8648b3d04c59a53383d47dd1b7c07703fdf"
        },
        "date": 1691741631051,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 319.829,
            "unit": "seconds",
            "range": 0.645961995785535
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "7f9ef8648b3d04c59a53383d47dd1b7c07703fdf",
          "message": "fix: Ensure count of node labels is correct if none exist",
          "timestamp": "2023-08-10T22:16:51-07:00",
          "tree_id": "a43aadcef4d33ab8eeee927946b7965737eba2a9",
          "url": "https://github.com/galacticusorg/galacticus/commit/7f9ef8648b3d04c59a53383d47dd1b7c07703fdf"
        },
        "date": 1691741639028,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.27601298855862,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.50668579880912,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.856109814363159,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.969225876387,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 70.7473742298182,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.35131210320676,
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
          "id": "ea48a439dc6908397670624edf80cac0e0f2c486",
          "message": "feat: Add the `rate()` method to the `mergerTreeBranchingProbabilityGnrlzdPrssSchchtr` class",
          "timestamp": "2023-08-14T20:55:05Z",
          "tree_id": "d13bcaa10158865edd50e753fddabd74c9438e16",
          "url": "https://github.com/galacticusorg/galacticus/commit/ea48a439dc6908397670624edf80cac0e0f2c486"
        },
        "date": 1692069745202,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 192.3,
            "unit": "seconds",
            "range": 0.905170702132917
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
          "id": "ea48a439dc6908397670624edf80cac0e0f2c486",
          "message": "feat: Add the `rate()` method to the `mergerTreeBranchingProbabilityGnrlzdPrssSchchtr` class",
          "timestamp": "2023-08-14T20:55:05Z",
          "tree_id": "d13bcaa10158865edd50e753fddabd74c9438e16",
          "url": "https://github.com/galacticusorg/galacticus/commit/ea48a439dc6908397670624edf80cac0e0f2c486"
        },
        "date": 1692069753384,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.27601298855862,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.50668579880912,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.856109814363159,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.969225876387,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 70.7473742298182,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.35131210320676,
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
          "id": "2081039481ead7c61f9b56dc67890ecf32755e16",
          "message": "fix(style): Formatting only",
          "timestamp": "2023-08-16T00:18:01Z",
          "tree_id": "3df04dc9fbd5ead8203330071042a65f7523f28a",
          "url": "https://github.com/galacticusorg/galacticus/commit/2081039481ead7c61f9b56dc67890ecf32755e16"
        },
        "date": 1692156163571,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 193.282,
            "unit": "seconds",
            "range": 0.0795210663877842
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
          "id": "2081039481ead7c61f9b56dc67890ecf32755e16",
          "message": "fix(style): Formatting only",
          "timestamp": "2023-08-16T00:18:01Z",
          "tree_id": "3df04dc9fbd5ead8203330071042a65f7523f28a",
          "url": "https://github.com/galacticusorg/galacticus/commit/2081039481ead7c61f9b56dc67890ecf32755e16"
        },
        "date": 1692156173348,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.44586931890625,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.64634609139956,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.900187219452953,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0654174765943,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 70.9525601138061,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.30750810796761,
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
          "id": "106f00ca438b3a4a4be3bc18ad12ba1e8dab3426",
          "message": "fix: Update cooling function test to match results from Cloudy v23.00",
          "timestamp": "2023-08-16T07:23:01-07:00",
          "tree_id": "4859a674c78a2353d94e87b714e2cb5bdbda9ee4",
          "url": "https://github.com/galacticusorg/galacticus/commit/106f00ca438b3a4a4be3bc18ad12ba1e8dab3426"
        },
        "date": 1692216465572,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 268.379,
            "unit": "seconds",
            "range": 0.963504488831343
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
          "id": "106f00ca438b3a4a4be3bc18ad12ba1e8dab3426",
          "message": "fix: Update cooling function test to match results from Cloudy v23.00",
          "timestamp": "2023-08-16T07:23:01-07:00",
          "tree_id": "4859a674c78a2353d94e87b714e2cb5bdbda9ee4",
          "url": "https://github.com/galacticusorg/galacticus/commit/106f00ca438b3a4a4be3bc18ad12ba1e8dab3426"
        },
        "date": 1692216474015,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.61900842182208,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.588824960662,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.799608411284678,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6785721019041,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 67.7384154401288,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.24778572576594,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c71a0f05613948ac93ffbf2d36742543f5dba130",
          "message": "Merge pull request #449 from galacticusorg/projectedHalos\n\nImplement projected mass profile calculations",
          "timestamp": "2023-08-17T02:17:29Z",
          "tree_id": "b0c21707056ee797cf7cd8b0c916bc2d3675cac9",
          "url": "https://github.com/galacticusorg/galacticus/commit/c71a0f05613948ac93ffbf2d36742543f5dba130"
        },
        "date": 1692249708879,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 174.75,
            "unit": "seconds",
            "range": 0.17840403583069
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c71a0f05613948ac93ffbf2d36742543f5dba130",
          "message": "Merge pull request #449 from galacticusorg/projectedHalos\n\nImplement projected mass profile calculations",
          "timestamp": "2023-08-17T02:17:29Z",
          "tree_id": "b0c21707056ee797cf7cd8b0c916bc2d3675cac9",
          "url": "https://github.com/galacticusorg/galacticus/commit/c71a0f05613948ac93ffbf2d36742543f5dba130"
        },
        "date": 1692249716490,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.43313684038849,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.0755070150941,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.02770619724301,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9047305636468,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 70.8471614676265,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.37390867528244,
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
          "id": "2e4d4ad34084440a26bcdfadc33601a4f16e53d7",
          "message": "fix: Remove OpenMP locks from `accretionDisksADAF` class\n\nThis are obsolete (and have been for a long time) as the tabulated solutions are now stored in each object of this class (and so are not shared).",
          "timestamp": "2023-08-21T15:29:58-07:00",
          "tree_id": "437fdaefcaddcf9efc219bb970f048b1e95bab2e",
          "url": "https://github.com/galacticusorg/galacticus/commit/2e4d4ad34084440a26bcdfadc33601a4f16e53d7"
        },
        "date": 1692668819083,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 237.131,
            "unit": "seconds",
            "range": 0.294276910407664
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
          "id": "2e4d4ad34084440a26bcdfadc33601a4f16e53d7",
          "message": "fix: Remove OpenMP locks from `accretionDisksADAF` class\n\nThis are obsolete (and have been for a long time) as the tabulated solutions are now stored in each object of this class (and so are not shared).",
          "timestamp": "2023-08-21T15:29:58-07:00",
          "tree_id": "437fdaefcaddcf9efc219bb970f048b1e95bab2e",
          "url": "https://github.com/galacticusorg/galacticus/commit/2e4d4ad34084440a26bcdfadc33601a4f16e53d7"
        },
        "date": 1692668827913,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.61900842182208,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.588824960662,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.799608411284678,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6785721019041,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 67.7384154401288,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.24778572576594,
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
          "id": "6ab982fafd44fd21a6631ae92238191f4ff50985",
          "message": "fix: Correct variable attributes",
          "timestamp": "2023-08-22T13:39:43-07:00",
          "tree_id": "2ad7dad4ffc88565fb13151fd06ab6462ab25040",
          "url": "https://github.com/galacticusorg/galacticus/commit/6ab982fafd44fd21a6631ae92238191f4ff50985"
        },
        "date": 1692757991521,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 237.43,
            "unit": "seconds",
            "range": 0.255968748092625
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
          "id": "6ab982fafd44fd21a6631ae92238191f4ff50985",
          "message": "fix: Correct variable attributes",
          "timestamp": "2023-08-22T13:39:43-07:00",
          "tree_id": "2ad7dad4ffc88565fb13151fd06ab6462ab25040",
          "url": "https://github.com/galacticusorg/galacticus/commit/6ab982fafd44fd21a6631ae92238191f4ff50985"
        },
        "date": 1692757998969,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.39228673805545,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.77557720723732,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.853704837402535,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1731965439773,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 70.1000719336354,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.29999257785656,
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
          "id": "7cd7b6df772e1a6bde6d74da6fb8fc894db2bee6",
          "message": "fix: Correct method descriptions\n\nIn several places \"cm\" was missing the \"m\".",
          "timestamp": "2023-08-23T07:39:03-07:00",
          "tree_id": "86aaa1b252e0b8d7d66fb1be787200c73580a158",
          "url": "https://github.com/galacticusorg/galacticus/commit/7cd7b6df772e1a6bde6d74da6fb8fc894db2bee6"
        },
        "date": 1692812647202,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 208.291,
            "unit": "seconds",
            "range": 0.214240285657083
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
          "id": "7cd7b6df772e1a6bde6d74da6fb8fc894db2bee6",
          "message": "fix: Correct method descriptions\n\nIn several places \"cm\" was missing the \"m\".",
          "timestamp": "2023-08-23T07:39:03-07:00",
          "tree_id": "86aaa1b252e0b8d7d66fb1be787200c73580a158",
          "url": "https://github.com/galacticusorg/galacticus/commit/7cd7b6df772e1a6bde6d74da6fb8fc894db2bee6"
        },
        "date": 1692812654391,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.42795701008481,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.3236158898131,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.840619438436014,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0612166943218,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 69.8922100635251,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.27384044800987,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "3d09e973a1c5355740ebeccc1e20eb5b71517979",
          "message": "Merge pull request #451 from galacticusorg/singleStep\n\nAdd a `mergerTreeBuildContoller` class that takes just a single step in the tree",
          "timestamp": "2023-08-24T05:25:38Z",
          "tree_id": "ccc351bea3a2dd0c65b1479f4de1ea85c2748c4d",
          "url": "https://github.com/galacticusorg/galacticus/commit/3d09e973a1c5355740ebeccc1e20eb5b71517979"
        },
        "date": 1692865896540,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 233.743,
            "unit": "seconds",
            "range": 0.251097789713961
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "3d09e973a1c5355740ebeccc1e20eb5b71517979",
          "message": "Merge pull request #451 from galacticusorg/singleStep\n\nAdd a `mergerTreeBuildContoller` class that takes just a single step in the tree",
          "timestamp": "2023-08-24T05:25:38Z",
          "tree_id": "ccc351bea3a2dd0c65b1479f4de1ea85c2748c4d",
          "url": "https://github.com/galacticusorg/galacticus/commit/3d09e973a1c5355740ebeccc1e20eb5b71517979"
        },
        "date": 1692865905786,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.51072909978044,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.4519243103872,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.827524153182893,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9474091034214,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 68.988705041606,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.27661123302472,
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
          "id": "c28a88f038736a156970f43dd0eb41037ff9f80d",
          "message": "fix: Add highlighted \"WARNING\" on mismatched cosmological parameters",
          "timestamp": "2023-08-24T15:17:02Z",
          "tree_id": "3e47956d4d064e5b4a07ed5ca59a1313553fb19e",
          "url": "https://github.com/galacticusorg/galacticus/commit/c28a88f038736a156970f43dd0eb41037ff9f80d"
        },
        "date": 1692912118942,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 228.513,
            "unit": "seconds",
            "range": 1.08628545972038
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
          "id": "c28a88f038736a156970f43dd0eb41037ff9f80d",
          "message": "fix: Add highlighted \"WARNING\" on mismatched cosmological parameters",
          "timestamp": "2023-08-24T15:17:02Z",
          "tree_id": "3e47956d4d064e5b4a07ed5ca59a1313553fb19e",
          "url": "https://github.com/galacticusorg/galacticus/commit/c28a88f038736a156970f43dd0eb41037ff9f80d"
        },
        "date": 1692912126773,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.55178149493395,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.833480318933,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.727384505571875,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1054424083548,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 68.3982979636022,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.26439024807642,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "4c020b1e878370a7337ed241db559f46d7d8e287",
          "message": "Merge pull request #452 from galacticusorg/schneiderCatch\n\nCatch numerical failures in the `darkMatterProfileConcentrationSchneider2015` class",
          "timestamp": "2023-08-24T22:00:40Z",
          "tree_id": "26538619dbd650fda193be4a388c37db9a479753",
          "url": "https://github.com/galacticusorg/galacticus/commit/4c020b1e878370a7337ed241db559f46d7d8e287"
        },
        "date": 1692925342653,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 235.509,
            "unit": "seconds",
            "range": 2.44606109899149
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "4c020b1e878370a7337ed241db559f46d7d8e287",
          "message": "Merge pull request #452 from galacticusorg/schneiderCatch\n\nCatch numerical failures in the `darkMatterProfileConcentrationSchneider2015` class",
          "timestamp": "2023-08-24T22:00:40Z",
          "tree_id": "26538619dbd650fda193be4a388c37db9a479753",
          "url": "https://github.com/galacticusorg/galacticus/commit/4c020b1e878370a7337ed241db559f46d7d8e287"
        },
        "date": 1692925350944,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.55178149493395,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.833480318933,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.727384505571875,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1054424083548,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 68.3982979636022,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.26439024807642,
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
          "id": "03df2130a7ba1dc0e7fec6f3af0a5a0e62d0eb54",
          "message": "feat: Add support for reading/writing of HDF5 3D ragged arrays",
          "timestamp": "2023-08-25T15:08:25-07:00",
          "tree_id": "dcfddf102aa3e28040cca905b62f75ea22ceb702",
          "url": "https://github.com/galacticusorg/galacticus/commit/03df2130a7ba1dc0e7fec6f3af0a5a0e62d0eb54"
        },
        "date": 1693017639119,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 204.975,
            "unit": "seconds",
            "range": 0.236255158669681
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
          "id": "03df2130a7ba1dc0e7fec6f3af0a5a0e62d0eb54",
          "message": "feat: Add support for reading/writing of HDF5 3D ragged arrays",
          "timestamp": "2023-08-25T15:08:25-07:00",
          "tree_id": "dcfddf102aa3e28040cca905b62f75ea22ceb702",
          "url": "https://github.com/galacticusorg/galacticus/commit/03df2130a7ba1dc0e7fec6f3af0a5a0e62d0eb54"
        },
        "date": 1693017647097,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.66546592943343,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.37085191838936,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.901394625736313,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7754860306449,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 69.2646162731899,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.9446588628862,
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
          "id": "0eb84762c6daca60c08640e1bc1b1bb211ec131a",
          "message": "feat: Add node and tree indices to debug log output",
          "timestamp": "2023-08-28T23:54:08Z",
          "tree_id": "5ea548ab11fae8562973877cb37b4c3572e70c40",
          "url": "https://github.com/galacticusorg/galacticus/commit/0eb84762c6daca60c08640e1bc1b1bb211ec131a"
        },
        "date": 1693286086238,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 234.963,
            "unit": "seconds",
            "range": 0.306806290677302
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
          "id": "0eb84762c6daca60c08640e1bc1b1bb211ec131a",
          "message": "feat: Add node and tree indices to debug log output",
          "timestamp": "2023-08-28T23:54:08Z",
          "tree_id": "5ea548ab11fae8562973877cb37b4c3572e70c40",
          "url": "https://github.com/galacticusorg/galacticus/commit/0eb84762c6daca60c08640e1bc1b1bb211ec131a"
        },
        "date": 1693286093356,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.66546592943343,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.37085191838936,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.901394625736313,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.7754860306449,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 69.2646162731899,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.9446588628862,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "ca41d8f3a9ec90d7b02161cec967e118949cd44f",
          "message": "Merge pull request #456 from galacticusorg/completionStatus\n\nAllow clean exit on task failure",
          "timestamp": "2023-08-31T14:30:38Z",
          "tree_id": "fc45d38eb8ef468fe8fca2d80e259965f51d3f8b",
          "url": "https://github.com/galacticusorg/galacticus/commit/ca41d8f3a9ec90d7b02161cec967e118949cd44f"
        },
        "date": 1693514173266,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 282.166,
            "unit": "seconds",
            "range": 3.90820424236017
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "ca41d8f3a9ec90d7b02161cec967e118949cd44f",
          "message": "Merge pull request #456 from galacticusorg/completionStatus\n\nAllow clean exit on task failure",
          "timestamp": "2023-08-31T14:30:38Z",
          "tree_id": "fc45d38eb8ef468fe8fca2d80e259965f51d3f8b",
          "url": "https://github.com/galacticusorg/galacticus/commit/ca41d8f3a9ec90d7b02161cec967e118949cd44f"
        },
        "date": 1693514182135,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.36326745327789,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.39830552400072,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.15503965592705,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0983850857773,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 68.5681473709282,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.7801728110126,
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
          "id": "db8121c69e6a0867d6a5d7d9bfa79f41efb001ab",
          "message": "fix: Use the deploy app to allow comments on PRs from external collaborators",
          "timestamp": "2023-09-01T18:11:16Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/db8121c69e6a0867d6a5d7d9bfa79f41efb001ab"
        },
        "date": 1693605236751,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 195.442,
            "unit": "seconds",
            "range": 0.692867664131626
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
          "id": "db8121c69e6a0867d6a5d7d9bfa79f41efb001ab",
          "message": "fix: Use the deploy app to allow comments on PRs from external collaborators",
          "timestamp": "2023-09-01T18:11:16Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/db8121c69e6a0867d6a5d7d9bfa79f41efb001ab"
        },
        "date": 1693605246281,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.34010807231743,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.05089950515783,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.995003358346423,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2806285130524,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 93.4045161107519,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.29736523344343,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "94ded2ded9c80fed323bfe76c2075977081ab394",
          "message": "Merge pull request #453 from galacticusorg/optimization\n\nBaryonic physics optimizations",
          "timestamp": "2023-09-05T14:41:44Z",
          "tree_id": "61e3a167309c1df2e1a176e2ec762919862e6ec0",
          "url": "https://github.com/galacticusorg/galacticus/commit/94ded2ded9c80fed323bfe76c2075977081ab394"
        },
        "date": 1693997325126,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 146.723,
            "unit": "seconds",
            "range": 0.211565828998024
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "94ded2ded9c80fed323bfe76c2075977081ab394",
          "message": "Merge pull request #453 from galacticusorg/optimization\n\nBaryonic physics optimizations",
          "timestamp": "2023-09-05T14:41:44Z",
          "tree_id": "61e3a167309c1df2e1a176e2ec762919862e6ec0",
          "url": "https://github.com/galacticusorg/galacticus/commit/94ded2ded9c80fed323bfe76c2075977081ab394"
        },
        "date": 1693997332830,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.582563003824,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.94176907891703,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.925672297704321,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2309698659553,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 72.0108989507493,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.20979088741696,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "92dfec05a324db5178bccd47285f1b4899ea4cbd",
          "message": "Merge pull request #464 from galacticusorg/sigmaTruncatedPower\n\nMake integration of $\\sigma(M)$ more robust for models with truncated power spectra",
          "timestamp": "2023-09-08T01:57:11Z",
          "tree_id": "ca3022418aa99d027443e8d34a61e8ab6c5ce4c1",
          "url": "https://github.com/galacticusorg/galacticus/commit/92dfec05a324db5178bccd47285f1b4899ea4cbd"
        },
        "date": 1694161693575,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 215.158,
            "unit": "seconds",
            "range": 0.285001052631738
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "92dfec05a324db5178bccd47285f1b4899ea4cbd",
          "message": "Merge pull request #464 from galacticusorg/sigmaTruncatedPower\n\nMake integration of $\\sigma(M)$ more robust for models with truncated power spectra",
          "timestamp": "2023-09-08T01:57:11Z",
          "tree_id": "ca3022418aa99d027443e8d34a61e8ab6c5ce4c1",
          "url": "https://github.com/galacticusorg/galacticus/commit/92dfec05a324db5178bccd47285f1b4899ea4cbd"
        },
        "date": 1694161702617,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.68914877952759,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.9172094282908,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.00886518625528,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2241451172097,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 71.4592815741363,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.23629461084522,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2e41e19b6199e6b4db4d9755326f6a5ee4f80e55",
          "message": "Merge pull request #463 from galacticusorg/roundOffFix\n\nAllow to tolerate round-off errors in merger tree branching integrals",
          "timestamp": "2023-09-08T14:33:45Z",
          "tree_id": "096fc4774a16a0df7d09de52ea767bf336602713",
          "url": "https://github.com/galacticusorg/galacticus/commit/2e41e19b6199e6b4db4d9755326f6a5ee4f80e55"
        },
        "date": 1694196574313,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 207.16,
            "unit": "seconds",
            "range": 0.179955550060414
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2e41e19b6199e6b4db4d9755326f6a5ee4f80e55",
          "message": "Merge pull request #463 from galacticusorg/roundOffFix\n\nAllow to tolerate round-off errors in merger tree branching integrals",
          "timestamp": "2023-09-08T14:33:45Z",
          "tree_id": "096fc4774a16a0df7d09de52ea767bf336602713",
          "url": "https://github.com/galacticusorg/galacticus/commit/2e41e19b6199e6b4db4d9755326f6a5ee4f80e55"
        },
        "date": 1694196582308,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.17454457340027,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.46153940840019,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.06095641891005,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.252661075156,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 70.991499069039,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 9.77241950199985,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "219789d07e9fd373d3150fb691d04a9659f34f9f",
          "message": "Merge pull request #465 from galacticusorg/fspsArm64\n\nFix build of FSPS on Apple M1/M2 chips",
          "timestamp": "2023-09-09T02:39:08Z",
          "tree_id": "230759199a9f1f1d41b627a1320054d4ce3bc165",
          "url": "https://github.com/galacticusorg/galacticus/commit/219789d07e9fd373d3150fb691d04a9659f34f9f"
        },
        "date": 1694252793962,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 168.812,
            "unit": "seconds",
            "range": 0.136724540589684
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "219789d07e9fd373d3150fb691d04a9659f34f9f",
          "message": "Merge pull request #465 from galacticusorg/fspsArm64\n\nFix build of FSPS on Apple M1/M2 chips",
          "timestamp": "2023-09-09T02:39:08Z",
          "tree_id": "230759199a9f1f1d41b627a1320054d4ce3bc165",
          "url": "https://github.com/galacticusorg/galacticus/commit/219789d07e9fd373d3150fb691d04a9659f34f9f"
        },
        "date": 1694252800929,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.17812392253161,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.59176895706195,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.3520004131317,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3002661797178,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 72.8505659807896,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.7651717039636,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "051be396b0a919add84d46a5b85b0254c7b43136",
          "message": "Merge pull request #466 from galacticusorg/optimization\n\nFurther optimization of baryonic physics models",
          "timestamp": "2023-09-09T16:47:29Z",
          "tree_id": "ac23c59150536c00e9d5ae80e78c411db8cac8b5",
          "url": "https://github.com/galacticusorg/galacticus/commit/051be396b0a919add84d46a5b85b0254c7b43136"
        },
        "date": 1694288571929,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 230.926,
            "unit": "seconds",
            "range": 0.272738702788789
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "051be396b0a919add84d46a5b85b0254c7b43136",
          "message": "Merge pull request #466 from galacticusorg/optimization\n\nFurther optimization of baryonic physics models",
          "timestamp": "2023-09-09T16:47:29Z",
          "tree_id": "ac23c59150536c00e9d5ae80e78c411db8cac8b5",
          "url": "https://github.com/galacticusorg/galacticus/commit/051be396b0a919add84d46a5b85b0254c7b43136"
        },
        "date": 1694288579521,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.49781059521228,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.95039293803183,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.09414061187602,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.4514447057735,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 94.8235487034022,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.25604926997866,
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
          "id": "3c616563c904639fb0e8ea211b6a1d5173ce1518",
          "message": "fix: Correct method names",
          "timestamp": "2023-09-14T04:19:12Z",
          "tree_id": "7e1d4cb2c352a6c6ca097763944a6c69a4d3c68c",
          "url": "https://github.com/galacticusorg/galacticus/commit/3c616563c904639fb0e8ea211b6a1d5173ce1518"
        },
        "date": 1694686074202,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 206.73,
            "unit": "seconds",
            "range": 0.344220859332119
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
          "id": "3c616563c904639fb0e8ea211b6a1d5173ce1518",
          "message": "fix: Correct method names",
          "timestamp": "2023-09-14T04:19:12Z",
          "tree_id": "7e1d4cb2c352a6c6ca097763944a6c69a4d3c68c",
          "url": "https://github.com/galacticusorg/galacticus/commit/3c616563c904639fb0e8ea211b6a1d5173ce1518"
        },
        "date": 1694686082748,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.62502284588273,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.89998486676598,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.906402852432514,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.243141878245,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 76.3364866012826,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.32204863767182,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "90736f9fff69acdc6c5e893519c274fb2f0260df",
          "message": "Merge pull request #470 from galacticusorg/arrayPropertyColumns\n\nAllow `nodePropertyExtractorArray` objects to specify numerical column values",
          "timestamp": "2023-09-14T20:51:32Z",
          "tree_id": "f112f521e66397511c8a2b2c1d252c2b5415125f",
          "url": "https://github.com/galacticusorg/galacticus/commit/90736f9fff69acdc6c5e893519c274fb2f0260df"
        },
        "date": 1694735802549,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 266.26,
            "unit": "seconds",
            "range": 1.7358663543032
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "90736f9fff69acdc6c5e893519c274fb2f0260df",
          "message": "Merge pull request #470 from galacticusorg/arrayPropertyColumns\n\nAllow `nodePropertyExtractorArray` objects to specify numerical column values",
          "timestamp": "2023-09-14T20:51:32Z",
          "tree_id": "f112f521e66397511c8a2b2c1d252c2b5415125f",
          "url": "https://github.com/galacticusorg/galacticus/commit/90736f9fff69acdc6c5e893519c274fb2f0260df"
        },
        "date": 1694735811596,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.57944474007854,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.0930259767054,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.926652115549488,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2724648032126,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 94.9665517712776,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.24127618049523,
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
          "id": "c02ffb64269e7c03a3471b6a8761d7833668507c",
          "message": "fix: Correct enclosed mass in dimensionfull Hernquist profiles",
          "timestamp": "2023-09-15T18:33:55Z",
          "tree_id": "9c722a35a7ee8662efccdf66ce4dfdc656a54ed0",
          "url": "https://github.com/galacticusorg/galacticus/commit/c02ffb64269e7c03a3471b6a8761d7833668507c"
        },
        "date": 1694813410007,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 247.297,
            "unit": "seconds",
            "range": 0.453488809123873
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
          "id": "c02ffb64269e7c03a3471b6a8761d7833668507c",
          "message": "fix: Correct enclosed mass in dimensionfull Hernquist profiles",
          "timestamp": "2023-09-15T18:33:55Z",
          "tree_id": "9c722a35a7ee8662efccdf66ce4dfdc656a54ed0",
          "url": "https://github.com/galacticusorg/galacticus/commit/c02ffb64269e7c03a3471b6a8761d7833668507c"
        },
        "date": 1694813417955,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.62502284588273,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.89998486676598,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.906402852432514,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.243141878245,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 76.3364866012826,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.32204863767182,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2c11769768762c96694def43cef28134e9856e9f",
          "message": "Merge pull request #475 from galacticusorg/singleStepLabels\n\nAdd labeling of progenitor halo origins when using the `singleStep` build controller",
          "timestamp": "2023-09-18T23:02:27Z",
          "tree_id": "724e53c9d7489384be374b7f6b6c4c3c5e09e3c3",
          "url": "https://github.com/galacticusorg/galacticus/commit/2c11769768762c96694def43cef28134e9856e9f"
        },
        "date": 1695091774588,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 214.385,
            "unit": "seconds",
            "range": 0.737146186315819
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2c11769768762c96694def43cef28134e9856e9f",
          "message": "Merge pull request #475 from galacticusorg/singleStepLabels\n\nAdd labeling of progenitor halo origins when using the `singleStep` build controller",
          "timestamp": "2023-09-18T23:02:27Z",
          "tree_id": "724e53c9d7489384be374b7f6b6c4c3c5e09e3c3",
          "url": "https://github.com/galacticusorg/galacticus/commit/2c11769768762c96694def43cef28134e9856e9f"
        },
        "date": 1695091783443,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.51945661377284,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.83484408288879,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.0769978886767,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2270217525507,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 71.5883879563574,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.25224901182828,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "40d12d74dfc67ddd41f69404da885368cdc235a6",
          "message": "Merge pull request #476 from galacticusorg/optimization\n\nFurther minor optimizations",
          "timestamp": "2023-09-19T02:55:55Z",
          "tree_id": "b53625ec22535b8ae414099fdd31bee401d6ac8e",
          "url": "https://github.com/galacticusorg/galacticus/commit/40d12d74dfc67ddd41f69404da885368cdc235a6"
        },
        "date": 1695103182774,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 261.856,
            "unit": "seconds",
            "range": 0.36417633091625
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "40d12d74dfc67ddd41f69404da885368cdc235a6",
          "message": "Merge pull request #476 from galacticusorg/optimization\n\nFurther minor optimizations",
          "timestamp": "2023-09-19T02:55:55Z",
          "tree_id": "b53625ec22535b8ae414099fdd31bee401d6ac8e",
          "url": "https://github.com/galacticusorg/galacticus/commit/40d12d74dfc67ddd41f69404da885368cdc235a6"
        },
        "date": 1695103191179,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.57944474007854,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.0930259767054,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.926652115549488,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2724648032126,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 94.9665517712776,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.24127618049523,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "736fa7d2b150b43bc23cdb67ef92b5366df06cee",
          "message": "Merge pull request #477 from galacticusorg/optimization\n\nGenerate custom equality operator for each enumeration",
          "timestamp": "2023-09-20T23:11:39Z",
          "tree_id": "bd9ad831db6bd6e1f9169848199172d24112190c",
          "url": "https://github.com/galacticusorg/galacticus/commit/736fa7d2b150b43bc23cdb67ef92b5366df06cee"
        },
        "date": 1695263403110,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 250.36,
            "unit": "seconds",
            "range": 0.198982411278064
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "736fa7d2b150b43bc23cdb67ef92b5366df06cee",
          "message": "Merge pull request #477 from galacticusorg/optimization\n\nGenerate custom equality operator for each enumeration",
          "timestamp": "2023-09-20T23:11:39Z",
          "tree_id": "bd9ad831db6bd6e1f9169848199172d24112190c",
          "url": "https://github.com/galacticusorg/galacticus/commit/736fa7d2b150b43bc23cdb67ef92b5366df06cee"
        },
        "date": 1695263410540,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.36461935373672,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.7151401559154,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.16008929828541,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9466564713043,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 69.3818874703465,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 9.9418741682224,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "7929a8eaaa3b36fee72cb343abe04644cbf3bc22",
          "message": "Merge pull request #478 from galacticusorg/fixNBodyOperator\n\nFix nBody operator",
          "timestamp": "2023-09-21T14:19:42Z",
          "tree_id": "a6dc53c9df6c0035b2555999c50aaa6393d44688",
          "url": "https://github.com/galacticusorg/galacticus/commit/7929a8eaaa3b36fee72cb343abe04644cbf3bc22"
        },
        "date": 1695316883941,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 236.779,
            "unit": "seconds",
            "range": 1.89409685602436
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "7929a8eaaa3b36fee72cb343abe04644cbf3bc22",
          "message": "Merge pull request #478 from galacticusorg/fixNBodyOperator\n\nFix nBody operator",
          "timestamp": "2023-09-21T14:19:42Z",
          "tree_id": "a6dc53c9df6c0035b2555999c50aaa6393d44688",
          "url": "https://github.com/galacticusorg/galacticus/commit/7929a8eaaa3b36fee72cb343abe04644cbf3bc22"
        },
        "date": 1695316894400,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.08604144839345,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.6182604290444,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.924168114205167,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6485312397504,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 66.0367708534659,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 9.63167271522589,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "54c7f9118c0af37e5841f387e790c6049de6dcc7",
          "message": "Merge pull request #479 from galacticusorg/noninstantFix\n\nCopy element index in to OpenMP parallel region",
          "timestamp": "2023-09-26T01:11:06Z",
          "tree_id": "c26dc4bde88c09d9cde120ca2bfadd947a0e07b8",
          "url": "https://github.com/galacticusorg/galacticus/commit/54c7f9118c0af37e5841f387e790c6049de6dcc7"
        },
        "date": 1695712619182,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 206.104,
            "unit": "seconds",
            "range": 0.124436329105092
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "54c7f9118c0af37e5841f387e790c6049de6dcc7",
          "message": "Merge pull request #479 from galacticusorg/noninstantFix\n\nCopy element index in to OpenMP parallel region",
          "timestamp": "2023-09-26T01:11:06Z",
          "tree_id": "c26dc4bde88c09d9cde120ca2bfadd947a0e07b8",
          "url": "https://github.com/galacticusorg/galacticus/commit/54c7f9118c0af37e5841f387e790c6049de6dcc7"
        },
        "date": 1695712627402,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.08016796083313,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.5375848481544,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.939221212504722,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6397508540005,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 66.9889662968008,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 9.59927968737974,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "ceb37911054a39d4fffe7855ba9050094ddab31f",
          "message": "Merge pull request #480 from galacticusorg/sourceDigestFix\n\nEnsure source digests are updated",
          "timestamp": "2023-09-26T13:16:37Z",
          "tree_id": "f9663af3d536968c3c5881e78426ef2979a120df",
          "url": "https://github.com/galacticusorg/galacticus/commit/ceb37911054a39d4fffe7855ba9050094ddab31f"
        },
        "date": 1695745010689,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 209.643,
            "unit": "seconds",
            "range": 0.453449115117419
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "ceb37911054a39d4fffe7855ba9050094ddab31f",
          "message": "Merge pull request #480 from galacticusorg/sourceDigestFix\n\nEnsure source digests are updated",
          "timestamp": "2023-09-26T13:16:37Z",
          "tree_id": "f9663af3d536968c3c5881e78426ef2979a120df",
          "url": "https://github.com/galacticusorg/galacticus/commit/ceb37911054a39d4fffe7855ba9050094ddab31f"
        },
        "date": 1695745020206,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 2.9594186480383,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.69418556042769,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.29144025280822,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9049345411609,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 67.8964708995721,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.6562052146004,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "456ee65559ae2d66240a0338a9207d59a749d28e",
          "message": "Merge pull request #482 from galacticusorg/bryanNorman1998GeneralCosmology\n\nAllow the `virialDensityContrastBryanNorman1998` class to work for arbitrary cosmologies",
          "timestamp": "2023-09-27T23:37:00Z",
          "tree_id": "6123b72aff0f3c5e0099a26b1c6499cecaacd08e",
          "url": "https://github.com/galacticusorg/galacticus/commit/456ee65559ae2d66240a0338a9207d59a749d28e"
        },
        "date": 1695869201233,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 228.507,
            "unit": "seconds",
            "range": 0.234230869016282
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "456ee65559ae2d66240a0338a9207d59a749d28e",
          "message": "Merge pull request #482 from galacticusorg/bryanNorman1998GeneralCosmology\n\nAllow the `virialDensityContrastBryanNorman1998` class to work for arbitrary cosmologies",
          "timestamp": "2023-09-27T23:37:00Z",
          "tree_id": "6123b72aff0f3c5e0099a26b1c6499cecaacd08e",
          "url": "https://github.com/galacticusorg/galacticus/commit/456ee65559ae2d66240a0338a9207d59a749d28e"
        },
        "date": 1695869209790,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.26508433888745,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.3101047357318,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.931769699117319,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6269493362323,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 65.3346334051298,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 9.79492758041714,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "563a05278a954caba9f81b617d8bece9800c5cb4",
          "message": "Merge pull request #483 from cgannonucm/tidally_truncated_improvement\n\nImprovements to the Radial Range that Subhalo Density Profiles are Tabulated when Fitting TNFW Profile",
          "timestamp": "2023-09-28T19:11:30Z",
          "tree_id": "869f3c09f051317a0b634f10f933ea66c1032e7f",
          "url": "https://github.com/galacticusorg/galacticus/commit/563a05278a954caba9f81b617d8bece9800c5cb4"
        },
        "date": 1695941957987,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 242.19,
            "unit": "seconds",
            "range": 0.457628670430354
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "563a05278a954caba9f81b617d8bece9800c5cb4",
          "message": "Merge pull request #483 from cgannonucm/tidally_truncated_improvement\n\nImprovements to the Radial Range that Subhalo Density Profiles are Tabulated when Fitting TNFW Profile",
          "timestamp": "2023-09-28T19:11:30Z",
          "tree_id": "869f3c09f051317a0b634f10f933ea66c1032e7f",
          "url": "https://github.com/galacticusorg/galacticus/commit/563a05278a954caba9f81b617d8bece9800c5cb4"
        },
        "date": 1695941966017,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 2.95014513235564,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.66026481258829,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.1512787904706,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3150443734923,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 72.1724829562187,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.7592151656849,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "82f64062a3b79f48ea8ee6fb5282436b0c155e8f",
          "message": "Merge pull request #487 from galacticusorg/brownianBridgeFix\n\nFix mathematical implementation of Brownian bridge solutions",
          "timestamp": "2023-10-03T04:56:09Z",
          "tree_id": "47de2ca5c923cedf8aeeb044c9ea346fdf77bfa0",
          "url": "https://github.com/galacticusorg/galacticus/commit/82f64062a3b79f48ea8ee6fb5282436b0c155e8f"
        },
        "date": 1696320726513,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 214.053,
            "unit": "seconds",
            "range": 0.2104331247647
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "82f64062a3b79f48ea8ee6fb5282436b0c155e8f",
          "message": "Merge pull request #487 from galacticusorg/brownianBridgeFix\n\nFix mathematical implementation of Brownian bridge solutions",
          "timestamp": "2023-10-03T04:56:09Z",
          "tree_id": "47de2ca5c923cedf8aeeb044c9ea346fdf77bfa0",
          "url": "https://github.com/galacticusorg/galacticus/commit/82f64062a3b79f48ea8ee6fb5282436b0c155e8f"
        },
        "date": 1696320734109,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 2.99585087369573,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.42399594307135,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.13606998291136,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2014044103602,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 74.4316183630359,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.6609219635255,
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
          "id": "7452d78a236da9cf40c3383eeb0b81f6100669f9",
          "message": "fix: Remove unused variable",
          "timestamp": "2023-10-03T14:34:37Z",
          "tree_id": "59998fcf7f59e58f41e70068930b4c8c04085fd5",
          "url": "https://github.com/galacticusorg/galacticus/commit/7452d78a236da9cf40c3383eeb0b81f6100669f9"
        },
        "date": 1696368219045,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 218.324,
            "unit": "seconds",
            "range": 0.376417852922867
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
          "id": "7452d78a236da9cf40c3383eeb0b81f6100669f9",
          "message": "fix: Remove unused variable",
          "timestamp": "2023-10-03T14:34:37Z",
          "tree_id": "59998fcf7f59e58f41e70068930b4c8c04085fd5",
          "url": "https://github.com/galacticusorg/galacticus/commit/7452d78a236da9cf40c3383eeb0b81f6100669f9"
        },
        "date": 1696368226662,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 2.99585087369573,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.42399594307135,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.13606998291136,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2014044103602,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 74.4316183630359,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.6609219635255,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "017bfd229dbee969fbeb24f6ca05e52cc054c55b",
          "message": "Merge pull request #488 from galacticusorg/noninstantTests\n\nAdd functionality useful for testing the non-instantaneous feedback model",
          "timestamp": "2023-10-04T04:12:14Z",
          "tree_id": "149548972a67effcde5cd1b7b851cdbfe0ad3847",
          "url": "https://github.com/galacticusorg/galacticus/commit/017bfd229dbee969fbeb24f6ca05e52cc054c55b"
        },
        "date": 1696404655956,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 308.341,
            "unit": "seconds",
            "range": 1.13894727709364
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "017bfd229dbee969fbeb24f6ca05e52cc054c55b",
          "message": "Merge pull request #488 from galacticusorg/noninstantTests\n\nAdd functionality useful for testing the non-instantaneous feedback model",
          "timestamp": "2023-10-04T04:12:14Z",
          "tree_id": "149548972a67effcde5cd1b7b851cdbfe0ad3847",
          "url": "https://github.com/galacticusorg/galacticus/commit/017bfd229dbee969fbeb24f6ca05e52cc054c55b"
        },
        "date": 1696404664234,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.40710430543489,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.07923576824673,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.06931789280255,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1982598706471,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 77.3349617019032,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.33281194103275,
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
          "id": "554b3cffba687d80957ed7348f8d5d791b8f807e",
          "message": "feat: Add a feedback outflow class that limits to some maximum mass loading factor\n\nThe limit is approach asymptotically using a `tanh()` function.",
          "timestamp": "2023-10-04T17:07:40Z",
          "tree_id": "9133eee86d152c5720417d784448286a665fbaa5",
          "url": "https://github.com/galacticusorg/galacticus/commit/554b3cffba687d80957ed7348f8d5d791b8f807e"
        },
        "date": 1696450555115,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 258.755,
            "unit": "seconds",
            "range": 0.130569904652507
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
          "id": "554b3cffba687d80957ed7348f8d5d791b8f807e",
          "message": "feat: Add a feedback outflow class that limits to some maximum mass loading factor\n\nThe limit is approach asymptotically using a `tanh()` function.",
          "timestamp": "2023-10-04T17:07:40Z",
          "tree_id": "9133eee86d152c5720417d784448286a665fbaa5",
          "url": "https://github.com/galacticusorg/galacticus/commit/554b3cffba687d80957ed7348f8d5d791b8f807e"
        },
        "date": 1696450563133,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.70885609720588,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 9.21600081932495,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.824276083210193,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3086093549023,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 73.3381313395941,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 1.27646139745279,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d6ca028811ec022b51a2ae1c0b21d3538ac56a06",
          "message": "Merge pull request #490 from galacticusorg/cgmCoolingFunction\n\nAdd a `nodePropertyExtractor` to output the cooling function in the CGM",
          "timestamp": "2023-10-05T14:30:53Z",
          "tree_id": "5ed79355033c2e09935deb3ef720d1bca878712d",
          "url": "https://github.com/galacticusorg/galacticus/commit/d6ca028811ec022b51a2ae1c0b21d3538ac56a06"
        },
        "date": 1696536321771,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 193.685,
            "unit": "seconds",
            "range": 0.272507798052126
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d6ca028811ec022b51a2ae1c0b21d3538ac56a06",
          "message": "Merge pull request #490 from galacticusorg/cgmCoolingFunction\n\nAdd a `nodePropertyExtractor` to output the cooling function in the CGM",
          "timestamp": "2023-10-05T14:30:53Z",
          "tree_id": "5ed79355033c2e09935deb3ef720d1bca878712d",
          "url": "https://github.com/galacticusorg/galacticus/commit/d6ca028811ec022b51a2ae1c0b21d3538ac56a06"
        },
        "date": 1696536329538,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 2.9594186480383,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.69418556042769,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.29144025280822,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9049345411609,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 67.8964708995721,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.6562052146004,
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
          "id": "7f2cf528e741255fb62f08f1db1a2ae0cce0a131",
          "message": "feat: Add functionality to output the (infimum of the) excursion along each branch\n\nAdds a `nodeOperator` and `nodePropertyExtractor` to output these.",
          "timestamp": "2023-10-06T18:53:55Z",
          "tree_id": "8b71a16284f768e8a225f4d38a676f11cf7ec5ae",
          "url": "https://github.com/galacticusorg/galacticus/commit/7f2cf528e741255fb62f08f1db1a2ae0cce0a131"
        },
        "date": 1696646072863,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 204.815,
            "unit": "seconds",
            "range": 0.206558708361339
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
          "id": "7f2cf528e741255fb62f08f1db1a2ae0cce0a131",
          "message": "feat: Add functionality to output the (infimum of the) excursion along each branch\n\nAdds a `nodeOperator` and `nodePropertyExtractor` to output these.",
          "timestamp": "2023-10-06T18:53:55Z",
          "tree_id": "8b71a16284f768e8a225f4d38a676f11cf7ec5ae",
          "url": "https://github.com/galacticusorg/galacticus/commit/7f2cf528e741255fb62f08f1db1a2ae0cce0a131"
        },
        "date": 1696646079872,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.04384072964081,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.3107827908493,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.11629035160799,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.5024823914336,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 72.9704556594375,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.4067870776505,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "4a1abebecd4065bde5191c8d37f61bc3c59e80dd",
          "message": "Merge pull request #492 from galacticusorg/memoryLeakFixes\n\nFix memory leaks",
          "timestamp": "2023-10-07T14:52:07Z",
          "tree_id": "b693cd8daca4d7bb18cb81b746e6eb13ce4b55c8",
          "url": "https://github.com/galacticusorg/galacticus/commit/4a1abebecd4065bde5191c8d37f61bc3c59e80dd"
        },
        "date": 1696702253996,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 219.335,
            "unit": "seconds",
            "range": 0.179756780125228
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "4a1abebecd4065bde5191c8d37f61bc3c59e80dd",
          "message": "Merge pull request #492 from galacticusorg/memoryLeakFixes\n\nFix memory leaks",
          "timestamp": "2023-10-07T14:52:07Z",
          "tree_id": "b693cd8daca4d7bb18cb81b746e6eb13ce4b55c8",
          "url": "https://github.com/galacticusorg/galacticus/commit/4a1abebecd4065bde5191c8d37f61bc3c59e80dd"
        },
        "date": 1696702261974,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.36461935373672,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.7151401559154,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.16008929828541,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9466564713043,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 69.3818874703465,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 9.9418741682224,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a29225c6ef9c5f69db7be1d84076f6272e863686",
          "message": "Merge pull request #493 from sachiwee/master\n\nChanges to accretion.halo.cold_mode.F90 and cooling.cooling_function.molecular_hydrogen_Galli_Palla.F90",
          "timestamp": "2023-10-10T12:46:37Z",
          "tree_id": "e512ab8f80b76a75a9c79caf1413ccd688e3310c",
          "url": "https://github.com/galacticusorg/galacticus/commit/a29225c6ef9c5f69db7be1d84076f6272e863686"
        },
        "date": 1696954591436,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 281.081,
            "unit": "seconds",
            "range": 1.213199447741
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a29225c6ef9c5f69db7be1d84076f6272e863686",
          "message": "Merge pull request #493 from sachiwee/master\n\nChanges to accretion.halo.cold_mode.F90 and cooling.cooling_function.molecular_hydrogen_Galli_Palla.F90",
          "timestamp": "2023-10-10T12:46:37Z",
          "tree_id": "e512ab8f80b76a75a9c79caf1413ccd688e3310c",
          "url": "https://github.com/galacticusorg/galacticus/commit/a29225c6ef9c5f69db7be1d84076f6272e863686"
        },
        "date": 1696954600476,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 2.95014513235564,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.66026481258829,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.1512787904706,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3150443734923,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 72.1724829562187,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.7592151656849,
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
          "id": "58aa67136b6a741ddc2363338b604055acad476d",
          "message": "fix: Update URL",
          "timestamp": "2023-10-11T06:47:47-07:00",
          "tree_id": "49bd2f41506894b7081486951d899277fd45c118",
          "url": "https://github.com/galacticusorg/galacticus/commit/58aa67136b6a741ddc2363338b604055acad476d"
        },
        "date": 1697043309012,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 194.698,
            "unit": "seconds",
            "range": 0.12188355097913
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
          "id": "58aa67136b6a741ddc2363338b604055acad476d",
          "message": "fix: Update URL",
          "timestamp": "2023-10-11T06:47:47-07:00",
          "tree_id": "49bd2f41506894b7081486951d899277fd45c118",
          "url": "https://github.com/galacticusorg/galacticus/commit/58aa67136b6a741ddc2363338b604055acad476d"
        },
        "date": 1697043317108,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 2.99585087369573,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 8.42399594307135,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.13606998291136,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2014044103602,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 74.4316183630359,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 10.6609219635255,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6ea5818074995f564179acdc3d43561338cbd8c7",
          "message": "Merge pull request #497 from galacticusorg/oldCloudy\n\nAllow download of Cloudy from `old/` subdirectory if main directory fails",
          "timestamp": "2023-10-19T00:13:35Z",
          "tree_id": "83aa356f85ec6cfbe93a07b919f5b4db8dd08100",
          "url": "https://github.com/galacticusorg/galacticus/commit/6ea5818074995f564179acdc3d43561338cbd8c7"
        },
        "date": 1697731205471,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 226.97,
            "unit": "seconds",
            "range": 0.0771880819869763
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6ea5818074995f564179acdc3d43561338cbd8c7",
          "message": "Merge pull request #497 from galacticusorg/oldCloudy\n\nAllow download of Cloudy from `old/` subdirectory if main directory fails",
          "timestamp": "2023-10-19T00:13:35Z",
          "tree_id": "83aa356f85ec6cfbe93a07b919f5b4db8dd08100",
          "url": "https://github.com/galacticusorg/galacticus/commit/6ea5818074995f564179acdc3d43561338cbd8c7"
        },
        "date": 1697731213347,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.07440025406794,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 10.5338743554984,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.950760056531201,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.604316850772,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 66.2044976862593,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 9.54165136000034,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "434e3b1a99e6580da8e60b4cae63b6c3f8a5db14",
          "message": "Merge pull request #495 from galacticusorg/optimization\n\nFurther optimization improvements",
          "timestamp": "2023-10-19T16:03:13Z",
          "tree_id": "c3abd7d712eb06c310a40d5ac9782f57693bf3d2",
          "url": "https://github.com/galacticusorg/galacticus/commit/434e3b1a99e6580da8e60b4cae63b6c3f8a5db14"
        },
        "date": 1697742763542,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 234.309,
            "unit": "seconds",
            "range": 0.622494096359002
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "434e3b1a99e6580da8e60b4cae63b6c3f8a5db14",
          "message": "Merge pull request #495 from galacticusorg/optimization\n\nFurther optimization improvements",
          "timestamp": "2023-10-19T16:03:13Z",
          "tree_id": "c3abd7d712eb06c310a40d5ac9782f57693bf3d2",
          "url": "https://github.com/galacticusorg/galacticus/commit/434e3b1a99e6580da8e60b4cae63b6c3f8a5db14"
        },
        "date": 1697742771081,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.23841396448377,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 3.9858051834278,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.19313032958238,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0809392929276,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 86.0986188134057,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.59177985821157,
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
          "id": "d3ea3e0ccd0311908a925b7933e1adb1cb06ceea",
          "message": "fix: Catch unphysical helium ionizing luminosities when computing emission line luminosities",
          "timestamp": "2023-10-24T09:36:14-07:00",
          "tree_id": "2f1e0cd50c64a76f3b651b80df4109e3ca54e7fc",
          "url": "https://github.com/galacticusorg/galacticus/commit/d3ea3e0ccd0311908a925b7933e1adb1cb06ceea"
        },
        "date": 1698177440537,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 227.608,
            "unit": "seconds",
            "range": 0.280331232650763
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
          "id": "d3ea3e0ccd0311908a925b7933e1adb1cb06ceea",
          "message": "fix: Catch unphysical helium ionizing luminosities when computing emission line luminosities",
          "timestamp": "2023-10-24T09:36:14-07:00",
          "tree_id": "2f1e0cd50c64a76f3b651b80df4109e3ca54e7fc",
          "url": "https://github.com/galacticusorg/galacticus/commit/d3ea3e0ccd0311908a925b7933e1adb1cb06ceea"
        },
        "date": 1698177448293,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.78832415645585,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.53260821574734,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.981611288705017,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0949507391994,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 85.0013535775173,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.48823833888572,
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
          "id": "a09c88e66393e68cbb03588ca637e41f5468fd02",
          "message": "feat: Update to use the new standard stellar properties compilation file",
          "timestamp": "2023-10-25T09:49:56-07:00",
          "tree_id": "ccd6ce8b65c9969b45d8d3f9181f717a656f9eba",
          "url": "https://github.com/galacticusorg/galacticus/commit/a09c88e66393e68cbb03588ca637e41f5468fd02"
        },
        "date": 1698265230368,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 232.965,
            "unit": "seconds",
            "range": 2.97990108896267
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
          "id": "a09c88e66393e68cbb03588ca637e41f5468fd02",
          "message": "feat: Update to use the new standard stellar properties compilation file",
          "timestamp": "2023-10-25T09:49:56-07:00",
          "tree_id": "ccd6ce8b65c9969b45d8d3f9181f717a656f9eba",
          "url": "https://github.com/galacticusorg/galacticus/commit/a09c88e66393e68cbb03588ca637e41f5468fd02"
        },
        "date": 1698265237616,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.75852871325266,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.58723611000534,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.898354170783437,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9605682277289,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 88.6527098733839,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.50287193001514,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "4ddb03adabf45de3bb5c89c7cb693bc68a546c2d",
          "message": "Merge pull request #501 from galacticusorg/sfhNodeDataOutput\n\nRefactor star formation histories to output to the `nodeData` group",
          "timestamp": "2023-10-29T14:35:51Z",
          "tree_id": "ea57bcf4d22d255988a2928904548c0550e2bd55",
          "url": "https://github.com/galacticusorg/galacticus/commit/4ddb03adabf45de3bb5c89c7cb693bc68a546c2d"
        },
        "date": 1698607997736,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 229.987,
            "unit": "seconds",
            "range": 0.349024497703918
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "4ddb03adabf45de3bb5c89c7cb693bc68a546c2d",
          "message": "Merge pull request #501 from galacticusorg/sfhNodeDataOutput\n\nRefactor star formation histories to output to the `nodeData` group",
          "timestamp": "2023-10-29T14:35:51Z",
          "tree_id": "ea57bcf4d22d255988a2928904548c0550e2bd55",
          "url": "https://github.com/galacticusorg/galacticus/commit/4ddb03adabf45de3bb5c89c7cb693bc68a546c2d"
        },
        "date": 1698608007285,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.2710797331703,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.33755238090893,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.17002349570196,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.04828315788,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 87.1728345534998,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.7879671075113,
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
          "id": "f9e67cd1dde359d0864c04813883ea8618ad1bc4",
          "message": "Merge pull request #502 from galacticusorg/dependabot/github_actions/tj-actions/changed-files-40\n\nbuild(deps): bump tj-actions/changed-files from 39 to 40",
          "timestamp": "2023-10-30T14:28:00Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/f9e67cd1dde359d0864c04813883ea8618ad1bc4"
        },
        "date": 1699476974752,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 145.618,
            "unit": "seconds",
            "range": 0.15230101772584
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
          "id": "f9e67cd1dde359d0864c04813883ea8618ad1bc4",
          "message": "Merge pull request #502 from galacticusorg/dependabot/github_actions/tj-actions/changed-files-40\n\nbuild(deps): bump tj-actions/changed-files from 39 to 40",
          "timestamp": "2023-10-30T14:28:00Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/f9e67cd1dde359d0864c04813883ea8618ad1bc4"
        },
        "date": 1699476983168,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.28178078755925,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.22702360097909,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.02674833069628,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1705993153765,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 86.3783487358013,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.7241140957634,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d9ea36be8fbdd4c4cf57cbcbb1ba10913aea780a",
          "message": "Merge pull request #503 from galacticusorg/perl38RegexFix\n\nFix broken regex behavior",
          "timestamp": "2023-11-09T20:53:29Z",
          "tree_id": "4a6ea660a42cf563a5c859aadda1a70c1f3869f2",
          "url": "https://github.com/galacticusorg/galacticus/commit/d9ea36be8fbdd4c4cf57cbcbb1ba10913aea780a"
        },
        "date": 1699572573604,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 271.328,
            "unit": "seconds",
            "range": 0.873920820209423
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d9ea36be8fbdd4c4cf57cbcbb1ba10913aea780a",
          "message": "Merge pull request #503 from galacticusorg/perl38RegexFix\n\nFix broken regex behavior",
          "timestamp": "2023-11-09T20:53:29Z",
          "tree_id": "4a6ea660a42cf563a5c859aadda1a70c1f3869f2",
          "url": "https://github.com/galacticusorg/galacticus/commit/d9ea36be8fbdd4c4cf57cbcbb1ba10913aea780a"
        },
        "date": 1699572582289,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.17508832542725,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.24824129370786,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.14229266604171,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1381068077174,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 121.211449594275,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.7561817326953,
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
          "id": "f98787d2300cfa787d4b55ee24f4d972812bd32b",
          "message": "feat: Check that `satellite` component `position` and `velocity` are gettable",
          "timestamp": "2023-11-13T18:18:34-08:00",
          "tree_id": "aef7b94cceda9cbbfd61f452ca6fc5b1c7cd34b3",
          "url": "https://github.com/galacticusorg/galacticus/commit/f98787d2300cfa787d4b55ee24f4d972812bd32b"
        },
        "date": 1699937461648,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 180.552,
            "unit": "seconds",
            "range": 0.266780059225807
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
          "id": "f98787d2300cfa787d4b55ee24f4d972812bd32b",
          "message": "feat: Check that `satellite` component `position` and `velocity` are gettable",
          "timestamp": "2023-11-13T18:18:34-08:00",
          "tree_id": "aef7b94cceda9cbbfd61f452ca6fc5b1c7cd34b3",
          "url": "https://github.com/galacticusorg/galacticus/commit/f98787d2300cfa787d4b55ee24f4d972812bd32b"
        },
        "date": 1699937467751,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.00893403619692,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.50114915636612,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.71009698372859,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.442930589529,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 285.450796641662,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 2.6943288509175,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "150e6e81fa667e06ab5d1db655c8a6bfea150199",
          "message": "Merge pull request #505 from galacticusorg/smbhPropertyOutput\n\nRefactor SMBH property output such that these are output to the `nodeData` group",
          "timestamp": "2023-11-17T22:22:56Z",
          "tree_id": "c1538aa88d1a2df5ea90792719d3442ad0b70e20",
          "url": "https://github.com/galacticusorg/galacticus/commit/150e6e81fa667e06ab5d1db655c8a6bfea150199"
        },
        "date": 1700268260687,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 143.042,
            "unit": "seconds",
            "range": 0.253691939171379
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "150e6e81fa667e06ab5d1db655c8a6bfea150199",
          "message": "Merge pull request #505 from galacticusorg/smbhPropertyOutput\n\nRefactor SMBH property output such that these are output to the `nodeData` group",
          "timestamp": "2023-11-17T22:22:56Z",
          "tree_id": "c1538aa88d1a2df5ea90792719d3442ad0b70e20",
          "url": "https://github.com/galacticusorg/galacticus/commit/150e6e81fa667e06ab5d1db655c8a6bfea150199"
        },
        "date": 1700268266612,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.15502294153761,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.92231083845279,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.593299212582026,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3430161971355,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 189.348730201804,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.0537302264587,
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
          "id": "f411bee51de1c75457b28077a09e74035ae015e9",
          "message": "feat: Remove debugging code",
          "timestamp": "2023-11-27T11:07:33-08:00",
          "tree_id": "e9c36b8e1721df0527ea44537aa2cf20ad755b62",
          "url": "https://github.com/galacticusorg/galacticus/commit/f411bee51de1c75457b28077a09e74035ae015e9"
        },
        "date": 1701124467904,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 170.338,
            "unit": "seconds",
            "range": 0.150477905353973
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
          "id": "f411bee51de1c75457b28077a09e74035ae015e9",
          "message": "feat: Remove debugging code",
          "timestamp": "2023-11-27T11:07:33-08:00",
          "tree_id": "e9c36b8e1721df0527ea44537aa2cf20ad755b62",
          "url": "https://github.com/galacticusorg/galacticus/commit/f411bee51de1c75457b28077a09e74035ae015e9"
        },
        "date": 1701124474290,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.47407780665755,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.37163159732478,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.04491724118473,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2331416934929,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 123.478000963063,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.8462548279889,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "60145e5c93a38a9103221c0ec80bf246656c0410",
          "message": "Merge pull request #507 from galacticusorg/openMPCriticals\n\nRemove some OpenMP critical sections",
          "timestamp": "2023-11-28T06:39:05Z",
          "tree_id": "5e51b48d4c733f98b344864c6d99eccd1adbcef1",
          "url": "https://github.com/galacticusorg/galacticus/commit/60145e5c93a38a9103221c0ec80bf246656c0410"
        },
        "date": 1701162181201,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 146.428,
            "unit": "seconds",
            "range": 0.261819785348983
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "60145e5c93a38a9103221c0ec80bf246656c0410",
          "message": "Merge pull request #507 from galacticusorg/openMPCriticals\n\nRemove some OpenMP critical sections",
          "timestamp": "2023-11-28T06:39:05Z",
          "tree_id": "5e51b48d4c733f98b344864c6d99eccd1adbcef1",
          "url": "https://github.com/galacticusorg/galacticus/commit/60145e5c93a38a9103221c0ec80bf246656c0410"
        },
        "date": 1701162187225,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.36067321870433,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.30446937878868,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.20066350382758,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1787597467278,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 121.354079112756,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.63914028137636,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "4ad1ed01f7180b8579aef9a691cb1d2e794a3463",
          "message": "Merge pull request #508 from galacticusorg/openMPCriticals\n\nRemove obsoleted OpenMP `critical` sections",
          "timestamp": "2023-11-28T21:47:23Z",
          "tree_id": "49482750bb25967d8e9369a7e1f944397a4521c2",
          "url": "https://github.com/galacticusorg/galacticus/commit/4ad1ed01f7180b8579aef9a691cb1d2e794a3463"
        },
        "date": 1701216781634,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 132.089,
            "unit": "seconds",
            "range": 0.387253018065046
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "4ad1ed01f7180b8579aef9a691cb1d2e794a3463",
          "message": "Merge pull request #508 from galacticusorg/openMPCriticals\n\nRemove obsoleted OpenMP `critical` sections",
          "timestamp": "2023-11-28T21:47:23Z",
          "tree_id": "49482750bb25967d8e9369a7e1f944397a4521c2",
          "url": "https://github.com/galacticusorg/galacticus/commit/4ad1ed01f7180b8579aef9a691cb1d2e794a3463"
        },
        "date": 1701216788094,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.52046857235995,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.70791184520522,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.01747541502993,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.247789925333,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 114.094186816649,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.49237622660232,
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
          "id": "3c8c9887321041b96ff0eb07f87090ec375bc042",
          "message": "fix: Revert incorrectly committed code",
          "timestamp": "2023-12-04T13:19:10-08:00",
          "tree_id": "8f1b01cd479346b1dea9a0db53dfd4747b137a17",
          "url": "https://github.com/galacticusorg/galacticus/commit/3c8c9887321041b96ff0eb07f87090ec375bc042"
        },
        "date": 1701733339128,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 195.966,
            "unit": "seconds",
            "range": 0.108519122733185
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
          "id": "3c8c9887321041b96ff0eb07f87090ec375bc042",
          "message": "fix: Revert incorrectly committed code",
          "timestamp": "2023-12-04T13:19:10-08:00",
          "tree_id": "8f1b01cd479346b1dea9a0db53dfd4747b137a17",
          "url": "https://github.com/galacticusorg/galacticus/commit/3c8c9887321041b96ff0eb07f87090ec375bc042"
        },
        "date": 1701733346813,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.50721527650958,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.21314268450191,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.953381440203268,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3133476958804,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 121.58464806153,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 7.65174618830223,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "1b4ef685df687dd6710fb5fdbdca360fcb841464",
          "message": "Merge pull request #510 from galacticusorg/agesFix\n\nCorrect indexing of meta-properties when computing stellar mass-weighted ages of disk and spheroid.",
          "timestamp": "2023-12-07T02:45:19Z",
          "tree_id": "e517388542a5d8e2423d2734ecdd285dbdc673b0",
          "url": "https://github.com/galacticusorg/galacticus/commit/1b4ef685df687dd6710fb5fdbdca360fcb841464"
        },
        "date": 1701925340287,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 175.287,
            "unit": "seconds",
            "range": 0.431476650585494
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "1b4ef685df687dd6710fb5fdbdca360fcb841464",
          "message": "Merge pull request #510 from galacticusorg/agesFix\n\nCorrect indexing of meta-properties when computing stellar mass-weighted ages of disk and spheroid.",
          "timestamp": "2023-12-07T02:45:19Z",
          "tree_id": "e517388542a5d8e2423d2734ecdd285dbdc673b0",
          "url": "https://github.com/galacticusorg/galacticus/commit/1b4ef685df687dd6710fb5fdbdca360fcb841464"
        },
        "date": 1701925348065,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.69133084011664,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.80172344441855,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.971517247365487,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0694849258449,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 123.663130274836,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.49489246491281,
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
          "id": "303504a414be4a532d36d27139be8ad3eadc6ded",
          "message": "fix: Merge branch 'agesFix'",
          "timestamp": "2023-12-12T08:12:47-08:00",
          "tree_id": "a1e796c7988ac1a0613cbad11bd9db6658e95d71",
          "url": "https://github.com/galacticusorg/galacticus/commit/303504a414be4a532d36d27139be8ad3eadc6ded"
        },
        "date": 1702406249893,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 149.374,
            "unit": "seconds",
            "range": 0.16238965484094
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
          "id": "303504a414be4a532d36d27139be8ad3eadc6ded",
          "message": "fix: Merge branch 'agesFix'",
          "timestamp": "2023-12-12T08:12:47-08:00",
          "tree_id": "a1e796c7988ac1a0613cbad11bd9db6658e95d71",
          "url": "https://github.com/galacticusorg/galacticus/commit/303504a414be4a532d36d27139be8ad3eadc6ded"
        },
        "date": 1702406257690,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.89258976759755,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.3952786701625,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.13560727215337,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9368755777895,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 113.965733473357,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.55556967366192,
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
          "id": "075e1183f4ca07ec9125497e850e9e918b97629c",
          "message": "fix: Update repo for codespaces image",
          "timestamp": "2023-12-13T10:30:27-08:00",
          "tree_id": "05128dc37177b9cb445c74cfe7ad63027c0941b3",
          "url": "https://github.com/galacticusorg/galacticus/commit/075e1183f4ca07ec9125497e850e9e918b97629c"
        },
        "date": 1702500853010,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 201.295,
            "unit": "seconds",
            "range": 0.182286861840554
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
          "id": "075e1183f4ca07ec9125497e850e9e918b97629c",
          "message": "fix: Update repo for codespaces image",
          "timestamp": "2023-12-13T10:30:27-08:00",
          "tree_id": "05128dc37177b9cb445c74cfe7ad63027c0941b3",
          "url": "https://github.com/galacticusorg/galacticus/commit/075e1183f4ca07ec9125497e850e9e918b97629c"
        },
        "date": 1702500860922,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.75082163360567,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.65265078716472,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.08959552614465,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9942329300209,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 115.595764854576,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.50779548803281,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "59f25e73014154ae31a6a2d797af229c7dbafd95",
          "message": "Merge pull request #511 from galacticusorg/fixOutputRank2ExtendSegFault\n\nFix a segfault associated with output of `rank?VarLen` datasets",
          "timestamp": "2023-12-15T02:07:31Z",
          "tree_id": "1807e9134e6c99eac3c8c7018039dc50c941a2df",
          "url": "https://github.com/galacticusorg/galacticus/commit/59f25e73014154ae31a6a2d797af229c7dbafd95"
        },
        "date": 1702617836717,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 161.878,
            "unit": "seconds",
            "range": 0.299619091515065
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "59f25e73014154ae31a6a2d797af229c7dbafd95",
          "message": "Merge pull request #511 from galacticusorg/fixOutputRank2ExtendSegFault\n\nFix a segfault associated with output of `rank?VarLen` datasets",
          "timestamp": "2023-12-15T02:07:31Z",
          "tree_id": "1807e9134e6c99eac3c8c7018039dc50c941a2df",
          "url": "https://github.com/galacticusorg/galacticus/commit/59f25e73014154ae31a6a2d797af229c7dbafd95"
        },
        "date": 1702617842770,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.69304732027374,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.61835421577017,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.14678653298399,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9058656291976,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 86.84331487735,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.50142452819427,
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
          "id": "238952f386539838344b47298704274566207a0c",
          "message": "Merge pull request #513 from galacticusorg/dependabot/github_actions/actions/upload-artifact-4\n\nbuild(deps): bump actions/upload-artifact from 3 to 4",
          "timestamp": "2023-12-18T03:57:38Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/238952f386539838344b47298704274566207a0c"
        },
        "date": 1702881671917,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 159.105,
            "unit": "seconds",
            "range": 0.125005999862403
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
          "id": "238952f386539838344b47298704274566207a0c",
          "message": "Merge pull request #513 from galacticusorg/dependabot/github_actions/actions/upload-artifact-4\n\nbuild(deps): bump actions/upload-artifact from 3 to 4",
          "timestamp": "2023-12-18T03:57:38Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/238952f386539838344b47298704274566207a0c"
        },
        "date": 1702881678627,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.15830176744519,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.10894600847223,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.649096979135215,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.4216941877242,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 188.517772924876,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 2.55626076502453,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "f36b106bdbda4fe2eb2649b12c2ae6bcad732189",
          "message": "Merge pull request #515 from galacticusorg/directivesErrorMessage\n\nAdd more useful error reporting when failing to parse source code directives",
          "timestamp": "2023-12-19T16:08:42Z",
          "tree_id": "1c92c7ad381ef098f2383587e35deb33a1b70ea1",
          "url": "https://github.com/galacticusorg/galacticus/commit/f36b106bdbda4fe2eb2649b12c2ae6bcad732189"
        },
        "date": 1703010704245,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 149.057,
            "unit": "seconds",
            "range": 0.349934422428436
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "f36b106bdbda4fe2eb2649b12c2ae6bcad732189",
          "message": "Merge pull request #515 from galacticusorg/directivesErrorMessage\n\nAdd more useful error reporting when failing to parse source code directives",
          "timestamp": "2023-12-19T16:08:42Z",
          "tree_id": "1c92c7ad381ef098f2383587e35deb33a1b70ea1",
          "url": "https://github.com/galacticusorg/galacticus/commit/f36b106bdbda4fe2eb2649b12c2ae6bcad732189"
        },
        "date": 1703010711931,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.03869615841866,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.80568382899929,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.66699871731044,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2159328108994,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 188.43241494481,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 2.53470107013292,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9a43b429c981765a16cde8b7e826d206d4d4c677",
          "message": "Merge pull request #516 from galacticusorg/adaptiveSFHLengths\n\nAdd missing `update` method to `starFormationHistoryAdsptive` class",
          "timestamp": "2023-12-21T05:48:51Z",
          "tree_id": "a292640f1ebe2735a04808208c19d1d76558e699",
          "url": "https://github.com/galacticusorg/galacticus/commit/9a43b429c981765a16cde8b7e826d206d4d4c677"
        },
        "date": 1703145973120,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 156.421,
            "unit": "seconds",
            "range": 0.137385952702897
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9a43b429c981765a16cde8b7e826d206d4d4c677",
          "message": "Merge pull request #516 from galacticusorg/adaptiveSFHLengths\n\nAdd missing `update` method to `starFormationHistoryAdsptive` class",
          "timestamp": "2023-12-21T05:48:51Z",
          "tree_id": "a292640f1ebe2735a04808208c19d1d76558e699",
          "url": "https://github.com/galacticusorg/galacticus/commit/9a43b429c981765a16cde8b7e826d206d4d4c677"
        },
        "date": 1703145979197,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.42213147422155,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.09345651150396,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.06049405980213,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0227743355522,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 120.34076655521,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 7.63492270806651,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "f167a854c58407da4412c3812d4994c0394bdb73",
          "message": "Merge pull request #521 from galacticusorg/fixAdaptiveSFHs\n\nFix adapative star formation history lengths and times",
          "timestamp": "2024-01-05T23:32:53Z",
          "tree_id": "35b78c2092837d0d11ad0cc78d912f2ed5319589",
          "url": "https://github.com/galacticusorg/galacticus/commit/f167a854c58407da4412c3812d4994c0394bdb73"
        },
        "date": 1704508903102,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 168.567,
            "unit": "seconds",
            "range": 0.200110219625845
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "f167a854c58407da4412c3812d4994c0394bdb73",
          "message": "Merge pull request #521 from galacticusorg/fixAdaptiveSFHs\n\nFix adapative star formation history lengths and times",
          "timestamp": "2024-01-05T23:32:53Z",
          "tree_id": "35b78c2092837d0d11ad0cc78d912f2ed5319589",
          "url": "https://github.com/galacticusorg/galacticus/commit/f167a854c58407da4412c3812d4994c0394bdb73"
        },
        "date": 1704508909320,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.36677198217072,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.7418641768436,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.14978263673784,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2183475521532,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 122.551893993565,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.7430169956265,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "33690072eff234ab616b7d59e781d60d84799539",
          "message": "Merge pull request #527 from galacticusorg/featMagnitudes\n\nAdd `nodePropertyExtractor`s to output absolute and apparent magnitudes",
          "timestamp": "2024-01-11T21:37:54Z",
          "tree_id": "5b3c95ad0d45637f63ad1becb71b4464b50c4532",
          "url": "https://github.com/galacticusorg/galacticus/commit/33690072eff234ab616b7d59e781d60d84799539"
        },
        "date": 1705037111090,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 164.91,
            "unit": "seconds",
            "range": 0.159974998046446
          }
        ]
      },
      {
        "commit": {
          "author": {
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
          "id": "33690072eff234ab616b7d59e781d60d84799539",
          "message": "Merge pull request #527 from galacticusorg/featMagnitudes\n\nAdd `nodePropertyExtractor`s to output absolute and apparent magnitudes",
          "timestamp": "2024-01-11T21:37:54Z",
          "tree_id": "5b3c95ad0d45637f63ad1becb71b4464b50c4532",
          "url": "https://github.com/galacticusorg/galacticus/commit/33690072eff234ab616b7d59e781d60d84799539"
        },
        "date": 1705037118901,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 2.93522857317178,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.29214981516082,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.877265577391712,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.6108078765042,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 375.380644563075,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 4.0703176705196,
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
          "id": "a34f4337ce846a73c1c3923bbda018546425c129",
          "message": "fix: Ensure variables are deallocated before potential re-allocation",
          "timestamp": "2024-01-17T14:52:00-08:00",
          "tree_id": "5af2d75c2b5c804f4ee9b071e88a3640ed350e30",
          "url": "https://github.com/galacticusorg/galacticus/commit/a34f4337ce846a73c1c3923bbda018546425c129"
        },
        "date": 1705552402730,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 166.493,
            "unit": "seconds",
            "range": 0.300599567530987
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
          "id": "a34f4337ce846a73c1c3923bbda018546425c129",
          "message": "fix: Ensure variables are deallocated before potential re-allocation",
          "timestamp": "2024-01-17T14:52:00-08:00",
          "tree_id": "5af2d75c2b5c804f4ee9b071e88a3640ed350e30",
          "url": "https://github.com/galacticusorg/galacticus/commit/a34f4337ce846a73c1c3923bbda018546425c129"
        },
        "date": 1705552408765,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.13227358249588,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.05511220140172,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.07387765728136,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3230197314207,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 122.779370258702,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.6484955931643,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "38eb56e63706ddadc6d5d8041c73d21137c8cf6a",
          "message": "Merge pull request #532 from galacticusorg/fixZeroStellarHistories\n\nEnsure stellar histories are zeroed when zeroing stellar mass due to ODE solver failure",
          "timestamp": "2024-01-18T15:20:03Z",
          "tree_id": "052cbacc96846aa8dd513d1abd956fec1d04b0f9",
          "url": "https://github.com/galacticusorg/galacticus/commit/38eb56e63706ddadc6d5d8041c73d21137c8cf6a"
        },
        "date": 1705607256175,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 141.753,
            "unit": "seconds",
            "range": 0.163878308509935
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "38eb56e63706ddadc6d5d8041c73d21137c8cf6a",
          "message": "Merge pull request #532 from galacticusorg/fixZeroStellarHistories\n\nEnsure stellar histories are zeroed when zeroing stellar mass due to ODE solver failure",
          "timestamp": "2024-01-18T15:20:03Z",
          "tree_id": "052cbacc96846aa8dd513d1abd956fec1d04b0f9",
          "url": "https://github.com/galacticusorg/galacticus/commit/38eb56e63706ddadc6d5d8041c73d21137c8cf6a"
        },
        "date": 1705607262354,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.23581063527455,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 5.12778440851184,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.665784165451604,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3949227267269,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 190.530964778531,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 2.52518522818629,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "738907272e8bedcc3d1c8a688abdd98d93bf258d",
          "message": "Merge pull request #529 from galacticusorg/featReuseableWorkflows\n\nMake use of reuseable workflows to simplify GitHub Actions workflow files",
          "timestamp": "2024-01-20T20:06:59Z",
          "tree_id": "146720a7231dc2c5d78c9c8b9783bbfbc4242b64",
          "url": "https://github.com/galacticusorg/galacticus/commit/738907272e8bedcc3d1c8a688abdd98d93bf258d"
        },
        "date": 1705792183271,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 164.417,
            "unit": "seconds",
            "range": 1.35993312335617
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "738907272e8bedcc3d1c8a688abdd98d93bf258d",
          "message": "Merge pull request #529 from galacticusorg/featReuseableWorkflows\n\nMake use of reuseable workflows to simplify GitHub Actions workflow files",
          "timestamp": "2024-01-20T20:06:59Z",
          "tree_id": "146720a7231dc2c5d78c9c8b9783bbfbc4242b64",
          "url": "https://github.com/galacticusorg/galacticus/commit/738907272e8bedcc3d1c8a688abdd98d93bf258d"
        },
        "date": 1705792190720,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.20890533289581,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.37256417745761,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.03038660268652,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0849085638188,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 122.411696654471,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.7563384932922,
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
          "id": "738907272e8bedcc3d1c8a688abdd98d93bf258d",
          "message": "Merge pull request #529 from galacticusorg/featReuseableWorkflows\n\nMake use of reuseable workflows to simplify GitHub Actions workflow files",
          "timestamp": "2024-01-20T20:06:59Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/738907272e8bedcc3d1c8a688abdd98d93bf258d"
        },
        "date": 1705871798969,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 168.784,
            "unit": "seconds",
            "range": 0.351557107734947
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
          "id": "738907272e8bedcc3d1c8a688abdd98d93bf258d",
          "message": "Merge pull request #529 from galacticusorg/featReuseableWorkflows\n\nMake use of reuseable workflows to simplify GitHub Actions workflow files",
          "timestamp": "2024-01-20T20:06:59Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/738907272e8bedcc3d1c8a688abdd98d93bf258d"
        },
        "date": 1705871805691,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.28742480387074,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.25572450949796,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.08946166008693,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2690341514682,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 122.076318656751,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.6431024061727,
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
          "id": "a7a30ec3c7ae4ec5df033cab4bcaeb0dee7899fb",
          "message": "fix: Merge branch 'featReuseableWorkflows'",
          "timestamp": "2024-01-21T22:47:59Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/a7a30ec3c7ae4ec5df033cab4bcaeb0dee7899fb"
        },
        "date": 1705891064675,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 168.653,
            "unit": "seconds",
            "range": 0.184309793555423
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
          "id": "a7a30ec3c7ae4ec5df033cab4bcaeb0dee7899fb",
          "message": "fix: Merge branch 'featReuseableWorkflows'",
          "timestamp": "2024-01-21T22:47:59Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/a7a30ec3c7ae4ec5df033cab4bcaeb0dee7899fb"
        },
        "date": 1705891071277,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.41185719503308,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.52306990951796,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.94243865656416,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9475534493175,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 85.4576585727402,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.66246714595877,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "426d47542279966670f16b52132911ffe1777496",
          "message": "Merge pull request #535 from galacticusorg/featFormationHalos\n\nAdd a `stellarFeedbackOutflows` class that uses the formation halo properties",
          "timestamp": "2024-01-26T06:09:08Z",
          "tree_id": "f93b656986590f5805ea143b9bfe76e911ec9ff4",
          "url": "https://github.com/galacticusorg/galacticus/commit/426d47542279966670f16b52132911ffe1777496"
        },
        "date": 1706260196258,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 169.002,
            "unit": "seconds",
            "range": 0.249250075225373
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "426d47542279966670f16b52132911ffe1777496",
          "message": "Merge pull request #535 from galacticusorg/featFormationHalos\n\nAdd a `stellarFeedbackOutflows` class that uses the formation halo properties",
          "timestamp": "2024-01-26T06:09:08Z",
          "tree_id": "f93b656986590f5805ea143b9bfe76e911ec9ff4",
          "url": "https://github.com/galacticusorg/galacticus/commit/426d47542279966670f16b52132911ffe1777496"
        },
        "date": 1706260202700,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.46467025534774,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.35378778486725,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.02650498039008,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0547159709972,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 121.877864036476,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.60133249781566,
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
          "id": "d277e62793a5e0160ff21b58bc3f4ae51e878bd4",
          "message": "fix: Merge branch 'master' of github.com:galacticusorg/galacticus",
          "timestamp": "2024-01-26T06:58:04-08:00",
          "tree_id": "790ade19acce7fdaf5f04fcfb390aeec709aac2c",
          "url": "https://github.com/galacticusorg/galacticus/commit/d277e62793a5e0160ff21b58bc3f4ae51e878bd4"
        },
        "date": 1706291824261,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 155.675,
            "unit": "seconds",
            "range": 0.151256404822848
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
          "id": "d277e62793a5e0160ff21b58bc3f4ae51e878bd4",
          "message": "fix: Merge branch 'master' of github.com:galacticusorg/galacticus",
          "timestamp": "2024-01-26T06:58:04-08:00",
          "tree_id": "790ade19acce7fdaf5f04fcfb390aeec709aac2c",
          "url": "https://github.com/galacticusorg/galacticus/commit/d277e62793a5e0160ff21b58bc3f4ae51e878bd4"
        },
        "date": 1706291831776,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.26620252383314,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.17167772936361,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.04092742075961,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 23.9937426253605,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 85.7336818051231,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.6486979306609,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "74c5285047ac77fa4398d9fe601da30c6427f868",
          "message": "Merge pull request #544 from galacticusorg/featPowerSpectrumNormalization\n\nAdd functionality to normalize the power spectrum via the $A_\\mathrm{s}$ parameter",
          "timestamp": "2024-01-27T02:45:59Z",
          "tree_id": "70c701f7ffed6f95139fc6eda794733aff851317",
          "url": "https://github.com/galacticusorg/galacticus/commit/74c5285047ac77fa4398d9fe601da30c6427f868"
        },
        "date": 1706334500607,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 180.959,
            "unit": "seconds",
            "range": 0.172113044249081
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "74c5285047ac77fa4398d9fe601da30c6427f868",
          "message": "Merge pull request #544 from galacticusorg/featPowerSpectrumNormalization\n\nAdd functionality to normalize the power spectrum via the $A_\\mathrm{s}$ parameter",
          "timestamp": "2024-01-27T02:45:59Z",
          "tree_id": "70c701f7ffed6f95139fc6eda794733aff851317",
          "url": "https://github.com/galacticusorg/galacticus/commit/74c5285047ac77fa4398d9fe601da30c6427f868"
        },
        "date": 1706334507657,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.35970234732318,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.66351561558961,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.18156218572975,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1798793944412,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 85.9539318000147,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.61110260930108,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c8bcf117d43a0d03576a6fa1ca4bb43be11abd09",
          "message": "Merge pull request #546 from galacticusorg/fixSubhaloPromotion\n\nFix subhalo promotion",
          "timestamp": "2024-01-30T05:36:55Z",
          "tree_id": "fbcda48fe6b54707b05014d67d062b3566d8af2a",
          "url": "https://github.com/galacticusorg/galacticus/commit/c8bcf117d43a0d03576a6fa1ca4bb43be11abd09"
        },
        "date": 1706604150971,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 201.677,
            "unit": "seconds",
            "range": 0.341783703533199
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c8bcf117d43a0d03576a6fa1ca4bb43be11abd09",
          "message": "Merge pull request #546 from galacticusorg/fixSubhaloPromotion\n\nFix subhalo promotion",
          "timestamp": "2024-01-30T05:36:55Z",
          "tree_id": "fbcda48fe6b54707b05014d67d062b3566d8af2a",
          "url": "https://github.com/galacticusorg/galacticus/commit/c8bcf117d43a0d03576a6fa1ca4bb43be11abd09"
        },
        "date": 1706604157404,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.26193195193211,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.16271737994937,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.30264711836685,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2512607570903,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 136.171519166732,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.59422789325731,
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
          "id": "970b8cf5439ec2ea0b846c9c8faa9735e170e6e9",
          "message": "fix: Ensure a host tree is set for unresolved halos in the `darkMatterProfileScaleRadiusJohnson2021` class\n\nThis is needed to allow access to random number generators attached to trees (for example, by a concentration class that wants to include some random scatter).",
          "timestamp": "2024-02-01T14:40:47-08:00",
          "tree_id": "eb4984ebdba10a91844d9e70bf96a10cc9d8ef6f",
          "url": "https://github.com/galacticusorg/galacticus/commit/970b8cf5439ec2ea0b846c9c8faa9735e170e6e9"
        },
        "date": 1706838425686,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 200.129,
            "unit": "seconds",
            "range": 0.492754401300353
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
          "id": "970b8cf5439ec2ea0b846c9c8faa9735e170e6e9",
          "message": "fix: Ensure a host tree is set for unresolved halos in the `darkMatterProfileScaleRadiusJohnson2021` class\n\nThis is needed to allow access to random number generators attached to trees (for example, by a concentration class that wants to include some random scatter).",
          "timestamp": "2024-02-01T14:40:47-08:00",
          "tree_id": "eb4984ebdba10a91844d9e70bf96a10cc9d8ef6f",
          "url": "https://github.com/galacticusorg/galacticus/commit/970b8cf5439ec2ea0b846c9c8faa9735e170e6e9"
        },
        "date": 1706838431663,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.37821951497771,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.73029549517395,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.18035674057112,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1405483143522,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 121.994173265711,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.61735255527312,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "cf5713b6826e52c776078beb867990d97e1d28b7",
          "message": "Merge pull request #547 from galacticusorg/fixStaticBuilds\n\nFix static builds of FSPS and CLASS",
          "timestamp": "2024-02-02T01:51:47Z",
          "tree_id": "75079517fe7039af285f228acf638968855e1059",
          "url": "https://github.com/galacticusorg/galacticus/commit/cf5713b6826e52c776078beb867990d97e1d28b7"
        },
        "date": 1706849546862,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 168.954,
            "unit": "seconds",
            "range": 0.242215606433206
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "cf5713b6826e52c776078beb867990d97e1d28b7",
          "message": "Merge pull request #547 from galacticusorg/fixStaticBuilds\n\nFix static builds of FSPS and CLASS",
          "timestamp": "2024-02-02T01:51:47Z",
          "tree_id": "75079517fe7039af285f228acf638968855e1059",
          "url": "https://github.com/galacticusorg/galacticus/commit/cf5713b6826e52c776078beb867990d97e1d28b7"
        },
        "date": 1706849554742,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.20351259769759,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.53064755776253,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.17299481904183,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1579471695706,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 88.1225835765995,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 14.000613786469,
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
          "distinct": false,
          "id": "b8d767f0c5dbece37ae2b52ace6d337d81f359ab",
          "message": "feat: Provide a more useful error message for missing power spectra files",
          "timestamp": "2024-02-06T08:18:17-08:00",
          "tree_id": "3ae7e1270216c6e85d40e63e4b51e2fe61e8f594",
          "url": "https://github.com/galacticusorg/galacticus/commit/b8d767f0c5dbece37ae2b52ace6d337d81f359ab"
        },
        "date": 1707261303539,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 154.724,
            "unit": "seconds",
            "range": 0.209962853858555
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
          "distinct": false,
          "id": "b8d767f0c5dbece37ae2b52ace6d337d81f359ab",
          "message": "feat: Provide a more useful error message for missing power spectra files",
          "timestamp": "2024-02-06T08:18:17-08:00",
          "tree_id": "3ae7e1270216c6e85d40e63e4b51e2fe61e8f594",
          "url": "https://github.com/galacticusorg/galacticus/commit/b8d767f0c5dbece37ae2b52ace6d337d81f359ab"
        },
        "date": 1707261310384,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.40045633572226,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.60680959270197,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.14221365309542,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1747898320541,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 86.8241400845919,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.64327462280099,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6590a605538a5f97ff3cd888539564acc1cc1d99",
          "message": "Merge pull request #552 from galacticusorg/fixPowerSpectrumErrorReporting\n\nFix power spectrum error reporting",
          "timestamp": "2024-02-07T02:14:17Z",
          "tree_id": "8c7393031e79e70b596ff31d5ade12adfc81ad68",
          "url": "https://github.com/galacticusorg/galacticus/commit/6590a605538a5f97ff3cd888539564acc1cc1d99"
        },
        "date": 1707282917358,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 161.907,
            "unit": "seconds",
            "range": 0.143471599975541
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6590a605538a5f97ff3cd888539564acc1cc1d99",
          "message": "Merge pull request #552 from galacticusorg/fixPowerSpectrumErrorReporting\n\nFix power spectrum error reporting",
          "timestamp": "2024-02-07T02:14:17Z",
          "tree_id": "8c7393031e79e70b596ff31d5ade12adfc81ad68",
          "url": "https://github.com/galacticusorg/galacticus/commit/6590a605538a5f97ff3cd888539564acc1cc1d99"
        },
        "date": 1707282923488,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.22960812992491,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.37642692126392,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.13694501516069,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.461918730439,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 87.2042952858633,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.9852144208426,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d0c330dce5fc9b6d455e455762c7e42167e50ac0",
          "message": "Merge pull request #553 from galacticusorg/fixHDF5VarLenStringRead\n\nAllow reading of variable length string types from HDF5 attributes",
          "timestamp": "2024-02-07T23:44:39Z",
          "tree_id": "e2ad3ae35eb6b1173f6eb7ebd44473d103daff08",
          "url": "https://github.com/galacticusorg/galacticus/commit/d0c330dce5fc9b6d455e455762c7e42167e50ac0"
        },
        "date": 1707360461370,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 156.741,
            "unit": "seconds",
            "range": 0.099754197904666
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d0c330dce5fc9b6d455e455762c7e42167e50ac0",
          "message": "Merge pull request #553 from galacticusorg/fixHDF5VarLenStringRead\n\nAllow reading of variable length string types from HDF5 attributes",
          "timestamp": "2024-02-07T23:44:39Z",
          "tree_id": "e2ad3ae35eb6b1173f6eb7ebd44473d103daff08",
          "url": "https://github.com/galacticusorg/galacticus/commit/d0c330dce5fc9b6d455e455762c7e42167e50ac0"
        },
        "date": 1707360468054,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.20692948724079,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.20306719736092,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.13123676748436,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0120062795959,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 122.194220953562,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.57814698511502,
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
          "id": "ed6bffd6daef4250ac13458f1dfc7b0f10d4de32",
          "message": "feat: Improve documentation",
          "timestamp": "2024-02-08T09:08:07-08:00",
          "tree_id": "1f949b67e972f8182efc4a51cc818fc2707c352d",
          "url": "https://github.com/galacticusorg/galacticus/commit/ed6bffd6daef4250ac13458f1dfc7b0f10d4de32"
        },
        "date": 1707423015105,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 159.085,
            "unit": "seconds",
            "range": 0.18163287147394
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
          "id": "ed6bffd6daef4250ac13458f1dfc7b0f10d4de32",
          "message": "feat: Improve documentation",
          "timestamp": "2024-02-08T09:08:07-08:00",
          "tree_id": "1f949b67e972f8182efc4a51cc818fc2707c352d",
          "url": "https://github.com/galacticusorg/galacticus/commit/ed6bffd6daef4250ac13458f1dfc7b0f10d4de32"
        },
        "date": 1707423020874,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.04320293481181,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 3.9552490996047,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.0598854107372,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.3402707952584,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 122.012368225686,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.6124974390149,
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
          "id": "ed1160a3906155fe885404761ec066763ace0832",
          "message": "fix: Fix bug in tablulation of times for the `powerSpectrumPrimordialTransferredFile` class\n\nA missing index caused all times in the array to be identical, breaking interpolation.\n\nAlso makes warnings about mismatched cosmological parameters follow new conventions (i.e. magenta \"WARNING:\" included at the start of the message).",
          "timestamp": "2024-02-09T07:52:01-08:00",
          "tree_id": "464b9e14277cbe56ad0b333d2cafb7fa7390de98",
          "url": "https://github.com/galacticusorg/galacticus/commit/ed1160a3906155fe885404761ec066763ace0832"
        },
        "date": 1707504871830,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 171.133,
            "unit": "seconds",
            "range": 0.204235403394854
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
          "id": "ed1160a3906155fe885404761ec066763ace0832",
          "message": "fix: Fix bug in tablulation of times for the `powerSpectrumPrimordialTransferredFile` class\n\nA missing index caused all times in the array to be identical, breaking interpolation.\n\nAlso makes warnings about mismatched cosmological parameters follow new conventions (i.e. magenta \"WARNING:\" included at the start of the message).",
          "timestamp": "2024-02-09T07:52:01-08:00",
          "tree_id": "464b9e14277cbe56ad0b333d2cafb7fa7390de98",
          "url": "https://github.com/galacticusorg/galacticus/commit/ed1160a3906155fe885404761ec066763ace0832"
        },
        "date": 1707504879596,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.26951075316755,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.20111441249431,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.09285926390837,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2717199523266,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 121.960607996758,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.6296644557524,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ec4a80d5bb4fa2fe1060bb01f2ac021e1eba8ce5",
          "message": "Merge pull request #555 from galacticusorg/validatePonosV\n\nAdd validation against the [PonosV](https://ui.adsabs.harvard.edu/abs/2016ApJ...824..144F) simulation measurement of the surface density of subhalos",
          "timestamp": "2024-02-14T21:50:31Z",
          "tree_id": "4d3cdf3f89ce955df3cfa4f4a189a95d9c2baa92",
          "url": "https://github.com/galacticusorg/galacticus/commit/ec4a80d5bb4fa2fe1060bb01f2ac021e1eba8ce5"
        },
        "date": 1707960183566,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 163.253,
            "unit": "seconds",
            "range": 0.137448535825044
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ec4a80d5bb4fa2fe1060bb01f2ac021e1eba8ce5",
          "message": "Merge pull request #555 from galacticusorg/validatePonosV\n\nAdd validation against the [PonosV](https://ui.adsabs.harvard.edu/abs/2016ApJ...824..144F) simulation measurement of the surface density of subhalos",
          "timestamp": "2024-02-14T21:50:31Z",
          "tree_id": "4d3cdf3f89ce955df3cfa4f4a189a95d9c2baa92",
          "url": "https://github.com/galacticusorg/galacticus/commit/ec4a80d5bb4fa2fe1060bb01f2ac021e1eba8ce5"
        },
        "date": 1707960189700,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.2555799784572,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.17264021557409,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.24835764547476,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2899336407244,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 122.321596023826,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.58870206342175,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5c3ea4a16d359eba4474d28f20d054c8be983b0e",
          "message": "Merge pull request #558 from galacticusorg/fixVarStringToVarString\n\nAllow read of an HDF5 variable length string to a `varying_string` object",
          "timestamp": "2024-02-19T22:29:35Z",
          "tree_id": "d2c78a681995907b0f6f65163637e969a64aade3",
          "url": "https://github.com/galacticusorg/galacticus/commit/5c3ea4a16d359eba4474d28f20d054c8be983b0e"
        },
        "date": 1708392862178,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 157.105,
            "unit": "seconds",
            "range": 0.141125830384223
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5c3ea4a16d359eba4474d28f20d054c8be983b0e",
          "message": "Merge pull request #558 from galacticusorg/fixVarStringToVarString\n\nAllow read of an HDF5 variable length string to a `varying_string` object",
          "timestamp": "2024-02-19T22:29:35Z",
          "tree_id": "d2c78a681995907b0f6f65163637e969a64aade3",
          "url": "https://github.com/galacticusorg/galacticus/commit/5c3ea4a16d359eba4474d28f20d054c8be983b0e"
        },
        "date": 1708392868006,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.08238439935854,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.13173302742336,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.0280534819961,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.263055274722,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 120.883085056726,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.5858199650133,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "f81750a3a21fe16e9847484b47d33825b7bc153f",
          "message": "Merge pull request #559 from galacticusorg/fixDeprecatePowerLaw\n\nDeprecate the (non-implemented) `powerLaw` extrapolation method for tables",
          "timestamp": "2024-02-21T03:38:39Z",
          "tree_id": "73a6a36c1e5bc1950812ffa01c625abeaa5a2d71",
          "url": "https://github.com/galacticusorg/galacticus/commit/f81750a3a21fe16e9847484b47d33825b7bc153f"
        },
        "date": 1708498531759,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 170.29,
            "unit": "seconds",
            "range": 0.187781788253771
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "f81750a3a21fe16e9847484b47d33825b7bc153f",
          "message": "Merge pull request #559 from galacticusorg/fixDeprecatePowerLaw\n\nDeprecate the (non-implemented) `powerLaw` extrapolation method for tables",
          "timestamp": "2024-02-21T03:38:39Z",
          "tree_id": "73a6a36c1e5bc1950812ffa01c625abeaa5a2d71",
          "url": "https://github.com/galacticusorg/galacticus/commit/f81750a3a21fe16e9847484b47d33825b7bc153f"
        },
        "date": 1708498538709,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.38954720826258,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.11136780236522,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.0006769282573,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1443349777288,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 87.1409483178553,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.60569780010348,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2347f1be84d73ea4f54355a8032e19069005df36",
          "message": "Merge pull request #566 from galacticusorg/fixVelocityOrbitalOutput\n\nMake orbital velocities relative to the eventual host halo when tracking pre-infall orbits",
          "timestamp": "2024-02-23T16:27:57Z",
          "tree_id": "e99409d613819615520400e14ac22da29f68ea4f",
          "url": "https://github.com/galacticusorg/galacticus/commit/2347f1be84d73ea4f54355a8032e19069005df36"
        },
        "date": 1708716762268,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 142.058,
            "unit": "seconds",
            "range": 0.152287885268986
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2347f1be84d73ea4f54355a8032e19069005df36",
          "message": "Merge pull request #566 from galacticusorg/fixVelocityOrbitalOutput\n\nMake orbital velocities relative to the eventual host halo when tracking pre-infall orbits",
          "timestamp": "2024-02-23T16:27:57Z",
          "tree_id": "e99409d613819615520400e14ac22da29f68ea4f",
          "url": "https://github.com/galacticusorg/galacticus/commit/2347f1be84d73ea4f54355a8032e19069005df36"
        },
        "date": 1708716769470,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.56085250951693,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.65938581251365,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.06899776645061,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0994077811737,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 121.307681128465,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.63620570693815,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "fe19e1a2fd3a716ef1ef35aed4e29eaf049dbf9d",
          "message": "Merge pull request #567 from galacticusorg/featOutputSelectorTimes\n\nAllow specifying output selection via time in addition to redshift",
          "timestamp": "2024-02-27T18:11:12Z",
          "tree_id": "6f70047c8290171cd9eb1c3d8b80811457587cd9",
          "url": "https://github.com/galacticusorg/galacticus/commit/fe19e1a2fd3a716ef1ef35aed4e29eaf049dbf9d"
        },
        "date": 1709078514133,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 150.987,
            "unit": "seconds",
            "range": 0.0646846194942114
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "fe19e1a2fd3a716ef1ef35aed4e29eaf049dbf9d",
          "message": "Merge pull request #567 from galacticusorg/featOutputSelectorTimes\n\nAllow specifying output selection via time in addition to redshift",
          "timestamp": "2024-02-27T18:11:12Z",
          "tree_id": "6f70047c8290171cd9eb1c3d8b80811457587cd9",
          "url": "https://github.com/galacticusorg/galacticus/commit/fe19e1a2fd3a716ef1ef35aed4e29eaf049dbf9d"
        },
        "date": 1709078520355,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.46448861772289,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.78393121510301,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.19170350679505,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2325982115047,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 86.2368241760886,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.9380560162071,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c851a21a3844ce74c1c35fb4f767761e802e6798",
          "message": "Merge pull request #568 from galacticusorg/featRunTimeFileDependencies\n\nAdd runtime file dependencies to `functionClass` `descriptor` methods",
          "timestamp": "2024-02-28T04:05:28Z",
          "tree_id": "403b915da36aa6ddf0c127bb18c7f568bd796bb7",
          "url": "https://github.com/galacticusorg/galacticus/commit/c851a21a3844ce74c1c35fb4f767761e802e6798"
        },
        "date": 1709104161494,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 187.181,
            "unit": "seconds",
            "range": 0.186426661184909
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c851a21a3844ce74c1c35fb4f767761e802e6798",
          "message": "Merge pull request #568 from galacticusorg/featRunTimeFileDependencies\n\nAdd runtime file dependencies to `functionClass` `descriptor` methods",
          "timestamp": "2024-02-28T04:05:28Z",
          "tree_id": "403b915da36aa6ddf0c127bb18c7f568bd796bb7",
          "url": "https://github.com/galacticusorg/galacticus/commit/c851a21a3844ce74c1c35fb4f767761e802e6798"
        },
        "date": 1709104167558,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.21357390954854,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.68559372591964,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.05964304157915,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1381155021015,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 86.839809127165,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.9854366925232,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5be0e4369d4aa4cd34b3ff88341fdb2ccb243dda",
          "message": "Merge pull request #569 from galacticusorg/featDetectDuplicatePropertyNames\n\nAdd detection of duplicated property names in output",
          "timestamp": "2024-02-29T04:39:01Z",
          "tree_id": "29cb5d856d39db670739d75acdb53c9c9ed0bfb5",
          "url": "https://github.com/galacticusorg/galacticus/commit/5be0e4369d4aa4cd34b3ff88341fdb2ccb243dda"
        },
        "date": 1709192655811,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 184.399,
            "unit": "seconds",
            "range": 0.365891923933079
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5be0e4369d4aa4cd34b3ff88341fdb2ccb243dda",
          "message": "Merge pull request #569 from galacticusorg/featDetectDuplicatePropertyNames\n\nAdd detection of duplicated property names in output",
          "timestamp": "2024-02-29T04:39:01Z",
          "tree_id": "29cb5d856d39db670739d75acdb53c9c9ed0bfb5",
          "url": "https://github.com/galacticusorg/galacticus/commit/5be0e4369d4aa4cd34b3ff88341fdb2ccb243dda"
        },
        "date": 1709192662436,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.19202534317258,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.33622987334921,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.08475433833191,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0661142130105,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 119.92808964816,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.7804872457243,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ba8109553e88d7c84ab0d67ab459aa8fd692cc27",
          "message": "Merge pull request #570 from galacticusorg/featFlexibleTypeIaDTDs\n\nRefactor the supernovae type Ia delay time distribution class",
          "timestamp": "2024-02-29T23:20:21Z",
          "tree_id": "3689bf360207b35fef58f17fbf48ede2170a9fa2",
          "url": "https://github.com/galacticusorg/galacticus/commit/ba8109553e88d7c84ab0d67ab459aa8fd692cc27"
        },
        "date": 1709261413760,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 156.324,
            "unit": "seconds",
            "range": 0.152270811384712
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ba8109553e88d7c84ab0d67ab459aa8fd692cc27",
          "message": "Merge pull request #570 from galacticusorg/featFlexibleTypeIaDTDs\n\nRefactor the supernovae type Ia delay time distribution class",
          "timestamp": "2024-02-29T23:20:21Z",
          "tree_id": "3689bf360207b35fef58f17fbf48ede2170a9fa2",
          "url": "https://github.com/galacticusorg/galacticus/commit/ba8109553e88d7c84ab0d67ab459aa8fd692cc27"
        },
        "date": 1709261419793,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.33736537407723,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.24195273876811,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.2238428991598,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0759399626732,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 88.0506295675713,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.59897416912105,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9a381319ea2196b3205966de35e5489e97f28118",
          "message": "Merge pull request #571 from galacticusorg/fixOptimizeBaryonicLSS\n\nOptimize baryonic spherical collapse solvers",
          "timestamp": "2024-03-01T20:22:50Z",
          "tree_id": "dc0476744365d235537e8ac69282d5dfceb7b4b9",
          "url": "https://github.com/galacticusorg/galacticus/commit/9a381319ea2196b3205966de35e5489e97f28118"
        },
        "date": 1709335739322,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 176.246,
            "unit": "seconds",
            "range": 0.206572989521184
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9a381319ea2196b3205966de35e5489e97f28118",
          "message": "Merge pull request #571 from galacticusorg/fixOptimizeBaryonicLSS\n\nOptimize baryonic spherical collapse solvers",
          "timestamp": "2024-03-01T20:22:50Z",
          "tree_id": "dc0476744365d235537e8ac69282d5dfceb7b4b9",
          "url": "https://github.com/galacticusorg/galacticus/commit/9a381319ea2196b3205966de35e5489e97f28118"
        },
        "date": 1709335745612,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.27172295785737,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.6341074608852,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.0764487771902,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0073040890109,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 121.074090792656,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.8017208728635,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d01251743649ae4335a060127395563cfd86d513",
          "message": "Merge pull request #573 from galacticusorg/fixOptimizeBaryonicLSS\n\nUse exclusive lock to avoid conflicts between processes",
          "timestamp": "2024-03-02T14:57:44Z",
          "tree_id": "85126c35f19e3f00a0347a88f8e1184e72a579d4",
          "url": "https://github.com/galacticusorg/galacticus/commit/d01251743649ae4335a060127395563cfd86d513"
        },
        "date": 1709402622034,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 152.54,
            "unit": "seconds",
            "range": 0.115195486020354
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d01251743649ae4335a060127395563cfd86d513",
          "message": "Merge pull request #573 from galacticusorg/fixOptimizeBaryonicLSS\n\nUse exclusive lock to avoid conflicts between processes",
          "timestamp": "2024-03-02T14:57:44Z",
          "tree_id": "85126c35f19e3f00a0347a88f8e1184e72a579d4",
          "url": "https://github.com/galacticusorg/galacticus/commit/d01251743649ae4335a060127395563cfd86d513"
        },
        "date": 1709402629541,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.38692333875035,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.12908611838645,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.989746819263137,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1525117191923,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 121.110580060666,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.8376386025071,
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
          "id": "6e64a03714bd1daeb9c975e5b9da4475cc06b0de",
          "message": "fix: Merge branch 'master' of github.com:galacticusorg/galacticus",
          "timestamp": "2024-03-06T10:58:09-08:00",
          "tree_id": "d6160e0734d63dc0a2fbe3642c7e9a796e38879c",
          "url": "https://github.com/galacticusorg/galacticus/commit/6e64a03714bd1daeb9c975e5b9da4475cc06b0de"
        },
        "date": 1709762328745,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 154.492,
            "unit": "seconds",
            "range": 0.139196264315283
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
          "id": "6e64a03714bd1daeb9c975e5b9da4475cc06b0de",
          "message": "fix: Merge branch 'master' of github.com:galacticusorg/galacticus",
          "timestamp": "2024-03-06T10:58:09-08:00",
          "tree_id": "d6160e0734d63dc0a2fbe3642c7e9a796e38879c",
          "url": "https://github.com/galacticusorg/galacticus/commit/6e64a03714bd1daeb9c975e5b9da4475cc06b0de"
        },
        "date": 1709762335926,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.08779492414224,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.11530259840329,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.16610887880892,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1539277972285,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 119.476649590127,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.6292371122148,
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
          "id": "c20456b4883f6ac5c9ad4d71f19afbb988029474",
          "message": "fix: Fix erroneously joined line",
          "timestamp": "2024-03-07T18:06:57Z",
          "tree_id": "0c2781c9fc1798b1bd600e15582862d043f0c906",
          "url": "https://github.com/galacticusorg/galacticus/commit/c20456b4883f6ac5c9ad4d71f19afbb988029474"
        },
        "date": 1709847120761,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 150.922,
            "unit": "seconds",
            "range": 0.548725432252944
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
          "id": "c20456b4883f6ac5c9ad4d71f19afbb988029474",
          "message": "fix: Fix erroneously joined line",
          "timestamp": "2024-03-07T18:06:57Z",
          "tree_id": "0c2781c9fc1798b1bd600e15582862d043f0c906",
          "url": "https://github.com/galacticusorg/galacticus/commit/c20456b4883f6ac5c9ad4d71f19afbb988029474"
        },
        "date": 1709847127951,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.11230747236508,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.27222287933723,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.23027926754141,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.0851332216068,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 120.479535263536,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.7627012674009,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "05526e1732a8eea0cf597557561ecd170cccbfc7",
          "message": "Merge pull request #574 from galacticusorg/featRamPressurePosition\n\nAdd a ram pressure force class which uses relative positions",
          "timestamp": "2024-03-08T04:04:16Z",
          "tree_id": "b6e25e2eeca8b4854f2df2c45009fdc0c0b8969b",
          "url": "https://github.com/galacticusorg/galacticus/commit/05526e1732a8eea0cf597557561ecd170cccbfc7"
        },
        "date": 1709881821731,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 191.995,
            "unit": "seconds",
            "range": 0.115232373927156
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "05526e1732a8eea0cf597557561ecd170cccbfc7",
          "message": "Merge pull request #574 from galacticusorg/featRamPressurePosition\n\nAdd a ram pressure force class which uses relative positions",
          "timestamp": "2024-03-08T04:04:16Z",
          "tree_id": "b6e25e2eeca8b4854f2df2c45009fdc0c0b8969b",
          "url": "https://github.com/galacticusorg/galacticus/commit/05526e1732a8eea0cf597557561ecd170cccbfc7"
        },
        "date": 1709881827664,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.34451503587741,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.79941684483739,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.979153593179563,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1028944895755,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 121.486602715398,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.6138228434182,
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
          "id": "cd809cf990fd32e4d6950fa725537232b4d2c208",
          "message": "feat: Add a merger tree operator to prune branches below a mass threshold and after the final output time",
          "timestamp": "2024-03-08T21:26:02Z",
          "tree_id": "8df396ca2f6f923003c691d78bb79edefea0ee22",
          "url": "https://github.com/galacticusorg/galacticus/commit/cd809cf990fd32e4d6950fa725537232b4d2c208"
        },
        "date": 1709944147291,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 144.013,
            "unit": "seconds",
            "range": 0.222719779095201
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
          "id": "cd809cf990fd32e4d6950fa725537232b4d2c208",
          "message": "feat: Add a merger tree operator to prune branches below a mass threshold and after the final output time",
          "timestamp": "2024-03-08T21:26:02Z",
          "tree_id": "8df396ca2f6f923003c691d78bb79edefea0ee22",
          "url": "https://github.com/galacticusorg/galacticus/commit/cd809cf990fd32e4d6950fa725537232b4d2c208"
        },
        "date": 1709944152858,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.34854287140891,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.34313405639497,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 0.984581866982154,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.1964085284445,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 122.254539630805,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.7982978939612,
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
          "id": "ae7f003c7ac6b69eef57fcf64694e92db51867e9",
          "message": "fix(style): Formatting only",
          "timestamp": "2024-03-11T14:47:15Z",
          "tree_id": "d54d75d380105fc61c7a962beb1f6854020a8686",
          "url": "https://github.com/galacticusorg/galacticus/commit/ae7f003c7ac6b69eef57fcf64694e92db51867e9"
        },
        "date": 1710179695829,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 201.665,
            "unit": "seconds",
            "range": 0.162291404575567
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
          "id": "ae7f003c7ac6b69eef57fcf64694e92db51867e9",
          "message": "fix(style): Formatting only",
          "timestamp": "2024-03-11T14:47:15Z",
          "tree_id": "d54d75d380105fc61c7a962beb1f6854020a8686",
          "url": "https://github.com/galacticusorg/galacticus/commit/ae7f003c7ac6b69eef57fcf64694e92db51867e9"
        },
        "date": 1710179703788,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.32266607057031,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.04080569577638,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.18142641227944,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2043948220685,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 120.097377734953,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.6315608659821,
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
          "id": "bb772e3b66b51f8c2e5081f5488ac8b90202d224",
          "message": "fix: Merge branch 'master' of github.com:galacticusorg/galacticus",
          "timestamp": "2024-03-11T20:41:09Z",
          "tree_id": "04f11eb7f1b836d3ffde9b1e7461006477b3af96",
          "url": "https://github.com/galacticusorg/galacticus/commit/bb772e3b66b51f8c2e5081f5488ac8b90202d224"
        },
        "date": 1710201798496,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 203.172,
            "unit": "seconds",
            "range": 0.318963320775826
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
          "id": "bb772e3b66b51f8c2e5081f5488ac8b90202d224",
          "message": "fix: Merge branch 'master' of github.com:galacticusorg/galacticus",
          "timestamp": "2024-03-11T20:41:09Z",
          "tree_id": "04f11eb7f1b836d3ffde9b1e7461006477b3af96",
          "url": "https://github.com/galacticusorg/galacticus/commit/bb772e3b66b51f8c2e5081f5488ac8b90202d224"
        },
        "date": 1710201805009,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.52607060882327,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.82404317654253,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.0695279901752,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2814549458896,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 122.141089155875,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.8576064687717,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "3ce8807226eb572cd2e69f27dff0e01ddefd9efd",
          "message": "Merge pull request #576 from galacticusorg/featTransferFunctionHalfModeSlope\n\nAdd a new transfer function class",
          "timestamp": "2024-03-12T04:50:30Z",
          "tree_id": "7a1d0cad97095dc74297e74aad36509622d468a8",
          "url": "https://github.com/galacticusorg/galacticus/commit/3ce8807226eb572cd2e69f27dff0e01ddefd9efd"
        },
        "date": 1710229937117,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 161.902,
            "unit": "seconds",
            "range": 0.194225641976331
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "3ce8807226eb572cd2e69f27dff0e01ddefd9efd",
          "message": "Merge pull request #576 from galacticusorg/featTransferFunctionHalfModeSlope\n\nAdd a new transfer function class",
          "timestamp": "2024-03-12T04:50:30Z",
          "tree_id": "7a1d0cad97095dc74297e74aad36509622d468a8",
          "url": "https://github.com/galacticusorg/galacticus/commit/3ce8807226eb572cd2e69f27dff0e01ddefd9efd"
        },
        "date": 1710229942920,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 3.39392315272724,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 4.11678503120485,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 1.1828077139497,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 24.2069788642423,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 121.088258786043,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 8.57906746448903,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "97efa132bd8578fd83bf04303d8c060c210aba3f",
          "message": "Merge pull request #583 from galacticusorg/orbitCalibration\n\nOrbit calibration",
          "timestamp": "2024-03-16T02:37:27Z",
          "tree_id": "d39ce6a0ac68f4b860f6e0f316e01a987fc7d42a",
          "url": "https://github.com/galacticusorg/galacticus/commit/97efa132bd8578fd83bf04303d8c060c210aba3f"
        },
        "date": 1710570806835,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 396.953,
            "unit": "seconds",
            "range": 0.172075855353461
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "97efa132bd8578fd83bf04303d8c060c210aba3f",
          "message": "Merge pull request #583 from galacticusorg/orbitCalibration\n\nOrbit calibration",
          "timestamp": "2024-03-16T02:37:27Z",
          "tree_id": "d39ce6a0ac68f4b860f6e0f316e01a987fc7d42a",
          "url": "https://github.com/galacticusorg/galacticus/commit/97efa132bd8578fd83bf04303d8c060c210aba3f"
        },
        "date": 1710570814826,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.88757547390109,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 59.588105071952,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 7.24064184333617,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.1627521915264,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 95.4272339234419,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.0112889620655,
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
          "id": "4372afab9f39d21a430619ea863635c6ee61f9df",
          "message": "feat: Validate using distribution of slopes",
          "timestamp": "2024-03-16T15:11:30Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/4372afab9f39d21a430619ea863635c6ee61f9df"
        },
        "date": 1710616349593,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 409.705,
            "unit": "seconds",
            "range": 0.24837974957316
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
          "id": "4372afab9f39d21a430619ea863635c6ee61f9df",
          "message": "feat: Validate using distribution of slopes",
          "timestamp": "2024-03-16T15:11:30Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/4372afab9f39d21a430619ea863635c6ee61f9df"
        },
        "date": 1710616356568,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.88175405144041,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 59.9009530765768,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 7.22447011572113,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.0828027589331,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 266.178929763363,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 11.9563492745274,
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
          "id": "8f25c914acd1feade805f306d219c5c1c21a9ac2",
          "message": "fix(style): Formatting only",
          "timestamp": "2024-03-18T16:55:21Z",
          "tree_id": "60c66da2d03823376c364c3f5d18224a1800dbe9",
          "url": "https://github.com/galacticusorg/galacticus/commit/8f25c914acd1feade805f306d219c5c1c21a9ac2"
        },
        "date": 1710795461080,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 431.101,
            "unit": "seconds",
            "range": 0.18074540105038
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
          "id": "8f25c914acd1feade805f306d219c5c1c21a9ac2",
          "message": "fix(style): Formatting only",
          "timestamp": "2024-03-18T16:55:21Z",
          "tree_id": "60c66da2d03823376c364c3f5d18224a1800dbe9",
          "url": "https://github.com/galacticusorg/galacticus/commit/8f25c914acd1feade805f306d219c5c1c21a9ac2"
        },
        "date": 1710795467731,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.93515494723007,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 59.6209319419666,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 7.21794014356623,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.2463789835063,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 96.7083361730445,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.0083753197272,
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
          "id": "7de5e4e6e3b6537ea986f05817e4ac54b8ae1604",
          "message": "fix: Catch alternate URLs for updated bibliography entries",
          "timestamp": "2024-03-19T23:53:15Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/7de5e4e6e3b6537ea986f05817e4ac54b8ae1604"
        },
        "date": 1710906928422,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 459.464,
            "unit": "seconds",
            "range": 0.28465839177874
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
          "id": "7de5e4e6e3b6537ea986f05817e4ac54b8ae1604",
          "message": "fix: Catch alternate URLs for updated bibliography entries",
          "timestamp": "2024-03-19T23:53:15Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/7de5e4e6e3b6537ea986f05817e4ac54b8ae1604"
        },
        "date": 1710906936224,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.87130688524518,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 60.7474045336967,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 7.29234262983404,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.1950253708527,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 272.026043626266,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 11.9517585759031,
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
          "id": "3ac102f22a46dd05613771c69899f2ca5b2b36c1",
          "message": "fix: Correct LaTeX journal reference commands",
          "timestamp": "2024-03-20T05:22:40Z",
          "tree_id": "5c3357ad8a4da33c36e1a354162d2bac4a0fa9b0",
          "url": "https://github.com/galacticusorg/galacticus/commit/3ac102f22a46dd05613771c69899f2ca5b2b36c1"
        },
        "date": 1710926315646,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 426.974,
            "unit": "seconds",
            "range": 0.401467806925042
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
          "id": "3ac102f22a46dd05613771c69899f2ca5b2b36c1",
          "message": "fix: Correct LaTeX journal reference commands",
          "timestamp": "2024-03-20T05:22:40Z",
          "tree_id": "5c3357ad8a4da33c36e1a354162d2bac4a0fa9b0",
          "url": "https://github.com/galacticusorg/galacticus/commit/3ac102f22a46dd05613771c69899f2ca5b2b36c1"
        },
        "date": 1710926321510,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 4.82801340690533,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 60.7872387045732,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 7.27422847946757,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.0716444601247,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 264.193357071225,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.0664898955798,
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
          "id": "759d2ea43ba78b72c1d2df543762181dbf715d67",
          "message": "fix: Correct file name for artifact upload",
          "timestamp": "2024-03-25T20:22:32Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/759d2ea43ba78b72c1d2df543762181dbf715d67"
        },
        "date": 1711413298959,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 525.739,
            "unit": "seconds",
            "range": 0.320501014060502
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
          "id": "759d2ea43ba78b72c1d2df543762181dbf715d67",
          "message": "fix: Correct file name for artifact upload",
          "timestamp": "2024-03-25T20:22:32Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/759d2ea43ba78b72c1d2df543762181dbf715d67"
        },
        "date": 1711413306464,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 5.61863847209972,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 60.4343372654477,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 6.73079520058824,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.5959359341889,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 84.1427837362023,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.9463374297978,
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
          "id": "c4658de9a148da9c2519051f7e2f5a43b2273a3d",
          "message": "feat: Add a title for the build profile page",
          "timestamp": "2024-03-26T07:55:11-07:00",
          "tree_id": "4670ece7f17d8544c7eede5a64105f14dc756336",
          "url": "https://github.com/galacticusorg/galacticus/commit/c4658de9a148da9c2519051f7e2f5a43b2273a3d"
        },
        "date": 1711479308322,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 410.576,
            "unit": "seconds",
            "range": 0.410317438098608
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
          "id": "c4658de9a148da9c2519051f7e2f5a43b2273a3d",
          "message": "feat: Add a title for the build profile page",
          "timestamp": "2024-03-26T07:55:11-07:00",
          "tree_id": "4670ece7f17d8544c7eede5a64105f14dc756336",
          "url": "https://github.com/galacticusorg/galacticus/commit/c4658de9a148da9c2519051f7e2f5a43b2273a3d"
        },
        "date": 1711479315917,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 5.66003342832257,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 59.5761158781408,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 6.64261990115022,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.6021795749817,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 85.8953515418612,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.9478967805266,
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
          "id": "5d8af8db7fbe60b96437d70092c152863e6c147d",
          "message": "fix(style): Formatting only",
          "timestamp": "2024-03-27T14:55:16Z",
          "tree_id": "9e5fe3de420be62d873888df77348eb8bd30bc85",
          "url": "https://github.com/galacticusorg/galacticus/commit/5d8af8db7fbe60b96437d70092c152863e6c147d"
        },
        "date": 1711586082606,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 496.664,
            "unit": "seconds",
            "range": 0.28647582793639
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
          "id": "5d8af8db7fbe60b96437d70092c152863e6c147d",
          "message": "fix(style): Formatting only",
          "timestamp": "2024-03-27T14:55:16Z",
          "tree_id": "9e5fe3de420be62d873888df77348eb8bd30bc85",
          "url": "https://github.com/galacticusorg/galacticus/commit/5d8af8db7fbe60b96437d70092c152863e6c147d"
        },
        "date": 1711586091652,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 5.5890471969813,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 59.6386065136656,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 6.61214189156761,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.5970475972473,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 68.9227199659161,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.0281805844055,
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
          "id": "a2807179bd76a02f1216608a7eb7e1b8308fa993",
          "message": "feat: Add helpful error message when failing to parse a fully-specified merger tree file",
          "timestamp": "2024-03-28T18:08:28-07:00",
          "tree_id": "c24c8b3b79152a3757845f20bf8c28a030bb1339",
          "url": "https://github.com/galacticusorg/galacticus/commit/a2807179bd76a02f1216608a7eb7e1b8308fa993"
        },
        "date": 1711690701347,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 431.406,
            "unit": "seconds",
            "range": 0.299550329658257
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
          "id": "a2807179bd76a02f1216608a7eb7e1b8308fa993",
          "message": "feat: Add helpful error message when failing to parse a fully-specified merger tree file",
          "timestamp": "2024-03-28T18:08:28-07:00",
          "tree_id": "c24c8b3b79152a3757845f20bf8c28a030bb1339",
          "url": "https://github.com/galacticusorg/galacticus/commit/a2807179bd76a02f1216608a7eb7e1b8308fa993"
        },
        "date": 1711690709541,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 5.63614659958942,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 59.0812750794803,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 6.74837119706251,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.599825527664,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 85.3178600409857,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.9665121303074,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e34b78a68f5847cf11983f1fd794f3ffe739f1d2",
          "message": "Merge pull request #595 from galacticusorg/fixUniverseEvolution\n\nAvoid deadlocking when processing universe events",
          "timestamp": "2024-03-30T15:21:22Z",
          "tree_id": "727519e6bc024b421e20b21916b2c3378f206ecc",
          "url": "https://github.com/galacticusorg/galacticus/commit/e34b78a68f5847cf11983f1fd794f3ffe739f1d2"
        },
        "date": 1711826452131,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 468.077,
            "unit": "seconds",
            "range": 0.290764681483689
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e34b78a68f5847cf11983f1fd794f3ffe739f1d2",
          "message": "Merge pull request #595 from galacticusorg/fixUniverseEvolution\n\nAvoid deadlocking when processing universe events",
          "timestamp": "2024-03-30T15:21:22Z",
          "tree_id": "727519e6bc024b421e20b21916b2c3378f206ecc",
          "url": "https://github.com/galacticusorg/galacticus/commit/e34b78a68f5847cf11983f1fd794f3ffe739f1d2"
        },
        "date": 1711826460105,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 5.6013943982179,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 60.2057566977968,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 6.71395685092119,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.5822707435313,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 85.4066598840991,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.9374714729122,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9a84f8498ff587af48f4d1d4d90a77622207bab8",
          "message": "Merge pull request #596 from galacticusorg/fixEnumerationDescribe\n\nFix generation of enumeration description functions",
          "timestamp": "2024-03-31T16:08:06Z",
          "tree_id": "8b2a44d5cbd9026cfbe313325560e8d709609125",
          "url": "https://github.com/galacticusorg/galacticus/commit/9a84f8498ff587af48f4d1d4d90a77622207bab8"
        },
        "date": 1711915685670,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 468.475,
            "unit": "seconds",
            "range": 0.326898302221995
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9a84f8498ff587af48f4d1d4d90a77622207bab8",
          "message": "Merge pull request #596 from galacticusorg/fixEnumerationDescribe\n\nFix generation of enumeration description functions",
          "timestamp": "2024-03-31T16:08:06Z",
          "tree_id": "8b2a44d5cbd9026cfbe313325560e8d709609125",
          "url": "https://github.com/galacticusorg/galacticus/commit/9a84f8498ff587af48f4d1d4d90a77622207bab8"
        },
        "date": 1711915695697,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 5.63075698315696,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 61.0438827724562,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 6.76476859420295,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.672815042201,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 68.1697712361321,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.9446590364978,
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
          "id": "ad57e23c33a9b4cde3ab619dee4a5c61c166c613",
          "message": "fix: Ensure objects are allocated before attempting to deallocate",
          "timestamp": "2024-04-02T11:44:35-07:00",
          "tree_id": "3048857c251de0375e874fac214fb0d2605cd580",
          "url": "https://github.com/galacticusorg/galacticus/commit/ad57e23c33a9b4cde3ab619dee4a5c61c166c613"
        },
        "date": 1712107823829,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 433.314,
            "unit": "seconds",
            "range": 0.287033795916627
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
          "id": "ad57e23c33a9b4cde3ab619dee4a5c61c166c613",
          "message": "fix: Ensure objects are allocated before attempting to deallocate",
          "timestamp": "2024-04-02T11:44:35-07:00",
          "tree_id": "3048857c251de0375e874fac214fb0d2605cd580",
          "url": "https://github.com/galacticusorg/galacticus/commit/ad57e23c33a9b4cde3ab619dee4a5c61c166c613"
        },
        "date": 1712107833932,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 6.18724334517588,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 56.2697199508759,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 6.24532850846259,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.1572973439512,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 54.7197392954377,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 15.8033053792017,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "49bb79e525b43544af3c0342b3199b8e3f656141",
          "message": "Merge pull request #598 from galacticusorg/fixJiang2014OrbitsFileLock\n\nFix file locking in the `virialOrbitJiang2014` class",
          "timestamp": "2024-04-05T02:44:00Z",
          "tree_id": "9dec9b16322946f41101f690f9b298053c8dd390",
          "url": "https://github.com/galacticusorg/galacticus/commit/49bb79e525b43544af3c0342b3199b8e3f656141"
        },
        "date": 1712299619224,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 470.891,
            "unit": "seconds",
            "range": 0.589042358411205
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "49bb79e525b43544af3c0342b3199b8e3f656141",
          "message": "Merge pull request #598 from galacticusorg/fixJiang2014OrbitsFileLock\n\nFix file locking in the `virialOrbitJiang2014` class",
          "timestamp": "2024-04-05T02:44:00Z",
          "tree_id": "9dec9b16322946f41101f690f9b298053c8dd390",
          "url": "https://github.com/galacticusorg/galacticus/commit/49bb79e525b43544af3c0342b3199b8e3f656141"
        },
        "date": 1712299627214,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 5.59713916554942,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 59.4099318111514,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 6.56957983139979,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.6493928402953,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 83.9809295598752,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.9101607593046,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "466c04b2b741a22d51551deb51143752b68dbcd7",
          "message": "Merge pull request #599 from galacticusorg/featDocClassDefaults\n\nSpecify class defaults in the documentation",
          "timestamp": "2024-04-10T04:55:20Z",
          "tree_id": "78ab2e2e8819507f92d5046852c62cfaf7f8dda5",
          "url": "https://github.com/galacticusorg/galacticus/commit/466c04b2b741a22d51551deb51143752b68dbcd7"
        },
        "date": 1712739052707,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 353.438,
            "unit": "seconds",
            "range": 0.200906943633796
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "466c04b2b741a22d51551deb51143752b68dbcd7",
          "message": "Merge pull request #599 from galacticusorg/featDocClassDefaults\n\nSpecify class defaults in the documentation",
          "timestamp": "2024-04-10T04:55:20Z",
          "tree_id": "78ab2e2e8819507f92d5046852c62cfaf7f8dda5",
          "url": "https://github.com/galacticusorg/galacticus/commit/466c04b2b741a22d51551deb51143752b68dbcd7"
        },
        "date": 1712739065598,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 5.5414360458036,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 60.2451611449146,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 6.71486520749867,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.5962878225013,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 84.6655862962841,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.9305507730765,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ace5ce0714ef1f3e486e4cf3d7d2a643b828e999",
          "message": "Merge pull request #600 from galacticusorg/featValidateBaryonicSuppression\n\nAdd a validation test for baryonic suppression of structure formation",
          "timestamp": "2024-04-13T16:36:37Z",
          "tree_id": "68fdbd599d3956a2242b5a837e1cbc29b6dc8bbb",
          "url": "https://github.com/galacticusorg/galacticus/commit/ace5ce0714ef1f3e486e4cf3d7d2a643b828e999"
        },
        "date": 1713044265430,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 381.857,
            "unit": "seconds",
            "range": 0.298868031068483
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ace5ce0714ef1f3e486e4cf3d7d2a643b828e999",
          "message": "Merge pull request #600 from galacticusorg/featValidateBaryonicSuppression\n\nAdd a validation test for baryonic suppression of structure formation",
          "timestamp": "2024-04-13T16:36:37Z",
          "tree_id": "68fdbd599d3956a2242b5a837e1cbc29b6dc8bbb",
          "url": "https://github.com/galacticusorg/galacticus/commit/ace5ce0714ef1f3e486e4cf3d7d2a643b828e999"
        },
        "date": 1713044272927,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 5.61693001797611,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 59.8123268568159,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 6.51504443460744,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.6201371984835,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 86.1381916808366,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 13.004689344697,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5c8dd069abdd6e76f85bd796c98848db2d3b86b7",
          "message": "Merge pull request #604 from galacticusorg/fixGeometryDatasets\n\nUpdate the URL for SDSS DR7 geometry data file",
          "timestamp": "2024-04-19T23:10:12Z",
          "tree_id": "32c13fd5d1c88e3e5c03a4ddebb6045f6cb2a491",
          "url": "https://github.com/galacticusorg/galacticus/commit/5c8dd069abdd6e76f85bd796c98848db2d3b86b7"
        },
        "date": 1713586334195,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 459.418,
            "unit": "seconds",
            "range": 0.258355568924718
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5c8dd069abdd6e76f85bd796c98848db2d3b86b7",
          "message": "Merge pull request #604 from galacticusorg/fixGeometryDatasets\n\nUpdate the URL for SDSS DR7 geometry data file",
          "timestamp": "2024-04-19T23:10:12Z",
          "tree_id": "32c13fd5d1c88e3e5c03a4ddebb6045f6cb2a491",
          "url": "https://github.com/galacticusorg/galacticus/commit/5c8dd069abdd6e76f85bd796c98848db2d3b86b7"
        },
        "date": 1713586342231,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 5.62649304891523,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 59.1555655934733,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 6.52423625555536,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.6151349430879,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 85.3153696189279,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.9471137503404,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c986aaafcb98f6b42ff0a61bce57444185c5b940",
          "message": "Merge pull request #605 from galacticusorg/featConcentrationTests\n\nExpand the tests of concentration models",
          "timestamp": "2024-04-25T14:13:09Z",
          "tree_id": "6546ebdbf1715cbf79460a1c573613bfdb4ac007",
          "url": "https://github.com/galacticusorg/galacticus/commit/c986aaafcb98f6b42ff0a61bce57444185c5b940"
        },
        "date": 1714072542628,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 351.195,
            "unit": "seconds",
            "range": 0.179411538085769
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c986aaafcb98f6b42ff0a61bce57444185c5b940",
          "message": "Merge pull request #605 from galacticusorg/featConcentrationTests\n\nExpand the tests of concentration models",
          "timestamp": "2024-04-25T14:13:09Z",
          "tree_id": "6546ebdbf1715cbf79460a1c573613bfdb4ac007",
          "url": "https://github.com/galacticusorg/galacticus/commit/c986aaafcb98f6b42ff0a61bce57444185c5b940"
        },
        "date": 1714072549734,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 5.66736064296011,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 60.0498933129212,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 6.62713864029349,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.6010571059684,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 66.5542084520898,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.90265907595,
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
          "id": "7e34697042dac0564a9235b9762ebd29a0420bb9",
          "message": "feat: Improve help message",
          "timestamp": "2024-04-30T14:23:10Z",
          "tree_id": "be372c3aad223084dbad67f2429e8835bf2adcfd",
          "url": "https://github.com/galacticusorg/galacticus/commit/7e34697042dac0564a9235b9762ebd29a0420bb9"
        },
        "date": 1714515227063,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 411.043,
            "unit": "seconds",
            "range": 0.187483599289771
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
          "id": "7e34697042dac0564a9235b9762ebd29a0420bb9",
          "message": "feat: Improve help message",
          "timestamp": "2024-04-30T14:23:10Z",
          "tree_id": "be372c3aad223084dbad67f2429e8835bf2adcfd",
          "url": "https://github.com/galacticusorg/galacticus/commit/7e34697042dac0564a9235b9762ebd29a0420bb9"
        },
        "date": 1714515233901,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 6.21765469232009,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 56.2398731093509,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 6.28963387714261,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.1862060820028,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 82.2694946998547,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 15.873914066425,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9e6a229062512e24a763699e077ff6fd844f3cd8",
          "message": "Merge pull request #607 from galacticusorg/fixTruncatedPowerTreeBuild\n\nFix truncated power tree build",
          "timestamp": "2024-05-03T12:56:09Z",
          "tree_id": "f20baa8b4200480f055a3d954784e6dcdac45320",
          "url": "https://github.com/galacticusorg/galacticus/commit/9e6a229062512e24a763699e077ff6fd844f3cd8"
        },
        "date": 1714759110260,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 354.715,
            "unit": "seconds",
            "range": 0.331714485660458
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9e6a229062512e24a763699e077ff6fd844f3cd8",
          "message": "Merge pull request #607 from galacticusorg/fixTruncatedPowerTreeBuild\n\nFix truncated power tree build",
          "timestamp": "2024-05-03T12:56:09Z",
          "tree_id": "f20baa8b4200480f055a3d954784e6dcdac45320",
          "url": "https://github.com/galacticusorg/galacticus/commit/9e6a229062512e24a763699e077ff6fd844f3cd8"
        },
        "date": 1714759117693,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Likelihood - localGroupMassMetallicityRelation",
            "value": 5.63281687588957,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassSizeRelation",
            "value": 59.2919290028983,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupMassVelocityDispersionRelation",
            "value": 6.65173005418985,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupOccupationFraction",
            "value": 16.6226151604118,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassFunction",
            "value": 84.351349161197,
            "unit": "-logℒ"
          },
          {
            "name": "Milky Way model - Likelihood - localGroupStellarMassHaloMassRelation",
            "value": 12.9216345820983,
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
          "id": "8b2d9cc16747e1474f9541a994a8e22d6cfc5f4d",
          "message": "fix: Ensure variables are always initialized",
          "timestamp": "2024-05-03T16:39:46-07:00",
          "tree_id": "f36ddf8576b86ae3483f146ecce13142970abacf",
          "url": "https://github.com/galacticusorg/galacticus/commit/8b2d9cc16747e1474f9541a994a8e22d6cfc5f4d"
        },
        "date": 1714798061256,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Milky Way model - Wall Time",
            "value": 366.157,
            "unit": "seconds",
            "range": 0.298141073983872
          }
        ]
      }
    ]
  }
}