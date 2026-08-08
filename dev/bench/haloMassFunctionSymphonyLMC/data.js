window.BENCHMARK_DATA = {
  "lastUpdate": 1786200373637,
  "repoUrl": "https://github.com/galacticusorg/galacticus",
  "entries": {
    "Halo mass function validation (Symphony LMC)": [
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
          "id": "7b8382082e295748f204b4a71e38bc1f6d18f1aa",
          "message": "ci: replace obsoleted test with new test",
          "timestamp": "2026-07-20T14:43:48-07:00",
          "tree_id": "f58631c70bbedd47f7f0513c0aca211629521180",
          "url": "https://github.com/galacticusorg/galacticus/commit/7b8382082e295748f204b4a71e38bc1f6d18f1aa"
        },
        "date": 1784610254404,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony LMC CDM resolutionX1 z=0.000 (38 realizations)",
            "value": 2002.6205948962652,
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
          "id": "7ccee793c100d43f06a33d7ab4cdf7a59b8536c8",
          "message": "ci: update gh-pages analysis thresholds\n\nAdd halo mass function metrics (COZMIC, MDPL, Symphony) and refresh\ncurrent/threshold values from the latest analysis run.",
          "timestamp": "2026-07-21T15:57:10-07:00",
          "tree_id": "7dc49e42433465665ece78a3f1bc36f3bb74c749",
          "url": "https://github.com/galacticusorg/galacticus/commit/7ccee793c100d43f06a33d7ab4cdf7a59b8536c8"
        },
        "date": 1784693892432,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony LMC CDM resolutionX1 z=0.000 (38 realizations)",
            "value": 2002.6212187274411,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6cd4316ed766428efed5b6330b379ec7c155abeb",
          "message": "Merge pull request #1306 from galacticusorg/satelliteFixesOnMaster\n\nfix: two dynamical friction bugs affecting subhalo counts, and the associated parameter migrations",
          "timestamp": "2026-07-29T15:07:04Z",
          "tree_id": "8d2a60073905259b45b1211dcc08447c8543e9f0",
          "url": "https://github.com/galacticusorg/galacticus/commit/6cd4316ed766428efed5b6330b379ec7c155abeb"
        },
        "date": 1785361919688,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony LMC CDM resolutionX1 z=0.000 (38 realizations)",
            "value": 2002.6211381265223,
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
          "id": "ae9cc618d0868833023d94cfd8288f42b5cffbe2",
          "message": "feat(ssps): add fixed stellar spectrum metallicity option for SSP tables\n\nAllows the ionizing stellar spectrum to be held fixed at a chosen\nmetallicity while the HII region gas metallicity continues to vary\nacross the grid, so that the effects of gas-phase abundances can be\nisolated from those of the stellar spectrum. The ionizing photon rate\nper unit stellar mass, the Cloudy normalization, and the validation\ncheck all follow the fixed spectrum, and the requested and adopted\nmetallicities are recorded in the output file.\n\nAlso gives `--abundanceAdjust` an explicit empty-list default, since\n`action='append'` otherwise leaves it as `None`.",
          "timestamp": "2026-08-03T07:58:19-07:00",
          "tree_id": "46ce2b23f221f6406a2b99cbf53fad18a17105bf",
          "url": "https://github.com/galacticusorg/galacticus/commit/ae9cc618d0868833023d94cfd8288f42b5cffbe2"
        },
        "date": 1785795478025,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony LMC CDM resolutionX1 z=0.000 (38 realizations)",
            "value": 2002.6211381265223,
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
          "id": "ab14075604f5713e98ac599a6c718606383c49dc",
          "message": "docs(physics): document mass distributions and potential zero points\n\nReplace the outdated \"galactic structure functions\" description with an\naccount of the `massDistribution` interface which replaced it, covering\ncomponent/mass type selection, the dark matter mass convention, and the\nfact that potential zero points are now set per mass distribution rather\nthan globally offset to the virial radius.",
          "timestamp": "2026-08-03T18:26:15-07:00",
          "tree_id": "28f74df2ad5bef46d25cae998246d165cc0306a0",
          "url": "https://github.com/galacticusorg/galacticus/commit/ab14075604f5713e98ac599a6c718606383c49dc"
        },
        "date": 1785934631010,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony LMC CDM resolutionX1 z=0.000 (38 realizations)",
            "value": 2002.6212179965205,
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
          "id": "64807e850821afdd9d0516390e1dd5bcecd70943",
          "message": "chore(ci): update gh-pages analysis thresholds for COZMIC models",
          "timestamp": "2026-08-06T08:10:40-07:00",
          "tree_id": "898a7a3024115bd46f58dd3d6344acfc739439ee",
          "url": "https://github.com/galacticusorg/galacticus/commit/64807e850821afdd9d0516390e1dd5bcecd70943"
        },
        "date": 1786131043679,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony LMC CDM resolutionX1 z=0.000 (38 realizations)",
            "value": 2002.6205998228145,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6703a1e9545f194745f78ccd1e4dc9384ac1ce3d",
          "message": "Merge pull request #1335 from galacticusorg/pr/farahi-pinning\n\nfeat(excursionSets): pin and reuse the Farahi first crossing tabulations",
          "timestamp": "2026-08-08T05:00:47Z",
          "tree_id": "657368afa0a13bcd1b1ffa9d85374a28c66839d4",
          "url": "https://github.com/galacticusorg/galacticus/commit/6703a1e9545f194745f78ccd1e4dc9384ac1ce3d"
        },
        "date": 1786200372547,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony LMC CDM resolutionX1 z=0.000 (38 realizations)",
            "value": 2002.6211442665935,
            "unit": "-logℒ"
          }
        ]
      }
    ]
  }
}