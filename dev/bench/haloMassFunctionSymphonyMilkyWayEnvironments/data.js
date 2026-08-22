window.BENCHMARK_DATA = {
  "lastUpdate": 1787441977720,
  "repoUrl": "https://github.com/galacticusorg/galacticus",
  "entries": {
    "Halo mass function validation (Symphony Milky Way environments)": [
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
        "date": 1784610277050,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1671.526349571755,
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
        "date": 1784693923245,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1671.5264492852175,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1785361947581,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1671.526271944108,
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
        "date": 1785795508190,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1671.526349571755,
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
        "date": 1785934660948,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1671.5263802516524,
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
        "date": 1786131069799,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1671.5264646779194,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
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
        "date": 1786200396339,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1671.5264646779194,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9e04f6f256cbf0ca9dcb567b478c773848082e61",
          "message": "Merge pull request #1347 from galacticusorg/fix/ethos-extended-slope-analytic-derivative\n\nfix(powerSpectrum): evaluate the ETHOS extended slope analytically",
          "timestamp": "2026-08-09T21:58:23Z",
          "tree_id": "00edc85b0c464f692594354d58b57a99179054be",
          "url": "https://github.com/galacticusorg/galacticus/commit/9e04f6f256cbf0ca9dcb567b478c773848082e61"
        },
        "date": 1786334912755,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1672.1924955489462,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "f62df0c3b1d868c052eedef79c20488955a738c9",
          "message": "Merge pull request #1366 from galacticusorg/update/bibliography\n\nfix: Update bibliography records",
          "timestamp": "2026-08-12T22:36:04Z",
          "tree_id": "b8f41bf1c8420cb8858a024e36483c41327a1f7d",
          "url": "https://github.com/galacticusorg/galacticus/commit/f62df0c3b1d868c052eedef79c20488955a738c9"
        },
        "date": 1786601576062,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1672.1925378606736,
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
          "id": "ea506749aa29eefea02d26d50af0db3050d203b9",
          "message": "chore(parameters): migrate parameter files to latest revision",
          "timestamp": "2026-08-13T08:11:28-07:00",
          "tree_id": "55c98049637be9f342b359850646090e1c40a53c",
          "url": "https://github.com/galacticusorg/galacticus/commit/ea506749aa29eefea02d26d50af0db3050d203b9"
        },
        "date": 1786658483740,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1672.1925385009288,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "0f36c1956a768da25d994e061f65826a137b2d48",
          "message": "Merge pull request #1371 from galacticusorg/fix/warn-quadratic-append\n\nfix(error): make recording of warnings independent of the number already issued",
          "timestamp": "2026-08-14T14:45:48Z",
          "tree_id": "f55822ad989c715b241aca7f1a6d704d5636a77f",
          "url": "https://github.com/galacticusorg/galacticus/commit/0f36c1956a768da25d994e061f65826a137b2d48"
        },
        "date": 1786741023678,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1672.192548607109,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9b6b0b129609b20637577dad446f2cdc9013b325",
          "message": "Merge pull request #1380 from galacticusorg/chore/gh-pages-thresholds-2026-08-14\n\nRe-anchor gh-pages analysis thresholds after the 2026-08-14 run",
          "timestamp": "2026-08-15T15:04:23Z",
          "tree_id": "2622d555d1bb8ae1967c3cf5362eafb5424908f2",
          "url": "https://github.com/galacticusorg/galacticus/commit/9b6b0b129609b20637577dad446f2cdc9013b325"
        },
        "date": 1786826941485,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1672.1925682155456,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ad3a4b3b14ecdd61e790a46a8db65373d9510e17",
          "message": "Merge pull request #1394 from galacticusorg/fix/derived-type-opener-regex\n\nfix(build): recognize `bind(c)` and spaced `extends(...)` in derived-type openers",
          "timestamp": "2026-08-18T18:06:32Z",
          "tree_id": "85b1feffe7df956af19e61e18157f4a7b47879ca",
          "url": "https://github.com/galacticusorg/galacticus/commit/ad3a4b3b14ecdd61e790a46a8db65373d9510e17"
        },
        "date": 1787093261019,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1672.192556479688,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "544daf1fc5eafad397a6143d479501f1129ee2a3",
          "message": "Merge pull request #1398 from galacticusorg/fix/preprocessor-generator-prerequisites\n\nfix(build): re-preprocess when the preprocessor's own Python sources change",
          "timestamp": "2026-08-18T23:38:46Z",
          "tree_id": "0514f3253b97556c90ddaa92e9002a55d383bd42",
          "url": "https://github.com/galacticusorg/galacticus/commit/544daf1fc5eafad397a6143d479501f1129ee2a3"
        },
        "date": 1787117918895,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1672.192473896027,
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
          "id": "0f386049a46fce7fce1c459a8c53cd79fa5e8c0b",
          "message": "chore(ci): raise the COZMIC WDM 3keV X1 benchmark threshold to 750%",
          "timestamp": "2026-08-20T06:10:12-07:00",
          "tree_id": "ca1c1aabceb526f2ff58815cc8c884590bfd43c5",
          "url": "https://github.com/galacticusorg/galacticus/commit/0f386049a46fce7fce1c459a8c53cd79fa5e8c0b"
        },
        "date": 1787250773553,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1672.1925682155456,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5e69da91e1dfdf8f745686f542e068deefbced61",
          "message": "Merge pull request #1413 from galacticusorg/feature/error-report-signal-handlers\n\nCall registered signal handlers on a deliberate fatal error",
          "timestamp": "2026-08-21T19:55:27Z",
          "tree_id": "804b1a1033550465b02702b2355fb9d935f37c6f",
          "url": "https://github.com/galacticusorg/galacticus/commit/5e69da91e1dfdf8f745686f542e068deefbced61"
        },
        "date": 1787363905898,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1672.1925682155456,
            "unit": "-logℒ"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "abensonca@gmail.com",
            "name": "Andrew Benson",
            "username": "abensonca"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "0ab8e753dde81ca158c013508a9f49681f00bd41",
          "message": "Merge pull request #1415 from galacticusorg/update/bibliography\n\nfix: Update bibliography records",
          "timestamp": "2026-08-22T18:51:24Z",
          "tree_id": "f23116617f68f078e71a6259f52e33a056aedef4",
          "url": "https://github.com/galacticusorg/galacticus/commit/0ab8e753dde81ca158c013508a9f49681f00bd41"
        },
        "date": 1787441976594,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Halo mass function - Likelihood - Symphony MilkyWay CDM resolutionX1 z=0.000 environments",
            "value": 1672.1925348030409,
            "unit": "-logℒ"
          }
        ]
      }
    ]
  }
}