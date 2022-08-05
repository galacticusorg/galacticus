window.BENCHMARK_DATA = {
  "lastUpdate": 1659740461968,
  "repoUrl": "https://github.com/galacticusorg/galacticus",
  "entries": {
    "Benchmark": [
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
          "id": "3c857a66d3bf1eeccf294e6b0759544d8a044b91",
          "message": "fix: Correct push task in benchmark workflow",
          "timestamp": "2022-08-05T18:38:22Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/3c857a66d3bf1eeccf294e6b0759544d8a044b91"
        },
        "date": 1659728112089,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Wall Time",
            "value": 41.239,
            "unit": "seconds",
            "range": 1.90165953840311
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
          "id": "9cdf8f4b349912c271b50e166dcca29c7f80102f",
          "message": "feat: Add a validation workflow\n\nIncludes validation of the subhalo orbital evolution model",
          "timestamp": "2022-08-05T21:43:47Z",
          "url": "https://github.com/galacticusorg/galacticus/commit/9cdf8f4b349912c271b50e166dcca29c7f80102f"
        },
        "date": 1659740460494,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloMassFunction",
            "value": 59.6374034521766,
            "unit": "|logℒ|"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloRadialDistribution",
            "value": 25.9438445610862,
            "unit": "|logℒ|"
          },
          {
            "name": "Dark Matter Only Subhalos - Likelihood - subhaloVelocityMaximumMean",
            "value": 1.79769313486232e+292,
            "unit": "|logℒ|"
          }
        ]
      }
    ]
  }
}