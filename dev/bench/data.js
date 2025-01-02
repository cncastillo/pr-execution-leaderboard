window.BENCHMARK_DATA = {
  "lastUpdate": 1735852846236,
  "repoUrl": "https://github.com/Stockless/pr-execution-leaderboard",
  "entries": {
    "Julia benchmark result": [
      {
        "commit": {
          "author": {
            "name": "Stockless",
            "username": "Stockless"
          },
          "committer": {
            "name": "Stockless",
            "username": "Stockless"
          },
          "id": "b0c4aea62650ba6bd370eb538174ef6cea4accdb",
          "message": "test1",
          "timestamp": "2025-01-02T21:19:50Z",
          "url": "https://github.com/Stockless/pr-execution-leaderboard/pull/8/commits/b0c4aea62650ba6bd370eb538174ef6cea4accdb"
        },
        "date": 1735852845813,
        "tool": "julia",
        "benches": [
          {
            "name": "function1/Stockless",
            "value": 9267,
            "unit": "ns",
            "extra": "gctime=727.8853190329343\nmemory=0\nallocs=0\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "function2/Stockless",
            "value": 1.563,
            "unit": "ns",
            "extra": "gctime=0.32621599943593943\nmemory=0\nallocs=0\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1000,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}