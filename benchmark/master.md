# Benchmark Report for *StableDQMC*

## Job Properties
* Time of benchmark: 22 May 2019 - 13:16
* Package commit: dirty
* Julia commit: 80516c
* Julia command flags: `-O3,--project=.`
* Environment variables: `OMP_NUM_THREADS => 4` `MKL_NUM_THREADS => 4` `JULIA_NUM_THREADS => 1`

## Results
Below is a table of this job's results, obtained by running the benchmarks.
The values listed in the `ID` column have the structure `[parent_group, child_group, ..., key]`, and can be used to
index into the BaseBenchmarks suite to retrieve the corresponding benchmarks.
The percentages accompanying time and memory values in the below table are noise tolerances. The "true"
time/memory value for a given benchmark is expected to fall within this percentage of the reported value.
An empty cell means that the value was zero.

| ID                                                                    | time            | GC time | memory          | allocations |
|-----------------------------------------------------------------------|----------------:|--------:|----------------:|------------:|
| `["trigonometry", "circular", "(\"cos\", 0.0)"]`                      |   3.399 ns (5%) |         |                 |             |
| `["trigonometry", "circular", "(\"cos\", π = 3.1415926535897...)"]`   |   0.001 ns (5%) |         |                 |             |
| `["trigonometry", "circular", "(\"sin\", 0.0)"]`                      |   3.300 ns (5%) |         |                 |             |
| `["trigonometry", "circular", "(\"sin\", π = 3.1415926535897...)"]`   |   0.001 ns (5%) |         |                 |             |
| `["trigonometry", "circular", "(\"tan\", 0.0)"]`                      |   3.299 ns (5%) |         |                 |             |
| `["trigonometry", "circular", "(\"tan\", π = 3.1415926535897...)"]`   |   0.001 ns (5%) |         |                 |             |
| `["trigonometry", "hyperbolic", "(\"cos\", 0.0)"]`                    |   3.299 ns (5%) |         |                 |             |
| `["trigonometry", "hyperbolic", "(\"cos\", π = 3.1415926535897...)"]` |   0.001 ns (5%) |         |                 |             |
| `["trigonometry", "hyperbolic", "(\"sin\", 0.0)"]`                    |   3.299 ns (5%) |         |                 |             |
| `["trigonometry", "hyperbolic", "(\"sin\", π = 3.1415926535897...)"]` |   0.001 ns (5%) |         |                 |             |
| `["trigonometry", "hyperbolic", "(\"tan\", 0.0)"]`                    |   3.299 ns (5%) |         |                 |             |
| `["trigonometry", "hyperbolic", "(\"tan\", π = 3.1415926535897...)"]` |   0.001 ns (5%) |         |                 |             |
| `["utf8", "join"]`                                                    | 156.431 ms (5%) |         | 156.27 MiB (1%) |          21 |
| `["utf8", "replace"]`                                                 | 144.701 μs (5%) |         |  12.09 KiB (1%) |           7 |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["trigonometry", "circular"]`
- `["trigonometry", "hyperbolic"]`
- `["utf8"]`

## Julia versioninfo
```
Julia Version 1.1.0
Commit 80516ca202 (2019-01-21 21:24 UTC)
Platform Info:
  OS: Windows (x86_64-w64-mingw32)
      Microsoft Windows [Version 10.0.17763.475]
  CPU: Intel(R) Core(TM) i5-6600 CPU @ 3.30GHz: 
              speed         user         nice          sys         idle          irq
       #1  3312 MHz  185570453            0    108983312    1011484359      4243656  ticks
       #2  3312 MHz  162304734            0     91267093    1052466078      1758203  ticks
       #3  3312 MHz  198191343            0     97118328    1010728234      2118921  ticks
       #4  3312 MHz  202511906            0    100509109    1003016875      1892953  ticks
       
  Memory: 31.858673095703125 GB (18060.64453125 MB free)
  Uptime: 1.529048e6 sec
  Load Avg:  0.0  0.0  0.0
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-6.0.1 (ORCJIT, skylake)
```