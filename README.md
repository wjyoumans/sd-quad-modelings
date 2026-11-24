# sd-quad-modelings

Code for the experiments of *An Improved Quadratic Modeling for Exact
Syndrome Decoding*.

## Requirements

Magma (the paper used V2.28-27) and Python 3.8 or later.

## Usage

```
magma -b check.m                                        # a few seconds
magma -b inst:=instances/SD_32_16_0 model:=CSA run.m    # a single job
python3 run.py --nmax 24                                # a quick run, under a minute
python3 run.py                                          # everything, about 8 hours
```

`run.py` writes the log of every job to `logs/` and the results to
`results/A.csv` and `results/B.csv`. Its options are `--cap` (time limit per
job, default 300 s), `--mem` (memory limit, default 40 GiB), `--cpu` (the
core to pin the jobs to on Linux, default 0) and `--nmax` (largest n, default
120); `python3 run.py A` runs only regime A. Jobs that already have a log are
skipped, so an interrupted run can be resumed.

## Experiments

Instances are planted: H is a uniformly random full-rank (n−k) x n matrix over
F_2, e a uniformly random vector of weight t = t_GV(n, k), the largest
t <= n/2 with binom(n, t) <= 2^(n−k), and s = He. An instance depends only on
(n, k, seed).

- Regime A: k = n/2.
- Regime B: k = floor(n/4).

Both use n = 16, 18, 20, ... and seeds 0 to 4. Each job has a limit of 300 s
and 40 GiB. A modeling's median time at n is infinite once the first seed, or
three of the five, hit a limit; its other seeds at that n are then skipped. A
modeling stops after two consecutive n with an infinite median.

## Results

Each row of `results/<regime>.csv` (and `paper/<regime>.csv`) is one job.
The paper reports, for each n and modeling, the median `gb_cpu_s` over the
seeds, the range of `soldeg` over the seeds, and the largest matrix of the finished
seed with the median number of entries (the lower median for an even number of seeds).

## License

MIT, see `LICENSE`.

## AI disclosure
Claude Opus 5.5 was used to fix issues in the human written code for the CCMMP2 and CSA modelings and help extend support to the RST and LMT modelings.
Claude was also used for general correctness checks, cleanup, and for preparing the final version of this repository, including generating `check.m`, `generate.py`, `instance.m`, `run.m`, and `run.py`.
All code was reviewed, edited, tested, and verified by the authors.
