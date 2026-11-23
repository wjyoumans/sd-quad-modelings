#!/usr/bin/env python3
"""Run the experiments and collect the results into results/<regime>.csv.

Usage: run.py [A] [B] [--cap 300] [--mem 40] [--cpu 0] [--nmax 120]

Regimes: planted instances with t = t_GV(n, k) (see generate.py), k = n/2 in
regime A and k = floor(n/4) in regime B, n = 16, 18, ..., nmax, and five
instances (seeds 0-4) per n. Jobs run one at a time, each in its own Magma
process (run.m), pinned to the CPU --cpu where the system supports it (Linux).
The output of each job goes to logs/<regime>/<model>_<n>_<seed>.log, with the
line KILLED appended if the job hit the time limit. A job whose log exists is
not run again, so delete the logs of the killed jobs to rerun them with a
larger --cap.

Stopping rule: a job hits a limit if its log has no RESULT line. A modeling
has an infinite median time at n once the first seed, or three of the five
seeds, hit a limit; its other seeds at that n are then skipped. A modeling
stops after two consecutive n with an infinite median.

Columns of results/<regime>.csv, one row per job:
  regime, n, k, t, seed, model
  status     ok, timeout (hit the time limit), memout (hit the memory limit)
             or error
  gb_cpu_s   CPU time of GroebnerBasis in seconds (only if ok)
  gb_real_s  real time of GroebnerBasis in seconds (only if ok)
  soldeg     solving degree: the largest step degree of Magma's F4, or the
             largest degree of the input polynomials if F4 ran no steps;
             for a job that hit a limit, the largest step degree it reached
  rows, cols, nnz, density
             the F4 matrix with the most entries (rows times columns),
             including those of the initial interreduction: its size, its
             number of nonzero entries, and nnz / (rows cols)
  nsols      number of solutions found (only if ok)
  ok         whether they all solve the instance and include the planted
             solution (only if ok)
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent
REGIMES = ["A", "B"]
MODELS = ["CCMMP2", "RST", "RCI", "CSA", "CHAM", "QFCHAM"]
SEEDS = range(5)
COLUMNS = ["regime", "n", "k", "t", "seed", "model", "status", "gb_cpu_s", "gb_real_s",
           "soldeg", "rows", "cols", "nnz", "density", "nsols", "ok"]

STEP_RE = re.compile(r"step degree: (\d+)")
MATRIX_RE = re.compile(r"= (\d+) rows / (\d+) columns")
DENSITY_RE = re.compile(r"^Density: .* total: (\d+)")
FIELDS_RE = re.compile(r"(\w+)=(\S+)")


# Running the jobs

def finished(log):
    """Whether the job of the log finished, i.e. printed its RESULT line."""
    return any(line.startswith("RESULT ") for line in log.read_text(errors="replace").splitlines())


def run_job(regime, model, n, seed, args):
    """Runs the job unless it already has a log, and returns whether it finished."""
    k = n // 2 if regime == "A" else n // 4
    inst = f"instances/SD_{n}_{k}_{seed}"
    if not (REPO / inst).exists():
        subprocess.run([sys.executable, "generate.py", str(n), str(k), str(seed)],
                       cwd=REPO, check=True)
    log = REPO / "logs" / regime / f"{model}_{n}_{seed}.log"
    if not log.exists():
        print(regime, model, f"n={n}", f"seed={seed}", flush=True)
        cmd = ["magma", "-b", f"inst:={inst}", f"model:={model}", f"mem:={args.mem}", "run.m"]
        try:
            with log.open("w") as f:
                subprocess.run(cmd, cwd=REPO, stdout=f, stderr=subprocess.STDOUT,
                               timeout=args.cap)
        except subprocess.TimeoutExpired:
            with log.open("a") as f:
                f.write("KILLED\n")
        except BaseException:
            # Interrupted (e.g. Ctrl-C): drop the partial log so the job runs again.
            log.unlink()
            raise
    return finished(log)


def run_regime(regime, args):
    (REPO / "logs" / regime).mkdir(parents=True, exist_ok=True)
    for model in MODELS:
        streak = 0  # consecutive n with an infinite median
        for n in range(16, args.nmax + 1, 2):
            limits = 0
            infinite = False
            for seed in SEEDS:
                if not run_job(regime, model, n, seed, args):
                    limits += 1
                    if seed == 0 or limits == 3:
                        infinite = True
                        break
            streak = streak + 1 if infinite else 0
            if streak == 2:
                break


# Collecting the results

def largest_matrix(matrices):
    """The (rows, cols, nnz) with the most entries, or None."""
    matrices = [m for m in matrices if m[2] is not None]
    return max(matrices, key=lambda m: m[0] * m[1]) if matrices else None


def row(regime, job, result, status, steps, maxdeg, matrix):
    """One CSV row. job: n, k, t, seed, model; result: the RESULT fields or
    None; steps: the F4 step degrees; matrix: (rows, cols, nnz) or None."""
    if steps:
        soldeg = max(steps)
    else:
        soldeg = maxdeg if status == "ok" else ""
    r = dict(regime=regime, **job, status=status, soldeg=soldeg)
    if matrix:
        rows, cols, nnz = matrix
        r.update(rows=rows, cols=cols, nnz=nnz, density=f"{nnz / (rows * cols):.4g}")
    if result:
        r.update(gb_cpu_s=f"{float(result['gb_cpu']):.3f}",
                 gb_real_s=f"{float(result['gb_real']):.3f}",
                 nsols=result["nsols"], ok=result["ok"])
    return r


def parse_log(regime, path):
    model, n, seed = path.stem.rsplit("_", 2)
    job = dict(n=int(n), k="", t="", seed=int(seed), model=model)
    result, maxdeg, steps, matrices = None, None, [], []
    text = path.read_text(errors="replace")
    for line in text.splitlines():
        if (m := STEP_RE.search(line)):
            steps.append(int(m.group(1)))
        elif (m := MATRIX_RE.search(line)):
            matrices.append([int(m.group(1)), int(m.group(2)), None])
        elif (m := DENSITY_RE.search(line)) and matrices and matrices[-1][2] is None:
            matrices[-1][2] = int(m.group(1))
        elif line.startswith("JOB "):
            fields = dict(FIELDS_RE.findall(line))
            job.update(k=int(fields["k"]), t=int(fields["t"]))
            maxdeg = int(fields["maxdeg"])
        elif line.startswith("RESULT "):
            result = dict(FIELDS_RE.findall(line))
    if result:
        status = "ok"
    elif "memory limit" in text.lower():
        status = "memout"
    elif text.rstrip().endswith("KILLED"):
        status = "timeout"
    else:
        status = "error"
    return row(regime, job, result, status, steps, maxdeg, largest_matrix(matrices))


def write_csv(path, rows):
    rows = sorted(rows, key=lambda r: (r["n"], MODELS.index(r["model"]), r["seed"]))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, COLUMNS)
        w.writeheader()
        w.writerows(rows)


def collect():
    (REPO / "results").mkdir(exist_ok=True)
    for logdir in sorted((REPO / "logs").glob("*")):
        rows = [parse_log(logdir.name, p) for p in logdir.glob("*.log")]
        write_csv(REPO / "results" / f"{logdir.name}.csv", rows)
        for r in rows:
            if r["status"] == "error" or r.get("ok") == "false":
                print(f"check {logdir.name} {r['model']} n={r['n']} seed={r['seed']}: "
                      f"status {r['status']}, ok {r.get('ok')}")
        print(f"{len(rows)} jobs -> results/{logdir.name}.csv")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("regimes", nargs="*", help="A, B or both (default: both)")
    ap.add_argument("--cap", type=float, default=300, help="time limit per job, seconds")
    ap.add_argument("--mem", type=int, default=40, help="memory limit per job, GiB")
    ap.add_argument("--cpu", type=int, default=0, help="the CPU to pin the jobs to (Linux)")
    ap.add_argument("--nmax", type=int, default=120, help="largest n")
    args = ap.parse_args()
    regimes = args.regimes or REGIMES
    for regime in regimes:
        if regime not in REGIMES:
            ap.error(f"unknown regime {regime}")

    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {args.cpu})  # inherited by the Magma processes
    for regime in regimes:
        run_regime(regime, args)
    collect()


if __name__ == "__main__":
    main()
