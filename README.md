# BLISlab: A Sandbox for Optimizing GEMM

BLISlab is a small, self-contained “lab” codebase for learning and experimenting with GEMM (matrix multiply)
optimization techniques. The repository is organized as a set of steps (`step1` … `step5`) where each step is
its own buildable project with:

- a reference GEMM implementation
- one or more optimized implementations (micro-kernels, assembly kernels, etc.)
- a simple benchmark/test driver that reports throughput (GFLOPS)

If you are looking for the lab handout, see `tutorial.pdf`.

## Repo Layout

- `step1`–`step4`: DGEMM (double-precision) variants and progressively more optimized kernels.
- `step5/single`: SGEMM (single-precision) kernels for x86 (uses `-mavx`).
- `step5/arm`: SGEMM kernels for ARMv7 (uses NEON flags like `-mfpu=neon -march=armv7-a`).
- `misc/results`: helper scripts and `.m` files for collecting results and plotting (MATLAB/Octave-style output).
- `misc/examples`: small standalone examples (not part of the step builds).
- `Untitled Diagram.drawio`: a diagram file (if you use draw.io).

## Quickstart (Run Any Step)

Each step is designed to be built and run from within its directory.

```bash
cd step1            # or: step2 step3 step4 step5/single step5/arm
source sourceme.sh  # sets BLISlab build/runtime environment
make
cd test
./run_bl_dgemm.sh   # step1–step4
# ./run_bl_sgemm.sh # step5
```

The `run_bl_*` scripts print MATLAB/Octave-friendly output (a matrix with columns like `m n k MY_GFLOPS REF_GFLOPS`).
You can redirect that output to a file and post-process it with scripts under `misc/results`.

## Configuration Notes

Most configuration is done via environment variables exported in each step’s `sourceme.sh`:

- `BLISLAB_USE_INTEL`: `true` uses Intel toolchain settings (`icc/icpc`); `false` uses GNU (`gcc/g++`).
- `BLISLAB_USE_BLAS`: when `true`, the reference path links against an external BLAS (often OpenBLAS/BLIS).
- `COMPILER_OPT_LEVEL`: `O0`/`O1`/`O2`/`O3`.
- `BLAS_DIR`: path to your BLAS install (only used when `BLISLAB_USE_BLAS=true` for GNU builds).
- `OMP_NUM_THREADS`, `BLISLAB_IC_NT`, `KMP_AFFINITY`: threading controls used by the benchmark scripts.

Important: the default `BLAS_DIR` values in `sourceme.sh` are machine-specific placeholders (e.g. `/u/jianyu/...`,
`/home/ubuntu/...`). If you enable `BLISLAB_USE_BLAS=true`, update `BLAS_DIR` to match your system.

## Build Outputs

After `make`, you should see:

- `lib/libblislab.a` and `lib/libblislab.so`
- test executables under `test/` (for example `test_bl_dgemm.x` or `test_bl_sgemm.x`)

## Common Gotchas

- x86 steps use `-mavx`; run them on an AVX-capable CPU.
- `step5/arm` uses ARMv7/NEON compiler flags; it is intended for ARM targets (or cross-compilation with a suitable toolchain).
- Some Intel-oriented scripts set `DYLD_LIBRARY_PATH` (macOS) or assume Intel runtime libraries; adjust as needed for your environment.
