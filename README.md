# Cardiac Electrophysiology Simulation

## What This Code Does

This is a 2D cardiac tissue simulation that models Ca-V dynamics in the presence of triggered calcium waves. The code specifically models cells with a poorly developed t-tubule system, where calcium influx through L-type calcium channels (LCC) is largely confined to the cell periphery.

**Key Model Features:**

- **Two-zone cellular structure**: The cell is divided into a boundary zone (containing junctional RyR clusters near the membrane) and an interior zone (containing non-junctional clusters)

- **Five calcium compartments**: Tracks calcium concentrations in boundary cytosol, boundary SR, interior cytosol, junctional SR (JSR), and network SR (NSR)

- **Stochastic calcium sparks**: Models calcium release using probabilistic spark recruitment at both boundary and interior sites

- **Graded release mechanism**: Boundary spark recruitment depends on LCC current, open probability, and SR calcium load with Hill function sensitivity

- **Wave propagation**: Interior recruitment includes baseline leak, boundary-triggered waves (when spark fraction exceeds threshold), and regenerative wave propagation with dynamic gating

- **Spatial coupling**: 2D tissue with electrical diffusion between neighboring cells through gap junctions

The simulation runs on a user-defined grid of cells (default 10×10) and outputs voltage and calcium time series data from a pacing protocol with 5 stimulus beats.

## Files Required

- **`wave2.f`** - Main simulation program
- **`currentv1.f`** - Ion channel subroutines
- **`currentsca1.f`** - Calcium handling subroutines  
- **`binom.f`** - Random number generation for calcium sparks

## How to Run

### 1. Compile
```bash
ifx wave2.f currentv1.f currentsca1.f binom.f -o cardiac_sim
```

**Alternative compiler:**
```bash
gfortran wave2.f currentv1.f currentsca1.f binom.f -o cardiac_sim
```

### 2. Execute
```bash
./cardiac_sim
```

**Runtime:** ~1-2 minutes

## Output Files

The simulation generates three data files with time series data:
- **`vx.dat`** - Membrane voltage (mV) vs time (ms)
- **`cbx.dat`** - Boundary calcium concentration (μM) vs time (ms)
- **`cix.dat`** - Interior calcium concentration (μM) vs time (ms)

Each file contains two columns: `time  value`
