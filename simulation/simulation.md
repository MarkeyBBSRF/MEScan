## MEScan paper simulation

### 1. Generate simulation data

The simulation data were generated using functions in `generate_data` folder.
  1. `SimulationFunction.R` simulation data not including highly mutated gene.
  2. `SimulationFunctionAddTP53.R` simulation data with one highly mutated gene.

### 2. Tests using simulation data.

We tested different algorithms using the simulation data created in step 1. The source code for testing each algorithm can be found in their individual folders.
  * MEScan
  * Dendrix
  * CoMEt
  * WExT
  * TiMEx
  * MESGA

### 3. Plot figures

Paper figures were created using functions in `plot_simulation`.

