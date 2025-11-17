First, we create a containerized conda environment: 
```bash
module load hpc-container-wrapper
conda-containerize new --prefix env_neutral_sim conda_env.yaml
```
To load the environment, execute: 
```bash
export PATH="$PWD/env_neutral_sim/bin:$PATH" 
```

To execute a single simulation: 
```
time ./pbs_sim.py config.yaml 1111 test.npy
```
