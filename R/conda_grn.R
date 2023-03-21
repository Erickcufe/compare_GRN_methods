library(reticulate)

py_config()
use_python("/Users/erickCuevas/opt/miniconda3/envs/grn_env/bin/python", required = TRUE)

# Activate the Conda environment
use_condaenv("grn_env", required = TRUE)

# Load necessary Python libraries
scanpy <- import("scanpy")
pyscenic <- import("pyscenic")
arboreto <- import("arboreto")
igraph <- import("igraph")
matplotlib <- import("matplotlib")
