# ActiveSwingAsymmetricWalkingIncreaseTrunkKinematicVariability.jl

This repository contains the Julia code, Jupyter notebook, and dataset used in the study "Active Arm Swing and Asymmetric Walking Leads to Increased Variability in Trunk Kinematics in Young Adults" by Mezher et al..

## Instructions

To run this analysis on your computer, both Julia and Jupyter Notebook must be available. A version of Julia appropriate for your OS can be downloaded from [the Julia website](https://julialang.org/downloads/), and Jupyter can be installed from within Julia (in the REPL) with `] add IJulia`; alternate instructions for installing Jupyter can be found on the [IJulia github](https://github.com/JuliaLang/IJulia.jl) or the [Jupyter homepage](https://jupyter.org/install) (not recommended).

From within the main repository directory, start Julia and then start Jupyter in the Julia REPL (`using IJulia; notebook(;dir=pwd())`, or if using a system Jupyter installation, start Jupyter from your favorite available shell (e.g. Powershell on Windows, bash on any \*nix variant, etc.). In Jupyter, open the `Analysis.ipynb` notebook. Running all cells will reproduce the results and sole figure of the above mentioned paper.

## Description of data

The `data` directory contains the demographics and raw data. Each `.mat` file contains gait events, including\*:

- `LFO`/`RFO` (Left/right foot lift-off)
- `LFC`/`RFC` (Left/right foot contact)

Data signals found in each `.mat` file are as follows:

- `FP1`/`FP2` (left and right force plates, respectively)
- `LFootPos`/`RFootPos` (Left/right foot COM position)
- `LFootLinVel`/`RFootLinVel` (Left/right foot linear velocity)
- `LFootVelwrtPelvis`/`RFootVelwrtPelvis` (Left/right foot linear velocity, relative to the pelvis)
- `LFootAngle`/`RFootAngle` (Left/right foot angle, relative to the lab reference frame)
- `LFootAngVel`/`RFootAngVel` (Left/right foot angular velocity, relative to the lab reference frame)
- `LShoulder`/`RShoulder` (Left/right shoulder angle, extracted in the order SAGITTAL-FRONTAL-CORONAL)
- `LHip`/`RHip` (Left/right hip angle, extracted in the order SAGITTAL-FRONTAL-CORONAL)
- `TrunkPos` (Trunk position)
- `TrunkAngle` (Trunk angle, relative to the lab reference frame)
- `TrunkLinVel`/`TrunkAngVel` (Trunk linear/angular velocity, relative to the lab reference frame)
- `COG` (Whole-body COM/COG)
- `COG_Velocity` (Whole-body COM/COG linear velocity)
- `ModelAngMmntm` (Whole-body angular momentum--WBAM)

Units for all data are standard Visual3D units.

Output variables in the `results.csv` file:

- `left_steplength`/`right_steplength` (Left/right average step length)
- `SD_left_steplength`/`SD_right_steplength` (Left/right step length standard deviation)
- `left_steptime`/`right_steptime` (Left/right average step time)
- `SD_left_steptime`/`SD_right_steptime` (Left/right step time standard deviation)
- `stepwidth`/`SD_stepwidth` (Average step width and step width standard deviation)
- `lvmean`/`avmean` (Mean of the average linear and angular velocity for each stride)
- `lvstd`/`avstd`  (Mean of the standard deviation of linear and angular velocity for each stride)
- `lvmax`/`avmax`  (Mean of the maximum (peak) linear and angular velocity for each stride)
- `WBAM_mean`/`WBAM_std`  (Mean of the average and standard deviation of WBAM for each stride)

\*Variables named `TRST`, `TREN`, `SLST`, `SLEN` are present, but empty, and should be ignored.

