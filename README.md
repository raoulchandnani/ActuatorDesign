# Actuator Design

Set of codes to design, analyse and optimize a 5-panel hinged structure with thermal/ SMA actuation.

**Abaqus scripts**

Code to generate and run model for a 5-panel hinged structure with conformal surface.

*5Panel_Thermal* - Runs abaqus model of 5 panel hinged assembly using thermal actuation

*5Panel_SMA* - Runs abaqus model of 5 panel hinged assembly using SMA actuation

*5Panel_Thermal_OML* - Runs abaqus model of 5 panel hinged assembly using thermal actuation with conformal surface as an input

*5Panel_Thermal_MultipleDevices* - Runs abaqus model of multiple device 5 Panel hinged assembly using thermal actuation with conformal surface as an input

*TorqueTube* -  Runs abaqus model of SMA torque tube analysis

**Data processing scripts**

Code to post-process data 

*Clustering* - Example code to cluster aero data into regions

*Post_Process_DOE* - Code to post-process DOE results ( for scripts in DOE & optimization scripts)

*Post_Process_3D* - Code to post-process multiple device results

**DOE & optimization scripts**

Code for DOE and optimization runs. Each folder contains an abaqus script and an 'optimization_script'. The 'optimization_script' has variables that control the design space and is to be run in a python shell. This script will call the python script and generate required results.

*3 panel- DOE* - Runs DOE for 3 panel case with 1 thermally actuated torque tube

*4 panel- DOE* - Runs DOE for 4 panel case with 2 thermally actuated torque tubes

*5 panel- DOE/optimization* - Runs DOE/optimization for symmetric 5 panel case with 2 thermally actuated torque tubes

**Python Kinematics Solver**

Code to calculate kinemtics of 5-Panel hinged structure using Pyslvs and find design variables that are closest to target curve.

*optimization_script* - Generates inputs (Design Variables) and calls 5Panel_Kinematics.

*5Panel_Kinematics* - Calculates geometry parameters and calls Pyslvs; saves output in txt file.

*PP_DOE* - Post processes results from tx file and plots best output from DOE.

**Library Dependencies**

*Abaqus scripts*

- Post_P_Script
- math
- os
- Aero_Data_I (text file containing aero data)

*optimization_script/5Panel_Kinematics*

- Pyslvs ( setup availabe in 'Pyslvs-UI-master' )
- subprocess
- math
- numpy
- scipy
- matplotlib
