**Actuator Design code**

Set of codes to design, analyse and optimize a 5-Panel hinged structure with SMA actuation.

**Abaqus scripts**

Code to generate and run model for a 5-Panel hinged structure with elastomer.
*5Panel* - Runs abaqus model of 5 Panel hinged assembly with elastomer.
*TorqueTube* -  Runs abaqus model of SMA torque tube analysiss

**Python Kinematics Solver**

Code to calculate kinemtics of 5-Panel hinged structure using Pyslvs and find design variables that are closest to target curve.
*optimization_script* - Generates inputs (Design Variables) and calls 5Panel_Kinematics.
*5Panel_Kinematics* - Calculates geometry parameters and calls Pyslvs; saves output in txt file.
*PP_DOE* - Post processes results from tx file and plots best output from DOE.

**Library Dependencies**

*5Panel*

- Post_P_Script
- math
- os

*optimization_script/5Panel_Kinematics*

- Pyslvs
- subprocess
- math
- numpy
- scipy
- matplotlib
