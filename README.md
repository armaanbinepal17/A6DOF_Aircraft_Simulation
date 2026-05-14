 6-DOF Aircraft Flight Dynamics Simulation

Overview
A six degree-of-freedom rigid body dynamics simulator 
built in MATLAB to model autonomous aircraft behaviour 
for a UAV competition. The simulation models all six 
axes of motion — three translational and three 
rotational — to predict flight characteristics and 
validate control system design before physical testing.

System Architecture
The simulator implements the full nonlinear equations 
of motion for a rigid aircraft in body-frame 
coordinates, integrated using MATLAB's ODE45 solver 
with tight relative and absolute tolerances 
(1e-6, 1e-9).

Core components:
- 6-DOF equations of motion in body frame
- Aerodynamic force and moment model with 
stability derivatives
- Euler angle kinematics for attitude propagation
- Rotation matrix transformation from body to 
inertial frame
- Visualisation suite with 8 subplot analysis panels

Technical Highlights
- Models full inertia tensor including cross-product 
term Ixz for realistic mass distribution
- Implements nonlinear aerodynamic coefficients 
including angle of attack dependent drag and 
lift, sideslip dependent side force, and full 
stability derivative moment model
- Angular momentum equations solved using inverse 
inertia tensor to correctly handle gyroscopic 
coupling terms
- 3D trajectory visualisation with energy 
conservation analysis to validate simulation 
integrity

Aircraft Parameters
- Mass: 18 kg
- Wingspan: 0.5 m
- Roll/Pitch/Yaw inertia: 0.4/0.6/1.2 kg·m²

Tools and Technologies
- MATLAB with ODE45 numerical integration
- SolidWorks for physical aircraft geometry 
(separate repo)

Results
Simulation successfully models stable flight 
from initial conditions of 100 m/s forward 
velocity at 1000 m altitude. Energy conservation 
analysis confirms simulation integrity across 
30 second flight envelope.

What I'd Do Differently
I would integrate the simulation directly with 
a Simulink control system model to enable 
hardware-in-the-loop testing, and add wind 
disturbance modelling to validate control 
robustness under realistic atmospheric conditions.
