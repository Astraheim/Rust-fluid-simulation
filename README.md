# Rust-fluid-simulation
A fluid simulation with rust for my MP2I TIPE school project. 

The goal is to code a real-time fluid simulation for a school project.
I've chosen an Eulerian simulation to allow for future expansion into an "infinite" world.
For solving advection of densities I'm using a shitty solution for now but planing on adding an Cubic Interpolated Propagation - Conservative Semi-Lagrangian (CIP-CSL) method.
Here, I plan to use WGPU to optimize my code as much as possible.
Rayon is also being used to parallelize the calculations.
By now, I'm using minifb to visualise my simulation, but a said, I'll switch to wgpu.
