// Last update to rendu_code_tex: 2025-05-02
// last modif: 2025-05-09

/*
FR :
    Ce fichier centralise tous les paramètres de la simulation : affichage, grille, fluides, méthodes numériques, etc.
    Il permet de configurer facilement la fenêtre (taille, adaptation à la grille), le comportement du fluide (viscosité, flux, vortex),
    le choix des algorithmes (projection, advection), ainsi que les sources de densité et de vitesse.

    L’objectif est de rendre la simulation modulaire, expérimentable, et facilement ajustable sans toucher au cœur du moteur.


ENG :
    This file centralizes all parameters for the simulation : display, grid, fluid properties, numerical methods, etc.
    It allows configuring the window (size, grid adaptation), fluid behavior (viscosity, flows, vortex),
    algorithm choices (projection, advection), and sources of density and velocity.

    The goal is to make the simulation modular, easy to experiment with, and adjustable without touching core engine logic.
*/
use crate::grid2::Vector22;

// Window parameters
pub const ADAPT_TO_WINDOW: bool = true;  // Adapt the window size to the grid size
pub const WINDOW_HEIGHT: usize = if ADAPT_TO_WINDOW == true { (N as usize +2) * DY as usize } else { 720 };  // Height of the window
pub const WINDOW_WIDTH: usize = if ADAPT_TO_WINDOW == true { (N as usize +2) * DX as usize } else { 720 };   // Width of the window



// Simulation parameters
pub const GRID : &str = "1"; // "1" for grid, "2" for grid2
pub const SIM_STEPS: usize = 5000; // Potentially the number of simulation steps
pub const PRINT_FORCES: bool = false; // Print forces
pub const CIP_CSL4 : bool = false; // (For now have to keep it on false) Use CIP-CSL4 method for density advection
pub const DENS_ADV_FAC: f32 = 0.1; // Factor for density advection
pub const VEL_STEP: &str = "2"; // "1" for vel_step, "2" for vel2_step, "cip_csl4" for vel_step_cip_csl4
pub const PROJECT : &str = "1"; // "1" for project, "2" for project2
pub const KARMAN_VORTEX: bool = false; // Use Karman vortex (at least tries to)
pub const CENTER_SOURCE: bool = false; // Use a circular source in the center of the grid
pub const CENTER_SOURCE_TYPE: bool = false; // "false" for a circular output, "true" for a radial output (like a fan)
pub const CENTER_SOURCE_RADIUS: f32 = N/40.0; // Radius of the center source
pub const CENTER_SOURCE_DENSITY: f32 = 0.15; // Density of the center source
pub const CENTER_SOURCE_VELOCITY: f32 = 0.5; // Velocity of the center source




// Flow parameters
pub const AIR_FLOW: bool = true; // Simulate air flow
pub const FLOW_DIRECTION: &str = "right"; //"left","right","up","down" // Direction of the flow (for now only right is properly implemented)
pub const FLOW_SPACE: usize = 4; // Space between two rows of flow
pub const FLOW_DENSITY: f32 = 25.0; // Density of the flow
pub const FLOW_VELOCITY: f32 = if AIR_FLOW == true {1.0} else { 0.0 }; // Velocity of the flow
pub const DRAW_VELOCITY_VECTORS: bool = false; // Draw velocity vectors
pub const VECTOR_SIZE_FACTOR : f32 = 8.0 ; // Factor for the size of the velocity vector
pub const INFLOW_VELOCITY: f32 = FLOW_VELOCITY; // Always force inflow velocity to be equal to the flow velocity
pub const PAINT_VORTICITY: bool = false; // Paint vorticity


// Grid parameters
pub const EXT_BORDER: bool = false;  // enable or disable external borders
pub const N: f32 = 510.0; // There are special conditions for the size of the grid: when using Morton encoding, the grid size must be a power of 2, then subtract 2. // IT MUST BE AN INTEGER \\
pub const DX: f32 = 2.0; // Size of a cell (in pixel), horizontal. // IT MUST BE AN INTEGER \\
pub const DY: f32 = 2.0; // Size of a cell (in pixel), vertical. // IT MUST BE AN INTEGER \\
pub const SIZE: f32 = (N + 2.0) * (N + 2.0); // IT WILL BE AN INTEGER \\




// Physical parameters
pub const DT: f32 = 1.0/600.0; // Time step (in seconds)
pub fn gravity() -> Vector22 {
    let gravity: Vector22 = Vector22::new(0.0, 0.0); // Gravity (in m/s^)
    gravity
}

// Fluid parameters
pub const VISCOSITY: f32 = 0.000; // A bit of viscosity for Karman vortex


// Log parameters
pub const LOG: bool = false; // Log simulation data



