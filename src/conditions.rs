use std::arch::aarch64::float32x2_t;

// Window parameters
pub const ADAPT_TO_WINDOW: bool = true;
pub const WINDOW_HEIGHT: usize = if ADAPT_TO_WINDOW { (N as usize +2) * DY as usize } else { 720 };
pub const WINDOW_WIDTH: usize = if ADAPT_TO_WINDOW { (N as usize +2) * DX as usize } else { 720 };



// Simulation parameters
pub const CIP_CSL4 : bool = false; // (For now have to keep it on false) Use CIP-CSL4 method for density advection
pub const DENS_ADV_FAC: f32 = 0.1;



// Flow parameters
pub const AIR_FLOW: bool = false; // Simulate air flow
pub const FLOW_DIRECTION: &str = "right"; //"left","right","up","down"  // Direction of the flow
pub const FLOW_SPACE: usize = 3; // Space between two rows of flow
pub const FLOW_DENSITY: f32 = 20.0; // Density of the flow
pub const FLOW_VELOCITY: f32 = 0.7; // Velocity of the flow  (Interesting values : 1.0, 0.25)



// Grid parameters
pub const EXT_BORDER: bool = true;  // have to keep it true for now
pub const N: f32 = 254.0; // There are special conditions for the size of the grid : when using Morton encoding, the grid size must be a power of 2, then subtract 2. // IT MUST BE INTEGER \\
pub const DX: f32 = 2.0; // Size of a cell (in pixel), horizontal. // IT MUST BE INTEGER \\
pub const DY: f32 = 2.0; // Size of a cell (in pixel), vertical. // IT MUST BE INTEGER \\
pub const SIZE: f32 = (N + 2.0) * (N + 2.0); // IT WILL BE INTEGER \\




// Physical parameters
pub const DT: f32 = 1.0/60.0;



// Fluid parameters
// todo : solve viscosity problems
pub const VISCOSITY: f32 = 0.0;

