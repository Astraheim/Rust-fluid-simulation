// Window parameters
pub const ADAPT_TO_WINDOW: bool = true;
pub const WINDOW_HEIGHT: usize = if ADAPT_TO_WINDOW { (N as usize +2) * DY as usize } else { 720 };
pub const WINDOW_WIDTH: usize = if ADAPT_TO_WINDOW { (N as usize +2) * DX as usize } else { 720 };



// Simulation parameters
pub const SIM_STEPS: usize = 5000;
pub const CIP_CSL4 : bool = false; // (For now have to keep it on false) Use CIP-CSL4 method for density advection
pub const DENS_ADV_FAC: f32 = 0.1; // Use Morton encoding for grid cells
pub const VEL_STEP: &str = "2"; // "1" for vel_step, "2" for vel2_step, "cip_csl4" for vel_step_cip_csl4
pub const PROJECT : &str = "1"; // "1" for project, "2" for project2



// Flow parameters
pub const AIR_FLOW: bool = false; // Simulate air flow
pub const FLOW_DIRECTION: &str = "right"; //"left","right","up","down"  // Direction of the flow
pub const FLOW_SPACE: usize = 1; // Space between two rows of flow
pub const FLOW_DENSITY: f32 = 25.0; // Density of the flow
pub const FLOW_VELOCITY: f32 = if AIR_FLOW {0.1} else { 0.0 }; // Velocity of the flow
pub const DRAW_VELOCITY_VECTORS: bool = true; // Draw velocity vectors
pub const VECTOR_SIZE_FACTOR : f32 = 8.0 ; // Factor for the size of the velocity vector
pub const INFLOW_VELOCITY: f32 = FLOW_VELOCITY; // Always force inflow velocity to be equal to the flow velocity
pub const PAINT_VORTICITY: bool = false; // Paint vorticity


// Grid parameters
pub const EXT_BORDER: bool = false;  // have to keep it true for now
pub const N: f32 = 254.0; // There are special conditions for the size of the grid : when using Morton encoding, the grid size must be a power of 2, then subtract 2. // IT MUST BE INTEGER \\
pub const DX: f32 = 3.0; // Size of a cell (in pixel), horizontal. // IT MUST BE INTEGER \\
pub const DY: f32 = 3.0; // Size of a cell (in pixel), vertical. // IT MUST BE INTEGER \\
pub const SIZE: f32 = (N + 2.0) * (N + 2.0); // IT WILL BE INTEGER \\




// Physical parameters
pub const DT: f32 = 1.0/60.0;



// Fluid parameters
// todo : solve viscosity problems
pub const VISCOSITY: f32 = 0.00; // A bit of viscosity for Karman vortex


