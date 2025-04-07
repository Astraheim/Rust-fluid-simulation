// Window parameters
pub const WINDOW_HEIGHT: usize = 612;
pub const WINDOW_WIDTH: usize = 612;


// Simulation parameters
pub const _SIM_STEPS: usize = 100;


// Grid parameters
pub const EXT_BORDER: bool = true;
pub const N: usize = 300; // Number of inside cells (non walls ones), if EXT_BORDER is enabled, you must subtract 2 to N
pub const SIZE: usize = (N + 2) * (N + 2);




// Physical parameters
pub const DT: f32 = 1.0/60.0;


// Fluid parameters
// todo : solve viscosity problems
pub const VISCOSITY: f32 = 0.0;

