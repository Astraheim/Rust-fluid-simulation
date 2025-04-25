use crate::conditions::*;
use crate::grid::*;
use minifb::{Window, WindowOptions};
use std::time::Instant;
use crate::pressure_computation::*;

// Vorticity calculation for colourisation (safe bounds)
fn vorticity(grid: &Grid, i: usize, j: usize) -> f32 {
    let dx = DX as f32;
    let dy = DY as f32;

    let im = i.saturating_sub(1).min(N as usize);
    let ip = (i + 1).min(N as usize);
    let jm = j.saturating_sub(1).min(N as usize);
    let jp = (j + 1).min(N as usize);

    let idx_up    = grid.to_index(i, jp);
    let idx_down  = grid.to_index(i, jm);
    let idx_left  = grid.to_index(im, j);
    let idx_right = grid.to_index(ip, j);

    let du_dy = (grid.cells[idx_up].velocity.x - grid.cells[idx_down].velocity.x) / (2.0 * dy);
    let dv_dx = (grid.cells[idx_right].velocity.y - grid.cells[idx_left].velocity.y) / (2.0 * dx);

    dv_dx - du_dy
}

// Draw a line using Bresenham's algorithm
fn draw_line(buffer: &mut [u32], width: usize, x0: usize, y0: usize, x1: usize, y1: usize, color: u32) {
    let dx = (x1 as isize - x0 as isize).abs();
    let dy = -(y1 as isize - y0 as isize).abs();
    let mut err = dx + dy;
    let mut x = x0 as isize;
    let mut y = y0 as isize;
    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };

    loop {
        if x >= 0 && x < width as isize && y >= 0 && y < width as isize {
            buffer[(y as usize) * width + (x as usize)] = color;
        }
        if x == x1 as isize && y == y1 as isize { break; }
        let e2 = 2 * err;
        if e2 >= dy {
            err += dy;
            x += sx;
        }
        if e2 <= dx {
            err += dx;
            y += sy;
        }
    }
}

// Map density to color when not painting vorticity
fn density_color(density: f32) -> u32 {
    if density <= 20.0 {
        let intensity = (255.0 * (1.0 - density / 20.0)).clamp(0.0, 255.0) as u32;
        (intensity << 16) | (intensity << 8) | 0xFF
    } else {
        let excess = (density - 20.0).clamp(0.0, 20.0);
        let red_intensity = (255.0 * (excess / 20.0)).clamp(0.0, 255.0) as u32;
        (red_intensity << 16) | 0xFF
    }
}

// Launch the simulation and rendering
pub fn run_simulation(grid: &mut Grid, mut step: i32) {
    let start = Instant::now();
    let width = WINDOW_WIDTH;
    let height = WINDOW_HEIGHT;
    let mut buffer: Vec<u32> = vec![0; width * height];

    let mut window = Window::new(
        "Grid Simulation",
        width,
        height,
        WindowOptions::default(),
    )
        .unwrap_or_else(|e| panic!("{}", e));

    println!("Starting simulation");

    let n_max = N as usize + 1;
    while window.is_open() {
        // Clear buffer
        buffer.fill(0xFFFFFFFF);

        // Draw cells
        for j in 1..n_max {
            for i in 1..n_max {
                let idx = grid.to_index(i, j);
                let x0 = i * DX as usize;
                let y0 = j * DY as usize;

                // Draw walls in black
                if grid.cells[idx].wall {
                    for dy in 0..DY as usize {
                        for dx in 0..DX as usize {
                            buffer[(y0 + dy) * width + (x0 + dx)] = 0x000000;
                        }
                    }
                    continue;
                }
                // Skip empty density
                if grid.cells[idx].density == 0.0 {
                    continue;
                }

                // Choose colour based on mode
                let colour = if PAINT_VORTICITY {
                    let vort = vorticity(grid, i, j);
                    let iv = ((vort.abs().min(5.0)) * 25.0) as u32;
                    if vort > 0.0 {
                        (iv << 8) | 0x0000FF
                    } else {
                        (iv << 16) | 0xFF0000
                    }
                } else {
                    density_color(grid.cells[idx].density)
                };

                // Fill cell block
                for dy in 0..DY as usize {
                    for dx in 0..DX as usize {
                        buffer[(y0 + dy) * width + (x0 + dx)] = colour;
                    }
                }

                // Draw velocity vector
                if DRAW_VELOCITY_VECTORS {
                    let vx = grid.cells[idx].velocity.x;
                    let vy = grid.cells[idx].velocity.y;
                    let cx = x0 + (DX as usize) / 2;
                    let cy = y0 + (DY as usize) / 2;
                    let mag = (vx*vx + vy*vy).sqrt().max(1e-5);
                    let scale = VECTOR_SIZE_FACTOR / mag;
                    let ex = (cx as f32 + vx * scale) as usize;
                    let ey = (cy as f32 + vy * scale) as usize;
                    draw_line(&mut buffer, width, cx, cy, ex, ey, 0x000000);
                }
            }
        }

        // Update window
        window
            .update_with_buffer(&buffer, width, height)
            .expect("Error while updating the window");

        step += 1;

        // Flow injection or custom source
        if AIR_FLOW {
            let hole_pos: Vec<usize> = (1..=N as usize)
                .filter(|&x| x % FLOW_SPACE == 0)
                .collect();
            grid.initialize_wind_tunnel(FLOW_DENSITY, &hole_pos);
        } else if KARMAN_VORTEX {
            for i in 0..=50 {
                Grid::cell_init(grid, 5, 127+i-25, 0.3, 0.0, 20.0);

                // for karman vortex
                Grid::velocity_init(grid, 10, 127+i-25, 0.9, 0.0);
            }
        } else if CENTER_SOURCE {
            grid.center_source(CENTER_SOURCE_RADIUS, CENTER_SOURCE_DENSITY, CENTER_SOURCE_VELOCITY, CENTER_SOURCE_TYPE);
        }


        // Perform simulation step
        if VEL_STEP == "1" {
            grid.vel_step();
        } else if VEL_STEP == "2" {
            grid.vel2_step();
        } else if VEL_STEP == "cip_csl4" {
            grid.vel_step_cip_csl4()
        } else {
            panic!("Invalid VEL_STEP value");
        }


        if PRINT_FORCES {
            grid.print_object_forces()
        }

        // Optionally print timing
        if step == 100  {
            let elapsed = start.elapsed();
            println!("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n Temps pour 100 vel_step: {:?} \n $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ", elapsed);
        }
        // Optionally print timing
        if step == SIM_STEPS as i32 {
            let elapsed = start.elapsed();
            println!("§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ \n Completed {} steps in {:?}\n§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§", SIM_STEPS, elapsed);
        }
    }
}

