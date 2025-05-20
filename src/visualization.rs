use crate::conditions::*;
use crate::grid::*;
use crate::grid2::*;
use minifb::{Window, WindowOptions};
use std::time::Instant;
use crate::conditions;
use crate::grid2::*;

/*
FR :
    Ce fichier contient la boucle principale de la simulation ainsi que les fonctions associées au rendu graphique, à l’interaction souris, 
    à la coloration des cellules et au tracé de vecteurs. 
    La fonction `run_simulation` gère le rendu de la grille à chaque pas de temps, l’injection de fluide, le tracé des murs à la souris, 
    et les différentes sources (flux d’air, vortex de Kármán, source centrale).
    Elle appelle aussi les étapes de simulation selon le mode choisi ('vel_step', `vel2_step', `cip_csl4').
    
    Fonctions principales :
    - `run_simulation` : boucle temps, dessin, interaction, injection, calcul
    - `vorticity` : calcule la vorticity locale pour la coloration
    - `density_color` : mappe une densité scalaire à une couleur RGB
    - `bresenham_line` / `draw_line_with_color` : algorithmes pour dessiner des lignes discrètes

    Interactions utilisateur :
    - Clic-gauche (maintenu) : dessine des murs à la souris sur la grille
    - Touche P : met la simulation en pause ou la redémarre

ENG :
    This file contains the main simulation loop along with graphical rendering, mouse interaction, 
    cell color mapping, and vector drawing functions.
    The `run_simulation` function handles grid rendering at each timestep, fluid injection, 
    mouse wall drawing, and various sources (air flow, Kármán vortex, central source).
    It also calls the appropriate velocity step method based on the selected mode ('vel_step', `vel2_step', `cip_csl4').

    Key functions :
    - `run_simulation' : time loop, rendering, interaction, injection, physics
    - `vorticity' : computes local vorticity for coloring
    - `density_color' : maps scalar density to an RGB color
    - `bresenham_line` / `draw_line_with_color' : discrete line drawing algorithms

    User interactions :
    - Left-click (drag) : draws walls interactively on the grid
    - P key : pauses/resumes the simulation
*/



/// Vorticity calculation for colourization (safe bounds)
fn vorticity(grid: &Grid, i: usize, j: usize) -> f32 {
    let dx = DX;
    let dy = DY;

    let im = i.saturating_sub(1).min(N as usize);
    let ip = (i + 1).min(N as usize);
    let jm = j.saturating_sub(1).min(N as usize);
    let jp = (j + 1).min(N as usize);

    let idx_up    = grid.to_index(i, jp);
    let idx_down  = grid.to_index(i, jm);
    let idx_left  = grid.to_index(im, j);
    let idx_right = grid.to_index(ip, j);

    let du_dy = (grid.cells[idx_up].velocity_x - grid.cells[idx_down].velocity_x) / (2.0 * dy);
    let dv_dx = (grid.cells[idx_right].velocity_y - grid.cells[idx_left].velocity_y) / (2.0 * dx);

    dv_dx - du_dy
}

/// Draw a line using Bresenham's algorithm
pub fn draw_line_with_color(buffer: &mut [u32], width: usize, x0: usize, y0: usize, x1: usize, y1: usize, color: u32) {
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

/// Map density to color when not painting vorticity
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

/// Function to compute a route between two points using the Bresenham algorithm
pub fn bresenham_line(x0: usize, y0: usize, x1: usize, y1: usize) -> Vec<(usize, usize)> {
    let mut points = Vec::new();

    let dx = (x1 as isize - x0 as isize).abs();
    let dy = (y1 as isize - y0 as isize).abs();

    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };

    let mut err = dx - dy;
    let mut x = x0 as isize;
    let mut y = y0 as isize;

    loop {
        points.push((x as usize, y as usize));

        if x == x1 as isize && y == y1 as isize { break; }

        let e2 = 2 * err;
        if e2 > -dy {
            err -= dy;
            x += sx;
        }
        if e2 < dx {
            err += dx;
            y += sy;
        }
    }

    points
}

/// Launch the simulation and rendering and handle user input
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

    // Enable key repeat for a better user experience
    window.set_key_repeat_delay(0.25);
    window.set_key_repeat_rate(0.05);

    // Variables to track mouse state
    let mut mouse_down = false;
    let mut last_grid_pos: Option<(usize, usize)> = None;
    let mut last_mouse_time = Instant::now();

    // State variable for pausing
    let mut paused = false;

    println!("Starting simulation");

    let n_max = N as usize + 1;
    while window.is_open() {
        // Check for a pause key (P)
        if window.is_key_pressed(minifb::Key::P, minifb::KeyRepeat::No) {
            paused = !paused;
            println!("Simulation {}", if paused { "PAUSED, PRESS 'P' TO RESUME "} else { "RUNNING, PRESS 'P' TO PAUSE" });
        }

        // Handle mouse events
        let mouse_pos = window.get_mouse_pos(minifb::MouseMode::Discard);
        let mouse_pressed = window.get_mouse_down(minifb::MouseButton::Left);

        // Update mouse state
        if mouse_pressed {
            let now = Instant::now();

            if !mouse_down {
                // First detected clic
                mouse_down = true;
                // Reset timer at the beginning of another clic
                last_mouse_time = now;
            }

            // Convert mouse position in grid coordinates
            if let Some((mouse_x, mouse_y)) = mouse_pos {
                let grid_x = (mouse_x as usize / DX as usize).clamp(1, N as usize);
                let grid_y = (mouse_y as usize / DY as usize).clamp(1, N as usize);
                let current_pos = (grid_x, grid_y);

                // If there is a change in position since last clic or first clic
                if last_grid_pos.is_none() || last_grid_pos.unwrap() != current_pos {
                    // Draw a wall on this position
                    let idx = grid.to_index(grid_x, grid_y);
                    grid.cells[idx].wall = true;

                    // Handle mouse drag and draw wall il between detected clics
                    let time_elapsed = now.duration_since(last_mouse_time);
                    if let Some((last_x, last_y)) = last_grid_pos {
                        if time_elapsed.as_secs_f32() < 1.0 {
                            // Uses Bresenham's algorithm to draw a line between the last and current position
                            let line_points = bresenham_line(last_x, last_y, grid_x, grid_y);
                            for (x, y) in line_points {
                                if x >= 1 && x <= N as usize && y >= 1 && y <= N as usize {
                                    let idx = grid.to_index(x, y);
                                    grid.cells[idx].wall = true;
                                }
                            }
                        }
                    }

                    last_grid_pos = Some(current_pos);
                    last_mouse_time = now;
                }
            }
        } else {
            // Mouse is released
            if mouse_down {
                mouse_down = false;
                last_grid_pos = None;
            }
        }

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

                // Choose color based on mode
                let colour = if PAINT_VORTICITY == true{
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

                // Draw a velocity vector
                if DRAW_VELOCITY_VECTORS == true{
                    let vx = grid.cells[idx].velocity_x;
                    let vy = grid.cells[idx].velocity_y;
                    let cx = x0 + (DX as usize) / 2;
                    let cy = y0 + (DY as usize) / 2;
                    let mag = (vx*vx + vy*vy).sqrt().max(1e-5);
                    let scale = VECTOR_SIZE_FACTOR / mag;
                    let ex = (cx as f32 + vx * scale) as usize;
                    let ey = (cy as f32 + vy * scale) as usize;
                    draw_line_with_color(&mut buffer, width, cx, cy, ex, ey, 0x000000);
                }
            }
        }

        // Update window
        window
            .update_with_buffer(&buffer, width, height)
            .expect("Error while updating the window");

        // If simulation is paused, don't update physic simulation
        if paused {
            continue;
        }

        step += 1;

        // Flow injection or custom source
        if AIR_FLOW ==true {
            let hole_pos: Vec<usize> = (1..=N as usize)
                .filter(|&x| x % FLOW_SPACE == 0)
                .collect();
            grid.initialize_wind_tunnel(FLOW_DENSITY, &hole_pos);
        } else if KARMAN_VORTEX == true{
            if GRID != "1"{
                panic!("Karman vortex only works with grid 1");
            }
            for i in 0..=50 {
                Grid::cell_init(grid, 5, 127+i-25, 0.3, 0.0, 20.0);

                // for karman vortex
                Grid::velocity_init(grid, 10, 127+i-25, 0.9, 0.0);
            }
        } else if CENTER_SOURCE == true{
            grid.center_source(CENTER_SOURCE_RADIUS, CENTER_SOURCE_DENSITY, CENTER_SOURCE_VELOCITY, CENTER_SOURCE_TYPE);
        } else {
            continue;
        }


        // Perform simulation step
        if GRID == "1" {
            if VEL_STEP == "1" {
                grid.vel_step();
            } else if VEL_STEP == "2" {
                grid.vel2_step();
            } else if VEL_STEP == "cip_csl4" {
                grid.vel_step_cip_csl4()
            } else {
                panic!("Invalid VEL_STEP value");
            }
        } else {
            panic!("Invalid GRID value");
        }

        if PRINT_FORCES == true{
            if GRID != "1"{
                panic!("PRINT_FORCES only works with grid 1");
            }
            grid.print_object_forces()
        }

        // Optionally print timing (this one especially for performance purposes)
        if step == 100  {
            let elapsed = start.elapsed();
            println!("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n Temps pour 100 vel_step: {:?} \n $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ", elapsed);
        }
        // Optionally print timing (to see time elapsed after SIM_STEPS steps)
        if step == SIM_STEPS as i32 {
            let elapsed = start.elapsed();
            println!("§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ \n Completed {} steps in {:?}\n§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§", SIM_STEPS, elapsed);
        }
    }
}


pub fn run_simulation2(grid2: &mut Grid2) {
    let start = Instant::now();
    let width = WINDOW_WIDTH;
    let height = WINDOW_HEIGHT;
    let mut buffer: Vec<u32> = vec![0; width * height];

    let mut window = Window::new(
        "Grid Simulation - Morton Staggered Grid",
        width,
        height,
        WindowOptions::default(),
    ).unwrap_or_else(|e| panic!("{}", e));

    // Enable key repeat for a better user experience
    window.set_key_repeat_delay(0.25);
    window.set_key_repeat_rate(0.05);

    // Variables to track mouse state
    let mut mouse_down = false;
    let mut last_grid_pos: Option<(usize, usize)> = None;
    let mut last_mouse_time = Instant::now();

    // State variable for pausing
    let mut paused = false;
    // Step counter
    let mut step = 0;

    println!("Starting simulation with Morton Staggered Grid");

    while window.is_open() {
        let gravity = gravity();
        // Check for a pause key (P)
        if window.is_key_pressed(minifb::Key::P, minifb::KeyRepeat::No) {
            paused = !paused;
            println!("Simulation {}", if paused { "PAUSED, PRESS 'P' TO RESUME "} else { "RUNNING, PRESS 'P' TO PAUSE" });
        }
        if window.is_key_pressed(minifb::Key::V, minifb::KeyRepeat::Yes){
            grid2.vel_step(DT, gravity);
        }
        if window.is_key_pressed(minifb::Key::D, minifb::KeyRepeat::Yes){
            grid2.advect_density(DT);
        }
        if window.is_key_pressed(minifb::Key::I, minifb::KeyRepeat::Yes){
            grid2.solve_incompressibility(DT, 20, 1000.0, 1.9);
        }
        if window.is_key_pressed(minifb::Key::A, minifb::KeyRepeat::Yes){
            grid2.advect_velocity(DT);
        }

        // Handle mouse events for wall creation
        let mouse_pos = window.get_mouse_pos(minifb::MouseMode::Discard);
        let mouse_pressed = window.get_mouse_down(minifb::MouseButton::Left);

        // Update mouse state
        if mouse_pressed {
            let now = Instant::now();

            if !mouse_down {
                // First detected click
                mouse_down = true;
                // Reset timer at the beginning of another click
                last_mouse_time = now;
            }

            // Convert mouse position to grid coordinates
            if let Some((mouse_x, mouse_y)) = mouse_pos {
                let cell_size = grid2.cell_size;
                let grid_x = (mouse_x as f32 / cell_size).floor() as usize;
                let grid_y = (mouse_y as f32 / cell_size).floor() as usize;
                let current_pos = (grid_x, grid_y);

                // If there is a change in position since last click or first click
                if last_grid_pos.is_none() || last_grid_pos.unwrap() != current_pos {
                    // Check if position is within grid bounds
                    if grid2.in_bounds(grid_x, grid_y) {
                        // Set the cell as solid/wall
                        if let Some(cell) = grid2.try_cell_mut_idx(grid_x, grid_y) {
                            cell.set_solid(true);
                        }

                        // Handle mouse drag and draw wall between detected clicks
                        let time_elapsed = now.duration_since(last_mouse_time);
                        if let Some((last_x, last_y)) = last_grid_pos {
                            if time_elapsed.as_secs_f32() < 1.0 {
                                // Uses Bresenham's algorithm to draw a line between the last and current position
                                let line_points = bresenham_line(last_x, last_y, grid_x, grid_y);
                                for (x, y) in line_points {
                                    if grid2.in_bounds(x, y) {
                                        if let Some(cell) = grid2.try_cell_mut_idx(x, y) {
                                            cell.set_solid(true);
                                        }
                                    }
                                }
                            }
                        }

                        last_grid_pos = Some(current_pos);
                        last_mouse_time = now;
                    }
                }
            }
        } else {
            // Mouse is released
            if mouse_down {
                mouse_down = false;
                last_grid_pos = None;
            }
        }

        // Clear buffer
        buffer.fill(0xFFFFFFFF);

        // Draw cells - using Morton order iteration
        for (i, j, idx) in grid2.iter_morton() {
            // Skip cells outside the simulation area (for optimization)
            if !grid2.in_simulation_bounds(i, j) {
                continue;
            }

            let cell = &grid2.cells[idx];
            let x0 = i * grid2.cell_size as usize;
            let y0 = j * grid2.cell_size as usize;
            let cell_size = grid2.cell_size as usize;

            // Draw walls in black
            if cell.is_solid() {
                for dy in 0..cell_size {
                    for dx in 0..cell_size {
                        if y0 + dy < height && x0 + dx < width {
                            buffer[(y0 + dy) * width + (x0 + dx)] = 0x000000;
                        }
                    }
                }
                continue;
            }

            // Skip empty density
            if cell.density.back == 0.0 {
                continue;
            }

            // Choose color based on mode
            let colour = if PAINT_VORTICITY == true {
                let vort = calculate_vorticity2(grid2, i, j);
                let iv = ((vort.abs().min(5.0)) * 25.0) as u32;
                if vort > 0.0 {
                    (iv << 8) | 0x0000FF  // Blue for positive vorticity
                } else {
                    (iv << 16) | 0xFF0000  // Red for negative vorticity
                }
            } else {
                density_color(cell.density.back)
            };

            // Fill cell block
            for dy in 0..cell_size {
                for dx in 0..cell_size {
                    if y0 + dy < height && x0 + dx < width {
                        buffer[(y0 + dy) * width + (x0 + dx)] = colour;
                    }
                }
            }

            // Draw velocity vectors
            if DRAW_VELOCITY_VECTORS == true {
                let vx = cell.velocity.back.x;
                let vy = cell.velocity.back.y;
                let cx = x0 + cell_size / 2;
                let cy = y0 + cell_size / 2;
                let mag = (vx*vx + vy*vy).sqrt().max(1e-5);
                let scale = VECTOR_SIZE_FACTOR / mag;
                let ex = (cx as f32 + vx * scale) as usize;
                let ey = (cy as f32 + vy * scale) as usize;

                // Ensure points are within window bounds
                if cx < width && cy < height && ex < width && ey < height {
                    draw_line_with_color(&mut buffer, width, cx, cy, ex, ey, 0x000000);
                }
            }
        }

        // Update window
        window
            .update_with_buffer(&buffer, width, height)
            .expect("Error while updating the window");

        // If simulation is paused, don't update physic simulation
        if paused {
            continue;
        }

        step += 1;

        // Apply sources based on settings
        if AIR_FLOW == true{
            let _cell_size = grid2.cell_size;
            let hole_pos: Vec<usize> = (1..=N as usize)
                .filter(|&x| x % FLOW_SPACE == 0)
                .collect();
            
            grid2.initialize_wind_tunnel(FLOW_DENSITY, &*hole_pos, FLOW_VELOCITY);

            /*// Apply flow at inlet positions
            for &i in &flow_pos {
                if grid2.in_simulation_bounds(i, 1) {
                    if let Some(cell) = grid2.try_cell_mut_idx(i, 1) {
                        if !cell.is_solid() {
                            cell.density.back = FLOW_DENSITY;
                            cell.velocity.back.x = FLOW_VELOCITY;
                            cell.velocity.back.y = 0.0;
                        }
                    }
                }
            }*/
        } else if KARMAN_VORTEX == true{
            // Position for Karman vortex shedding (need to adjust for Grid2)
            let center_j = grid2.height / 2;

            // Create an obstacle if not already present
            if step == 1 {
                let cx = grid2.width as f32 * 0.3;
                let cy = grid2.height as f32 * 0.5;
                let radius = grid2.cell_size * 5.0;
                grid2.add_circle_obstacle(cx, cy, radius, None);
            }

            // Add inlet flow
            for j in (center_j-25)..(center_j+25) {
                if grid2.in_simulation_bounds(5, j) {
                    if let Some(cell) = grid2.try_cell_mut_idx(5, j) {
                        if !cell.is_solid() {
                            cell.density.back = 0.3;
                            cell.velocity.back.x = 0.9;
                            cell.velocity.back.y = 0.0;
                        }
                    }
                }
            }
        } else if CENTER_SOURCE == true{
            // Add central source
            let center_i = grid2.width / 2;
            let center_j = grid2.height / 2;
            let radius = CENTER_SOURCE_RADIUS as usize;

            // Apply a central source
            for j in (center_j-radius)..(center_j+radius) {
                for i in (center_i-radius)..(center_i+radius) {
                    if grid2.in_simulation_bounds(i, j) {
                        // Calculate distance from a center
                        let di = i as f32 - center_i as f32;
                        let dj = j as f32 - center_j as f32;
                        let dist_sq = di*di + dj*dj;

                        if dist_sq <= (radius as f32 * radius as f32) {
                            if let Some(cell) = grid2.try_cell_mut_idx(i, j) {
                                if !cell.is_solid() {
                                    cell.density.back = CENTER_SOURCE_DENSITY;

                                    // Different source types
                                    match CENTER_SOURCE_TYPE { //"false" for a circular output, "true" for a radial output (like a fan)
                                        true => {
                                            let dist = dist_sq.sqrt().max(0.1);
                                            cell.velocity.back.x = di / dist * CENTER_SOURCE_VELOCITY;
                                            cell.velocity.back.y = dj / dist * CENTER_SOURCE_VELOCITY;
                                        },
                                        false => {
                                            cell.velocity.back.x = -dj / radius as f32 * CENTER_SOURCE_VELOCITY;
                                            cell.velocity.back.y = di / radius as f32 * CENTER_SOURCE_VELOCITY;
                                        } 
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Perform a simulation step with Grid2

        //grid2.cell_init(100,100,0.0,0.0,30.0);
        
        //grid2.vel_step(DT, gravity);

        /*grid2.apply_forces(DT, gravity);
        grid2.solve_incompressibility(DT, 20, 1000.0, 1.9);
        grid2.advect_velocity(DT);
        grid2.advect_density(DT);
        grid2.apply_boundary_conditions();*/

        // Optionally print forces for an obstacle
        if PRINT_FORCES == true {
            // Calculate and print forces on obstacles
            grid2.print_object_forces();
        }

        // Optionally print timing (for performance analysis)
        if step == 100 {
            let elapsed = start.elapsed();
            println!("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n Time for 100 steps: {:?} \n $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ", elapsed);
        }

        if step == SIM_STEPS as i32 {
            let elapsed = start.elapsed();
            println!("§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ \n Completed {} steps in {:?}\n§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§", SIM_STEPS, elapsed);
        }
    }
}

// Helper function to calculate vorticity at a specific cell in Grid2
fn calculate_vorticity2(grid: &Grid2, i: usize, j: usize) -> f32 {
    // Skip boundary cells
    if !grid.in_simulation_bounds(i, j) {
        return 0.0;
    }

    // Check if we have all required neighboring cells
    if i == 0 || j == 0 || i >= grid.width-1 || j >= grid.height-1 {
        return 0.0;
    }

    // Get velocity at neighboring cells
    let vx_up = grid.cell_idx(i, j+1).velocity.back.x;
    let vx_down = grid.cell_idx(i, j-1).velocity.back.x;
    let vy_right = grid.cell_idx(i+1, j).velocity.back.y;
    let vy_left = grid.cell_idx(i-1, j).velocity.back.y;

    // Calculate vorticity (curl of velocity field)
    let dvy_dx = (vy_right - vy_left) / (2.0 * grid.cell_size);
    let dvx_dy = (vx_up - vx_down) / (2.0 * grid.cell_size);

    return dvy_dx - dvx_dy;
}
