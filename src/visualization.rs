use crate::conditions::*;
use crate::grid::*;
use minifb::{Window, WindowOptions};
use std::time::{Duration, Instant};

// Minifb window parameters
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

fn force_to_color(force: f32) -> u32 {
    let intensity = (force.clamp(0.0, 255.0)) as u32;
    (intensity << 8) | 0x000000 // Color intensity in green
}

// Launch the simulation and call draw
pub fn run_simulation(grid: &mut Grid, mut step: i32) {
    let start2 = Instant::now();
    let width = WINDOW_WIDTH;
    let height = WINDOW_HEIGHT;
    let mut buffer: Vec<u32> = vec![0; width * height];

    let mut window = Window::new(
        "Grid Simulation",
        width,
        height,
        WindowOptions::default(),
    )
        .unwrap_or_else(|e| {
            panic!("{}", e);
        });

    //Steps
    println!("Steps: {}", step);

    while window.is_open() {
        // Erase the buffer
        for pixel in buffer.iter_mut() {
            *pixel = 0xFFFFFFFF; // White background
        }

        // Draw the grid
        for j in 0..=(N + 1.0) as usize {
            for i in 0..=(N + 1.0) as usize {
                let idx = grid.to_index(i, j);
                let color = if grid.cells[idx].wall {
                    force_to_color(0.0)
                } else {
                    if grid.cells[idx].density == 0.0 {
                        continue;
                    }
                    let density = grid.cells[idx].density;
                    let pressure = grid.cells[idx].pressure;
                    if density <= 20.0 {
                        let intensity = (255.0 * (1.0 - density / 20.0)).clamp(0.0, 255.0) as u32;
                        (intensity << 16) | (intensity << 8) | 0xFF // Bleu profond pour faible densité
                    } else {
                        let excess = (density - 20.0).clamp(0.0, 20.0);
                        let red_intensity = (255.0 * (excess / 20.0)).clamp(0.0, 255.0) as u32;
                        (red_intensity << 16) | 0xFF // Tendance au rouge pour densité élevée
                    }
                };

                let x = i * 2;
                let y = j * 2;
                for dy in 0..2 {
                    for dx in 0..2 {
                        let buffer_idx = (y + dy) * width + (x + dx);
                        if buffer_idx < buffer.len() && x + dx < width && y + dy < height {
                            buffer[buffer_idx] = color;
                        }
                    }
                }

                /*
                // Draw velocity vectors
                if !grid.cells[idx].wall {
                    let vx = grid.cells[idx].velocity.x;
                    let vy = grid.cells[idx].velocity.y;
                    let center_x = x + 10;
                    let center_y = y + 10;
                    let end_x = (center_x as f32 + vx * 20.0) as usize;
                    let end_y = (center_y as f32 + vy * 20.0) as usize;
                    draw_line(&mut buffer, width, center_x, center_y, end_x, end_y, 0x000000);
                }
                */
            }
        }


        // Update the window with the buffer
        window
            .update_with_buffer(&buffer, width, height)
            .expect("Error while updating the window");


        step += 1;




        // AIRFLOW Gestion

        if AIR_FLOW {
            let hole_pos: Vec<usize> = (1..=(N as usize)).filter(|&x| x % FLOW_SPACE == 0).collect();
            grid.initialize_wind_tunnel(FLOW_DENSITY, &hole_pos);
            // Pause
            //std::thread::sleep(std::time::Duration::from_millis(1000));
        } else {
           // if step < 10 {
                Grid::cell_init(grid, 5, 100, 8.0, 0.0, 350.0);
           // }else {
                //Pause
                //std::thread::sleep(std::time::Duration::from_millis(1000));
           // }
        }




        // println!("Global density {:2}", grid.total_density());
        if step == 100 {
            let duration2 = start2.elapsed();
            println!("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n Temps pour 100 vel_step: {:?} \n $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ", duration2);
        }


        // Perform the grid update step
        grid.vel_step();

    }
}
