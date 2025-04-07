use crate::conditions::*;
use crate::grid::*;
use minifb::{Window, WindowOptions};
use std::time::Instant;


// minifb function to draw lines
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
    (intensity << 8) | 0x000000 // Walls greener with force
}



//function that run simulation and calls draw
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
    // Number of steps
    println!("Steps: {}", step);

    while window.is_open() {
        // Clear the buffer
        for i in buffer.iter_mut() {
            *i = 0xFFFFFFFF; // White background
        }

        // Draw the grid
        for j in 0..=N + 1 {
            for i in 0..=N + 1 {
                let idx = grid.index(i, j);
                let color =
                    if grid.cells[idx].wall {
                        force_to_color(0.0)
                } else {
                    if grid.cells[idx].density == 0.0 {
                        continue;
                    }
                    let density = grid.cells[idx].density;
                    let pressure = grid.cells[idx].pressure;
                    if density <= 20.0 {
                        let intensity = (255.0 * (1.0 - density / 20.0)).clamp(0.0, 255.0) as u32;
                        (intensity << 16) | (intensity << 8) | 0xFF // Deeper blue with density before 20

                    } else {
                        let excess = (density - 20.0).clamp(0.0, 20.0); // Higher than 20 density
                        let red_intensity = (255.0 * (excess / 20.0)).clamp(0.0, 255.0) as u32;
                        (red_intensity << 16) | 0xFF // Becomes more red after density of 20
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

                    // Draw line from center to end point
                    draw_line(&mut buffer, width, center_x, center_y, end_x, end_y, 0x000000); // the velocity vector
                }
*/

            }
        }

        // Update the window
        window.update_with_buffer(&buffer, width, height).expect("Error while updating the window");


        // Update the grid
        grid.vel_step();


        // Pause for observation
        //std::thread::sleep(std::time::Duration::from_millis(50));


        step += 1;

        if step < 375 || step > 500 {
            Grid::cell_init(grid, 5, 100, 5.0, 0.0, 300.0);
        }






        /* if step % 30 == 0 {
             Grid::cell_init(grid, 5, 5, 1.5, 1.5, 2000.0);
             Grid::cell_init(grid, 3, 35, 1.5, 1.5, 2000.0);
         }*/

        /* if step == 100 {
            Grid::cell_init(grid, 5, 20,5.0, 1.0, 10000.0);
        }*/

        //println!("Global density {:2}", grid.total_density());
        if step == 100 {
            let duration2 = start2.elapsed();
            println!("Temps pour 100 vel_step: {:?}", duration2);
        }

    }
}
