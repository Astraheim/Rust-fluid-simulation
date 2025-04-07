use crate::conditions::*;
use std::vec::Vec;
use rayon::prelude::*;
use std::fmt::Write;


#[derive(Clone, Copy, Debug, Default)]
pub struct Vector2 {
    pub x : f32,
    pub y: f32,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Cell {
    pub velocity: Vector2,
    pub pressure: f32,
    pub density: f32,
    pub wall: bool,
}

#[derive(Clone, Debug)]
pub struct Grid {
    pub cells: Vec<Cell>,
}

impl Vector2 {
    pub fn magnitude(&self) -> f32 {
        (self.x.powi(2) + self.y.powi(2)).sqrt()
    }
}


impl Grid {

    // Creates a new grid with borders set as walls if EXT_BORDER is enabled
    pub fn new() -> Self {
        let mut grid = Self {
            cells: vec![Cell { velocity: Vector2 { x: 0.0, y: 0.0 }, pressure: 0.0, density: 0.0, wall: false }; SIZE],
        };

        if EXT_BORDER {
            for i in 0..=N + 1 {
                let left_wall_index = grid.index(i, 0);
                let right_wall_index = grid.index(i, N + 1);
                grid.cells[left_wall_index].wall = true; // Left wall
                grid.cells[right_wall_index].wall = true; // Right wall
            }
            for j in 0..=N + 1 {
                let top_wall_index = grid.index(0, j);
                let bottom_wall_index = grid.index(N + 1, j);
                grid.cells[top_wall_index].wall = true; // Top wall
                grid.cells[bottom_wall_index].wall = true; // Bottom wall
            }
        }

        grid
    }


    // Returns the index of a cell in the grid based on its coordinates
    pub fn index(&self, i: usize, j: usize) -> usize {
        // todo: switch to morton encoding
         i + (N + 2) * j
    }


    pub fn get_pos(index: usize) -> (usize, usize) {
        let y = index / (N + 2);
        let x = index % (N + 2);
        (x, y)
    }


    // Adds a source to the grid cells' pressure
    pub fn add_source(&mut self, source: &[f32], dt: f32) {
        self.cells.par_iter_mut().enumerate().for_each(|(i, cell)| {
            cell.pressure += dt * source[i];
        });
    }


    // Diffuses the density of the grid cells
    pub fn diffuse(&mut self, diff: f32, dt: f32) {
        let a = dt * diff * (N * N) as f32;
        let mut new_density = vec![0.0; self.cells.len()];

        for _ in 0..20 {
            new_density.par_iter_mut().zip(self.cells.par_iter()).enumerate().for_each(|(i, (new_cell, cell))| {
                let (x, y) = Self::get_pos(i);
                if cell.wall { return; }

                let mut sum = 0.0;
                let mut count = 0;

                for &(di, dj) in &[(-1, 0), (1, 0), (0, -1), (0, 1)] {
                    let ni = (x as isize + di) as usize;
                    let nj = (y as isize + dj) as usize;
                    let n_idx = self.index(ni, nj);

                    if !self.cells[n_idx].wall {
                        sum += self.cells[n_idx].density;
                        count += 1;
                    }
                }

                if count > 0 {
                    *new_cell = (cell.density + a * sum) / (1.0 + a * count as f32);
                }
            });

            for (cell, &new_dens) in self.cells.iter_mut().zip(new_density.iter()) {
                cell.density = new_dens;
            }
        }
    }


    // Advects the velocity of the grid cells

    pub fn advect_velocity(&mut self, dt: f32) {
        let dt0 = dt * N as f32;
        let mut new_velocities = vec![[0.0; 2]; self.cells.len()]; // [vx, vy]
        new_velocities.par_iter_mut().enumerate().for_each(|(idx, new_vel)| {
                let (i, j) = Self::get_pos(idx);
                if self.cells[idx].wall {
                    return;
                }

                let cell = &self.cells[idx];
                let x = (i as f32) - dt0 * cell.velocity.x;
                let y = (j as f32) - dt0 * cell.velocity.y;

                let x = x.clamp(0.5, N as f32 + 0.5);
                let y = y.clamp(0.5, N as f32 + 0.5);

                let i0 = x.floor() as usize;
                let i1 = i0 + 1;
                let j0 = y.floor() as usize;
                let j1 = j0 + 1;

                let s1 = x - i0 as f32;
                let s0 = 1.0 - s1;
                let t1 = y - j0 as f32;
                let t0 = 1.0 - t1;

                let idx00 = self.index(i0, j0);
                let idx01 = self.index(i0, j1);
                let idx10 = self.index(i1, j0);
                let idx11 = self.index(i1, j1);

                let vx = s0 * (t0 * self.cells[idx00].velocity.x + t1 * self.cells[idx01].velocity.x)
                    + s1 * (t0 * self.cells[idx10].velocity.x + t1 * self.cells[idx11].velocity.x);

                let vy = s0 * (t0 * self.cells[idx00].velocity.y + t1 * self.cells[idx01].velocity.y)
                    + s1 * (t0 * self.cells[idx10].velocity.y + t1 * self.cells[idx11].velocity.y);

                new_vel[0] = vx;
                new_vel[1] = vy;
            });

        // Apply the new velocities to the cells
        self.cells.iter_mut().zip(new_velocities.iter()).for_each(|(cell, &vel)| {
                cell.velocity.x = vel[0];
                cell.velocity.y = vel[1];
            });
    }


    // Advects the density of the grid cells
    pub fn advect_density(&mut self, dt: f32) {
        let dt0 = dt * N as f32;
        let total_density_before: f32 = self.cells.par_iter().map(|cell| cell.density).sum();

        let new_cells: Vec<Cell> = self.cells.par_iter().enumerate().map(|(idx, cell)| {
            let i = idx % (N + 2);
            let j = idx / (N + 2);

            if cell.wall || i == 0 || i == N + 1 || j == 0 || j == N + 1 {
                return *cell;  // We don't modify walls
            }

            let x = (i as f32) - dt0 * cell.velocity.x;
            let y = (j as f32) - dt0 * cell.velocity.y;

            let i0 = x.floor() as usize;
            let i1 = i0 + 1;
            let j0 = y.floor() as usize;
            let j1 = j0 + 1;

            let s1 = x - i0 as f32;
            let s0 = 1.0 - s1;
            let t1 = y - j0 as f32;
            let t0 = 1.0 - t1;

            let density = s0 * (t0 * self.cells[self.index(i0, j0)].density
                + t1 * self.cells[self.index(i0, j1)].density)
                + s1 * (t0 * self.cells[self.index(i1, j0)].density
                + t1 * self.cells[self.index(i1, j1)].density);

            Cell { density, ..*cell }
        }).collect();

        // Cheat condition to remove density loss
        // todo : change this shitty solution
        let total_density_after: f32 = new_cells.par_iter().map(|cell| cell.density).sum();
        let correction_factor = if total_density_after > 0.0 {
            total_density_before / total_density_after
        } else {
            1.0
        };

        self.cells.par_iter_mut().zip(new_cells).for_each(|(cell, new_cell)| {
            cell.density = new_cell.density * correction_factor;
            cell.velocity = new_cell.velocity;
        });
    }


    // Projects the velocity field to make it incompressible
    pub fn project(&mut self) {
        let h = 1.0 / N as f32;
        let mut pressure = vec![0.0; SIZE];
        let mut div = vec![0.0; SIZE];

        // Clone cells for parallel access
        let cells = self.cells.clone();
        let index = |i: usize, j: usize| -> usize {
            j * (N + 2) + i
        };

        // Divergence computation
        div.par_iter_mut().enumerate().for_each(|(idx, d)| {
            let (i, j) = Self::get_pos(idx);
            if i == 0 || i > N || j == 0 || j > N || cells[idx].wall {
                *d = 0.0;
            } else {
                *d = -0.5 * h * (
                    cells[index(i + 1, j)].velocity.x - cells[index(i - 1, j)].velocity.x
                        + cells[index(i, j + 1)].velocity.y - cells[index(i, j - 1)].velocity.y
                );
            }
        });


        // Pressure computation
        let tolerance = 1e-5;
        for _ in 0..20 {
            let (new_pressure, local_errors): (Vec<f32>, Vec<f32>) = (0..SIZE).into_par_iter()
                .map(|idx| {
                    let (i, j) = Self::get_pos(idx);
                    if i == 0 || i > N || j == 0 || j > N || cells[idx].wall {
                        return (pressure[idx], 0.0);
                    }

                    let p_new = (div[idx]
                        + pressure[index(i - 1, j)]
                        + pressure[index(i + 1, j)]
                        + pressure[index(i, j - 1)]
                        + pressure[index(i, j + 1)]) / 4.0;

                    let err = (p_new - pressure[idx]).abs();
                    (p_new, err)
                })
                .unzip();

            let max_error = local_errors.into_par_iter().reduce(|| 0.0, f32::max);
            pressure.copy_from_slice(&new_pressure);

            if max_error < tolerance {
                break;
            }
        }


        // Updating velocities
        self.cells.par_iter_mut().enumerate().for_each(|(idx, cell)| {
            let (i, j) = Self::get_pos(idx);
            if i == 0 || i > N || j == 0 || j > N || cell.wall {
                return;
            }

            cell.velocity.x -= 0.5 * (pressure[index(i + 1, j)] - pressure[index(i - 1, j)]) / h;
            cell.velocity.y -= 0.5 * (pressure[index(i, j + 1)] - pressure[index(i, j - 1)]) / h;
        });
    }


    // Computes the wall forces for each cell
    pub fn compute_wall_forces(&self) -> Vec<Vector2> {
        self.cells.par_iter().enumerate().map(|(index, cell)| {
                if cell.wall {
                    let (x, y) = Self::get_pos(index);
                    let mut force = Vector2 { x: 0.0, y: 0.0 };
                    for &(di, dj) in &[(-1, 0), (1, 0), (0, -1), (0, 1)] {
                        let ni = (x as isize + di) as usize;
                        let nj = (y as isize + dj) as usize;
                        let neighbor_idx = self.index(ni, nj);
                        if neighbor_idx < self.cells.len() && !self.cells[neighbor_idx].wall {
                            let neighbor = &self.cells[neighbor_idx];
                            let pressure_force = neighbor.pressure;
                            force.x += pressure_force * di as f32;
                            force.y += pressure_force * dj as f32;
                        }
                    }
                    force
                } else {
                    Vector2 { x: 0.0, y: 0.0 }
                }
            })
            .collect()
    }


    // Initializes the velocity of a specific cell
    pub fn cell_init(&mut self, line: usize, column: usize, vx: f32, vy: f32, density: f32) {
        if line <= N + 1 && column <= N + 1 {
            let idx = self.index(line, column);
            if !self.cells[idx].wall {
                self.cells[idx].velocity = Vector2 { x: vx, y: vy };
                self.cells[idx].density = density;
            } else {
                println!("Impossible to modify a wall!");
            }
        } else {
            println!("Error: indices out of bounds!");
        }
    }

    // Initializes the wall of a specific cell
    pub fn wall_init(&mut self, line: usize, column: usize, wall: bool) {
        if line <= N + 1 && column <= N + 1 {
            let idx = self.index(column, line);
            if !self.cells[idx].wall {
                self.cells[idx].wall = wall;

        } else {
            println!("Impossible to modify a wall!");
        }
    } else {
    println ! ("Error: indices out of bounds!");}
    }


    // Initializes the velocity of a specific cell
    pub fn velocity_init(&mut self, line: usize, column: usize, vx: f32, vy: f32) {
        if line <= N + 1 && column <= N + 1 {
            let idx = self.index(line, column);
            if !self.cells[idx].wall {
                self.cells[idx].velocity = Vector2 { x: vx, y: vy };
            } else {
                println!("Impossible to modify a wall!");
            }
        } else {
            println!("Error: indices out of bounds!");
        }
    }


    // Prints the grid with walls
    pub fn print_grid_walls(&self) {
        let lines: Vec<String> = (0..=N + 1).into_par_iter().map(|j| {
            let mut line = String::with_capacity((N + 2) * 3);
            for i in 0..=N + 1 {
                let idx = self.index(i, j);
                line.push_str(if self.cells[idx].wall { " █ " } else { " · " });
            }
            line.push('\n');
            line
        }).collect();

        for line in lines {
            print!("{line}");
        }
    }


    // Prints the complete grid with velocities
    pub fn print_grid_velocity(&self) {
        let lines: Vec<String> = (0..=N + 1).into_par_iter().map(|j| {
            let mut line = String::with_capacity((N + 2) * 18);
            for i in 0..=N + 1 {
                let idx = self.index(i, j);
                if self.cells[idx].wall {
                    line.push_str("   | █████ |   ");
                } else {
                    let v = &self.cells[idx].velocity;
                    line.push_str(&format!("|{:.3}, {:.3}| ", v.x, v.y));
                }
            }
            line.push('\n');
            line
        }).collect();

        for line in lines {
            print!("{}", line);
        }
    }


    // Print the complete grid with densities
    pub fn print_grid_density(&self) {
        let output: String = (0..=N + 1).into_par_iter().map(|j| {
                let mut line = String::with_capacity((N + 2) * 10);
                for i in 0..=N + 1 {
                    let idx = self.index(i, j);
                    if self.cells[idx].wall {
                        line.push_str("   |█████|   ");
                    } else {
                        let _ = write!(line, "  |{:5.2}|  ", self.cells[idx].density);
                    }
                }
                line.push('\n');
                line
            })
            .collect();

        print!("{output}");
    }


    // Perform the total density in the grid
    pub fn total_density(&self) -> f32 {
        let total_density: f32 = self.cells.iter().map(|cell| cell.density).sum();
        //let mut total_density = 0.0;
        //for cell in &self.cells {
        //    total_density += cell.density;
        //}
        total_density
    }


    // Make a circle of walls in the grid
    pub fn circle (&mut self, center_x : isize, center_y : isize, radius: f32) {
        let center_x = center_x as usize;
        let center_y = center_y as usize;
        for i in 0..=N {
            for j in 0..=N {
                let dx = i as isize - center_x as isize;
                let dy = j as isize - center_y as isize;
                if dx * dx + dy * dy <= (radius * radius) as isize {
                    self.wall_init(j, i, true);
                }
            }
        }
    }


    // Performs a velocity step including diffusion, advection, and projection
    pub fn vel_step(&mut self) {
        self.diffuse(VISCOSITY, DT);
        // println!("Total density after diffuse {:2}", self.total_density());
        self.project();
        // println!("Total density after project {:2}", self.total_density());
        self.advect_velocity(DT);
        // println!("Total density after advect velocity {:2}", self.total_density());
        self.advect_density(DT);
        // println!("Total density after apres advect density {:2}", self.total_density());
        self.project();
    }


}

