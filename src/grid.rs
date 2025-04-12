use crate::conditions::*;
use std::vec::Vec;
use rayon::prelude::*;
use std::fmt::Write;
use crate::cip_csl4_v7_2d::*;

#[derive(Clone, Copy, Debug, Default)]
pub struct Vector2 {
    pub x: f32,
    pub y: f32,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Cell {
    pub velocity: Vector2,
    pub pressure: f32,
    pub density: f32,
    pub density_x: f32,
    pub density_y: f32,
    pub density_xx: f32,
    pub density_yy: f32,
    pub density_xy: f32,
    pub wall: bool,
}

#[derive(Clone, Debug)]
pub struct Grid {
    pub cells: Vec<Cell>,
}

impl Vector2 {

    // Return vector magnitude
    pub fn magnitude(&self) -> f32 {
        (self.x.powi(2) + self.y.powi(2)).sqrt()
    }
}


// Encode Morton coordinates (x, y) into a single index
fn morton_encode(x: usize, y: usize) -> usize {
    part1by1(x) | (part1by1(y) << 1)
}

// Decode Morton index into coordinates (x, y)
fn morton_decode(z: usize) -> (usize, usize) {
    (compact1by1(z), compact1by1(z >> 1))
}

// Interleave bits of x and y (morton_encode)
fn part1by1(mut n: usize) -> usize {
    n &= 0x0000ffff;
    n = (n | (n << 8)) & 0x00ff00ff;
    n = (n | (n << 4)) & 0x0f0f0f0f;
    n = (n | (n << 2)) & 0x33333333;
    n = (n | (n << 1)) & 0x55555555;
    n
}

// Compact bits of n into a single number (morton_decode)
fn compact1by1(mut n: usize) -> usize {
    n &= 0x55555555;
    n = (n | (n >> 1)) & 0x33333333;
    n = (n | (n >> 2)) & 0x0f0f0f0f;
    n = (n | (n >> 4)) & 0x00ff00ff;
    n = (n | (n >> 8)) & 0x0000ffff;
    n
}

impl Grid {

    // Encoding to morton index
    pub fn to_index(&self, i: usize, j: usize) -> usize {
        morton_encode(i, j)
    }

    // Return the coordinates (i, j) from the morton index
    pub fn decode_index(&self, idx: usize) -> (usize, usize) {
        morton_decode(idx)
    }

    // Return Some(idx) if the index is valid, otherwise None
    pub fn try_index(&self, i: usize, j: usize) -> Option<usize> {
        let idx = self.to_index(i, j);
        if idx < self.cells.len() {
            Some(idx)
        } else {
            None
        }
    }

    // Itère sur toutes les cellules dans l'ordre de Morton.
    pub fn iter_morton(&self) -> impl Iterator<Item = (usize, usize, usize)> + '_ {
        self.cells
            .iter()
            .enumerate()
            .map(move |(idx, _)| {
                let (i, j) = self.decode_index(idx);
                (i, j, idx)
            })
    }

    // Version parallèle avec Rayon.
    pub fn iter_morton_par(&self) -> impl ParallelIterator<Item = (usize, usize, usize)> + '_ {
        use rayon::prelude::*;
        self.cells
            .par_iter()
            .enumerate()
            .map(move |(idx, _)| {
                let (i, j) = self.decode_index(idx);
                (i, j, idx)
            })
    }

    // Create a new grid with the specified size (with wall if EXT_BORDER is enabled)
    pub fn new() -> Self {
        let mut grid = Self {
            cells: vec![Cell::default(); SIZE as usize],
        };

        if EXT_BORDER {
            for i in 0..=(N + 1.0) as usize {
                for &j in &[0, (N + 1.0) as usize] {
                    if let Some(idx) = grid.try_index(i, j) {
                        grid.cells[idx].wall = true;
                    }
                }
            }
            for j in 0..=(N + 1.0) as usize {
                for &i in &[0, (N + 1.0) as usize] {
                    if let Some(idx) = grid.try_index(i, j) {
                        grid.cells[idx].wall = true;
                    }
                }
            }
        }

        grid
    }

    // Adds a source to the pressure of the cells
    pub fn add_source(&mut self, source: &[f32], dt: f32) {
        self.cells.par_iter_mut().enumerate().for_each(|(i, cell)| {
            cell.pressure += dt * source[i];
        });
    }

    // Diffuse density in the grid
    pub fn diffuse(&mut self, diff: f32, dt: f32) {
        let a = dt * diff * (N * N) as f32;
        let mut new_density = vec![0.0; self.cells.len()];

        for _ in 0..20 {
            new_density.par_iter_mut()
                .zip(self.cells.par_iter())
                .enumerate()
                .for_each(|(idx, (new_cell, cell))| {
                    let (x, y) = morton_decode(idx);
                    if cell.wall {
                        *new_cell = cell.density; // Les murs conservent leur densité
                        return;
                    }

                    let mut sum = 0.0;
                    let mut count = 0;

                    for &(di, dj) in &[(-1, 0), (1, 0), (0, -1), (0, 1)] {
                        let ni = (x as isize + di) as usize;
                        let nj = (y as isize + dj) as usize;
                        if ni > 0 && ni <= N as usize && nj > 0 && nj <= N as usize {
                            let n_idx = morton_encode(ni, nj);
                            if !self.cells[n_idx].wall {
                                sum += self.cells[n_idx].density;
                                count += 1;
                            } else {
                                // Réflexion ou contribution des murs
                                sum += cell.density;
                                count += 1;
                            }
                        }
                    }

                    if count > 0 {
                        *new_cell = (cell.density + a * sum) / (1.0 + a * count as f32);
                    }
                });

            // Mise à jour des densités
            for (cell, &new_dens) in self.cells.iter_mut().zip(new_density.iter()) {
                cell.density = new_dens;
            }
        }
    }

    // Advect the velocity of the cells in the grid
    pub fn advect_velocity(&mut self, dt: f32) {
        let dt0 = dt * N as f32;
        let mut new_velocities = vec![[0.0; 2]; self.cells.len()];
        new_velocities.par_iter_mut().enumerate().for_each(|(idx, new_vel)| {
            let (i, j) = morton_decode(idx);
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

            let idx00 = morton_encode(i0, j0);
            let idx01 = morton_encode(i0, j1);
            let idx10 = morton_encode(i1, j0);
            let idx11 = morton_encode(i1, j1);

            let vx = s0 * (t0 * self.cells[idx00].velocity.x + t1 * self.cells[idx01].velocity.x)
                + s1 * (t0 * self.cells[idx10].velocity.x + t1 * self.cells[idx11].velocity.x);

            let vy = s0 * (t0 * self.cells[idx00].velocity.y + t1 * self.cells[idx01].velocity.y)
                + s1 * (t0 * self.cells[idx10].velocity.y + t1 * self.cells[idx11].velocity.y);

            new_vel[0] = vx;
            new_vel[1] = vy;
        });

        self.cells.iter_mut().zip(new_velocities.iter()).for_each(|(cell, &vel)| {
            cell.velocity.x = vel[0];
            cell.velocity.y = vel[1];
        });
    }



    // Advect density in the grid
    pub fn advect_density(&mut self, dt: f32) {
        let dt0 = dt * N;

        let new_cells: Vec<Cell> = self.cells.par_iter().enumerate().map(|(idx, cell)| {
            let (i, j) = morton_decode(idx);
            if cell.wall || i == 0 || i == (N + 1.0) as usize || j == 0 || j == (N + 1.0) as usize {
                return *cell;
            }

            let x = (i as f32) - dt0 * cell.velocity.x;
            let y = (j as f32) - dt0 * cell.velocity.y;
            let x = x.clamp(1.0, N);
            let y = y.clamp(1.0, N);

            let i0 = x.floor().clamp(1.0, N as f32) as usize;
            let i1 = (i0 + 1).min((N + 1.0) as usize);
            let j0 = y.floor().clamp(1.0, N as f32) as usize;
            let j1 = (j0 + 1).min((N + 1.0) as usize);

            let s1 = x - i0 as f32;
            let s0 = 1.0 - s1;
            let t1 = y - j0 as f32;
            let t0 = 1.0 - t1;

            let get_density = |i: usize, j: usize| -> f32 {
                if i <= (N + 1.0) as usize && j <= (N + 1.0) as usize {
                    let idx = morton_encode(i, j);
                    if !self.cells[idx].wall {
                        return self.cells[idx].density;
                    }
                }
                0.0 // on ignore les murs ou bordures en considérant 0
            };

            let density =
                s0 * (t0 * get_density(i0, j0) + t1 * get_density(i0, j1)) +
                    s1 * (t0 * get_density(i1, j0) + t1 * get_density(i1, j1));

            Cell { density: density.clamp(0.0, 255.0), ..*cell }
        }).collect();

        self.cells.par_iter_mut().zip(new_cells).for_each(|(cell, new_cell)| {
            cell.density = new_cell.density;
        });
    }

    // Project the velocity field to ensure incompressibility
    pub fn project(&mut self) {
        let h = 1.0 / N as f32;
        let mut pressure = vec![0.0; SIZE as usize];
        let mut div = vec![0.0; SIZE as usize];

        // Compute divergence
        div.par_iter_mut().enumerate().for_each(|(idx, d)| {
            let (i, j) = morton_decode(idx);
            if i == 0 || i > N as usize || j == 0 || j > N as usize|| self.cells[idx].wall {
                *d = 0.0;
            } else {
                *d = -0.5 * h * (
                    self.cells[morton_encode(i + 1, j)].velocity.x - self.cells[morton_encode(i - 1, j)].velocity.x +
                        self.cells[morton_encode(i, j + 1)].velocity.y - self.cells[morton_encode(i, j - 1)].velocity.y
                );
            }
        });

        let tolerance = 1e-5;
        for _ in 0..20 {
            let (new_pressure, local_errors): (Vec<f32>, Vec<f32>) = (0..SIZE as usize).into_par_iter()
                .map(|idx| {
                    let (i, j) = morton_decode(idx);
                    if i == 0 || i > N as usize || j == 0 || j > N as usize|| self.cells[idx].wall {
                        return (pressure[idx], 0.0);
                    }
                    let p_new = (div[idx]
                        + pressure[morton_encode(i - 1, j)]
                        + pressure[morton_encode(i + 1, j)]
                        + pressure[morton_encode(i, j - 1)]
                        + pressure[morton_encode(i, j + 1)]) / 4.0;
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

        // Update velocities
        self.cells.par_iter_mut().enumerate().for_each(|(idx, cell)| {
            let (i, j) = morton_decode(idx);
            if i == 0 || i > N as usize|| j == 0 || j > N as usize || cell.wall {
                return;
            }
            cell.velocity.x -= 0.5 * (pressure[morton_encode(i + 1, j)] - pressure[morton_encode(i - 1, j)]) / h;
            cell.velocity.y -= 0.5 * (pressure[morton_encode(i, j + 1)] - pressure[morton_encode(i, j - 1)]) / h;
        });
    }


    pub fn advect_density_cip_csl4(&mut self, dt: f32) {
        use crate::cip_csl4_v7_2d::advect_cip_csl4_2d;
        advect_cip_csl4_2d( self, dt, DX, DY, (N + 2.0) as usize);
    }



    // Compute wall forces applied to the cells
    pub fn compute_wall_forces(&self) -> Vec<Vector2> {
        self.cells.par_iter().enumerate().map(|(idx, cell)| {
            if cell.wall {
                let (x, y) = morton_decode(idx);
                let mut force = Vector2 { x: 0.0, y: 0.0 };
                for &(di, dj) in &[(-1, 0), (1, 0), (0, -1), (0, 1)] {
                    let ni = (x as isize + di) as usize;
                    let nj = (y as isize + dj) as usize;
                    let neighbor_idx = morton_encode(ni, nj);
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
        }).collect()
    }

    // Initialize the grid with a given velocity and density
    pub fn cell_init(&mut self, line: usize, column: usize, vx: f32, vy: f32, density: f32) {
        if line <= (N + 1.0) as usize && column <= (N + 1.0) as usize {
            let idx = self.to_index(line, column);
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

    // Initialize the wall character of a cell
    pub fn wall_init(&mut self, line: usize, column: usize, wall: bool) {
        if line <= (N + 1.0) as usize && column <= (N + 1.0) as usize {
            let idx = self.to_index(column, line);
            if !self.cells[idx].wall {
                self.cells[idx].wall = wall;
            } else {
                println!("Impossible to modify a wall!");
            }
        } else {
            println!("Error: indices out of bounds!");
        }
    }

    // Initialize the velocity of a cell
    pub fn velocity_init(&mut self, line: usize, column: usize, vx: f32, vy: f32) {
        if line <= (N + 1.0) as usize && column <= (N + 1.0) as usize {
            let idx = self.to_index(line, column);
            if !self.cells[idx].wall {
                self.cells[idx].velocity = Vector2 { x: vx, y: vy };
            } else {
                println!("Impossible to modify a wall!");
            }
        } else {
            println!("Error: indices out of bounds!");
        }
    }

    // Print the grid with walls
    pub fn print_grid_walls(&self) {
        let lines: Vec<String> = (0..=(N + 1.0) as usize).into_par_iter().map(|j| {
            let mut line = String::with_capacity(((N + 2.0) as usize) * 3);
            for i in 0..=(N + 1.0) as usize {
                let idx = self.to_index(i, j);
                line.push_str(if self.cells[idx].wall { " █ " } else { " · " });
            }
            line.push('\n');
            line
        }).collect();

        for line in lines {
            print!("{line}");
        }
    }

    // Print the grid with velocities
    pub fn print_grid_velocity(&self) {
        let lines: Vec<String> = (0..=(N + 1.0) as usize).into_par_iter().map(|j| {
            let mut line = String::with_capacity(((N + 2.0) as usize) * 18);
            for i in 0..=(N + 1.0) as usize {
                let idx = self.to_index(i, j);
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
            print!("{line}");
        }
    }

    // Print the grid with densities
    pub fn print_grid_density(&self) {
        let output: String = (0..=(N + 1.0) as usize).into_par_iter().map(|j| {
            let mut line = String::with_capacity(((N + 2.0) as usize) * 10);
            for i in 0..=(N + 1.0) as usize {
                let idx = self.to_index(i, j);
                if self.cells[idx].wall {
                    line.push_str("   |█████|   ");
                } else {
                    let _ = write!(line, "  |{:5.2}|  ", self.cells[idx].density);
                }
            }
            line.push('\n');
            line
        }).collect();

        print!("{output}");
    }

    // Compute the total density of the grid
    pub fn total_density(&self) -> f32 {
        self.cells.iter().map(|cell| cell.density).sum()
    }

    // Create a circle in the grid
    pub fn circle(&mut self, center_x: isize, center_y: isize, radius: f32) {
        let center_x = center_x as usize;
        let center_y = center_y as usize;
        for i in 0..= N as usize{
            for j in 0..= N as usize {
                let dx = i as isize - center_x as isize;
                let dy = j as isize - center_y as isize;
                if dx * dx + dy * dy <= (radius * radius) as isize {
                    self.wall_init(j, i, true);
                }
            }
        }
    }

    pub fn initialize_wind_tunnel(&mut self, density: f32, hole_positions: &[usize]) {
        let left_wall = 1;
        let box_width = (N as usize / 20).max(2); // Largeur de la boîte (5% de N, minimum 2)
        let right_wall = left_wall + box_width;

        // Créer les murs extérieurs de la boîte
        for j in 1..=N as usize {
            let idx_left = self.to_index(left_wall, j);
            if !self.cells[idx_left].wall {
                self.wall_init(j, left_wall, true); // Mur gauche
            }

            if !hole_positions.contains(&j) {
                let idx_right = self.to_index(right_wall, j);
                if !self.cells[idx_right].wall {
                    self.wall_init(j, right_wall, true); // Mur droit, sauf aux trous
                }
            }
        }

        // Ajouter une densité élevée dans la boîte
        for i in (left_wall + 2)..=right_wall - 5 {
            for j in 2..=N as usize - 1 {
                let idx = self.to_index(i, j);
                if !self.cells[idx].wall { // Vérifie que ce n'est pas un mur
                    self.cell_init(i, j, FLOW_VELOCITY, 0.0, density);
                }
            }
        }
    }


    // Perform a velocity step in the simulation
    pub fn vel_step(&mut self) {
        self.diffuse(VISCOSITY, DT);
        //println!("Total density after diffuse {:2}", self.total_density());
        self.project();
        //println!("Total density after project {:2}", self.total_density());
        self.advect_velocity(DT);
        //println!("Total density after advect velocity {:2}", self.total_density());
        if CIP_CSL4 {
            self.advect_density_cip_csl4(DT);
        } else {
            self.advect_density(DT);
        }
        println!("Total density after apres advect density {:2}", self.total_density());
        self.project();
    }
}

