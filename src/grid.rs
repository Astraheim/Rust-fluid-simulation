use crate::conditions::*;
use std::vec::Vec;
use rayon::prelude::*;
use std::fmt::Write;
use std::ops::Mul;
use crate::pressure_computation::*;

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



#[derive(Clone, Debug)]
pub struct ObjectForce {
    pub id: usize,
    pub center_of_mass: Vector2,
    pub total_force: Vector2,
    pub torque: f32,        // Moment de force
    pub cell_count: usize,  // Nombre de cellules dans l'objet
}



impl Mul<f32> for Vector2 {
    type Output = Vector2;

    fn mul(self, rhs: f32) -> Vector2 {
        Vector2 {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl Mul<Vector2> for f32 {
    type Output = Vector2;

    fn mul(self, rhs: Vector2) -> Vector2 {
        Vector2 {
            x: self * rhs.x,
            y: self * rhs.y,
        }
    }
}



impl Vector2 {

    // Return vector magnitude
    pub fn magnitude(&self) -> f32 {
        (self.x.powi(2) + self.y.powi(2)).sqrt()
    }

    pub fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }

    pub fn length(&self) -> f32 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    pub fn normalize(&self) -> Self {
        let len = self.length();
        if len == 0.0 {
            Self::new(0.0, 0.0)
        } else {
            Self::new(self.x / len, self.y / len)
        }
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

    // Check if the coordinates (i, j) are within the grid bounds
    pub fn in_bounds(&self, i: usize, j: usize) -> bool {
        i >= 1 && i <= N as usize && j >= 1 && j <= N as usize
    }

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

    // Calcul de la moyenne de la vitesse verticale v entre deux cellules adjacentes (cellules voisines)
    fn avg_v(&self, i: usize, j: usize) -> f32 {
        let n = N + 1.0 ;
        0.5 * (self.cells[i * n as usize + j].velocity.y + self.cells[i * n as usize + j + 1].velocity.y)
    }

    // Calcul de la moyenne de la vitesse horizontale u entre deux cellules adjacentes (cellules voisines)
    fn avg_u(&self, i: usize, j: usize) -> f32 {
        let n = N;
        0.5 * (self.cells[i * n as usize + j].velocity.x + self.cells[(i + 1) * n as usize + j].velocity.x)
    }



    // Advect density in the grid
    pub fn advect_density(&mut self, dt: f32) {
        let dt0 = dt * N as f32;

        // Save total density before advection
        // let total_density_before = self.total_density();



        // New densities after advection
        let new_densities: Vec<f32> = self.cells.par_iter().enumerate().map(|(idx, cell)| {
            let (i, j) = morton_decode(idx);

            // Preserve density for walls and boundaries
            if cell.wall || i == 0 || i >= N as usize || j == 0 || j >= N as usize {
                return cell.density;
            }

            // Backtracking: where did the fluid particle come from
            let x = i as f32 - dt0 * cell.velocity.x;
            let y = j as f32 - dt0 * cell.velocity.y;

            // Clamping with reflection at boundaries (but not limiting the actual density values)
            let x = x.clamp(0.5, N as f32 + 0.5);
            let y = y.clamp(0.5, N as f32 + 0.5);

            // Indices and weights for bilinear interpolation
            let i0 = x.floor() as usize;
            let i1 = i0 + 1;
            let j0 = y.floor() as usize;
            let j1 = j0 + 1;

            let s1 = x - i0 as f32;
            let s0 = 1.0 - s1;
            let t1 = y - j0 as f32;
            let t0 = 1.0 - t1;

            // Improved function to get density with reflection at walls
            let get_density = |i: usize, j: usize| -> f32 {
                if i < self.cells.len() && j < self.cells.len() {
                    let idx = morton_encode(i, j);
                    if idx < self.cells.len() {
                        if !self.cells[idx].wall {
                            return self.cells[idx].density;
                        } else {
                            // Reflection: use density of current cell
                            return cell.density;
                        }
                    }
                }
                // For out-of-bounds cells, use current cell's density
                cell.density
            };

            // Bilinear interpolation
            let density = s0 * (t0 * get_density(i0, j0) + t1 * get_density(i0, j1)) +
                s1 * (t0 * get_density(i1, j0) + t1 * get_density(i1, j1));

            density.max(0.0) // Avoid negative values but no upper limit
        }).collect();

        // Update densities
        for (i, &new_density) in new_densities.iter().enumerate() {
            self.cells[i].density = new_density;
        }

        /*// Mass conservation correction (optional)
        let total_density_after = self.total_density();
        if total_density_after > 0.0 && total_density_before > 0.0 {
            let correction_factor = total_density_before / total_density_after;
                // Apply correction if difference is significan
            if (correction_factor - 1.0).abs() > 1e-6 {
                for cell in self.cells.iter_mut().filter(|c| !c.wall) {
                    cell.density *= correction_factor;
                }
            }
        }*/

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



    /*pub fn project2(&mut self, num_iters: usize, dt: f32, over_relaxation: f32) {
        let cp = 1.0 * (1.0 / N as f32) / dt; // densité * h / dt, ici densité = 1
        let h = 1.0 / N as f32;

        for _ in 0..num_iters {
            let mut updates = Vec::new();

            for (i, j, idx) in self.iter_morton() {
                if !self.in_bounds(i, j) || self.cells[idx].wall {
                    continue;
                }

                let sx0 = if self.in_bounds(i - 1, j) && !self.cells[self.to_index(i - 1, j)].wall { 1.0 } else { 0.0 };
                let sx1 = if self.in_bounds(i + 1, j) && !self.cells[self.to_index(i + 1, j)].wall { 1.0 } else { 0.0 };
                let sy0 = if self.in_bounds(i, j - 1) && !self.cells[self.to_index(i, j - 1)].wall { 1.0 } else { 0.0 };
                let sy1 = if self.in_bounds(i, j + 1) && !self.cells[self.to_index(i, j + 1)].wall { 1.0 } else { 0.0 };

                let s = sx0 + sx1 + sy0 + sy1;
                if s == 0.0 {
                    continue;
                }

                let u_r = self.cells[self.to_index(i + 1, j)].velocity.x;
                let u_l = self.cells[self.to_index(i, j)].velocity.x;
                let v_t = self.cells[self.to_index(i, j + 1)].velocity.y;
                let v_b = self.cells[self.to_index(i, j)].velocity.y;

                let div = u_r - u_l + v_t - v_b;

                let mut p = -div / s;
                p *= over_relaxation;
                let dp = cp * p;

                updates.push((i, j, idx, dp, sx0, sx1, sy0, sy1));
            }

            // Application des mises à jour (en dehors de l'itération de lecture)
            for (i, j, idx, dp, sx0, sx1, sy0, sy1) in &updates {
                self.cells[*idx].pressure += *dp;

                if *sx0 != 0.0 {
                    self.cells[self.to_index(i - 1, *j)].velocity.x -= sx0 * *dp;
                }
                if *sx1 != 0.0 {
                    self.cells[self.to_index(i + 1, *j)].velocity.x += sx1 * *dp;
                }
                if *sy0 != 0.0 {
                    self.cells[self.to_index(*i, j - 1)].velocity.y -= sy0 * *dp;
                }
                if *sy1 != 0.0 {
                    self.cells[self.to_index(*i, j + 1)].velocity.y += sy1 * *dp;
                }
            }
        }
    }*/



    // Impose u=0, v=0 in obstacles and walls, constant inflow on the left, outflow (∂/∂x=0) on the right and on the top/bottom
    pub fn apply_boundary_conditions(&mut self) {
        let n = N as usize;
        let len = self.cells.len();
        // 1) On clone l’état actuel des vitesses
        let old_vels: Vec<Vector2> = self.cells.iter().map(|c| c.velocity).collect();
        // 2) On prépare un buffer pour les nouvelles vitesses
        let mut new_vels = Vec::with_capacity(len);

        for idx in 0..len {
            let (i, j) = morton_decode(idx);
            let cell = &self.cells[idx];
            let v = if cell.wall {
                // no‑slip
                Vector2::default()
            } else if i == 1 {
                // inflow
                Vector2 { x: INFLOW_VELOCITY, y: 0.0 }
            } else if i == n {
                // outflow : on reprend la vitesse à gauche
                old_vels[morton_encode(i - 1, j)]
            } else if j == 1 {
                // top : ∂u/∂y = 0 ⇒ on prend la cellule du bas
                old_vels[morton_encode(i, j + 1)]
            } else if j == n {
                // bottom : ∂u/∂y = 0 ⇒ on prend la cellule du haut
                old_vels[morton_encode(i, j - 1)]
            } else {
                // pas de BC particulière
                old_vels[idx]
            };
            new_vels.push(v);
        }

        // 3) On écrit les nouvelles vitesses
        for (cell, &v) in self.cells.iter_mut().zip(new_vels.iter()) {
            cell.velocity = v;
        }
    }




    pub fn extrapolate(&mut self) {
        // On suppose que N est la taille logique (hors-bord), donc les bords sont 0 et N+1
        let n = N as usize;

        // Extrapolation verticale pour u (velocity.x)
        for i in 1..=n {
            let y0 = 0;
            let y1 = 1;
            if let (Some(i0), Some(i1)) = (self.try_index(i, y0), self.try_index(i, y1)) {
                if !self.cells[i0].wall && !self.cells[i1].wall {
                    self.cells[i0].velocity.x = self.cells[i1].velocity.x;
                }
            }

            let y0 = n + 1;
            let y1 = n;
            if let (Some(i0), Some(i1)) = (self.try_index(i, y0), self.try_index(i, y1)) {
                if !self.cells[i0].wall && !self.cells[i1].wall {
                    self.cells[i0].velocity.x = self.cells[i1].velocity.x;
                }
            }
        }

        // Extrapolation horizontale pour v (velocity.y)
        for j in 1..=n {
            let x0 = 0;
            let x1 = 1;
            if let (Some(i0), Some(i1)) = (self.try_index(x0, j), self.try_index(x1, j)) {
                if !self.cells[i0].wall && !self.cells[i1].wall {
                    self.cells[i0].velocity.y = self.cells[i1].velocity.y;
                }
            }

            let x0 = n + 1;
            let x1 = n;
            if let (Some(i0), Some(i1)) = (self.try_index(x0, j), self.try_index(x1, j)) {
                if !self.cells[i0].wall && !self.cells[i1].wall {
                    self.cells[i0].velocity.y = self.cells[i1].velocity.y;
                }
            }
        }
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


    // Ajoute une source circulaire de densité et de vitesse au centre de la grille
    pub fn center_source(&mut self, radius: f32, density_value: f32, velocity_magnitude: f32, is_radial: bool) {
        // Calculer le centre de la grille
        let center_x = (N / 2.0) as usize;
        let center_y = (N / 2.0) as usize;

        // Parcourir toutes les cellules dans un carré englobant le cercle
        let r_int = radius.floor() as usize;
        let start_i = if center_x > r_int { center_x - r_int } else { 1 };
        let end_i = (center_x + r_int).min(N as usize);
        let start_j = if center_y > r_int { center_y - r_int } else { 1 };
        let end_j = (center_y + r_int).min(N as usize);

        for i in start_i..=end_i {
            for j in start_j..=end_j {
                // Calculer la distance au centre
                let dx = i as f32 - center_x as f32;
                let dy = j as f32 - center_y as f32;
                let distance_squared = dx * dx + dy * dy;

                // Si la cellule est dans le cercle et n'est pas un mur
                if distance_squared <= radius * radius {
                    let idx = self.to_index(i, j);
                    if !self.cells[idx].wall {
                        // Ajouter de la densité
                        self.cells[idx].density += density_value;

                        // Calculer la direction de la vitesse
                        if is_radial {
                            // Vitesse radiale (vers l'extérieur depuis le centre)
                            if distance_squared > 0.0 {
                                let distance = distance_squared.sqrt();
                                let dir_x = dx / distance;
                                let dir_y = dy / distance;

                                // Ajouter la vitesse (diminuant avec la distance du centre)
                                let factor = 1.0 - (distance / radius); // Plus fort au centre
                                self.cells[idx].velocity.x += dir_x * velocity_magnitude * factor;
                                self.cells[idx].velocity.y += dir_y * velocity_magnitude * factor;
                            }
                        } else {
                            // Vitesse circulaire (tourbillon)
                            if distance_squared > 0.0 {
                                let distance = distance_squared.sqrt();
                                // Direction perpendiculaire (pour rotation)
                                let dir_x = -dy / distance;
                                let dir_y = dx / distance;

                                // Vitesse augmentant avec la distance jusqu'à un certain point
                                let factor = (distance / radius) * (1.0 - distance / radius) * 4.0; // Max au milieu
                                self.cells[idx].velocity.x += dir_x * velocity_magnitude * factor;
                                self.cells[idx].velocity.y += dir_y * velocity_magnitude * factor;
                            }
                        }
                    }
                }
            }
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
        let box_width = (N as usize / 30).max(2); // Largeur de la boîte (5% de N, minimum 2)
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


    pub fn setup_karman_vortex(&mut self) {
        use crate::conditions::*;
        let radius = (N as usize / 6) as f32 -10.5 ; // Rayon de l'obstacle
        let center = ((N as usize + 2) / 6, (N as usize + 2) / 2);

        self.circle(center.0 as isize, center.1 as isize, radius as f32);

        // Écoulement horizontal de gauche à droite avec densité
        for j in 1..(N as usize + 1) {
            for i in 1..15 {
                let idx = self.to_index(i, j);
                if !self.cells[idx].wall {
                    self.cells[idx].density = FLOW_DENSITY;
                    self.cells[idx].velocity = Vector2::new(FLOW_VELOCITY, 0.0);
                }
            }
        }

        println!("✔ Kármán vortex setup complete: central obstacle + inflow.");
    }


    // Perform a velocity step in the simulation
    pub fn vel_step(&mut self) {
        //self.diffuse(VISCOSITY, DT);
        //println!("Total density after diffuse {:2}", self.total_density());
        if PROJECT == "1"{
            self.project();
        } else if PROJECT == "2"{
            //self.project2(80, DT, 1.9);
        }
        //println!("Total density after project {:2}", self.total_density());
        self.advect_velocity(DT);
        //println!("Total density after advect velocity {:2}", self.total_density());
        self.advect_density(DT);
        //println!("Total density after apres advect density {:2}", self.total_density());
        self.apply_boundary_conditions();
    }

    pub fn vel2_step(&mut self) {
        if PROJECT == "1"{
            self.project();
        } else if PROJECT == "2"{
            //self.project2(80, DT, 1.9);
        }
        //println!("Total density after project {:2}", self.total_density());
        self.extrapolate();
        //println!("Total density after extrapolate {:2}", self.total_density());
        self.advect_velocity(DT);
        //println!("Total density after advect velocity {:2}", self.total_density());
        self.advect_density(DT);
        //println!("Total density after apres advect density {:2}", self.total_density());
    }


    pub fn vel_step_cip_csl4(&mut self) {
        // Diffusion si nécessaire
        // self.diffuse(VISCOSITY, DT);

        // Projection pour assurer l'incompressibilité
        if PROJECT == "1"{
            self.project();
        } else if PROJECT == "2" {
            //self.project2(80, DT, 1.9);
        }

        // Advection de la vitesse
        self.advect_velocity(DT);

        // Utiliser la méthode CIP-CSL4 au lieu de l'advection standard
        Self::advect_density_cip_csl4(self, DT);

        // Appliquer les conditions aux limites
        self.apply_boundary_conditions();
    }


}



