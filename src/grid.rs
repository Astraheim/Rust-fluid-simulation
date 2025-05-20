use crate::conditions::*;
use rayon::prelude::*;
use std::fmt::Write;
use std::ops::Mul;
use std::vec::Vec;

#[derive(Clone, Copy, Debug, Default)]
pub struct Vector2 {
    pub x: f32,
    pub y: f32,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Cell {
    // Staggered grid: vitesses stockées sur les interfaces
    pub velocity_x: f32, // vitesse horizontale (stockée sur interfaces E/O)
    pub velocity_y: f32, // vitesse verticale (stockée sur interfaces N/S)
    pub pressure: f32,   // pression au centre
    pub density: f32,    // densité au centre
    pub density_x: f32,
    pub density_y: f32,
    pub density_xx: f32,
    pub density_yy: f32,
    pub density_xy: f32,
    pub wall: bool,
}

// Modifier Grid pour accéder aux vitesses correctement
impl Grid {
    // Nouvelles méthodes pour accéder aux vitesses selon leur position

    /// Obtenir la vitesse horizontale u(i,j) (située sur la face entre (i-1,j) et (i,j))
    pub fn get_u(&self, i: usize, j: usize) -> f32 {
        if i > 0 && i <= (N + 1.0) as usize && j <= (N + 1.0) as usize {
            let idx = self.to_index(i, j);
            self.cells[idx].velocity_x
        } else {
            0.0 // Hors limites
        }
    }

    /// Obtenir la vitesse verticale v(i,j) (située sur la face entre (i,j-1) et (i,j))
    pub fn get_v(&self, i: usize, j: usize) -> f32 {
        if i <= (N + 1.0) as usize && j > 0 && j <= (N + 1.0) as usize {
            let idx = self.to_index(i, j);
            self.cells[idx].velocity_y
        } else {
            0.0 // Hors limites
        }
    }

    /// Obtenir un vecteur de vitesse au centre de la cellule (i,j) par interpolation
    pub fn get_velocity_at_center(&self, i: usize, j: usize) -> Vector2 {
        let u_left = self.get_u(i, j);
        let u_right = self.get_u(i + 1, j);
        let v_bottom = self.get_v(i, j);
        let v_top = self.get_v(i, j + 1);

        Vector2 {
            x: 0.5 * (u_left + u_right),
            y: 0.5 * (v_bottom + v_top),
        }
    }

    /// Définir la vitesse horizontale u(i,j)
    pub fn set_u(&mut self, i: usize, j: usize, value: f32) {
        if i > 0 && i <= (N + 1.0) as usize && j <= (N + 1.0) as usize {
            let idx = self.to_index(i, j);
            self.cells[idx].velocity_x = value;
        }
    }

    /// Définir la vitesse verticale v(i,j)
    pub fn set_v(&mut self, i: usize, j: usize, value: f32) {
        if i <= (N + 1.0) as usize && j > 0 && j <= (N + 1.0) as usize {
            let idx = self.to_index(i, j);
            self.cells[idx].velocity_y = value;
        }
    }
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
    pub torque: f32,
    pub cell_count: usize,
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
    /// Return vector magnitude
    pub fn magnitude(&self) -> f32 {
        (self.x.powi(2) + self.y.powi(2)).sqrt()
    }

    /// Create a new vector
    pub fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }

    /// Exactly the same as magnitude, returns vector magnitude
    pub fn length(&self) -> f32 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    /// Create a new vector with the same direction but normalized
    pub fn normalize(&self) -> Self {
        let len = self.length();
        if len == 0.0 {
            Self::new(0.0, 0.0)
        } else {
            Self::new(self.x / len, self.y / len)
        }
    }
}

/// Encode Morton coordinates (x, y) into a single index
fn morton_encode(x: usize, y: usize) -> usize {
    part1by1(x) | (part1by1(y) << 1)
}

/// Decode Morton index into coordinates (x, y)
fn morton_decode(z: usize) -> (usize, usize) {
    (compact1by1(z), compact1by1(z >> 1))
}

/// Interleave bits of x and y (morton_encode)
fn part1by1(mut n: usize) -> usize {
    n &= 0x0000ffff;
    n = (n | (n << 8)) & 0x00ff00ff;
    n = (n | (n << 4)) & 0x0f0f0f0f;
    n = (n | (n << 2)) & 0x33333333;
    n = (n | (n << 1)) & 0x55555555;
    n
}

/// Compact bits of n into a single number (morton_decode)
fn compact1by1(mut n: usize) -> usize {
    n &= 0x55555555;
    n = (n | (n >> 1)) & 0x33333333;
    n = (n | (n >> 2)) & 0x0f0f0f0f;
    n = (n | (n >> 4)) & 0x00ff00ff;
    n = (n | (n >> 8)) & 0x0000ffff;
    n
}

impl Grid {
    /// Check if the coordinates (i, j) are within the grid bounds
    pub fn in_bounds(&self, i: usize, j: usize) -> bool {
        i >= 1 && i <= N as usize && j >= 1 && j <= N as usize
    }

    /// Encoding to morton index
    pub fn to_index(&self, i: usize, j: usize) -> usize {
        morton_encode(i, j)
    }

    /// Return the coordinates (i, j) from the morton index
    pub fn decode_index(&self, idx: usize) -> (usize, usize) {
        morton_decode(idx)
    }

    /// Return Some(idx) if the index is valid, otherwise None
    pub fn try_index(&self, i: usize, j: usize) -> Option<usize> {
        let idx = self.to_index(i, j);
        if idx < self.cells.len() {
            Some(idx)
        } else {
            None
        }
    }

    /// Iterate on every cell in Morton order
    pub fn iter_morton(&self) -> impl Iterator<Item = (usize, usize, usize)> + '_ {
        self.cells.iter().enumerate().map(move |(idx, _)| {
            let (i, j) = self.decode_index(idx);
            (i, j, idx)
        })
    }

    /// Create a new grid with the specified size (with wall if EXT_BORDER is enabled)
    pub fn new() -> Self {
        let mut grid = Self {
            cells: vec![Cell::default(); SIZE as usize],
        };

        if EXT_BORDER == true {
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

    /// Adds a source to the pressure of the cells
    pub fn add_source(&mut self, source: &[f32], dt: f32) {
        self.cells.par_iter_mut().enumerate().for_each(|(i, cell)| {
            cell.pressure += dt * source[i];
        });
    }

    /// Diffuse density in the grid
    pub fn diffuse(&mut self, diff: f32, dt: f32) {
        let a = dt * diff * (N * N);
        let mut new_density = vec![0.0; self.cells.len()];

        for _ in 0..20 {
            new_density
                .par_iter_mut()
                .zip(self.cells.par_iter())
                .enumerate()
                .for_each(|(idx, (new_cell, cell))| {
                    let (x, y) = morton_decode(idx);
                    if cell.wall {
                        *new_cell = cell.density; // Conservation of density in walls
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
                                sum += cell.density;
                                count += 1;
                            }
                        }
                    }

                    if count > 0 {
                        *new_cell = (cell.density + a * sum) / (1.0 + a * count as f32);
                    }
                });

            // Update densities
            for (cell, &new_dens) in self.cells.iter_mut().zip(new_density.iter()) {
                cell.density = new_dens;
            }
        }
    }

    /// Advect the velocity of the cells in the grid
    // Modification de advect_velocity pour une staggered grid
    pub fn advect_velocity(&mut self, dt: f32) {
        let dt0 = dt * N;
        let mut new_u = vec![0.0; self.cells.len()];
        let mut new_v = vec![0.0; self.cells.len()];

        // Advection de la vitesse horizontale u
        for i in 1..=(N as usize) {
            for j in 1..=(N as usize) {
                let idx = self.to_index(i, j);
                if self.cells[idx].wall {
                    new_u[idx] = 0.0;
                    continue;
                }

                // Position de u(i,j) sur la grille
                let pos_x = i as f32; // face entre (i-1,j) et (i,j)
                let pos_y = j as f32 + 0.5; // centré verticalement

                // Calculer la vitesse interpolée au point (pos_x, pos_y)
                let u_center = self.get_u(i, j);

                // Déterminer les vitesses moyennes aux centres des cellules adjacentes
                let v_nw = self.get_v(i - 1, j);
                let v_ne = self.get_v(i, j);
                let v_sw = self.get_v(i - 1, j + 1);
                let v_se = self.get_v(i, j + 1);
                let v_center = 0.25 * (v_nw + v_ne + v_sw + v_se);

                // Backtracking
                let x = pos_x - dt0 * u_center;
                let y = pos_y - dt0 * v_center;

                // Limiter aux frontières
                let x = x.clamp(0.5, N + 0.5);
                let y = y.clamp(0.5, N + 0.5);

                // Interpolation bilinéaire
                let i0 = x.floor() as usize;
                let i1 = i0 + 1;
                let j0 = y.floor() as usize;
                let j1 = j0 + 1;

                let s1 = x - i0 as f32;
                let s0 = 1.0 - s1;
                let t1 = y - j0 as f32;
                let t0 = 1.0 - t1;

                // Interpoler la vitesse u
                let sample_u = |i: usize, j: usize| -> f32 {
                    if i > 0 && i <= N as usize && j > 0 && j < N as usize {
                        let idx = self.to_index(i, j);
                        if !self.cells[idx].wall {
                            return self.cells[idx].velocity_x;
                        }
                    }
                    0.0 // Paroi ou hors limites
                };

                new_u[idx] = s0 * (t0 * sample_u(i0, j0) + t1 * sample_u(i0, j1))
                    + s1 * (t0 * sample_u(i1, j0) + t1 * sample_u(i1, j1));
            }
        }

        // Advection de la vitesse verticale v
        for i in 1..=(N as usize) {
            for j in 1..=(N as usize) {
                let idx = self.to_index(i, j);
                if self.cells[idx].wall {
                    new_v[idx] = 0.0;
                    continue;
                }

                // Position de v(i,j) sur la grille
                let pos_x = i as f32 + 0.5; // centré horizontalement
                let pos_y = j as f32; // face entre (i,j-1) et (i,j)

                // Calculer la vitesse interpolée au point (pos_x, pos_y)
                let v_center = self.get_v(i, j);

                // Déterminer les vitesses moyennes aux centres des cellules adjacentes
                let u_nw = self.get_u(i, j - 1);
                let u_ne = self.get_u(i + 1, j - 1);
                let u_sw = self.get_u(i, j);
                let u_se = self.get_u(i + 1, j);
                let u_center = 0.25 * (u_nw + u_ne + u_sw + u_se);

                // Backtracking
                let x = pos_x - dt0 * u_center;
                let y = pos_y - dt0 * v_center;

                // Limiter aux frontières
                let x = x.clamp(0.5, N + 0.5);
                let y = y.clamp(0.5, N + 0.5);

                // Interpolation bilinéaire
                let i0 = x.floor() as usize;
                let i1 = i0 + 1;
                let j0 = y.floor() as usize;
                let j1 = j0 + 1;

                let s1 = x - i0 as f32;
                let s0 = 1.0 - s1;
                let t1 = y - j0 as f32;
                let t0 = 1.0 - t1;

                // Interpoler la vitesse v
                let sample_v = |i: usize, j: usize| -> f32 {
                    if i > 0 && i < N as usize && j > 0 && j <= N as usize {
                        let idx = self.to_index(i, j);
                        if !self.cells[idx].wall {
                            return self.cells[idx].velocity_y;
                        }
                    }
                    0.0 // Paroi ou hors limites
                };

                new_v[idx] = s0 * (t0 * sample_v(i0, j0) + t1 * sample_v(i0, j1))
                    + s1 * (t0 * sample_v(i1, j0) + t1 * sample_v(i1, j1));
            }
        }

        // Mettre à jour les vitesses
        for idx in 0..self.cells.len() {
            self.cells[idx].velocity_x = new_u[idx];
            self.cells[idx].velocity_y = new_v[idx];
        }
    }

    /// Process average vertical velocity
    fn avg_v(&self, i: usize, j: usize) -> f32 {
        let n = N + 1.0;
        0.5 * (self.cells[i * n as usize + j].velocity_y
            + self.cells[i * n as usize + j + 1].velocity_y)
    }

    /// Process average horizontal velocity
    fn avg_u(&self, i: usize, j: usize) -> f32 {
        let n = N;
        0.5 * (self.cells[i * n as usize + j].velocity_x
            + self.cells[(i + 1) * n as usize + j].velocity_x)
    }

    /// Advect density in the grid
    // Modification de advect_density pour une staggered grid
    pub fn advect_density(&mut self, dt: f32) {
        let dt0 = dt * N;

        // Nouvelles densités après advection
        let new_densities: Vec<f32> = self
            .cells
            .par_iter()
            .enumerate()
            .map(|(idx, cell)| {
                let (i, j) = morton_decode(idx);

                // Préserver la densité pour les murs et les frontières
                if cell.wall || i == 0 || i > N as usize || j == 0 || j > N as usize {
                    return cell.density;
                }

                // Obtenir les vitesses au centre de la cellule (i,j) par interpolation
                let u_left = self.get_u(i, j);
                let u_right = self.get_u(i + 1, j);
                let v_bottom = self.get_v(i, j);
                let v_top = self.get_v(i, j + 1);

                let u_center = 0.5 * (u_left + u_right);
                let v_center = 0.5 * (v_bottom + v_top);

                // Backtracking: d'où venait la particule de fluide
                let x = i as f32 - dt0 * u_center;
                let y = j as f32 - dt0 * v_center;

                // Limiter avec réflexion aux frontières
                let x = x.clamp(0.5, N + 0.5);
                let y = y.clamp(0.5, N + 0.5);

                // Indices et poids pour l'interpolation bilinéaire
                let i0 = x.floor() as usize;
                let i1 = i0 + 1;
                let j0 = y.floor() as usize;
                let j1 = j0 + 1;

                let s1 = x - i0 as f32;
                let s0 = 1.0 - s1;
                let t1 = y - j0 as f32;
                let t0 = 1.0 - t1;

                // Fonction améliorée pour obtenir la densité avec réflexion aux murs
                let get_density = |i: usize, j: usize| -> f32 {
                    if i <= N as usize && j <= N as usize {
                        let idx = self.to_index(i, j);
                        if idx < self.cells.len() {
                            return if !self.cells[idx].wall {
                                self.cells[idx].density
                            } else {
                                // Réflexion: utiliser la densité de la cellule actuelle
                                cell.density
                            };
                        }
                    }
                    // Pour les cellules hors limites, utiliser la densité de la cellule actuelle
                    cell.density
                };

                // Interpolation bilinéaire
                let density = s0 * (t0 * get_density(i0, j0) + t1 * get_density(i0, j1))
                    + s1 * (t0 * get_density(i1, j0) + t1 * get_density(i1, j1));

                density.max(0.0) // Éviter les valeurs négatives
            })
            .collect();

        // Mettre à jour les densités
        for (i, &new_density) in new_densities.iter().enumerate() {
            self.cells[i].density = new_density;
        }
    }

    /// Project the velocity field to ensure incompressibility
    // Modification de project pour une staggered grid
    pub fn project(&mut self) {
        let h = 1.0 / N;
        let mut pressure = vec![0.0; SIZE as usize];
        let mut div = vec![0.0; SIZE as usize];

        // Calculer la divergence
        for i in 1..=(N as usize) {
            for j in 1..=(N as usize) {
                let idx = self.to_index(i, j);
                if self.cells[idx].wall {
                    div[idx] = 0.0;
                    continue;
                }

                // Calculer div = ∂u/∂x + ∂v/∂y
                // Pour une staggered grid, la divergence utilise directement les vitesses aux interfaces
                let u_right = self.get_u(i + 1, j); // u à droite de la cellule (i,j)
                let u_left = self.get_u(i, j); // u à gauche de la cellule (i,j)
                let v_top = self.get_v(i, j + 1); // v en haut de la cellule (i,j)
                let v_bottom = self.get_v(i, j); // v en bas de la cellule (i,j)

                div[idx] = -h * ((u_right - u_left) + (v_top - v_bottom));
            }
        }

        // Résoudre l'équation de Poisson pour la pression
        let tolerance = 1e-5;
        for _ in 0..20 {
            let mut max_error: f32 = 0.0;

            for i in 1..=(N as usize) {
                for j in 1..=(N as usize) {
                    let idx = self.to_index(i, j);
                    if self.cells[idx].wall {
                        continue;
                    }

                    // Calculer la nouvelle pression en utilisant Gauss-Seidel
                    let p_left = if i > 1 && !self.cells[self.to_index(i - 1, j)].wall {
                        pressure[self.to_index(i - 1, j)]
                    } else {
                        pressure[idx] // Réflexion à la paroi
                    };

                    let p_right = if i < N as usize && !self.cells[self.to_index(i + 1, j)].wall {
                        pressure[self.to_index(i + 1, j)]
                    } else {
                        pressure[idx] // Réflexion à la paroi
                    };

                    let p_bottom = if j > 1 && !self.cells[self.to_index(i, j - 1)].wall {
                        pressure[self.to_index(i, j - 1)]
                    } else {
                        pressure[idx] // Réflexion à la paroi
                    };

                    let p_top = if j < N as usize && !self.cells[self.to_index(i, j + 1)].wall {
                        pressure[self.to_index(i, j + 1)]
                    } else {
                        pressure[idx] // Réflexion à la paroi
                    };

                    let p_new = (div[idx] + p_left + p_right + p_bottom + p_top) / 4.0;
                    let err = (p_new - pressure[idx]).abs();
                    max_error = max_error.max(err);
                    pressure[idx] = p_new;
                }
            }

            if max_error < tolerance {
                break;
            }
        }

        // Appliquer le gradient de pression pour corriger les vitesses
        for i in 1..=(N as usize + 1) {
            for j in 1..=(N as usize) {
                if i <= N as usize {
                    // Ne pas aller au-delà de la dernière cellule horizontale
                    let idx = self.to_index(i, j);
                    if !self.cells[idx].wall {
                        // Corriger u(i,j) (vitesse horizontale à la face de gauche de la cellule)
                        let p_left = if i > 1 {
                            pressure[self.to_index(i - 1, j)]
                        } else {
                            pressure[idx] // Réflexion à la frontière
                        };
                        let p_right = pressure[idx];

                        // Gradient de pression: ∂p/∂x
                        self.cells[idx].velocity_x -= (p_right - p_left) / h;
                    }
                }

                if j <= N as usize {
                    // Ne pas aller au-delà de la dernière cellule verticale
                    let idx = self.to_index(i, j);
                    if !self.cells[idx].wall {
                        // Corriger v(i,j) (vitesse verticale à la face du bas de la cellule)
                        let p_bottom = if j > 1 {
                            pressure[self.to_index(i, j - 1)]
                        } else {
                            pressure[idx] // Réflexion à la frontière
                        };
                        let p_top = pressure[idx];

                        // Gradient de pression: ∂p/∂y
                        self.cells[idx].velocity_y -= (p_top - p_bottom) / h;
                    }
                }
            }
        }
    }

    /// Project the velocity field to ensure incompressibility (alternative method)
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

            // Updates
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

    /// Impose u=0, v=0 in obstacles and walls, constant inflow on the left, outflow (∂/∂x=0) on the right and on the top/bottom
    // Modification de apply_boundary_conditions pour une staggered grid
    pub fn apply_boundary_conditions(&mut self) {
        let n = N as usize;

        // Conditions aux limites pour la vitesse horizontale u
        for j in 1..=n {
            // Entrée à gauche (inflow)
            self.set_u(1, j, INFLOW_VELOCITY);

            // Sortie à droite (outflow): ∂u/∂x = 0
            self.set_u(n + 1, j, self.get_u(n, j));

            // Parois (no-slip)
            for i in 1..=n + 1 {
                let idx = self.to_index(i, j);
                if self.cells[idx].wall {
                    self.set_u(i, j, 0.0);
                }
            }
        }

        // Conditions aux limites pour la vitesse verticale v
        for i in 1..=n {
            // Haut: ∂v/∂y = 0
            self.set_v(i, 1, self.get_v(i, 2));

            // Bas: ∂v/∂y = 0
            self.set_v(i, n + 1, self.get_v(i, n));

            // Parois (no-slip)
            for j in 1..=n + 1 {
                let idx = self.to_index(i, j);
                if self.cells[idx].wall {
                    self.set_v(i, j, 0.0);
                }
            }
        }

        // Conditions spéciales aux coins

        // Coin supérieur gauche
        if !self.cells[self.to_index(1, 1)].wall {
            self.set_u(1, 1, INFLOW_VELOCITY);
            self.set_v(1, 1, 0.0);
        }

        // Coin supérieur droit
        if !self.cells[self.to_index(n, 1)].wall {
            self.set_u(n + 1, 1, self.get_u(n, 1));
            self.set_v(n, 1, 0.0);
        }

        // Coin inférieur gauche
        if !self.cells[self.to_index(1, n)].wall {
            self.set_u(1, n, INFLOW_VELOCITY);
            self.set_v(1, n + 1, 0.0);
        }

        // Coin inférieur droit
        if !self.cells[self.to_index(n, n)].wall {
            self.set_u(n + 1, n, self.get_u(n, n));
            self.set_v(n, n + 1, 0.0);
        }
    }

    /// Extrapolate the velocity field at the boundaries
    pub fn extrapolate(&mut self) {
        let n = N as usize;

        // Vertical u extrapolation
        for i in 1..=n {
            let y0 = 0;
            let y1 = 1;
            if let (Some(i0), Some(i1)) = (self.try_index(i, y0), self.try_index(i, y1)) {
                if !self.cells[i0].wall && !self.cells[i1].wall {
                    self.cells[i0].velocity_x = self.cells[i1].velocity_x;
                }
            }

            let y0 = n + 1;
            let y1 = n;
            if let (Some(i0), Some(i1)) = (self.try_index(i, y0), self.try_index(i, y1)) {
                if !self.cells[i0].wall && !self.cells[i1].wall {
                    self.cells[i0].velocity_x = self.cells[i1].velocity_x;
                }
            }
        }

        // Horizontal v extrapolation
        for j in 1..=n {
            let x0 = 0;
            let x1 = 1;
            if let (Some(i0), Some(i1)) = (self.try_index(x0, j), self.try_index(x1, j)) {
                if !self.cells[i0].wall && !self.cells[i1].wall {
                    self.cells[i0].velocity_y = self.cells[i1].velocity_y;
                }
            }

            let x0 = n + 1;
            let x1 = n;
            if let (Some(i0), Some(i1)) = (self.try_index(x0, j), self.try_index(x1, j)) {
                if !self.cells[i0].wall && !self.cells[i1].wall {
                    self.cells[i0].velocity_y = self.cells[i1].velocity_y;
                }
            }
        }
    }

    /// Initialize the grid with a given velocity and density
    pub fn cell_init(&mut self, line: usize, column: usize, vx: f32, vy: f32, density: f32) {
        if line <= (N + 1.0) as usize && column <= (N + 1.0) as usize {
            let idx = self.to_index(line, column);
            if !self.cells[idx].wall {
                self.cells[idx].velocity_x = vx;
                self.cells[idx].velocity_y = vy;
                self.cells[idx].density = density;
            } else {
                println!("Impossible to modify a wall!");
            }
        } else {
            println!("Error: indices out of bounds!");
        }
    }

    /// Initialize the wall character of a cell
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

    /// Initialize the velocity of a cell
    pub fn velocity_init(&mut self, line: usize, column: usize, vx: f32, vy: f32) {
        if line <= (N + 1.0) as usize && column <= (N + 1.0) as usize {
            let idx = self.to_index(line, column);
            if !self.cells[idx].wall {
                self.cells[idx].velocity_x = vx;
                self.cells[idx].velocity_y = vy;
            } else {
                println!("Impossible to modify a wall!");
            }
        } else {
            println!("Error: indices out of bounds!");
        }
    }

    /// Add a circular density and velocity source in the middle of the grid
    pub fn center_source(
        &mut self,
        radius: f32,
        density_value: f32,
        velocity_magnitude: f32,
        is_radial: bool,
    ) {
        // Process the center of the grid
        let center_x = (N / 2.0) as usize;
        let center_y = (N / 2.0) as usize;

        // Scan every cell in the inscribed square of the circle
        let r_int = radius.floor() as usize;
        let start_i = if center_x > r_int {
            center_x - r_int
        } else {
            1
        };
        let end_i = (center_x + r_int).min(N as usize);
        let start_j = if center_y > r_int {
            center_y - r_int
        } else {
            1
        };
        let end_j = (center_y + r_int).min(N as usize);

        for i in start_i..=end_i {
            for j in start_j..=end_j {
                // Compute the distance from the center
                let dx = i as f32 - center_x as f32;
                let dy = j as f32 - center_y as f32;
                let distance_squared = dx * dx + dy * dy;

                // If the current cell is not a wall
                if distance_squared <= radius * radius {
                    let idx = self.to_index(i, j);
                    if !self.cells[idx].wall {
                        // Add density
                        self.cells[idx].density += density_value;

                        // Compute velocity direction
                        if is_radial {
                            // Initialize radial velocity
                            if distance_squared > 0.0 {
                                let distance = distance_squared.sqrt();
                                let dir_x = dx / distance;
                                let dir_y = dy / distance;

                                // Add velocity in the direction of the center
                                let factor = 1.0 - (distance / radius); // Plus fort au centre
                                self.cells[idx].velocity_x += dir_x * velocity_magnitude * factor;
                                self.cells[idx].velocity_y += dir_y * velocity_magnitude * factor;
                            }
                        } else {
                            // Circular velocity (whirlwind)
                            if distance_squared > 0.0 {
                                let distance = distance_squared.sqrt();
                                let dir_x = -dy / distance;
                                let dir_y = dx / distance;

                                // Velocity increases with distance from the center up to a certain point
                                let factor = (distance / radius) * (1.0 - distance / radius) * 4.0;
                                self.cells[idx].velocity_x += dir_x * velocity_magnitude * factor;
                                self.cells[idx].velocity_y += dir_y * velocity_magnitude * factor;
                            }
                        }
                    }
                }
            }
        }
    }

    /// Print the grid with walls
    pub fn print_grid_walls(&self) {
        let lines: Vec<String> = (0..=(N + 1.0) as usize)
            .into_par_iter()
            .map(|j| {
                let mut line = String::with_capacity(((N + 2.0) as usize) * 3);
                for i in 0..=(N + 1.0) as usize {
                    let idx = self.to_index(i, j);
                    line.push_str(if self.cells[idx].wall {
                        " █ "
                    } else {
                        " · "
                    });
                }
                line.push('\n');
                line
            })
            .collect();

        for line in lines {
            print!("{line}");
        }
    }

    /// Print the grid with velocities
    pub fn print_grid_velocity(&self) {
        let lines: Vec<String> = (0..=(N + 1.0) as usize)
            .into_par_iter()
            .map(|j| {
                let mut line = String::with_capacity(((N + 2.0) as usize) * 18);
                for i in 0..=(N + 1.0) as usize {
                    let idx = self.to_index(i, j);
                    if self.cells[idx].wall {
                        line.push_str("   | █████ |   ");
                    } else {
                        let vx = &self.cells[idx].velocity_x;
                        let vy = &self.cells[idx].velocity_y;
                        line.push_str(&format!("|{:.3}, {:.3}| ", vx, vy));
                    }
                }
                line.push('\n');
                line
            })
            .collect();

        for line in lines {
            print!("{line}");
        }
    }

    /// Print the grid with densities
    pub fn print_grid_density(&self) {
        let output: String = (0..=(N + 1.0) as usize)
            .into_par_iter()
            .map(|j| {
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
            })
            .collect();

        print!("{output}");
    }

    /// Compute the total density of the grid
    pub fn total_density(&self) -> f32 {
        self.cells.iter().map(|cell| cell.density).sum()
    }

    /// Create a circle in the grid
    pub fn circle(&mut self, center_x: isize, center_y: isize, radius: f32) {
        let center_x = center_x as usize;
        let center_y = center_y as usize;
        for i in 0..=N as usize {
            for j in 0..=N as usize {
                let dx = i as isize - center_x as isize;
                let dy = j as isize - center_y as isize;
                if dx * dx + dy * dy <= (radius * radius) as isize {
                    self.wall_init(j, i, true);
                }
            }
        }
    }

    /// Create a wind tunnel like setup to test as in wind tunnel experiments
    pub fn initialize_wind_tunnel(&mut self, density: f32, hole_positions: &[usize]) {
        let left_wall = 1;
        let box_width = (N as usize / 30).max(2); // Width of the box
        let right_wall = left_wall + box_width;

        // Create the walls
        for j in 1..=N as usize {
            let idx_left = self.to_index(left_wall, j);
            if !self.cells[idx_left].wall {
                self.wall_init(j, left_wall, true); // Left wall
            }

            if !hole_positions.contains(&j) {
                let idx_right = self.to_index(right_wall, j);
                if !self.cells[idx_right].wall {
                    self.wall_init(j, right_wall, true); // Right wall, except for holes
                }
            }
        }

        // Add density in the box
        for i in (left_wall + 2)..=right_wall - 5 {
            for j in 2..=N as usize - 1 {
                let idx = self.to_index(i, j);
                if !self.cells[idx].wall {
                    // Ensure this is not a wall
                    self.cell_init(i, j, FLOW_VELOCITY, 0.0, density);
                }
            }
        }
    }

    /// Attempt to create karman vortex formation favorable conditions
    pub fn setup_karman_vortex(&mut self) {
        use crate::conditions::*;
        let radius = (N as usize / 6) as f32 - 10.5; // Obstacle radius
        let center = ((N as usize + 2) / 6, (N as usize + 2) / 2);

        self.circle(center.0 as isize, center.1 as isize, radius);

        // Horizontal flow
        for j in 1..(N as usize + 1) {
            for i in 1..15 {
                let idx = self.to_index(i, j);
                if !self.cells[idx].wall {
                    self.cells[idx].density = FLOW_DENSITY;
                    self.cells[idx].velocity_x = FLOW_VELOCITY
                }
            }
        }

        println!("✔ Kármán vortex setup complete: central obstacle + inflow.");
    }

    /// Perform a step in the simulation
    pub fn vel_step(&mut self) {
        //self.diffuse(VISCOSITY, DT);
        //println!("Total density after diffuse {:2}", self.total_density());
        if PROJECT == "1" {
            self.project();
        } else if PROJECT == "2" {
            //self.project2(80, DT, 1.9);
        }
        //println!("Total density after project {:2}", self.total_density());
        self.advect_velocity(DT);
        //println!("Total density after advect velocity {:2}", self.total_density());
        self.advect_density(DT);
        //println!("Total density after apres advect density {:2}", self.total_density());
        //self.apply_boundary_conditions();
    }

    /// Perform a step in the simulation with another method
    pub fn vel2_step(&mut self) {
        if PROJECT == "1" {
            self.project();
        } else if PROJECT == "2" {
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

    /// Perform a step in the simulation with cip-csl4 method (not working)
    pub fn vel_step_cip_csl4(&mut self) {
        // Diffusion only if needed
        // self.diffuse(VISCOSITY, DT);

        // Projection to ensure incompressibility
        if PROJECT == "1" {
            self.project();
        } else if PROJECT == "2" {
            //self.project2(80, DT, 1.9);
        }

        // Velocity advection
        self.advect_velocity(DT);

        // Use CIP-CSL4 method over other ones
        Self::advect_density_cip_csl4(self, DT);

        // Apply boundary conditions
        self.apply_boundary_conditions();
    }
}
