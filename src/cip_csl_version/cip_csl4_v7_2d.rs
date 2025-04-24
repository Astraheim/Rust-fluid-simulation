use rayon::prelude::*;
use crate::{Grid, Cell};

/// Construit la base monomiale cubique sur 4x4 points centrés.
fn build_monomial_matrix() -> [[f32; 16]; 16] {
    let mut mat = [[0.0; 16]; 16];
    let mut row = 0;
    for j in -1..=2 {
        for i in -1..=2 {
            let mut col = 0;
            for q in 0..4 {
                for p in 0..4 {
                    mat[row][col] = (i as f32).powi(p) * (j as f32).powi(q);
                    col += 1;
                }
            }
            row += 1;
        }
    }
    mat
}

/// Inversion naïve d'une matrice 16×16
fn invert_matrix_16x16(matrix: [[f32; 16]; 16]) -> [[f32; 16]; 16] {
    let mut a = matrix;
    let mut inv = [[0.0; 16]; 16];
    for i in 0..16 { inv[i][i] = 1.0; }

    for i in 0..16 {
        let pivot = a[i][i];
        for j in 0..16 {
            a[i][j] /= pivot;
            inv[i][j] /= pivot;
        }
        for k in 0..16 {
            if k != i {
                let factor = a[k][i];
                for j in 0..16 {
                    a[k][j] -= factor * a[i][j];
                    inv[k][j] -= factor * inv[i][j];
                }
            }
        }
    }

    inv
}

/// Produit scalaire entre deux tableaux de 16 éléments.
fn dot16(a: &[f32; 16], b: &[f32; 16]) -> f32 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Retourne le vecteur de base monomiale à un point (x, y)
fn monomial_basis_vector(x: f32, y: f32) -> [f32; 16] {
    let mut basis = [0.0; 16];
    let mut col = 0;
    for q in 0..4 {
        for p in 0..4 {
            basis[col] = x.powi(p) * y.powi(q);
            col += 1;
        }
    }
    basis
}

/// Évalue le polynôme 2D en un point à partir des coefficients
fn eval_poly(coeffs: &[f32; 16], x: f32, y: f32) -> f32 {
    let basis = monomial_basis_vector(x, y);
    dot16(coeffs, &basis)
}

/// Intégrale du polynôme sur [-1/2, 1/2]^2 pour avoir un flux conservatif
fn integrate_poly(coeffs: &[f32; 16]) -> f32 {
    let mut sum = 0.0;
    for q in 0..4 {
        for p in 0..4 {
            let idx = q * 4 + p;
            sum += coeffs[idx] / ((p + 1) as f32 * (q + 1) as f32)
                * (0.5f32).powi((p + q + 2) as i32) * 4.0;
        }
    }
    sum
}

/// Advection conservatrice CIP–CSL4 avec reconstruction bicubique
pub fn advect_cip_csl4_2d(grid: &mut Grid, dt: f32, dx: f32, dy: f32, grid_size: usize) {
    let monomial_matrix = build_monomial_matrix();
    let inv_basis = invert_matrix_16x16(monomial_matrix);

    let new_cells: Vec<Cell> = (0..grid_size * grid_size)
        .into_par_iter()
        .map(|flat_idx| {
            let (i, j) = grid.decode_index(flat_idx);
            let mut out_cell = grid.cells[flat_idx]; // Copie initiale

            if i < 2 || i >= grid_size - 2 || j < 2 || j >= grid_size - 2 || out_cell.wall {
                return out_cell;
            }

            let idx = flat_idx;
            let cell = grid.cells[idx];
            let alpha_x = cell.velocity.x * dt / dx;
            let alpha_y = cell.velocity.y * dt / dy;

            let x_back = (i as f32 + 0.5) - alpha_x;
            let y_back = (j as f32 + 0.5) - alpha_y;

            let i0 = x_back.floor() as isize;
            let j0 = y_back.floor() as isize;
            let frac_x = x_back - i0 as f32;
            let frac_y = y_back - j0 as f32;
            let local_x = x_back - i0 as f32;
            let local_y = y_back - j0 as f32;

            if i0 < 1 || j0 < 1 || i0 + 2 >= (grid_size - 1) as isize || j0 + 2 >= (grid_size - 1) as isize {
                return out_cell;
            }

            let mut stencil = [0.0f32; 16];
            let mut k = 0;
            for dj in -1..=2 {
                for di in -1..=2 {
                    let ni = (i0 + di) as usize;
                    let nj = (j0 + dj) as usize;
                    let nidx = grid.to_index(ni, nj);
                    let neighbor = grid.cells[nidx];
                    if neighbor.wall {
                        return out_cell;
                    }
                    stencil[k] = neighbor.density;
                    k += 1;
                }
            }

            let mut coeffs = [0.0f32; 16];
            for row in 0..16 {
                coeffs[row] = dot16(&inv_basis[row], &stencil);
            }

            let flux = eval_poly(&coeffs, local_x, local_y);

            let center_vals = [
                stencil[5], stencil[6], // ligne 2
                stencil[9], stencil[10] // ligne 3
            ];
            let (min_val, max_val) = center_vals.iter().fold((center_vals[0], center_vals[0]), |(min, max), &v| {
                (min.min(v), max.max(v))
            });
            out_cell.density = flux.clamp(min_val, max_val);
            out_cell
        })
        .collect();

    grid.cells = new_cells;
}
