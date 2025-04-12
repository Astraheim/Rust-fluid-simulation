use rayon::prelude::*;

/// Retourne l'inverse pré-calculé de la matrice A pour un stencil 5 points en 1D.
/// (Même fonction que dans la version 1D.)
fn inv_A_const() -> [[f32; 5]; 5] {
    // Les valeurs ici doivent correspondre aux moments I_n(k) pour k = -2,...,2,
    // et leur inversion. Les valeurs présentées précédemment sont un exemple.
    let A = [
        [1.0,   -2.0,    4.08333, -8.5,    18.0125],
        [1.0,   -1.0,    1.08333, -1.25,   1.5125],
        [1.0,    0.0,    0.08333,  0.0,    0.0125],
        [1.0,    1.0,    1.08333,  1.25,   1.5125],
        [1.0,    2.0,    4.08333,  8.5,    18.0125],
    ];
    invert_matrix_5x5(A)
}

/// Inverse d'une matrice 5x5 par Pivot de Gauss.
fn invert_matrix_5x5(A: [[f32; 5]; 5]) -> [[f32; 5]; 5] {
    let mut M = A;
    let mut inv = [[0.0f32; 5]; 5];
    for i in 0..5 { inv[i][i] = 1.0; }
    for i in 0..5 {
        let pivot = M[i][i];
        for j in 0..5 {
            M[i][j] /= pivot;
            inv[i][j] /= pivot;
        }
        for k in 0..5 {
            if k != i {
                let factor = M[k][i];
                for j in 0..5 {
                    M[k][j] -= factor * M[i][j];
                    inv[k][j] -= factor * inv[i][j];
                }
            }
        }
    }
    inv
}

use rayon::prelude::*;
use crate::{Grid, Cell};

pub fn advect_cip_csl4_2d(grid: &mut Grid, dt: f32, dx: f32, dy: f32, grid_size: usize) {
    let inv_A = inv_A_const();
    let mut new_cells = grid.cells.clone();

    for j in 2..(grid_size - 2) {
        for i in 2..(grid_size - 2) {
            let idx = i + j * grid_size;
            let cell = grid.cells[idx];

            if cell.wall {
                continue;
            }

            let alpha_x = cell.velocity.x * dt / dx;
            let alpha_y = cell.velocity.y * dt / dy;

            let x_back = i as f32 - alpha_x;
            let y_back = j as f32 - alpha_y;
            let i_center = x_back.floor() as isize;
            let j_center = y_back.floor() as isize;
            let frac_x = x_back - i_center as f32;
            let frac_y = y_back - j_center as f32;

            if i_center < 2 || j_center < 2 || (i_center + 2) as usize >= grid_size || (j_center + 2) as usize >= grid_size {
                continue;
            }

            let mut b = [[0.0f32; 5]; 5];
            let mut skip = false;
            for r in 0..5 {
                for s in 0..5 {
                    let ii = (i_center - 2 + s as isize) as usize;
                    let jj = (j_center - 2 + r as isize) as usize;
                    let n_idx = ii + jj * grid_size;
                    let neighbor = grid.cells[n_idx];
                    if neighbor.wall {
                        skip = true;
                        break;
                    }
                    b[r][s] = neighbor.density;
                }
                if skip { break; }
            }
            if skip { continue; }

            let mut a = [[0.0f32; 5]; 5];
            for p in 0..5 {
                for q in 0..5 {
                    for r in 0..5 {
                        for s in 0..5 {
                            a[p][q] += inv_A[r][p] * inv_A[s][q] * b[r][s];
                        }
                    }
                }
            }

            let F = |x: f32, y: f32| -> f32 {
                let mut sum = 0.0;
                for p in 0..5 {
                    for q in 0..5 {
                        sum += a[p][q] * x.powi((p + 1) as i32) / (p as f32 + 1.0)
                            * y.powi((q + 1) as i32) / (q as f32 + 1.0);
                    }
                }
                sum
            };

            let flux = F(0.5 + frac_x, 0.5 + frac_y)
                - F(-0.5 + frac_x, 0.5 + frac_y)
                - F(0.5 + frac_x, -0.5 + frac_y)
                + F(-0.5 + frac_x, -0.5 + frac_y);

            let mut local_min = b[0][0];
            let mut local_max = b[0][0];
            for r in 0..5 {
                for s in 0..5 {
                    local_min = local_min.min(b[r][s]);
                    local_max = local_max.max(b[r][s]);
                }
            }

            new_cells[idx].density = flux.clamp(local_min, local_max);
        }
    }

    grid.cells.copy_from_slice(&new_cells);
}




fn monomial_basis_vector_4(x: f32, y: f32) -> [f32; 16] {
    let x_powers = [1.0, x, x * x, x * x * x];
    let y_powers = [1.0, y, y * y, y * y * y];

    let mut basis = [0.0; 16];
    for i in 0..4 {
        for j in 0..4 {
            basis[i * 4 + j] = x_powers[i] * y_powers[j];
        }
    }
    basis
}



fn monomial_basis_4(x: f32, y: f32) -> [[f32; 4]; 4] {
    let x_powers = [1.0, x, x * x, x * x * x];
    let y_powers = [1.0, y, y * y, y * y * y];

    let mut basis = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            basis[i][j] = x_powers[i] * y_powers[j];
        }
    }
    basis
}


fn evaluate_polynomial(coeffs: &[f32; 16], x: f32, y: f32) -> f32 {
    let monomials = monomial_basis_vector_4(x, y);
    monomials.iter().zip(coeffs.iter()).map(|(m, c)| m * c).sum()
}
