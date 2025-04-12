use crate::grid::{Grid};
use rayon::prelude::*;
use crate::conditions::*;

/// Effectue une étape d'advection CIP — CSL4 sur la densité.
///
/// Cette fonction réalise une interpolation tensorielle sur un stencil 5×5
/// afin d'obtenir une interpolation d'ordre 4 en chaque dimension.
///
/// Note : Cette version utilise une interpolation Lagrange classique et n'assure
/// pas encore la conservation stricte de la masse par intégration.
pub fn advect_cip_csl4(grid: &mut Grid, dt: f32) {
    // On considère que la taille de cellule est 1 (DX = 1.0)
    let grid_size = (N + 2.0) as usize; // taille globale (nombre de cellules en x = N+2, idem en y)

    // Pour chaque cellule de la grille, on calcule la nouvelle densité.
    let new_density: Vec<f32> = (0..grid.cells.len()).into_par_iter().map(|idx| {
        let (i, j) = grid.decode_index(idx);

        // Pour les murs, on conserve la densité
        if grid.cells[idx].wall {
            return grid.cells[idx].density;
        }

        // Backtracking : position amont
        let vel = grid.cells[idx].velocity;
        let x_back = i as f32 - vel.x * dt;
        let y_back = j as f32 - vel.y * dt;

        // On choisit un stencil de 5 points en chaque dimension.
        // Les nœuds théoriques seront aux positions -2, -1, 0, 1, 2.
        // On définit le centre du stencil comme floor(x_back) et floor(y_back).
        let i_center = x_back.floor() as isize;
        let j_center = y_back.floor() as isize;
        let frac_x = x_back - i_center as f32;
        let frac_y = y_back - j_center as f32;

        // Pour chaque des 5 lignes autour de j_center (de j_center-2 à j_center+2),
        // on réalise une interpolation 1D en x.
        let mut interp_rows = [0.0f32; 5];
        for m in 0..5 {
            let j_idx = j_center - 2 + m as isize;
            let mut fvals = [0.0f32; 5];
            for n in 0..5 {
                let i_idx = i_center - 2 + n as isize;
                // Vérifie les limites (si hors de la grille, on utilise la valeur du point courant)
                if i_idx < 0 || j_idx < 0 || (i_idx as usize) >= grid_size || (j_idx as usize) >= grid_size {
                    fvals[n] = grid.cells[idx].density;
                } else {
                    let neighbor_idx = grid.to_index(i_idx as usize, j_idx as usize);
                    // Si le voisin est un mur, on utilise la densité de la cellule courante
                    fvals[n] = if grid.cells[neighbor_idx].wall {
                        grid.cells[idx].density
                    } else {
                        grid.cells[neighbor_idx].density
                    };
                }
            }
            // Interpolation 1D en x sur le stencil de 5 points
            interp_rows[m] = cip_csl4_1d(frac_x, &fvals);
        }
        // Interpolation 1D en y sur les 5 lignes obtenues
        let final_density = cip_csl4_1d(frac_y, &interp_rows);
        final_density.max(0.0)
    }).collect();

    // Mise à jour de la grille (pour les cellules non-murs)
    for (idx, &density) in new_density.iter().enumerate() {
        if !grid.cells[idx].wall {
            grid.cells[idx].density = density;
        }
    }
}

/// Effectue une interpolation 1D d'ordre 4 (basée sur 5 points) par Lagrange.
///
/// - `x` est la coordonnée fractionnaire relative aux nœuds (les nœuds sont placés à -2, -1, 0, 1, 2).
/// - `f` est le tableau des valeurs aux 5 nœuds.
fn cip_csl4_1d(x: f32, f: &[f32; 5]) -> f32 {
    let nodes = [-2.0, -1.0, 0.0, 1.0, 2.0];
    let mut result = 0.0;
    for i in 0..5 {
        let mut L = 1.0;
        for j in 0..5 {
            if i != j {
                L *= (x - nodes[j]) / (nodes[i] - nodes[j]);
            }
        }
        result += f[i] * L;
    }
    result
}
