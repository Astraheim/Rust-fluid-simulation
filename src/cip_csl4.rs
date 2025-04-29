use crate::conditions::*;
use crate::grid::{Grid, Vector2};


/*
FR :
    Fonctionnalités clés de l'implémentation CIP-CSL4

    Calcul des gradients : La méthode nécessite les dérivées premières et secondes de la densité pour construire une interpolation précise.
    Conservation de la masse : Nous calculons les masses des cellules et les flux traversant les faces pour conserver la masse totale.
    Coefficients d'interpolation : Les coefficients a, b, c pour l'interpolation polynomiale d'ordre 4 sont calculés à partir des valeurs de densité, des gradients et de la masse.
    Traçage des caractéristiques : Nous suivons les particules fluides en remontant le temps (backtracking) pour déterminer d'où vient chaque valeur.
    Mise à jour des dérivées : Les dérivées sont également advectées pour maintenir la cohérence de l'ordre élevé du schéma.

    Améliorations possibles :

    L'implémentation actuelle pourrait être optimisée davantage en utilisant Rayon pour la parallélisation.
    Le calcul des gradients pourrait être amélioré avec des méthodes plus précises (WENO, ENO).
    La gestion des conditions aux limites pourrait être affinée pour mieux traiter les cas particuliers.
    Un limiteur de pente pourrait être ajouté pour éviter les oscillations dans les régions à fort gradient.


ENG :

    Key Features of the CIP-CSL4 Implementation

    Gradient Computation : The method requires first and second derivatives of the density to construct accurate interpolation.
    Mass Conservation : We compute the mass of each cell and the fluxes across cell faces to conserve total mass.
    Interpolation Coefficients : The coefficients a, b, c for the 4th-order polynomial interpolation are calculated from density values, gradients, and mass.
    Characteristic Tracing : Fluid particles are tracked backward in time (backtracking) to determine the origin of each value.
    Derivative Update : Derivatives are also advected to maintain the high-order consistency of the scheme.

    Possible Improvements :

    The current implementation could be further optimized using Rayon for parallelization.
    Gradient computation could be enhanced using more accurate methods (WENO, ENO).
    Boundary condition handling could be refined to better address edge cases.
    A slope limiter could be added to prevent oscillations in regions with strong gradients.


*/


impl Grid {
    // CIP-CSL4 method for density advection
    pub fn advect_density_cip_csl4(&mut self, dt: f32) {
        let dt0 = dt * N;

        // Temporary variables to store intermediate results
        let mut a_coeffs = vec![Vector2::default(); self.cells.len()];
        let mut b_coeffs = vec![Vector2::default(); self.cells.len()];
        let mut c_coeffs = vec![Vector2::default(); self.cells.len()];
        let mut dr = vec![0.0; self.cells.len()]; // Flux across cell faces
        let mut cell_masses = vec![0.0; self.cells.len()]; // Cell masses

        // Compute the initial masses of the cells
        // For density, mass is density * cell volume (here 1.0 / (N*N))
        let cell_volume = 1.0 / (N * N);

        // Initialize cell masses
        for idx in 0..self.cells.len() {
            cell_masses[idx] = self.cells[idx].density * cell_volume;
        }

        // Compute density gradients (1st order)
        self.compute_density_gradients();

        // Collect necessary data for coefficient computation
        let mut cells_to_process = Vec::new();
        for (i, j, idx) in self.iter_morton() {
            if !self.in_bounds(i, j) || self.cells[idx].wall {
                continue;
            }

            // Determine the "upwind" cell based on velocity direction
            let vx = self.cells[idx].velocity.x;
            let vy = self.cells[idx].velocity.y;

            let sign_x = if vx >= 0.0 { 1 } else { -1 };
            let sign_y = if vy >= 0.0 { 1 } else { -1 };

            let i_up = (i as isize - sign_x) as usize;
            let j_up = (j as isize - sign_y) as usize;

            // Check if the upwind cell is within bounds
            if !self.in_bounds(i_up, j_up) {
                continue;
            }

            let idx_up = self.to_index(i_up, j_up);
            if self.cells[idx_up].wall {
                continue;
            }

            cells_to_process.push((i, j, idx, i_up, j_up, idx_up));
        }

        // Compute coefficients and fluxes
        for &(i, j, idx, i_up, j_up, idx_up) in &cells_to_process {
            // Retrieve cell data
            let vx = self.cells[idx].velocity.x;
            let vy = self.cells[idx].velocity.y;
            let u_i = self.cells[idx].density;
            let u_up = self.cells[idx_up].density;
            let du_i_x = self.cells[idx].density_x;
            let du_i_y = self.cells[idx].density_y;
            let du_up_x = self.cells[idx_up].density_x;
            let du_up_y = self.cells[idx_up].density_y;
            let r_c = cell_masses[idx];

            // Spatial displacement (characteristic)
            let e_x = -vx * dt0;
            let e_y = -vy * dt0;

            // Square, cube, etc., of the displacement
            let e2_x = e_x * e_x;
            let e3_x = e_x * e2_x;
            let e4_x = e_x * e3_x;
            let e5_x = e_x * e4_x;

            let e2_y = e_y * e_y;
            let e3_y = e_y * e2_y;
            let e4_y = e_y * e3_y;
            let e5_y = e_y * e4_y;

            // Distance to the upwind cell
            let dx_i = 1.0 / N; // Distance between cells
            let dx5 = dx_i.powi(5);
            let dx4 = dx_i.powi(4);
            let dx3 = dx_i.powi(3);
            let dx2 = dx_i.powi(2);

            let v_sgn_x = if vx >= 0.0 { 1.0 } else { -1.0 };
            let v_sgn_y = if vy >= 0.0 { 1.0 } else { -1.0 };

            // Compute coefficients for the X direction
            a_coeffs[idx].x = (-2.5 / dx5) * (6.0 * (u_up + u_i) * dx_i
                - (du_up_x - du_i_x) * dx2
                + 12.0 * v_sgn_x * r_c);

            b_coeffs[idx].x = (4.0 / dx4) * ((7.0 * u_up + 8.0 * u_i) * dx_i
                - (du_up_x - 1.5 * du_i_x) * dx2
                + 15.0 * v_sgn_x * r_c);

            c_coeffs[idx].x = (-1.5 / dx3) * (4.0 * (2.0 * u_up + 3.0 * u_i) * dx_i
                - (du_up_x - 3.0 * du_i_x) * dx2
                + 20.0 * v_sgn_x * r_c);

            // Compute coefficients for the Y direction
            a_coeffs[idx].y = (-2.5 / dx5) * (6.0 * (u_up + u_i) * dx_i
                - (du_up_y - du_i_y) * dx2
                + 12.0 * v_sgn_y * r_c);

            b_coeffs[idx].y = (4.0 / dx4) * ((7.0 * u_up + 8.0 * u_i) * dx_i
                - (du_up_y - 1.5 * du_i_y) * dx2
                + 15.0 * v_sgn_y * r_c);

            c_coeffs[idx].y = (-1.5 / dx3) * (4.0 * (2.0 * u_up + 3.0 * u_i) * dx_i
                - (du_up_y - 3.0 * du_i_y) * dx2
                + 20.0 * v_sgn_y * r_c);

            // Compute flux across faces (X component)
            dr[idx] = -1.0 * (
                e5_x * a_coeffs[idx].x / 5.0 +
                    e4_x * b_coeffs[idx].x / 4.0 +
                    e3_x * c_coeffs[idx].x / 3.0 +
                    e2_x * du_i_x / 2.0 +
                    e_x * u_i
            );

            // Add Y component to the flux
            dr[idx] += -1.0 * (
                e5_y * a_coeffs[idx].y / 5.0 +
                    e4_y * b_coeffs[idx].y / 4.0 +
                    e3_y * c_coeffs[idx].y / 3.0 +
                    e2_y * du_i_y / 2.0 +
                    e_y * u_i
            );
        }

        // Collect information for mass updates
        let mut mass_updates = Vec::new();
        for (i, j, idx) in self.iter_morton() {
            if !self.in_bounds(i, j) || self.cells[idx].wall {
                continue;
            }

            // Compute indices of adjacent cells
            let idx_right = if i < N as usize { self.to_index(i + 1, j) } else { idx };
            let idx_top = if j < N as usize { self.to_index(i, j + 1) } else { idx };

            let mut delta_mass = 0.0;

            // Update mass using incoming and outgoing fluxes
            if idx_right < self.cells.len() && !self.cells[idx_right].wall {
                delta_mass += dr[idx] - dr[idx_right];
            }

            if idx_top < self.cells.len() && !self.cells[idx_top].wall {
                delta_mass += dr[idx] - dr[idx_top];
            }

            mass_updates.push((idx, delta_mass));
        }

        // Apply mass updates
        for (idx, delta_mass) in mass_updates {
            cell_masses[idx] += delta_mass;
        }

        // Collect data for density and derivative updates
        let mut density_updates = Vec::new();
        for (i, j, idx) in self.iter_morton() {
            if !self.in_bounds(i, j) || self.cells[idx].wall {
                continue;
            }

            // Spatial displacement (characteristic)
            let vx = self.cells[idx].velocity.x;
            let vy = self.cells[idx].velocity.y;
            let e_x = -vx * dt0;
            let e_y = -vy * dt0;

            // Powers of the displacement
            let e2_x = e_x * e_x;
            let e3_x = e_x * e2_x;

            let e2_y = e_y * e_y;
            let e3_y = e_y * e2_y;
            
            // Compute new derivatives
            let new_density_x = 4.0 * a_coeffs[idx].x * e3_x +
                3.0 * b_coeffs[idx].x * e2_x +
                2.0 * c_coeffs[idx].x * e_x +
                self.cells[idx].density_x;

            let new_density_y = 4.0 * a_coeffs[idx].y * e3_y +
                3.0 * b_coeffs[idx].y * e2_y +
                2.0 * c_coeffs[idx].y * e_y +
                self.cells[idx].density_y;

            // Convert mass to density
            let new_density = cell_masses[idx] / cell_volume;

            density_updates.push((idx, new_density, new_density_x, new_density_y));
        }

        // Apply density and derivative updates
        for (idx, new_density, new_density_x, new_density_y) in density_updates {
            self.cells[idx].density = new_density;
            self.cells[idx].density_x = new_density_x;
            self.cells[idx].density_y = new_density_y;
        }

        // Apply boundary conditions for density
        self.apply_density_boundary_conditions();
    }

    // Compute density gradients
    pub fn compute_density_gradients(&mut self) {
        // Use centered finite differences to compute derivatives
        let n = N as usize;
        let h = 1.0 / N;

        // Structure to store all updates
        struct GradientUpdate {
            idx: usize,
            density_x: f32,
            density_y: f32,
            density_xx: Option<f32>,
            density_yy: Option<f32>,
            density_xy: Option<f32>,
        }

        // Collect information for updates
        let mut updates = Vec::new();

        // First loop to compute values
        for (i, j, idx) in self.iter_morton() {
            if !self.in_bounds(i, j) || self.cells[idx].wall {
                continue;
            }

            let mut update = GradientUpdate {
                idx,
                density_x: 0.0,
                density_y: 0.0,
                density_xx: None,
                density_yy: None,
                density_xy: None,
            };

            // Retrieve the density of the current cell
            let density_current = self.cells[idx].density;

            // Centered finite differences for dx
            if i > 1 && i < n {
                let density_left = self.cells[self.to_index(i - 1, j)].density;
                let density_right = self.cells[self.to_index(i + 1, j)].density;
                update.density_x = (density_right - density_left) / (2.0 * h);
            } else if i == 1 {
                // Forward difference for the left edge
                let density_right = self.cells[self.to_index(i + 1, j)].density;
                update.density_x = (density_right - density_current) / h;
            } else if i == n {
                // Backward difference for the right edge
                let density_left = self.cells[self.to_index(i - 1, j)].density;
                update.density_x = (density_current - density_left) / h;
            }

            // Centered finite differences for dy
            if j > 1 && j < n {
                let density_bottom = self.cells[self.to_index(i, j - 1)].density;
                let density_top = self.cells[self.to_index(i, j + 1)].density;
                update.density_y = (density_top - density_bottom) / (2.0 * h);
            } else if j == 1 {
                // Forward difference for the bottom edge
                let density_top = self.cells[self.to_index(i, j + 1)].density;
                update.density_y = (density_top - density_current) / h;
            } else if j == n {
                // Backward difference for the top edge
                let density_bottom = self.cells[self.to_index(i, j - 1)].density;
                update.density_y = (density_current - density_bottom) / h;
            }

            // Compute second derivatives (for better accuracy)
            if i > 1 && i < n && j > 1 && j < n {
                let density_left = self.cells[self.to_index(i - 1, j)].density;
                let density_right = self.cells[self.to_index(i + 1, j)].density;
                let density_bottom = self.cells[self.to_index(i, j - 1)].density;
                let density_top = self.cells[self.to_index(i, j + 1)].density;

                // Second derivatives
                update.density_xx = Some((density_right - 2.0 * density_current + density_left) / (h * h));
                update.density_yy = Some((density_top - 2.0 * density_current + density_bottom) / (h * h));

                // Cross derivative
                let density_top_right = self.cells[self.to_index(i + 1, j + 1)].density;
                let density_top_left = self.cells[self.to_index(i - 1, j + 1)].density;
                let density_bottom_right = self.cells[self.to_index(i + 1, j - 1)].density;
                let density_bottom_left = self.cells[self.to_index(i - 1, j - 1)].density;

                update.density_xy = Some((density_top_right - density_top_left - density_bottom_right + density_bottom_left) / (4.0 * h * h));
            }

            updates.push(update);
        }

        // Apply all updates
        for update in updates {
            self.cells[update.idx].density_x = update.density_x;
            self.cells[update.idx].density_y = update.density_y;

            if let Some(density_xx) = update.density_xx {
                self.cells[update.idx].density_xx = density_xx;
            }

            if let Some(density_yy) = update.density_yy {
                self.cells[update.idx].density_yy = density_yy;
            }

            if let Some(density_xy) = update.density_xy {
                self.cells[update.idx].density_xy = density_xy;
            }
        }
    }

    // Apply boundary conditions for density
    pub fn apply_density_boundary_conditions(&mut self) {
        let n = N as usize;

        // Structure to store updates
        struct BoundaryUpdate {
            idx: usize,
            density: f32,
            density_x: Option<f32>,
            density_y: Option<f32>,
        }

        let mut updates = Vec::new();

        // Boundary conditions for horizontal edges
        for i in 1..=n {
            // Bottom edge
            let idx_bottom = self.to_index(i, 1);
            let idx_above = self.to_index(i, 2);
            if !self.cells[idx_bottom].wall {
                updates.push(BoundaryUpdate {
                    idx: idx_bottom,
                    density: self.cells[idx_above].density,
                    density_x: None,
                    density_y: Some(self.cells[idx_above].density_y),
                });
            }

            // Top edge
            let idx_top = self.to_index(i, n);
            let idx_below = self.to_index(i, n - 1);
            if !self.cells[idx_top].wall {
                updates.push(BoundaryUpdate {
                    idx: idx_top,
                    density: self.cells[idx_below].density,
                    density_x: None,
                    density_y: Some(self.cells[idx_below].density_y),
                });
            }
        }

        // Boundary conditions for vertical edges
        for j in 1..=n {
            // Left edge - inlet condition
            let idx_left = self.to_index(1, j);
            if !self.cells[idx_left].wall {
                // Constant density at the inlet
                updates.push(BoundaryUpdate {
                    idx: idx_left,
                    density: FLOW_DENSITY,
                    density_x: Some(0.0),
                    density_y: None,
                });
            }

            // Right edge - outlet condition
            let idx_right = self.to_index(n, j);
            let idx_before = self.to_index(n - 1, j);
            if !self.cells[idx_right].wall {
                updates.push(BoundaryUpdate {
                    idx: idx_right,
                    density: self.cells[idx_before].density,
                    density_x: Some(self.cells[idx_before].density_x),
                    density_y: None,
                });
            }
        }

        // Collect information for wall conditions
        let mut wall_updates = Vec::new();
        for (i, j, idx) in self.iter_morton() {
            if !self.in_bounds(i, j) || !self.cells[idx].wall {
                continue;
            }

            // For each wall, look at neighboring non-wall cells
            let mut sum_density = 0.0;
            let mut count = 0;

            for &(di, dj) in &[(-1, 0), (1, 0), (0, -1), (0, 1)] {
                let ni = (i as isize + di) as usize;
                let nj = (j as isize + dj) as usize;

                if self.in_bounds(ni, nj) {
                    let n_idx = self.to_index(ni, nj);
                    if !self.cells[n_idx].wall {
                        sum_density += self.cells[n_idx].density;
                        count += 1;
                    }
                }
            }

            // Reflect density for walls
            if count > 0 {
                wall_updates.push(BoundaryUpdate {
                    idx,
                    density: sum_density / count as f32,
                    density_x: Some(0.0),
                    density_y: Some(0.0),
                });
            }
        }

        // Add wall updates to the global list
        updates.extend(wall_updates);

        // Apply all updates
        for update in updates {
            self.cells[update.idx].density = update.density;

            if let Some(density_x) = update.density_x {
                self.cells[update.idx].density_x = density_x;
            }

            if let Some(density_y) = update.density_y {
                self.cells[update.idx].density_y = density_y;
            }
        }
    }
}
