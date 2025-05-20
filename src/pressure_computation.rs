use crate::conditions::*;
use crate::grid::{Grid, Vector2, ObjectForce};
// Last update to rendu_code_tex: 2025-05-02
// last modif: 2025-05-02

/*
FR :
    Fonctionnalités clés de l'analyse des forces sur objets solides (murs)

    Détection des objets : Les cellules marquées comme "wall" sont regroupées en objets connexes via une exploration en profondeur (DFS).
    Calcul des forces murales : Pour chaque cellule murale, on évalue les forces dues à la pression du fluide, la traînée (drag) et le cisaillement (shear).
    Agrégation par objet : Les forces sont ensuite regroupées par objet pour obtenir la force totale, le centre de masse et le moment de force (torque).
    Analyse directionnelle : L'orientation principale de la force résultante est estimée, avec une indication directionnelle (ex : ↑, ↘, ←).
    Visualisation : Les résultats sont affichés dans la console pour aider au diagnostic des interactions fluide-objet.

    Améliorations possibles :

    L'accumulation du couple (torque) pourrait tenir compte de la distance physique réelle plutôt que des indices de cellule.
    Une visualisation graphique des forces serait utile pour le debugging et la validation.


ENG :
    Key Features of Solid Object Force Analysis (walls)

    Object Detection : Cells marked as "wall" are grouped into connected objects using depth first search (DFS).
    Wall Force Computation : Each wall cell is subject to pressure, drag, and shear forces from adjacent fluid cells.
    Object-wise Aggregation : Forces are aggregated per object to compute total force, center of mass, and torque.
    Directional Analysis : The main direction of the resulting force is estimated and reported (e.g., ↑, ↘, ←).
    Visualization : Results are printed to the console to support debugging of fluid object interactions.

    Possible Improvements :

    Torque accumulation could use actual physical distances instead of cell indices.
    Graphical visualization of the forces would help with debugging and validation.
*/



impl Grid {


    /// Process the forces of pressure caused by fluid on the walls of the grid
    pub fn compute_wall_forces(&self) -> Vec<Vector2> {
        let h = 1.0 / N; // Size of a cell
        let mut forces = vec![Vector2::default(); self.cells.len()];

        for (i, j, idx) in self.iter_morton() {
            if !self.cells[idx].wall {
                continue;
            }

            // Directions around wall cells
            let neighbors = [
                (i.wrapping_sub(1), j, Vector2 { x: -1.0, y: 0.0 }), // left
                (i + 1, j, Vector2 { x: 1.0, y: 0.0 }),              // right
                (i, j.wrapping_sub(1), Vector2 { x: 0.0, y: -1.0 }), // bottom
                (i, j + 1, Vector2 { x: 0.0, y: 1.0 }),              // top
            ];

            for &(ni, nj, normal) in &neighbors {
                if let Some(nidx) = self.try_index(ni, nj) {
                    if !self.cells[nidx].wall {
                        // Pressure force
                        let p = self.cells[nidx].pressure;
                        let vx = &self.cells[nidx].velocity_x;
                        let vy = &self.cells[nidx].velocity_y;
                        let v_normal = vx * normal.x + vy * normal.y;
                        // Tangential velocity
                        let v_tang_x = vx - v_normal * normal.x;
                        let v_tang_y = vy - v_normal * normal.y;

                        // Pression force (normal to the wall)
                        let f_pressure = -p * normal;

                        // Drag force (opposite to the velocity)
                        let drag_coef = 0.5; // Coefficient à ajuster
                        let f_drag_x = drag_coef * v_normal.abs() * v_normal * normal.x;
                        let f_drag_y = drag_coef * v_normal.abs() * v_normal * normal.y;

                        // Wind shear force
                        let visc_coef = VISCOSITY;
                        let f_shear_x = visc_coef * v_tang_x;
                        let f_shear_y = visc_coef * v_tang_y;

                        // Total force
                        forces[idx].x += (f_pressure.x + f_drag_x + f_shear_x) * h;
                        forces[idx].y += (f_pressure.y + f_drag_y + f_shear_y) * h;
                    }
                }
            }
        }

        forces
    }

    /// Objects identification
    pub fn identify_objects(&self) -> Vec<usize> {
        let mut object_ids = vec![0; self.cells.len()];
        let mut current_id = 1;
        let mut stack = Vec::new();

        for (i, j, idx) in self.iter_morton() {
            if self.cells[idx].wall && object_ids[idx] == 0 {
                // New object found
                object_ids[idx] = current_id;
                stack.push((i, j));

                // Exploration of the object using Depth-First Search (DFS)
                while let Some((ci, cj)) = stack.pop() {
                    let neighbors = [
                        (ci.wrapping_sub(1), cj), // left
                        (ci + 1, cj),             // right
                        (ci, cj.wrapping_sub(1)), // down
                        (ci, cj + 1),             // top
                    ];

                    for &(ni, nj) in &neighbors {
                        if let Some(nidx) = self.try_index(ni, nj) {
                            if self.cells[nidx].wall && object_ids[nidx] == 0 {
                                object_ids[nidx] = current_id;
                                stack.push((ni, nj));
                            }
                        }
                    }
                }

                current_id += 1;
            }
        }

        object_ids
    }

    /// Computes total forces on each object
    pub fn compute_object_forces(&self) -> Vec<ObjectForce> {
        let cell_forces = self.compute_wall_forces();
        let object_ids = self.identify_objects();
        let max_id = *object_ids.iter().max().unwrap_or(&0);

        if max_id == 0 {
            return Vec::new();
        }

        // Initialization of objects accumulators
        let mut objects = Vec::with_capacity(max_id);
        for id in 1..=max_id {
            objects.push(ObjectForce {
                id,
                center_of_mass: Vector2::default(),
                total_force: Vector2::default(),
                torque: 0.0,
                cell_count: 0,
            });
        }

        // For each wall cell, we accumulate the forces and the center of mass
        for (i, j, idx) in self.iter_morton() {
            if self.cells[idx].wall {
                let obj_id = object_ids[idx];
                if obj_id > 0 {
                    let obj_idx = obj_id - 1;

                    // Accumulation in a center of mass
                    objects[obj_idx].center_of_mass.x += i as f32;
                    objects[obj_idx].center_of_mass.y += j as f32;
                    objects[obj_idx].cell_count += 1;

                    // Forces accumulation
                    objects[obj_idx].total_force.x += cell_forces[idx].x;
                    objects[obj_idx].total_force.y += cell_forces[idx].y;
                }
            }
        }

        // Final computation of the center of mass and torque
        for obj in &mut objects {
            if obj.cell_count > 0 {
                // Process center of mass
                obj.center_of_mass.x /= obj.cell_count as f32;
                obj.center_of_mass.y /= obj.cell_count as f32;

                // Recalculating torque
                for (i, j, idx) in self.iter_morton() {
                    if self.cells[idx].wall && object_ids[idx] == obj.id {
                        // Vecteur du centre de masse à la cellule
                        let r_x = i as f32 - obj.center_of_mass.x;
                        let r_y = j as f32 - obj.center_of_mass.y;

                        // 2D Vectorial product : r × F = r_x*F_y - r_y*F_x
                        obj.torque += r_x * cell_forces[idx].y - r_y * cell_forces[idx].x;
                    }
                }
            }
        }

        objects
    }




    /// Identifies objects, compute forces and print them
    pub fn print_object_forces(&self) {
        let objects = self.compute_object_forces();

        if objects.is_empty() {
            println!("Aucun objet détecté.");
            return;
        }

        println!("\n=== Forces sur les objets ===");
        for obj in &objects {
            let force_magnitude = obj.total_force.magnitude();

            println!("Objet #{} :", obj.id);
            println!("  - Nombre de cellules: {}", obj.cell_count);
            println!("  - Centre de masse: ({:.2}, {:.2})", obj.center_of_mass.x, obj.center_of_mass.y);
            println!("  - Force totale: ({:.4}, {:.4}) [magnitude: {:.4}]",
                     obj.total_force.x, obj.total_force.y, force_magnitude);
            println!("  - Moment de force: {:.4}", obj.torque);

            // Principal direction of the force
            if force_magnitude > 0.001 {
                let direction = obj.total_force.normalize();
                println!("  - Direction de la force: ({:.2}, {:.2})", direction.x, direction.y);

                // Approximation of the angle (only for visualization purposes)
                let angle = direction.y.atan2(direction.x) * 180.0 / std::f32::consts::PI;
                let direction_desc = match angle {
                    a if a > -22.5 && a <= 22.5 => "→ (droite)",
                    a if a > 22.5 && a <= 67.5 => "↗ (haut-droite)",
                    a if a > 67.5 && a <= 112.5 => "↑ (haut)",
                    a if a > 112.5 && a <= 157.5 => "↖ (haut-gauche)",
                    a if a > 157.5 || a <= -157.5 => "← (gauche)",
                    a if a > -157.5 && a <= -112.5 => "↙ (bas-gauche)",
                    a if a > -112.5 && a <= -67.5 => "↓ (bas)",
                    _ => "↘ (bas-droite)",
                };
                println!("  - Orientation: {} ({:.1}°)", direction_desc, angle);
            } else {
                println!("  - Force négligeable");
            }
            println!();
        }
        println!("===================\n \n");
    }



}