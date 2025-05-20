//use crate::conditions::*;
use crate::grid2::{Grid2, Vector22,Cell2Type};
// Last update to rendu_code_tex: 2025-05-02
// last modif: 2025-05-13

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



// Structure pour représenter les forces sur un objet solide
#[derive(Clone, Debug)]
pub struct ObjectForce {
    pub id: usize,               // Identifiant unique de l'objet
    pub center_of_mass: Vector22, // Centre de masse de l'objet
    pub total_force: Vector22,    // Force totale appliquée sur l'objet
    pub torque: f32,             // Moment de force (couple)
    pub cell_count: usize,       // Nombre de cellules composant l'objet
}

impl Grid2 {
    /// Calcule les forces de pression causées par le fluide sur les murs de la grille
    pub fn compute_wall_forces(&self) -> Vec<Vector22> {
        let h = self.cell_size; // Taille d'une cellule
        let mut forces = vec![Vector22::zeros(); self.cells.len()];

        for (i, j, idx) in self.iter_morton() {
            if self.cells[idx].cell_type != Cell2Type::Solid {
                continue;
            }

            // Directions autour des cellules solides
            let neighbors = [
                (i.wrapping_sub(1), j, Vector22::new(-1.0, 0.0)), // gauche
                (i + 1, j, Vector22::new(1.0, 0.0)),             // droite
                (i, j.wrapping_sub(1), Vector22::new(0.0, -1.0)), // bas
                (i, j + 1, Vector22::new(0.0, 1.0)),             // haut
            ];

            for &(ni, nj, normal) in &neighbors {
                if let Some(nidx) = self.try_idx(ni, nj) {
                    if self.cells[nidx].cell_type == Cell2Type::Fluid {
                        // Force de pression
                        let p = self.cells[nidx].pressure;
                        let v = &self.cells[nidx].velocity.back;
                        let v_normal = v.x * normal.x + v.y * normal.y;

                        // Vitesse tangentielle
                        let v_tang_x = v.x - v_normal * normal.x;
                        let v_tang_y = v.y - v_normal * normal.y;

                        // Force de pression (normale au mur)
                        let f_pressure = Vector22::new(-p * normal.x, -p * normal.y);

                        // Force de traînée (opposée à la vitesse)
                        let drag_coef = 0.5; // Coefficient à ajuster
                        let f_drag_x = drag_coef * v_normal.abs() * v_normal * normal.x;
                        let f_drag_y = drag_coef * v_normal.abs() * v_normal * normal.y;

                        // Force de cisaillement (shear)
                        // On suppose VISCOSITY = 0.1 car non défini dans le code fourni
                        let visc_coef = 0.1;
                        let f_shear_x = visc_coef * v_tang_x;
                        let f_shear_y = visc_coef * v_tang_y;

                        // Force totale
                        forces[idx].x += (f_pressure.x + f_drag_x + f_shear_x) * h;
                        forces[idx].y += (f_pressure.y + f_drag_y + f_shear_y) * h;
                    }
                }
            }
        }

        forces
    }

    /// Identification des objets solides connectés
    pub fn identify_objects(&self) -> Vec<usize> {
        let mut object_ids = vec![0; self.cells.len()];
        let mut current_id = 1;
        let mut stack = Vec::new();

        for (i, j, idx) in self.iter_morton() {
            if self.cells[idx].cell_type == Cell2Type::Solid && object_ids[idx] == 0 {
                // Nouvel objet trouvé
                object_ids[idx] = current_id;
                stack.push((i, j));

                // Exploration de l'objet par parcours en profondeur (DFS)
                while let Some((ci, cj)) = stack.pop() {
                    let neighbors = [
                        (ci.wrapping_sub(1), cj), // gauche
                        (ci + 1, cj),            // droite
                        (ci, cj.wrapping_sub(1)), // bas
                        (ci, cj + 1),            // haut
                    ];

                    for &(ni, nj) in &neighbors {
                        if let Some(nidx) = self.try_idx(ni, nj) {
                            if self.cells[nidx].cell_type == Cell2Type::Solid && object_ids[nidx] == 0 {
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

    /// Calcule les forces totales sur chaque objet
    pub fn compute_object_forces(&self) -> Vec<ObjectForce> {
        let cell_forces = self.compute_wall_forces();
        let object_ids = self.identify_objects();
        let max_id = *object_ids.iter().max().unwrap_or(&0);

        if max_id == 0 {
            return Vec::new();
        }

        // Initialisation des accumulateurs d'objets
        let mut objects = Vec::with_capacity(max_id);
        for id in 1..=max_id {
            objects.push(ObjectForce {
                id,
                center_of_mass: Vector22::zeros(),
                total_force: Vector22::zeros(),
                torque: 0.0,
                cell_count: 0,
            });
        }

        // Pour chaque cellule solide, on accumule les forces et le centre de masse
        for (i, j, idx) in self.iter_morton() {
            if self.cells[idx].cell_type == Cell2Type::Solid {
                let obj_id = object_ids[idx];
                if obj_id > 0 {
                    let obj_idx = obj_id - 1;

                    // Accumulation du centre de masse
                    objects[obj_idx].center_of_mass.x += i as f32;
                    objects[obj_idx].center_of_mass.y += j as f32;
                    objects[obj_idx].cell_count += 1;

                    // Accumulation des forces
                    objects[obj_idx].total_force.x += cell_forces[idx].x;
                    objects[obj_idx].total_force.y += cell_forces[idx].y;
                }
            }
        }

        // Calcul final du centre de masse et du couple
        for obj in &mut objects {
            if obj.cell_count > 0 {
                // Traitement du centre de masse
                obj.center_of_mass.x /= obj.cell_count as f32;
                obj.center_of_mass.y /= obj.cell_count as f32;

                // Recalcul du couple (torque)
                for (i, j, idx) in self.iter_morton() {
                    if self.cells[idx].cell_type == Cell2Type::Solid && object_ids[idx] == obj.id {
                        // Vecteur du centre de masse à la cellule
                        let r_x = i as f32 - obj.center_of_mass.x;
                        let r_y = j as f32 - obj.center_of_mass.y;

                        // Produit vectoriel 2D : r × F = r_x*F_y - r_y*F_x
                        obj.torque += r_x * cell_forces[idx].y - r_y * cell_forces[idx].x;
                    }
                }
            }
        }

        objects
    }

    /// Calcule la magnitude d'un vecteur Vector22
    pub fn magnitude(v: &Vector22) -> f32 {
        (v.x * v.x + v.y * v.y).sqrt()
    }

    /// Normalise un vecteur Vector22
    pub fn normalize(v: &Vector22) -> Vector22 {
        let mag = Self::magnitude(v);
        if mag > 1e-6 {
            Vector22 {
                x: v.x / mag,
                y: v.y / mag,
            }
        } else {
            Vector22::zeros()
        }
    }

    /// Identifie les objets, calcule les forces et les affiche
    pub fn print_object_forces(&self) {
        let objects = self.compute_object_forces();

        if objects.is_empty() {
            println!("Aucun objet détecté.");
            return;
        }

        println!("\n=== Forces sur les objets ===");
        for obj in &objects {
            let force_magnitude = Self::magnitude(&obj.total_force);

            println!("Objet #{} :", obj.id);
            println!("  - Nombre de cellules: {}", obj.cell_count);
            println!("  - Centre de masse: ({:.2}, {:.2})", obj.center_of_mass.x, obj.center_of_mass.y);
            println!("  - Force totale: ({:.4}, {:.4}) [magnitude: {:.4}]",
                     obj.total_force.x, obj.total_force.y, force_magnitude);
            println!("  - Moment de force: {:.4}", obj.torque);

            // Direction principale de la force
            if force_magnitude > 0.001 {
                let direction = Self::normalize(&obj.total_force);
                println!("  - Direction de la force: ({:.2}, {:.2})", direction.x, direction.y);

                // Approximation de l'angle (à des fins de visualisation uniquement)
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