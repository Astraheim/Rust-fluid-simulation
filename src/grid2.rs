use crate::conditions::*;
use rayon::prelude::*;

#[derive(Clone, Copy, Debug, Default)]
pub struct Vector22 {
    pub x: f32,
    pub y: f32,
}

impl Vector22 {
    pub fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }

    pub fn zeros() -> Self {
        Self { x: 0.0, y: 0.0 }
    }

    // Méthodes utiles pour manipulation vectorielle
    pub fn add(&self, other: &Vector22) -> Vector22 {
        Vector22 {
            x: self.x + other.x,
            y: self.y + other.y
        }
    }

    pub fn scale(&self, factor: f32) -> Vector22 {
        Vector22 {
            x: self.x * factor,
            y: self.y * factor
        }
    }

    pub fn norm_squared(&self) -> f32 {
        self.x * self.x + self.y * self.y
    }
}

// Double buffer pour les champs qui nécessitent un échange front/back
#[derive(Clone, Copy, Debug, Default)]
pub struct DoubleBuffer<T: Copy + Default> {
    pub front: T,
    pub back: T,
}

impl<T: Copy + Default> DoubleBuffer<T> {
    pub fn new(val: T) -> Self {
        Self {
            front: val,
            back: val
        }
    }

    pub fn swap(&mut self) {
        std::mem::swap(&mut self.front, &mut self.back);
    }
}

// Type de cellule
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Cell2Type {
    Fluid,
    Solid,
}

impl Default for Cell2Type {
    fn default() -> Self {
        Cell2Type::Fluid
    }
}

// Cellule modifiée pour une grille décalée
#[derive(Clone, Debug, Default)]
pub struct Cell2 {
    pub velocity: DoubleBuffer<Vector22>,  // Vélocités sur les arêtes
    pub pressure: f32,                    // Pression au centre
    pub density: DoubleBuffer<f32>,       // Densité au centre (fumée/colorant)
    pub div: f32,                         // Divergence calculée
    pub cell_type: Cell2Type,             // Type de cellule (fluide ou solide)

    // Facteurs de relaxation pour la résolution de l'incompressibilité
    pub s_nbs: [Vector22; 2],              // Facteurs s pour les voisins [négatifs, positifs]
    pub s_tot_inv: f32,                   // Inverse de la somme des facteurs s
}

impl Cell2 {
    pub fn new() -> Self {
        Self {
            velocity: DoubleBuffer::default(),
            pressure: 0.0,
            density: DoubleBuffer::default(),
            div: 0.0,
            cell_type: Cell2Type::Fluid,
            s_nbs: [Vector22::zeros(), Vector22::zeros()],
            s_tot_inv: 0.0,
        }
    }

    // Pour convertir l'ancien "wall" en "cell_type"
    pub fn set_solid(&mut self, is_solid: bool) {
        self.cell_type = if is_solid { Cell2Type::Solid } else { Cell2Type::Fluid };
    }

    pub fn is_solid(&self) -> bool {
        self.cell_type == Cell2Type::Solid
    }
}

// Fonctions pour le codage/décodage Morton
// Ces tables de lookup peuvent être pré-calculées pour de meilleures performances
fn morton_encode(x: usize, y: usize) -> usize {
    part1by1(x) | (part1by1(y) << 1)
}

fn morton_decode(z: usize) -> (usize, usize) {
    (compact1by1(z), compact1by1(z >> 1))
}

fn part1by1(mut n: usize) -> usize {
    n &= 0x0000ffff;
    n = (n | (n << 8)) & 0x00ff00ff;
    n = (n | (n << 4)) & 0x0f0f0f0f;
    n = (n | (n << 2)) & 0x33333333;
    n = (n | (n << 1)) & 0x55555555;
    n
}

fn compact1by1(mut n: usize) -> usize {
    n &= 0x55555555;
    n = (n | (n >> 1)) & 0x33333333;
    n = (n | (n >> 2)) & 0x0f0f0f0f;
    n = (n | (n >> 4)) & 0x00ff00ff;
    n = (n | (n >> 8)) & 0x0000ffff;
    n
}

// Grille avec gestion du codage Morton
#[derive(Clone, Debug)]
pub struct Grid2 {
    pub cells: Vec<Cell2>,
    pub width: usize,          // Nombre de cellules en largeur
    pub height: usize,         // Nombre de cellules en hauteur
    pub cell_size: f32,        // Taille d'une cellule

    // Offsets pour les vélocités sur la grille décalée
    pub offsets: [Vector22; 2],

    // Cache pour accélérer les opérations
    idx_cache: Vec<(usize, usize)>,  // Cache pour les décodages Morton
}

impl Grid2 {
    pub fn new(width: usize, height: usize, cell_size: f32) -> Self {
        // Assurer que width et height sont des puissances de 2 pour tolérer le codage Morton
        let next_pow2_width = width.next_power_of_two();
        let next_pow2_height = height.next_power_of_two();

        // Pré-allouer le tableau avec la taille nécessaire pour la courbe Morton
        let total_cells = next_pow2_width * next_pow2_height;

        let half_cell = cell_size * 0.5;

        // Initialiser les cellules
        let cells = vec![Cell2::new(); total_cells];

        // Pré-calculer le cache des décodages Morton
        let mut idx_cache = vec![(0, 0); total_cells];
        for idx in 0..total_cells {
            idx_cache[idx] = morton_decode(idx);
        }

        let mut grid = Self {
            cells,
            width: next_pow2_width,
            height: next_pow2_height,
            cell_size,
            // Offset pour u (x) et v (y) sur la grille décalée
            offsets: [
                Vector22::new(0.0, half_cell),  // Offset u (x). Arêtes verticales
                Vector22::new(half_cell, 0.0),  // Offset v (y). Arêtes horizontales
            ],
            idx_cache,
        };

        // Configuration des murs aux frontières
        grid.set_boundary_walls();

        grid
    }

    // Configuration des murs aux frontières
    fn set_boundary_walls(&mut self) {
        // Murs en haut et en bas
        for i in 0..self.width {
            self.cell_mut_idx(i, 0).set_solid(true);
            self.cell_mut_idx(i, self.height - 1).set_solid(true);
        }

        // Murs à gauche et à droite
        for j in 0..self.height {
            self.cell_mut_idx(0, j).set_solid(true);
            self.cell_mut_idx(self.width - 1, j).set_solid(true);
        }
    }

    // Conversion index (i,j) vers index Morton
    pub fn idx(&self, i: usize, j: usize) -> usize {
        morton_encode(i, j)
    }

    // Décodage d'un index Morton vers coordonnées (i, j) utilise le cache
    pub fn decode_idx(&self, idx: usize) -> (usize, usize) {
        if idx < self.idx_cache.len() {
            self.idx_cache[idx]
        } else {
            morton_decode(idx)  // Fallback si hors du cache
        }
    }

    // Vérifier si l'index est dans les limites de la grille
    pub fn in_bounds(&self, i: usize, j: usize) -> bool {
        i < self.width && j < self.height
    }

    // Vérifier si l'index est dans les limites de la grille de simulation (sans les cellules fantômes)
    pub fn in_simulation_bounds(&self, i: usize, j: usize) -> bool {
        i > 0 && i < self.width - 1 && j > 0 && j < self.height - 1
    }

    // Accès à une cellule par index (i, j) avec vérification
    pub fn cell_idx(&self, i: usize, j: usize) -> &Cell2 {
        debug_assert!(self.in_bounds(i, j), "Cell index out of bounds: ({}, {})", i, j);
        let idx = self.idx(i, j);
        &self.cells[idx]
    }

    pub fn cell_mut_idx(&mut self, i: usize, j: usize) -> &mut Cell2 {
        debug_assert!(self.in_bounds(i, j), "Cell index out of bounds: ({}, {})", i, j);
        let idx = self.idx(i, j);
        &mut self.cells[idx]
    }

    // Retourne Some(idx) si l'index est valide, sinon None - plus sûr
    pub fn try_idx(&self, i: usize, j: usize) -> Option<usize> {
        if self.in_bounds(i, j) {
            Some(self.idx(i, j))
        } else {
            None
        }
    }

    // Version sécurisée pour obtenir une référence à une cellule
    pub fn try_cell_idx(&self, i: usize, j: usize) -> Option<&Cell2> {
        self.try_idx(i, j).map(|idx| &self.cells[idx])
    }

    pub fn try_cell_mut_idx(&mut self, i: usize, j: usize) -> Option<&mut Cell2> {
        if let Some(idx) = self.try_idx(i, j) {
            Some(&mut self.cells[idx])
        } else {
            None
        }
    }

    // Itérer sur toutes les cellules en ordre Morton — plus efficace avec le cache
    pub fn iter_morton(&self) -> impl Iterator<Item = (usize, usize, usize)> + '_ {
        self.cells
            .iter()
            .enumerate()
            .map(move |(idx, _)| {
                let (i, j) = self.decode_idx(idx);
                (i, j, idx)
            })
    }

    // Itérer uniquement sur les cellules de simulation (évite les bords)
    pub fn iter_simulation_cells(&self) -> impl Iterator<Item = (usize, usize, usize)> + '_ {
        self.iter_morton()
            .filter(move |&(i, j, _)| self.in_simulation_bounds(i, j))
    }

    // Ajouter un obstacle circulaire
    pub fn add_circle_obstacle(&mut self, center_x: f32, center_y: f32, radius: f32, vel: Option<Vector22>) {
        let velocity = vel.unwrap_or(Vector22::zeros());

        // Parcourir toutes les cellules de la grille de simulation
        for j in 1..self.height-1 {
            for i in 1..self.width-1 {
                let cell_center_x = (i as f32 + 0.5) * self.cell_size;
                let cell_center_y = (j as f32 + 0.5) * self.cell_size;

                let dx = cell_center_x - center_x;
                let dy = cell_center_y - center_y;
                let dist_squared = dx*dx + dy*dy;

                if dist_squared <= radius*radius {
                    let cell = self.cell_mut_idx(i, j);
                    cell.set_solid(true);
                    cell.velocity.back = velocity;
                }
            }
        }
    }

    // Obtenir les indices des voisins d'une cellule [négatifs, positifs] pour chaque direction
    // avec vérification des limites
    fn get_neighbor_indices(&self, i: usize, j: usize) -> [[Option<usize>; 2]; 2] {
        [
            [
                if i > 0 { Some(self.idx(i-1, j)) } else { None },
                if j > 0 { Some(self.idx(i, j-1)) } else { None }
            ],  // Voisins négatifs (gauche, haut)
            [
                if i + 1 < self.width { Some(self.idx(i+1, j)) } else { None },
                if j + 1 < self.height { Some(self.idx(i, j+1)) } else { None }
            ],  // Voisins positifs (droite, bas)
        ]
    }

    // === ÉTAPE PRINCIPALE DE SIMULATION ===
    // Amélioration de la fonction vel_step avec plus de diagnostics
    pub fn vel_step(&mut self, dt: f32, gravity: Vector22) {
        if LOG == true {
            // Vérifier l'état initial
            println!("----- ÉTAT INITIAL -----");
            self.check_velocity_noise();
            println!("Total densité initial: {:.2}", self.total_density());

            // 1. Appliquer les forces externes (comme la gravité)
            self.apply_forces(dt, gravity);
            println!("----- APRÈS FORCES -----");
            self.check_velocity_noise();
            println!("Total densité après forces externes: {:.2}", self.total_density());

            // 2. Résoudre l'incompressibilité
            self.solve_incompressibility(dt, 20, 1000.0, 1.9);
            println!("----- APRÈS INCOMPRESSIBILITÉ -----");
            self.check_velocity_noise();
            println!("Total densité après solve incompressibility: {:.2}", self.total_density());

            // 3. Advection
            self.advect_velocity(dt);
            println!("----- APRÈS ADVECTION VITESSE -----");
            self.check_velocity_noise();
            println!("Total densité après advect velocity: {:.2}", self.total_density());

            self.advect_density(dt);
            println!("----- APRÈS ADVECTION DENSITÉ -----");
            self.check_velocity_noise();
            println!("Total densité après advect density: {:.2}", self.total_density());

            // 4. Appliquer les conditions aux limites
            self.apply_boundary_conditions();
            println!("----- APRÈS CONDITIONS AUX LIMITES -----");
            self.check_velocity_noise();
            println!("Total densité après apply boundary conditions: {:.2}", self.total_density());
        } else if LOG == false {
            // Version sans logs

            // 1. Appliquer les forces externes (comme la gravité)
            self.apply_forces(dt, gravity);

            // 2. Résoudre l'incompressibilité
            self.solve_incompressibility(dt, 20, 1000.0, 1.9);

            // 3. Advection
            self.advect_velocity(dt);
            self.advect_density(dt);

            // 4. Appliquer les conditions aux limites
            self.apply_boundary_conditions();
        }
    }

    // === FORCES EXTERNES ===
    pub fn apply_forces(&mut self, dt: f32, gravity: Vector22) {
        // Version Morton — itérer sur les cellules en ordre Morton
        // Collecter d'abord les indices des cellules fluides
        let fluid_indices: Vec<_> = self.iter_simulation_cells()
            .filter(|&(_, _, idx)| self.cells[idx].cell_type == Cell2Type::Fluid)
            .map(|(_, _, idx)| idx)
            .collect();

        // Appliquer ensuite les modifications
        for &idx in &fluid_indices {
            self.cells[idx].velocity.back.y += gravity.y * dt;
            self.cells[idx].velocity.back.x += gravity.x * dt;
        }
    }

    // === RÉSOLUTION DE L'INCOMPRESSIBILITÉ ===
    // Améliorons aussi le solveur d'incompressibilité pour garantir une meilleure convergence
    pub fn solve_incompressibility(&mut self, dt: f32, iterations: usize, density: f32, relaxation: f32) {
        let cp = density * self.cell_size / dt;

        // Fonction pour déterminer si une cellule contribue à la pression (1.0 pour fluide, 0.0 pour solide)
        let s_factor = |cell: &Cell2| -> f32 {
            if cell.cell_type == Cell2Type::Solid { 0.0 } else { 1.0 }
        };

        // Calculer les facteurs s pour tous les voisins — version sécurisée
        for j in 1..self.height-1 {
            for i in 1..self.width-1 {
                let idx = self.idx(i, j);

                // Vérifier si la cellule actuelle est solide
                if self.cells[idx].cell_type == Cell2Type::Solid {
                    continue;
                }

                // Vérifier et calculer les indices des voisins de manière sûre
                let has_neg_i = i > 0;
                let has_neg_j = j > 0;
                let has_pos_i = i + 1 < self.width;
                let has_pos_j = j + 1 < self.height;

                let neg_i = if has_neg_i { self.idx(i-1, j) } else { idx };
                let neg_j = if has_neg_j { self.idx(i, j-1) } else { idx };
                let pos_i = if has_pos_i { self.idx(i+1, j) } else { idx };
                let pos_j = if has_pos_j { self.idx(i, j+1) } else { idx };

                // Calculer tous les facteurs s d'abord
                let cell_s = s_factor(&self.cells[idx]);
                let s_neg_i = if has_neg_i { s_factor(&self.cells[neg_i]) } else { 0.0 };
                let s_neg_j = if has_neg_j { s_factor(&self.cells[neg_j]) } else { 0.0 };
                let s_pos_i = if has_pos_i { s_factor(&self.cells[pos_i]) } else { 0.0 };
                let s_pos_j = if has_pos_j { s_factor(&self.cells[pos_j]) } else { 0.0 };

                // Mettre à jour les facteurs s pour tous les voisins de manière sécurisée
                if has_neg_i {
                    self.cells[neg_i].s_nbs[1].x = cell_s;
                }
                if has_neg_j {
                    self.cells[neg_j].s_nbs[1].y = cell_s;
                }

                self.cells[idx].s_nbs[1].x = s_pos_i;
                self.cells[idx].s_nbs[1].y = s_pos_j;
                self.cells[idx].s_nbs[0].x = s_neg_i;
                self.cells[idx].s_nbs[0].y = s_neg_j;
            }
        }

        // Calculer la somme inverse des facteurs s pour chaque cellule
        let cell_updates: Vec<_> = self.iter_simulation_cells()
            .filter(|&(_, _, idx)| self.cells[idx].cell_type != Cell2Type::Solid)
            .map(|(_, _, idx)| {
                let sum = self.cells[idx].s_nbs[0].x + self.cells[idx].s_nbs[0].y +
                    self.cells[idx].s_nbs[1].x + self.cells[idx].s_nbs[1].y;
                (idx, sum)
            })
            .collect();

        for (idx, sum) in cell_updates {
            self.cells[idx].pressure = 0.0;
            self.cells[idx].s_tot_inv = if sum > 0.0 { 1.0 / sum } else { 0.0 };
        }

        // CORRECTION: Initialisation de max_div avec 0.0 au lieu de std::f32::MAX
        let mut max_div = 0.0;

        // Vérification initiale de la divergence
        let mut initial_div = 0.0;
        for j in 1..self.height-1 {
            for i in 1..self.width-1 {
                let idx = self.idx(i, j);

                // Ignorer les cellules solides
                if self.cells[idx].cell_type == Cell2Type::Solid {
                    continue;
                }

                // Vérifier les limites
                if i + 1 >= self.width || j + 1 >= self.height {
                    continue;
                }

                let pos_i = self.idx(i+1, j);
                let pos_j = self.idx(i, j+1);

                // Calculer la divergence
                let vel_x = self.cells[idx].velocity.back.x;
                let vel_y = self.cells[idx].velocity.back.y;
                let vel_pos_i = self.cells[pos_i].velocity.back.x;
                let vel_pos_j = self.cells[pos_j].velocity.back.y;

                let div_x = vel_pos_i - vel_x;
                let div_y = vel_pos_j - vel_y;
                let div = div_x + div_y;

                // Si la divergence initiale est non nulle, nous devrions itérer
                if div.abs() > initial_div {
                    initial_div = div.abs();
                }
            }
        }

        // Si la divergence initiale est déjà très petite, pas besoin d'itérer
        let skip_iterations = initial_div < 1e-6;
        if skip_iterations {
            if LOG == true {
                println!("Divergence initiale très faible ({:.6}), itérations ignorées", initial_div);
            }
            max_div = initial_div;
        } else {
            // Effectuer les itérations de résolution de pression
            for iter in 0..iterations {
                let mut current_max_div : f32 = 0.0;

                for j in 1..self.height-1 {
                    for i in 1..self.width-1 {
                        let idx = self.idx(i, j);

                        // Vérifions d'abord le type de cellule
                        if self.cells[idx].cell_type == Cell2Type::Solid {
                            continue;
                        }

                        // Vérifier et calculer les indices des voisins de manière sûre
                        let has_pos_i = i + 1 < self.width;
                        let has_pos_j = j + 1 < self.height;

                        if !has_pos_i || !has_pos_j {
                            continue;  // Ignorer les cellules au bord
                        }

                        // Voisins positifs (droite, bas)
                        let pos_i = self.idx(i+1, j);
                        let pos_j = self.idx(i, j+1);

                        // Calculer la divergence
                        let vel_x = self.cells[idx].velocity.back.x;
                        let vel_y = self.cells[idx].velocity.back.y;
                        let vel_pos_i = self.cells[pos_i].velocity.back.x;
                        let vel_pos_j = self.cells[pos_j].velocity.back.y;
                        let s_nbs_0 = self.cells[idx].s_nbs[0];
                        let s_nbs_1 = self.cells[idx].s_nbs[1];
                        let s_tot_inv = self.cells[idx].s_tot_inv;

                        let div_x = vel_pos_i - vel_x;
                        let div_y = vel_pos_j - vel_y;
                        let div = div_x + div_y;
                        let div_normed = div * s_tot_inv;

                        // Mise à jour de la pression
                        self.cells[idx].pressure -= cp * div_normed;

                        // AMÉLIORATION: Suivi de la divergence max
                        current_max_div = current_max_div.max(div_normed.abs());

                        // Mise à jour des vitesses — version sécurisée
                        self.cells[idx].velocity.back.x += relaxation * s_nbs_0.x * div_normed;
                        self.cells[idx].velocity.back.y += relaxation * s_nbs_0.y * div_normed;

                        // Vérifier à nouveau avant d'accéder aux cellules voisines
                        if self.cells[pos_i].cell_type != Cell2Type::Solid {
                            self.cells[pos_i].velocity.back.x -= relaxation * s_nbs_1.x * div_normed;
                        }

                        if self.cells[pos_j].cell_type != Cell2Type::Solid {
                            self.cells[pos_j].velocity.back.y -= relaxation * s_nbs_1.y * div_normed;
                        }
                    }
                }

                // AMÉLIORATION: Vérification de convergence
                if LOG == true {
                    if iter % 5 == 0 {  // Affichage tous les 5 itérations pour éviter de spammer
                        println!("Itération {}: divergence max = {:.6}", iter, current_max_div);
                    }
                }

                // CORRECTION: Mise à jour de max_div
                max_div = current_max_div;

                // Sortir si la convergence est suffisante
                
                if current_max_div < 1e-5 {
                    if LOG == true {
                        println!("Convergence atteinte après {} itérations (divergence = {:.6})", iter + 1, current_max_div);
                    }
                    break;
                }
            }
        }
        if LOG == true {
            println!("Divergence finale après {} itérations: {:.6}", iterations, max_div);
        }
    }

    // === ADVECTION DE LA VITESSE ===
    pub fn advect_velocity(&mut self, dt: f32) {
        // Copier les vitesses actuelles dans le front buffer
        for cell in self.cells.iter_mut() {
            cell.velocity.front = cell.velocity.back;
        }

        // Advection pour chaque cellule de la grille de simulation
        for j in 1..self.height-1 {
            for i in 1..self.width-1 {
                let idx = self.idx(i, j);

                if self.cells[idx].cell_type == Cell2Type::Solid {
                    continue;
                }

                // Advection de la vitesse en x (u) — avec vérification
                if i > 0 && self.cells[self.idx(i-1, j)].cell_type != Cell2Type::Solid {
                    let pos_x = (i as f32) * self.cell_size;
                    let pos_y = (j as f32 + 0.5) * self.cell_size;

                    // Échantillonnage bilinéaire de la vitesse à cette position
                    let vel = self.sample_velocity_at(pos_x, pos_y);

                    // Position de départ de la particule
                    let start_x = pos_x - dt * vel.x;
                    let start_y = pos_y - dt * vel.y;

                    // Vitesse à la position de départ
                    self.cells[idx].velocity.front.x = self.sample_velocity_component_at(start_x, start_y, 0);
                }

                // Advection de la vitesse en y (v) — avec vérification
                if j > 0 && self.cells[self.idx(i, j-1)].cell_type != Cell2Type::Solid {
                    let pos_x = (i as f32 + 0.5) * self.cell_size;
                    let pos_y = (j as f32) * self.cell_size;

                    // Échantillonnage bilinéaire de la vitesse à cette position
                    let vel = self.sample_velocity_at(pos_x, pos_y);

                    // Position de départ de la particule
                    let start_x = pos_x - dt * vel.x;
                    let start_y = pos_y - dt * vel.y;

                    // Vitesse à la position de départ
                    self.cells[idx].velocity.front.y = self.sample_velocity_component_at(start_x, start_y, 1);
                }
            }
        }

        // Échanger les buffers de vitesse
        for cell in self.cells.iter_mut() {
            cell.velocity.swap();
        }
    }

    // === ADVECTION DE LA DENSITÉ ===
    pub fn advect_density(&mut self, dt: f32) {
        // Copier les densités actuelles dans le front buffer
        for cell in self.cells.iter_mut() {
            cell.density.front = cell.density.back;
        }

        // Advection pour chaque cellule de la grille de simulation
        for j in 1..self.height-1 {
            for i in 1..self.width-1 {
                let idx = self.idx(i, j);

                if self.cells[idx].cell_type == Cell2Type::Solid {
                    continue;
                }

                // Position au centre de la cellule
                let pos_x = (i as f32 + 0.5) * self.cell_size;
                let pos_y = (j as f32 + 0.5) * self.cell_size;

                // Échantillonnage bilinéaire de la vitesse à cette position
                let vel = self.sample_velocity_at(pos_x, pos_y);

                // Position de départ de la particule (méthode de "backtracing")
                let start_x = pos_x - dt * vel.x;
                let start_y = pos_y - dt * vel.y;

                // Ajouter des logs de diagnostic pour les cellules contenant de la densité
                if LOG == true {
                    if self.cells[idx].density.back > 0.1 {
                        println!(
                            "Advection densité Cellule ({}, {}): Densité: {:.3}, Vitesse: ({:.3}, {:.3}), Trace vers: ({:.3}, {:.3})",
                            i, j,
                            self.cells[idx].density.back,
                            vel.x, vel.y,
                            start_x / self.cell_size, start_y / self.cell_size
                        );
                    }
                }

                // Densité à la position de départ
                self.cells[idx].density.front = self.sample_density_at(start_x, start_y);
            }
        }

        // Échanger les buffers de densité
        for cell in self.cells.iter_mut() {
            cell.density.swap();
        }
    }

    pub fn advect_density2(&mut self, dt: f32) {
        let dt0 = dt * self.width as f32; // Utiliser width comme facteur de normalisation (équivalent à N)

        // Collecter toutes les nouvelles densités en parallèle
        let new_densities: Vec<f32> = self.cells.par_iter().enumerate().map(|(idx, cell)| {
            // Décodage de l'index
            let i = idx % self.width;
            let j = idx / self.width;

            // Préserver la densité pour les cellules solides et les frontières
            if cell.cell_type == Cell2Type::Solid ||
                i == 0 || i >= self.width - 1 ||
                j == 0 || j >= self.height - 1 {
                return cell.density.back;
            }

            // Position au centre de la cellule
            let pos_x = (i as f32 + 0.5) * self.cell_size;
            let pos_y = (j as f32 + 0.5) * self.cell_size;

            // Échantillonnage du vecteur vitesse à cette position
            let vel = self.sample_velocity_at(pos_x, pos_y);

            // Backtracking: d'où venait la particule de fluide
            let start_x = pos_x - dt0 * vel.x;
            let start_y = pos_y - dt0 * vel.y;

            // Limiter avec réflexion aux frontières
            let start_x = start_x.max(0.5 * self.cell_size).min((self.width as f32 - 0.5) * self.cell_size);
            let start_y = start_y.max(0.5 * self.cell_size).min((self.height as f32 - 0.5) * self.cell_size);

            // Ajout de logs de diagnostic pour les cellules contenant de la densité (si activé)
            if LOG == true && cell.density.back > 0.1 {
                println!(
                    "Advection densité Cellule ({}, {}): Densité: {:.3}, Vitesse: ({:.3}, {:.3}), Trace vers: ({:.3}, {:.3})",
                    i, j,
                    cell.density.back,
                    vel.x, vel.y,
                    start_x / self.cell_size, start_y / self.cell_size
                );
            }

            // Fonction modifiée pour échantillonner la densité avec réflexion aux murs
            let sample_density_with_reflection = |x: f32, y: f32| -> f32 {
                // Convertir en indices de grille
                let grid_x = x / self.cell_size;
                let grid_y = y / self.cell_size;

                // Calculer les indices et poids pour l'interpolation
                let i0 = grid_x.floor() as usize;
                let j0 = grid_y.floor() as usize;

                // Limiter les indices aux bornes du domaine
                let i0 = i0.clamp(0, self.width - 2);
                let j0 = j0.clamp(0, self.height - 2);
                let i1 = (i0 + 1).min(self.width - 1);
                let j1 = (j0 + 1).min(self.height - 1);

                // Calculer les facteurs d'interpolation
                let s1 = (grid_x - i0 as f32).clamp(0.0, 1.0);
                let s0 = 1.0 - s1;
                let t1 = (grid_y - j0 as f32).clamp(0.0, 1.0);
                let t0 = 1.0 - t1;

                // Fonction pour obtenir la densité avec réflexion aux murs
                let get_density = |i: usize, j: usize| -> f32 {
                    let idx = self.idx(i, j);
                    if self.cells[idx].cell_type == Cell2Type::Solid {
                        // Réflexion: utiliser la densité de la cellule actuelle
                        cell.density.back
                    } else {
                        self.cells[idx].density.back
                    }
                };

                // Interpolation bilinéaire
                s0 * (t0 * get_density(i0, j0) + t1 * get_density(i0, j1)) +
                    s1 * (t0 * get_density(i1, j0) + t1 * get_density(i1, j1))
            };

            // Échantillonner la densité à la position de départ avec réflexion
            let new_density = sample_density_with_reflection(start_x, start_y);

            // Éviter les valeurs négatives
            new_density.max(0.0)
        }).collect();

        // Mettre à jour les densités en une seule passe
        for (idx, &new_density) in new_densities.iter().enumerate() {
            self.cells[idx].density.front = new_density;
        }

        // Échanger les buffers de densité
        for cell in self.cells.iter_mut() {
            cell.density.swap();
        }
    }



    // === CONDITIONS AUX LIMITES ===
    pub fn apply_boundary_conditions(&mut self) {
        // Pour chaque cellule sur les bords (sans les coins)

        // Bord gauche (entrée)
        for j in 1..self.height-1 {
            let idx = self.idx(1, j);
            let cell = &mut self.cells[idx];
            cell.velocity.back.x = 0.0;
            cell.velocity.back.y = 0.0;
        }

        // Bord droit (sortie libre)
        for j in 1..self.height-1 {
            let i = self.width - 2;
            let idx = self.idx(i, j);

            // Vérifier si la cellule à gauche est disponible
            if i > 0 && self.in_bounds(i-1, j) {
                let idx_left = self.idx(i-1, j);
                self.cells[idx].velocity.back.x = self.cells[idx_left].velocity.back.x;
            }
        }

        // Bords haut et bas (pas de glissement)
        for i in 1..self.width-1 {
            // Bord haut
            if self.in_bounds(i, 1) {
                let idx_top = self.idx(i, 1);
                self.cells[idx_top].velocity.back.y = 0.0;
            }

            // Bord bas
            if self.in_bounds(i, self.height-2) {
                let idx_bottom = self.idx(i, self.height-2);
                self.cells[idx_bottom].velocity.back.y = 0.0;
            }
        }

        // S'assurer que les cellules solides ont une vitesse nulle
        let solid_indices: Vec<_> = self.iter_morton()
            .filter_map(|(_i, _j, idx)| {
                if self.cells[idx].cell_type == Cell2Type::Solid {
                    Some(idx)
                } else {
                    None
                }
            })
            .collect();

        for &idx in solid_indices.iter() {
            self.cells[idx].velocity.back = Vector22::zeros();
        }

        // Extrapolation des vitesses aux bords
        self.extrapolate_velocity();
    }


    // Extrapolation des vitesses aux bords
    // Amélioration de l'extrapolation des vitesses aux bords pour éviter la propagation d'erreurs
    fn extrapolate_velocity(&mut self) {
        // MODIFICATION: Vérification des valeurs avant extrapolation
        let mut has_non_zero = false;

        // Extrapolation horizontale (pour u)
        for i in 1..self.width-1 {
            // Copier la vitesse u du rang 1 vers le rang 0
            let vel_x = self.cells[self.idx(i, 1)].velocity.back.x;
            if vel_x.abs() > 1e-6 {
                has_non_zero = true;
            }
            let idx = self.idx(i, 0);
            self.cells[idx].velocity.back.x = vel_x;

            // Copier la vitesse u du rang height-2 vers le rang height-1
            let j_last = self.height - 1;
            let j_before_last = self.height - 2;
            let vel_x_last = self.cells[self.idx(i, j_before_last)].velocity.back.x;
            if vel_x_last.abs() > 1e-6 {
                has_non_zero = true;
            }
            let idx = self.idx(i, j_last);
            self.cells[idx].velocity.back.x = vel_x_last;
        }

        // Extrapolation verticale (pour v)
        for j in 1..self.height-1 {
            // Copier la vitesse v du rang 1 vers le rang 0
            let vel_y = self.cells[self.idx(1, j)].velocity.back.y;
            if vel_y.abs() > 1e-6 {
                has_non_zero = true;
            }
            let idx = self.idx(0, j);
            self.cells[idx].velocity.back.y = vel_y;

            // Copier la vitesse v du rang width-2 vers le rang width-1
            let i_last = self.width - 1;
            let i_before_last = self.width - 2;
            let vel_y_last = self.cells[self.idx(i_before_last, j)].velocity.back.y;
            if vel_y_last.abs() > 1e-6 {
                has_non_zero = true;
            }
            let idx = self.idx(i_last, j);
            self.cells[idx].velocity.back.y = vel_y_last;
        }

        if has_non_zero && LOG == true {
            println!("Attention: Extrapolation de vitesses non nulles aux bords");
        }
    }

    // === FONCTIONS D'ÉCHANTILLONNAGE ===

    // Échantillonnage bilinéaire d'une composante de vitesse à une position
    // Échantillonnage bilinéaire d'une composante de vitesse à une position
    fn sample_velocity_component_at(&self, x: f32, y: f32, component: usize) -> f32 {
        // Limiter la position à l'intérieur de la grille avec une marge
        let x = x.max(0.0).min((self.width as f32) * self.cell_size);
        let y = y.max(0.0).min((self.height as f32) * self.cell_size);

        // Pour tenir compte de la grille décalée, ajuster les coordonnées en fonction de la composante
        let (px, py) = match component {
            0 => (x, y - self.offsets[0].y), // Pour u (composante x)
            1 => (x - self.offsets[1].x, y), // Pour v (composante y)
            _ => panic!("Composant invalide")
        };

        // Position dans la grille (en unités de cellule)
        let grid_x = px / self.cell_size;
        let grid_y = py / self.cell_size;

        // Cellule contenant la position avec limites sécurisées
        let i0 = grid_x.floor() as usize;
        let j0 = grid_y.floor() as usize;

        // Assurer que les indices sont dans les limites
        let i0 = i0.clamp(0, self.width - 2);
        let j0 = j0.clamp(0, self.height - 2);
        let i1 = (i0 + 1).min(self.width - 1);
        let j1 = (j0 + 1).min(self.height - 1);

        // Recalculer les facteurs d'interpolation basés sur les indices ajustés
        let s1 = (grid_x - i0 as f32).clamp(0.0, 1.0);
        let s0 = 1.0 - s1;
        let t1 = (grid_y - j0 as f32).clamp(0.0, 1.0);
        let t0 = 1.0 - t1;

        // Échantillonner la composante de vitesse de manière sécurisée
        let v00 = self.cells[self.idx(i0, j0)].velocity.back[component];
        let v01 = self.cells[self.idx(i0, j1)].velocity.back[component];
        let v10 = self.cells[self.idx(i1, j0)].velocity.back[component];
        let v11 = self.cells[self.idx(i1, j1)].velocity.back[component];

        // Interpolation bilinéaire
        s0 * (t0 * v00 + t1 * v01) + s1 * (t0 * v10 + t1 * v11)
    }

    // Échantillonnage bilinéaire du vecteur vitesse complet à une position
    fn sample_velocity_at(&self, x: f32, y: f32) -> Vector22 {
        Vector22 {
            x: self.sample_velocity_component_at(x, y, 0),
            y: self.sample_velocity_component_at(x, y, 1),
        }
    }

    // Échantillonnage bilinéaire de la densité à une position
    fn sample_density_at(&self, x: f32, y: f32) -> f32 {
        // Limiter la position à l'intérieur de la grille
        let x = x.max(0.0).min((self.width as f32) * self.cell_size);
        let y = y.max(0.0).min((self.height as f32) * self.cell_size);

        // Position dans la grille (en unités de cellule)
        let grid_x = x / self.cell_size;
        let grid_y = y / self.cell_size;

        // Cellule contenant la position avec limites sécurisées
        let i0 = grid_x.floor() as usize;
        let j0 = grid_y.floor() as usize;

        // Assurer que les indices sont dans les limites
        let i0 = i0.clamp(0, self.width - 2);
        let j0 = j0.clamp(0, self.height - 2);
        let i1 = (i0 + 1).min(self.width - 1);
        let j1 = (j0 + 1).min(self.height - 1);

        // Recalculer les facteurs d'interpolation
        let s1 = (grid_x - i0 as f32).clamp(0.0, 1.0);
        let s0 = 1.0 - s1;
        let t1 = (grid_y - j0 as f32).clamp(0.0, 1.0);
        let t0 = 1.0 - t1;

        // Échantillonner la densité de manière sécurisée
        let d00 = self.cells[self.idx(i0, j0)].density.back;
        let d01 = self.cells[self.idx(i0, j1)].density.back;
        let d10 = self.cells[self.idx(i1, j0)].density.back;
        let d11 = self.cells[self.idx(i1, j1)].density.back;

        // Interpolation bilinéaire
        s0 * (t0 * d00 + t1 * d01) + s1 * (t0 * d10 + t1 * d11)
    }

    // === FONCTIONS UTILITAIRES ===

    pub fn test_density_advection(&mut self, dt: f32, steps: usize) {
        // Ajouter une source de densité au centre
        let center_i = self.width / 2;
        let center_j = self.height / 2;

        // Ajouter de la densité dans une zone 3x3
        for di in -1..=1 {
            for dj in -1..=1 {
                let i = (center_i as isize + di) as usize;
                let j = (center_j as isize + dj) as usize;
                self.add_density_source(i, j, 1.0);
            }
        }

        // Ajouter une vitesse constante (par exemple, nord-est)
        for i in 1..self.width-1 {
            for j in 1..self.height-1 {
                let idx = self.idx(i, j);
                if self.cells[idx].cell_type == Cell2Type::Fluid {
                    // Vitesse constante vers la droite et le bas
                    self.cells[idx].velocity.back.x = 0.5;  // Vitesse vers la droite
                    self.cells[idx].velocity.back.y = 0.5;  // Vitesse vers le bas
                }
            }
        }

        println!("=== TEST ADVECTION DE DENSITÉ ===");
        println!("Configuration initiale:");
        self.print_grid_density();

        for step in 1..=steps {
            // Advection de la densité seulement
            self.advect_density(dt);

            // Afficher l'état
            println!("\nApres {} pas de temps:", step);
            self.print_grid_density();

            // Calculer le centre de masse de la densité
            let mut total_density = 0.0;
            let mut weighted_x = 0.0;
            let mut weighted_y = 0.0;

            for j in 1..self.height-1 {
                for i in 1..self.width-1 {
                    let idx = self.idx(i, j);
                    let d = self.cells[idx].density.back;
                    total_density += d;
                    weighted_x += i as f32 * d;
                    weighted_y += j as f32 * d;
                }
            }

            if total_density > 0.0 {
                let center_of_mass_x = weighted_x / total_density;
                let center_of_mass_y = weighted_y / total_density;
                println!("Centre de masse: ({:.2}, {:.2})", center_of_mass_x, center_of_mass_y);
            }
        }
    }

    pub fn check_velocity_noise(&self) -> (f32, f32, f32) {
        let mut sum_vel_x = 0.0;
        let mut sum_vel_y = 0.0;
        let mut count = 0;

        // Calculer la vitesse moyenne dans toute la grille
        for (_, _, idx) in self.iter_simulation_cells() {
            if self.cells[idx].cell_type == Cell2Type::Fluid {
                sum_vel_x += self.cells[idx].velocity.back.x;
                sum_vel_y += self.cells[idx].velocity.back.y;
                count += 1;
            }
        }

        let avg_vel_x = if count > 0 { sum_vel_x / count as f32 } else { 0.0 };
        let avg_vel_y = if count > 0 { sum_vel_y / count as f32 } else { 0.0 };
        let avg_magnitude = (avg_vel_x * avg_vel_x + avg_vel_y * avg_vel_y).sqrt();

        println!("Vitesse résiduelle moyenne: ({:.6}, {:.6}), magnitude: {:.6}",
                 avg_vel_x, avg_vel_y, avg_magnitude);

        (avg_vel_x, avg_vel_y, avg_magnitude)
    }



    // Ajouter une source de densité
    pub fn add_density_source(&mut self, i: usize, j: usize, amount: f32) {
        if self.in_simulation_bounds(i, j) {
            let idx = self.idx(i, j);
            if self.cells[idx].cell_type == Cell2Type::Fluid {
                self.cells[idx].density.back += amount;
            }
        }
    }

    // Calculer la densité totale dans la grille (pour vérification de conservation)
    pub fn total_density(&self) -> f32 {
        self.cells.iter()
            .filter(|cell| cell.cell_type == Cell2Type::Fluid)
            .map(|cell| cell.density.back)
            .sum()
    }
   

    /// Initialiser une cellule avec une vitesse et densité données
    pub fn cell_init(&mut self, i: usize, j: usize, vx: f32, vy: f32, density: f32) {
        if self.in_bounds(i, j) {
            let idx = self.idx(i, j);
            if self.cells[idx].cell_type == Cell2Type::Fluid {
                self.cells[idx].velocity.back = Vector22::new(vx, vy);
                self.cells[idx].density.back = density;
            } else {
                println!("Impossible de modifier une cellule solide!");
            }
        } else {
            println!("Erreur: indices hors limites!");
        }
    }

    /// Initialiser le caractère solide d'une cellule
    pub fn wall_init(&mut self, i: usize, j: usize, is_solid: bool) {
        if self.in_bounds(i, j) {
            let idx = self.idx(i, j);
            self.cells[idx].set_solid(is_solid);
        } else {
            println!("Erreur: indices hors limites!");
        }
    }

    /// Initialiser la vitesse d'une cellule
    pub fn velocity_init(&mut self, i: usize, j: usize, vx: f32, vy: f32) {
        if self.in_bounds(i, j) {
            let idx = self.idx(i, j);
            if self.cells[idx].cell_type == Cell2Type::Fluid {
                self.cells[idx].velocity.back = Vector22::new(vx, vy);
            } else {
                println!("Impossible de modifier une cellule solide!");
            }
        } else {
            println!("Erreur: indices hors limites!");
        }
    }

    /// Ajouter une source circulaire de densité et de vitesse au centre de la grille
    pub fn center_source(&mut self, radius: f32, density_value: f32, velocity_magnitude: f32, is_radial: bool) {
        // Calculer le centre de la grille
        let center_x = (self.width as f32 / 2.0) as usize;
        let center_y = (self.height as f32 / 2.0) as usize;

        // Déterminer la zone à scanner (carré inscrit dans le cercle)
        let r_int = radius.floor() as usize;
        let start_i = if center_x > r_int { center_x - r_int } else { 1 };
        let end_i = (center_x + r_int).min(self.width - 1);
        let start_j = if center_y > r_int { center_y - r_int } else { 1 };
        let end_j = (center_y + r_int).min(self.height - 1);

        for i in start_i..=end_i {
            for j in start_j..=end_j {
                // Calculer la distance au centre
                let dx = i as f32 - center_x as f32;
                let dy = j as f32 - center_y as f32;
                let distance_squared = dx * dx + dy * dy;

                // Si la cellule courante n'est pas une cellule solide et est dans le rayon
                if distance_squared <= radius * radius {
                    let idx = self.idx(i, j);
                    if self.cells[idx].cell_type == Cell2Type::Fluid {
                        // Ajouter de la densité
                        self.cells[idx].density.back += density_value;

                        // Calculer la direction de la vitesse
                        if is_radial {
                            // Initialiser une vitesse radiale
                            if distance_squared > 0.0 {
                                let distance = distance_squared.sqrt();
                                let dir_x = dx / distance;
                                let dir_y = dy / distance;

                                // Ajouter une vitesse dans la direction du centre
                                let factor = 1.0 - (distance / radius); // Plus fort au centre
                                self.cells[idx].velocity.back.x += dir_x * velocity_magnitude * factor;
                                self.cells[idx].velocity.back.y += dir_y * velocity_magnitude * factor;
                            }
                        } else {
                            // Vitesse circulaire (tourbillon)
                            if distance_squared > 0.0 {
                                let distance = distance_squared.sqrt();
                                let dir_x = -dy / distance;
                                let dir_y = dx / distance;

                                // Vitesse qui augmente avec la distance jusqu'à un certain point
                                let factor = (distance / radius) * (1.0 - distance / radius) * 4.0;
                                self.cells[idx].velocity.back.x += dir_x * velocity_magnitude * factor;
                                self.cells[idx].velocity.back.y += dir_y * velocity_magnitude * factor;
                            }
                        }
                    }
                }
            }
        }
    }

    /// Afficher la grille avec les murs/cellules solides
    pub fn print_grid_walls(&self) {
        for j in 0..self.height {
            let mut line = String::with_capacity(self.width * 3);
            for i in 0..self.width {
                let is_solid = self.cell_idx(i, j).is_solid();
                line.push_str(if is_solid { " █ " } else { " · " });
            }
            println!("{}", line);
        }
    }

    /// Afficher la grille avec les vélocités
    pub fn print_grid_velocity(&self) {
        for j in 0..self.height {
            let mut line = String::with_capacity(self.width * 18);
            for i in 0..self.width {
                let cell = self.cell_idx(i, j);
                if cell.is_solid() {
                    line.push_str("   | █████ |   ");
                } else {
                    let v = &cell.velocity.back;
                    line.push_str(&format!("|{:.3}, {:.3}| ", v.x, v.y));
                }
            }
            println!("{}", line);
        }
    }

    /// Afficher la grille avec les densités
    pub fn print_grid_density(&self) {
        for j in 0..self.height {
            let mut line = String::with_capacity(self.width * 10);
            for i in 0..self.width {
                let cell = self.cell_idx(i, j);
                if cell.is_solid() {
                    line.push_str("   |█████|   ");
                } else {
                    line.push_str(&format!("  |{:5.2}|  ", cell.density.back));
                }
            }
            println!("{}", line);
        }
    }

    /// Créer un cercle dans la grille
    pub fn circle(&mut self, center_x: isize, center_y: isize, radius: f32) {
        let center_x = center_x as usize;
        let center_y = center_y as usize;

        for i in 1..self.width-1 {
            for j in 1..self.height-1 {
                let dx = i as isize - center_x as isize;
                let dy = j as isize - center_y as isize;
                if dx * dx + dy * dy <= (radius * radius) as isize {
                    self.wall_init(i, j, true);
                }
            }
        }
    }

    /// Créer une configuration type tunnel de vent
    pub fn initialize_wind_tunnel(&mut self, density: f32, hole_positions: &[usize], flow_velocity: f32) {
        let left_wall = 1;
        let box_width = (self.width / 30).max(2); // Largeur de la boîte
        let right_wall = left_wall + box_width;

        // Créer les murs
        for j in 1..self.height-1 {
            // Mur gauche
            self.wall_init(left_wall, j, true);

            // Mur droit, sauf pour les trous
            if !hole_positions.contains(&j) {
                self.wall_init(right_wall, j, true);
            }
        }

        // Ajouter de la densité dans la boîte
        for i in (left_wall + 1)..right_wall {
            for j in 1..self.height-1 {
                if self.cells[self.idx(i, j)].cell_type == Cell2Type::Fluid {
                    self.cell_init(i, j, flow_velocity, 0.0, density);
                }
            }
        }

        // Initialiser les vitesses aux sorties
        for &j in hole_positions {
            // Initialiser plusieurs cellules à droite de chaque trou pour créer un jet
            for i in right_wall..(right_wall + 5).min(self.width - 1) {
                if self.cells[self.idx(i, j)].cell_type == Cell2Type::Fluid {
                    self.cell_init(i, j, flow_velocity * 2.0, 0.0, density * 0.5);
                }
            }
        }
    }


    /// Configuration pour favoriser la formation de tourbillons de Kármán
    pub fn setup_karman_vortex(&mut self, flow_density: f32, flow_velocity: f32) {
        let radius = (self.width / 6) as f32 - 10.5; // Rayon de l'obstacle
        let center = ((self.width + 2) / 6, (self.height + 2) / 2);

        // Créer l'obstacle circulaire
        self.circle(center.0 as isize, center.1 as isize, radius);

        // Flux horizontal
        for j in 1..self.height-1 {
            for i in 1..15 {
                let idx = self.idx(i, j);
                if self.cells[idx].cell_type == Cell2Type::Fluid {
                    self.cells[idx].density.back = flow_density;
                    self.cells[idx].velocity.back = Vector22::new(flow_velocity, 0.0);
                }
            }
        }

        println!("✔ Configuration de tourbillons de Kármán terminée: obstacle central + flux entrant.");
    }

}

// Extension: Helper pour Vector22
impl std::ops::Index<usize> for Vector22 {
    type Output = f32;

    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            _ => panic!("Index out of bounds for Vector22"),
        }
    }
}

impl std::ops::IndexMut<usize> for Vector22 {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            _ => panic!("Index out of bounds for Vector22"),
        }
    }
}