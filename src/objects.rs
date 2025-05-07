use crate::conditions::*;
use crate::grid::{Grid};
use rand::{Rng};
use rand::distr::Uniform;
use rand::distr::Distribution;
// Last update to rendu_code_tex: 2025-05-02
// last modif: 2025-05-02

/*
FR :
    Placement aléatoire d'obstacles dans la grille

    Génération procédurale : Cette section permet de générer automatiquement des formes solides (murs) dans la grille de simulation.
    Types de formes : Carré, rectangle, triangle plein, triangle en contour (outline) et cercle.
    Contrôle de position : Les formes sont placées à l'intérieur de la grille avec une marge pour éviter les débordements.
    Dessin précis : L'algorithme de Bresenham est utilisé pour tracer des lignes nettes dans le cas de contours triangulaires.
    Objectif : Faciliter la génération de scènes variées pour tester l’interaction fluide/obstacle.

    Améliorations possibles :

    Ajouter une rotation aléatoire pour les formes non circulaires.
    Introduire des formes plus complexes (polygones, ellipses).
    Ajouter une gestion de la superposition ou des collisions lors du placement.
    Permettre des motifs périodiques ou déterministes pour certaines études.

ENG :
    Random Placement of Obstacles in the Grid

    Procedural Generation : This section enables automatic generation of solid shapes (walls) in the fluid simulation grid.
    Shape Types : Square, rectangle, filled triangle, triangle outline, and circle (via external function).
    Position Control : Shapes are placed within grid bounds, avoiding overflow via border margins.
    Precise Drawing : Bresenham's algorithm is used to draw clean line segments for triangle outlines.
    Purpose : Allows generating divers test scenes to study fluid — obstacle interaction.

    Possible Improvements :

    Add random rotation for non-circular shapes.
    Introduce more complex shapes (polygons, ellipses).
    Handle overlaps or collisions between placed objects.
    Enable periodic or deterministic patterns for specific simulations.
*/



impl Grid {
    /// Place randomly different forms in the grid
    pub fn place_random_objects(&mut self, num_objects: usize, min_size: f32, max_size: f32) {
        let mut rng = rand::rng();
        let shape_dist = Uniform::new(0, 4).unwrap();

        // Setup limits for movement
        let border_margin = (max_size * 1.5) as usize;
        let min_pos = border_margin;
        let max_pos_x = (N as usize).saturating_sub(border_margin);
        let max_pos_y = (N as usize).saturating_sub(border_margin);

        for _ in 0..num_objects {
            let shape_type = shape_dist.sample(&mut rng);
            let size = rng.random_range(min_size..=max_size);
            let x = rng.random_range(min_pos..=max_pos_x) as isize;
            let y = rng.random_range(min_pos..=max_pos_y) as isize;

            match shape_type {
                0 => self.circle(x, y, size),
                1 => self.square(x, y, size),
                2 => self.triangle(x, y, size),
                3 => {
                    let width = size;
                    let height = rng.random_range(min_size..=max_size);
                    self.rectangle(x, y, width, height);
                },
                _ => unreachable!(),
            }
        }

        println!("✔ Placed {} random objects in the grid.", num_objects);
    }


    /// Create a square in the grid
    pub fn square(&mut self, center_x: isize, center_y: isize, side_length: f32) {
        let half_side = (side_length / 2.0) as isize;
        let start_x = center_x - half_side;
        let end_x = center_x + half_side;
        let start_y = center_y - half_side;
        let end_y = center_y + half_side;

        for i in start_x..=end_x {
            for j in start_y..=end_y {
                if i >= 0 && i <= N as isize && j >= 0 && j <= N as isize {
                    self.wall_init(j as usize, i as usize, true);
                }
            }
        }
    }


    /// Create a rectangle in the grid
    pub fn rectangle(&mut self, center_x: isize, center_y: isize, width: f32, height: f32) {
        let half_width = (width / 2.0) as isize;
        let half_height = (height / 2.0) as isize;
        let start_x = center_x - half_width;
        let end_x = center_x + half_width;
        let start_y = center_y - half_height;
        let end_y = center_y + half_height;

        for i in start_x..=end_x {
            for j in start_y..=end_y {
                if i >= 0 && i <= N as isize && j >= 0 && j <= N as isize {
                    self.wall_init(j as usize, i as usize, true);
                }
            }
        }
    }


    /// Create a triangle in the grid
    pub fn triangle(&mut self, center_x: isize, center_y: isize, size: f32) {
        let height = size as isize;
        let half_base = (size / 2.0) as isize;

        // Draw equilateral triangle
        for dy in 0..height {
            let width_at_y = (2 * dy * half_base) / height;
            let start_x = center_x - width_at_y;
            let end_x = center_x + width_at_y;

            for x in start_x..=end_x {
                let y = center_y + height - dy; 
                if x >= 0 && x <= N as isize && y >= 0 && y <= N as isize {
                    self.wall_init(y as usize, x as usize, true);
                }
            }
        }
    }


    /// Draw only the outline of a triangle in the grid
    pub fn triangle_outline(&mut self, center_x: isize, center_y: isize, size: f32) {
        let height = size as isize;
        let half_base = (size / 2.0) as isize;

        // Three corners of the triangle
        let top = (center_x, center_y - height/2);
        let bottom_left = (center_x - half_base, center_y + height/2);
        let bottom_right = (center_x + half_base, center_y + height/2);

        // Draw the three sides
        self.draw_line(top.0, top.1, bottom_left.0, bottom_left.1);
        self.draw_line(bottom_left.0, bottom_left.1, bottom_right.0, bottom_right.1);
        self.draw_line(bottom_right.0, bottom_right.1, top.0, top.1);
    }


    /// Another utility using bresenham algorithm to draw lines
    fn draw_line(&mut self, x0: isize, y0: isize, x1: isize, y1: isize) {
        let mut x0 = x0;
        let mut y0 = y0;
        let dx = (x1 - x0).abs();
        let sx = if x0 < x1 { 1 } else { -1 };
        let dy = -(y1 - y0).abs();
        let sy = if y0 < y1 { 1 } else { -1 };
        let mut err = dx + dy;

        loop {
            if x0 >= 0 && x0 <= N as isize && y0 >= 0 && y0 <= N as isize {
                self.wall_init(y0 as usize, x0 as usize, true);
            }

            if x0 == x1 && y0 == y1 {
                break;
            }

            let e2 = 2 * err;
            if e2 >= dy {
                if x0 == x1 {
                    break;
                }
                err += dy;
                x0 += sx;
            }
            if e2 <= dx {
                if y0 == y1 {
                    break;
                }
                err += dx;
                y0 += sy;
            }
        }
    }
}