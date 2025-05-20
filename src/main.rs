use crate::conditions::*;
use crate::grid::*;
use crate::grid2::*;
use crate::visualization::{run_simulation, run_simulation2};
use std::time::Instant;
pub mod cip_csl4;
pub mod conditions;
pub mod grid;
mod grid2;
mod objects;
mod pressure_computation;
mod pressure_computation2;
pub mod visualization;

fn main() {
    let start = Instant::now();

    if GRID == "1" {
        //Grid::print_grid_walls(&grid);
        // Grid::print_grid_velocity(&grid);
        // Grid::print_grid_density(&grid);

        let mut grid = Grid::new();

        Grid::place_random_objects(&mut grid, 9, 15.0, 20.0);
        //Grid::circle(&mut grid, 75,120,30.9);

        run_simulation(&mut grid, 0);

        println!("Global density {:2}", grid.total_density());
    } else if GRID == "2" {
        let mut grid2 = Grid2::new(N as usize, N as usize, DX);

        //let gravity = Vector22::new(0.0, -0.5); // Gravité vers le bas

        //grid2.add_circle_obstacle(N/2.0, N/2.0, N/8.0, None);
        //grid2.add_circle_obstacle(300.0, 300.0, 50.5, None);

        for i in 0..3 {
            //grid2.wall_init(200, 200+i, true);
            grid2.cell_init(210, 200 + i, 1.0, 0.0, 2000.0);
        }

        println!("Démarrage de la simulation...");
        println!(
            "Utilisation d'une grille {}x{} avec codage Morton",
            grid2.width, grid2.height
        );

        let _zero = Vector22::new(0.0, 0.0);
        println!("Densité avant run_simulation2 : {}", grid2.total_density());

        /*for _ in 0..100{
            grid2.vel_step(DT, _zero);
        }*/

        run_simulation2(&mut grid2);
    }

    let duration = start.elapsed();
    println!(
        " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n Temps pour vel_step: {:?} \n $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$",
        duration
    );
}
