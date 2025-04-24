use crate::grid::*;
use plotters::prelude::*;
use std::time::Instant;
use conditions::*;
use crate::visualization::run_simulation;

pub mod grid;
pub mod conditions;
pub mod visualization;
pub mod cip_csl4;


fn main() {
    let step: i32 = 0;
    let mut grid = Grid::new();
    let start = Instant::now();

    //Grid::print_grid_walls(&grid);
    // Grid::print_grid_velocity(&grid);
    // Grid::print_grid_density(&grid);


    /*let d: usize = 60;
    for i in 1..=N as usize {
        grid.wall_init(d, i, true);
    }  // Upper wall
    for i in 1..=N as usize {
        grid.wall_init(N as usize - d, i, true);
    }  // Lower wall
*/


    //grid.setup_karman_vortex(); // <<< ici on prépare l'expérience

    Grid::cell_init(&mut grid, N as usize/5, N as usize/2, 0.0, 0.0, 500.0);
    Grid::circle(&mut grid, (N as usize / 4) as isize, (N as usize / 2) as isize, (N as usize / 10) as f32 -6.5);
    run_simulation(&mut grid, 0);








    println!("Global density {:2}", grid.total_density());

    let duration = start.elapsed();
    println!(" $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n Temps pour vel_step: {:?} \n $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$", duration);
}

