use crate::grid::*;
use plotters::prelude::*;
use std::time::Instant;
use conditions::*;
use crate::visualization::run_simulation;

mod grid;
mod conditions;
mod visualization;
mod cip_csl4_v7_2d;

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


    Grid::circle(&mut grid, 110, 100, 14.5);

    //Grid::cell_init(&mut grid, 100, 100, 1.0, 1.0, 5000.0);


    println!("Global density {:2}", grid.total_density());

    run_simulation(&mut grid, step);

    let duration = start.elapsed();
    println!(" $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n Temps pour vel_step: {:?} \n $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$", duration);

}

