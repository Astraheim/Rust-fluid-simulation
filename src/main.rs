use crate::grid::*;
use plotters::prelude::*;
use std::time::Instant;
use crate::visualization::run_simulation;

pub mod grid;
pub mod conditions;
pub mod visualization;
pub mod cip_csl4;
mod pressure_computation;
mod objects;

fn main() {
    let mut grid = Grid::new();
    let start = Instant::now();

    //Grid::print_grid_walls(&grid);
    // Grid::print_grid_velocity(&grid);
    // Grid::print_grid_density(&grid);


    //Grid::place_random_objects(&mut grid, 9, 15.0, 40.0);
    run_simulation(&mut grid, 0);
    

    println!("Global density {:2}", grid.total_density());

    let duration = start.elapsed();
    println!(" $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n Temps pour vel_step: {:?} \n $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$", duration);
}

