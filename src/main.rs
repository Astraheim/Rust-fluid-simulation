use crate::grid::*;
use plotters::prelude::*;
use std::time::Instant;
use conditions::*;
use crate::visualization::run_simulation;

mod grid;
mod conditions;
mod visualization;



fn main() {
    let step: i32 = 0;
    let mut grid = Grid::new();
    let start = Instant::now();

    //Grid::print_grid_walls(&grid);
    // Grid::print_grid_velocity(&grid);
    // Grid::print_grid_density(&grid);

    /*
    for i in 15..=25 {
        grid.wall_init(10, i, true); }  // Upper wall
    for i in 15..=25 {
        grid.wall_init(30, i, true); }  // Lower wall
     */

    let d: usize = 50;
    for i in 0..=N {
        grid.wall_init(d, i, true);
    }  // Upper wall
    for i in 0..=N {
        grid.wall_init(N- d, i, true);
    }  // Lower wall


    Grid::circle(&mut grid, 150, 100, 14.0);


    Grid::cell_init(&mut grid, 15, 20, 2.0,0.01,5000.0);

/*
    for i in 0..100 {
        step+=1;
        println!("Step n°{}", step);
        Grid::vel_step(&mut grid);
    }*/

    println!("Global density {:2}", grid.total_density());

    run_simulation(&mut grid, step);


    let duration = start.elapsed();
    println!("Temps pour vel_step: {:?}", duration);
}
