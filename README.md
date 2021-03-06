# Particle Filter
## Self-Driving Car Engineer Nanodegree Program

Implementation of a two-dimensional particle filter in C++. The particle filter is given a map and some initial localization information (analogous to what a GPS would provide). At each time step the filter receives observation and control data. 

## Visualization of result
Number of particles: 20  

![visualization](pf_visualization.gif)  
Matlab script is `visualization.m`

## Landmark association
I am using simple nearest-neighbor association to assign a landmark to a measurement

## Resampling strategy
Resampling is carried out with replacement. It is simply based on a uniformly generated random number between 0 and 1 - with the combined particle's weights adding up to 1 as well. The list of particles is traversed until the sum of all weights is equal or greater than the random number.

## Running the Code
Once you have this repository on your machine, `cd` into the repository's root directory and run the following commands from the command line:

```
> ./clean.sh
> ./build.sh
> ./run.sh
```
## Inputs to the Particle Filter
You can find the inputs to the particle filter in the `data` directory. 

#### The Map*
`map_data.txt` includes the position of landmarks (in meters) on an arbitrary Cartesian coordinate system. Each row has three columns
1. x position
2. y position
3. landmark id

> * Map data provided by 3D Mapping Solutions GmbH.

#### Control Data
`control_data.txt` contains rows of control data. Each row corresponds to the control data for the corresponding time step. The two columns represent
1. vehicle speed (in meters per second)
2. vehicle yaw rate (in radians per second)

#### Observation Data
The `observation` directory includes around 2000 files. Each file is numbered according to the timestep in which that observation takes place. 

These files contain observation data for all "observable" landmarks. Here observable means the landmark is sufficiently close to the vehicle. Each row in these files corresponds to a single landmark. The two columns represent:
1. x distance to the landmark in meters (right is positive) RELATIVE TO THE VEHICLE. 
2. y distance to the landmark in meters (forward is positive) RELATIVE TO THE VEHICLE.
