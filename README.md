# Particle Filter: Kidnapped Vehicle Project

[//]: # (Image References)
[image1]: ./doc/intro.png  "intro"
[image2]: ./doc/Filter_algo.png  "algo"
[image3]: ./doc/result.png  "result"

## Project Introduction
A robot has been kidnapped and transported to a new location! Luckily it has a map of this location, a (noisy) GPS estimate of its initial location, and lots of (noisy) sensor and control data. 

In this project, I implemented a 2 dimensional particle filter in C++. The particle filter is given a map and some initial localization information (analogous to what a GPS would provide). At each time step the filter will also get observation and control data, and will update the filter based on this information.  In the project simulator, the result looks like this:

![alt text][image1]

## [Rubric Points](https://review.udacity.com/#!/rubrics/747/view)

Here I will consider the rubric points individually.  

### Accuracy: Does your particle filter localize the vehicle to within the desired accuracy?

Yes, the error values never exceed the specified bounds.  The final error values are:

| Measurement |  Error  |
|:-----------:|:-------:|
|      x      |  0.113  |
|      y      |  0.109  |
|     yaw     |  0.004  |

### Performance: Does your particle run within the specified time of 100 seconds?

Yes, the run completes in about 50 seconds.

### General: Does your code use a particle filter to localize the robot?

Yes, the localization is visualized in the simulator as a blue circle.  The green lines indicate sensor measurements to the various landmarks, and the blue lines indicate the observations from the best predicted particle in the filter. 

## Running the Code
This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

The main program can be built and run by doing the following from the project top directory.

1. mkdir build && cd build
2. cmake .. && make
3. ./particle_filter

## Implementing the Particle Filter

The overall particle filter algorithm looks like this:

![alt text][image2]

The sections below will go through each of the 4 steps in the algorithm, which can be seen above in green boxes. 

### Initialization 

In the initialization step (see the "init()" function), the "ParticleFilter" class takes in the initial GPS measurements, and sets the initial states of all particles to those values +/- some Gaussian noise using the provided standard deviation and the "normal_distribution" class from the C++ standard library.  I chose to use 100 particles as this gave me good accuracy with reasonable performance.  

### Prediction Step 

In the prediction step (see the "prediction()" function), the "ParticleFilter" class takes in the velocity and yaw rate of the vehicle, and uses this to predict the new x, y, and theta values using trigonometry and our "bicycle" motion model.  It again adds some Gaussian noise to the measurements using the provided standard deviation and the "normal_distribution" class.

### Update Step 

The update step (see the "updateWeights() function) is by far the most complex step in the process.  For each particle, the algorithm preforms the following:
1. Gather a list of predicted measurements by filtering the list of all landmark objects to only include those within the provide sensor range.
2. Transform the input set of observations from the car from the Vehicle's coordinate system into the Map's coordinate system.  
3. Associate each of the observations with a predicted measurement using the nearest neighbor algorithm
4. Calculate the new weight of the particle to be the product of each observations [Multivariate-Gaussian probability density](https://en.wikipedia.org/wiki/Multivariate_normal_distribution).

### Resample 

The resampling step can be found in the "resample()" function of the "ParticleFilter" class.  In the course materials on the resample step, the instructor used a "resampling wheel" to generate the new set of particles from the assigned weights.  However, as was suggested in the provided source code, the C++ standard library has a "discrete_distribution" class that simplifies this operation.  Thus, all I needed to do was extract the weights of each particle into a separate vector, and pass that vector into a "discrete_distribution" class, which I then use to select indexes in the previous list of particles.

## Results

Here's a [link to a video recording of my final result](./project_recording.mp4).  

Here is a snapshot at the end confirming the success criteria has been met:

![alt text][image3]

## Discussion

In the course materials, Sebastian described particle filters as "easy to code".  In the end, the concepts were much easier to understand than the previously discussed Kalman filters (mostly due to the lack of linear algebra), but the amount of guidance provided for creating the particle filter was also quite sparse.  The update step in particular was lacking in detail and required me to piece together details from the python and C++ based lessons.  I was most confused initial about the sensor range input, but the "Explanation of Project Code" section (which I initially skipped as the shown source code was outdated) finally gave some clues about what this was used for.  This helped clarify how to generate the initial set of "predicted" measurements. After that, the rest was just implementing equations for transforming the coordinate spaces, and the Multivariate-Gaussian probability density.  

The amount of loops, and nested loops over all the data made it clear why this needed to be done in C++.  I'm sure there were some optimizations possible by using more efficient data types, but with 100 particles, I saw good results with mostly smooth runtime performance. With more time, I would have liked to experiment with more particles and see if accuracy and performance could be improved with additional tweaks.  