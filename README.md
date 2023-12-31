# Pose-From-Conics

## Fabricated data
``generateConicData()`` gives the randomly generated data that could be feed into the methods that create constraints.

## Methods
This folder contains helper functions and several functions that creates the constrains.
- createPointConstraints
- createConicConstraints
- createPolarityConstraints
- createTangencyConstraints
- degPtsConstraints: Used to find the dimension and degree for the point problems.

## Experiements
To run experiments with points problems:

In ``findDimDeg.m2``, Call the ``degPtsConstraints(R,L)`` to find the dimension and degree of the point problems created. ``L`` is a list {n1,n2} that contains the number of points on the conics, where n1 is the number of points on the first image conic, and n2 is the number of points on the second image conic.

## Todo

### New Methods:
- Generic Number of Solutions:
    - degTanConstraints: similar to degPtsConstraints (input: R (polynomial ring), L (list of number of tangents), output: dimension and degree).
    - degPtsTanConstraints: compute the dimension and degree for the mixed point-tangent problems (input: R, L(list of number of points), S (list of number of tangents)).
    - degSpecialConstraints: functions for any other possible minimal problems.
- Homotopy Methods:
    - ``conicMonodromy.m2`` (detailed codes for implementing monodromySolve)
  
### Experiments
- Find the degrees from different tangent problems and mixed point-tangent problems.
- In ``homotopyExperiment.m2``, apply the homotopy method to the minimal problems, report the number of solutions to each of the problems, mention possible symmetries.
- Using ``GAP``, compute the Galois group for result from the homotopy method, try to find out any non-trivial groups.
- Experiments with noisy inputs.