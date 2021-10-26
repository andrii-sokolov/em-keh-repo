# Modelling of Electromechanical Coupling in KEH system.

This software is a part of a PhD project **"Investigation and Optimisation of Kinetic  Energy Harvesters with Nonlinear Electromechanical Coupling Mechanisms"**. It contains both models for the reconstruction of e.m.f. and the electromagnetic mechanical force acting on the permanent magnet that moves close to the coil. In addition it have an optimization procedure for the resonance curve fitting and parameters reconstruction of the system. 

All results were published: 

Modelling and Verification of Nonlinear Electromechanical Coupling in Micro-Scale Kinetic Electromagnetic Energy Harvesters
doi.org/10.1109/TCSI.2019.2938421

## Requirements 

1. gcc >= 7.8 
2. make >= 4.3
3. boost >= 1.65
4. gnuplot >= 5.2

## Usage 

There are two folders named as KEHs used in the corresponding paper (or thesis). To run the modelling of self-consistent electromechanical force model one needs to type:

```
make magnet
```
and 
```
./mag_field
```

To run the optimization procedure: 
```
make 
```

To clean all temporary files:
```
make clear
```
