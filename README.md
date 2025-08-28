# particle-fortran

## Introduction

Welcome to my particle simulation project. This project is a contiuation of a couple of my previous projects: one a C++ particle simulation I wrote in graduate school, and the other a CUDA fluid simulation. I found that I wanted to explore better paralelization than what was achived in my CUDA fluid simulation project. To modernize and improve parallelization, I have decided to swap to a more-portable compute library with OpenMP and to write multi-device code using MPI. I also wanted to update my previous particle simulation with better output rendering, which I will write in python. The result will be a simulated dust cloud of particles orbiting the sun at a radius between Venus and Earth. The simulation is an N-body problem of asteroid-sized particles which will be allowed to collide with each other via perfectly inelasitic collisions in order to form larger object.

The primary learning objectives of this project are in working with Fortran, OpenMP, and MPI to yeild parallel code. I will then gather some performance benchmarks on CPU, single-GPU, and multi-GPU implementations.

## A CPU particle Simulation

To begin, I want to write a simulation that runs purely in the CPU