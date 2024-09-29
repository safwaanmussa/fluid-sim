<!-- PROJECT LOGO -->
<br />
<div align="center">
  <h3 align="center">2D Fluid Simulation</h3>

  <p align="center">
    A real-time fluid dynamics simulation using C and Raylib
    <br />
    <a href="https://github.com/safwaanmussa/fluid-sim"><strong>Explore the code »</strong></a>
    <br />
    <br />
    <a href="https://github.com/safwaanmussa/fluid-sim">View Demo</a>
    ·
    <a href="https://github.com/safwaanmussa/fluid-sim/issues">Report Bug</a>
  </p>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project">About The Project</a></li>
    <li><a href="#built-with">Built With</a></li>
    <li><a href="#getting-started">Getting Started</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#key-features">Key Features</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About The Project

This project is a 2D fluid simulation based on the Navier-Stokes equations, implemented in C using the Raylib graphics library. It demonstrates the application of advanced numerical methods and computational fluid dynamics in a real-time, interactive environment.

<!-- BUILT WITH -->
## Built With

* C
* [Raylib](https://www.raylib.com/)
* GCC Compiler

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running, follow these steps:

1. Ensure you have GCC and Raylib installed on your system. Installation details can be found [here](https://www.raylib.com/). Note that if using the local MinGW GCC compiler bundled with Raylib, the path\to\raylib\w64devkit\bin filepath must be added to the PATH environment variable.
2. Clone the repo
   ```sh
   git clone https://github.com/safwaanmussa/fluid-sim.git
   ```
3. Compile the project
   ```sh
   make
   ```
4. Run the simulation. For Windows this is main.exe in the root directory.

<!-- USAGE -->
## Usage

Once the simulation is running, you can interact with the fluid using your mouse:

- Click and drag to add fluid and velocity
- Observe how the fluid moves and interacts with itself

<!-- KEY FEATURES -->
## Key Features

- Real-time simulation of fluid dynamics
- Interactive fluid manipulation via mouse input
- Implementation of the Navier-Stokes equations for incompressible flow
- Utilization of a multigrid solver for efficient pressure projection
- Visualization of fluid density using color mapping

<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

This project was inspired by and built upon the work of several researchers and educators in the field of computational fluid dynamics:

- Stam, J. (1999). Stable Fluids. In Proceedings of the 26th annual conference on Computer graphics and interactive techniques (SIGGRAPH '99). [Link](https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/ns.pdf)
- Harris, M. (2004). Fast Fluid Dynamics Simulation on the GPU. In GPU Gems: Programming Techniques, Tips, and Tricks for Real-Time Graphics. Addison-Wesley Professional. [Link](https://developer.nvidia.com/gpugems/gpugems/part-vi-beyond-triangles/chapter-38-fast-fluid-dynamics-simulation-gpu)
- The Coding Train. (2017). Coding Challenge #132: Fluid Simulation. [Youtube](https://www.youtube.com/watch?v=alhpH6ECFvQ)


<!-- CONTACT -->
## Contact

Your Name - safwaanmussa28@gmail.com

Project Link: [https://github.com/safwaanmussa/fluid-sim](https://github.com/safwaanmussa/fluid-sim)

<p align="right">(<a href="#readme-top">back to top</a>)</p>
