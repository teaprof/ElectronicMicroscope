The aim of this project is devoted to simulation of the Single Shot Difraction
microscope. The aim of the simulation is to find the best waveguide
geometry which leads to maximum compression the electron cloud.


Write here how single shot difraction microscope operates.


# How to run

1. Clone the project using one of the following commands:
    ```
    git clone https://github.com/teaprof/Optimum-design-of-Single-Shot-Difraction-Microscope.git
    ```
    or
    ```
    git clone git@github.com:teaprof/Optimum-design-of-Single-Shot-Difraction-Microscope.git
    ```
2. Start MATLAB and cd into the project root dir.
3. Initialize MATLAB path with the command:
    ```
    addpath(genpath('.'))
    ```
4. Compile getEBClassicCoder.m to produce getEBClassicCoder_mex.mex file. This can be done by the command:
    ```
        codegen getEBClassicCoder
    ```
5. Run optimMain:
    ```
    optimMain
    ```

# Waveguide solver


```
drawEM
```

![Alt text](figures/TM_00_01_mode_EzEr.png)

![Alt text](figures/TM_00_01_mode_EzHphi.png)

![Alt text](figures/TM_00_01_mode_3d.png)


# Physics

In this project the motion of the electrons affected by laser electromagnetic
field is simulated. The motion is integrated using Runge-Kutta ode45 integrator.
For right-hand side three alternative approaches are implemented.

1. `rightSideNoInteraction` - no interaction between electrons and classical EM field
2. `rightSideClassic` - with instant interaction between electrons and with classical EM field
3. `rightSideRetarded` - With retarded interaction between electrosna and with classical EM field

There is no special mechanism to switch between these rhs function. To select another physics
the appropriate blocks of the code should be commented or uncommented (see [simElectrons.m](./dynamics/sumElectrons.m)).


## rightSideNoInteraction
This is relatively simple algorithm - just take the waveguide solution for EM field and calculate 
$\vec{F_1} = \vec{E}\cdot q +q [\vec{E} \times \vec{B}]$,
where $E$ and $B$ are the electostatic and magnetic fields at the point where the electron is located.


## rightSideClassic

This rhs part is similar to previous one except that $\vec{E}$ and $\vec{B}$ contains additional terms $E_{ij}$ and $B_{ij}$. $\vec{E}_{ij}$ is the Couloumb electric field of charge $i$
at point where charge $j$ is located and $\vec{B}_{ij}$ is the magnetic field defined by Biot-Savart law for the same pair of charged particles.


Since the computational complexity of this rhs call is $O(n^2)$ where $n$ is the number of electrons, the following implementations
are proposed (see [rightSideClassic.m](./dynamics/rightSide/rightSideClassic.m)):
1. `[E, B] = getFieldClassic(q, r, v)` - uses nested for-loops and scalar computations. This implementation is well tested, compared to known analytical solutions and used
    to test other implementations.
2. `[E, B] = getFieldClassicVectorized(q, r, v)` - inner loop is replaced by vectorized computations
3. `[E, B] = getFieldClassicGPU(q, r, v)` is the GPU port (using MATLAB GPU capabilities)
4. `[E, B] = getFieldClassicCoderMex(q, r, v)` can be compiled to binary code (known as mex file) and used by default.

To compile `getFieldClassicCoderMex` just run `codegen getEBClassicCoder` from the MATLAB command line.

Switching between these implementations can be done by commenting or uncommenting blocks in [rightSideClassic.m](./dynamics/rightSide/rightSideClassic.m).

The following picture shows the performance of different implementations. On my ASUS Z170 Pro with i7-6700K and EVGA GTX1080 mex-based implementation
is slightly faster than others.

![benchEBClassic](./figures/benchEBclassic.png)

This picture can be reproduced on your system by typing a command
````
benchEBClassicVariants
````

## rightSideRetarded

Write here what is retarded potential

The difference between this and previous implementation is that this implementation uses retarded potentials while calculating $\vec{E_{ij}}$ and $\vec{B_{ij}}$.
The retarded potential of the charge $i$ can be calculated if we know the previous locations of this particle: we should find time $t$ such that the light emited at this point of time
reaches ... 
The algorithm is slow since it has following steps
1. On each major step of ode45 Runge-Kutta solver the new point is added to the list of previous locations. 
2. Finding roots of the function defined by interpolant

retarded -- very slow


### SimMain

`SimMain` is the function that just simulates the motion of the electrons in the waveguide. It produces avi-file with animation. To run this simulation just type
```
simMain
```

### optimMain

This is a main function which optimizes the waveguide geometry. It uses modified multistart method to solve the optimization problem: first it generates a lot of
points and calculates the target function value at each point, then it takes points with best values and runs the local solver from each of them.

The laser should compress the electron cloud and Couloumb interaction between electrons will grow. This leads to stretching of the cloud and decreasing
of the Couloumb forces. So the length of the electron cloud should oscilate if the waveguide is endless. To find optimal solution one should cut off the
waveguide at the position where the length of the cloud reaches the minimum. 

In this project we simulate the flight of electrons during some constant time. The target function (objective function) is the minimum length
of the electron cloud during the simulation time span.

To run optimization process, just type

```
optimMain
```


### Where results are located

See folders `results` and `figures`



