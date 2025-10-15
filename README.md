# Optimum design of the Single shot diffraction electronic microscope

This project simulates a Single Shot Diffraction microscope. The aim is to 
find the best waveguide geometry which produces the maximum compression of the electron cloud.

This project is inspired by the following work:

1. **Wong, Jie & Fallahi, Arya & Kärtner, Franz.** (2013). *Compact electron acceleration and bunch compression in THz waveguides.* Optics Express. 21. 9792-9806. 10.1364/OE.21.009792. 
([see here](https://www.researchgate.net/publication/236266731\_Compact\_electron\_acceleration\_and\_bunch\_compression\_in\_THz\_waveguides))

A *Single-Shot Diffraction Electronic Microscope* is an advanced imaging system 
designed to capture ultrafast structural information about microscopic objects — 
such as biological samples or nanostructures — using a single ultrashort 
pulse. One of the ways to form a compact electronic pulse is using of THz waveguide to decrease
the thickness of the electron cloud. This project simulates such system.

Figures below shows the motion of electrons in the laser EM field in the waveguide.
![motion](docs/Fig01\_compressed.gif)
![motion](docs/Fig04\_compressed.gif)

The next figure shows the thickness of the electron cloud. At 11 ps it reaches the minimum.
![cloud length along z axis](docs/cloud\_diamz.png)

# How to run

1. Clone the project using one of the following commands:
    ```
    git clone https://github.com/teaprof/Optimum-design-of-Single-Shot-Diffraction-Microscope.git
    ```
    or
    ```
    git clone git@github.com:teaprof/Optimum-design-of-Single-Shot-Diffraction-Microscope.git
    ```
2. Start MATLAB and cd into the project root dir.
3. Initialize MATLAB path with the command:
    ```
    addpath(genpath('.'))
    ```
4. Compile getEBClassicCoder.m to produce getEBClassicCoder\_mex.mex file. This can be done by the command:
    ```
    codegen getEBClassicCoder
    ```
5. Run optimMain:
    ```
    optimMain
    ```


# Algorithm

## Waveguide solver

Finding the electromagnetic field in a waveguide is a classic problem in electromagnetics — 
and essential for designing optical or microwave systems, including diffraction microscopes.

Here we consider the cylindrical (circular) waveguide which consists of two layers: vacuum inner layer and dielectrical outer layer.
In this case we know the analytical solution which gives us $\mathbf{E}$ and $\mathbf{B}$ (or more precisely, $\mathbf{H}$) as functions
of polar coordinates and time. We consider only $TM$ modes ($TE$ modes are not considered).

The following figures shows the different modes of EM-waves in the waveguide.

![TM 01\_01\_mode: Ez, Er](figures/TM\_00\_01\_mode\_EzEr.png)

![TM 01\_01\_mode: Ez, Hphi](figures/TM\_00\_01\_mode\_EzHphi.png)

![TM 01\_01\_mode: 3d](figures/TM\_00\_01\_mode\_3d.png)

These figures can be reproduced by the command 
```
drawEM
```

## Physics

This project simulates the motion of the electrons affected by laser electromagnetic
field. The motion is integrated using Runge-Kutta `ode45` integrator.
For right-hand side (rhs) three alternative approaches are implemented.

1. `rightSideNoInteraction` - no interaction between electrons and classical EM field;
2. `rightSideClassic` - instant interaction between electrons and classical EM field;
3. `rightSideRetarded` - retarded interaction between electrons and classical EM field.

There is no special mechanism to switch between these rhs functions is provided. To select another physics
the appropriate blocks of the code should be commented or uncommented (see [simElectrons.m](./dynamics/sumElectrons.m)).


### rightSideNoInteraction

This is relatively simple algorithm - just take the waveguide solution for EM field and calculate 
$\mathbf{F\_1} = \mathbf{E}\cdot q +q \left[\mathbf{E} \times \mathbf{B}\right]$,
where $\mathbf{E}$ and $\mathbf{B}$ are the electric and magnetic fields at the point where the electron is located.


### rightSideClassic

This rhs part is similar to previous one except that $\mathbf{E}$ and $\mathbf{B}$ contain additional terms $\mathbf{E}\_{other}$ and $\mathbf{B}\_{other}$.

The first one is 

$\mathbf{E}\_{other} = \sum\_{i, i\ne j}\mathbf{E}\_{ij}$ , 

where $\mathbf{E}\_{ij}$ is the Coulomb electric field of charge $i$ at point where charge $j$ is located.

The second additional term is 

$\mathbf{B}\_{other} = \sum\_{i, i\ne j}\mathbf{B}\_{ij}$, 

where $\mathbf{B}\_{ij}$ corresponds to the magnetic field defined by Biot-Savart law for the same pair of charged particles.

This rhs is used by default.

Since the computational complexity of this rhs call is $O(n^2)$ where $n$ is the number of electrons, the following implementations
are proposed (see [rightSideClassic.m](./dynamics/rightSide/rightSideClassic.m)):
1. `[E, B] = getFieldClassic(q, r, v)` - uses nested for-loops and scalar computations. This implementation is well tested, compared to known analytical solutions and used
    to test other implementations;
2. `[E, B] = getFieldClassicVectorized(q, r, v)` - inner loop is replaced by vectorized computations;
3. `[E, B] = getFieldClassicGPU(q, r, v)` is the GPU port (using MATLAB GPU capabilities);
4. `[E, B] = getFieldClassicCoderMex(q, r, v)` can be compiled to binary code (known as mex file) and used by default.

To compile `getFieldClassicCoderMex` just run `codegen getEBClassicCoder` from the MATLAB command line.

Switching between these implementations can be done by commenting or uncommenting blocks in [rightSideClassic.m](./dynamics/rightSide/rightSideClassic.m).

The following picture shows the performance of different implementations on my *ASUS Z170 Pro with i7-6700K and EVGA GTX1080*. It is shown that mex-based implementation is significantly faster than others so this implementation is chosen by default.

![benchEBClassic](./figures/benchEBclassic.png)

This picture can be reproduced on your system by typing a command
````
benchEBClassicVariants
````

### rightSideRetarded

The difference between this and previous implementation is that this implementation uses retarded potentials while calculating $\mathbf{E}\_{ij}$ and $\mathbf{B}\_{ij}$.
They are called retarded because they depend on the state of the sources at an earlier ("retarded") time — the time when the electromagnetic influence left the source and started traveling toward the observation point.

So, to find $\mathbf{E}\_{ij}$ and $\mathbf{B}\_{ij}$ we must find $t\_r$ such that

$t\_r = t - \frac{\left|\mathbf{r}\_i(t\_r) - \mathbf{r}\_j(t)\right|}{c}$

holds, where $\mathbf{r}\_i\left(t\_r\right)$ is the position of the $i$-th charged
particle at time $t\_r$, $\mathbf{r}\_j\left(t\right)$ is the position of $j$-th particle at time $t$, $t$ is the current time, $c$ is the speed of light.

This equation implies that we know $\mathbf{r}\_i\left(t\_r\right)$ function. This function is implemented by the linear interpolation algorithm that uses the list of all previous locations of $i$-th particle.

The resulting algorithm is very slow due to two reasons: on each step of Runge-Kutta it finds the zeros of the function defined as the interpolant and this function should be updated after each major step of Runge-Kutta solver.

By default, the retarded potentials are not used.


# Numerical solvers

## SimMain

`SimMain` is the function that just simulates the motion of the electrons in the waveguide. It produces avi-file with animation. To run this simulation just type
```
simMain
```

## optimMain

This is a main function which optimizes the waveguide geometry. It uses modified multistart method to solve the optimization problem: first it randomly generates a lot of
points and calculates the target function value at each point, then it takes points with best values and runs the local solver from each of them.

The target function (objective function) is the minimum length of the electron cloud during the entire simulation time span.

[comment]: The thickness of the electron cloud should oscillate if the waveguide is endless. To find optimal length of the waveguide one should cut off it at the position where the thickness of the cloud reaches its minimum. So, the length of the waveguide could be excluded from sampled variables. This explains our choice of the target function.

To run optimization process, just type
```
optimMain
```

### Where the results are saved?

See folders `optimresults` (generated by `optimMain`) and `figures` (generated by other scripts).


# End notes

The present document is partially generated by ChatGPT.

The project was in active phase in 2021-2022.

