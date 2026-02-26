# Kettera - A numerical solver for the time-dependent Schrödinger equation

Kettera (from the finnish word ketterä, meaning agile or nimble and from 'ket', the name commonly used for quantum state vectors) is 
a numerical solver for the time-dependent Schrödinger equation (TDSE) written in modern Fortran. It also comes equipped with
scripts for generating plots and animations based on the data produced by the solver. **Kettera is intented to be used on
GNU/Linux systems and no guarantees for functionality are provided for users on other operating systems.**

## Features

Kettera can handle multiple different initial conditions and potentials to simulate many different quantum systems.
Initial conditons can also be read in from an external file. Kettera implements the Crank - Nicolson method as well
as the split-step Fourier method for solving the TDSE. Kettera also comes equipped with two python scripts 'ketplot'
and 'ketanim' which can be used to plot physical quantities as a function of time and for creating animations
that showcase the dynamics of the simulated system. In addition, there is a bash script called 'ketract' bundled with
kettera, which can be used to extract a single frame from an output file to be used as an initial condition for a
simulation.

### Initial conditions

The initial state of the system can be provided in two ways:
1. As normal input parameters describing the type, initial offset and initial momentum of the wave function
2. As an external file that contains the shape of the wave function at time t = 0. This file can be generated
from a previous simulation using the bash script 'ketract' that is bundled with kettera

### Potentials

The user can decide what kind of environment the wavefunction is in by providing the parameters for a potential function.
The user can control the type, offset and strength of the potential function or choose to set the potential to zero
to simulate a free particle.

### Visualization

After a simulation is done, some physical quantities describing the system (eg. total energy, wavefunction norm squared, etc.)
can be plotted as a function of time using the python script 'ketplot'. High quality animations of the simulation can also be
generated using the python script 'ketanim'. Ketanim supports plotting single frames and animating the whole simulation as
a .gif or .mov file.

## External Dependencies

Kettera uses the [FFTW](https://www.fftw.org/) library to calculate Fourier transforms and the python scripts require
[Numpy](https://numpy.org/) as well as [Matplotlib](https://matplotlib.org/) to function. In addition, 'ketract' contains
bash-specific syntax meaning that bash is also an implicit requirement.

## Installation

**NOTE**: This section assumes you have configured the external dependencies listed in the previous section to work on your system.

TBD

## Usage

To invoke Kettera, simply type `kettera` into a terminal. This should produce some simple usage information about kettera.
Providing the option `--help` will provide more detailed usafe information about Kettera and providing the option `--version`
will output Kettera's current version.

### Running simulations

TBD

### Using external input files and ketract

TBD

### Using ketanim and ketplot

TBS

## License

Kettera is provided under [version 3 of the GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Author

Kettera is written by Juuso Kaarela in 2026


