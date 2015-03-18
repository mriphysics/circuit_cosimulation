# circuit-cosimulation

### A MATLAB implementation of circuit co-simulation

Circuit co-simulation is a method for interactively adapting EM simulations of radio-frequency devices. The basic principal is to replace lumped components by ports when running a simulation. For example a 2-port device with 10 lumped components would be simulated as a 12 port device. The advantage is that from this larger data-set, a reduced 2-port model can be deterimined by assigning any impedance to the lumped elements without running another full EM simulation. This enables the behaviour of the system with different electrical properties to be seen very quickly. The Matlab code in this repo is for performing this combination. The methods are fully described in this [publication](http://dx.doi.org/10.1002/mrm.25504): please cite this in any work resulting from use of our code.

#### Worked Example:
A worked example of an [8-port device](http://dx.doi.org/10.1002/mrm.21294) with 120 lumped elements is provided. Simply run the script **runme.m**. In order to run this code you will need to download binary data for the 128-port S-matrix generated using CST microwave studio from [here](http://bit.ly/1H1sJLO) *(note this is a ~150Mb Matlab mat file)*. 

### Optimisation
It is possible to use co-simulation to compute optimal lumped element values for tuning, matching and decoupling of RF devices. This is again described in our [publication](http://dx.doi.org/10.1002/mrm.25504) and an example is included in **runme.m**.

The optimisation example uses Matlab in-built function `fminsearch`. This is liable to fall into local minima, and a better option is to use the 'SOMA'  algorithm which you can download from the creator's [website](http://www.ft.utb.cz/people/zelinka/soma/).

### Limitations


1. It is not yet possible to attach different networks (i.e. matching networks) to the driving ports and then drive the system externally in the circuit domain. Other contributors are of course welcome to add this functionality: it would simply involve a 2-port network relation between the driving port and the port connecting to the other network.

2. For the lumped element ports, only the capacitance values are considered as part of the optimisation but it would only involve minor editing to also include the resistance or inductance.

3. It is assumed that all the impedances for a given lumped element are in series so please edit this if this is not the case.

4. Currently the full transient LRC behaviour is not fully taken into account so be careful if putting in very high capacitances as the large time constant is not considered; hence leading to inaccurate results.

