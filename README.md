# NeuronBuilder2

Joe Snider, 10/29/2017

This is the c++ source for a an algorithm to grow neurons.

TODO: Summary of algorithm.

## Installation

Compile with the handy-dandy comp2.sh script, or just use gcc. It requires boost and gsl. Otherwise, it is pretty standard C++ (99).

## Running

Run it as

.\neuronbuilder2 <remove radius> 

and it will create an arbor and dump segments to stdout (pipe it wherever).

The remove radius parameter controls the scaling exponent of the arbor: If Rr is big, then the resulting arbors are pitiful little lines (dimension 1). If Rr is small (i.e. 0), then the arbors fill up the available space (dimension 3). In between, the arbors are statistically self-similar. If you create a few hundred arbors and log-log plot their volume (e.g. convex hull volume) versus total length, then you will get a straight(ish) line with a slope that is the fractal dimension. The fractal
dimension is non-trivial around values of Rr=1.


