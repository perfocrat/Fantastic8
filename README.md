# Equilateral Triangular Infinite-Well Module

## Overview

The Equilateral Triangular Infinite-Well Module is part of the Fantastic8 Physics Project. It models a quantum particle confined within an equilateral triangular potential well.

The module uses an analytical wavefunction model to calculate the energy coefficient and evaluate the wavefunction inside the triangular region. It also generates visual representations of the wavefunction, probability density, and energy level.

## Features

The module validates user inputs, constructs a triangular mesh, calculates six wavevectors, evaluates the wavefunction, and generates three graphical outputs.

The wavefunction graph represents the spatial variation of the wavefunction, while the probability-density graph shows the corresponding probability distribution. The energy-level graph displays the selected quantum state's energy relative to the reference energy.

## Input Parameters

The module accepts the side length of the equilateral triangle and two quantum numbers, m and n.

The side length must be a positive, finite number. The current implementation supports integer quantum numbers satisfying 1 ≤ n ≤ m ≤ 5.

Invalid inputs are rejected before calculations begin.

## Energy Model

The energy of the selected quantum state depends on the triangle's side length and the quantum numbers. The energy coefficient is determined by the expression m² + mn + n².

The module reports the energy coefficient relative to a reference energy that depends on the triangle's side length, the particle's mass, and Planck's reduced constant.

## Graphical Outputs

The module generates three graphical outputs:

1. Wavefunction plot, showing the spatial variation and sign of the wavefunction.
2. Probability-density plot, showing the squared magnitude of the wavefunction.
3. Energy-level plot, showing the energy of the selected quantum state relative to the reference energy.

The graphs are saved as PNG files in the specified output directory.

## Validation and Testing

The module has been tested across all 15 supported quantum-number combinations. Testing has covered input validation, boundary behaviour, wavevector magnitudes, energy consistency, numerical normalization, and graph generation.

Numerical normalization checks showed close agreement with the normalization formula used by the implementation. These tests provide numerical evidence of consistency, but they do not independently prove that the complete physical model is theoretically correct.

## Integration

The module is designed to be called by the Fantastic8 user interface when the equilateral triangular shape is selected. It returns calculation metadata and generates graphical outputs for the selected state.

The current implementation saves the graphs but does not return the sampled wavefunction and probability-density arrays. Additional integration work may be needed if the user interface requires interactive plots or direct access to numerical data.

## Limitations

The current implementation supports only the sine-form wavefunction through its main execution function. The supported quantum-number range is limited to 1 ≤ n ≤ m ≤ 5.

The module's numerical tests do not replace an independent theoretical derivation of the complete triangular-well model. Further theoretical review is required before claiming full physical validation.