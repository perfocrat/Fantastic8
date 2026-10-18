# Triangle Backend Testing and Validation

## 1. Overview

This module tests and validates the equilateral triangular infinite-well backend used in the Fantastic8 Physics Project.

The backend calculates quantum wavefunctions, probability densities, and energy levels for an equilateral triangular infinite potential well. The testing module checks whether the implementation behaves as expected for valid and invalid inputs, whether the wavefunction satisfies the boundary conditions, whether the normalization is correct, and whether the expected graphical outputs are generated.

The testing framework includes automated tests and numerical validation using integration.

## 2. Objectives of Testing

The tests are designed to verify the following:

- Input validation for side length and quantum numbers.
- Correct handling of valid and invalid quantum-number combinations.
- Consistency of the wavevector calculations.
- Satisfaction of the triangular well's boundary conditions.
- Normalization of the wavefunction.
- Consistency of calculated energy levels.
- Successful generation of wavefunction, probability-density, and energy-level plots.
- Numerical verification of the normalization constant using independent numerical integration.

## 3. Files Used

The testing module uses the following files:

- `trianglebackend.py` — Main triangular infinite-well backend containing the mathematical calculations and plotting functions.
- `test_triangle.py` — Automated tests for the backend's input handling, calculations, and output generation.
- `verify_triangle_normalization.py` — Numerical validation of the wavefunction normalization constant using Gauss–Legendre quadrature.
- `README.md` — Documentation explaining the testing procedure and expected results.

The backend requires NumPy and Matplotlib. The normalization verification script also uses NumPy's numerical functionality.

## 4. Requirements

Python 3 and the following Python libraries are required:

- NumPy
- Matplotlib

The `math` and `pathlib` modules are part of Python's standard library.

Install the required libraries using:

`python3 -m pip install numpy matplotlib`

If the project already has a virtual environment, activate it before installing the dependencies.

On macOS, run the following command from the project directory:

`source .venv/bin/activate`

## 5. How to Run the Tests

### Step 1: Open the project directory

Open Terminal and navigate to the directory containing the testing files and `trianglebackend.py`.

### Step 2: Activate the virtual environment

If the project uses a virtual environment named `.venv`, activate it using:

`source .venv/bin/activate`

### Step 3: Run the automated tests

Execute:

`python3 test_triangle.py`

This runs the automated tests implemented in the testing script. The tests check the supported inputs, invalid inputs, calculation consistency, and graph generation according to the assertions implemented in the file.

A successful run should finish without unexpected exceptions or failed assertions. The exact printed output depends on the messages included in the testing script.

### Step 4: Run the normalization verification

Execute:

`python3 verify_triangle_normalization.py`

This script numerically integrates the squared, unnormalized wavefunction over the triangular domain and uses the result to infer the normalization constant. It then compares the inferred value with the normalization constant used by the backend.

The numerical integration uses Gauss–Legendre quadrature. Increasing the quadrature order improves the numerical approximation, subject to floating-point precision and implementation details.

### Step 5: Check the generated plots

After running the backend or the corresponding graph-generation tests, inspect the output directory specified in the script. The backend's default output directory is `generated_plots`.

## 6. Tests Performed

### 6.1 Input Validation

The backend validates the side length and quantum numbers before performing calculations.

The tests cover valid positive finite side lengths, supported integer quantum numbers, invalid values, and quantum-number combinations that violate the backend's restrictions.

The supported quantum-number range is:

- \(1 \leq n \leq m \leq 5\)

This gives 15 valid combinations.

The backend rejects combinations where \(m<n\), as well as invalid quantum numbers and invalid side lengths. Boolean values are also rejected where integer quantum numbers are expected.

**Expected result:** Valid inputs are accepted, while invalid inputs are rejected through input validation.

### 6.2 Wavevector Consistency

The wavevector function calculates the wavevectors used to construct the triangular wavefunction.

The tests check the consistency of the calculated components with the mathematical expressions implemented in the backend.

**Expected result:** The wavevectors follow the implemented equations and maintain the expected relationships between their components.

### 6.3 Boundary Condition Verification

For an infinite potential well, the wavefunction must vanish along the boundaries of the allowed region.

Numerical tests evaluate the wavefunction at points on the triangular boundaries and measure the residual values.

**Expected result:** The boundary residuals are approximately \(10^{-15}\), which is consistent with floating-point numerical precision for the tested points.

### 6.4 Wavefunction Normalization

A normalized wavefunction must satisfy:

\[
\int_{\text{triangle}} |\psi(x,y)|^2\,dA = 1
\]

The normalization is checked numerically over the triangular domain for all 15 supported quantum-number combinations.

The verification script also independently estimates the normalization constant by integrating the squared raw wavefunction and comparing the inferred constant against the value used by the backend.

**Expected result:** The numerical normalization integral approaches 1, and the independently inferred normalization constant agrees with the implemented value to within numerical precision.

### 6.5 Energy-Level Consistency

The backend calculates energy levels using the triangular well's implemented energy model.

The tests check consistency between the calculated energy and the quantum-number expression used by the backend.

**Expected result:** The energy calculations are consistent with the implemented analytical formula.

This is an internal consistency check, not an independent proof that the complete physical model is correct.

### 6.6 Graph Generation

The backend generates plots for the selected quantum state.

The testing process checks that the backend can execute and generate the expected graphical outputs.

**Expected result:** The requested plot files are generated successfully when the backend completes without errors.

## 7. Expected Output Files

For a test input with side length \(a=5\), quantum numbers \(m=2\) and \(n=1\), the backend is expected to generate the following files in its output directory:

- `wavefunction_m2_n1_s.png` — Visualizes the spatial wavefunction for the selected sine state.
- `probability_density_m2_n1_s.png` — Visualizes the probability density, \(|\psi(x,y)|^2\).
- `energy_level_m2_n1.png` — Visualizes the calculated energy level.

The exact output directory depends on the `output_dir` argument supplied to `run_triangle`. By default, the backend uses `generated_plots`.

The wavefunction plot represents the signed spatial amplitude, whereas the probability-density plot represents the non-negative squared magnitude. The energy-level plot shows the energy associated with the selected quantum numbers according to the implemented model.

## 8. Numerical Validation Results

The following results were obtained during the validation of the current implementation:

- All 15 supported quantum-number combinations were tested for numerical normalization.
- Ten expected invalid-input cases were tested.
- No unexpected failures were observed in the automated test batch.
- Boundary-condition residuals were approximately \(10^{-15}\).
- Numerical integration on a fine mesh produced normalization values ranging from approximately 0.99974301 to 0.99998972.
- Gauss–Legendre quadrature at order 160 produced a maximum relative difference of approximately \(3.64 \times 10^{-15}\) between the independently inferred normalization constant and the implemented normalization constant.

These values describe the tests performed on the current implementation. They are not guarantees for every possible input, numerical environment, or future code modification.

## 9. Understanding Test Failures

If a test fails, check the following:

- Confirm that `trianglebackend.py` is available in the expected directory.
- Ensure NumPy and Matplotlib are installed in the active Python environment.
- Verify that the quantum numbers satisfy the backend's input restrictions.
- Check that the output directory is writable.
- Confirm that the expected plot filenames match those generated by the current backend.
- If a normalization test fails, inspect the integration domain, quadrature order, wavevector expressions, normalization constant, and numerical tolerances.
- If a boundary test fails, check the wavefunction expression and the coordinates used to evaluate the triangular boundaries.

## 10. Limitations

The backend uses analytical expressions for the triangular infinite-well wavefunction and energy rather than numerically solving the Schrödinger equation from scratch.

The numerical tests verify implementation consistency and selected mathematical properties. In particular, energy consistency tests reuse the backend's analytical model and therefore do not independently establish the physical correctness of that model.

The normalization verification independently evaluates the integral numerically, but it still uses the wavefunction and wavevector expressions supplied by the backend. A complete theoretical validation would additionally require a rigorous derivation of the wavefunction, boundary conditions, normalization, and energy spectrum.

The test results apply to the implementation and cases tested. They should not be interpreted as proof that every aspect of the physical model has been independently verified.

## 11. Integration with Fantastic8

This testing module supports the integration of the triangular backend into the Fantastic8 application.

The frontend can supply the selected shape, quantum numbers, and shape-specific parameters. The integration layer can then call the corresponding backend and display its calculated results and plots.

Before integration, the backend should be tested independently to confirm that valid inputs are handled correctly, invalid inputs are rejected, calculations complete successfully, and the expected plots are generated.

For integration into the shared application, the backend's return values and plot paths should be exposed consistently so the frontend can retrieve the numerical results and display the generated plots.

## 12. Conclusion

The testing module provides automated checks and numerical validation for the equilateral triangular infinite-well backend. It covers input validation, wavevector consistency, boundary conditions, normalization, energy calculations, and graph generation.

The results support the numerical consistency of the current implementation for the tested cases while identifying the distinction between software validation and complete theoretical verification of the underlying quantum-mechanical model.
