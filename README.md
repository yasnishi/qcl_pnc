## A quantum circuit learning model using particle-number-conserving state
This repository contains the source code used to produce the results presented in "Applications of quantum circuit learning model using particle-number-conserving state on quantum chemical calculations". The source code can reproduce the inferred potential energy surfaces of H2 molecule shown in Figs 5, 7 and 11 and H3 molecule in Fig. 8 of the above-mentioned paper.

## Requirements

The following packages is needed:

```bash
pip install qiskit
pip install pyscf
pip install openfermion
pip install tqdm
```

We have confirmed that this code works with qiskit==0.45.0.

## How to use

To use the scripts, simply set the input data in notebooks. However, only n=4 or 6 (number of qubits) trained quantum circuit model  with D=4 (PQC depth) are available. Optimized parameters are stored in the params directory.

## License

Please see the license file.


