'''
  pqc

'''

import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Pauli

# ansatz
from sub.ansatz import pqc_ansatz_PNC
from sub.ansatz import pqc_ansatz_HE

# expval
from sub.utils import expval

# qiskit primitives
from qiskit.primitives import Estimator


# Quantum surrogate circuit using particle number conserving ansatz
def QNNcircuit_PNC(num_qubits, X, reps_pqc, reps, params):
    X_min = 0.45
    X_max = 3.0

    qc = QuantumCircuit(num_qubits)

    ## normalization of X
    _X = 2.0 * (X - X_min) / (X_max - X_min) + 0.5

    sigma_x = [0, 0, _X]

    # HF initial state
    for i in range(num_qubits // 2):
        qc.x(2 * i)

    block_params_num = (num_qubits - 1) * reps_pqc - 1

    ## First Layer -------------------------------------------
    layer_num = 0
    params0 = [sigma_x[2]] + list(
        params[int(block_params_num * layer_num):int(block_params_num *
                                                     (layer_num + 1))])
    pqc = pqc_ansatz_PNC(qc, params0, reps_pqc)
    ## -------------------------------------------------------

    op0 = Pauli('IIIZ')
    sigma_z = expval(op0, pqc)

    ## Second Layer -------------------------------------------
    qc2 = QuantumCircuit(num_qubits)

    # HF initial state
    for i in range(num_qubits // 2):
        qc2.x(2 * i)

    layer_num = 1
    params1 = [sigma_z[0]] + list(
        params[int(block_params_num * layer_num):int(block_params_num *
                                                     (layer_num + 1))])
    pqc2 = pqc_ansatz_PNC(qc2, params1, reps_pqc)
    ## -------------------------------------------------------

    return pqc2


# Quantum surrogate circuit using particle number conserving ansatz
# witout intermediate measurements
def QNNcircuit_PNC_linear(num_qubits, X, reps_pqc, reps, params):
    X_min = 0.45
    X_max = 3.0

    qc = QuantumCircuit(num_qubits)

    ## normalization of X
    _X = 2.0 * (X - X_min) / (X_max - X_min) + 0.5

    sigma_x = [0, 0, _X]

    # HF initial state
    for i in range(num_qubits // 2):
        qc.x(2 * i)

    block_params_num = ((num_qubits - 1) * reps_pqc - 1) * reps + (reps - 1)

    ## First Layer -------------------------------------------
    layer_num = 0
    params0 = [sigma_x[2]] + list(
        params[int(block_params_num * layer_num):int(block_params_num *
                                                     (layer_num + 1))])
    pqc = pqc_ansatz_PNC(qc, params0, reps_pqc * reps)
    ## -------------------------------------------------------

    return pqc


# model for Symmetry Preserving
def QNNcircuit_PNC_with_g(num_qubits, X, reps_pqc, reps, params, weight=1.0):
    X_min = 0.45
    X_max = 3.0

    _X = 2.0 * (X - X_min) / (X_max - X_min) + 0.5

    block_params_num = (num_qubits - 1) * reps_pqc - 1

    sigma_z = funcg(_X)

    ## Second Layer -------------------------------------------
    qc2 = QuantumCircuit(num_qubits)

    # HF initial state
    for i in range(num_qubits // 2):
        qc2.x(2 * i)

    layer_num = 0
    params1 = [sigma_z] + list(
        params[int(block_params_num * layer_num):int(block_params_num *
                                                     (layer_num + 1))])
    pqc2 = pqc_ansatz_PNC(qc2, params1, reps_pqc)
    ## -------------------------------------------------------

    return pqc2


def QNNcircuit_PNC_with_h(num_qubits, X, reps_pqc, reps, params, weight=1.0):
    X_min = 0.45
    X_max = 3.0

    _X = 2.0 * (X - X_min) / (X_max - X_min) + 0.5

    block_params_num = (num_qubits - 1) * reps_pqc - 1

    sigma_z = funch(_X)

    ## Second Layer -------------------------------------------
    qc2 = QuantumCircuit(num_qubits)

    # HF initial state
    for i in range(num_qubits // 2):
        qc2.x(2 * i)

    layer_num = 0
    params1 = [sigma_z] + list(
        params[int(block_params_num * layer_num):int(block_params_num *
                                                     (layer_num + 1))])
    pqc2 = pqc_ansatz_PNC(qc2, params1, reps_pqc)
    ## -------------------------------------------------------

    return pqc2


def funcg(x):

    c0 = 0.1
    c1 = 0.65
    c2 = 0.8
    c3 = 0.65
    c4 = 2.0

    return c0 + c1 * np.sin(c2 * (x - c3))**c4


def funch(x):

    c0 = 0.0
    c1 = 0.5
    c2 = 0.8
    c3 = 0.4
    c4 = 2.0

    return c0 + c1 * np.sin(c2 * (x - c3))**c4


# Quantum surrogate circuit using particle number conserving ansatz
def QNNcircuit_PNC_H3(num_qubits, X, reps_pqc, reps, params):

    qc = QuantumCircuit(num_qubits)

    ## normalization of X
    sigma_x = [X[0], X[1], X[2], 0]

    # HF initial state
    for i in range(num_qubits // 2):
        qc.x(2 * i)

    block_params_num = (num_qubits - 1) * reps_pqc - 1

    ## First Layer -------------------------------------------
    layer_num = 0
    params0 = [sigma_x[0]] + list(
        params[int(block_params_num * layer_num):int(block_params_num *
                                                     (layer_num + 1))])
    pqc = pqc_ansatz_PNC(qc, params0, reps_pqc)
    ## -------------------------------------------------------

    op0 = Pauli('IIIIIZ')
    observable = [op0]

    estimator = Estimator()
    job = estimator.run([pqc for i in range(len(observable))], observable)
    result = job.result()

    sigma_z = result.values

    ## Second Layer -------------------------------------------
    qc2 = QuantumCircuit(num_qubits)

    # HF initial state
    for i in range(num_qubits // 2):
        qc2.x(2 * i)

    layer_num = 1
    params1 = [sigma_z[0]] + list(
        params[int(block_params_num * layer_num):int(block_params_num *
                                                     (layer_num + 1))])
    pqc2 = pqc_ansatz_PNC(qc2, params1, reps_pqc)
    ## -------------------------------------------------------

    return pqc2


# model
def QNNcircuit_HE(num_qubits, X, reps_pqc, reps, params, weight=1.0):
    '''
    Hardware Efficient ansatz (Xia model)
    ref) R.Xia and S. Kais, Hybrid quantum-classical neural network
    for calculating ground state energies of molecules,” Entropy 22, 828 (2020)
    '''
    X_min = 0.45
    X_max = 3.0

    qc = QuantumCircuit(num_qubits)

    ## normalization of X
    _X = 2.0 * (X - X_min) / (X_max - X_min) + 0.5

    sigma_x = [0, 0, _X]

    ## encoding
    for i in range(num_qubits):
        qc.h(i)
        qc.ry(sigma_x[2], i)  # one-encoding

    block_params_num = num_qubits * reps_pqc  # HE PQC

    ## First Layer -------------------------------------------
    layer_num = 0
    pqc = pqc_ansatz_HE(
        qc, params[int(block_params_num * layer_num):int(block_params_num *
                                                         (layer_num + 1))],
        reps_pqc)
    ## -------------------------------------------------------

    op0 = Pauli('IIIZ')
    op1 = Pauli('IIZI')
    op2 = Pauli('IZII')
    op3 = Pauli('ZIII')
    observable = [op0, op1, op2, op3]

    estimator = Estimator()
    job = estimator.run([pqc for i in range(len(observable))], observable)
    result = job.result()

    sigma_z = result.values * weight

    ## Second Layer -------------------------------------------
    qc2 = QuantumCircuit(num_qubits)

    ## encoding
    for i in range(num_qubits):
        qc2.h(i)
        qc2.ry(sigma_z[i], i)

    layer_num = 1
    pqc2 = pqc_ansatz_HE(
        qc2, params[int(block_params_num * layer_num):int(block_params_num *
                                                          (layer_num + 1))],
        reps_pqc)
    ## -------------------------------------------------------

    return pqc2


def QNNcircuit_HE_H3(num_qubits, X, reps_pqc, reps, params, weight=1.0):
    '''
    Hardware Efficient ansatz (Xia model)
    ref) R.Xia and S. Kais, Hybrid quantum-classical neural network
    for calculating ground state energies of molecules,” Entropy 22, 828 (2020) 
    '''
    qc = QuantumCircuit(num_qubits)

    ## normalization of X
    sigma_x = [X[0], X[1], X[2], 0]

    ## encoding
    for i in range(num_qubits):
        qc.h(i)
        qc.ry(sigma_x[0], i)  # one-encoding

    block_params_num = num_qubits * reps_pqc  # HE PQC

    ## First Layer -------------------------------------------
    layer_num = 0
    pqc = pqc_ansatz_HE(
        qc, params[int(block_params_num * layer_num):int(block_params_num *
                                                         (layer_num + 1))],
        reps_pqc)
    ## -------------------------------------------------------

    op0 = Pauli('IIIIIZ')
    op1 = Pauli('IIIIZI')
    op2 = Pauli('IIIZII')
    op3 = Pauli('IIZIII')
    op4 = Pauli('IZIIII')
    op5 = Pauli('ZIIIII')
    observable = [op0, op1, op2, op3, op4, op5]

    estimator = Estimator()
    job = estimator.run([pqc for i in range(len(observable))], observable)
    result = job.result()

    sigma_z = result.values * weight

    ## Second Layer -------------------------------------------
    qc2 = QuantumCircuit(num_qubits)

    ## encoding
    for i in range(num_qubits):
        qc2.h(i)
        qc2.ry(sigma_z[i], i)

    layer_num = 1
    pqc2 = pqc_ansatz_HE(
        qc2, params[int(block_params_num * layer_num):int(block_params_num *
                                                          (layer_num + 1))],
        reps_pqc)
    ## -------------------------------------------------------

    return pqc2
