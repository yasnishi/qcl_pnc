'''
  utilities

'''

# qiskit
from qiskit import QuantumCircuit
from qiskit.primitives import Estimator
from qiskit.quantum_info import SparsePauliOp

# openfermion
from openfermion.chem import MolecularData
from openfermion.transforms import get_fermion_operator, jordan_wigner, bravyi_kitaev
from openfermion.linalg import get_sparse_operator
from openfermion.ops import FermionOperator
from openfermionpyscf import run_pyscf, generate_molecular_hamiltonian

import numpy as np

#unit
from ase.units import Bohr


class Mole():
    def __init__(self,
                 basis='sto3g',
                 multiplicity=None,
                 charge=None,
                 geometry=None):

        self.basis = basis
        self.multiplicity = multiplicity
        self.charge = charge
        self.geometry = geometry

    def __str__(self):
        return str({
            "basis": self.basis,
            "multiplicity": self.multiplicity,
            "charge": self.charge,
            "geometry": self.geometry,
        })

    def build_structure_for_openfermion(self, atoms, coords):
        self.atoms = atoms
        self.x = coords
        n_atom = len(self.atoms)
        geoms = []
        for i in range(n_atom):
            geoms.append([self.atoms[i], tuple(self.x[i * 3:i * 3 + 3])])

        self.geometry = geoms


def get_qubit_hamiltonian_of_H2(X):

    mol = Mole()
    mol.multiplicity = 1
    mol.charge = 0

    atoms = ['H', 'H']

    H1 = [0, 0, 0]
    H2 = [0, 0, X]
    x = H1 + H2

    mol.build_structure_for_openfermion(atoms, x)
    qubit_hamiltonian, num_qubits = H_from_openfermion(mol)

    return qubit_hamiltonian, mol


def get_qubit_hamiltonian_of_H3(X, R):

    mol = Mole()
    mol.multiplicity = 2
    mol.charge = 0

    atoms = ['H', 'H', 'H']

    H1 = [0, 0, 0]
    H2 = [0, 0, X]
    H3 = [0, 0, R * Bohr]
    x = H1 + H2 + H3

    mol.build_structure_for_openfermion(atoms, x)
    qubit_hamiltonian, num_qubits = H_from_openfermion(mol)

    return qubit_hamiltonian, mol


def H_from_openfermion(mol):

    atoms = mol.atoms
    #define constants
    basis = mol.basis  #"sto-3g"  #basis set
    multiplicity = mol.multiplicity  #spin multiplicity
    charge = mol.charge  #total charge for the molecule
    geometry = mol.geometry

    molecule = MolecularData(geometry,
                             basis,
                             multiplicity,
                             charge,
                             filename='./mol')
    molecule = run_pyscf(molecule, run_scf=True, run_fci=False)
    hamil = jordan_wigner(
        get_fermion_operator(molecule.get_molecular_hamiltonian()))

    pauli_str = list(hamil.terms.keys())
    pauli_coeffs = list(hamil.terms.values())

    # num_qubits
    temp = []
    for i in range(len(hamil.terms)):
        temp.append(np.max(len(pauli_str[i])))
    num_qubit = np.max(temp)

    # qubit hamiltonian

    hamil_str = []
    hamil_coeffs = []

    for i in range(len(hamil.terms)):
        hamil_coeffs.append(pauli_coeffs[i])

    for i in range(len(hamil.terms)):
        if pauli_str[i] == ():
            op = ''
            for j in range(num_qubit):
                op = op + 'I'
            hamil_str.append(op)
        else:
            temp_op_index = []
            temp_op = []
            for j in range(len(pauli_str[i])):
                temp_op_index.append(pauli_str[i][j][0])
                temp_op.append(pauli_str[i][j][1])

            temp_op_index2 = [0 for i in range(num_qubit)]
            for k in range(num_qubit):
                if k in temp_op_index:
                    for k2 in range(len(temp_op_index)):
                        if k == temp_op_index[k2]:
                            temp_op_index2[k] = k2
                else:
                    temp_op_index2[k] = num_qubit

            op = ''
            for k in range(num_qubit):
                if temp_op_index2[k] == num_qubit:
                    op = op + 'I'
                else:
                    op = op + temp_op[temp_op_index2[k]]

            hamil_str.append(op)

    # convert qiskit hamiltonian from openfermion
    hamil = SparsePauliOp(hamil_str, coeffs=hamil_coeffs)
    return hamil, num_qubit


def expval(op, qc0):

    estimator = Estimator()
    job = estimator.run(qc0, op)
    result = job.result()

    return result.values
