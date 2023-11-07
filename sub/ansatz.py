'''
  ansatz
'''


# particle number conserving ansatz
def pqc_ansatz_PNC(qc, params, reps_pqc):
    num_qubits = qc.num_qubits

    for i in range(0, reps_pqc):

        i_num = num_qubits // 2
        for j in range(i_num):
            qc.cx(j * 2, j * 2 + 1)
            qc.ry(-params[(qc.num_qubits - 1) * i + j], j * 2)
            qc.cx(j * 2 + 1, j * 2)
            qc.ry(params[(qc.num_qubits - 1) * i + j], j * 2)
            qc.cx(j * 2, j * 2 + 1)

        for k in range(i_num - 1):
            qc.cx(k * 2 + 1, k * 2 + 2)
            qc.ry(-params[(qc.num_qubits - 1) * i + k + i_num], k * 2 + 1)
            qc.cx(k * 2 + 2, k * 2 + 1)
            qc.ry(params[(qc.num_qubits - 1) * i + k + i_num], k * 2 + 1)
            qc.cx(k * 2 + 1, k * 2 + 2)

    return qc


# hardware-efficient ansatz
def pqc_ansatz_HE(qc, params, reps_pqc):
    num_qubits = qc.num_qubits

    for i in range(0, reps_pqc):
        for j in range(num_qubits - 1):
            if j % 2 == 0:
                qc.cx(j, j + 1)
            else:
                pass
        for j in range(num_qubits - 1):
            if j % 2 == 1:
                qc.cx(j, j + 1)
            else:
                pass
        for j in range(num_qubits):
            qc.ry(params[qc.num_qubits * i + j], j)

        qc.barrier()

    return qc
