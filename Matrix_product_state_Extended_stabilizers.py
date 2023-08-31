"""Utility functions for generating random Clifford circuits."""

import numpy as np
import qiskit
from qiskit.circuit import ClassicalRegister, QuantumCircuit, CircuitInstruction
from qiskit.circuit.library import standard_gates
from random import randrange


def random_clifford_circuit(
    num_qubits, depth, T_gate=0, measure=False, seed=None, draw=False
):
    #qiskit.quantum_info.Clifford

    if num_qubits == 0:
        return QuantumCircuit()

    gates_1q = [
        # (Gate class, number of qubits)
        #(standard_gates.IGate, 1),
        (standard_gates.XGate, 1),
        (standard_gates.HGate, 1),
        (standard_gates.SGate, 1),
        (standard_gates.SdgGate, 1),
        (standard_gates.YGate, 1),
        (standard_gates.ZGate, 1),
    ]
    gates_2q = [
        (standard_gates.CXGate, 2),
        (standard_gates.CZGate, 2),
        (standard_gates.SwapGate, 2),

    ]

    gates = gates_1q.copy()
    gates.extend(gates_2q)

    gates = np.array(
        gates, dtype=[("class", object), ("num_qubits", np.int64)]
    )
    gates_1q = np.array(gates_1q, dtype=gates.dtype)

    qc = QuantumCircuit(num_qubits)

    if measure:
        cr = ClassicalRegister(num_qubits, "c")
        qc.add_register(cr)

    if seed is None:
        seed = np.random.randint(0, np.iinfo(np.int32).max)
    rng = np.random.default_rng(seed)

    qubits = np.array(qc.qubits, dtype=object, copy=True)

    # Apply arbitrary random operations in layers across all qubits.
    for _ in range(depth):
        # We generate all the randomness for the layer in one go, to avoid many separate calls to
        # the randomisation routines, which can be fairly slow.

        # This reliably draws too much randomness, but it's less expensive than looping over more
        # calls to the rng. After, trim it down by finding the point when we've used all the qubits.
        gate_specs = rng.choice(gates, size=len(qubits))
        cumulative_qubits = np.cumsum(gate_specs["num_qubits"], dtype=np.int64)
        # Efficiently find the point in the list where the total gates would use as many as
        # possible of, but not more than, the number of qubits in the layer.  If there's slack, fill
        # it with 1q gates.
        max_index = np.searchsorted(cumulative_qubits, num_qubits, side="right")
        gate_specs = gate_specs[:max_index]
        slack = num_qubits - cumulative_qubits[max_index - 1]
        if slack:
            gate_specs = np.hstack((gate_specs, rng.choice(gates_1q, size=slack)))


        # For efficiency in the Python loop, this uses Numpy vectorisation to pre-calculate the
        # indices into the lists of qubits and parameters for every gate, and then suitably
        # randomises those lists.
        q_indices = np.empty(len(gate_specs) + 1, dtype=np.int64)
        q_indices[0] = 0
        np.cumsum(gate_specs["num_qubits"], out=q_indices[1:])
        rng.shuffle(qubits)

        
        for gate, q_start, q_end in zip(
            gate_specs["class"], q_indices[:-1], q_indices[1:]
        ):
            operation = gate()
            qc._append(CircuitInstruction(operation=operation, qubits=qubits[q_start:q_end]))

    # insert here t gates
    qc_data = []
    for instruction, qargs, cargs in qc:
        qc_data.append((instruction, qargs, cargs))

    while T_gate!=0:
        # substitute gates with t gate
        idx=randrange(len(qc_data))
        #print(qc_data[idx])
        #print(qc_data[idx][0].num_qubits)

        if qc_data[idx][0].num_qubits==1:
            if qc_data[idx][0].name=='t':
                continue
            qc_data[idx]=(qiskit.circuit.library.TGate(),qc_data[idx][1],[])

        if qc_data[idx][0].num_qubits==2:
            #print(type(qc_data[idx][1]))
            #print(qc_data[idx][1][0])
            args=qc_data[idx][1]
            qc_data[idx]=(qiskit.circuit.library.TGate(),[args[0]],[])
            qc_data.append((rng.choice(gates_1q)[0](),[args[1]],[]))

        T_gate=T_gate-1

    qc.data = qc_data


    if measure:
        qc.measure(qc.qubits, cr)
    if draw:
        print(qc.draw())

    return qc



#######################################

ch = input("Select simulation method:\n1. extended_stabilizer\n2. matrix_product_state\n3. exit\n")

if ch=="1":
    method="extended_stabilizer" #(23 qb: t=12 -> 97 s; t=0 -> 6 s)
elif ch=="2":
    method="matrix_product_state" # (23 qb: t=12 -> 146 s; t=0 -> 600 s (number double gateis higher))
else:
    exit()


simulator = qiskit.Aer.get_backend("aer_simulator_"+method) 


with open(method+".dat","w") as file:

    line = "qubits\tdepth\tT_gates\ttime\n"
    file.write(line)

    for qb in range(16,21,4):
        for dp in range(5,51,15):
            for tg in range(0,11,2):
                for _ in range(5):

                    print("progress: ",qb,"/30 qubit\t",dp,"/50 depth\t",tg,"/20 T Gate")
                    file.write(str(qb)+"\t"+str(dp)+"\t"+str(tg)+"\t")

                    circ = random_clifford_circuit(qb,dp, T_gate=tg, measure=True)
                    transpiled = qiskit.transpile(circ,simulator)
                    result = simulator.run(transpiled).result()

                    file.write(str(result.metadata['time_taken_execute'])+"\n")


