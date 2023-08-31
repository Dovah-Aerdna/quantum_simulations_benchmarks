import qiskit
from qiskit_aer import AerSimulator
from qiskit import Aer
from qiskit.circuit.random import random_circuit

# ok for cpu and gpu
# qubit = 28 
# depth = 30

# too much
# qubit = 29
# depth = 3
cpu_time=0.
gpu_time=0.

for qubit in range(22,29,2):
        for depth in range(5,101,5):
                        for _ in range(3):
                                
                                circ = random_circuit(qubit, depth, measure=True)

                                aersim = AerSimulator()

                                result_ideal = qiskit.execute(circ, aersim).result()
                                cpu_time = result_ideal.time_taken

                                #counts_ideal = result_ideal.get_counts(0)
                                #print('random circuit, qubits: ',qubit,' depth: ',depth)

                                #print('time cpu: ', result_ideal.time_taken)

                                #for bk in Aer.backends():
                                #    print(bk)


                                # simulator_gpu = Aer.get_backend('aer_simulator')
                                # simulator_gpu.set_options(device='GPU')


                                gpu_sim = Aer.get_backend("aer_simulator_statevector_gpu")

                                circ = qiskit.transpile(circ, gpu_sim)

                                result = gpu_sim.run(circ).result()
                                gpu_time = result.time_taken
                                #counts = result.get_counts(circ)

                                #print('time gpu: ',result.time_taken)
                                print(qubit,'\t',depth,'\t',cpu_time,'\t',gpu_time)
                        