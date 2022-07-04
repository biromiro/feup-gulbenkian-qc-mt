import math
import numpy as np
from qiskit import Aer
from qiskit.utils import QuantumInstance
from qiskit.algorithms import Shor
from qiskit.providers.aer import AerError
from qiskit import IBMQ

# print available qiskit simulators
IBMQ.load_account()  # Load account from disk
provider = IBMQ.get_provider(hub='ibm-q')
# print(provider.backends())
backend = provider.get_backend('ibmq_qasm_simulator')


def add_prime_mr(q1, p, a):
    print(f"q1 = {q1}, p = {p}")
    q2 = q1 * p
    quantum_instance1 = QuantumInstance(backend, shots=1024)
    quantum_instance2 = QuantumInstance(backend, shots=1024)
    shor1 = Shor(quantum_instance=quantum_instance1)
    shor2 = Shor(quantum_instance=quantum_instance2)

    result1 = shor1.factor(q1)
    # print(f"{result1}")
    result2 = shor2.factor(q2)
    # print(f"{result2}")

    resultList1 = [] if result1.factors == [] else result1.factors[0]
    resultList2 = [] if result2.factors == [] else result2.factors[0]

    difference = list(set(resultList2) - set(resultList1))

    return (True, '') if (difference == [] or difference == [p]) else (False, (result1, result2))


primes_list = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
               41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

for prime in primes_list:
    for q in range(15, 50, 2):
        if (q in primes_list or q < prime):
            continue
        result = ''
        for a in [2] + primes_list:
            if (a > q or result != ''):
                break
            val, result = add_prime_mr(q, prime, a)
        print(result)

        # if not val:
        #    print(
        #        f"The metamorphic test did not succeed. The primes of {q} were {result[0]}, and the primes of {q}*{prime} were {result[1]}.")
        # else:
        #    print(f"The metamorphic test succeeded.")
