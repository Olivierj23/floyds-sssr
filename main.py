import itertools
import time

def path_combination(pathlist_1, pathlist_2):
    return_list = []
    for path1 in pathlist_1:
        for path2 in pathlist_2:
            combined_path = path1.union(path2)
            return_list.append(combined_path)

    return return_list

def make_pid_matrices(adjacency_dict):
    # intialization
    distance_matrix = {}
    pid_e = {}
    pid_e_prime = {}
    seen = set()
    index_to_edge_map = {}
    for start_node in adjacency_dict:
        distance_matrix[start_node] = {}
        pid_e[start_node] = {}
        pid_e_prime[start_node] = {}
        for end_node in adjacency_dict[start_node]:
            if adjacency_dict[start_node][end_node] == 1:
                distance_matrix[start_node][end_node] = 1
                pid_e[start_node][end_node] = [set([(start_node, end_node)])]
            else:
                distance_matrix[start_node][end_node] = float("inf")
                pid_e[start_node][end_node] = []

            if (start_node, end_node) not in seen:
                index_to_edge_map[len(index_to_edge_map)] = [(start_node, end_node), (end_node, start_node)]
                seen.add((start_node, end_node))
                seen.add((end_node, start_node))


            pid_e_prime[start_node][end_node] = []

    # floyd warshall

    n = len(adjacency_dict)

    for k in adjacency_dict:
        for i in adjacency_dict:
            for j in adjacency_dict:
                if distance_matrix[i][j] > distance_matrix[i][k] + distance_matrix[k][j]:
                    if distance_matrix[i][j] == distance_matrix[i][k] + distance_matrix[k][j] + 1:
                        pid_e_prime[i][j] = pid_e[i][j]

                    else:
                        pid_e_prime[i][j] = []

                    distance_matrix[i][j] = distance_matrix[i][k] + distance_matrix[k][j]
                    pid_e[i][j] = path_combination(pid_e[i][k], pid_e[k][j])

                elif distance_matrix[i][j] == distance_matrix[i][k] + distance_matrix[k][j]:
                    pid_e[i][j].extend(path_combination(pid_e[i][k], pid_e[k][j]))

                elif distance_matrix[i][j] == distance_matrix[i][k] + distance_matrix[k][j] - 1:
                    pid_e_prime[i][j].extend(path_combination(pid_e[i][k], pid_e[k][j]))
                else:
                    pass

    for i in distance_matrix:
        for j in distance_matrix:
            if i == j:
                distance_matrix[i][j] = 0
                pid_e[i][j] = []
                pid_e_prime[i][j] = []

    return distance_matrix, pid_e, pid_e_prime, index_to_edge_map

def find_theoretical_sssr(adjacency_dict):
    m = 0
    n = len(adjacency_dict)
    for i in adjacency_dict:
        for j in adjacency_dict:
            if adjacency_dict[i][j] == 1:
                m += 0.5

    return m - n + 1



def make_cset(distance_matrix, pid_e, pid_e_prime):
    n = len(distance_matrix)

    c_set = []

    for i in distance_matrix:
        for j in distance_matrix:
            if distance_matrix[i][j] == 0 or (len(pid_e[i][j]) == 1 and len(pid_e_prime[i][j]) == 0):
                continue
            else:
                if len(pid_e_prime[i][j]) != 0:
                    c_num = 2 * (distance_matrix[i][j] + 0.5)
                    c_set.append([c_num, pid_e[i][j], pid_e_prime[i][j]])
                if len(pid_e[i][j]) > 1:
                    c_num = 2 * (distance_matrix[i][j])
                    c_set.append([c_num, pid_e[i][j], pid_e_prime[i][j]])

    c_set = sorted(c_set)

    return c_set

def is_cycle(path1, path2):
    for edge in path1:
        if edge in path2 or edge[::-1] in path2:
            return False

    return True

def xor(cycle_1, cycle_2):
    return_set = set()

    for edge in cycle_1:
        if not edge in cycle_2 and not edge[::-1] in cycle_2:
            return_set.add(edge)

    for edge in cycle_2:
        if not edge in cycle_1 and not edge[::-1] in cycle_1:
            return_set.add(edge)

    return return_set

def is_good_cycle(cycle, c_sssr):
    cycle_1 = cycle
    for length in range(len(c_sssr) + 1):
        for subset in itertools.combinations(c_sssr, length):
            target_xor = cycle_1
            for cycle_3 in subset:
                target_xor = xor(target_xor, cycle_3)

            if len(target_xor) == 0:
                return False

    return True


def find_sssr(c_set, n_sssr):
    n_ringidx = 0
    c_sssr = []

    for i in range(len(c_set)):
        if int(c_set[i][0]) % 2 == 1:
            for j in range(len(c_set[i][2])):
                if is_cycle(c_set[i][1][0], c_set[i][2][j]):
                    cycle = c_set[i][1][0].union(c_set[i][2][j])
                    if is_good_cycle(cycle, c_sssr):
                        c_sssr.append(cycle)
                        n_ringidx += 1
                        print(len(c_sssr))
                    if n_ringidx == n_sssr:
                        return c_sssr

        elif int(c_set[i][0]) % 2 == 0:
            for j in range(len(c_set[i][1]) - 1):
                if is_cycle(c_set[i][1][j], c_set[i][1][j + 1]):
                    cycle = c_set[i][1][j].union(c_set[i][1][j + 1])
                    if is_good_cycle(cycle, c_sssr):
                        c_sssr.append(cycle)
                        n_ringidx += 1
                        print(len(c_sssr))
                    if n_ringidx == n_sssr:
                        return c_sssr

def floyds_sssr(adjacency_dict):
    start = time.time()
    distance_matrix, pid_e, pid_e_prime, _ = make_pid_matrices(adjacency_dict)

    print(time.time()-start)
    start = time.time()

    c_set = make_cset(distance_matrix, pid_e, pid_e_prime)

    print(time.time() - start)
    start = time.time()

    n_sssr = find_theoretical_sssr(adjacency_dict)
    result = find_sssr(c_set, n_sssr)

    print(time.time() - start)

    return result

def xor_bit(cycle_1, cycle_2):
    return bin(int(cycle_1, 2) ^ int(cycle_2, 2))[2:]

def is_good_cycle_bit(cycle_bit, c_sssr_bit):
    cycle_1 = cycle_bit
    for length in range(len(c_sssr_bit) + 1):
        for subset in itertools.combinations(c_sssr_bit, length):
            target_xor = cycle_1
            for cycle_3 in subset:
                target_xor = xor_bit(target_xor, cycle_3)

            if target_xor == "0":
                return False

    return True

def converter(cycle, index_to_edge_map):
    return_list = []
    for i in range(len(index_to_edge_map)):
        edge, reverse_edge = index_to_edge_map[i]
        if edge in cycle or reverse_edge in cycle:
            return_list.append("1")
        else:
            return_list.append("0")

    return "".join(return_list)


def find_sssr_bit(c_set, n_sssr, index_to_edge_map):
    n_ringidx = 0
    c_sssr = []
    c_sssr_bit = []

    for i in range(len(c_set)):
        if int(c_set[i][0]) % 2 == 1:
            for j in range(len(c_set[i][2])):
                if is_cycle(c_set[i][1][0], c_set[i][2][j]):
                    cycle = c_set[i][1][0].union(c_set[i][2][j])
                    cycle_bit = converter(cycle, index_to_edge_map)
                    if is_good_cycle_bit(cycle_bit, c_sssr_bit):
                        c_sssr_bit.append(cycle_bit)
                        c_sssr.append(cycle)
                        n_ringidx += 1
                        print(len(c_sssr))
                    if n_ringidx == n_sssr:
                        return c_sssr

        elif int(c_set[i][0]) % 2 == 0:
            for j in range(len(c_set[i][1]) - 1):
                if is_cycle(c_set[i][1][j], c_set[i][1][j + 1]):
                    cycle = c_set[i][1][j].union(c_set[i][1][j + 1])
                    cycle_bit = converter(cycle, index_to_edge_map)
                    if is_good_cycle_bit(cycle_bit, c_sssr_bit):
                        c_sssr_bit.append(cycle_bit)
                        c_sssr.append(cycle)
                        n_ringidx += 1
                        print(len(c_sssr))
                    if n_ringidx == n_sssr:
                        return c_sssr

def floyds_sssr_bit(adjacency_dict):
    start = time.time()
    distance_matrix, pid_e, pid_e_prime, index_to_edge_map = make_pid_matrices(adjacency_dict)

    print(time.time()-start)
    start = time.time()

    c_set = make_cset(distance_matrix, pid_e, pid_e_prime)

    print(time.time() - start)
    start = time.time()

    n_sssr = find_theoretical_sssr(adjacency_dict)
    result = find_sssr_bit(c_set, n_sssr, index_to_edge_map)

    print(time.time() - start)

    return result

adjacency_dict_1 = {0: {0: 0, 1: 1, 2: 0, 3: 1},
                  1: {0: 1, 1: 0, 2: 1, 3: 0},
                  2: {0: 1, 1: 0, 2: 0, 3: 1},
                  3: {0: 1, 1: 0, 2: 0, 3: 0}
                }

adjacency_dict_2 = {0: {0: 0, 1: 1, 2: 0, 3: 0, 4: 0},
                    1: {0: 0, 1: 0, 2: 1, 3: 0, 4: 0},
                    2: {0: 0, 1: 1, 2: 0, 3: 0, 4: 0},
                    3: {0: 0, 1: 0, 2: 1, 3: 0, 4: 1},
                    4: {0: 1, 1: 1, 2: 0, 3: 0, 4: 0},
                }

adjacency_dict_3 = {"a": {"a": 0, "b": 1, "c": 0, "d": 0, "e": 0, "f": 0, "g": 1, "h": 1},
                    "b": {"a": 1, "b": 0, "c": 1, "d": 0, "e": 0, "f": 1, "g": 1, "h": 1},
                    "c": {"a": 0, "b": 1, "c": 0, "d": 0, "e": 1, "f": 0, "g": 0, "h": 1},
                    "d": {"a": 0, "b": 0, "c": 0, "d": 0, "e": 1, "f": 0, "g": 0, "h": 0},
                    "e": {"a": 0, "b": 0, "c": 1, "d": 1, "e": 0, "f": 1, "g": 0, "h": 0},
                    "f": {"a": 0, "b": 1, "c": 0, "d": 0, "e": 1, "f": 0, "g": 1, "h": 0},
                    "g": {"a": 1, "b": 1, "c": 0, "d": 0, "e": 0, "f": 1, "g": 0, "h": 0},
                    "h": {"a": 1, "b": 1, "c": 1, "d": 0, "e": 0, "f": 0, "g": 0, "h": 0}
                }

adjacency_dict_4 = {1: {1: 0, 2: 1, 3: 0, 4: 0, 5: 0, 6: 0, 7: 1 },
                    2: {1: 1, 2: 0, 3: 1, 4: 0, 5: 0, 6: 1, 7: 0 },
                    3: {1: 0, 2: 1, 3: 0, 4: 1, 5: 0, 6: 0, 7: 0 },
                    4: {1: 0, 2: 0, 3: 1, 4: 0, 5: 1, 6: 0, 7: 0 },
                    5: {1: 0, 2: 0, 3: 0, 4: 1, 5: 0, 6: 1, 7: 0 },
                    6: {1: 0, 2: 1, 3: 0, 4: 0, 5: 1, 6: 0, 7: 1 },
                    7: {1: 1, 2: 0, 3: 0, 4: 0, 5: 0, 6: 1, 7: 0 },
                }

adjacency_dict_B = {1: {1: 0, 2: 1, 3: 1, 4: 0, 5: 1, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0},
                    2: {1: 1, 2: 0, 3: 0, 4: 1, 5: 0, 6: 0, 7: 0, 8: 0, 9: 1, 10: 0, 11: 0, 12: 0, 13: 0},
                    3: {1: 1, 2: 0, 3: 0, 4: 1, 5: 0, 6: 1, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0},
                    4: {1: 0, 2: 1, 3: 1, 4: 0, 5: 0, 6: 0, 7: 1, 8: 1, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0},
                    5: {1: 1, 2: 0, 3: 0, 4: 0, 5: 0, 6: 1, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 1, 13: 0},
                    6: {1: 0, 2: 0, 3: 1, 4: 0, 5: 1, 6: 0, 7: 0, 8: 0, 9: 0, 10: 1, 11: 0, 12: 0, 13: 0},
                    7: {1: 0, 2: 0, 3: 0, 4: 1, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 1, 11: 0, 12: 0, 13: 0},
                    8: {1: 0, 2: 0, 3: 0, 4: 1, 5: 0, 6: 0, 7: 0, 8: 0, 9: 1, 10: 0, 11: 1, 12: 0, 13: 0},
                    9: {1: 0, 2: 1, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 1, 9: 0, 10: 0, 11: 0, 12: 0, 13: 1},
                    10: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 1, 7: 1, 8: 0, 9: 0, 10: 0, 11: 1, 12: 1, 13: 0},
                    11: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 1, 9: 0, 10: 1, 11: 0, 12: 0, 13: 1},
                    12: {1: 0, 2: 0, 3: 0, 4: 0, 5: 1, 6: 0, 7: 0, 8: 0, 9: 0, 10: 1, 11: 0, 12: 0, 13: 1},
                    13: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 1, 10: 0, 11: 1, 12: 1, 13: 0}
                }

adjacency_dict_C = {1: {1: 0, 2: 1, 3: 1, 4: 0, 5: 1, 6: 1, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0},
                    2: {1: 1, 2: 0, 3: 1, 4: 1, 5: 0, 6: 0, 7: 0, 8: 0, 9: 1, 10: 0, 11: 0, 12: 0, 13: 0},
                    3: {1: 1, 2: 1, 3: 0, 4: 1, 5: 0, 6: 1, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0},
                    4: {1: 0, 2: 1, 3: 1, 4: 0, 5: 0, 6: 0, 7: 1, 8: 1, 9: 1, 10: 0, 11: 0, 12: 0, 13: 0},
                    5: {1: 1, 2: 0, 3: 0, 4: 0, 5: 0, 6: 1, 7: 0, 8: 0, 9: 0, 10: 1, 11: 0, 12: 1, 13: 0},
                    6: {1: 1, 2: 0, 3: 1, 4: 0, 5: 1, 6: 0, 7: 0, 8: 0, 9: 0, 10: 1, 11: 0, 12: 0, 13: 0},
                    7: {1: 0, 2: 0, 3: 0, 4: 1, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 1, 11: 0, 12: 0, 13: 0},
                    8: {1: 0, 2: 0, 3: 0, 4: 1, 5: 0, 6: 0, 7: 0, 8: 0, 9: 1, 10: 0, 11: 1, 12: 0, 13: 1},
                    9: {1: 0, 2: 1, 3: 0, 4: 1, 5: 0, 6: 0, 7: 0, 8: 1, 9: 0, 10: 0, 11: 0, 12: 0, 13: 1},
                    10: {1: 0, 2: 0, 3: 0, 4: 0, 5: 1, 6: 1, 7: 1, 8: 0, 9: 0, 10: 0, 11: 1, 12: 1, 13: 0},
                    11: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 1, 9: 0, 10: 1, 11: 0, 12: 1, 13: 1},
                    12: {1: 0, 2: 0, 3: 0, 4: 0, 5: 1, 6: 0, 7: 0, 8: 0, 9: 0, 10: 1, 11: 1, 12: 0, 13: 1},
                    13: {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 1, 9: 1, 10: 0, 11: 1, 12: 1, 13: 0}
                }



# print(floyds_sssr(adjacency_dict_C))
# for elem in floyds_sssr_bit(adjacency_dict_C):
#     print(len(elem))
#     print(elem)

def generate_adjacency_dict(adjacency_list):
    n = len(adjacency_list)
    adjacency_dict = {}
    for i in range(n):
        adjacency_dict[i] = {}
        for j in range(n):
            adjacency_dict[i][j] = adjacency_list[i][j]

    return adjacency_dict

from rdkit.Chem import MolFromSmiles, rdmolops
from rdkit.Chem import AllChem, AddHs

# SMILES for methanol
smi='CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C'

# SMILES to mol
mol=MolFromSmiles(smi)

# Add hydrogens

# Print atom ordering
print(f'Ordering={[atom.GetSymbol() for atom in mol.GetAtoms()]}')

# Print the adjacency matrix
A=rdmolops.GetAdjacencyMatrix(mol)
print(A)

start = time.time()
print(sorted([len(r) for r in rdmolops.GetSSSR(mol)]))
print(time.time() - start)

print([len(r) for r in floyds_sssr(generate_adjacency_dict(A))])












