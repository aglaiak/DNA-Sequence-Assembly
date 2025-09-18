
"""
author: Aglaia Kakoulidou
"""

import pandas as pd
import json
import networkx as nx
import matplotlib.pyplot as plt
import os


def read_csv(name: str) -> pd.DataFrame:
    """ Function that reads a file with a given name"""
    df = pd.read_csv(name, names=["SegmentNr", "Position", "A", "C", "G", "T"], header=None)
    return df


def _sequencer(df: pd.DataFrame) -> dict:

    """ Helper function that reads a pandas dataframe with format
     [SegmentNr, Position, A, C, G, T] and extracts all the sequences in 
      a dictionary with SegmentNr as keys and sequences as the values """

    sequences = {}

    for seg in df["SegmentNr"].unique():
        seg_df = df.loc[df["SegmentNr"] == seg]

        # temporary list that stores the sequences of each segment
        seq_temp = []

        for A, C, G, T in zip(seg_df["A"], seg_df["C"], seg_df["G"], seg_df["T"]):
            if A == 1:
                seq_temp.append("A")
            elif C == 1:
                seq_temp.append("C")
            elif G == 1:
                seq_temp.append("G")
            elif T == 1:
                seq_temp.append("T")            

        sequences[int(seg)] = "".join(seq_temp)
    return sequences


def clean_data(df: pd.DataFrame) -> pd.DataFrame:
    """ Function that cleans the data of a given dataframe """

    clean_df = df.copy()
    bad_idx = []    

    # remove duplicates where position and nucleotides are the same
    clean_df = df.drop_duplicates(subset=["SegmentNr", "Position", "A", "C", "G", "T"], keep="first")

    # sort the values based on segmentnr and position
    clean_df = (clean_df.sort_values(["SegmentNr", "Position"])
                .reset_index(drop=True))

    for seg in clean_df["SegmentNr"].unique():
        seg_df = clean_df.loc[clean_df["SegmentNr"] == seg]     

        m = seg_df["Position"].max()  # calculate the value m
        len_position = len(seg_df["Position"])  # calculate the length of each segment

        # dealing with missing positions in a segment
        if len_position != m:
            bad_idx.extend(seg_df.index.tolist())
            continue

        # check for position duplicates with different nucleotide entries
        duplicates = seg_df.duplicated(subset=["Position"], keep=False)
        if duplicates.any():
            bad_idx.extend(seg_df.index.tolist())
            continue

        # check for wrong position
        # if the sum of the columns A, C, G, T > 1 it means that the segment will be ignored
        if (seg_df[["A", "C", "G", "T"]].sum(axis=1) != 1).any():
            bad_idx.extend(seg_df.index.tolist())  
            continue

    clean_df = clean_df.drop(index=bad_idx).reset_index(drop=True)

    seen = set()

    # iterate over the keys and values of the sequences dict and
    # remove the clean_df segments with the same sequences

    seen = set()
    bad_seq = []

    for key, value in _sequencer(clean_df).items():
        if value in seen:
            bad_seq.extend(
                clean_df.loc[clean_df["SegmentNr"] == key].index.tolist()
            )
        else:
            seen.add(value)

    clean_df = clean_df.drop(index=bad_seq).reset_index(drop=True)
    clean_df = clean_df.sort_values(["SegmentNr", "Position"])

    return clean_df


def generate_sequences(df: pd.DataFrame) -> json:
    """ Function that generates all the sequences of a given dataframe
    and converts them to JSON format"""

    json_seq = json.dumps(_sequencer(df))
    return json_seq


def _find_kmers(sequence: str, k: int) -> list:
    """ Helper function that finds all the kmers of a given sequence"""
    result = []
    # base case
    if len(sequence) < k:
        return result

    # find the kmer of each given sequence
    kmer1 = sequence[0:k]
    result.append(kmer1)
    result.extend(_find_kmers(sequence[1:], k))

    return result


def construct_graph(json_data: json, k: int) -> nx.MultiGraph:
    """ Function that creates a Bruijn graph based a json file with 
    sequences and a k integer"""

    # extract the sequences of the json_data file and save them in a list
    # named raw_sequences

    seq_dict = json.loads(json_data)
    raw_sequences = [sequence for sequence in seq_dict.values()]

    # create a new list (of lists) that only has the kmers
    kmers = []
    for seq in raw_sequences:
        kmers_seq = _find_kmers(seq, k)
        kmers.append(kmers_seq)

    # separate the contents of the kmers list into (k-1)mers and save them into 
    # a dictionary with keys L and R

    edges = {"L": [], "R": []}
    for temp_kmers in kmers:
        for kmer in temp_kmers:
            left = kmer[:-1]
            right = kmer[1:]
            edges["L"].append(left)
            edges["R"].append(right)

    # initialize the graph object, iterate over the dictionary
    # and iteratively add the edges

    graph = nx.MultiDiGraph()
    for L, R in zip(edges["L"], edges["R"]):
        graph.add_edge(L, R)
    return graph


def plot_graph(graph: nx.MultiDiGraph, filename: str):
    """ Function that creates a file with the given graph object"""

    pos = nx.shell_layout(graph)
    plt.figure(figsize=(8, 6), constrained_layout=True)

    # adjusting the size of the nodes
    nx.draw_networkx_nodes(
        graph,
        pos,
        node_size=200,
        linewidths=0.5,
    )

    # adjusting the size of letters
    nx.draw_networkx_labels(
        graph,
        pos,
        font_size=4
    )

    nx.draw_networkx_edges(
        graph,
        pos,
        arrowstyle='-|>',
        arrowsize=2
    )

    plt.axis("off")
    plt.savefig(filename, dpi=300)
    plt.close()


def is_valid_graph(graph: nx.MultiDiGraph) -> bool:
    """ Function that checks if a graph object is Eulerian"""

    # initialize a list that collects all the non compliant vertices
    nc_vertices = {}

    for v in graph.nodes():
        degree_diff = graph.out_degree(v) - graph.in_degree(v)

        # complying to the in_degree and out_degree rule
        if degree_diff == 0:
            continue
        else:
            nc_vertices[v] = degree_diff

    # we check if the dictionary with the non-compliant vertices
    # has a length of 0 or 2
    length = len(nc_vertices)

    if length != 0 and length != 2:
        return False

    # if the length is 0 we start from any vertex
    if length == 0:
        for v in graph.nodes():
            start_vertex = v
            break

    # if the length is 2 we need to start from the one that has
    # the outgoing edge

    elif length == 2:
        if 1 not in nc_vertices.values() or -1 not in nc_vertices.values():
            return False
        else:
            for n, d in nc_vertices.items():
                # we ensure that the starting node is with the outgoing edge
                if d == 1:
                    start_vertex = n

    # check for connectivity
    c_nodes = []
    for v in graph.nodes():
        if graph.degree(v) > 0:
            c_nodes.append(v)

    # breadth first search algorithm to traverse all the edges
    if len(c_nodes) != 0:
        visited = []
        queue = []
        queue.append(start_vertex)
        while len(queue) != 0:
            v = queue.pop(0)
            if v not in visited:
                visited.append(v)
                neighbours = list(graph[v])
                for w in neighbours:
                    queue.append(w)

    # compare c_nodes with the ones visited with bfs algorithm
    final_check = []
    for v in c_nodes:
        if v not in visited:
            final_check.append(v)
    if len(final_check) > 0:
        return False

    return True


def _find_start(graph: nx.MultiDiGraph) -> int:

    """ Helper function that finds the starting vertex of a given
    graph by taking into account the outgoing and ingoing connectivity
    """
    nc_vertices = {}

    for v in graph.nodes():
        degree_diff = graph.out_degree(v) - graph.in_degree(v)

        # complying to the in_degree and out_degree rule
        if degree_diff == 0:
            continue
        else:
            nc_vertices[v] = degree_diff

    # we check if the dictionary with the non-compliant vertices
    # has a length of 0 or 2

    length = len(nc_vertices)

    if length != 0 and length != 2:
        return False

    # if the length is 0 we start from any vertex
    if length == 0:
        for v in graph.nodes():
            start_vertex = v
            break
    # if the length is 2 we need to start from the one that has the outgoing edge             
    elif length == 2:
        if 1 not in nc_vertices.values() or -1 not in nc_vertices.values():
            return False
        else:
            for n, d in nc_vertices.items():
                # we ensure that the starting node is with the outgoing edge
                if d == 1:
                    start_vertex = n

    return start_vertex


def _construct_euler_path(graph: nx.MultiDiGraph) -> list:

    """ Helper function that builds the euler path of a given de Bruijn 
    graph """

    nc_vertices = {}

    for v in graph.nodes():
        degree_diff = graph.out_degree(v) - graph.in_degree(v)

        # complying to the in_degree and out_degree rule
        if degree_diff == 0:
            continue
        else:
            nc_vertices[v] = degree_diff

    start_vertex = _find_start(graph)

    nbrs = {}
    for v in graph.nodes():
        nbrs[v] = list(graph[v])

    stack = [start_vertex]
    euler_path = []

    # traverse the graph and extract the euler path
    while len(stack) > 0:
        v = stack[-1]
        if len(nbrs[v]) > 0:
            w = nbrs[v].pop(0)
            stack.append(w)
        else:
            euler_path.append(stack.pop())

    # in case we are dealing with a closed graph

    euler_path.reverse()

    # dealing with the edge case where the graph is circular and ends at
    # the same vertex
    if len(nc_vertices) == 0:
        euler_path.pop()

    return euler_path


def construct_dna_sequence(graph: nx.MultiDiGraph) -> str:
    """ Function that constructs the sequence given an Eulerian graph """

    # string where we will store our output
    sequence = ''

    if is_valid_graph(graph) is True:

        euler_path = _construct_euler_path(graph)

        # construct the sequence based on the euler path
        if len(euler_path) > 0:
            sequence = euler_path[0]
            for kmer in euler_path[1:]:
                sequence += kmer[-1]
    return sequence


def save_output(s: str, filename: str) -> str:
    """ Function that saves the constructed sequence to .txt"""

    if s == '':
        print("The sequence cannot be constructed")
    else:
        with open(filename, 'w') as f:
            f.write(s)

# run the program through the command line


if __name__ == "__main__":

    argv = os.sys.argv
    input_file = argv[1]

    result = input_file.rstrip(".csv")
    result = result.split('_')

    # extract x and k
    x = int(result[1])
    k = int(result[2])

    # build the pipeline for running in command line
    df = read_csv(input_file)
    df_cleaned = clean_data(df)
    json_sequences = generate_sequences(df_cleaned)
    graph_object = construct_graph(json_sequences, k)
    filename = f"DNA_{x}.png"
    graph_image = plot_graph(graph_object, filename)
    sequence = construct_dna_sequence(graph_object)
    euler_path = _construct_euler_path(graph_object)
    if is_valid_graph(graph_object) is True:
        for i, kmer in enumerate(euler_path):
            if i < len(euler_path) - 1:
                print(f"{kmer} -", end=' ')
            else:
                print(kmer)
        save_output(sequence, f"DNA_{x}.txt")
