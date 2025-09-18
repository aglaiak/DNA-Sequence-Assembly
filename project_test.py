from project import clean_data, generate_sequences, construct_graph, is_valid_graph, construct_dna_sequence

from pytest import mark
import pandas as pd
import networkx as nx


@mark.parametrize(
    'dna_df, expected',
    [
        (
                pd.DataFrame(data=[
                    [1, 1, 1, 0, 0, 1],
                    [1, 2, 0, 0, 0, 1],
                    [2, 1, 1, 0, 0, 0],
                    [2, 2, 0, 1, 0, 0]],
                    columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
                pd.DataFrame(data=[
                    [2, 1, 1, 0, 0, 0],
                    [2, 2, 0, 1, 0, 0]],
                    columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']
                    ),
        ),   
        ### Test from the project example
        (        pd.DataFrame(data= [
                    [1, 1, 1, 0, 0, 0],
                    [1, 2, 0, 0, 0, 1],
                    [1, 3, 0, 0, 1, 0],
                    [1, 4, 1, 0, 0, 0],
                    [1, 4, 1, 0, 0, 0],
                    [1, 5, 1, 0, 0, 0],
                    [5, 6, 1, 0, 0, 0],
                    [5, 6, 0, 0, 0, 1],
                    [2, 1, 0, 1, 0, 0],
                    [2, 2, 0, 0, 0, 1],
                    [2, 3, 0, 0, 1, 0],
                    [2, 4, 1, 0, 0, 0],
                    [2, 5, 1, 0, 0, 0],
                    [2, 6, 0, 0, 0, 1],
                    [2, 7, 0, 0, 1, 0],
                    [2, 8, 1, 0, 0, 0],
                    [3, 1, 0, 1, 0, 0],
                    [3, 2, 0, 0, 0, 1],
                    [3, 3, 0, 0, 1, 0],
                    [3, 4, 1, 0, 0, 0],
                    [3, 5, 1, 0, 0, 0],
                    [3, 6, 0, 0, 0, 1],
                    [3, 7, 0, 0, 1, 0],
                    [3, 8, 1, 0, 0, 0],
                    [4, 1, 1, 0, 0, 1],
                    [4, 2, 0, 0, 0, 1],
                    [4, 3, 0, 0, 1, 0],
                    [4, 4, 1, 0, 0, 0],
                    [4, 4, 0, 0, 1, 0],
                    [5, 1, 1, 0, 0, 0],
                    [5, 3, 0, 0, 0, 1]],
                    columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']
                    ),
                pd.DataFrame([
                    [1, 1, 1, 0, 0, 0],
                    [1, 2, 0, 0, 0, 1],
                    [1, 3, 0, 0, 1, 0],
                    [1, 4, 1, 0, 0, 0],
                    [1, 5, 1, 0, 0, 0],
                    [2, 1, 0, 1, 0, 0],
                    [2, 2, 0, 0, 0, 1],
                    [2, 3, 0, 0, 1, 0],
                    [2, 4, 1, 0, 0, 0],
                    [2, 5, 1, 0, 0, 0],
                    [2, 6, 0, 0, 0, 1],
                    [2, 7, 0, 0, 1, 0],
                    [2, 8, 1, 0, 0, 0]],
                    columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T'])    
        ),

         (
                pd.DataFrame(data=[
                    [1, 1, 1, 1, 0, 1],
                    [1, 2, 0, 0, 0, 1],
                    [2, 1, 0, 0, 0 ,1]
                    ],
                    columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
                pd.DataFrame(data=
                    [[2, 1, 0, 0, 0 ,1]],
                    columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']
                    ),
        ),

    ],

)
def test_clean_data(dna_df: pd.DataFrame, expected: pd.DataFrame) -> None:
    assert clean_data(dna_df).equals(expected)


@mark.parametrize(
    'dna_df, expected_json_str',
    [ ( 
            # test case 1 with 1 segment
        pd.DataFrame(data=[
                    [1, 1, 1, 0, 0, 1],
                    [1, 2, 0, 0, 0, 1],
                    [1, 3, 1, 0, 0, 0],
                    [1, 4, 0, 1, 0, 0],
                    ],
                    columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
                '{"1": "ATAC"}'
        ),
            # test case 2 with multiple segments
    (pd.DataFrame(data=[
                    [1, 1, 1, 0, 0, 1],
                    [1, 2, 0, 0, 0, 1],
                    [1, 3, 1, 0, 0, 0],
                    [1, 4, 0, 1, 0, 0],
                    [4, 1, 0, 1, 0, 0],
                    [4, 2, 0, 1, 0, 0],
                    [4, 3, 0, 0, 1, 0],
                    [4, 4, 1, 0, 0, 0],
                    [6, 1, 0, 0, 0, 1],
                    [6, 2, 0, 0, 1, 0],
                    [6, 3, 1, 0, 0, 0],
                    [6, 4, 0, 0, 0, 1],
                    [6, 5, 1, 0, 0, 0],
                    [6, 6, 0, 0, 1, 0],
                    ],
                    columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']
                ),
                '{"1": "ATAC", "4": "CCGA", "6": "TGATAG"}'
    ),
        (
    pd.DataFrame(columns=['SegmentNr', 'Position', 'A', 'C', 'G', 'T']),
    '{}'
),
    ],
)


def test_generate_sequences(dna_df: pd.DataFrame, expected_json_str: str) -> None:
    assert (generate_sequences(dna_df) == expected_json_str)


@mark.parametrize(
    'json_data, k,  expected_edge_list',
    [
        # test case 1, with only one sequence
         (
            '{"1":"CCTGAACC"}',
            3,
            [
                ("CC","CT"), ("CT","TG"), ("TG","GA"),
                ("GA","AA"), ("AA","AC"), ("AC","CC"),
            ],
        ),

        (   # test case 2,  with 3 segments and k = 3
            '{"1":"AACTGC","2":"AATCC","3":"CCAAT"}',
            3,
            [
                ("AA","AC"), ("AC","CT"), ("CT","TG"), ("TG","GC"),
                ("AA","AT"), ("AT","TC"), ("TC","CC"),
                ("CC","CA"), ("CA","AA"), ("AA","AT"),
            ],
        ),
        (   # test case 3 with a k value higher than the length of the sequence
            '{"2":"AACCC"}',
            7,
            [],
        ), # test case 4 with 3 segments but one of which has a length below k
          (
            '{"1":"GGGCTA","7":"GGGGGT","10":"AA"}',
            5,
            [
                ("GGGC","GGCT"),
                ("GGCT","GCTA"),
                ("GGGG","GGGG"),
                ("GGGG","GGGT"),
            ],
        ),
    ])
def test_construct_graph(json_data: str, k: int, expected_edge_list: list) -> None:
    G = construct_graph(json_data, k)
    assert sorted(G.edges()) == sorted(expected_edge_list)

@mark.parametrize(
    'DNA_edge_list,  expected_validity',
    [
        (       # test case 1 with valid vertices
                [('ATTA', 'TTAC'), ('TTAC', 'TACT'), ('TACT', 'ACTC'), ('ACTC', 'ATTA')],
                True
        ), 

        (       # test case 2 with no valid vertices
                [('GG', 'AT'), ('AT', 'CT'), ('TT', 'TA')],
                False
        ),

        (       # test case 3 with 3 invalid vertices
                [('GG', 'AT'), ('GG', 'GT'), ('GT', 'TA')],
                False
        ),

        (       # test case 4 with one edge having two outgoing edges
                [('GG', 'GA'), ('GG', 'GT'), ('GA', 'AG'), ('GT', 'TG'),
                ('AG', 'GG'), ('TG', 'GG')],
                True
        ),
    ])
def test_is_valid_graph(DNA_edge_list: list, expected_validity: bool) -> None:
    debruijn_graph = nx.MultiDiGraph()
    for edge in DNA_edge_list:
        debruijn_graph.add_edge(edge[0], edge[1])

    assert is_valid_graph(debruijn_graph) is expected_validity


@mark.parametrize(
    'DNA_edge_list,  possible_dna_sequence',
    [       # test case 1
        (
                [('AAA', 'AAC'), ('AAC', 'ACA'), ('ACA', 'CAC')],
                ["AAACAC"]
        ),
        (   # test case 2
                [('TTA', 'TTA')],
                ['TTA']
        ),
        (   # test case 3 for circular graphs
                [('AA','TT'), ('TT','GG'), ('GG', 'CA'), ('CA', 'AA')],
                ['AATGA']
        ),
    ])
def test_construct_dna_sequence(DNA_edge_list: list, possible_dna_sequence) -> None:
    debruijn_graph = nx.MultiDiGraph()
    for edge in DNA_edge_list:
        debruijn_graph.add_edge(edge[0], edge[1])

    assert construct_dna_sequence(debruijn_graph) in possible_dna_sequence