# DNA-Sequence-Assembly
Python script that generates the original DNA sequence with the use of De Bruijn graphs

## Contents

1. [Introduction](#1-introduction)
2. [De Bruijn graphs](#2-de-bruijn-graphs)
3. [Euler Paths](#3-euler-paths)
4. [Script Implementation](#4-script-implementation)
5. [Test Cases](#5-test-cases)
6. [How to Use](#6-how-to-use)
7. [References](#7-references)
8. [License](#8-license)

## 1. Introduction
A nucleotide is the basic building block of DNA sequences. The four nucleobases that can appear
in a DNA sequence are A,C,G and T. Thus, a DNA sequence is a long chain (string) of nucleotides
such as (the following sequence will be referred to it as SQ1) TTAATTACTCACTGGCTAATTACTCACTGGGTCACTACGCACTG

In order to construct a DNA sequence of a given DNA, typically the sequence is multiplied
and segmented into different lengths. These smaller segments of the DNA undergo a chemical
process in order to know the order of nucleotides in each segment. The following sequences are
segments of SQ1:

- TTAATTA
- ATTACTC
- ACTCAC
- TCACTGGCTAA
- CTAATTACTCACTGG
- CTGGGT
- GGGTCACT
- CACTACGCACTG
  
These DNA sequences of the smaller segments typically overlap. Thus, they have to be
aligned or connected together in order to know the correct sequence. For example, the correct
alignment of the previous set of sequences relative to the original DNA sequence is as follows

<p align="center">
  <img width="371" height="180" alt="image" src="https://github.com/user-attachments/assets/a442ce88-4261-4495-91d1-e21d40777809" />
</p>
Without a reference DNA sequence, there is no algorithm that will get the correct sequence
for all the possible cases of segments. However, the overall script's logic is working under the assumption that the correct DNA sequence will be given from the given segments.
Main Objective of that project, is to implement a program that given overlapping
sequences of small DNA segments, it assembles these segments and returns the whole DNA
sequence.

## 2. De Bruijn graphs

In graph theory, the standard de Bruijn graph is the graph obtained by taking all strings over any finite alphabet of length ℓ
 as vertices, and adding edges between vertices that have an overlap of ℓ−1. In the following, we consider assembly using a slightly modified version of the standard de Bruijn graph from the L-spectrum of a genome.
Given the L-spectrum of a genome, we construct a de Bruijn graph as follows:

- Add a vertex for each (L-1)-mer in the L-spectrum.

- Add k edges between two (L-1)-mers if their overlap has length L-2 and the corresponding L-mer appears k times in the L-spectrum.

An example de Bruijn graph construction is shown below.

<p align="center">
  <img width="826" height="521" alt="image" src="https://github.com/user-attachments/assets/21bf7c45-d167-4e74-b7ef-4996a9505e06" />
</p>

Source: Pevsner, J. (2015). Bioinformatics and Functional Genomics (3rd ed.). John Wiley & Sons. 

## 3. Euler Paths

An Euler Path is a route through a connected graph such that each edge of the graph is used exactly once. This graphical property is important in reconstructing sequences from a de Bruijn Graph. The goal is to find a path that traverses every "edge" (k-mer) of the graph exactly once and then reconstruct the original, long DNA sequence. 

<p align = "center">
<img width="363" height="362" alt="image" src="https://github.com/user-attachments/assets/cd4456e3-8b00-4d05-97f5-668887638839" />
</p>

Image source: https://dragoncurvetutoring.org/graphtheory.html

## 4. Script Implementation

### 2.1 Function `read_csv(name)`

For that function, the `read_csv` function from the pandas library was used. An
alternative approach would have been to read each .csv file with the built-in
csv module by creating a reader object and creating a for loop to iterate over
each line separately. After collecting all the rows, object could get converted
to a pandas DataFrame. I find my approach much better as it doesn’t need an
extended amount of lines as its alternative in order to perform the same task.

The provided .csv file must have the following format in order to be parsed properly:
<p align="center">
<img width="320" height="125" alt="image" src="https://github.com/user-attachments/assets/3c665bff-ee87-43ed-84f2-39eea1aa8299" />
</p>

### 2.2 Function `clean_data(df)`

The `clean_data()` function first makes a copy of the input DataFrame, then
drops exact duplicates on the columns SegmentNr, Position, A, C, G, T and
sorts by segment number and position. It then builds a list of “bad” row indices
by scanning each segment for three error conditions: (1) missing positions (the
number of rows is less than the maximum position), (2) conflicting entries at the
same position (duplicate positions with different nucleotide hits), and (3) rows
where the nucleotide hits don’t sum to exactly one. All rows at those bad indices
are removed in bulk. Finally, it calls the `_sequencer` helper to reconstruct each
segment’s nucleotide string, identifies any segments whose sequences are similar
to another one, and drops those segments as well, returning a fully cleaned
DataFrame.

### 2.3 Function `generate_sequences(df)`

The `generate_sequences()` function uses a pandas DataFrame that contains
the sequence information, extracts the nucleotide sequences for each segment
using the helper function `_sequencer()`, and converts the resulting dictionary
of sequences to JSON. An alternative way would have been to manually read the
dictionary generated by the helper function and construct a string that mimicks
the JSON format.

### 2.4 Function `construct_graph(json data, k)`

To generate all k-mers for each segment, I wrote a recursive helper called
`find_kmers(sequence, k)`. It first checks the base case, if the remaining
sequence is shorter than k, it returns an empty list. Otherwise, it takes the
substring from position 0 to k, appends it to the result, and then calls itself on
the sequence starting at position 1. This recursion continues until the base case
is reached. I chose recursion here because it generates the k-mers in just a few
lines. Although an iterative for/while version would work too, it tends to be
harder to read, so I chose recursion instead.

### 2.5 Function `plot_graph(graph, filename)`

That function inputs a graph object and plots a graph with the help of networkx
package. After initializing the layout (shell), the function draws the nodes,
labels, and edges. When the function is called through command line, it saves
the figure under the prespecified filename as a .PNG file.

### 2.6 Function `is_valid_graph(graph)`

That function determines whether a directed multigraph meets the requirements
for being Eulerian by first computing, for each vertex, the difference between
its out-degree and in-degree and collecting any “non-compliant” vertices (those
whose difference is nonzero) into the dictionary nc_vertices. A requirement
for a graph to be Eulerian is that nc_vertices = 0, while it permits exactly two
entries with differences of +1 and −1. If nc_vertices is neither 0 nor 2, the
function immediately returns False. When there are no non-compliant vertices,
it selects the first node as the start vertex, and when there are exactly two, it
chooses the node with a difference +1 as the start_vertex. It then builds the
list c_nodes of all vertices with nonzero total degree (e.g., those participating
in at least one edge) and performs a breadth-first search from start_vertex,
recording each reachable node in visited. Finally, it checks that every node in
c_nodes appears in visited; if any are missing, it returns False, otherwise it
returns True, thereby adhering to all given rules for Eulerian graphs.

### 2.7 Function `construct_dna_sequence(graph)`

This function assembles the DNA sequence by traversing the Eulerian path of
a de Bruijn graph, using the helper `_construct_euler_path(graph)`. That
helper identifies the starting vertex by calling the helper `_find_start(graph)`,
which verifies the graph’s Eulerian conditions and returns the appropriate start_vertex.
After identifying the starting vertex, the function employs Hierholzer’s algo-
rithm, using a stack to ensure each edge is visited exactly once, and produces
the ordered list of k-mer nodes. Because vertices are appended in reverse order,
the result is reversed to restore the correct path. The main function first calls
the `is_valid_graph` function and, if that returns True, it then obtains the Eu-
lerian path. Finally, it concatenates the first k-mer and the last character of
each subsequent k-mer.

### 2.8 Function `save_output(s, filename)`

That function accepts the full DNA sequence (s) and saves a .txt file containing
the sequence locally.

## 5. Test Cases
`project_test.py` corresponds to all the unit tests that were used throughout the project to ensure robust performance of the script

## 6. How-to-Use
The program is composed from a single .py file. One goal in the future is to design a GUI and make the experience more user friendly. 
Download the `project.py` file and before you run it, ensure that you have Python (version 3) installed to your computer.

The program is runnable with the following command from the command line:

  `python project.py DNA_[x]_[k].csv`

where `x` is the number referring to which DNA this .csv file corresponds to and `k` the length of kmers which is needed to construct the de Bruijn graph.

## 7. References
1. https://dragoncurvetutoring.org/graphtheory.html
2. Pevsner, J. (2015). Bioinformatics and Functional Genomics (3rd ed.). John Wiley & Sons. 

## 8. License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).
