# DNA-Sequence-Assembly
Program that generates the original DNA sequence with De Bruijn graphs

<img width="600" height="378" alt="image" src="https://github.com/user-attachments/assets/fbe9682b-3956-4a64-9e73-69ff4cfe6375" />


# Objective & General Information
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

<img width="371" height="180" alt="image" src="https://github.com/user-attachments/assets/a442ce88-4261-4495-91d1-e21d40777809" />

Without a reference DNA sequence, there is no algorithm that will get the correct sequence
for all the possible cases of segments. However, we are going to add some assumptions under
which we guarantee to find the correct DNA sequence from the given segments.
Main Objective of that project, is toimplement a program that given overlapping
sequences of small DNA segments it assembles these segments and returns the whole DNA
sequence.

# De Bruijn graphs

In graph theory, the standard de Bruijn graph is the graph obtained by taking all strings over any finite alphabet of length ℓ
 as vertices, and adding edges between vertices that have an overlap of ℓ−1. In the following, we consider assembly using a slightly modified version of the standard de Bruijn graph from the L-spectrum of a genome.
Given the L-spectrum of a genome, we construct a de Bruijn graph as follows:

- Add a vertex for each (L-1)-mer in the L-spectrum.

- Add k edges between two (L-1)-mers if their overlap has length L-2 and the corresponding L-mer appears k times in the L-spectrum.

An example de Bruijn graph construction is shown below.




<img width="1042" height="548" alt="image" src="https://github.com/user-attachments/assets/76e5a4ff-c768-4f58-8548-949204083272" />
