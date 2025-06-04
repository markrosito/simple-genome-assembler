# Basic Genome Assembler - v1.0

A simple genome assembler written in Python as part of the final project for the course **BioQuants** at Politecnico di Torino.

## 📘 Description

This module reads DNA reads from a file, computes pairwise overlaps, and assembles the full genome sequence based on those overlaps.

The tool operates under two simplifying assumptions:
1. There are no sequencing errors (all overlaps are perfect matches).
2. No reads are nested within others (no full containment).

## 🚀 Features

- Read DNA fragments from a text file.
- Calculate all pairwise overlaps.
- Automatically determine the read order.
- Assemble the full genome sequence.
- Save the assembled genome to a file.

## 🧪 Example Input Format

Each line in the input file should have this format:

1:ATGCG... 
2:AGGCG... 
3:TGAAG...


## 🖥️ Requirements

- Python 3.6 or higher  
- No external libraries required (uses only Python standard library)

## 🧠 Example of Use

```python
from GenomeAssembler import *

fileName='genome_assembly-input.txt' 
reads=readDataFromFile(fileName)

print("\nMean length of sequences:")
print(meanLength(fileName), "\n") 

overlaps=getAllOverlaps(reads)
print("Overlaps: ") 
prettyPrint(overlaps)

name=findFirstRead(overlaps) 
order=findOrder(name,overlaps) 
genome=assembleGenome(order,reads,overlaps) 
print("Assembled genome: ")
print(genome, "\n") 

saveGenomeToFile(genome, 'genome.txt')