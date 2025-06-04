# Example usage of the GenomeAssembler module:

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
