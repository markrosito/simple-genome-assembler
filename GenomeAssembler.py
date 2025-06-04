#######################################################################################
# Basic Genome Assembler Module - v.1.0                                               #
# Author: s319706 Marco Rosella - Politecnico di Torino                               #
# Email: s319706@studenti.polito.it                                                   # 
# This module is part of the exam of the course "BioQuants" at Politecnico di Torino. #
# Date: 2025-05-15                                                                    #
# Python version: 3.6+                                                                #
#######################################################################################

"""Genome Assembler:
This simple module provides some functions to read genome sequences from a file,
calculate overlaps between sequences, and assemble a genome from those sequences.
It includes functions to read data, calculate mean sequence length,
find overlaps, determine the order of reads, and assemble the genome.
It also provides a function to save the assembled genome to a file.
For the the develop of this tool we make two important simplifying assumptions:
    1. There are no sequencing errors, so that overlaps between reads are
       perfect matches. 
    2. No reads are nested in other reads, implying that overlaps are always
       of this type:  
            XXXXXXXXXX 
                  XXXXXXXXXXX 
       and never of this type: 
            XXXXXXXXXXXXXX 
               XXXXXXXX
"""

def readDataFromFile( fileName ):
    """
    This function reads data from a file and returns the content in a dictionary.
    Each line in the file should contain a name and a sequence
    separated by a semicolon (:).
    The function returns a dictionary for each line, where the key
    is the name and the value is the corresponding sequence.
    
    Params:
        fileName (str): The name of the file to read from, in the form "fileName.txt".

    Returns: 
        data (dict): A dictionary containing names and sequences.
    """
    
    data = {}
    with open(fileName, 'r') as file:
        for line in file:
            name, seq = line.strip().split(':')
            data[name] = str(seq)
    return data

def meanLength(fileName):
    """
    This function calculates the mean length of the sequences in a .txt file,
    where each line contains a name and a sequence separated by a semicolon (:).
    
    Params: 
        fileName (str): The name of the file to read from, in the form "fileName.txt".

    Returns: 
        meanLength (float): The mean length of sequences.
    """
    data = readDataFromFile(fileName)
    totalLength = sum(len(seq) for seq in data.values())
    meanLength = totalLength / len(data)
    return meanLength

def getOverlap(left, right):
    """
    This function finds the overlapping characters between two strings.
    In the DNA case, it checks if the 3' end of the left read
    overlaps the 5' end of the right read.
    It checks all possible combinations and it returns the overlapping sequence
    in the form of a string, if any overlap exists.
    
    Params:
        left (str): The left sequence.
        right (str): The right sequence.
    
    Returns:
        overlap (str): The overlap sequence.
    """
    overlap = ''
    for i in range(min(len(left), len(right))-1):
        if right[:i+1] == left[-i-1:]:
            overlap = right[:i+1]
    return overlap

def getAllOverlaps(reads):
    """
    This function finds all overlaps between all the sequences in a
    dictionary in the form returned from the readDataFromFile() function.
    
    Params:
        reads (dict): A dictionary as produced by readDataFromFile().
    
    Returns:
        d (dict): A dictionary of dictionaries with overlaps length
    """
    d = {}
    for i in reads.keys():
        d[i] = {}
        for j in reads.keys():
            if i != j:
                d[i] = {**d[i], j : len(getOverlap(reads[i], reads[j]))}
    return d

def prettyPrint(overlaps):
    """
    This function prints in a matrix-like way the overlaps dictionary
    found by getAllOverlaps(). You can read the overlaps considering
    a row index as the first read and a column index as the second read.
    
    Params:
        overlaps (dict): The overlaps dictionary returned by getAllOverlaps().

    Returns:
        None: This function prints the overlaps in a formatted way.
    """
    names = list(overlaps.keys())
    names.sort()
    print(' ' * 3, end = ' ')
    for i in range(len(names)):
        print('%3s' % names[i], end = ' ')
    
    for i in names:
        print('\n%3s' % i, end = ' ')
        for j in names:
            if i == j:
                print('  -', end = ' ')
            else:
                print('%3d' % overlaps[i][j], end = ' ')
    print('\n')
        
def findFirstRead(overlaps):
    """
    This function finds the first read in the overlaps dictionary.
    If the overlap length is less than 3, we does not consider it significant.
    It returns the name of the first read of the genome sequence.
    
    Params:
        overlaps (dict): The overlaps dictionary to check,
            as produced by getAllOverlaps().

    Returns: 
        (str): The name of the first read.
    """
    for i in overlaps.keys():
        c = 0    
        for j in overlaps.keys():
            if i != j and overlaps[j][i] < 3:
                c = c+1
            if c == len(overlaps.keys()) - 1: 
                return str(i)
            
def findKeyForLargestValue(d):
    """
    Find the key for the largest value in a dictionary.
    
    Params:
        d (dict): The dictionary to check.
    
    Returns:
        (str): The key with the largest associated value.
    """
    return max(d, key=d.get)

def findOrder(name, overlaps):
    """
    This function finds the order of the reads based on overlaps dictionary
    as returned by getAllOverlaps().
    This function starts from the first read and finds the next read,
    which is the one that overlaps the most with the current read.
    It continues to find the next read until it finds a read
    that has already been added to the order or the overlap length
    is less than 3.
    
    Params:
        name (str): The name of the read to check.
        overlaps (dict): The overlaps dictionary to check.
    
    Returns:
        order (list): The order of reads.
    """
    order = [name]
    while True:
        nextRead = findKeyForLargestValue(overlaps[order[-1]])
        if nextRead in order or overlaps[order[-1]][nextRead] < 3:
            break
        order.append(nextRead)
    return order

def assembleGenome(readOrder, reads, overlaps):
    """
    This function assembles the genome based on the list with the reading order,
    as given by findOrder().
    It concatenates the reads in the order specified by readOrder,
    using the overlaps dictionary to determine where to cut the reads.
    
    Params:
        readOrder (list): The order of reads to assemble.
        reads (dict): The dictionary of reads.
        overlaps (dict): The overlaps dictionary.
    
    Returns:
        genome (str): The assembled genome.
    """
    genome = reads[readOrder[0]]
    for i in range(len(readOrder)-1):
        genome += reads[readOrder[i+1]][overlaps[readOrder[i]][readOrder[i+1]]:]
    return genome


def saveGenomeToFile(genome, fileName):
    """
    Save the assembled genome to a .txt file.
    
    Params:
        genome (str): The assembled genome.
        fileName (str): The name of the file to save to, in the form "fileName.txt".
    
    Returns:
        None: The function writes the genome to the specified file.
    """
    with open(fileName, 'w') as file:
        file.write(genome)



