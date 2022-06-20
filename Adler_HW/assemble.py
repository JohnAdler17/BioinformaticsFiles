"""
assembly.py assembles a De Bruijn Graph from an input file of reads
author: John Adler
"""
import sys

#This function reads in the arguments, or prints an error message if the
#incorrect number of arguments was given
#Parameters: None
#Returns: 2 values consisting of the filename and length of each k-mer    
def parse_arguments():
    if len(sys.argv) != 3:
        print("Incorrect number of parameters")
        sys.exit(-1)
    else:
        filename = sys.argv[1]
        lenkmer = int(sys.argv[2])

    return filename, lenkmer

def main():
    filename, lenkmer = parse_arguments()

    #I used these to hard code in the file and kmer length
    """
    filename = "sample_c12_r_50_e0.01.txt"
    lenkmer = 13
    """
    
    file = open(filename, "r+")
    ListVertices = [] #list that will hold the number of vertices in the De Bruijn Graph
    EdgeDictionary = {} #dictionary that will hold the number of edges in the De Bruijn Graph
    for line in file:
        line = line.strip()
        for i in range(len(line) - lenkmer + 1):
            kmer = line[i:i+lenkmer]
            left = kmer[:-1]
            right = kmer[1:]
            if left not in ListVertices:
                ListVertices.append(left)
                
            if right not in ListVertices:
                ListVertices.append(right)
                
            if kmer not in EdgeDictionary:
                EdgeDictionary[kmer] = 1

            if kmer in EdgeDictionary:
                EdgeDictionary[kmer] = EdgeDictionary[kmer] + 1
        
    """     
    for i in range(len(ListVertices)):
        print(ListVertices[i])
    """
    print("Number of vertices in graph: ", len(ListVertices))
    print("Number of edges in graph: ", len(EdgeDictionary))
if __name__ == '__main__':
    main()
