"""
assembly.py assembles a De Bruijn Graph from an input file of reads
author: John Adler
"""
import sys
import random

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
    
    #filename = "reads.txt"
    #lenkmer = 5
    
    
    file = open(filename, "r+")
    ListVertices = [] #list that will hold the number of vertices in the De Bruijn Graph
    EdgeDictionary = {} #dictionary that will hold the number of edges in the De Bruijn Graph

    VertexEdge = []
    Sequence = ''
    ListContigs = []

    for line in file:
        
        line = line.strip()
        #contigStart = random.randrange(2, len(line)- 2, 1)
        """
        contig = line
        ListContigs.append(contig)
        """
        
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
                VertexEdge.append(kmer)

            if kmer in EdgeDictionary:
                EdgeDictionary[kmer] = EdgeDictionary[kmer] + 1

    #print(len(VertexEdge))
    
    while (len(EdgeDictionary) != 0):
        #if the VertexEdge list of kmers is empty, then break
        if len(VertexEdge) == 0:
            break
        
        contig = '' #an empty string to contain the constructed contig
        randkmer = random.randint(0, abs(len(VertexEdge)-1)) #a value that picks a random kmer to start the contig with
        startkmer = VertexEdge[randkmer]
        
        #print(randkmer)
        #randcontiglength = random.randint(0, len(VertexEdge)-1)
        
        contig = startkmer #contig starts with random kmer
        VertexEdge.remove(startkmer) #removes start kmer from list

        #removes start kmer from edge dictionary if the value in dictionary is 0
        EdgeDictionary[startkmer] = EdgeDictionary[startkmer] - 1
        if EdgeDictionary[startkmer] == 0:
            EdgeDictionary.pop(startkmer)

        nextletterindex = 0
        for i in range(0, len(VertexEdge), 1):
            if nextletterindex == 4:
                    break
            
            if nextletterindex == 0:
                nextletter = 'A'
            elif nextletterindex == 1:
                nextletter = 'T'
            elif nextletterindex == 2:
                nextletter = 'C'
            elif nextletterindex == 3:
                nextletter = 'G'

            
            if contig[(i+1):]+nextletter not in VertexEdge:
                nextletterindex = nextletterindex + 1
                continue
            
            #print(contig)
            #print(contig[(i+1):])
            indexnextkmer = VertexEdge.index(contig[(i+1):]+nextletter)
            nextkmer = VertexEdge[indexnextkmer]
            VertexEdge.remove(nextkmer)
            EdgeDictionary[nextkmer] = EdgeDictionary[nextkmer] - 1
            if EdgeDictionary[nextkmer] == 0:
                EdgeDictionary.pop(nextkmer)
            #print(nextkmer)
            contig = contig + nextletter
            
            if nextletterindex == 4:
                nextletterindex = 0
        ListContigs.append(contig)
    
    for i in range(len(ListContigs)):
        Sequence = Sequence + ListContigs[i]
        #print(ListContigs[i])
    lenContigs = len(Sequence)/2

    ListLargeContigs = []

    newfile = open("contigs.txt", "w+")
    for i in range(len(ListContigs)):
        contig = ListContigs[i]
        newfile.write("%s\n" % contig)
    
    while (lenContigs != 0):
        if len(ListContigs) == 0:
            break
        largestcontig = ListContigs[0]
        for i in range(1, len(ListContigs)):
            if len(largestcontig) < len(ListContigs[i]):
                largestcontig = ListContigs[i]
        lenContigs = lenContigs - len(largestcontig)
        #if the current largest contig makes the lenContigs negative, then it continues the loop
        if lenContigs < 0:
            lenContigs = lenContigs + len(largestcontig)
            ListContigs.remove(largestcontig)
            continue
        ListContigs.remove(largestcontig)
        ListLargeContigs.append(largestcontig)

    """
    for i in range(len(ListLargeContigs)):
        print(ListLargeContigs[i])
    """
    
    N50 = len(ListLargeContigs[0])
    for i in range(1, len(ListLargeContigs)):
        if N50 > len(ListLargeContigs[i]):
            N50 = len(ListLargeContigs[i])

    L50 = len(ListLargeContigs)
        
    #print("Sequence: ", Sequence)
    """     
    for i in range(len(ListVertices)):
        print(ListVertices[i])
    """
    #print("Number of vertices in graph: ", len(ListVertices))
    #print("Number of edges in graph: ", len(EdgeDictionary))

    print("N50 Score: ", N50)
    print("L50 Score: ", L50)

    
if __name__ == '__main__':
    main()
