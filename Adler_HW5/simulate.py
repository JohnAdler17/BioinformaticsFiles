"""
So far, this program reads in the data from the command line arguments 
and if the correct number of arguments is given, it opens and reads in 
the sequence from the fasta file.

You will need to add the functionality to sample reads from the sequence,
add "noise" to the reads as described in the homework assignment,
and write the reads out to a file named reads.txt

@author: cwelsh1
"""
import sys
import random



#This function takes in the name of a FASTA file, and returns
#a single string that contains the entire DNA/RNA sequence found in the file
#Parameters: filename - a string that contains the name of a FASTA file
#Returns: a string
def read_fasta(filename):
    sequence = ""
    file = open(filename, "r")
    for line in file:
        line = line.rstrip()
        if line.startswith(">"):
            sequence = sequence
        else:
            for char in line:
                if char.upper() != "N" and char.upper() != " ":
                    sequence = sequence + char.upper()
                else:
                    sequence = sequence.lstrip()
                    sequence = sequence.rstrip()
    file.close()
    return sequence

#This function reads in the arguments, or prints an error message if the
#incorrect number of arguments was given
#Parameters: None
#Returns: 4 values consisting of the FASTA filename, coverage, length, 
#and error rate    
def parse_arguments():
    if len(sys.argv) != 5:
        print("Incorrect number of parameters")
        sys.exit(-1)
    else:
        filename = sys.argv[1]
        coverage = int(sys.argv[2])
        length = int(sys.argv[3])
        error = float(sys.argv[4])

    return filename, coverage, length, error

def main():
    '''
    Runs simulation program.
    '''
    filename, coverage, length, error = parse_arguments()

    #I used these to hard code in the file, coverage, length, and error rate
    """
    filename = "sample.fasta"
    coverage = 12
    length = 50
    error = 0.01
    """
    #seq = "ACCCTATTATTAAGGG"
    seq = read_fasta(filename)

    N = int((coverage * len(seq))/length)
    
    newfile = open("reads.txt", "w+")
    for i in range(N):
        j = random.randint(0, len(seq) - length)
        read = seq[j:j+length]
        #print(read, "started at position ", j)
        newread = ''
        for k in range(len(read)):
            errorprob = error * 100
            genrand = random.randrange(1, 100, 1)

            #checks every position k in every read for an error
            if errorprob >= genrand:
                #print("error generated at position ", k)
                randread = random.randrange(0, 2, 1)
                #if error exists, a replacement character is chosen depending
                #on the character that is replaced
                if read[k] == 'A':
                    
                    if randread == 0:
                        newread = newread + 'C'
                    elif randread == 1:
                        newread = newread + 'T'
                    elif randread == 2:
                        newread = newread + 'G'
                        
                elif read[k] == 'C':
                    
                    if randread == 0:
                        newread = newread + 'A'
                    elif randread == 1:
                        newread = newread + 'T'
                    elif randread == 2:
                        newread = newread + 'G'
                        
                elif read[k] == 'T':
                    
                    if randread == 0:
                        newread = newread + 'C'
                    elif randread == 1:
                        newread = newread + 'A'
                    elif randread == 2:
                        newread = newread + 'G'
                        
                elif read[k] == 'G':
                    
                    if randread == 0:
                        newread = newread + 'C'
                    elif randread == 1:
                        newread = newread + 'T'
                    elif randread == 2:
                        newread = newread + 'A'
            else:
                newread = newread + read[k]
        newfile.write("%s\n" % newread)
        #print("Newread: ", newread)
    newfile.close()
    
    
if __name__ == '__main__':
    main()
