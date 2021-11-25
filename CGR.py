import collections
from matplotlib import cm
import pylab
import math 
from urllib.request import urlopen

#called by chaos graph
#Counting K_mers
def count_kmers(sequence, k):
    #any new key introduced to the dictionary 
    # will have default value zero
    d = collections.defaultdict(int)
    #loop over number of possible kmers ((L - K) + 1) or (L -(K-1))
    #reporting number of occurence for each kmer 
    for i in range(len(sequence)-(k-1)):
        #introducing new keys , setting their values 
        #updating old keys' value
        d[sequence[i:i+k]] +=1
        #{"att":3,"acg":1,....}
    for key in d.keys():
        if "N" in key:
            del d[key]
    return d

#called by chaos graph
#Calculating Prob of each k-mer [agt:25%, etc... ]
def probabilities(kmer_count, k,squence):
     #any new key introduced to the dictionary 
     # will have default value 0.0
    kmer_prop = collections.defaultdict(float)
    squenceLength = len(squence)
    for key, value in kmer_count.items():
        kmer_prop[key] = float(value) / (squenceLength - k + 1)
    return kmer_prop


#called by chaos graph
#calculating CGR Arrary
def chaos_game_representation(probabilities, k):
    try:
        chaos = []
        array_size = int(math.sqrt(4**k))
        for i in range(array_size):
            chaos.append([0]*array_size)
        #initializing position of current location and determining the last location of the array
        maxx = array_size
        maxy = array_size
        posx = 1
        posy = 1
        for key, value in probabilities.items():
            for char in key:
                # in respect of
                #  A = upper left quadrant 
                #  C = lower left quadrant  then (maxY/2)
                #  G = upper right quadrant then (maxX/2)
                #  T = lower right quadrant then (maxX/2 or maxY/2)
                if char == "T":
                    posx += maxx / 2
                    posy += maxy / 2
                elif char == "C":
                    posy += maxy / 2
                elif char == "G":
                    posx += maxx / 2
                maxx = maxx / 2
                maxy /= 2
            #chaos array that is 2 dimentional and will be showed in graphical representation
            chaos[int(posy-1)][int(posx-1)] = value
            maxx = array_size
            maxy = array_size
            posx = 1
            posy = 1
        return chaos
    except:
        print("\nThe size of k-mer was too big or too small \nKmer is set to the default value 4\n")
        return 3




def graph(chaos_kx,kmerSize):
    pylab.figure(1)
    pylab.title('Chaos game representation for {}-mers'.format(kmerSize))
    pylab.imshow(chaos_kx, interpolation='nearest', origin="lower")
    pylab.axis("off")
    pylab.show()        
 
#called by chaosGfraph
def baseForChaosGraph(sequence,kmerSize):
    fx = count_kmers(sequence, kmerSize)
    fx_prob = probabilities(fx,kmerSize,sequence)
    chaos_kx = chaos_game_representation(fx_prob, kmerSize)
    if chaos_kx == 3:
        return 4
    else:
        graph(chaos_kx,kmerSize)
        
#called by the main program
def chaosGraph(data,modifiedKmer=None):
    if  not bool(modifiedKmer):
        iterator = 0
        try:
            numOfKmers=int(input("How many kmers would you like to represent ? \n"))
            while iterator<numOfKmers:
                kmerSize = int(input("Enter size of k-mer {} \n".format(iterator+1)))
                chaos_kx = baseForChaosGraph(data,kmerSize)
                #if the kemer set by the user was not suitble 
                #it's overwritten with 4 and passed again to the function
                if chaos_kx == 4:
                    baseForChaosGraph(data,4)
                iterator += 1
        except:
            print("number of kmers is undefined and set to one kmer with size 4\n")
            baseForChaosGraph(data,4)
    
    
        
    

#revers translation of protien back into dna 
def reverseTranslation(protienFile,filetype=None):     
    #reading the file (Either it's local or was created out of the url)
    with open(protienFile,"r") as f:
        fil = f.read()
        fil = check_file_type(fil,filetype) #fasta or fastq 
        table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }
        allKeys = ""
        i = 0
        while(i < len(fil)):
            for key, value in table.items():
                #if the value of the AA is found return its corrusponding key (triplet bases
                if fil[i] == value:
                    allKeys = allKeys + "," + key
                    #break after finding the first occurence of the amino acid 
                    #so it won't be reported more than once
                    break
            i += 1
        return allKeys


#url method
def check_file_type(file=None,type=None,chekUrlIsUsed=None):
        f=""

        if chekUrlIsUsed ==1:
            with open(file,"r") as DnaFile:
                file = DnaFile.read()
        if type == "A":
                 f="".join(file.split("\n")[1:])
        elif type == "Q":
                 seq_reads=[]
                 data_list = file.split("\n") # a list of every line in the file
                 for i in range(1,len(data_list),4):
                            seq_reads.append(data_list[i]) # only seq reads
                 f="".join(seq_reads) #as a full text(joint_seq_reads)
        return f
        
def URL_method():
    url_user=input('''\nEnter your file's URL 
    \t(it is recommended to search by seq accession number in the database )\n''')
    response = urlopen(url_user)
    file_from_url = response.read().decode("utf-8", "ignore")
    return File_method(file_from_url,"Called by URL")

#file method
def File_method(file_from_comp=None,x=None):
    errorhandler="wrong path or name of the file \n"
    protein_file_holder=None    
    file_from_user=file_from_comp
    chekUrlIsUsed=0 #false
    while(True):
        #your file is local then give its name or path
        if x != "Called by URL":
            file_from_user= input ("Enter the (name or the path) of the file  ")
        else:
            #not local then set the chekking varible to 1
            #  to be used in turning the string into file
            chekUrlIsUsed = 1 #true   

        #cheking the file_from_user is not empty
        if bool(file_from_user):
            input_data_type = input("Is your data DNA or Protien \n D for DNA \n P for protien \n ").upper()  #DNA or Protein
            input_file_type = input("Is your file FASTA OR FASTQ \n A for FASTA \n Q for FASTQ \n").upper() #FASTA OR FASTQ 
            if input_data_type == 'D':
                try:
                        if chekUrlIsUsed == 1 :
                            with open("newDNAFile",'w+') as DnaFl:
                                dfile=DnaFl.write(file_from_user) 
                            data =check_file_type('newDNAFile',input_file_type,chekUrlIsUsed)
                            chaosGraph(data)
                            break
                        else:
                                 with open(file_from_user,"r") as f:
                                         s1 = f.read()
                                         data =check_file_type(s1,input_file_type)
                                         chaosGraph(data)
                                         break
                except:
                    print(errorhandler)
            elif input_data_type == 'P':
                    try:
                        #if your file came of an url the first condition is performed
                        if chekUrlIsUsed == 1: # true 
                            # writing the data from the string into a file 
                            # to be able to deal with it as in the local file 
                            with open("newProtienFile",'w+') as protienFl:
                                pfile=protienFl.write(file_from_user)
                                print(pfile)
                            protein_file_holder=reverseTranslation("newProtienFile",input_file_type)
                        else:
                            protein_file_holder=reverseTranslation(file_from_user,input_file_type)     
                    except:
                        print(errorhandler)
                    data = "".join(protein_file_holder.split(","))
                    chaosGraph(data)
                    break
            else:
                print("\n File not found \n")
                continue







#Reading the file from the user and defining its type
try:
    fileOrUrl=input('''please enter F for file or U for url: ''').upper()
    if fileOrUrl=="F":
        File_method()
    elif fileOrUrl=="U":
        URL_method()
    else:
        print("\nNot recognized input")
except:
    print("Somthing went wrong")
