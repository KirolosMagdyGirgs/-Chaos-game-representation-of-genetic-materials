# The integer chaos game (iCGR) representation of DNA sequences: encoding and decoding
import math


# ------------------------------------------------------------------------------
# Helper function: get the nucleotide based on the signs of integers x and y
# ------------------------------------------------------------------------------
def getNucleotide(x, y):
    if (x > 0 and y > 0):
        nucleotide = 'A'
    elif (x > 0 and y < 0):
        nucleotide = 'G'
    elif (x < 0 and y > 0):
        nucleotide = 'T'
    elif (x < 0 and y < 0):
        nucleotide = 'C'
    else:
        nucleotide = 'N'

    return nucleotide


# ------------------------------------------------------------------------------
# Helper function: get iCGR vertices for given signs of integers x and y
# ------------------------------------------------------------------------------
def getCGRVertex(x, y):
    Nx = 0;
    Ny = 0;

    if (x > 0 and y > 0):
        Nx = 1
        Ny = 1
    elif (x > 0 and y < 0):
        Nx = 1
        Ny = -1
    elif (x < 0 and y > 0):
        Nx = -1
        Ny = 1
    elif (x < 0 and y < 0):
        Nx = -1
        Ny = -1
    else:
        Nx = 0
        Ny = 0

    return Nx, Ny


# ------------------------------------------------------------------------------
# iCGR Encoding: encode a DNA sequence into three integers (iCGR encoding)
# Input: a DNA sequence
# Outputs: three integers, x_n, y_n, and n
# ------------------------------------------------------------------------------
def encodeDNASequence(seq):
    A = [1, 1]
    T = [-1, 1]
    C = [-1, -1]
    G = [1, -1]
    a = 0
    b = 0
    x = []
    y = []
    n = len(seq)

    if seq[0] == 'A':
        a = int(A[0])
        b = int(A[1])
    elif seq[0] == 'T':
        a = int(T[0])
        b = int(T[1])
    elif seq[0] == 'C':
        a = int(C[0])
        b = int(C[1])
    else:
        a = int(G[0])
        b = int(G[1])

    x.append(a)
    y.append(b)

    for i in range(1, n):
        if seq[i] == 'A':
            a = int(x[i - 1]) + int(math.pow(2, i))
            b = int(y[i - 1]) + int(math.pow(2, i))
        elif seq[i] == 'T':
            a = int(x[i - 1]) - int(math.pow(2, i))
            b = int(y[i - 1]) + int(math.pow(2, i))
        elif seq[i] == 'C':
            a = int(x[i - 1]) - int(math.pow(2, i))
            b = int(y[i - 1]) - int(math.pow(2, i))
        else:
            a = int(x[i - 1]) + int(math.pow(2, i))
            b = int(y[i - 1]) - int(math.pow(2, i))

        x.append(a)
        y.append(b)

    x_n = int(x[n - 1])
    y_n = int(y[n - 1])

    return x_n, y_n, n


# ------------------------------------------------------------------------------
# iCGR Decoding: decode three integers to a DNA sequence
# Inputs: three integers, x_n, y_n, and n
# Outputs: the DNA sequence represented by the three integers
# ------------------------------------------------------------------------------
def decodeDNASequence(x_n, y_n, n):
    x = [0] * n
    y = [0] * n
    x[n - 1] = x_n
    y[n - 1] = y_n

    seq = []
    for i in range(n - 1, -1, -1):
        nt = getNucleotide(x[i], y[i])
        seq.insert(i - n + 1, nt)
        Nx, Ny = getCGRVertex(x[i], y[i])
        x[i - 1] = int(x[i]) - int(math.pow(2, i) * Nx)
        y[i - 1] = int(y[i]) - int(math.pow(2, i) * Ny)
    else:
        DNASeq = ''.join(seq)
        return DNASeq



# The use case of integer chaos game (iCGR) encoding and decoding DNA sequences

# Homo sapiens SMPX gene for alternative protein SMPX, isolate 15610, GenBank: HF583935.1
# https://www.ncbi.nlm.nih.gov/nuccore/HF583935.1?report=fasta
DNASeq = 'ATGGAGAATGATATGGCAATGTGCCTAACGATTTTGATGAAAAGTTTCCCAAGCTACTTCCTACAGTATT'
print('Original DNA sequence:', DNASeq)

# 1. Encoding DNA sequence by integer CGR (iCGR)
x_n, y_n, n = encodeDNASequence(DNASeq)
print('Encoded tri-integers:', x_n, y_n, n)

# 2. Decoding DNA sequence from tri-integers
DNAseq_decoded = decodeDNASequence(x_n, y_n, n)
print('Decoded DNA sequence:', DNAseq_decoded)

if DNASeq == DNAseq_decoded:
    print('Successfully encoding and decoding by iCGR.')
else:
    print('Unsuccessfully encoding and decoding.')