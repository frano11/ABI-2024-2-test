import csv
import pandas as pd


import Bio.Entrez as ez
ez.email = "c.knorr@cq-bildung.de"
from Bio import SeqIO

from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq

dna_h2b_seqs = [] # List to store fetched ID sequences

with ez.esearch(db="nuccore", term='"H2B"[Gene Name] AND "Mus musculus"[Organism]', retmax=20) as query:
    esearch11 = ez.read(query)
    print("=========== TOTAL EFETCH INFO OF IDs FROM NUCLEOTIDE - GENE: H2B / ORGANISM: MUS MUSCULUS -  =============")
    print(esearch11)
    print()
    print("=========== FETCHED IDs FROM IdList  =============")
    print(esearch11["IdList"])
    print()
    for ids in esearch11["IdList"]:
        with ez.efetch(db="nuccore", id=ids, rettype="fasta", retmode="text") as handle:
            efetch11 = SeqIO.read(handle, "fasta")
            dna_h2b_seqs.append(efetch11)
            print("=========== LOOPED EFETCH IDs FROM NUCLEOTIDE - GENE: H2B / ORGANISM: MUS MUSCULUS -  =============")
            print(efetch11)
            print()
            print("=========== USING .SEQ =============")
            print(efetch11.id)
            print(repr(efetch11.seq))
            print(len(efetch11))
            print()
print()
print("======== Visualization of the list dna_h2b_seqs (CHAOTIC!!) ===========")
print()
print(dna_h2b_seqs) # If we print the list as it is, we'll get a mass difficult to read.
print()
# print(dna_h2b_seqs.seq) # it gives an error
# To see the IDs and Seqs from the list, we must create a loop. See the next script right below:
print("========= USING A FOR LOOP FOR PRINTING THE CONTENT FROM LIST dna_h2b_seqs (ORGANIZED!!) ===========")
# The only way to get the info from all IDs in the list, we must use the list we created in the loop (dna_h2b_seq1), ...
# ... and make a for loop.
print()
for seqs in dna_h2b_seqs:
    #print(ID_)
    print(seqs)
    print()


print()
print("============ PLAYING WITH .ID and .SEQ ================")
print()
dna_h2b_id1 = (efetch11.id)
dna_h2b_seq1 = (efetch11.seq)
print(dna_h2b_id1)
print(dna_h2b_seq1)
print(repr(dna_h2b_seq1), "representation display")
print()
print(efetch11)
# Here, when the variable efetch11 is outside the loop can only show the info of 1 ID, which is the last one that overwrote the previous ID.
# The same if we use the list (dna_h2b_seq1) outside a loop. It will only show the last ID.

print()
dna = Seq(str(dna_h2b_seq1))  
print(repr(dna), "DNA 5' --> 3'")
print(repr(dna.complement()), "DNA 3' --> 5' complement")
print(repr(dna.transcribe()), "mRNA 5' --> 3'")
#print(dna.translate(), "Protein sequence with * (STOP)")
# to translate into a protein, the DNA/mRNA should be multiple of three, otherwise the tranlation won't be accurate.
#print(dna.translate(to_stop=True), "Protein sequence ")
#translated1 = dna.translate()  # translate can use either mRNA strand or DNA coding strand. There's no need to transcribe first.
print()
print(len(dna_h2b_seq1), "bp linear DNA or mRNA")
#print(len(translated1), "amino acids total")
print()
print()

print("================= DISPLAYING THE SEQ IDs AS DNA, COMPLEMENTARY AND  MRNA IN THE LIST dna_h2b_seqs USING A FOR LOOP =====================")
print()
for seqs in dna_h2b_seqs:  # dna_h2b_seqs --> is a list with the two SeqRecords
    seq = seqs.seq         # Access the sequence from SeqRecord object
    print(seqs, "DNA 5' --> 3' (DNA CODING STRAND)")
    print(repr(seq.complement()), "DNA 3' --> 5' complement (DNA TEMPLATE STRAND)" )
    print(repr(seq.transcribe()), "mRNA 5' --> 3'")
    #print((seq.translate(), "Protein"))  # REMEMBER that the DNA/mRNA sequence is not multiple of three, thus the protein seq is not accurate.
    print()
print()
print()
# print("================= SAVING DNA - COMP. DNA - MRNA TO A TXT FILE =====================")
# with open("240403_ex12-Pipeline_Nuccore1.txt", "w") as file:
#     for seqs in dna_h2b_seqs:  # dna_h2b_seqs --> is a list with the two SeqRecords
#         seq = seqs.seq  # Access the sequence from SeqRecord object
#         # Write the sequence, complement, and mRNA transcription to the file
#         file.write(str(seqs) + "DNA 5' --> 3'\n")
#         file.write(str(repr(seq.complement()) + "DNA 3' --> 5' complement \n"))
#         file.write(repr(seq.transcribe()) + " mRNA 5' --> 3'\n")
#         file.write("\n")  # Add a newline for separation
# print()
# 
# print("================= SAVING DNA - COMP. DNA - MRNA TO A CSV FILE USING IMPORT CSV =====================")
# # Open a CSV file for writing
# with open("240403_ex12-Pipeline_Nuccore1_CSV.csv", "w", newline="") as csvfile:
#     writer = csv.writer(csvfile)
#     # Write the header row
#     writer.writerow(["Sequence", "Complement", "mRNA"])
# 
#     # Iterate over each SeqRecord in dna_h2b_seqs
#     for seqs in dna_h2b_seqs:
#         seq = seqs.seq  # Access the sequence from SeqRecord object
#         # Write the sequence, complement, and mRNA transcription to the CSV file
#         writer.writerow([str(seqs), str(seq.complement()), repr(seq.transcribe())])
# 
# 
# print()
# 
# print("================= SAVING DNA - COMP. DNA - MRNA TO A CSV FILE USING IMPORT PANDAS =====================")
# # Create a list to store the data
# data = []
# 
# # Iterate over each SeqRecord in dna_h2b_seqs
# for seqs in dna_h2b_seqs:
#     seq = seqs.seq  # Access the sequence from SeqRecord object
#     # Append the sequence, complement, and mRNA transcription to the data list
#     data.append([str(seqs), str(seq.complement()), repr(seq.transcribe())])
# 
# # Create a DataFrame from the data
# df = pd.DataFrame(data, columns=["Sequence", "Complement", "mRNA"])
# 
# # Save the DataFrame to a CSV file
# df.to_csv("240403_ex12-Pipeline_Nuccore1_PANDAS.csv", index=False)

# print()
# print("=========== EGQUERY OF TERM=H2B =============")
# print()
# with ez.egquery(term="H2B") as handle:
#     egquery1 = ez.read(handle)
#     for row in egquery1["eGQueryResult"]:
#         print(row["DbName"], row["Count"])
# print()

print()
print("=========== ELINK OF IDs FROM H2B LIST =============")
print()
elink_list1 = esearch11["IdList"]
print(elink_list1)
print()
#elink_dict1[databases]={}
for ids in elink_list1:
    with ez.elink(db="gene", dbfrom="nuccore", id=ids) as query:
        elink1 = ez.read(query)
        print("From Nucleotide ID: ",elink1[0]["IdList"],"\nTo Gene IDs: ",elink1[0]["LinkSetDb"][0]["Link"])
        print()
    with ez.elink(db="protein", dbfrom="nuccore", id=ids) as query:
        elink1 = ez.read(query)
        print("From Nucleotide ID: ",elink1[0]["IdList"],"\nTo Protein IDs: ",elink1[0]["LinkSetDb"][0]["Link"])
        print()
    with ez.elink(db="pubmed", dbfrom="nuccore", id=ids) as query:
        elink1 = ez.read(query)
        print("From Nucleotide ID: ",elink1[0]["IdList"],"\nTo Pubmed IDs: ",elink1[0]["LinkSetDb"][0]["Link"])
        print()
#     with ez.elink(db="Biosample", dbfrom="nuccore", id=ids) as query:
#         elink1 = ez.read(query)
#         print("From Nucleotide ID: ",elink1[0]["IdList"],"\nTo BioSample IDs: ",elink1[0]["LinkSetDb"][0]["Link"])
#         print()
#     with ez.elink(db="SRA", dbfrom="nuccore", id=ids) as query:
#         elink1 = ez.read(query)
#         print("From Nucleotide ID: ",elink1[0]["IdList"],"\nTo SRA IDs: ",elink1[0]["LinkSetDb"][0]["Link"])
#         print()
    with ez.elink(db="pmc", dbfrom="nuccore", id=ids) as query:
        elink1 = ez.read(query)
        print("From Nucleotide ID: ",elink1[0]["IdList"],"\nTo PMC IDs: ",elink1[0]["LinkSetDb"][0]["Link"])
        print()

# print("========= SOME SEQUENCES HK2 - MUS MUSCULUS IMPORTED FROM FASTA FILE =============")
# print()
# for seqs in SeqIO.parse("ex9_nuccore_hk2_mouse_fasta", "fasta"):
#     print("ID: ",seqs.id)
#     print("Descrition: ",seqs.description)
#     print((seqs.seq), "DNA 5' --> 3'(DNA CODING STRAND)")
#     print((seqs.seq.complement()), "DNA 3' --> 5' complement (DNA TEMPLATE STRAND)" )
#     print((seqs.seq.transcribe()), "mRNA 5' --> 3'")
#     print(len(seqs),"bp linear DNA or mRNA")
#     seq_obj = (seqs)
#     gc_content = gc_fraction(seq_obj)
#     print("GC content:", gc_content*100,"%")
#     print()