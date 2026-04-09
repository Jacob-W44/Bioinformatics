# Set Working Directory to Folder with FASTA file in it.
setwd("C:/Users/jacob/Downloads/BioinformaticsClass/Bioinformatics/Midterm2_TakeHome")

# Load in necessary libraries to determine protein.
library(GenomicAlignments)
library(UniprotR)
library(protti)
library(r3dmol)
library(Biostrings)


# Read the specific alignment from FASTA file into DNA String set
homo.sapien <- readDNAStringSet("metazoa_alignment.gene.fasta") [[5]] # Homo sapiens were 5th alignment.

# There were gaps in the alignment, so this code removes the gaps to allow for translation.
homo.sapien <- gsub("[^ATGCatgc]", "", homo.sapien)

# This code sets the new alignment back into a string format after removing the spaces.
homo.sapien.string <- DNAString(homo.sapien)

# Translate the newly strung DNA alignment from previous code.
homo.sapien.prot <- Biostrings::translate(homo.sapien.string, if.fuzzy.codon = "solve")

# Convert amino acid translation into a string set to write as new FASTA file.
homo.sapien.prot <- AAStringSet(homo.sapien.prot)

# Create new FASTA file with AA sequence for searching on UniProt
Biostrings::writeXStringSet(homo.sapien.prot, filepath = "homo.sapien.prot", format = "fasta")

# For question 9, I used UniProt to compare my AA string set to potential proteins, 
# and the best match was gene POLG, which is a human DNA polymerase subunit (gamma-1)

# Gather information from UniProt to R studio using accession number of protein.
UniProt.hs <- fetch_uniprot("P54098")

# Gather Gene Ontology information from UniProt.
# The vector wasn't working so I just reused the accession number.
Prot.GO.Info <- GetProteinGOInfo("P54098") 

# Create Plot using the obtained Gene Ontology Information.
GO.Plot <- PlotGoInfo(Prot.GO.Info)
