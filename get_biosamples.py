from Bio import Entrez
import tqdm
import pandas as pd
import numpy as np

#from https://github.com/annacprice/ena-fastq-fetch


# Get SRA ID from GenBank accession
def get_sra(genbank_accession):
    """Access the sequence read archive accession for a GenBank consensus sequence"""
    Entrez.email = "theo@theo.io"
    handle = Entrez.efetch(db="nucleotide",
                           id=genbank_accession,
                           rettype="gb",
                           retmode="xml")
    record = Entrez.read(handle)
    handle.close()

    # Find DBLink to Sequence Read Archive
    # Iterate through GBXref_dbname and GBXref_id to find SRA
    for xref in record[0]["GBSeq_xrefs"]:
        print(xref)
        if xref["GBXref_dbname"] == "Sequence Read Archive":
            return xref["GBXref_id"]
    return None


def get_biosample_from_list(list_of_accessions):
    joined_list = ",".join(list_of_accessions)
    Entrez.email = "theo@theo.io"
    handle = Entrez.efetch(db="nucleotide",
                           id=joined_list,
                           rettype="gb",
                           retmode="xml")
    record = Entrez.read(handle)
    handle.close()
    results = {}

    # Find DBLink to Biosample
    for item in record:
        if "GBSeq_xrefs" in item:
            for xref in item["GBSeq_xrefs"]:
                if xref["GBXref_dbname"] == "BioSample":
                    results[item["GBSeq_locus"]] = xref["GBXref_id"]
                # print("success")
    return results


batch_size = 100
batch = []
all_results = {}

accessions_list = pd.read_csv("covspectrum.csv").genbankAccession.values

for accession in tqdm.tqdm(accessions_list):
    accession = accession.strip()
    batch.append(accession)
    if len(batch) == batch_size:
        print("s")
        results = None
        while results == None:
            try:
                results = get_biosample_from_list(batch)
            except Exception as e:
                print("There was a problem with", e, batch)
        print("e")
        all_results.update(results)
        batch = []

output_file = open("biosamples.txt", "wt")
for accession, sra in all_results.items():
    output_file.write(accession + "\t" + sra + "\n")
