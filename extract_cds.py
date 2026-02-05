from Bio import Entrez, SeqIO
import csv
import time

# REQUIRED: put your email (NCBI requirement)
Entrez.email = "zakariyyayevrasul@gmail.com"

# List of genome accessions
accessions = [
    "MZ501051","MZ501080","MZ501087","MZ501069","MZ501101","MZ501090",
    "MZ501079","MZ501059","MZ501099","MZ501077","MZ501085","MZ501053",
    "MZ501092","MZ501107","MZ501097","MZ501070","MZ501088","MZ501095",
    "MZ501057","MZ501068","MZ501071","MZ501091","MZ501112","MZ501104",
    "MZ501110","MZ501072","MZ501108","MZ501076","MZ501105","MZ501109",
    "MZ501058","MZ501075","MZ501074","MZ501103","MZ501113","MZ501096",
    "MZ501089","MZ501052","MZ501066","MZ501065","MZ501067","MZ501050",
    "MZ501106","MZ501046","MZ501098","MZ501056","MZ501047","MZ501054",
    "MZ501062","MZ501060","MZ501111","MZ501102","MZ501083","MZ501093",
    "MZ501082","MZ501048","MZ501094","MZ501073","MZ501061","MZ501100",
    "MZ501063","MZ501084","MZ501086","MZ501081","MZ501078","MZ501064",
    "MZ501055","MZ501049","PQ850619","PQ850615","PQ850621","PQ850603",
    "PQ850627","PQ850596","PQ850632","PQ850598","PQ850625","PQ850629",
    "PQ850597","PQ850630","PQ850612","PQ850604","PQ850611","PQ850624",
    "PQ850622","PQ850594","PQ850628","PQ850613","PQ850626","PQ850614",
    "PQ850610","PQ850593","PQ850617","PQ850607","PQ850623","PQ850606",
    "PQ850608","PQ850618","PQ850605","PQ850631","PQ850620","PQ850602",
    "PQ850595","PQ850601","PQ850599","NC_005859.1","MH751506.1",
    "NC_000866.4","MH550421.1","NC_001604","NC_008720.1","NC_001416",
    "CP099588","NC_005856.1","NC_001895.1"
]

# Output file
output_file = "phage_proteins.csv"

with open(output_file, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow([
        "Locus_tag",
        "Genome_accession",
        "Protein_accession",
        "Protein_name",
        "Protein_sequence",
        "Protein_length"
    ])

    for acc in accessions:
        print(f"Processing {acc}")
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        for feature in record.features:
            if feature.type == "CDS":
                locus_tag = feature.qualifiers.get("locus_tag", ["NA"])[0]
                protein_id = feature.qualifiers.get("protein_id", ["NA"])[0]
                product = feature.qualifiers.get("product", ["NA"])[0]

                if "translation" in feature.qualifiers:
                    protein_seq = feature.qualifiers["translation"][0]
                    protein_len = len(protein_seq)
                else:
                    protein_seq = "NA"
                    protein_len = "NA"

                writer.writerow([
                    locus_tag,
                    record.id,
                    protein_id,
                    product,
                    protein_seq,
                    protein_len
                ])

        time.sleep(0.4)  # be polite to NCBI servers

print("Done.")
