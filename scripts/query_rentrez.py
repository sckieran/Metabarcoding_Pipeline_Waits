# This is a sample Python script.

from Bio import Entrez
from Bio import SeqIO
import collections
import sys
import os
import pandas

prefix = sys.argv[1]
os.chdir(sys.argv[2])
genelist = open(sys.argv[3], "r")
genus_search = sys.argv[6]
retmax=sys.argv[5]
if not sys.argv[7]:
    print(
        "WARNING: not providing a valid NCBI API key will slow this step down. Please register for a free account with NCBI and provide a key and the associated email.")
else:
    Entrez.api_key = sys.argv[7]
if not sys.argv[8]:
    print(
        "WARNING: not providing an email address will cause this program to throw warnings, as NCBI requires an email address to use its e-Utilities.")
else:
    Entrez.email = sys.argv[8]
col1 = genelist.readline()
list_of_genes = col1.split()
print(list_of_genes)
gldf = pandas.read_csv(sys.argv[3], sep='\t')
for gene in list_of_genes:
    taxalist = open(sys.argv[4], "r")
    taxalist2 = open(sys.argv[4], "r")
    term_list = gldf[gene]
    term_list = term_list.dropna()
    gene_terms = ' OR '.join(term_list)
    outname = f"{prefix}_{gene}_taxalist_results.txt"
    out1 = open(outname, "w")
    headline = f"taxname\tavail_seq?\tnum_avail\ttaxid\n"
    out1.write(headline)
    out1.close()
    out = open(outname, "a")
    genus_yes = list()
    genus_all = list()
    x = 1
    taxlen = len(taxalist2.readlines())
    for line in taxalist:
        taxname = line.strip()
        taxiddict = {accession: id for key, value in []}
        genus = taxname.split()[0]
        species = taxname.split()[1]
        genus_all.append(genus)
        term_search = f"{taxname}[ORGN] AND ({gene_terms})"
        results = Entrez.esearch(db="nuccore", term=term_search, retmax=20)
        result = Entrez.read(results)
        print("Doing taxon", x, "of", taxlen, "species is", genus, species, "NCBI returns", len(result['IdList']),
              "result(s) for", gene)
        x += 1
        if (int(result['Count'])) <= 0:
            ln2 = f"{taxname}\tno\t0\tNA\n"
            out.write(ln2)
        else:
            genus_yes.append(genus)
            avail_seq = "yes"
            num_avail = int(result['Count'])
            ids = result['IdList']
            handle = Entrez.esummary(db="nuccore", id=ids, rettype="gb")
            records = Entrez.read(handle)
            #print(records)
            for record in records:
                accession = record['Id']
                taxids = int(record['TaxId'])
                taxiddict.update({accession: taxids})
            #print(len(list(set(list(taxiddict.values())))))
            if len(list(set(list(taxiddict.values())))) == 1:
                tid = int(records[0]['TaxId'])
                ln2 = f"{taxname}\t{avail_seq}\t{num_avail}\t{tid}\n"
                out.write(ln2)
                fname = f"{prefix}_{gene}_{genus}_{species}_seqs.fasta"
                fasta = open(fname, "w")
                fres = Entrez.efetch(db="nuccore", id=ids, rettype="fasta")
                for seq_record in SeqIO.parse(fres, "fasta"):
                    name1 = str(seq_record.id)
                    headstring = f"{name1} {genus} {species} taxid={tid}"
                    fasta_format_string = f">{headstring}\n%s\n" % seq_record.seq
                    fasta.write(fasta_format_string)
                fasta.close()
            else:
                # taxid_list = list(set(val for dic in taxiddict for val in taxiddict.values()))
                print("more than one taxid found in hits, removing species that are not matches")
                # taxid_dict = collections.Counter(taxiddict)
                # true_taxid = taxid_dict.most_common()
                new_ids = list()
                tid_term = f"{taxname}[SCIN]"
                h4 = Entrez.esearch(db="taxonomy", term=tid_term)
                r4 = Entrez.read(h4)
                true_taxid = int(r4['IdList'][0])
                #print("taxid_num is",true_taxid)
                for id in taxiddict.keys():
                    id2=taxiddict[id]
                    print(id2)
                    if id2 == true_taxid:
                        new_ids.append(id)
                    else:
                        h6 = Entrez.esummary(db="taxonomy", id=id2)
                        r6 = Entrez.read(h6)
                        taxgenus2 = str(r6[0]['Genus'])
                        taxspecies2 = str(r6[0]['Species'])
                        taxname2=f"{taxgenus2} {taxspecies2}"
                        #print(r6[0]['Rank'])
                        if r6[0]['Rank'] == "subspecies":
                            print("taxa's rank is subspecies. Checking if it is subspecies of correct species.")
                            #print(taxname2, "=", taxname)
                            if taxname2 == taxname:
                                new_ids.append(id)
                                print("all taxa are a subspecies of species in taxlist, including in reference database.")
                ln2 = f"{taxname}\t{avail_seq}\t{num_avail}\t{true_taxid}\n"
                out.write(ln2)
                fname = f"{prefix}_{gene}_{genus}_{species}_seqs.fasta"
                fasta = open(fname, "w")
                fres = Entrez.efetch(db="nuccore", id=new_ids, rettype="fasta")
                for seq_record in SeqIO.parse(fres, "fasta"):
                    name1 = str(seq_record.id)
                    headstring = f"{name1} {genus} {species} taxid={true_taxid}"
                    fasta_format_string = f">{headstring}\n%s\n" % seq_record.seq
                    fasta.write(fasta_format_string)
                fasta.close()
    genus_yes2 = collections.Counter(genus_yes)
    genus_all2 = collections.Counter(genus_all)
    missing_genus = list()
    print(genus_yes2)
    print(genus_all2)
    for taxon in genus_all2:
        if genus_yes2[taxon] == 0:
            missing_genus.append(taxon)
    if len(missing_genus) >= 1:
        print("There were one or more genera in your taxalist not found in NCBI. Missing genera are:", missing_genus)
        if genus_search == "TRUE":
            print(
                "Because genus_search is set to TRUE, we will search for any species in this genus at this gene and include them in your reference database.")
            genus_newyes = list()
            for line in missing_genus:
                taxname = line.strip()
                taxiddict = {accession: id for key, value in []}
                term_search = f"{taxname}[ORGN] AND ({gene_terms})"
                results = Entrez.esearch(db="nuccore", term=term_search, retmax=20)
                result = Entrez.read(results)
                print("genus is", genus, " sp. NCBI returns", len(result['IdList']), "result(s) for", gene)
                if (int(result['Count'])) <= 0:
                    ln2 = f"{taxname}\tno\t0\tNA\n"
                    out.write(ln2)
                    genus_newyes.append(genus)
                else:
                    avail_seq = "yes"
                    num_avail = int(result['Count'])
                    ids = result['IdList']
                    handle = Entrez.esummary(db="nuccore", id=ids, rettype="gb")
                    records = Entrez.read(handle)
                    # print(records)
                    tid = int(records[0]['TaxId'])
                    ln2 = f"{taxname}\t{avail_seq}\t{num_avail}\t{tid}\n"
                    out.write(ln2)
                    fname = f"{prefix}_{gene}_{taxname}_sp._seqs.fasta"
                    fasta = open(fname, "w")
                    fres = Entrez.efetch(db="nuccore", id=ids, rettype="fasta")
                    for seq_record in SeqIO.parse(fres, "fasta"):
                        name1 = str(seq_record.id)
                        headstring = f"{name1} {genus} {species} taxid={tid}"
                        fasta_format_string = f">{headstring}\n%s\n" % seq_record.seq
                        fasta.write(fasta_format_string)
        genus_newyes2 = collections.Counter(genus_newyes)
        genus_stillmissing = list()
        for taxon in missing_genus:
            if genus_newyes2[taxon] == 0:
                genus_stillmissing.append(taxon)
        print("there are", len(genus_stillmissing), "genera that have no sequences for any species at", gene,
              ". They are:", genus_stillmissing)
    taxalist.close()
    taxalist2.close()
    out.close()
