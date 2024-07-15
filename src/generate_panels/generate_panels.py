"""
Author: Charles Van Goethem
Affiliation: CHU Montpellier
Email: c (-) vangoethem (@) chu (-) montpellier (.) fr
Date: 2024-07-03
"""

import argparse
import json
import os
import re
import sys

from cyvcf2 import VCF
from liftover import get_lifter
import pandas as pd
import pybedtools

def common_member(a, b):
    """
    Checks if there are any common elements between two sets.

    Args:
        a (set): The first set.
        b (set): The second set.

    Returns:
        bool: True if there are common elements, False otherwise.
    """
    return len(set(a).intersection(set(b))) > 0


def parse_gtf(gtf_file, source=["BestRefSeq"], chromosomes=False, chunksize=5000):
    """
    Parses a GTF file and returns a pandas DataFrame containing information about genes, transcripts, and exons.

    Args:
        gtf_file (str): Path to the GTF file.
        source (list, optional): Source(s) to filter by (default: ["BestRefSeq"]).
        chromosomes (dict, optional): Optional dictionary to map chromosome names (default: None, keep original names).
        chunksize (int, optional): Number of lines to read at a time (default: 5000).

    Returns:
        pandas.DataFrame: A DataFrame containing parsed GTF data.
    """
    
    # Open the GTF file to count total lines
    with open(gtf_file, "r") as f:
        number_of_lines = sum(1 for _ in f)
    
    chunks = 0
    gene_info = []

    for chunk in pd.read_csv(gtf_file, sep='\t', names=['chrom', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attributes'], comment='#', chunksize=chunksize):
        # Convert start/end to integers, filter by sources and remapped chromosomes
        chunk["source"] = chunk["source"].apply(lambda row: row.split("%2C"))
        chunk = chunk[chunk['source'].apply(lambda x: common_member(x, source))]

        chunk = chunk[chunk['chrom'].isin(chromosomes)]
        if chromosomes:
            # chunk = chunk[chunk['chrom'].isin([])]
            chunk.loc[:, "chrom"] = chunk["chrom"].map(chromosomes)

        chunk.loc[:, "start"] = chunk["start"].astype(int)
        chunk.loc[:, "end"] = chunk["end"].astype(int)

        for index, row in chunk.iterrows():
            # Progress bar
            w = int((index+1)*20/number_of_lines)
            p = (index+1)*100/number_of_lines
            sys.stdout.write('\r')
            sys.stdout.write("[%-20s] %d%% Treat GTF" % ('='*w, p))
            sys.stdout.flush()

            # split attributes
            attributes = row['attributes'].split('"; ')
            attr = dict()
            pattern = r'(.*?) "(.*+)'
            for attribute in attributes:
                if attribute:
                    match = re.findall(pattern, attribute)
                    if match[0][0] in attr:
                        if type(attr[match[0][0]]) is list:
                            attr[match[0][0]].append(match[0][1])
                        else:
                            attr[match[0][0]] = [attr[match[0][0]], match[0][1]]
                    else:
                        attr[match[0][0]] = match[0][1]
                
            infos = row[:8]
            
            # Create an ID
            name = f'{attr["gene"]}-{row["type"]}'
            if "transcript_id" in attr and attr["transcript_id"]:
                name = f'{name}-{attr["transcript_id"]}'
            if "exon_number" in attr and attr["exon_number"]:
                name = f'{name}-{attr["exon_number"]}'
            infos["name"] = name
            
            gene_info.append(dict(infos, **attr))
        chunks+=chunksize
    sys.stdout.write('\r')
    sys.stdout.write("[%-20s] %d%% Treat GTF" % ('='*20, 100))
    sys.stdout.flush()
    sys.stdout.write("\n")

    return pd.DataFrame(gene_info)

def split_df(df, column, output):
    """
    Splits a pandas DataFrame into multiple DataFrames based on a unique column value.

    Args:
        df (pandas.DataFrame): The input DataFrame.
        column (str): The column name to split by.
        output (str): The base filename for the split DataFrames.

    Returns:
        dict: A dictionary where keys are the unique column values and values are the corresponding DataFrames.
    """
    elements = df[column].unique()
    df_splitted = dict()
    for element in elements:
        df_splitted[element] = df.loc[df[column] == element]
        df_splitted[element].to_csv(f'{output}.{element}.tsv', index=False, sep="\t")
    
    return df_splitted

def get_bed(df, types, padding, chr):
    """
    Creates a BED file from a pandas DataFrame containing gene annotations.

    Args:
        df (pandas.DataFrame): The input DataFrame with gene annotations.
        types (list): A list of feature types to include (e.g., ["CDS", "start_codon", "stop_codon"]).
        padding (int): The padding to apply around each feature.
        chr (dict): A dictionary mapping chromosome names to their UCSC style names (e.g., {"chr1": "1"}).

    Returns:
        pandas.DataFrame: A DataFrame containing BED records for the specified feature types.
    """
    subset_df = df[df['type'].isin(types)]
    
    subset_df.loc[:, "start"] = subset_df["start"] - padding
    subset_df.loc[:, "end"] = subset_df["end"] + padding

    if chr:
        subset_df.loc[:, "chrom"] = subset_df["chrom"].map(chr)

    bed_df = subset_df[["chrom", "start", "end", "name"]]
    
    return bed_df


def treat_variant(v, clnsig_matches=False, df=pd.DataFrame({"chrom": [],"start":[],"end":[],"name":[]}),max_len=50):
    """
    Processes a VCF variant record and updates a pandas DataFrame with relevant information.

    Args:
        v (cyvcf2.Variant): The VCF variant record to process.
        clnsig_matches (bool, optional): Whether to filter by CLNSIG terms (default: False).
        df (pandas.DataFrame, optional): The DataFrame to update with variant information (default: empty DataFrame).
        max_len (int, optional): The maximum length of variants to consider (default: 50).

    Returns:
        pandas.DataFrame: The updated DataFrame with variant information.
    """
    if not "CLNSIG" in dict(v.INFO):
        return df
    if not "CLNREVSTAT" in dict(v.INFO):
        return df
    if v.INFO["CLNREVSTAT"] == "no_classification_provided":
        return df
    try:
        length = max(len(v.REF), len(v.ALT[0]))
    except IndexError:
        length = 0
    if length > max_len:
        return df

    clnsig = v.INFO["CLNSIG"].replace("/",",").replace("|",",").split(",")
    if (set(clnsig_matches) & set(clnsig)):
        v.INFO["CLNSIG"]
        df2 = pd.Series({"chrom": v.CHROM, "start": v.POS,"end":v.POS + length,"name": f'Clinvar:{v.ID} ({length})'})
        df = pd.concat([df, df2.to_frame().T], ignore_index=True)
    return df



def clinvar_to_bed(clinvar, region_df=False, clnsig_matches=False, clnrevstat_matches=False, max_len=50):
    """
    Converts a DataFrame containing ClinVar variants to a BED file format.

    Args:
        df (pandas.DataFrame): The DataFrame containing ClinVar variant information.
        output_file (str): The filename for the output BED file.
    """
    df = pd.DataFrame({"chrom": [],"start":[],"end":[],"name":[]})

    vcf = VCF(clinvar)
    if isinstance(region_df, pd.DataFrame):
        for index, row in region_df.iterrows():
            for v in vcf(f'{row["chrom"]}:{row["start"]}-{row["end"]}'):
                df=treat_variant(v, clnsig_matches, df, max_len=max_len)
    else:
        for v in vcf():
            df=treat_variant(v, clnsig_matches, df)

    return df

def lovd_to_bed(filename, max_len):
    """
    Converts a LOVD variant file to a BED file format.

    Args:
        lovd_file (str): Path to the LOVD variant file.
        output_file (str): Path to the output BED file.
    """

    data = []
    converter = get_lifter('hg19', 'hg38', one_based=True)

    df= pd.read_csv(filename, sep='\t')

    df.loc[:, "position_g_start"] = df["position_g_start"].astype(int)
    df.loc[:, "position_g_end"] = df["position_g_end"].astype(int)

    
    for idx, row in df.iterrows():

        chrom = row["chromosome"]
        start = row["position_g_start"]
        end = row["position_g_end"]
        conv_start = converter[chrom][start][0][1]
        conv_end = converter[chrom][end][0][1] + 1
        name = row["URL"][16:]
        if (conv_end - conv_start) > max_len:
            continue

        data.append({
            "chrom":chrom,
            "start":conv_start,
            "end": conv_end,
            "name": name
        })

    return pd.DataFrame(data)


def main(args):
    """
    Main function to execute the script.

    Args:
        args: Command-line arguments parsed by argparse.
    """
    # Traitement des arguments
    if args.verbose:
        print("Mode verbose activé")

    df = pd.read_table(args.chromosomes, sep='\t')
    if args.without_chr:
        df['UCSC style name'] = df['UCSC style name'].apply(lambda x: str(x).replace("chr", ''))
    df=df[df['UCSC style name'].isin(args.usable_chromosomes)]
    chromosomes = dict(df[["RefSeq seq accession", "UCSC style name"]].to_dict('split')['data'])
    
    genes = []
    if args.gene:
        genes = args.gene
    if args.genes_list:
        genes.append([line.rstrip('\n') for line in open(args.genes_list, 'r')])
    genes = sorted(set(sum(genes, [])))

    gtf_df = parse_gtf(args.gtf, args.source, chromosomes)
    missing_genes = set(genes).difference(set(gtf_df["gene"]))
    keep_genes = set(genes).difference(missing_genes)

    filtered_df = gtf_df[(gtf_df['gene'].isin(genes)) & (gtf_df['type'].isin(args.type))]
    df_selected = get_bed(filtered_df, args.type, args.padding, chromosomes)
    df_selected.to_csv(f'{args.output}.selected.bed', index=False, sep="\t", header=None)
    
    missing_type = keep_genes.difference(set(filtered_df["gene"]))
    bed_selected = pybedtools.BedTool.from_dataframe(df_selected)

    filtered_gene_df = gtf_df[(gtf_df['gene'].isin(genes)) & (gtf_df['type'].isin(["gene"]))]
    filtered_gene_df = get_bed(filtered_gene_df, ["gene"], args.padding, chromosomes)
    filtered_gene_df.to_csv(f'{args.output}.gene.bed', index=False, sep="\t", header=None)
    
    if args.clinvar:
        df_clinvar = clinvar_to_bed(args.clinvar, region_df=gtf_df[(gtf_df['gene'].isin(genes)) & gtf_df['type'].isin(args.clinvar_type)], clnsig_matches=args.clinvar_clnsig, max_len=args.max_len_variants)
        df_clinvar.to_csv(f'{args.output}.clinvar.bed', index=False, sep="\t", header=None)

    
    bed_clinvar = pybedtools.BedTool.from_dataframe(df_clinvar)
    bed_clinvar_sub = bed_clinvar.subtract(bed_selected).to_dataframe()

    
    if args.lovd:
        df_lovd = lovd_to_bed(args.lovd, args.max_len_variants)
        df_lovd.to_csv(f'{args.output}.lovd.bed', index=False, sep="\t", header=None)
        bed_lovd = pybedtools.BedTool.from_dataframe(df_lovd)
        bed_lovd_sub = bed_lovd.subtract(bed_selected).to_dataframe()

    df = pd.concat([df_selected, bed_clinvar_sub, bed_lovd_sub], ignore_index=True)
    df.to_csv(f'{args.output}.merge.bed', index=False, sep="\t", header=None)

    if missing_genes or missing_type:
        print(f'Missing genes : {len(missing_genes | missing_type)}')
        for i in missing_genes:
            print(f'    - {i} (no genes found)')
        for i in missing_type:
            print(f'    - {i} (type {args.type} not found for this gene ; clinvar variants with {args.clinvar_clnsig} were selected)')

    return
    
def script():
    # Argument configuration
    parser = argparse.ArgumentParser(description="My Generic Python Script", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose mode")
    parser.add_argument("-o", "--output", help="Specify an output file for gene information", default=os.path.dirname(os.path.realpath(__file__)))
    parser.add_argument("--gtf", help="Specify the gtf file to use", required=True)
    parser.add_argument("-s", "--source", help='Specify the source(s) to use ("BestRefSeq", "Curated", "Gnomon", "RefSeq", "cmsearch", "tRNAscan-SE")', nargs="+", default=["BestRefSeq"])
    parser.add_argument("-t", "--type", help='Specify the sequence type to analyze ("CDS", "exon", "gene", "start_codon", "stop_codon", "transcript")', nargs="+", default=["CDS", "start_codon", "stop_codon"])

    parser.add_argument("-g", "--gene", help="Specify the gene(s)", action='append', nargs='*', required=False)
    parser.add_argument("--genes-list", help="Specify the gene(s)", required=False)

    parser.add_argument("--chromosomes", help="Specify the file with chromosome correspondences", required=False)
    parser.add_argument("--usable-chromosomes", help='Specify the chromosomes you want to use', nargs="+", default=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", 'X', 'Y'])
    parser.add_argument("-w", "--without-chr", help="Specify the chromosome format", action="store_false", default=True)

    parser.add_argument("-c", "--clinvar", help="Specify the ClinVar vcf file", required=False)
    parser.add_argument( "--clinvar-type", help='Specify the sequence type to check on ClinVar ("CDS", "exon", "gene", "start_codon", "stop_codon", "transcript")', nargs="+", default=["gene"])
    parser.add_argument("--clinvar-clnsig", help='Specify the significance levels for ClinVar ("Pathogenic", "Likely_pathogenic", "Pathogenic_low_penetrance", "Likely_pathogenic_low_penetrance", "Established_risk_allele", "Likely_risk_allele")', nargs="+", default=["Pathogenic", "Likely_pathogenic", "Pathogenic_low_penetrance", "Likely_pathogenic_low_penetrance", "Established_risk_allele", "Likely_risk_allele"])
    parser.add_argument("--lovd", help='Specify LOVD file', required=False)
    parser.add_argument("--max-len-variants", help='Specify the max length for clinvar variants', default=150)

    parser.add_argument("-p", "--padding", help='Specify the padding to apply', default=50)
    # Argument parsing
    args = parser.parse_args()

    if not (args.gene or args.genes_list):
        parser.error('No genes provided, add --gene (-g) or --genes-list')

    # Execution of the main function
    main(args)

if __name__ == "__main__":
    script()