"""
Author: Charles Van Goethem
Affiliation: CHU Montpellier
Email: c (-) vangoethem (@) chu (-) montpellier (.) fr
Date: 2024-07-03
"""

import argparse
import os
import re
import sys

from cyvcf2 import VCF
from liftover import get_lifter
import pandas as pd
import pybedtools
from extendLib.extArgparse import ListAppend


def common_member(first_list, second_list):
    """
    Checks if there are any common elements between two sets.

    Args:
        a (set): The first set.
        b (set): The second set.

    Returns:
        bool: True if there are common elements, False otherwise.
    """
    return len(set(first_list).intersection(set(second_list))) > 0


def parse_gtf(gtf_file, source=None, chromosomes=False, chunksize=5000):
    """
    Parses a GTF file and returns a pandas DataFrame containing information
    about genes, transcripts, and exons.

    Args:
        gtf_file (str): Path to the GTF file.
        source (list, optional): Source(s) to filter by
        (default: None).
        chromosomes (dict, optional): Optional dictionary to map chromosome
        names (default: None, keep original names).
        chunksize (int, optional): Number of lines to read at a time
        (default: 5000).

    Returns:
        pandas.DataFrame: A DataFrame containing parsed GTF data.
    """
    if source is None:
        source = []
        source.append("BestRefSeq")

    # Open the GTF file to count total lines
    with open(gtf_file, "r") as gtf_open:
        number_of_lines = sum(1 for _ in gtf_open)

    chunks = 0
    gene_info = []
    c_names = [
        'chrom', 'source', 'type', 'start',
        'end', 'score', 'strand', 'frame',
        'attributes']

    for chunk in pd.read_csv(gtf_file, sep='\t', names=c_names, comment='#', chunksize=chunksize):
        # Convert start/end to integers, filter by sources and remapped chromosomes
        chunk["source"] = chunk["source"].apply(
            lambda row: row.split("%2C"))
        chunk = chunk[chunk['source'].apply(
            lambda sources_gtf: common_member(sources_gtf, source))]

        if chromosomes:
            chunk = chunk[chunk['chrom'].isin(chromosomes)]
            # chunk = chunk[chunk['chrom'].isin([])]
            chunk.loc[:, "chrom"] = chunk["chrom"].map(chromosomes)

        chunk.loc[:, "start"] = chunk["start"].astype(int)
        chunk.loc[:, "end"] = chunk["end"].astype(int)

        for index, row in chunk.iterrows():
            # Progress bar
            steps = int((index+1)*20/number_of_lines)
            percent = (index+1)*100/number_of_lines
            sys.stdout.write('\r')
            sys.stdout.write("[%-20s] %d%% Treat GTF" % ('='*steps, percent))
            sys.stdout.flush()

            # split attributes
            attributes = row['attributes'].split('"; ')
            attr = dict()
            pattern = r'(.*?) "(.*+)'
            for attribute in attributes:
                if attribute:
                    match = re.findall(pattern, attribute)
                    if match[0][0] in attr:
                        if isinstance(attr[match[0][0]], list):
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
        chunks += chunksize
    sys.stdout.write('\r')
    sys.stdout.write("[%-20s] %d%% Treat GTF" % ('='*20, 100))
    sys.stdout.flush()
    sys.stdout.write("\n")

    return pd.DataFrame(gene_info)

def get_bed(dataframe, types, padding, chr):
    """
    Creates a BED file from a pandas DataFrame containing gene annotations.

    Args:
        df (pandas.DataFrame): The input DataFrame with gene annotations.
        types (list): A list of feature types to include
        (e.g., ["CDS", "start_codon", "stop_codon"]).
        padding (int): The padding to apply around each feature.
        chr (dict): A dictionary mapping chromosome names to their UCSC style
        names (e.g., {"chr1": "1"}).

    Returns:
        pandas.DataFrame: A DataFrame containing BED records for the specified
        feature types.
    """
    subset_df = dataframe[dataframe['type'].isin(types)]

    subset_df.loc[:, "start"] = subset_df["start"] - int(padding)
    subset_df.loc[:, "end"] = subset_df["end"] + int(padding)

    bed_df = subset_df[["chrom", "start", "end", "name"]]

    return bed_df


def treat_variant(variant, clnsig_matches=False, dataframe=pd.DataFrame({"chrom": [], "start":[], "end":[], "name":[]}), max_len=50):
    """
    Processes a VCF variant record and updates a pandas DataFrame with
    relevant information.

    Args:
        variant (cyvcf2.Variant): The VCF variant record to process.
        clnsig_matches (bool, optional): Whether to filter by CLNSIG terms
        (default: False).
        df (pandas.DataFrame, optional): The DataFrame to update with variant
        information (default: empty DataFrame).
        max_len (int, optional): The maximum length of variants to consider
        (default: 50).

    Returns:
        pandas.DataFrame: The updated DataFrame with variant information.
    """
    if "CLNSIG" not in dict(variant.INFO):
        return dataframe
    if "CLNREVSTAT" not in dict(variant.INFO):
        return dataframe
    if variant.INFO["CLNREVSTAT"] == "no_classification_provided":
        return dataframe
    try:
        length = max(len(variant.REF), len(variant.ALT[0]))
    except IndexError:
        length = 0
    if length > max_len:
        return dataframe

    clnsig = variant.INFO["CLNSIG"].replace("/", ",").replace("|", ",").split(",")
    if set(clnsig_matches) & set(clnsig):
        df2 = pd.Series({
            "chrom": variant.CHROM,
            "start": variant.POS,
            "end":variant.POS + length,
            "name": f'Clinvar:{variant.ID} ({length})'})
        dataframe = pd.concat([dataframe, df2.to_frame().T], ignore_index=True)
    return dataframe

def clinvar_to_bed(clinvar, region_df=False, clnsig_matches=False, max_len=50):
    """
    Converts a DataFrame containing ClinVar variants to a BED file format.

    Args:
        df (pandas.DataFrame): The DataFrame containing ClinVar variant
        information.
        output_file (str): The filename for the output BED file.
    """
    bed_df = pd.DataFrame({"chrom": [], "start":[], "end":[], "name":[]})

    vcf = VCF(clinvar)
    if isinstance(region_df, pd.DataFrame):
        for index, row in region_df.iterrows():
            for variant in vcf(f'{row["chrom"]}:{row["start"]}-{row["end"]}'):
                bed_df = treat_variant(variant, clnsig_matches, bed_df, max_len=max_len)
    else:
        for variant in vcf():
            bed_df = treat_variant(variant, clnsig_matches, bed_df)

    return bed_df

def lovd_to_bed(filename, max_len=50):
    """
    Converts a LOVD variant file to a BED file format.

    Args:
        lovd_file (str): Path to the LOVD variant file.
        output_file (str): Path to the output BED file.
    """

    data = []
    converter = get_lifter('hg19', 'hg38', one_based=True)

    lovd_df = pd.read_csv(filename, sep='\t')

    lovd_df.loc[:, "position_g_start"] = lovd_df["position_g_start"].astype(int)
    lovd_df.loc[:, "position_g_end"] = lovd_df["position_g_end"].astype(int)


    for idx, row in lovd_df.iterrows():
        chrom = row["chromosome"]
        start = row["position_g_start"]
        end = row["position_g_end"]
        name = row["URL"][16:]
        conv_start = converter[chrom][start]
        conv_end = converter[chrom][end]
        if len(conv_start) != 1 or len(conv_end) != 1:
            print(f"Multiple position for : '{name}' => please check manually")
            continue
        if (conv_end[0][1] - conv_start[0][1] + 1) > max_len:
            continue

        data.append({
            "chrom":chrom,
            "start":conv_start[0][1],
            "end": conv_end[0][1] + 1,
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

    if args.chromosomes:
        chr_df = pd.read_table(args.chromosomes, sep='\t')
        if args.without_chr:
            chr_df['UCSC style name'] = chr_df['UCSC style name'].apply(
                lambda x: str(x).replace("chr", ''))
        chr_df = chr_df[chr_df['UCSC style name'].isin(args.usable_chromosomes)]
        chr_names = ["RefSeq seq accession", "UCSC style name"]
        chromosomes = dict(chr_df[chr_names].to_dict('split')['data'])
    else:
        chromosomes = False

    genes = []
    if args.gene:
        genes = args.gene
    if args.genes_list:
        genes_list_file = [line.rstrip('\n') for line in open(args.genes_list, 'r')]
        genes += genes_list_file
    genes = sorted(set(genes))

    gtf_df = parse_gtf(args.gtf, args.source, chromosomes)
    
    missing_genes = set(genes).difference(set(gtf_df["gene"]))
    keep_genes = set(genes).difference(missing_genes)
    if genes:
        filtered_df = gtf_df[(gtf_df['gene'].isin(genes)) & (gtf_df['type'].isin(args.type))]
    else:
        filtered_df = gtf_df[(gtf_df['type'].isin(args.type))]
    df_selected = get_bed(filtered_df, args.type, args.padding, chromosomes)
    df_selected.to_csv(
        f'{args.output}.selected.bed',
        index=False, sep="\t", header=None)

    missing_type = keep_genes.difference(set(filtered_df["gene"]))
    bed_selected = pybedtools.BedTool.from_dataframe(df_selected)

    if genes:
        filtered_gene_df = gtf_df[(gtf_df['gene'].isin(genes)) & (gtf_df['type'].isin(["gene"]))]
    else:
        filtered_gene_df = gtf_df[(gtf_df['type'].isin(["gene"]))]
    filtered_gene_df = get_bed(
        filtered_gene_df, ["gene"],
        args.padding, chromosomes)
    filtered_gene_df.to_csv(
        f'{args.output}.gene.bed',
        index=False, sep="\t", header=None)
    
    merged_df = df_selected

    if args.clinvar:
        df_clinvar = clinvar_to_bed(
            args.clinvar,
            region_df=filtered_gene_df,
            clnsig_matches=args.clinvar_clnsig,
            max_len=args.max_len_variants)
        df_clinvar.to_csv(
            f'{args.output}.clinvar.bed',
            index=False, sep="\t", header=None)

        bed_clinvar = pybedtools.BedTool.from_dataframe(df_clinvar)
        bed_clinvar_sub = bed_clinvar.subtract(bed_selected).to_dataframe()
        merged_df = pd.concat([df_selected, bed_clinvar_sub], ignore_index=True)

    if args.lovd:
        df_lovd = lovd_to_bed(args.lovd, args.max_len_variants)
        df_lovd.to_csv(f'{args.output}.lovd.bed', index=False, sep="\t", header=None)
        bed_lovd = pybedtools.BedTool.from_dataframe(df_lovd)
        bed_lovd_sub = bed_lovd.subtract(bed_selected).to_dataframe()

        merged_df.to_csv(f'{args.output}.merge.bed', index=False, sep="\t", header=None)
        merged_df = pd.concat([df_selected, bed_lovd_sub], ignore_index=True)

    if missing_genes or missing_type:
        print(f'Missing genes : {len(missing_genes | missing_type)}')
        for i in sorted(missing_genes):
            print(f'    - {i} (gene not found)')
        for i in sorted(missing_type):
            print(f'    - {i} (type {args.type} not found for this gene ; '
                  'if files were provided variants from LOVD and/or ClinVar will be kept for this gene)')

def script():
    """
    Argument for script launching
    """
    # Argument configuration
    parser = argparse.ArgumentParser(
        description="My Generic Python Script",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose mode")
    parser.add_argument(
        "-o",
        "--output",
        help="Specify an output file for gene information",
        default=os.path.dirname(os.path.realpath(__file__)))
    parser.add_argument(
        "--gtf",
        help="Specify the gtf file to use", required=True)
    parser.add_argument(
        "-s",
        "--source",
        help='Specify the source(s) to use ("BestRefSeq", "Curated", '
        '"Gnomon", "RefSeq", "cmsearch", "tRNAscan-SE")',
        action=ListAppend,
        nargs="+",
        default=["BestRefSeq"])
    parser.add_argument(
        "-t",
        "--type",
        help='Specify the sequence type to analyze ("CDS", "exon", "gene", '
        '"start_codon", "stop_codon", "transcript")',
        action=ListAppend,
        nargs="+",
        default=["CDS", "start_codon", "stop_codon"])

    parser.add_argument(
        "-g",
        "--gene",
        help="Specify the gene(s)",
        action=ListAppend,
        nargs='*',
        required=False)
    parser.add_argument(
        "--genes-list",
        help="Specify the gene(s)",
        required=False)

    parser.add_argument(
        "--chromosomes",
        help="Specify the file with chromosome correspondences",
        required=False)
    parser.add_argument(
        "--usable-chromosomes",
        help='Specify the chromosomes you want to use',
        nargs="+",
        default=[
            "1", "2", "3", "4", "5", "6", "7", "8", "9",
            "10", "11", "12", "13", "14", "15", "16", "17",
            "18", "19", "20", "21", "22", 'X', 'Y'
            ])
    parser.add_argument(
        "-w",
        "--without-chr",
        help="Specify the chromosome format",
        action="store_false",
        default=True)

    parser.add_argument(
        "-c", "--clinvar",
        help="Specify the ClinVar vcf file",
        required=False)
    parser.add_argument(
        "--clinvar-type",
        help='Specify the sequence type to check on ClinVar ("CDS", "exon", '
             '"gene", "start_codon", "stop_codon", "transcript")',
        nargs="+",
        default=["gene"])
    parser.add_argument(
        "--clinvar-clnsig",
        help='Specify the significance levels for ClinVar ("Pathogenic", '
             '"Likely_pathogenic", "Pathogenic_low_penetrance", '
             '"Likely_pathogenic_low_penetrance", "Established_risk_allele", '
             '"Likely_risk_allele")',
        nargs="+",
        default=[
            "Pathogenic", "Likely_pathogenic", "Pathogenic_low_penetrance",
            "Likely_pathogenic_low_penetrance", "Established_risk_allele",
            "Likely_risk_allele"
            ])
    parser.add_argument(
        "--lovd",
        help='Specify LOVD file',
        required=False)
    parser.add_argument(
        "--max-len-variants",
        help='Specify the max length for clinvar variants',
        default=150)

    parser.add_argument(
        "-p",
        "--padding",
        help='Specify the padding to apply',
        default=50)

    # Argument parsing
    args = parser.parse_args()
    print(args.type)

    # Execution of the main function
    main(args)

if __name__ == "__main__":
    script()
