# Generate Panel

Here's a script for generating bed files to create templates for capture-targeted high-throughput sequencing.

## Usage

```
Generate a bed panel from a GTF file

options:
  -h, --help            show this help message and exit
  -v, --verbose         Enable verbose mode (default: False)
  -o OUTPUT, --output OUTPUT
                        Specify an output file for gene information (default: /Users/Charles/Documents/codes/generate_panels/src)
  --gtf GTF             Specify the gtf file to use (default: None)

GTF options:
  -s SOURCE [SOURCE ...], --source SOURCE [SOURCE ...]
                        Specify the source(s) to use ("BestRefSeq", "Curated", "Gnomon", "RefSeq", "cmsearch", "tRNAscan-SE") (default: ['BestRefSeq'])
  -t TYPE [TYPE ...], --type TYPE [TYPE ...]
                        Specify the sequence type to analyze ("CDS", "exon", "gene", "start_codon", "stop_codon", "transcript") (default: ['exon'])

Genes options:
  -g [GENE ...], --gene [GENE ...]
                        Specify the gene(s) (default: None)
  --genes-list GENES_LIST
                        Specify the gene(s) (default: None)

Chromosomes options:
  --chromosomes CHROMOSOMES
                        Specify the file with chromosome correspondences (default: None)
  --usable-chromosomes USABLE_CHROMOSOMES [USABLE_CHROMOSOMES ...]
                        Specify the chromosomes you want to use (default: ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y'])
  -w, --without-chr     Specify the chromosome format (default: True)

Variant (ClinVar/LOVD) options:
  -c CLINVAR, --clinvar CLINVAR
                        Specify the ClinVar vcf file (default: None)
  --clinvar-type CLINVAR_TYPE [CLINVAR_TYPE ...]
                        Specify the sequence type to check on ClinVar ("CDS", "exon", "gene", "start_codon", "stop_codon", "transcript") (default: ['gene'])
  --clinvar-clnsig CLINVAR_CLNSIG [CLINVAR_CLNSIG ...]
                        Specify the significance levels for ClinVar ("Pathogenic", "Likely_pathogenic", "Pathogenic_low_penetrance", "Likely_pathogenic_low_penetrance", "Established_risk_allele", "Likely_risk_allele") (default:
                        ['Pathogenic', 'Likely_pathogenic', 'Pathogenic_low_penetrance', 'Likely_pathogenic_low_penetrance', 'Established_risk_allele', 'Likely_risk_allele'])
  --lovd LOVD           Specify LOVD file (default: None)
  --max-len-variants MAX_LEN_VARIANTS
                        Specify the max length for clinvar variants (default: 150)

Other options:
  -p PADDING, --padding PADDING
                        Specify the padding to apply (default: 0)
```

### GTF file

> Recommended :
> download the GTF from NCBI (e.g. GRCh38 : https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)

### Example

- Generate exome
> `generate-panel --gtf /ncbi_dataset/ncbi_dataset/data/GCF_000001405.40/genomic.gtf `

- Generate panel of genes with only coding region with a pdding of 50

> `generate-panel --gtf /ncbi_dataset/ncbi_dataset/data/GCF_000001405.40/genomic.gtf -g TP53 TNF EGFR VEGFA APOE IL6 TGFB1 MTHFR ESR1 AKT1 -t CDS -p 50`