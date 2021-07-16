# VCF_annotation

## Part I: Import database, libraries and functions defined in "find_effect_variation.R"

## Part II: VCF file cleaning
Expand INFO and FORMAT to columns (Output: "vcf_extend.csv")

## Part III: VCF file annotation
#### 0. Initialize annotated dataframe: vcf_annotation
#### 1. Type of variation and their effect
**Score severity of variation effect**
High risk (3): Nonsense mutation (Stop codon gain), Stop codon lost, Start codon lost, Frameshift
Moderate risk (2): Missense mutation, Inframe insertion, Inframe deletion
Low risk (1): Intergenic/Intragenic variant
No risk (0): Silent mutation
**Annotate variation effect**
Type: snp
- If variation is in cds region:
  - NO: *Intergenic/Intragenic variant*
  - YES -> find codon -> if encoding amino acid is changed:
    - NO: *Silent mutation*
    - YES -> if this is start codon or stop codon:
      - NO: *Missense mutation*
      - YES: *Nonsense mutation* or *Stop codon lost* or *Start codon lost*
Type: insertion or deletion
- If variation is in cds region:
  - NO: *Intergenic/Intragenic variant*
  - YES -> if the insertion/deletion is multiple of 3:
    - NO: *Frameshift*
    - YES: *Inframe insertion* or *Inframe deletion*
Type: mnp (length = 2)
- If variation is in cds region:
  - NO: *Intergenic/Intragenic variant*
  - YES -> if two nucleotides are in the same codon:
    - YES: check codon as defined for snp
    - NO: check two codons individually as defined for snp, then select the most deleterious effect
Type: complex
- If variation is in cds region:
  - NO: *Intergenic/Intragenic variant*
  - YES -> if the length is changed in variant:
    - YES: follow steps defined for insertion/deletion
    - NO: check all related codons, then select the most deleterious effect

#### 2. Depth of sequence coverage at the site of variation
#### 3. Number of reads supporting the variant
#### 4. Percentage of reads supporting the variant versus those supporting reference reads
#### 5. Allele frequency of variant from ExAC API

## Part IV: save annotated VCF files
Output: vcf_annotation.csv
