#!/bin/bash
# download https://rdr.ucl.ac.uk/ndownloader/files/28624629 & extract
cd 'Nannini et al_mAb_1864084_PacBio data'
cp $(find -type f -name "CD160*fastq") ../0.Raw/

# Define output directory and input FASTQ files
res="1.Preprocessed"
files=("0.Raw/CD160_P0_CCS.fastq" \
       "0.Raw/CD160_P1_CCS.fastq" \
       "0.Raw/CD160_P2_CCS.fastq" \
       "0.Raw/CD160_P3_CCS.fastq")

# Define primers
primer5="GTCGTCTTTCCAGACGTTAGT"
primer3="CAGGAAACAGCTATGAC"

# Define size filter parameters
minlength=900
maxlength=1200

# Process each file
for fastq in "${files[@]}"; do
  basename=$(basename "$fastq" .fastq)
  
  # Step 1: basseto convert
  ./basseto convert --input "$fastq" --input-format fastq \
  --output "${res}/${basename}.filtered.fastq" --output-format fastq \
  --quality 0.99 --threads 8 
  
  # Step 2: cutadapt remove primers and reverse complement
  cutadapt -a "${primer5}...${primer3}" "${res}/${basename}.filtered.fastq" \
  --revcomp --untrimmed-output "${res}/${basename}.filtered.noadapt.fastq" \
  -o "${res}/${basename}.filtered.adapt.fastq" > "${res}/${basename}.cutadapt.out"
  
  # Step 3: cutadapt filter by size
  cutadapt "${res}/${basename}.filtered.adapt.fastq" \
  -m "${minlength}" -M "${maxlength}" \
  --too-short-output "${res}/${basename}.filtered.adapt.${minlength}.fastq" \
  --too-long-output "${res}/${basename}.filtered.adapt.${maxlength}.fastq" \
  -o "${res}/${basename}.filtered.adapt.sizetrim.fastq" > "${res}/${basename}.cutsize.out"
  
  # Step 4: convert to fasta
  basseto convert -i "${res}/${basename}.filtered.adapt.sizetrim.fastq" \
  -n fastq -o "${res}/${basename}.filtered.adapt.sizetrim.fasta" -m fasta
done