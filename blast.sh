#!/bin/bash
#SBATCH --job-name=blast_ulcerans            # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=16		                        # Number of cores per task 
#SBATCH --mem=80gb			                            # Total memory for job
#SBATCH --time=00-10:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/cgw47706/scratch/log.%j.out
#SBATCH --error=/scratch/cgw47706/scratch/log.%j.err
#SBATCH --mail-user=cgw47706@uga.edu                    # Where to send mail 
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set directory variables
OUTDIR="/scratch/cgw47706/ulcerans"
FNA="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/925/GCF_000013925.1_ASM1392v2/GCF_000013925.1_ASM1392v2_genomic.fna.gz"
FAA="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/925/GCF_000013925.1_ASM1392v2/GCF_000013925.1_ASM1392v2_protein.faa.gz"
H37RVfna="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz"
H37Rvfaa="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_protein.faa.gz"

if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

#load modules
module load BLAST+/2.14.1-gompi-2023a

cd $OUTDIR

#download ref seq GCF_000013925.1 (M. Ulcerans) from NCBI.
curl -s $FNA | gunzip -c > $OUTDIR/mulcerans.fna
curl -s $FAA | gunzip -c > $OUTDIR/mulcerans.faa

#download ref seq 
curl -s $H37RVfna | gunzip -c > $OUTDIR/H37RV.fna
curl -s $ | gunzip -c > $OUTDIR/H37Rv.faa

#check that downloads were successful before moving to blast.
if [[ ! -f $OUTDIR/mulcerans.fna || ! -f $OUTDIR/mulcerans.faa || ! -f $OUTDIR/H37RV.fna || ! -f $OUTDIR/H37Rv.faa ]]; then
    echo "Error: One or more files failed to download." >&2
    exit 1
fi


#make mtb database to compare m ulcerans against. & makes them run in parallel. 
makeblastdb -in $OUTDIR/H37RV.fna -dbtype nucl -out $OUTDIR/mtb_fna_db &
makeblastdb -in $OUTDIR/H37Rv.faa -dbtype prot -out $OUTDIR/mtb_protein_db &

#blast ref nucleotide sequence to mtb database. Look for high scoring pairs (HSPs) and percent identity to evaluate homology.
#Use BLAST viewer or Biopython for result visualization if needed. Outfmt 6 = tabular out put. outfmt 5 is an XML for more detailed records.
blastn -query $OUTDIR/mulcerans.fna -db $OUTDIR/mtb_fna_db -out $OUTDIR/blast_fna_results.out -outfmt 6 -num_threads 8 &

#blast ref protein sequence to mtb database. 
blastp -query $OUTDIR/H37Rv.faa -db $OUTDIR/mtb_protein_db -out $OUTDIR/blast_protein_results.out -outfmt 6 -evalue 1e-5 -num_threads 8 &

#wait until both processes are finished.
wait