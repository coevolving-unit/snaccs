# Create barcode list from 10X reference
echo "Creating RNA barcode list..."
if [ -f "/data/3M-february-2018.txt" ]; then
    cat /data/3M-february-2018.txt | \
    sed 's/-1/-RNA/g' > RNA.filtered.barcode.list.tsv
    echo "Barcode list created: RNA.filtered.barcode.list.tsv"
else
    echo "WARNING: 3M-february-2018.txt not found. Please download 10X barcode whitelist."
fi

# Create short sample names
if [ -f "bams.txt" ]; then
    cut -d'_' -f2 bams.txt > bams_short.txt
    echo "Created short sample names: bams_short.txt"
else
    echo "WARNING: bams.txt not found. Please create this file with full sample names."
fi

# Download reference genome
echo "Downloading X chromosome reference..."
FASTA_URL="http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz"

if [ ! -f "chromosome.X.fa" ]; then
    wget "$FASTA_URL"
    gunzip -c Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz > chromosome.X.fa
    rm Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
    echo "X chromosome reference downloaded: chromosome.X.fa"
else
    echo "X chromosome reference already exists: chromosome.X.fa"
fi

echo "Setup completed!"
