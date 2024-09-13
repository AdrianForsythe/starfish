#!/bin/bash
#SBATCH -A naiss2023-22-237
#SBATCH -p devcore
#SBATCH -n 2
#SBATCH -t 0-00:10:00

###################
# variables
###################
species="Pae"

threads=2

BASEDIR=/home/adrianf/systbio-project-folder/Paecilomyces_starfish

MISSING=1
MAXCOPY=5

###################
# environment
###################
cd $BASEDIR || exit

if [[ $USER == "adrianf" ]]; then
  source /sw/apps/conda/latest/rackham/bin/activate
  conda activate /crex/proj/naiss2023-23-128/nobackup/ADRIAN/bin/starfish-env || exit
  STARFISH=/home/adrianf/systbio-project-folder/starfish
  GENOME_DIR=$(find /home/adrianf/systbio-project-folder/Genomes/ -maxdepth 1 -name "$species*")

  module load bioinfo-tools samtools
else
  source "/home/adrian/anaconda3/bin/activate"
  conda activate "/home/adrian/anaconda3/envs/starfish" || exit
  STARFISH=/home/adrian/starfish
  GENOME_DIR=$(find /home/adrian/Genomes/"$species"*/)

fi

export PATH=$PATH:$STARFISH/
export PATH=$PATH:$STARFISH/CNEFinder/

STARFISH_DB=$STARFISH/database

###################
# prepare for starfish
###################
# clean up
rm -f "$GENOME_DIR"/*/*.header.* "$GENOME_DIR"/*/*.starfish_format.*

# make sure dir names don't contain '_'
# find /home/adrianf/systbio-project-folder/Genomes/Paecilomyces/ -type d -name '*_*' -exec sh -c 'mv "$0" "$(echo "$0" | sed "s/_//g")"' {} \;

for file in $(find "$GENOME_DIR"/*/*.fna)
  do
    header=$(echo $file | rev | cut -d "/" -f2 | rev)
    echo "adding "$header" to "$file"..."
    awk '{print $1}' "$file" | awk '{sub(/scaffold_/, "scaffold", $1)} 1' | sed "s/>/>$header_/g" > "${file%.fna}".header.fna
  done

for file in $(find "$GENOME_DIR"/*/*.gff3)
  do
    header=$(echo $file | rev | cut -d "/" -f2 | rev)
    echo "adding "$header" to "$file"..."
    awk '{sub(/scaffold_/, "scaffold", $1)} 1' $file | sed "/^#/! s/^/$header_/" | sed 's/jgi.p|Paevar1|//g; s/jgi.p|Paevar_HGY_1|//g' > ${file%.gff3}.header.gff3
  done

# create .txt files detailing the absolute path to each genome's gff3 and assembly:
# realpath "$GENOME_DIR"/*/*.header.fna | perl -pe 's/^(.+?([^\/]+?).fasta)$/\2\t\1/' > assemblies.txt
find "$GENOME_DIR"/*/*.header.fna > assemblies.txt
# realpath "$GENOME_DIR"/*/*.gff3 | perl -pe 's/^(.+?([^\/]+?).final.gff3)$/\2\t\1/' > gffs.txt
find "$GENOME_DIR"/*/*.header.gff3 > gffs.txt

# add columns containing genome code
cat assemblies.txt | rev | cut -d "/" -f2 | rev > genome-codes.txt

files="assemblies.txt
gffs.txt"

for i in $files
do
  echo "formatting "$i"..."
  paste genome-codes.txt $i > $i.temp
  mv $i.temp $i
done

echo "formatting gffs to starfish format..."
# will crash if there are extra delimeter characters in fasta headers or 1st column of gff
$STARFISH/starfish format -f assemblies.txt -g gffs.txt -s '_'

# concatenate all gff3 files into a single file (a useful shortcut for some analyses):
cat "$GENOME_DIR"/*/*starfish_format.gff3 > "$BASEDIR"/"$species".gff3

# concatenate all assembly files and make a blastn database:
rm -fr blastdb
mkdir blastdb

# make sure that this file now points to formated fasta's and gff's
sed -i 's/.fna/.starfish_format.fna/g' assemblies.txt
sed -i 's/.gff3/.starfish_format.gff3/g' gffs.txt

# concatenate
cat "$GENOME_DIR"/*/*.starfish_format.fna > blastdb/"$species".assemblies.fna
makeblastdb -in blastdb/"$species".assemblies.fna -out blastdb/"$species".assemblies -parse_seqids -dbtype nucl

# calculate %GC content across all genomes (useful for visualizing elements later):
$STARFISH/scripts/seq-gc.sh -Nbw 1000 blastdb/"$species".assemblies.fna > "$species".assemblies.gcContent_w1000.bed
rm -f blastdb/"$species".assemblies.fna

# parse the provided eggnog mapper annotations (NB the format of the output file has changed in more recent emapper versions):
cut -f1,12  $STARFISH/examples/ann/*emapper.annotations | grep -v  '#' | grep -v -P '\t-' | perl -pe 's/\t/\tEMAP\t/' | grep -vP '\tNA' > $STARFISH/examples/ann/gene2emap.txt

# retrieve the narrowest eggnog ortholog group per sequence and convert to mcl format:
cut -f1,10 $STARFISH/examples/ann/*emapper.annotations | grep -v '#' | perl -pe 's/^([^\s]+?)\t([^\|]+).+$/\1\t\2/' > $STARFISH/examples/ann/gene2og.txt

# Important! Note that the above command only works with older emapper annotation output, such as the output provided as part of this tutorial.
# In order to parse the output from more recent versions of emapper, you must use the following command to retrieve the narrowest eggnog ortholog group per sequence:
# cut -f1,5 $STARFISH/examples/ann/*emapper.annotations | grep -v '#' | perl -pe 's/(^.+?)\t.+,([^,]+)$/\1\t\2/' | perl -pe 's/@/\t/' > $STARFISH/examples/ann/gene2og.txt

# note that for this command to run, the second column cannot be empty
# convert to .mcl format:
$STARFISH/scripts/geneOG2mclFormat.pl -i $STARFISH/examples/ann/gene2og.txt -o $STARFISH/examples/ann/

# make sure to keep formatted copies of assembly and annotation files for later
rm -f "$GENOME_DIR"/*/*.header.fna "$GENOME_DIR"/*/*.header.gff3

mkdir $BASEDIR/ann/

# you can increase confidence in region homology by only looking at gene ortholog groups with low copy numbers missing from few genomes
$STARFISH/scripts/filterOG.pl \
		-O $STARFISH/examples/ann/gene2og.mcl \
		-a $MISSING \
		-c $MAXCOPY \
		-o $BASEDIR/ann/