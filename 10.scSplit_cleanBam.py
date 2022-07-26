#logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
#qsub -q short.q -b y -wd $logDir -j y -R y -pe smp 12 -V python /share/ScratchGeneral/nenbar/projects/CITE/scripts/10.scSplit_cleanBam.py        
import pysam
import csv

cluster_dict = {}
with open('/share/ScratchGeneral/nenbar/projects/CITE/raw_files/3921_3963_4066/HTO_fastqs/3921_3963_4066_whitelist.txt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    #skip header
    for row in csv_reader:
        cluster_dict[row[0]+"-1"] = 1

clusters = set(x for x in cluster_dict.values())


fin = pysam.AlignmentFile("/share/ScratchGeneral/nenbar/projects/CITE/raw_files/3921_3963_4066/outs/possorted_genome_bam.bam", "rb")

# open one bam file
fouts_dict = {}
fout = pysam.AlignmentFile("/share/ScratchGeneral/nenbar/projects/CITE/sc_split/3921_3963_4066_whitelist.bam", "wb", template = fin)
fouts_dict[1] = fout

for read in fin:
    tags = read.tags
    CB_list = [ x for x in tags if x[0] == "CB"]
    if CB_list:
        cell_barcode = CB_list[0][1]
    # the bam files may contain reads not in the final clustered barcodes
    # will be None if the barcode is not in the clusters.csv file
    else: 
        continue
    cluster_id = cluster_dict.get(cell_barcode)
    if cluster_id:
        fouts_dict[cluster_id].write(read)

## do not forget to close the files
fin.close()
for fout in fouts_dict.values():
    fout.close()
