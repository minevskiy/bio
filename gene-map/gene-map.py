# Copyright, 2014-2015, Ivan Minevskiy

'''
http://home.cc.umanitoba.ca/~psgendb/birchhomedir/BIRCHDEV/doc/NCBI/blast_formatter.txt

qlen means Query sequence length
slen means Subject sequence length
qstart means Start of alignment in query
qend means End of alignment in query
sstart means Start of alignment in subject
send means End of alignment in subject
evalue means Expect value
bitscore means Bit score
score means Raw score
length means Alignment length
pident means Percentage of identical matches
'''

import csv, sys
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

if len(sys.argv) == 0:
    print 'Usage: gene-map.py <subject file>'
    sys.exit(2)

subject_filename = sys.argv[1]

SUBJ_LEN = 45000
LEN_LIMIT = 1000

# BLAST column indices
QSEQID = 0
SSEQID = 1
PIDENT = 2
SSTART = 3
SEND = 4
EVALUE = 5
BITSCORE = 6
SCORE = 7
QSTART = 8
QEND = 9
QLEN = 10
LENGTH = 11

def prepare_values(row):
    row[PIDENT] = float(row[PIDENT])
    row[SSTART] = int(row[SSTART])
    row[SEND] = int(row[SEND])
    row[QSTART] = int(row[QSTART])
    row[QEND] = int(row[QEND])
    row[QLEN] = int(row[QLEN])
    row[EVALUE] = float(row[EVALUE])
    row[BITSCORE] = float(row[BITSCORE])
    row[SCORE] = int(row[SCORE])
    row[LENGTH] = int(row[LENGTH])
    return row

gd_diagram = GenomeDiagram.Diagram(fragments=1, x=0.1, yt=0.7, yb=0, start=0, end=SUBJ_LEN, tracklines=False)

# Scale track
gd_track_for_scale = GenomeDiagram.Track(scale=True, scale_ticks=True, 
              scale_largetick_interval=10000, scale_smalltick_interval=1000,
              scale_largetick_labels=True, scale_smalltick_labels=True,
              scale_fontangle=315)
gd_diagram.add_track(gd_track_for_scale, 1)

# Subject track
gd_feature_set = GenomeDiagram.FeatureSet()
subject_handle = open(subject_filename, "rU")
csvreader = csv.reader(subject_handle, delimiter='\t')

# the header row
start_index = 0
end_index = 0
strand_index = 0
product_index = 0

row = csvreader.next()
for i in xrange(0, len(row)):
    if row[i] == "start":
        start_index = i
    elif row[i] == "end":
        end_index = i
    elif row[i] == "strand":
        strand_index = i
    elif row[i] == "product":
        product_index = i

for row in csvreader:
    if len(row) < 4 or row[0] == "":
        continue
    start = int(row[start_index])
    end = int(row[end_index])
    strand = 1
    if row[strand_index] == "-":
        strand = -1
    name = row[product_index]
    
    subject_feature = SeqFeature(FeatureLocation(start, end), strand=strand)
    gd_feature_set.add_feature(subject_feature, color=colors.green, sigil="BIGARROW", 
                           arrowshaft_height=0.5, arrowhead_length=1,
                           label=True, label_size = 10, label_position="start",
                           label_strand=1, label_angle=35, name=name)
        
subject_handle.close()
gd_track = GenomeDiagram.Track()
gd_track.add_set(gd_feature_set)
gd_diagram.add_track(gd_track, 1)

gd_diagram.draw(format="linear", pagesize=(120*cm, 15*cm))
gd_diagram.write(subject_filename + ".svg", output="SVG", dpi=600)


