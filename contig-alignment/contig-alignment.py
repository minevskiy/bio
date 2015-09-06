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

import csv
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

result_handle = open("sample-input.txt")
subject_handle = open("sample-input.csv", "rU")

SUBJ_LEN = 34000
LEN_LIMIT = 500

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

gd_diagram = GenomeDiagram.Diagram(fragments=1, x=0.01, yt=0.05, yb=0, start=0, end=SUBJ_LEN, tracklines=False)

# Scale track
gd_track_for_scale = GenomeDiagram.Track(scale=True, scale_ticks=True, 
              scale_largetick_interval=10000, scale_smalltick_interval=1000,
              scale_largetick_labels=True, scale_smalltick_labels=True,
              scale_fontangle=315)
gd_diagram.add_track(gd_track_for_scale, 1)

# Query tracks
index = 1
csvreader = csv.reader(result_handle, delimiter='\t')
for row in csvreader:
    if len(row) == 0:
        continue
    row = prepare_values(row)
    if row[LENGTH] < LEN_LIMIT:
        continue

    if row[SSTART] > row[SEND]:
        # swap values
        row[SSTART],row[SEND] = row[SEND],row[SSTART]
        strand = -1
    else:
        strand = 1

    if row[QSTART] > row[QEND]:
        # swap values
        row[QSTART],row[END] = row[QEND],row[QSTART]
        
    # query feature coordinates
    q_feat_start = row[SSTART] - row[QSTART]
    q_feat_end = q_feat_start + row[QLEN]

    # trim query coordinates to fit inside the subject. Render an arrow on the end which does not fit.
    arrowhead_length_qstart=0
    arrowhead_length_qend=0
    if q_feat_start < 0:
        q_feat_start = 0
        arrowhead_length_qstart=0.5
    if q_feat_end > SUBJ_LEN:
        q_feat_end = SUBJ_END
        arrowhead_length_qend=0.5

    gd_feature_set = GenomeDiagram.FeatureSet()
    gd_track = GenomeDiagram.Track()

    # because a feature can only have one arrow, render two features for the query
    # add/subtract 200 to not cover the arrow ends
    query_feature_start = SeqFeature(FeatureLocation(q_feat_start, q_feat_end-200), strand=-1)
    gd_feature_set.add_feature(query_feature_start, color=colors.lightblue, sigil="BIGARROW", 
                           arrowshaft_height=1.0, arrowhead_length=arrowhead_length_qstart)

    query_feature_end = SeqFeature(FeatureLocation(q_feat_start+200, q_feat_end), strand=1)
    gd_feature_set.add_feature(query_feature_end, color=colors.lightblue, sigil="BIGARROW", 
                           arrowshaft_height=1.0, arrowhead_length=arrowhead_length_qend,
                           label=True, label_size = 10, label_position="end",
                           label_color=colors.white, label_angle=180, label_strand=-1,
                           name=row[QSEQID])
                 
    # render the alignmet feature with an arrow
    match_feature = SeqFeature(FeatureLocation(row[SSTART], row[SEND]),strand=strand)
    gd_feature_set.add_feature(match_feature, color=colors.blue, sigil="BIGARROW",
                           arrowshaft_height=1.0, arrowhead_length=0)
                           
    gd_track.add_set(gd_feature_set)
    gd_diagram.add_track(gd_track, index + 1)

    index = index + 1

result_handle.close()
print("Number of tracks: " + str(index-1))

# Subject track
gd_feature_set = GenomeDiagram.FeatureSet()
csvreader = csv.reader(subject_handle, delimiter=',')
# skip the header row
csvreader.next()
for row in csvreader:
    if len(row) == 0:
        continue
    start = int(row[1])
    end = int(row[2])
    strand = 1
    if row[3] == "-":
        strand = -1
    name = row[4]
    
    subject_feature = SeqFeature(FeatureLocation(start, end), strand=strand)
    gd_feature_set.add_feature(subject_feature, color=colors.green, sigil="BIGARROW", 
                           arrowshaft_height=0.5, arrowhead_length=1,
                           label=True, label_size = 10, label_position="start",
                           label_strand=1, label_angle=35, name=name)
        
subject_handle.close()
gd_track = GenomeDiagram.Track()
gd_track.add_set(gd_feature_set)
gd_diagram.add_track(gd_track, index + 1)

gd_diagram.draw(format="linear", pagesize='A2')
gd_diagram.write("bio"+str(LEN_LIMIT)+".svg", output="SVG", dpi=300)


