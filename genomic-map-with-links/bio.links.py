# Copyright, 2014-2015, Ivan Minevskiy

import csv, sys
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

if len(sys.argv) <= 1:
    print 'Usage: bio.links.py <mapping file.txt> <func and taxo table file.txt> <link file.txt>'
    sys.exit(2)

mapping_filename = sys.argv[1]
table_filename = sys.argv[2]
if len(sys.argv) == 4:
    links_filename = sys.argv[3]

SUBJ_LEN = 45000
# space between contigs
CONTIG_BUFFER = 1000

# default input table column indices
ID = 0
LEN = 1
START = 2
END = 3
CID = 4
CLEN = 5
DIR = 6

def parse_table_row(row):
    feature = dict()
    feature["id"] = row[ID]
    feature["start"] = int(row[START])
    feature["end"] = int(row[END])
    feature["len"] = int(row[LEN])
    feature["cid"] = row[CID]
    feature["clen"] = int(row[CLEN])
    feature["dir"] = 1
    if row[DIR] == "-":
        feature["dir"] = -1
    return feature

def is_continuation(name, next_name):
    parts1 = name.rpartition("_")
    parts2 = next_name.rpartition("_")
    if parts1[0] == "" or parts2[0] == "":
        return False
    try:
        index1 = int(parts1[2])
        index2 = int(parts2[2])
    except ValueError:
        return False
    # assume parts are listed in the right order, so don't check index values
    if parts1[0] == parts2[0]:
        return True

def next(name):
    parts = name.rpartition("_")
    if parts[0] == "":
        return ""
    index = int(parts[2])
    return parts[0] + "_" + str(index + 1)

def prev(name):
    parts = name.rpartition("_")
    if parts[0] == "":
        return ""
    index = int(parts[2])
    if index == 0:
        return ""
    return parts[0] + "_" + str(index - 1)
    
# ====================================================================
# Read name mapping
# ====================================================================

mapping_handle = open(mapping_filename, "rU")
csvreader = csv.reader(mapping_handle, delimiter='\t')
name_mapping = dict()
track_order = []
track_count = 0
prev_name = ""
for row in csvreader:
    if len(row) < 3 or row[0] == "":
        continue
    name_mapping[row[0]] = row[1]
    track_order.append(row[0])
    if not is_continuation(prev_name, row[1]):
        track_count = track_count + 1
    prev_name = row[1]
mapping_handle.close()

# ====================================================================
# Read func and taxo table. Stored in a map {contig.id : {subject.id : {subject data} } }
# ====================================================================

table_handle = open(table_filename, "rU")
csvreader = csv.reader(table_handle, delimiter='\t')
subjects = dict()
lengths = dict()
features = dict()
# skip header row and iterate over other rows
csvreader.next()
for row in csvreader:
    if len(row) < 7 or row[0] == "":
        continue
    feature = parse_table_row(row)
    features[feature["id"]] = feature
    if feature["cid"] not in subjects.keys():
        subjects[feature["cid"]] = dict()
        lengths[feature["cid"]] = feature["clen"]
    subject = subjects[feature["cid"]]
    subject[feature["id"]] = feature
table_handle.close()

# ====================================================================
# Build tracks. Order by name mapping.
# ====================================================================

tracks = dict()
track_list = []
prev_cid = ""
default_color = colors.green
color = default_color
offset = 0
for cid in track_order:
    name = name_mapping[cid]
    is_cont = prev_cid != "" and is_continuation(name_mapping[prev_cid], name)
    if is_cont:
        offset = offset + CONTIG_BUFFER + lengths[prev_cid]
    else:
        offset = 0
    if is_cont and color == default_color:
        color = colors.blue
    else:
        color = default_color
    
    gd_feature_set = GenomeDiagram.FeatureSet()
    show_label = True
    for id, feature in subjects[cid].iteritems():
        # update coordinates 
        feature["start"] = feature["start"] + offset
        feature["end"] = feature["end"] + offset
        # create SeqFeature
        gd_feature = SeqFeature(FeatureLocation(feature["start"], feature["end"]), strand=feature["dir"])
        gd_feature_set.add_feature(gd_feature, color=color, sigil="BIGARROW", 
                               arrowshaft_height=0.5, arrowhead_length=1, name=name,
                               label=show_label, label_size = 8, label_position="end",
                               label_strand=1, label_angle=0)
        show_label = False
    
    # if it a continuation of previous contig, reuse the previous track.
    if not is_cont:
        gd_track = GenomeDiagram.Track(name=cid)
        track_list.append(gd_track)
    # map track by cid regardless of continuation
    tracks[cid] = gd_track
    gd_track.add_set(gd_feature_set)
    prev_cid = cid

# ==========================================================================
# Build diagram, add tracks to it
# ==========================================================================

gd_diagram = GenomeDiagram.Diagram(fragments=1, x=0.01, yt=0.01, yb=0, start=0, end=SUBJ_LEN, 
    tracklines=False, track_size=0.2, fragment_size=1)
# Scale track
gd_track_for_scale = GenomeDiagram.Track(scale=True, scale_ticks=True, 
              scale_largetick_interval=10000, scale_smalltick_interval=1000,
              scale_largetick_labels=True, scale_smalltick_labels=True,
              scale_fontangle=315)
gd_diagram.add_track(gd_track_for_scale, 1)

track_index = 2
# Add tracks in reverse order because they are arranged bottom-up
for gd_track in reversed(track_list):
    gd_diagram.add_track(gd_track, track_index)
    track_index = track_index + 1
  
# ==========================================================================
# Read links file
# ==========================================================================
  
links = dict()
links_handle = open(links_filename, "rU")
csvreader = csv.reader(links_handle, delimiter='\t')
for row in csvreader:
    if len(row) < 2 or row[0] == "":
        continue
    # skip dups
    if row[0] + row[1] in links.keys() or row[1] + row[0] in links.keys():
        continue
    feat1 = features[row[0]]
    feat2 = features[row[1]]
    track1 = tracks[feat1["cid"]]
    track2 = tracks[feat2["cid"]]
    if track1 == track2:
        # print "Could not link " + row[0] + " and " + row[1] + " because they are on the same track"
        continue
    link = {"id": row[0] + row[1], "track1" : track1, "start1" : feat1["start"], "end1": feat1["end"], "fid1": feat1["id"],
                                   "track2" : track2, "start2" : feat2["start"], "end2": feat2["end"], "fid2": feat2["id"]}
    links[link["id"]] = link
'''
    if "NapDC1_B04" in name_mapping[link["feat1"]["cid"]] or "NapDC1_B04" in name_mapping[link["feat2"]["cid"]]:
        print str(link["feat1"]["start"]) + ", " + str(link["feat1"]["end"]) + ", " + str(link["feat2"]["start"]) + ", " + str(link["feat2"]["end"]) + ", " + link["feat1"]["id"] + ", " + link["feat2"]["id"]
''' 
links_handle.close()

# ==========================================================================
# Merge adjacent links
# ==========================================================================

merged_links = []

link_keys = links.keys()
for id, link in links.iteritems():
    if id in merged_links:
        continue
    
    # merge forward
    cur_link = link
    next_link = None
    while True:
        fid1 = cur_link["fid1"]
        fid2 = cur_link["fid2"]
        if next(fid1) + next(fid2) in link_keys:
            next_link = links[next(fid1) + next(fid2)]
            if next_link["id"] not in merged_links:
                link["end1"] = next_link["end1"]
                link["end2"] = next_link["end2"]
                '''
                if "NapTol_Functional_sc_0" in link["fid1"] or "NapTol_Functional_sc_0" in link["fid2"]:
                    print "looking at next: " + next_link["id"]
                    print "extended " + link["feat1"]["id"] + " to [" + str(link["feat1"]["start"]) + ", " + str(link["feat1"]["end"])  + "]"
                    print "extended " + link["feat2"]["id"] + " to [" + str(link["feat2"]["start"]) + ", " + str(link["feat2"]["end"])  + "]"
                    print "consumed " + next_link["id"]
                '''
                cur_link = next_link
                merged_links.append(next_link["id"])
                continue
        break

    # merge across
    cur_link = link
    next_link = None
    while True:
        fid1 = cur_link["fid1"]
        fid2 = cur_link["fid2"]
        if next(fid1) + prev(fid2) in link_keys:
            next_link = links[next(fid1) + prev(fid2)]
            if next_link["id"] not in merged_links:
                link["end1"] = next_link["end1"]
                link["start2"] = next_link["start2"]
                cur_link = next_link
                merged_links.append(next_link["id"])
                continue
        break

# ==========================================================================
# Build links
# ==========================================================================

for id, link in links.iteritems():
    if id in merged_links:
        continue
    feat1 = features[link["fid1"]]
    feat2 = features[link["fid2"]]
    '''
    if "NapDC1_B04" not in name_mapping[feat1["cid"]] and "NapDC1_B04" not in name_mapping[feat2["cid"]]:
        continue
    '''
    color = colors.Color(50, 50, 50, alpha=0)
    gd_link = CrossLink((link["track1"], link["start1"], link["end1"]),
                        (link["track2"], link["start2"], link["end2"]),
                         color, colors.lightgrey, feat1["dir"] != feat2["dir"])
    gd_diagram.cross_track_links.append(gd_link)

print "Number of links: " + str(len(gd_diagram.cross_track_links))
# ==========================================================================
# Write to file
# ==========================================================================

gd_diagram.draw(format="linear", pagesize="A1")
gd_diagram.write(mapping_filename + ".svg", output="SVG", dpi=600)


