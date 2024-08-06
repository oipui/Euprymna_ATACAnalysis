#ich würde ein dictionary für start und eins für stop machen und ein dictionary of lists mit allen synt_ids per chromosome
from collections import defaultdict #das brauchst du für das dictionary of lists
import re
import os
start_dict={}
stop_dict={}
id_dict=defaultdict(list) #list of chromosomes and synteny ID
list_of_chroms=[] #make a list of all the chromosomes so you can loop throug it later
syntenyID_dict={}

with open("/Users/pui/Documents/Lab_Vienna/ATAC/Cluster/EupSc_inCeph.clus") as ceph:
    for line in ceph.readlines():
        columns = line.rstrip("\n").split("\t")
        #print(columns)
        region = range(int(columns[1]),int(columns[2]))
        ID = columns[5]
        chrom=columns[0]
        start=columns[1]
        stop=columns[2]
        start_dict[ID]=int(start)
        stop_dict[ID]=int(stop)
        if ID not in id_dict[chrom]:
            id_dict[chrom].append(ID) #füge jede ID zu einer liste des jeweiligen chromosoms hinzu
        else:
            print("THIS ID IS IN HERE TWICE; WHY? %s ID %s id_dict[chrom] %s " %(ID, id_dict[chrom] ,line))
        if chrom not in list_of_chroms:
            list_of_chroms.append(chrom)
        #syntenyID_dict[region]= ID
        #print(syntenyID_dict)
#print(id_dict)
#print(start_dict)

#jetz hast du: ein dictionary mit den synteny IDs für jedes chromosome, ein dictionary mit synt_ID als key und dem start und ein dictionary mit der synt_id als key und dem stop

#again make dictionaries for the peak file


#you can also make a longname dict, so its easier to rename the chromosomes back to the long names

directory="/Users/pui/Documents/Lab_Vienna/ATAC/Peakfiles/Annotations/Genrich_peakfiles/Bedfiles_forHomerMotifsearch/"
for filename in os.listdir(directory):
    if filename.endswith("longnames.bed"):
        peak_start_dict = {}
        peak_stop_dict = {}
        peak_id_dict = defaultdict(list)
        longname_dict = {}
        with open(
                "/Users/pui/Documents/Lab_Vienna/ATAC/Peakfiles/Annotations/Genrich_peakfiles/Bedfiles_forHomerMotifsearch/"+filename) as bed:
            # newfile = open("/Users/pui/Documents/Lab_Vienna/ATAC/Peakfiles/Annotations/Genrich_peakfiles/Bedfiles_forHomerMotifsearch/Stage20_peaks_bothRep_all.longnames.bed2","w")
            for line in bed.readlines():
                columns = line.rstrip("\n").split("\t")
                # print(columns)
                start = int(columns[1])
                stop = int(columns[2])
                chrom_long = (columns[0])
                peak_id = columns[3]
                chrom = re.split("__", columns[0])[0]  # rename long chromosome names to chr
                chrom = chrom.replace("Lachesis_group", "chr")
                peak_start_dict[peak_id] = start
                peak_stop_dict[peak_id] = stop
                peak_id_dict[chrom].append(peak_id)
                if chrom not in longname_dict:
                    longname_dict[chrom] = chrom_long
        #newfile = open(
         #   "/Users/pui/Documents/Lab_Vienna/ATAC/Peakfiles/Annotations/Genrich_peakfiles/Bedfiles_forHomerMotifsearch/" + filename + "2",
          #  "w")
        for chrom in list_of_chroms:
            # print(chrom)
            for ID in id_dict[chrom]:
                # print(ID)
                # print(id_dict[chrom])
                for Peak_ID in peak_id_dict[chrom]:
                    if peak_start_dict[Peak_ID] and peak_stop_dict[Peak_ID] in range(start_dict[ID], stop_dict[ID]):
                        # print(longname_dict[chrom], peak_start_dict[Peak_ID], peak_stop_dict[Peak_ID], Peak_ID, start_dict[ID], stop_dict[ID],ID)
                        strand = "."
                        empty = ""
                        print(ID)
                        print(Peak_ID)
                       # newfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (longname_dict[chrom], peak_start_dict[Peak_ID], peak_stop_dict[Peak_ID], Peak_ID, empty, strand, ID))
                    else:
                        print("uh this peak is lost")
        #print(peak_start_dict)
        #print(peak_stop_dict)
        #print(peak_id_dict)

