import re
longname_dict={}
with open("/Users/pui/Documents/Lab_Vienna/ATAC/Peakfiles/new_all_chroms_uniq_forr.gff") as bed:
   # longnames = open("Chromosomes_longnames_dictionary","w")
    for lines in bed.readlines():
        columns = lines.rstrip("\n").split("\t")
        long=columns[0]
        chrom = re.split('__',columns[0])[0]
        if chrom not in longname_dict:
            longname_dict[chrom]= long
    #print(longname_dict['Lachesis_group0'])


import os

directory="/Users/pui/Documents/Lab_Vienna/ATAC/Peakfiles/Annotations/Genrich_peakfiles/Bedfiles_forHomerMotifsearch/"
for filename in os.listdir(directory):
    #print(filename)
    #print(filename)
    bed = open("/Users/pui/Documents/Lab_Vienna/ATAC/Peakfiles/Annotations/Genrich_peakfiles/Bedfiles_forHomerMotifsearch/"+ filename + ".longnames.bed",'w')
    with open("/Users/pui/Documents/Lab_Vienna/ATAC/Peakfiles/Annotations/Genrich_peakfiles/Bedfiles_forHomerMotifsearch/"+ filename) as peaks:
        #newfile = open("test_file","w")
        for lines in peaks.readlines():
            columns = lines.rstrip("\n").split("\t")
            # print(columns)
            shortname = columns[0]
            # print(shortname)
            longname = longname_dict[shortname]
            strand = "."
            # print(longname)
            bed.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (longname,columns[1],columns[2],columns[4],columns[3],strand))
            #print("%s\t%s\t%s\t%s\t%s\t%s\n" % (longname,columns[1],columns[2],columns[4],columns[3],strand))



