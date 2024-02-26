import csv
import os

filename = "output_file.csv"
f = open(filename, "r")
header = f.readline().strip().split() #Header to ignore
lines = f.readlines() #Get all lines (Besides header)
f.close()

with open("Link_Output.csv", "a", newline='') as o:
    csv_writer = csv.writer(o)
    for line in lines:
        line = line.strip().split()
        link = "https://img.jgi.doe.gov/cgi-bin/mer/main.cgi?section=MetaGeneDetail&page=genePageMainFaa&taxon_oid=" + line[0] + "&data_type=assembled&gene_oid=" + line[2].split(",")[0]
        if line[10] == "alcohol":
            mlist = [line[14], link]
        else:
            mlist = [line[10], link]


def linkgen():
    filename = "output_file.csv"
    f = open(filename, "r")
    header = f.readline().strip().split() #Header to ignore
    lines = f.readlines() #Get all lines (Besides header)
    f.close()

    with open("Link_Output.csv", "a", newline='') as o:
        csv_writer = csv.writer(o)
        csv_writer.writerow(["Gene seq", "Link to JGI page"])
        for line in lines:
            line = line.strip().split()
            link = "https://img.jgi.doe.gov/cgi-bin/mer/main.cgi?section=MetaGeneDetail&page=genePageMainFaa&taxon_oid=" + line[0] + "&data_type=assembled&gene_oid=" + line[2].split(",")[0]
            if line[10] == "alcohol":
                mlist = [line[14], link]
            else:
                mlist = [line[10], link]
            csv_writer.writerow(mlist)
    print("done")
