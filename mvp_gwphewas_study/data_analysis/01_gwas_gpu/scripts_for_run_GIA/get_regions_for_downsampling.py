#!/usr/bin/python
import sys, os, glob

regions_file = "/ccs/home/arodriguez/MVP.rel4.regions.txt"
regions = {}
fh = open(regions_file, "r")
for line in fh:
    (name, chr, start, end) = line.rstrip("\n").split("\t")
    regions[name] = [chr, start, end]

downsampling_variants_dir = "/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/downsampling_variant_lists/"
downsampling_files = glob.glob("%s/*downsampling_variants.txt" % downsampling_variants_dir)

for f in downsampling_files:
    fn = open(f, "r")
    outf = "%s.alex.txt" % f
    fo = open(outf, "w")
    for var in fn:
        (chr, pos, ref, allele) = var.rstrip("\n").split(":")
        for region in regions:
            if (int(regions[region][0]) == int(chr) and int(regions[region][1]) <= int(pos) and int(regions[region][2]) >= int(pos)):
                fo.write("%s\t%s\n" % (var.rstrip("\n"), region))
                break

