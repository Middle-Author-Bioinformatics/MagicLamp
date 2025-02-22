#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys
import time


# TODO: ADD CYTOCHROME 579 HMM
# TODO: ADD COLUMN WITH ORF STRAND


def main():
    def SUM(ls):
        count = 0
        for i in ls:
            count += float(i)
        return count

    def firstNum(string):
        outputNum = []
        for i in string:
            try:
                int(i)
                outputNum.append(i)
            except ValueError:
                break
        Num = "".join(outputNum)
        return Num

    def Strip(ls):
        outList = []
        for i in ls:
            gene = i.split("|")[0]
            outList.append(gene)
        return outList

    def unique(ls, ls2):
        unqlist = []
        for i in ls:
            if i not in unqlist and i in ls2:
                unqlist.append(i)
        return len(unqlist)

    def Unique(ls):
        unqList = []
        for i in ls:
            if i not in unqList:
                unqList.append(i)
        return unqList

    def Unique2(ls):
        unqList = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in unqList:
                unqList.append(hmm)
        return unqList

    def checkFe(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_reduction", "iron_oxidation"]:
                    count += 1
        return count

    def checkDFE1(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if hmm in ["DFE_0461", "DFE_0462", "DFE_0463", "DFE_0464", "DFE_0465"]:
                    count += 1
        return count

    def checkDFE2(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if hmm in ["DFE_0448", "DFE_0449", "DFE_0450", "DFE_0451"]:
                    count += 1
        return count

    def checkGACE(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if hmm in ["GACE_1843", "GACE_1844", "GACE_1845", "GACE_1846", "GACE_1847"]:
                    count += 1
        return count

    def check1(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_aquisition-siderophore_transport_potential", "iron_aquisition-heme_transport", "iron_aquisition-siderophore_transport"]:
                    count += 1
        return count

    def check1_2(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_aquisition-siderophore_synthesis"]:
                    count += 1
        return count

    def check2(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_aquisition-iron_transport", "iron_aquisition-heme_oxygenase"]:
                    count += 1
        return count

    def check3(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] in ["iron_aquisition-siderophore_synthesis"]:
                    count += 1
        return count

    def checkReg(ls):
        count = 0
        for i in ls:
            hmm = i.split("|")[0]
            if re.findall(r'regulation', geneToCatDict[hmm]):
                count += 1
        return count

    def checkMam(ls):
        count = 0
        uniqueLS = []
        for i in ls:
            hmm = i.split("|")[0]
            if hmm not in uniqueLS:
                uniqueLS.append(hmm)
                if geneToCatDict[hmm] == "magnetosome_formation":
                    count += 1
        return count

    def derep(ls):
        outLS = []
        for i in ls:
            if i not in outLS:
                outLS.append(i)
        return outLS

    def cluster(data, maxgap):
        '''Arrange data into groups where successive elements
           differ by no more than *maxgap*

            #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
            [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

            #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
            [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

        '''
        # data = sorted(data)
        data.sort(key=int)
        groups = [[data[0]]]
        for x in data[1:]:
            if abs(x - groups[-1][-1]) <= maxgap:
                groups[-1].append(x)
            else:
                groups.append([x])
        return groups

    def lastItem(ls):
        x = ''
        for i in ls:
            x = i
        return x

    def RemoveDuplicates(ls):
        empLS = []
        for i in ls:
            if i not in empLS:
                empLS.append(i)
            else:
                pass
        return empLS

    def allButTheLast(iterable, delim):
        x = ''
        length = len(iterable.split(delim))
        for i in range(0, length - 1):
            x += iterable.split(delim)[i]
            x += delim
        return x[0:len(x) - 1]

    def secondToLastItem(ls):
        x = ''
        for i in ls[0:len(ls) - 1]:
            x = i
        return x

    def pull(item, one, two):
        ls = []
        counter = 0
        for i in item:
            if counter == 0:
                if i != one:
                    pass
                else:
                    counter += 1
                    ls.append(i)
            else:
                if i != two:
                    ls.append(i)
                else:
                    ls.append(i)
                    counter = 0
        outstr = "".join(ls)
        return outstr

    def stabilityCounter(int):
        if len(str(int)) == 1:
            string = (str(0) + str(0) + str(0) + str(0) + str(int))
            return (string)
        if len(str(int)) == 2:
            string = (str(0) + str(0) + str(0) + str(int))
            return (string)
        if len(str(int)) == 3:
            string = (str(0) + str(0) + str(int))
            return (string)
        if len(str(int)) == 4:
            string = (str(0) + str(int))
            return (string)

    def replace(stringOrlist, list, item):
        emptyList = []
        for i in stringOrlist:
            if i not in list:
                emptyList.append(i)
            else:
                emptyList.append(item)
        outString = "".join(emptyList)
        return outString

    def remove(stringOrlist, list):
        emptyList = []
        for i in stringOrlist:
            if i not in list:
                emptyList.append(i)
            else:
                pass
        outString = "".join(emptyList)
        return outString

    def remove2(stringOrlist, list):
        emptyList = []
        for i in stringOrlist:
            if i not in list:
                emptyList.append(i)
            else:
                pass
        # outString = "".join(emptyList)
        return emptyList

    def removeLS(stringOrlist, list):
        emptyList = []
        for i in stringOrlist:
            if i not in list:
                emptyList.append(i)
            else:
                pass
        return emptyList

    def fasta(fasta_file):
        seq = ''
        header = ''
        Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in fasta_file:
            i = i.rstrip()
            if re.match(r'^>', i):
                if len(seq) > 0:
                    Dict[header] = seq
                    header = i[1:]
                    seq = ''
                else:
                    header = i[1:]
                    seq = ''
            else:
                seq += i
        Dict[header] = seq
        # print(count)
        return Dict

    def fastaRename(fasta_file):
        counter = 0
        seq = ''
        header = ''
        Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in fasta_file:
            i = i.rstrip()
            if re.match(r'^>', i):
                if len(seq) > 0:
                    Dict[header] = seq
                    header = i[1:]
                    header = header.split(" ")[0]
                    counter += 1
                    header = header + "_" + str(counter)
                    seq = ''
                else:
                    header = i[1:]
                    header = header.split(" ")[0]
                    counter += 1
                    header = header + "_" + str(counter)
                    seq = ''
            else:
                seq += i
        Dict[header] = seq
        # print(count)
        return Dict

    def filter(list, items):
        outLS = []
        for i in list:
            if i not in items:
                outLS.append(i)
        return outLS

    def delim(line):
        ls = []
        string = ''
        for i in line:
            if i != " ":
                string += i
            else:
                ls.append(string)
                string = ''
        ls = filter(ls, [""])
        return ls

    parser = argparse.ArgumentParser(
        prog="FeGenie.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        *******************************************************

        Developed by Arkadiy Garber and Nancy Merino;
        University of Southern California, Earth Sciences
        Please send comments and inquiries to arkadiyg@usc.edu

            )`-.--.  )\.---.     )\.-.    )\.---.   )\  )\  .'(   )\.---.  
            ) ,-._( (   ,-._(  ,' ,-,_)  (   ,-._( (  \, /  \  ) (   ,-._( 
            \ `-._   \  '-,   (  .   __   \  '-,    ) \ (   ) (   \  '-,   
             ) ,_(    ) ,-`    ) '._\ _)   ) ,-`   ( ( \ \  \  )   ) ,-`   
            (  \     (  ``-.  (  ,   (    (  ``-.   `.)/  )  ) \  (  ``-.  
             ).'      )..-.(   )/'._.'     )..-.(      '.(    )/   )..-.(                                                                                    
                                  %(?/////////&//%                                                
              .,,.                   (%((&@@@#/*.                      .,,.        
              .,,.                     @(((/&@@@#///**                  ...        
                                         #&((///////////////*/@                                
                                                             #*@.                             
                                      ()                   * )//*
                                      <^^>             *     (/*   .
                                     .-""-.                  *)
                          .---.    ."-....-"-._     _...---''`/. '
                         ( (`\ \ .'            ``-''    _.-"'`
                          \ \ \ : :.                 .-'
                           `\`.\: `:.             _.'
                           (  .'`.`            _.'
                            ``    `-..______.-'
                                      ):.  (
                                    ."-....-".
                                  .':.        `.
                                  "-..______..-"

        Image design: Nancy Merino (2018);
        ASCII art: https://manytools.org/hacker-tools/convert-images-to-ascii-art/
        https://ascii.co.uk/text
        *******************************************************
        '''))

    parser.add_argument('-bin_dir', type=str, help="directory of bins", default="NA")

    parser.add_argument('-bin_ext', type=str, help="extension for bins (do not include the period)", default="NA")

    parser.add_argument('-d', type=int, help="maximum distance between genes to be considered in a genomic \'cluster\'."
                                             "This number should be an integer and should reflect the maximum number of "
                                             "genes in between putative iron-related genes identified by the HMM database "
                                             "(default=5)", default=5)

    parser.add_argument('-ref', type=str, help="path to a reference protein database, which must be in FASTA format",
                        default="NA")

    parser.add_argument('-out', type=str, help="name output directory (default=fegenie_out)",
                        default="fegenie_out")

    parser.add_argument('-inflation', type=int, help="inflation factor for final gene category counts (default=1000)",
                        default=1000)

    parser.add_argument('-t', type=int, help="number of threads to use for DIAMOND BLAST and HMMSEARCH "
                                             "(default=1, max=16)", default=1)

    parser.add_argument('-bams', type=str, help="a tab-delimited file with two columns: first column has the genome or "
                                                "metagenome file names; second column has the corresponding BAM file "
                                                "(provide full path to the BAM file). Use this option if you have genomes "
                                                "that each have different BAM files associated with them. If you have a set "
                                                "of bins from a single metagenome sample and, thus, have only one BAM file, "
                                                " then use the \'-bam\' option. BAM files are only required if you would like to create "
                                                "a heatmap that summarizes the abundance of a certain gene that is based on "
                                                "read coverage, rather than gene counts.", default="NA")

    parser.add_argument('-which_bams', type=str, help="if you provided a tab-delimited file specifying multiple BAM files for "
                                                "your metagenome assemblies or bins/genomes, FeGenie will, by default, "
                                                "make the heatmap CSV and dotplot based on the average depth across all of BAM files. "
                                                "However, with this argument, you can specify which bam in that file that you "
                                                      "want FeGenie to use for the generation of a heatmap/dotplot. "
                                                      "For example, if only coverage from the first BAM file is desired, "
                                                      "then you can specify \'-which_bams 1\'. "
                                                      "For the third BAM file in the provided tab-delimited file, \'-which_bams 3\ should be specified'", default="average")

    parser.add_argument('-bam', type=str, help="BAM file. This option is only required if you would like to create "
                                                "a heatmap that summarizes the abundance of a certain gene that is based on "
                                                "read coverage, rather than gene counts. If you have more than one BAM file"
                                               "corresponding to different genomes that you are providing, please use the \'-bams\' "
                                               "argument to provide a tab-delimited file that denotes which BAM file (or files) belongs "
                                               "with which genome", default="NA")

    # parser.add_argument('-delim', type=str, help="delimiter that separates contig names from ORF names (provide this flag if you are "
    #                          "providing your own ORFs. Default delimiter for Prodigal-predicted ORFs is \'_\'", default="_")

    parser.add_argument('-contig_names', type=str, help="contig names in your provided FASTA files. Use this option"
                                                        "if you are providing gene calls in amino acid format (don't forget"
                                                        "to add the \'--orfs\' flag)", default="NA")

    parser.add_argument('-cat', type=str, help="comma-separated list of iron gene categories you'd like FeGenie to look for (default = all categories)", default="NA")

    parser.add_argument('--gbk', type=str, help="include this flag if your bins are in Genbank format", const=True,
                        nargs="?")

    parser.add_argument('--orfs', type=str,
                        help="include this flag if you are providing bins as open-reading frames or genes in FASTA amino-acid format",
                        const=True,
                        nargs="?")

    parser.add_argument('--skip', type=str,
                        help="skip the main part of the algorithm (ORF prediction and HMM searching) "
                             "and re-summarize previously produced results (for example, if you want to re-run using "
                             "the --norm flag, or providing a BAM file). All other flags/arguments need to "
                             "be provided (e.g. -bin_dir, -bin_ext, -out, etc.)", const=True, nargs="?")

    parser.add_argument('--meta', type=str,
                        help="include this flag if the provided contigs are from metagenomic/metatranscriptomic assemblies",
                        const=True, nargs="?")

    parser.add_argument('--norm', type=str,
                        help="include this flag if you would like the gene counts for each iron gene category to be normalized to "
                             "the number of predicted ORFs in each genome or metagenome. Without "
                             "normalization, FeGenie will create a heatmap-compatible "
                             "CSV output with raw gene counts. With normalization, FeGenie will create a "
                             "heatmap-compatible with \'normalized gene abundances\'", const=True, nargs="?")

    parser.add_argument('--all_results', type=str,
                        help="report all results, regardless of clustering patterns and operon structure", const=True, nargs="?")

    parser.add_argument('--heme', type=str,
                        help="find all genes with heme-binding motifs (CXXCH), and output them to a separate summary file", const=True, nargs="?")

    parser.add_argument('--hematite', type=str,
                        help="find all genes with hematite-binding motifs, and output them to a separate summary file", const=True, nargs="?")

    parser.add_argument('--makeplots', type=str,
                        help="include this flag if you would like FeGenie to make some figures from your data?. "
                             "To take advantage of this part of the pipeline, you will need to have Rscipt installed. It is a way for R to be called directly from the command line. "
                             "Please be sure to install all the required R packages as instrcuted in the FeGenie Wiki: "
                             "https://github.com/Arkadiy-Garber/FeGenie/wiki/Installation. "
                             "If you see error or warning messages associated with Rscript, you can still expect to "
                             "see the main output (CSV files) from FeGenie.", const=True, nargs="?")

    parser.add_argument('--nohup', type=str, help="include this flag if you are running FeGenie under \'nohup\', and would like to re-write a currently existing directory.", const=True,
                        nargs="?")


    # CHECKING FOR CONDA INSTALL
    os.system("echo ${iron_hmms} > HMMlib.txt")
    os.system("echo ${rscripts} > rscripts.txt")

    file = open("HMMlib.txt")
    HMMdir = ""
    for i in file:
        HMMdir = i.rstrip()

    bits = "/home/ark/MAB/bin/MagicLamp/hmms/iron/HMM-bitcutoffs.txt"
    rscriptDir = "/home/ark/MAB/bin/MagicLamp/rscripts/"

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_known_args()[0]

    # ************** Checking for the required arguments ******************* #
    cwd = os.getcwd()
    print("checking arguments")

    if args.bin_dir != "NA":
        binDir = args.bin_dir + "/"
        binDirLS = os.listdir(args.bin_dir)
        print(".")
    else:
        print("Looks like you did not provide a directory of genomes/bins or assemblies.")
        print("Exiting")
        raise SystemExit

    if args.bin_ext != "NA":
        print(".")
    else:
        print(
            'Looks like you did not provide an extension for your genomes/bins or assemblies, so FeGenie does not know'
            ' which files in the provided directory are fasta files that you would like analyzed.')
        print("Exiting")
        raise SystemExit


    print(".")
    os.system("mkdir -p %s" % args.out)
    os.system("mkdir -p %s/ORF_calls" % args.out)

    if lastItem(args.out) == "/":
        outDirectory = "%s" % args.out[0:len(args.out)-1]
        outDirectoryLS = os.listdir(outDirectory)
    else:
        outDirectory = "%s" % args.out
        outDirectoryLS = os.listdir("%s" % args.out)

    print("All required arguments provided!")
    print("")

    prodigal = 0
    # *************** MAKE NR A DIAMOND DB AND READ THE FILE INTO HASH MEMORY ************************ #
    if args.ref != "NA":
        try:
            testFile = open(args.ref + ".dmnd")

        except FileNotFoundError:
            print("Making diamond database out of provided reference file")
            os.system("diamond makedb --in %s -d %s" % (args.ref, args.ref))

    # *************** CALL ORFS FROM BINS AND READ THE ORFS INTO HASH MEMORY ************************ #
    BinDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    binCounter = 0
    for i in binDirLS:
        if lastItem(i.split(".")) == args.bin_ext and not re.match(r'\.', i):
            cell = i
            binCounter += 1
            if not args.gbk:

                if args.orfs:
                    testFile = open("%s/%s" % (binDir, i), "r")
                    for line in testFile:
                        if re.match(r'>', line):
                            if re.findall(r'\|]', line):
                                print("Looks like one of your fasta files has a header containing the character: \|")
                                print(
                                    "Unfortunately, this is a problem for FeGenie because it uses that character as delimiter to store important information.")
                                print("Please rename your FASTA file headers")
                                raise SystemExit

                else:
                    try:
                        testFile = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                        print("ORFS for %s found. Skipping Prodigal, and going with %s-proteins.faa" % (i, i))
                        for line in testFile:
                            if re.match(r'>', line):
                                if re.findall(r'\|]', line):
                                    print(
                                        "Looks like one of your fasta files has a header containing the character: \|")
                                    print(
                                        "Unfortunately, this is a problem for FeGenie because it uses that character as delimiter to store important information.")
                                    print("Please rename your FASTA file headers")
                                    raise SystemExit

                    except FileNotFoundError:
                        binFile = open("%s/%s" % (binDir, i), "r")
                        for line in binFile:
                            if re.match(r'>', line):
                                if re.findall(r'\|]', line):
                                    print("Looks like one of your fasta files has a header containing the character: \|")
                                    print(
                                        "Unfortunately, this is a problem for FeGenie because it uses that character as delimiter to store important information.")
                                    print("Please rename your FASTA file headers")
                                    raise SystemExit

                        prodigal = 1
                        print("Finding ORFs for " + cell)
                        if args.meta:
                            os.system("prodigal -i %s/%s -a %s/ORF_calls/%s-proteins.faa -o %s/ORF_calls/%s-prodigal.out -p meta -q" % (
                                binDir, i, outDirectory, i, outDirectory, i))
                        else:
                            os.system(
                                "prodigal -i %s/%s -a %s/ORF_calls/%s-proteins.faa -o %s/ORF_calls/%s-prodigal.out -q" % (
                                    binDir, i, outDirectory, i, outDirectory, i))
            else:
                os.system('gtt-genbank-to-AA-seqs -i %s/%s -o %s/%s.faa' % (binDir, i, outDirectory, i))

                faa = open("%s/%s.faa" % (outDirectory, i))
                faa = fasta(faa)

                gbkDict = defaultdict(list)
                counter = 0

                count = 0
                gbk = open("%s/%s" % (binDir, i))
                for gbkline in gbk:
                    ls = gbkline.rstrip()
                    if re.findall(r'/locus_tag', ls):
                        count += 1

                if count > 0:
                    gbk = open("%s/%s" % (binDir, i))
                    for gbkline in gbk:
                        ls = gbkline.rstrip()
                        if re.findall(r'LOCUS', ls):
                            locus = (ls)
                            locus = (locus.split("       ")[1])
                            locus = locus.split(" ")[0]
                        if re.findall(r'gene   ', ls):
                            gene = (ls)
                            gene = (gene.split("            ")[1])
                            start = (gene.split("..")[0])
                            end = (gene.split("..")[1])
                            start = remove(start, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                            end = remove(end, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                            altContigName = (locus + "_" + start + "_" + end)

                        if re.findall(r'/locus_tag', ls):
                            locusTag = (ls)
                            locusTag = (locusTag.split("=")[1])
                            locusTag = remove(locusTag, ["\""])
                            counter += 1

                        if counter > 0:
                            gbkDict[locus].append(locusTag)
                            counter = 0
                else:
                    # print(i)
                    gbk = open("%s/%s" % (binDir, i))
                    for gbkline in gbk:
                        ls = gbkline.rstrip()
                        if re.findall(r'LOCUS', ls):
                            locus = (ls)
                            locus = (locus.split("       ")[1])
                            locus = locus.split(" ")[0]
                        if re.findall(r'gene   ', ls):
                            gene = (ls)
                            gene = (gene.split("            ")[1])
                            start = (gene.split("..")[0])
                            end = (gene.split("..")[1])
                            start = remove(start, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                            end = remove(end, ["c", "o", "m", "p", "l", "e", "m", "e", "n", "t", "(", ")"])
                            altContigName = (locus + "_" + start + "_" + end)
                            counter += 1

                        if re.findall(r'/locus_tag', ls):
                            locusTag = (ls)
                            locusTag = (locusTag.split("=")[1])
                            locusTag = remove(locusTag, ["\""])

                        if counter > 0:
                            gbkDict[locus].append(altContigName)
                            counter = 0

                idxOut = open("%s/ORF_calls/%s-proteins.idx" % (outDirectory, i), "w")
                faaOut = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "w")

                for gbkkey1 in gbkDict.keys():
                    counter = 0
                    for gbkey2 in gbkDict[gbkkey1]:
                        counter += 1
                        if len(faa[gbkey2]) > 0:
                            newOrf = gbkkey1 + "_" + str(counter)
                            idxOut.write(gbkey2 + "," + newOrf + "\n")
                            faaOut.write(">" + newOrf + "\n")
                            faaOut.write(str(faa[gbkey2]) + "\n")

                idxOut.close()
                faaOut.close()

            if args.orfs:
                file = open("%s/%s" % (binDir, i))
            else:
                file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i))
            file = fasta(file)
            for j in file.keys():
                orf = j.split(" ")[0]
                BinDict[cell][orf] = file[j]

    if binCounter == 0:
        print("Did not detect any files in the provided directory (%s) matching the provided filename extension (%s). "
              "Please double-check the filenames and your command" % (args.bin_dir, args.bin_ext))
        raise SystemExit
    # ******************** READ BITSCORE CUT-OFFS INTO HASH MEMORY ****************************** #
    meta = open(bits, "r")
    metaDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in meta:
        ls = i.rstrip().split("\t")
        metaDict[ls[0]] = ls[1]

    # ******************* BEGINNING MAIN ALGORITHM **********************************))))
    HMMdir = "/home/ark/MAB/bin/MagicLamp/hmms/iron"
    if not args.skip:
        if args.cat == "NA":
            catList = []
            HMMdirLS = os.listdir(HMMdir)
            for FeCategory in HMMdirLS:
                if not re.match(r'\.', FeCategory) and FeCategory not in ["HMM-bitcutoffs.txt", "FeGenie-map.txt"]:
                    catList.append(FeCategory)
        else:
            categories = args.cat
            catList = categories.split(",")

        print("starting main pipeline...")
        HMMdirLS = os.listdir(HMMdir)
        for FeCategory in HMMdirLS:
            if not re.match(r'\.', FeCategory) and FeCategory not in ["HMM-bitcutoffs.txt", "FeGenie-map.txt"] and FeCategory in catList:
                print("")
                print("Looking for following iron-related functional category: " + FeCategory)
                hmmDir = "%s/%s/" % (HMMdir, FeCategory)
                hmmDirLS2 = os.listdir("%s/%s" % (HMMdir, FeCategory))

                HMMdict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: "EMPTY")))
                for i in binDirLS:  # ITERATION THROUGH EACH BIN IN A GIVEN DIRECTORY OF BINS
                    if lastItem(i.split(".")) == args.bin_ext:  # FILTERING OUT ANY NON-BIN-RELATED FILES
                        os.system(
                            "mkdir -p " + outDirectory + "/" + i + "-HMM")  # CREATING DIRECTORY, FOR EACH BIN, TO WHICH HMMSEARCH RESULTS WILL BE WRITTEN

                        count = 0
                        for hmm in hmmDirLS2:  # ITERATING THROUGH ALL THE HMM FILES IN THE HMM DIRECTORY
                            count += 1
                            perc = (count / len(hmmDirLS2)) * 100
                            sys.stdout.write("analyzing " + i + ": %d%%   \r" % (perc))
                            sys.stdout.flush()
                            if len(metaDict[hmm.split(".")[0]]) == 0:
                                bit = 0
                            else:
                                bit = metaDict[hmm.split(".")[0]]

                                if args.orfs:
                                    os.system(
                                        "hmmsearch --cpu %d -T %d --tblout %s/%s-HMM/%s.tblout -o %s/%s-HMM/%s.txt %s/%s %s/%s"
                                        % (int(args.t), float(bit), outDirectory, i, hmm, outDirectory, i, hmm, hmmDir, hmm, binDir, i)
                                    )
                                else:
                                    os.system(
                                        "hmmsearch --cpu %d -T %d --tblout %s/%s-HMM/%s.tblout -o %s/%s-HMM/%s.txt %s/%s %s/ORF_calls/%s-proteins.faa"
                                        % (int(args.t), float(bit), outDirectory, i, hmm, outDirectory, i, hmm, hmmDir, hmm, outDirectory, i)
                                    )

                                # REMOVING THE STANDARD OUTPUT FILE
                                os.system(
                                    "rm " + outDirectory + "/" + i + "-HMM/" + hmm + ".txt"
                                )

                                # READING IN THE HMMSEARCH RESULTS (TBLOUT) OUT FILE
                                try:
                                    hmmout = open(outDirectory + "/" + i + "-HMM/" + hmm + ".tblout", "r")
                                except FileNotFoundError:
                                    print("FeGenie cannot find the correct hmmsearch output files. "
                                          "If you provided gene or ORF-call sequences, "
                                          "please be sure to specify this in the command using the \'--orfs\' flag")

                                # COLLECTING SIGNIFICANT HMM HITS IN THE FILE
                                for line in hmmout:
                                    if not re.match(r'#', line):
                                        ls = delim(line)
                                        evalue = float(ls[4])
                                        bit = float(ls[5])
                                        orf = ls[0]
                                        if evalue < float(1E-1):  # FILTERING OUT BACKGROUND NOISE
                                            # LOADING HMM HIT INTO DICTIONARY, BUT ONLY IF THE ORF DID NOT HAVE ANY OTHER HMM HITS

                                            if orf not in HMMdict[i]:
                                                HMMdict[i][orf]["hmm"] = hmm
                                                HMMdict[i][orf]["evalue"] = evalue
                                                HMMdict[i][orf]["bit"] = bit
                                            else:
                                                # COMPARING HITS FROM DIFFERENT HMM FILES TO THE SAME ORF
                                                if bit > HMMdict[i][orf]["bit"]:
                                                    HMMdict[i][orf]["hmm"] = hmm
                                                    HMMdict[i][orf]["evalue"] = evalue
                                                    HMMdict[i][orf]["bit"] = bit

                        print("")

                out = open(outDirectory + "/%s-summary.csv" % (FeCategory), "w")
                out.write("cell" + "," + "ORF" + "," + "HMM" + "," + "evalue" + "," + "bitscore" + "\n")
                for key in HMMdict.keys():
                    for j in HMMdict[key]:
                        out.write(key + "," + j + "," + HMMdict[key][j]["hmm"] + "," +
                                  str(HMMdict[key][j]["evalue"]) + "," + str(HMMdict[key][j]["bit"]) + "\n")

                out.close()
                time.sleep(5)

        print("\n")
        print("Consolidating summary files into one master summary file")
        out = open(outDirectory + "/FinalSummary.csv", "w")
        out.write("category" + "," + "cell" + "," + "orf" + "," + "related_hmm" + "," + "HMM-bitscore" + "\n")

        resultsDir = os.listdir(outDirectory)
        for i in resultsDir:
            if lastItem(i.split("-")) == "summary.csv":
                result = open(outDirectory + "/" + i, "r")
                for j in result:
                    ls = j.rstrip().split(",")
                    cell = ls[0]
                    orf = ls[1]
                    hmm = ls[2]

                    bit = ls[4]

                    if cell != "cell":
                        out.write(i.split("-summary")[0] + "," + cell + "," + orf + "," + hmm + "," + str(bit) + "\n")

        out.close()
        time.sleep(5)

        # ****************************************** DEREPLICATION *********************************************************
        summary = open(outDirectory + "/FinalSummary.csv", "r")
        SummaryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
        for i in summary:
            ls = i.rstrip().split(",")
            if ls[0] != "category" and ls[0] != "FeGenie":
                if len(ls) > 0:
                    category = ls[0]
                    cell = ls[1]
                    orf = ls[2]
                    hmm = ls[3]
                    hmmBit = ls[4]

                    if cell not in SummaryDict.keys():
                        SummaryDict[cell][orf]["hmm"] = hmm
                        SummaryDict[cell][orf]["hmmBit"] = hmmBit
                        SummaryDict[cell][orf]["category"] = category

                    else:
                        if orf not in SummaryDict[cell]:
                            SummaryDict[cell][orf]["hmm"] = hmm
                            SummaryDict[cell][orf]["hmmBit"] = hmmBit
                            SummaryDict[cell][orf]["category"] = category

                        else:
                            if float(hmmBit) > float(SummaryDict[cell][orf]["hmmBit"]):
                                SummaryDict[cell][orf]["hmm"] = hmm
                                SummaryDict[cell][orf]["hmmBit"] = hmmBit
                                SummaryDict[cell][orf]["category"] = category

        # ************************** CLUSTERING OF ORFS BASED ON GENOMIC PROXIMITY *********************************
        if not args.orfs:
            print("Identifying genomic proximities and putative operons")
            CoordDict = defaultdict(lambda: defaultdict(list))
            orfNameDict = defaultdict(lambda: defaultdict(list))
            for i in SummaryDict.keys():
                if i != "category":
                    for j in SummaryDict[i]:
                        contig = allButTheLast(j, "_")
                        numOrf = lastItem(j.split("_"))
                        CoordDict[i][contig].append(int(numOrf))

            counter = 0
            print("Clustering ORFs...")
            print("")
            out = open(outDirectory + "/FinalSummary-dereplicated-clustered.csv", "w")
            for i in CoordDict.keys():
                print(".")
                for j in CoordDict[i]:
                    LS = (CoordDict[i][j])
                    clusters = (cluster(LS, args.d))
                    unknown = [[759,762,763,764,765], [5079,5080,5081]]
                    for k in clusters:
                        if len(RemoveDuplicates(k)) == 1:

                            orf = j + "_" + str(k[0])

                            out.write(SummaryDict[i][orf]["category"] + "," + i + "," + orf + "," + SummaryDict[i][orf]["hmm"] +
                                      "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(counter) + "\n")

                            out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                            counter += 1

                        else:
                            for l in RemoveDuplicates(k):

                                orf = j + "_" + str(l)

                                out.write(SummaryDict[i][orf]["category"] + "," + i + "," + orf + "," + SummaryDict[i][orf][
                                    "hmm"] + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(counter) + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                            counter += 1
            out.close()

        else:

            if args.contig_names != "NA":
                contigNameDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
                contigNames = open(args.contig_names)
                for orfname in contigNames:
                    orfnameLS = orfname.rstrip().split("\t")
                    contigNameDict[orfnameLS[0]]["contig"] = ls[1]
                    contigNameDict[orfnameLS[0]]["position"] = ls[2]

                CoordDict = defaultdict(lambda: defaultdict(list))
                orfNameDict = defaultdict(lambda: defaultdict(list))
                for i in SummaryDict.keys():
                    if i != "category":
                        for j in SummaryDict[i]:
                            # contigLS = contig.split(args.contig_names + args.delim)
                            numOrf = contigNameDict[j]["position"]
                            contig = contigNameDict[j]["contig"]
                            CoordDict[i][contig].append(int(numOrf))
                            orfNameDict[contig + "_" + numOrf] = j
                            CoordDict[i][contig].append(int(numOrf))

                counter = 0
                print("Clustering ORFs...")
                print("")
                out = open(outDirectory + "/FinalSummary-dereplicated-clustered.csv", "w")
                for i in CoordDict.keys():
                    print(".")
                    for j in CoordDict[i]:
                        LS = (CoordDict[i][j])
                        clusters = (cluster(LS, args.d))
                        for k in clusters:
                            if len(RemoveDuplicates(k)) == 1:

                                orf = orfNameDict[j + "_" + str(k[0])]

                                out.write(
                                    SummaryDict[i][orf]["category"] + "," + i + "," + orf + "," + SummaryDict[i][orf][
                                        "hmm"] +
                                    "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(counter) + "\n")

                                out.write(
                                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                                counter += 1

                            else:
                                for l in RemoveDuplicates(k):
                                    orf = orfNameDict[j + "_" + str(l)]

                                    out.write(SummaryDict[i][orf]["category"] + "," + i + "," + orf + "," +
                                              SummaryDict[i][orf][
                                                  "hmm"] + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(
                                        counter) + "\n")

                                out.write(
                                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                                counter += 1
                out.close()

            else:
                counter = 0
                out = open(outDirectory + "/FinalSummary-dereplicated-clustered.csv", "w")
                for i in SummaryDict.keys():
                    if i != "category":
                        for j in SummaryDict[i]:
                            orf = j
                            out.write(SummaryDict[i][orf]["category"] + "," + i + "," + orf + "," + SummaryDict[i][orf][
                                "hmm"] + "," + str(SummaryDict[i][orf]["hmmBit"]) + "," + str(counter) + "\n")
                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                            counter += 1


        time.sleep(5)

        # else:
        #     counter = 0
        #     for i in SummaryDict.keys():
        #         if i != "category":
        #             for j in SummaryDict[i]:
        #
        #                 out.write(SummaryDict[i][j]["category"] + "," + i + "," + j + "," + SummaryDict[i][j][
        #                     "hmm"] + "," + str(SummaryDict[i][j]["hmmBit"]) + "," + str(counter) + "\n")
        #
        #                 out.write(
        #                     "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
        #
        #                 counter += 1

        # ************************** BLAST-BASED METHODS/LOOKING FOR UNMODELED MARKERS ********************************
        thermincola = "%s/iron_reduction/non-aligned/TherJR_SLCs.faa" % HMMdir
        geobacter = "%s/iron_reduction/non-aligned/geobacter_PCCs.faa" % HMMdir

        print("Looking for Thermincola S-layer cytochromes and Geobacter-related porin-cytochrome operons")

        for i in binDirLS:
            if lastItem(i.split(".")) == args.bin_ext:
                if args.orfs:
                    os.system(
                        "makeblastdb -dbtype prot -in %s/%s -out %s/%s -logfile %s/makedbfile.txt" % (
                            binDir, i, binDir, i, binDir))
                    os.system("rm %s/makedbfile.txt" % binDir)

                    os.system(
                        "blastp -query %s -db %s/%s -num_threads %s -outfmt 6 -out %s/%s-thermincola.blast -evalue 1E-10"
                        % (thermincola, binDir, i, args.t, outDirectory, i))

                    os.system(
                        "blastp -query %s -db %s/%s -num_threads %s -outfmt 6 -out %s/%s-geobacter.blast -evalue 1E-10"
                        % (geobacter, binDir, i, args.t, outDirectory, i))
                    os.system("rm %s/%s.phr" % (binDir, i))
                    os.system("rm %s/%s.psq" % (binDir, i))
                    os.system("rm %s/%s.pin" % (binDir, i))

                else:
                    os.system(
                        "makeblastdb -dbtype prot -in %s/ORF_calls/%s-proteins.faa -out %s/ORF_calls/%s-proteins.faa -logfile %s/makedbfile.txt" % (
                            outDirectory, i, outDirectory, i, outDirectory))
                    os.system("rm %s/makedbfile.txt" % outDirectory)

                    os.system(
                        "blastp -query %s -db %s/ORF_calls/%s-proteins.faa -num_threads %s -outfmt 6 -out %s/%s-thermincola.blast -evalue 1E-10"
                        % (thermincola, outDirectory, i, args.t, outDirectory, i))

                    os.system(
                        "blastp -query %s -db %s/ORF_calls/%s-proteins.faa -num_threads %s -outfmt 6 -out %s/%s-geobacter.blast -evalue 1E-10"
                        % (geobacter, outDirectory, i, args.t, outDirectory, i))
                    os.system("rm %s/ORF_calls/%s-proteins.faa.phr" % (outDirectory, i))
                    os.system("rm %s/ORF_calls/%s-proteins.faa.psq" % (outDirectory, i))
                    os.system("rm %s/ORF_calls/%s-proteins.faa.pin" % (outDirectory, i))

        geoDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        geo = open("%s/iron_reduction/non-aligned/geobacter_PCCs.faa" % HMMdir)
        geo = fasta(geo)
        for i in geo.keys():
            id = i.split(" ")[0]
            type = (i.split(" ")[2])
            type = type[1:len(type) - 1]
            geoDict[id]["type"] = type
            geoDict[id]["header"] = i

        thermDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        therm = open("%s/iron_reduction/non-aligned/TherJR_SLCs.faa" % HMMdir)
        therm = fasta(therm)
        for i in therm.keys():
            id = i.split(" ")[0]
            thermDict[id]["header"] = i

        out = open(outDirectory + "/GeoThermin.csv", "w")
        for blastresult in os.listdir(outDirectory):
            blastDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            if re.findall(r'geobacter\.blast', blastresult):
                blast = open(outDirectory + "/" + blastresult, "r")

                for i in blast:
                    ls = i.rstrip().split("\t")
                    if ls[1] not in blastDict.keys():
                        blastDict[ls[1]]["query"] = ls[0]
                        blastDict[ls[1]]["e"] = ls[10]
                    else:
                        if float(ls[10]) < float(blastDict[ls[1]]["e"]):
                            blastDict[ls[1]]["query"] = ls[0]
                            blastDict[ls[1]]["e"] = ls[10]
                        else:
                            pass

                operonDict = defaultdict(list)
                for i in blastDict.keys():
                    orf = (i)
                    contig = allButTheLast(i, "_")
                    orfCall = lastItem(orf.split("_"))
                    id = (blastDict[i]["query"])
                    type = geoDict[id]["type"]
                    operonDict[contig].append(int(orfCall))

                count = 0
                operonDict2 = defaultdict(lambda: defaultdict(list))
                for i in operonDict.keys():
                    contig = (i)
                    clu = cluster(operonDict[i], 2)
                    for j in clu:
                        if len(j) > 1:
                            operon = (j)
                            count += 1
                            for k in operon:
                                orf = str(contig) + "_" + str(k)
                                id = blastDict[orf]["query"]
                                header = geoDict[id]["header"]
                                type = geoDict[id]["type"]
                                operonDict2["operon" + str(count)]["types"].append(type)
                                operonDict2["operon" + str(count)]["orfs"].append(orf)
                                operonDict2["operon" + str(count)]["headers"].append(header)

                for i in operonDict2.keys():
                    if "porin" in operonDict2[i]["types"] and (
                                    "pc" in operonDict2[i]["types"] or "omc" in operonDict2[i]["types"]):
                        genome = blastresult.split("-geobacter.blas")[0]
                        category = "iron_reduction"
                        for j in range(0, len(operonDict2[i]["types"])):
                            orf = (operonDict2[i]["orfs"][j])
                            evalue = blastDict[orf]["e"]
                            header = (operonDict2[i]["headers"][j])
                            out.write(category + "," + genome + "," + orf + "," + replace(header, [","],
                                                                                          ";") + "," + "evalue: " + str(
                                evalue) + "," + str(counter) + "\n")
                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                        counter += 1
                        # print("\n\n")

            if re.findall(r'thermincola\.blast', blastresult):
                blast = open(outDirectory + "/" + blastresult, "r")

                for i in blast:
                    ls = i.rstrip().split("\t")
                    if ls[1] not in blastDict.keys():
                        blastDict[ls[1]]["query"] = ls[0]
                        blastDict[ls[1]]["e"] = ls[10]
                    else:
                        if float(ls[10]) < float(blastDict[ls[1]]["e"]):
                            blastDict[ls[1]]["query"] = ls[0]
                            blastDict[ls[1]]["e"] = ls[10]
                        else:
                            pass

                operonDict = defaultdict(lambda: defaultdict(list))
                for i in blastDict.keys():
                    orf = (i)
                    contig = allButTheLast(i, "_")
                    id = (blastDict[i]["query"])
                    header = thermDict[id]["header"]
                    operonDict[blastresult]["headers"].append(header)
                    operonDict[blastresult]["orfs"].append(orf)

                count = 0
                operonDict2 = defaultdict(lambda: defaultdict(list))
                for i in operonDict.keys():
                    gene1 = "646797728 YP_003639887 TherJR_1122 cytochrome C family protein [Thermincola sp. JR: NC_014152]"
                    gene2 = "646799199 YP_003641333 TherJR_2595 hypothetical protein [Thermincola sp. JR: NC_014152]"
                    gene3 = "646796949 YP_003639120 TherJR_0333 hypothetical protein [Thermincola sp. JR: NC_014152]"
                    if gene1 in operonDict[i]["headers"] and gene2 in operonDict[i]["headers"] and gene3 in operonDict[i][
                        "headers"]:
                        genome = blastresult.split("-thermincola.blas")[0]
                        category = "iron_reduction"
                        for j in range(0, len(operonDict[i]["orfs"])):
                            orf = (operonDict[i]["orfs"][j])
                            header = (operonDict[i]["headers"][j])
                            evalue = blastDict[orf]["e"]
                            out.write(category + "," + genome + "," + orf + "," + replace(header, [","], ";") + "," + "evalue: " + str(
                                evalue) + "," + str(counter) + "\n")
                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                        counter += 1
        out.close()
        time.sleep(5)

        summary = open("%s/FinalSummary-dereplicated-clustered.csv" % outDirectory)
        out = open("%s/FinalSummary-dereplicated-clustered-blast.csv" % outDirectory, "w")
        DeRepDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in summary:
            ls = i.rstrip().split(",")
            unq = ls[1] + "|" + ls[2]
            DeRepDict[unq] = ls[0]
            out.write(i.rstrip() + "\n")

        blastHits = open("%s/GeoThermin.csv" % outDirectory)
        for i in blastHits:
            ls = i.rstrip().split(",")
            unq = ls[1] + "|" + ls[2]
            if unq not in DeRepDict.keys():
                out.write(i.rstrip() + "\n")

        out.close()
        time.sleep(5)
        # ****************************** FILTERING OUT LIKELY FALSE POSITIVES *************************************

        if not args.all_results and not args.orfs:

            print("Pre-processing of final outout file")
            # fleet = ["EetA.hmm", "EetB.hmm", "Ndh2.hmm", "FmnB.hmm", "FmnA.hmm", "DmkA.hmm", "DmkB.hmm", "PplA.hmm"]
            # mam = ["MamA.hmm", "MamB.hmm", "MamE.hmm", "MamK.hmm", "MamP.hmm", "MamM.hmm", "MamP.hmm", "MamQ.hmm", "MamI.hmm",
            #        "MamL.hmm", "MamO.hmm"]
            # foxabc = ["FoxA.hmm", "FoxB.hmm", "FoxC.hmm"]
            # foxeyz = ["FoxE.hmm", "FoxY.hmm", "FoxZ.hmm"]

            clusterDict = defaultdict(lambda: defaultdict(list))
            summary = open("%s/FinalSummary-dereplicated-clustered-blast.csv" % outDirectory, "r")
            for i in summary:
                if not re.match(r'#', i):
                    ls = i.rstrip().split(",")
                    clusterDict[ls[5]]["line"].append(ls)
                    clusterDict[ls[5]]["gene"].append(ls[3])
                    clusterDict[ls[5]]["category"].append(ls[0])

            out = open("%s/FinalSummary-dereplicated-clustered-blast-filtered.csv" % outDirectory, "w")
            for i in (clusterDict.keys()):
                ls = (clusterDict[i]["gene"])
                if "EetA.hmm" in ls or "EetB.hmm" in ls or "Ndh2.hmm" in ls or "FmnB.hmm" in ls or "FmnA.hmm" in ls or "DmkA.hmm" in ls or "DmkB.hmm" in ls or "PplA.hmm" in ls:
                    fleet = ["EetA.hmm", "EetB.hmm", "Ndh2.hmm", "FmnB.hmm", "FmnA.hmm", "DmkA.hmm", "DmkB.hmm", "PplA.hmm"]

                    if unique(ls, fleet) < 5:  # If there are less than 5 FLEET genes in the cluster
                        if len(remove2(ls, fleet)) < 1:  # If FLEET genes are the only ones in the cluster
                            pass
                        else:  # If there are other genes in the cluster that are not FLEET
                            for j in clusterDict[i]["line"]:
                                if j[3] not in fleet:  # avoiding the fleet genes
                                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:  # if there are 5 or more of the FLEET genes present within cluster
                        for j in clusterDict[i]["line"]:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                elif "MamA.hmm" in ls or "MamB.hmm" in ls or "MamE.hmm" in ls or "MamK.hmm" in ls or "MamM.hmm" in ls or "MamO.hmm" \
                        in ls or "MamP.hmm" in ls or "MamQ.hmm" in ls or "MamI.hmm" in ls or "MamL.hmm" in ls:
                    mam = ["MamA.hmm", "MamB.hmm", "MamE.hmm", "MamK.hmm", "MamP.hmm", "MamM.hmm", "MamP.hmm", "MamQ.hmm",
                           "MamI.hmm", "MamL.hmm", "MamO.hmm"]
                    if unique(ls, mam) < 5:
                        if len(remove2(ls, mam)) < 1:
                            pass
                        else:
                            for j in clusterDict[i]["line"]:
                                if j[3] not in mam:
                                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                elif "MtoA.hmm" in ls or "MtrA.hmm" in ls or "MtrC_TIGR03507.hmm" in ls or "MtrB_TIGR03509.hmm" in ls:
                    if "MtoA.hmm" in ls and "MtrB_TIGR03509.hmm" in ls and "MtrC_TIGR03507.hmm" not in ls:
                        for j in clusterDict[i]["line"]:
                            if j[3] in ["MtrB_TIGR03509.hmm", "MtoA.hmm", "CymA.hmm"]:
                                out.write(
                                    "iron_oxidation" + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                            else:
                                out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    elif "MtrA.hmm" in ls and "MtrB_TIGR03509.hmm" in ls:
                        for j in clusterDict[i]["line"]:
                            if j[3] in ["MtrA.hmm", "MtrB_TIGR03509.hmm"]:
                                out.write(
                                    "iron_reduction" + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                            else:
                                out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    elif "MtrC_TIGR03507.hmm" in ls:
                        for j in clusterDict[i]["line"]:
                            if j[3] in ["MtrA.hmm", "MtrB_TIGR03509.hmm"]:
                                out.write(
                                    "iron_reduction" + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                            else:
                                out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")


                    # elif "MtrB_TIGR03509.hmm" not in ls:
                    #     pass

                elif "FoxA.hmm" in ls or "FoxB.hmm" in ls or "FoxC.hmm" in ls:
                    foxabc = ["FoxA.hmm", "FoxB.hmm", "FoxC.hmm"]
                    if unique(ls, foxabc) < 2:
                        if len(remove2(ls, foxabc)) < 1:
                            pass

                        else:
                            for j in clusterDict[i]["line"]:
                                if j[3] not in foxabc:
                                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                    else:

                        for j in clusterDict[i]["line"]:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                elif "FoxE.hmm" in ls or "FoxY.hmm" in ls or "FoxZ.hmm" in ls:
                    foxeyz = ["FoxE.hmm", "FoxY.hmm", "FoxZ.hmm"]
                    if "FoxE.hmm" not in ls:
                        if len(remove2(ls, foxeyz)) < 1:
                            pass

                        else:
                            for j in clusterDict[i]["line"]:
                                if j[3] not in foxeyz:
                                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")
                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                elif "DFE_0448.hmm" in ls or "DFE_0449.hmm" in ls or "DFE_0450.hmm" in ls or "DFE_0451.hmm" in ls:
                    DFE1 = ["DFE_0448.hmm", "DFE_0449.hmm", "DFE_0450.hmm", "DFE_0451.hmm"]

                    if unique(ls, DFE1) < 3:
                        if len(remove2(ls, DFE1)) < 1:
                            pass

                        else:
                            for j in clusterDict[i]["line"]:
                                if j[3] not in DFE1:
                                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                    else:

                        for j in clusterDict[i]["line"]:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                elif "DFE_0461.hmm" in ls or "DFE_0462.hmm" in ls or "DFE_0463.hmm" in ls or "DFE_0464.hmm" in ls or "DFE_0465.hmm" in ls:
                    DFE2 = ["DFE_0461.hmm", "DFE_0462.hmm", "DFE_0463.hmm", "DFE_0464.hmm", "DFE_0465"]

                    if unique(ls, DFE2) < 3:
                        if len(remove2(ls, DFE2)) < 1:
                            pass

                        else:
                            for j in clusterDict[i]["line"]:
                                if j[3] not in DFE2:
                                    out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                    else:

                        for j in clusterDict[i]["line"]:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                elif "Cyc1.hmm" in ls:
                    if "Cyc2_repCluster3.hmm" not in ls and "Cyc2_repCluster2.hmm" not in ls and "Cyc2_repCluster1.hmm" not in ls:
                        pass

                    else:
                        for j in clusterDict[i]["line"]:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                elif "CymA.hmm" in ls:
                    if "MtrB_TIGR03509.hmm" not in ls and "MtrA.hmm" not in ls and "MtoA.hmm" not in ls and "MtrC_TIGR03507.hmm" not in ls:
                        pass

                elif "iron_aquisition-siderophore_synthesis" in clusterDict[i]["category"] or \
                                "iron_aquisition-siderophore_transport_potential" in clusterDict[i]["category"] or \
                                "iron_aquisition-siderophore_transport" in clusterDict[i]["category"] or \
                        "iron_aquisition-iron_transport" in clusterDict[i][
                            "category"] or "iron_aquisition-heme_transport" in clusterDict[i]["category"]:

                    if len(Unique(ls)) > 1:
                        for j in clusterDict[i]["line"]:
                            out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                        out.write(
                            "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                    else:
                        if "FutA1-iron_ABC_transporter_iron-binding-rep.hmm" in ls or "FutA2-iron_ABC_transporter_iron-binding-rep.hmm" in ls \
                                or "FutC-iron_ABC_transporter_ATPase-rep.hmm" in ls or "LbtU-LvtA-PiuA-PirA-RhtA.hmm" in ls or "LbtU-LbtB-legiobactin_receptor.hmm" in ls \
                                or "LbtU_LbtB-legiobactin_receptor_2.hmm" in ls or "IroC-salmochelin_transport-rep.hmm" in ls or "LbtU-LbtB-legiobactin_receptor.hmm" in ls:
                            for j in clusterDict[i]["line"]:
                                out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                            out.write(
                                "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

                        else:
                            pass

                else:
                    linels = (clusterDict[i]["line"])
                    for j in linels:
                        out.write(j[0] + "," + j[1] + "," + j[2] + "," + j[3] + "," + j[4] + "," + j[5] + "\n")

                    out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")

            out.close()
            time.sleep(5)

            # REMOVING FILES
            os.system("rm %s/GeoThermin.csv" % outDirectory)
            os.system("rm %s/*summary*" % outDirectory)
            os.system("rm %s/FinalSummary-dereplicated-clustered-blast.csv" % outDirectory)
            os.system("rm %s/*blast" % outDirectory)
            os.system("rm %s/FinalSummary.csv" % outDirectory)
            os.system("rm %s/FinalSummary-dereplicated-clustered.csv" % outDirectory)
            os.system("mv %s/FinalSummary-dereplicated-clustered-blast-filtered.csv %s/FeGenie-summary.csv" % (outDirectory, outDirectory))

            # OPTIONAL CROSS-VALIDATION AGAINST REFERENCE DATABASE
            if args.ref != "NA":
                print("")
                print("Performing Diamond BLASTx search of putative iron genes against reference database")
                summary = open("%s/FeGenie-summary.csv" % outDirectory, "r")
                out = open("%s/FeGenie-summary.fasta" % outDirectory, "w")
                for i in summary:
                    if not re.match(r'#', i):
                        ls = (i.rstrip().split(","))
                        seq = (BinDict[ls[1]][ls[2]])
                        header = (">" + ls[1] + "|" + ls[2])
                        out.write(header + "\n")
                        out.write(seq + "\n")
                out.close()
                time.sleep(5)

                os.system("diamond blastp --db %s.dmnd --out "
                          "%s/FeGenie-summary.dmndout --max-target-seqs 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle "
                          "--threads %s --query %s/FeGenie-summary.fasta --quiet" % (args.ref, outDirectory, str(args.t), outDirectory))

                dmndblast = open("%s/FeGenie-summary.dmndout" % outDirectory)
                dmndblastDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'No_Hits')))
                for hit in dmndblast:
                    lst = hit.rstrip().split("\t")
                    evalue = lst[10]
                    cell = lst[0].split("|")[0]
                    orf = lst[0].split("|")[1]
                    target = lst[12]
                    target = replace(target, [","], ";")
                    dmndblastDict[cell][orf]["e"] = evalue
                    dmndblastDict[cell][orf]["target"] = target

            summary = open("%s/FeGenie-summary.csv" % outDirectory, "r")
            out = open("%s/FeGenie-summary-blasthits.csv" % outDirectory, "w")
            if args.ref != "NA":
                out.write(
                    "category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "cluster" + "," + "heme_binding_motifs" + "," + "top_blast_hit" + "," + "blast_hit_evalue" + "," + "protein_sequence" + "\n")
            else:
                out.write(
                    "category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "cluster" + "," + "heme_binding_motifs" + "," + "protein_sequence" + "\n")


            # SUMMARIZING CROSS-REFERENCE RESULTS AND COUNTING HEME-BINDING MOTIFS
            print("Counting heme-binding motifs")
            aromatics = ["F", "Y", "W", "H"]
            counter = 1
            for i in summary:
                if not re.match(r'#', i):
                    ls = i.rstrip().split(",")
                    seq = BinDict[ls[1]][ls[2]]

                    hemes = len(re.findall(r'C(..)CH', seq)) + len(re.findall(r'C(...)CH', seq)) \
                            + len(re.findall(r'C(....)CH', seq)) + len(re.findall(r'C(..............)CH', seq)) \
                            + len(re.findall(r'C(...............)CH', seq))

                    if args.ref != "NA":
                        blasthit = dmndblastDict[ls[1]][ls[2]]["target"]
                        e = dmndblastDict[ls[1]][ls[2]]["e"]
                        try:
                            out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + str(ls[4]) + "," + str(
                                metaDict[ls[3].split(".")[0]]) + "," + str(counter) + "," + str(
                                hemes) + "," + blasthit + "," + str(
                                e) + "," + seq + "\n")
                        except TypeError:
                            out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + ls[4] + "," + str(counter) + "," + str(
                                    hemes) + "," + blasthit + "," + str(e) + "," + seq + "\n")

                    else:
                        try:
                            out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + str(ls[4]) + "," + str(
                                metaDict[ls[3].split(".")[0]]) + "," + str(counter) + "," + str(hemes) + "," + seq + "\n")

                        except TypeError:
                            out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + ls[4] + "," + str(counter) + "," + str(
                                    hemes) + "," + seq + "\n")

                else:
                    counter += 1
                    out.write(i)

            out.close()
            time.sleep(5)

        else:
            print("Pre-processing of final outout file")

            summary = open("%s/FinalSummary-dereplicated-clustered-blast.csv" % outDirectory, "r")
            out = open("%s/FinalSummary-dereplicated-clustered-blast-filtered.csv" % outDirectory, "w")
            for i in summary:
                out.write(i.rstrip() + "\n")

            out.close()

            time.sleep(5)
            # REMOVING FILES
            os.system("rm %s/GeoThermin.csv" % outDirectory)
            os.system("rm %s/*summary*" % outDirectory)
            os.system("rm %s/FinalSummary-dereplicated-clustered-blast.csv" % outDirectory)
            os.system("rm %s/*blast" % outDirectory)
            os.system("rm %s/FinalSummary.csv" % outDirectory)
            os.system("rm %s/FinalSummary-dereplicated-clustered.csv" % outDirectory)
            os.system("mv %s/FinalSummary-dereplicated-clustered-blast-filtered.csv %s/FeGenie-summary.csv" % (
            outDirectory, outDirectory))

            # OPTIONAL CROSS-VALIDATION AGAINST REFERENCE DATABASE
            if args.ref != "NA":
                print("")
                print("Performing Diamond BLASTx search of putative iron genes against reference database")
                summary = open("%s/FeGenie-summary.csv" % outDirectory, "r")
                out = open("%s/FeGenie-summary.fasta" % outDirectory, "w")
                for i in summary:
                    if not re.match(r'#', i):
                        ls = (i.rstrip().split(","))
                        seq = (BinDict[ls[1]][ls[2]])
                        header = (">" + ls[1] + "|" + ls[2])
                        out.write(header + "\n")
                        out.write(seq + "\n")
                out.close()
                time.sleep(5)

                os.system("diamond blastp --db %s.dmnd --out "
                          "%s/FeGenie-summary.dmndout --max-target-seqs 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle "
                          "--threads %s --query %s/FeGenie-summary.fasta --quiet" % (
                          args.ref, outDirectory, str(args.t), outDirectory))

                dmndblast = open("%s/FeGenie-summary.dmndout" % outDirectory)
                dmndblastDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'No_Hits')))
                for hit in dmndblast:
                    lst = hit.rstrip().split("\t")
                    evalue = lst[10]
                    cell = lst[0].split("|")[0]
                    orf = lst[0].split("|")[1]
                    target = lst[12]
                    target = replace(target, [","], ";")
                    dmndblastDict[cell][orf]["e"] = evalue
                    dmndblastDict[cell][orf]["target"] = target

            time.sleep(5)
            summary = open("%s/FeGenie-summary.csv" % outDirectory, "r")
            out = open("%s/FeGenie-summary-blasthits.csv" % outDirectory, "w")
            if args.ref != "NA":
                out.write(
                    "category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "cluster" + "," + "heme_c_binding_motifs" + "," + "top_blast_hit" + "," + "blast_hit_evalue" + "," + "protein_sequence" + "\n")
            else:
                out.write(
                    "category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "cluster" + "," + "heme_c_binding_motifs" + "," + "protein_sequence" + "\n")

            # SUMMARIZING CROSS-REFERENCE RESULTS AND COUNTING HEME-BINDING MOTIFS
            print("Counting c-type heme-binding motifs")
            counter = 1
            for i in summary:
                if not re.match(r'#', i):
                    ls = i.rstrip().split(",")
                    seq = BinDict[ls[1]][ls[2]]
                    hemes = len(re.findall(r'C(..)CH', seq)) + len(re.findall(r'C(...)CH', seq)) \
                            + len(re.findall(r'C(....)CH', seq)) + len(re.findall(r'C(..............)CH', seq)) \
                            + len(re.findall(r'C(...............)CH', seq))

                    if args.ref != "NA":
                        blasthit = dmndblastDict[ls[1]][ls[2]]["target"]
                        e = dmndblastDict[ls[1]][ls[2]]["e"]
                        try:
                            out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + str(ls[4]) + "," + str(
                                metaDict[ls[3].split(".")[0]]) + "," + str(counter) + "," + str(
                                hemes) + "," + blasthit + "," + str(
                                e) + "," + seq + "\n")
                        except TypeError:
                            out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + ls[4] + "," + str(
                                counter) + "," + str(
                                hemes) + "," + blasthit + "," + str(e) + "," + seq + "\n")

                    else:
                        try:
                            out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + str(ls[4]) + "," + str(
                                metaDict[ls[3].split(".")[0]]) + "," + str(counter) + "," + str(hemes) + "," + seq + "\n")

                        except TypeError:
                            out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + ls[4] + "," + str(
                                counter) + "," + str(
                                hemes) + "," + seq + "\n")

                else:
                    counter += 1
                    out.write(i)

            out.close()

        time.sleep(5)
        # REMOVING FILES
        if args.ref != "NA":
            os.system("rm %s/FeGenie-summary.dmndout" % outDirectory)
            os.system("rm %s/FeGenie-summary.fasta" % outDirectory)

        os.system("rm %s/FeGenie-summary.csv" % outDirectory)
        os.system("mv %s/FeGenie-summary-blasthits.csv %s/FeGenie-summary.csv" % (outDirectory, outDirectory))

        # FILTERING OUT FALSE POSITIVES FOR SIDEROPHORE GENES

        if not args.all_results:

            print("Final processing of output\n")
            MAP = open(HMMdir + "/FeGenie-map.txt", "r")
            mapDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in MAP:
                ls = i.rstrip().split("\t")
                mapDict[ls[0]] = ls[1]

            out = open(outDirectory + "/FeGenie-summary-fixed.csv", "w")
            geneToCatDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            memoryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
            clusterDict = defaultdict(list)
            infile = open(outDirectory + "/FeGenie-summary.csv")
            for i in infile:
                if not re.match(r'#', i):
                    ls = i.rstrip().split(",")
                    if ls[6] != "cluster":
                        if not re.findall(r'defaultdict', ls[5]):
                            clu = ls[6]
                            cat = ls[0]
                            dataset = ls[1]
                            orf = ls[2]
                            hmm = allButTheLast(ls[3], ".")
                            clusterDict[clu].append(hmm + "|" + dataset + "|" + orf)
                            geneToCatDict[hmm] = cat
                            hmm = allButTheLast(ls[3], ".")
                            memoryDict[dataset][orf]["cat"] = ls[0]
                            memoryDict[dataset][orf]["gene"] = ls[3]
                            memoryDict[dataset][orf]["bit"] = ls[4]
                            memoryDict[dataset][orf]["cutoff"] = ls[5]
                            memoryDict[dataset][orf]["clu"] = clu
                            memoryDict[dataset][orf]["heme"] = ls[7]
                            memoryDict[dataset][orf]["seq"] = ls[8]
                            if args.ref != "NA":
                                memoryDict[dataset][orf]["blastHit"] = ls[8]
                                memoryDict[dataset][orf]["blastEval"] = ls[9]
                                memoryDict[dataset][orf]["seq"] = ls[10]
                        else:
                            cat = ls[0]
                            dataset = ls[1]
                            orf = ls[2]
                            clu = ls[7]
                            hmm = ls[3]
                            memoryDict[dataset][orf]["cat"] = ls[0]
                            memoryDict[dataset][orf]["gene"] = ls[3]
                            memoryDict[dataset][orf]["bit"] = ls[4]
                            memoryDict[dataset][orf]["cutoff"] = "evalue-cutoff: 1E-10"
                            memoryDict[dataset][orf]["clu"] = ls[6]
                            memoryDict[dataset][orf]["heme"] = ls[7]
                            memoryDict[dataset][orf]["seq"] = ls[8]
                            if args.ref != "NA":
                                memoryDict[dataset][orf]["blastHit"] = ls[8]
                                memoryDict[dataset][orf]["blastEval"] = ls[9]
                                memoryDict[dataset][orf]["seq"] = ls[10]

                            geneToCatDict[hmm] = cat
                            clusterDict[ls[7]].append(hmm + "|" + dataset + "|" + orf)
                    else:
                        out.write(i.rstrip())
            print("\n\n\n\n")
            for i in clusterDict.keys():
                out.write("#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                for j in clusterDict[i]:
                    hmm = j.split("|")[0]
                    dataset = j.split("|")[1]
                    orf = j.split("|")[2]
                    cat = memoryDict[dataset][orf]["cat"]

                    if cat in ["iron_aquisition-siderophore_transport_potential", "iron_aquisition-siderophore_transport", "iron_aquisition-heme_transport"]:
                        if len(Unique2(clusterDict[i])) < 2:
                            break
                        elif check1(clusterDict[i]) < 2:
                            pass
                        else:
                            if args.ref != "NA":
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + mapDict[hmm] + "," +
                                    memoryDict[dataset][orf]["bit"] + "," + memoryDict[dataset][orf]["cutoff"] + "," +
                                    memoryDict[dataset][orf]["clu"] + "," + memoryDict[dataset][orf]["heme"] + "," +
                                    memoryDict[dataset][orf]["blastHit"] + "," + memoryDict[dataset][orf]["blastEval"] + "," +
                                    memoryDict[dataset][orf]["seq"] + "\n")
                            else:
                                out.write(memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + mapDict[hmm] + "," +
                                          memoryDict[dataset][orf]["bit"] + "," + memoryDict[dataset][orf]["cutoff"] + "," +
                                          memoryDict[dataset][orf]["clu"] + "," + memoryDict[dataset][orf]["heme"] + "," +
                                          memoryDict[dataset][orf]["seq"] + "\n")

                    elif cat in ["iron_aquisition-siderophore_synthesis"]:
                        if len(Unique2(clusterDict[i])) < 3:
                            break
                        elif check1_2(clusterDict[i]) < 3:
                            pass
                        else:
                            if args.ref != "NA":
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + mapDict[hmm] + "," +
                                    memoryDict[dataset][orf]["bit"] + "," + memoryDict[dataset][orf]["cutoff"] + "," +
                                    memoryDict[dataset][orf]["clu"] + "," + memoryDict[dataset][orf]["heme"] + "," +
                                    memoryDict[dataset][orf]["blastHit"] + "," + memoryDict[dataset][orf]["blastEval"] + "," +
                                    memoryDict[dataset][orf]["seq"] + "\n")
                            else:

                                out.write(memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + mapDict[hmm] + "," +
                                          memoryDict[dataset][orf]["bit"] + "," + memoryDict[dataset][orf]["cutoff"] + "," +
                                          memoryDict[dataset][orf]["clu"] + "," + memoryDict[dataset][orf]["heme"] + "," +
                                          memoryDict[dataset][orf]["seq"] + "\n")

                    elif cat in ["iron_aquisition-iron_transport", "iron_aquisition-heme_oxygenase"]:
                        if len(Unique2(clusterDict[i])) < 2:
                            break
                        elif check2(clusterDict[i]) < 2:
                            pass
                        else:
                            if args.ref != "NA":
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + mapDict[hmm] + "," +
                                    memoryDict[dataset][orf]["bit"] + "," + memoryDict[dataset][orf]["cutoff"] + "," +
                                    memoryDict[dataset][orf]["clu"] + "," + memoryDict[dataset][orf]["heme"] + "," +
                                    memoryDict[dataset][orf]["blastHit"] + "," + memoryDict[dataset][orf]["blastEval"] + "," +
                                    memoryDict[dataset][orf]["seq"] + "\n")
                            else:
                                out.write(memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + mapDict[hmm] + "," +
                                          memoryDict[dataset][orf]["bit"] + "," + memoryDict[dataset][orf]["cutoff"] + "," +
                                          memoryDict[dataset][orf]["clu"] + "," + memoryDict[dataset][orf]["heme"] + "," +
                                          memoryDict[dataset][orf]["seq"] + "\n")

                    elif cat == "iron_gene_regulation":
                        if checkReg(clusterDict[i]) < 1:
                            pass
                        else:
                            if args.ref != "NA":
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                    memoryDict[dataset][orf]["bit"] + "," +
                                    memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                    memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] + "," +
                                    memoryDict[dataset][orf]["blastEval"] + "," + memoryDict[dataset][orf]["seq"] + "\n")
                            else:
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                    memoryDict[dataset][orf]["bit"] + "," +
                                    memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                    memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["seq"] + "\n")

                    elif cat == "magnetosome_formation":
                        if checkMam(clusterDict[i]) < 5:
                            pass
                        else:
                            if args.ref != "NA":
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                    memoryDict[dataset][orf]["bit"] + "," +
                                    memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                    memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] + "," +
                                    memoryDict[dataset][orf]["blastEval"] + "," + memoryDict[dataset][orf]["seq"] + "\n")
                            else:
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                    memoryDict[dataset][orf]["bit"] + "," +
                                    memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                    memoryDict[dataset][orf]["heme"] + "," +
                                    memoryDict[dataset][orf]["seq"] + "\n")

                    elif cat == "iron_oxidation":
                        if hmm == "Cyc1":
                            if checkFe(clusterDict[i]) < 2:
                                pass
                            else:
                                if args.ref != "NA":
                                    out.write(
                                        memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] + "," +
                                        memoryDict[dataset][orf]["blastEval"] + "," + memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")

                        elif hmm in ["Cyc2_repCluster1", "Cyc2_repCluster2", "Cyc2_repCluster3"]:

                            seq = memoryDict[dataset][orf]["seq"]

                            numhemes = len(re.findall(r'C(..)CH', seq)) + len(re.findall(r'C(...)CH', seq)) \
                            + len(re.findall(r'C(....)CH', seq)) + len(re.findall(r'C(..............)CH', seq)) \
                            + len(re.findall(r'C(...............)CH', seq))
                            if args.ref == "NA":
                                if len(seq) >= 365:
                                    if numhemes > 0:
                                        out.write(
                                            memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                            memoryDict[dataset][orf]["bit"] + "," +
                                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                            memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["seq"] + "\n")
                            else:
                                if len(seq) >= 365:
                                    if numhemes > 0:
                                        out.write(
                                            memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                            memoryDict[dataset][orf]["bit"] + "," +
                                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                            memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf][
                                                "blastHit"] + "," +
                                            memoryDict[dataset][orf]["blastEval"] + "," + memoryDict[dataset][orf][
                                                "seq"] + "\n")


                        elif hmm in ["MtoA", "MtrA", "MtrB_TIGR03509", "MtrC_TIGR03507"]:
                            operon = clusterDict[i]
                            operon = Strip(operon)
                            if "MtrB_TIGR03509" in operon and "MtrC_TIGR03507" in operon and "MtrA" in operon:
                                if args.ref != "NA":
                                    out.write(
                                        "iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        "iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")

                            elif "MtoA" in operon and "MtrB_TIGR03509" in operon and "MtrC_TIGR03507" in operon:
                                if args.ref != "NA":
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")

                            elif "MtrB_TIGR03509" in operon and "MtrA" in operon:
                                if args.ref != "NA":
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")

                            elif "MtoA" in operon and "MtrB_TIGR03509" in operon:
                                if args.ref != "NA":
                                    out.write(
                                        "possible_iron_oxidation_and_possible_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        "possible_iron_oxidation_and_possible_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")

                            elif "MtrC_TIGR03507" in operon and "MtrB_TIGR03509" in operon:
                                if args.ref != "NA":
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")

                            elif "MtrC_TIGR03507" in operon:
                                if args.ref != "NA":
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")

                            else:
                                pass

                        else:
                            if args.ref != "NA":
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                    memoryDict[dataset][orf]["bit"] + "," +
                                    memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                    memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                    memoryDict[dataset][orf]["blastEval"] + "," +
                                    memoryDict[dataset][orf]["seq"] + "\n")
                            else:
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                    memoryDict[dataset][orf]["bit"] + "," +
                                    memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                    memoryDict[dataset][orf]["heme"] + "," +
                                    memoryDict[dataset][orf]["seq"] + "\n")

                    elif cat == "iron_reduction":
                        if hmm in ["t4ap"]:
                            aromatics = ["F", "Y", "W", "H"]
                            seq = memoryDict[dataset][orf]["seq"]
                            aromaticAAs = 0
                            aromaticFreeGap = 0
                            gapThreshold = 0
                            for aa in seq:
                                if aa in aromatics:
                                    if aromaticFreeGap > 35:
                                        gapThreshold += 1
                                    aromaticAAs += 1
                                    aromaticFreeGap = 0
                                else:
                                    aromaticFreeGap += 1

                            percAromatic = aromaticAAs / len(seq)
                            if percAromatic > 0.097 and gapThreshold == 0 and \
                                    re.findall(r'[FYWH](......................)[FYWH](..)[FYWH](....)[FYWH](.................)[FYWH][FYWH](.....)[FYWH]', seq):
                                if args.ref != "NA":
                                    out.write(
                                        memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")

                        elif hmm in ["GACE_1843", "GACE_1844", "GACE_1845", "GACE_1846", "GACE_1847"]:
                            if checkGACE(clusterDict[i]) < 3:
                                pass
                            else:
                                if args.ref != "NA":
                                    out.write(
                                        memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")


                        elif hmm in ["DFE_0465", "DFE_0464", "DFE_0463", "DFE_0462", "DFE_0461"]:
                            if checkDFE1(clusterDict[i]) < 3:
                                pass
                            else:
                                if args.ref != "NA":
                                    out.write(
                                        memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")

                        elif hmm in ["DFE_0451", "DFE_0450", "DFE_0449", "DFE_0448"]:
                            if checkDFE2(clusterDict[i]) < 3:
                                pass
                            else:
                                if args.ref != "NA":
                                    out.write(
                                        memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")

                        elif hmm in ["MtrC_TIGR03507", "MtrA", "MtrB_TIGR03509", "MtoA"]:
                            operon = clusterDict[i]
                            operon = Strip(operon)
                            if "MtrB_TIGR03509" in operon and "MtrC_TIGR03507" in operon and "MtrA" in operon:
                                if args.ref != "NA":
                                    out.write(
                                        "iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        "iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")

                            elif "MtoA" in operon and "MtrB_TIGR03509" in operon and "MtrC_TIGR03507" in operon:
                                if args.ref != "NA":
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")

                            elif "MtrB_TIGR03509" in operon and "MtrA" in operon:
                                if args.ref != "NA":
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:

                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                    # print("writing iron_reduction_or_oxidation")

                            elif "MtoA" in operon and "MtrB_TIGR03509" in operon:
                                if args.ref != "NA":
                                    out.write(
                                        "possible_iron_oxidation_and_possible_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        "possible_iron_oxidation_and_possible_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")

                            elif "MtrC_TIGR03507" in operon and "MtrB_TIGR03509" in operon:
                                if args.ref != "NA":
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")

                            elif "MtrC_TIGR03507" in operon:
                                if args.ref != "NA":
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                        memoryDict[dataset][orf]["blastEval"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                                else:
                                    out.write(
                                        "probable_iron_reduction" + "," + dataset + "," + orf + "," + hmm + "," +
                                        memoryDict[dataset][orf]["bit"] + "," +
                                        memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                        memoryDict[dataset][orf]["heme"] + "," +
                                        memoryDict[dataset][orf]["seq"] + "\n")
                            else:
                                pass

                        else:
                            if args.ref != "NA":
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                    memoryDict[dataset][orf]["bit"] + "," +
                                    memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                    memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                    memoryDict[dataset][orf]["blastEval"] + "," +
                                    memoryDict[dataset][orf]["seq"] + "\n")
                            else:
                                out.write(
                                    memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                    memoryDict[dataset][orf]["bit"] + "," +
                                    memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                    memoryDict[dataset][orf]["heme"] + "," +
                                    memoryDict[dataset][orf]["seq"] + "\n")

                    else:
                        if args.ref != "NA":
                            out.write(
                                memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                memoryDict[dataset][orf]["bit"] + "," +
                                memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                                memoryDict[dataset][orf]["blastEval"] + "," +
                                memoryDict[dataset][orf]["seq"] + "\n")
                        else:
                            out.write(
                                memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                                memoryDict[dataset][orf]["bit"] + "," +
                                memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                                memoryDict[dataset][orf]["heme"] + "," +
                                memoryDict[dataset][orf]["seq"] + "\n")

            out.close()

            time.sleep(5)
            os.system("mv %s/FeGenie-summary-fixed.csv %s/FeGenie-summary.csv" % (outDirectory, outDirectory))
        else:
            print("Final processing of output\n")
            MAP = open(HMMdir + "/FeGenie-map.txt", "r")
            mapDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in MAP:
                ls = i.rstrip().split("\t")
                mapDict[ls[0]] = ls[1]

            out = open(outDirectory + "/FeGenie-summary-fixed.csv", "w")
            geneToCatDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            memoryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
            clusterDict = defaultdict(list)
            infile = open(outDirectory + "/FeGenie-summary.csv")
            for i in infile:
                if not re.match(r'#', i):
                    ls = i.rstrip().split(",")
                    if ls[6] != "cluster":

                        if not re.findall(r'defaultdict', ls[5]):
                            clu = ls[6]
                            cat = ls[0]
                            dataset = ls[1]
                            orf = ls[2]
                            hmm = allButTheLast(ls[3], ".")
                            clusterDict[clu].append(hmm + "|" + dataset + "|" + orf)
                            geneToCatDict[hmm] = cat
                            hmm = allButTheLast(ls[3], ".")
                            memoryDict[dataset][orf]["cat"] = ls[0]
                            memoryDict[dataset][orf]["gene"] = ls[3]
                            memoryDict[dataset][orf]["bit"] = ls[4]
                            memoryDict[dataset][orf]["cutoff"] = ls[5]
                            memoryDict[dataset][orf]["clu"] = clu
                            memoryDict[dataset][orf]["heme"] = ls[7]
                            memoryDict[dataset][orf]["seq"] = ls[8]
                            if args.ref != "NA":
                                memoryDict[dataset][orf]["blastHit"] = ls[8]
                                memoryDict[dataset][orf]["blastEval"] = ls[9]
                                memoryDict[dataset][orf]["seq"] = ls[10]
                        else:
                            cat = ls[0]
                            dataset = ls[1]
                            orf = ls[2]
                            clu = ls[7]
                            hmm = ls[3]
                            memoryDict[dataset][orf]["cat"] = ls[0]
                            memoryDict[dataset][orf]["gene"] = ls[3]
                            memoryDict[dataset][orf]["bit"] = ls[4]
                            memoryDict[dataset][orf]["cutoff"] = "evalue-cutoff: 1E-10"
                            memoryDict[dataset][orf]["clu"] = ls[6]
                            memoryDict[dataset][orf]["heme"] = ls[7]
                            memoryDict[dataset][orf]["seq"] = ls[8]
                            if args.ref != "NA":
                                memoryDict[dataset][orf]["blastHit"] = ls[8]
                                memoryDict[dataset][orf]["blastEval"] = ls[9]
                                memoryDict[dataset][orf]["seq"] = ls[10]


                            geneToCatDict[hmm] = cat
                            clusterDict[ls[7]].append(hmm + "|" + dataset + "|" + orf)
                    else:
                        out.write(i.rstrip())

            for i in clusterDict.keys():
                out.write(
                    "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "," + "#" + "\n")
                for j in clusterDict[i]:
                    hmm = j.split("|")[0]
                    dataset = j.split("|")[1]
                    orf = j.split("|")[2]
                    cat = memoryDict[dataset][orf]["cat"]

                    if args.ref != "NA":
                        out.write(
                            memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," + memoryDict[dataset][orf]["blastHit"] +
                            memoryDict[dataset][orf]["blastEval"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")
                    else:
                        out.write(
                            memoryDict[dataset][orf]["cat"] + "," + dataset + "," + orf + "," + hmm + "," +
                            memoryDict[dataset][orf]["bit"] + "," +
                            memoryDict[dataset][orf]["cutoff"] + "," + memoryDict[dataset][orf]["clu"] + "," +
                            memoryDict[dataset][orf]["heme"] + "," +
                            memoryDict[dataset][orf]["seq"] + "\n")

            out.close()
            time.sleep(5)
            if not args.orfs:
                os.system("mv %s/FeGenie-summary-fixed.csv %s/FeGenie-summary.csv" % (outDirectory, outDirectory))

        # ****************************** PRE-FINAL ALTERATION OF THE OUTPUT FILE ***************************************
        clu = 0
        summaryDict = defaultdict(list)
        summary = open(outDirectory + "/FeGenie-summary.csv")

        out = open(outDirectory + "/FeGenie-summary-altered.csv", "w")
        if args.ref != "NA":
            out.write(
                "category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "clusterID" + "," + "heme_c_binding_motifs" + "," + "heme_b_binding_motifs" + "," + "hematite_binding_motifs" + "," + "top_blast_hit" + "," + "blast_hit_evalue" + "," + "protein_sequence" + "\n")
        else:
            out.write(
                "category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "clusterID" + "," + "heme_c_binding_motifs" + "," + "heme_b_binding_motifs" + "," + "hematite_binding_motifs" + "," + "protein_sequence" + "\n")

        for i in summary:
            if re.search(r'#', i):
                clu += 1
            else:
                summaryDict[clu].append(i.rstrip())

        if args.gbk:
            idxDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for idxfile in binDirLS:
                if lastItem(idxfile.split(".")) == "idx":
                    idxfileopen = open("%s/%s" % (binDir, idxfile))
                    for idxline in idxfileopen:
                        ls = idxline.rstrip().split(",")
                        newOrf = ls[1]
                        oldOrf = ls[0]
                        idxDict[newOrf] = oldOrf

            print("\n")
            for i in summaryDict.keys():
                if len(summaryDict[i]) > 0:
                    if args.ref != "NA":
                        for j in summaryDict[i]:
                            ls = j.split(",")
                            seq = ls[10]
                            hemeb = len(re.findall(r'G(.)[HR](.)C[PLAV]G', seq))
                            hbm = len(re.findall(r'[ST][AVILMFYWH][ST]P[ST]', seq))
                            out.write(ls[0] + "," + ls[1] + "," + str(idxDict[ls[2]]) + "," + ls[3] + "," + ls[4] + "," + ls[
                                    5] + "," + ls[6] + "," + ls[7] + "," + str(hemeb) + "," + str(hbm) + "," + ls[8] + "," + ls[9] + "," + ls[10] + "\n")
                        out.write("##########################################################################################################################################################################################################\n")
                    else:
                        for j in summaryDict[i]:
                            ls = j.split(",")
                            seq = ls[8]
                            hemeb = len(re.findall(r'G(.)[HR](.)C[PLAV]G', seq))
                            hbm = len(re.findall(r'[ST][AVILMFYWH][ST]P[ST]', seq))
                            out.write(
                                ls[0] + "," + ls[1] + "," + str(idxDict[ls[2]]) + "," + ls[3] + "," + ls[4] + "," + ls[
                                    5] + "," + ls[6] + "," + ls[7] + "," + str(hemeb) + "," + str(hbm) + "," + ls[8] + "\n")
                        out.write(
                            "#####################################################################################################"
                            "#####################################################################################################\n")
        else:
            for i in summaryDict.keys():
                if len(summaryDict[i]) > 0:
                    if args.ref != "NA":
                        for j in summaryDict[i]:
                            ls = j.split(",")
                            seq = ls[8]
                            hemeb = len(re.findall(r'G(.)[HR](.)C[PLAV]G', seq))
                            hbm = len(re.findall(r'[STC][AVILMFYWH][ST]P[ST]', seq))
                            out.write(
                                ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + ls[4] + "," + ls[
                                    5] + "," + ls[6] + "," + ls[7] + "," + str(hemeb) + "," + str(hbm) + "," + ls[
                                    8] + "," + ls[9] + "," + ls[10] + "\n")
                        out.write(
                            "#####################################################################################################"
                            "#####################################################################################################\n")
                    else:
                        for j in summaryDict[i]:
                            ls = j.split(",")
                            seq = ls[8]
                            hemeb = len(re.findall(r'G(.)[HR](.)C[PLAV]G', seq))
                            hbm = len(re.findall(r'[STC][AVILMFYWH][ST]P[ST]', seq))
                            out.write(
                                ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + ls[4] + "," + ls[
                                    5] + "," + ls[6] + "," + ls[7] + "," + str(hemeb) + "," + str(hbm) + "," + ls[8] + "\n")
                        out.write(
                            "#####################################################################################################"
                            "#####################################################################################################\n")

        out.close()

        time.sleep(5)
        os.system("mv %s/FeGenie-summary-altered.csv %s/FeGenie-geneSummary-clusters.csv" % (outDirectory, outDirectory))
        os.system("rm %s/FeGenie-summary.csv" % outDirectory)

        # ****************************** REMOVING #'S ***************************************
        summary = open("%s/FeGenie-geneSummary-clusters.csv" % outDirectory, "r")
        out = open("%s/FeGenie-geneSummary.csv" % outDirectory, "w")
        if args.ref != "NA":
            out.write(
                "category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "clusterID" + "," + "heme_c_binding_motifs" + "," + "top_blast_hit" + "," + "blast_hit_evalue" + "," + "protein_sequence" + "\n")
        else:
            out.write(
                "category" + "," + "genome/assembly" + "," + "orf" + "," + "HMM" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "clusterID" + "," + "heme_c_binding_motifs" + "," + "protein_sequence" + "\n")

        for i in summary:
            if not re.search(r'#', i):
                out.write(i.rstrip() + "\n")

        out.close()

        time.sleep(5)
        try:
            hmmout = os.listdir("%s/HMM_results" % outDirectory)
            os.system("rm -rf %s/HMM_results/*" % outDirectory)
            os.system("mv %s/*-HMM %s/HMM_results/" % (outDirectory, outDirectory))
        except FileNotFoundError:
            os.system("mkdir %s/HMM_results" % outDirectory)
            os.system("mv %s/*-HMM %s/HMM_results/" % (outDirectory, outDirectory))

        # if prodigal == 1:
        #     os.system("rm %s/ORF_calls/*-prodigal.out" % outDirectory)

        os.system("rm -rf %s/makedbfile.txt.perf" % outDirectory)

        if args.heme:
            print("Looking for heme-binding motifs")
            out = open("%s/FeGenie-hemeProteins.csv" % outDirectory, "w")
            out.write("genome/assembly,orf,heme_c_binding_motifs,heme_b_binding_motifs,seq\n")
            for i in binDirLS:
                if lastItem(i.split(".")) == args.bin_ext:  # FILTERING OUT ANY NON-BIN-RELATED FILES
                    if args.orfs:
                        BIN = open("%s/%s" % (binDir, i))
                        BIN = fasta(BIN)
                        for j in BIN.keys():
                            seq = BIN[j]
                            hemec = len(re.findall(r'C(..)CH', seq)) + len(re.findall(r'C(...)CH', seq)) \
                                    + len(re.findall(r'C(....)CH', seq)) + len(re.findall(r'C(..............)CH', seq)) \
                                    + len(re.findall(r'C(...............)CH', seq))
                            hemeb = len(re.findall(r'G(.)[HR](.)C[PLAV]G', seq))
                            if hemec > 0 or hemeb > 0:
                                out.write(i + "," + j + "," + str(hemec) + "," + str(hemeb) + "," + seq + "\n")

                    else:
                        BIN = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i))
                        BIN = fasta(BIN)
                        for j in BIN.keys():
                            seq = BIN[j]
                            hemec = len(re.findall(r'C(..)CH', seq)) + len(re.findall(r'C(...)CH', seq)) \
                                    + len(re.findall(r'C(....)CH', seq)) + len(re.findall(r'C(..............)CH', seq)) \
                                    + len(re.findall(r'C(...............)CH', seq))
                            hemeb = len(re.findall(r'G(.)[HR](.)C[PLAV]G', seq))
                            if hemec > 0 or hemeb > 0:
                                out.write(i + "," + j + "," + str(hemec) + "," + str(hemeb) + "," + seq + "\n")
            out.close()

        if args.hematite:
            print("Looking for hematite-binding motifs")
            positive = ["R", "H", "K"]
            negative = ["D", "E"]
            polar = ["S", "T", "N", "Q"]
            special = ["C", "G", "P"]
            hydrophobic = ["A", "V", "I", "L", "M", "F", "Y", "W"]

            out = open("%s/FeGenie-hematiteProteins.csv" % outDirectory, "w")
            out.write("genome/assembly,orf,hematite_binding_motifs,seq\n")
            for i in binDirLS:
                if lastItem(i.split(".")) == args.bin_ext:  # FILTERING OUT ANY NON-BIN-RELATED FILES
                    if args.orfs:
                        BIN = open("%s/%s" % (binDir, i))
                        BIN = fasta(BIN)
                        for j in BIN.keys():
                            seq = BIN[j]
                            hbm = len(re.findall(r'[STC][AVILMFYWH][ST]P[ST]', seq))
                            if hbm > 0:
                                out.write(i + "," + j + "," + str(hbm) + "," + seq + "\n")

                    else:
                        BIN = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i))
                        BIN = fasta(BIN)
                        for j in BIN.keys():
                            seq = BIN[j]
                            hbm = len(re.findall(r'[STC][AVILMFYWH][ST]P[ST]', seq))
                            if hbm > 0:
                                out.write(i + "," + j + "," + str(hbm) + "," + seq + "\n")
            out.close()

        time.sleep(5)

        print("Writen summary to file: %s/FeGenie-geneSummary-clusters.csv for visual inspection" % outDirectory)
        print("Writen summary to file: %s/FeGenie-geneSummary.csv for downstream parsing and analyses" % outDirectory)
        # ****************************** CREATING A HEATMAP-COMPATIBLE CSV FILE *************************************
        print("Writing heatmap-formatted output file: %s/FeGenie-heatmap-data.csv\n" % outDirectory)

        # GENE-COUNTS BASED ABUNDANCE

        if args.bam == "NA" and args.bams == "NA":
            cats = ["iron_aquisition-iron_transport", "iron_aquisition-heme_transport", "iron_aquisition-heme_oxygenase",
                    "iron_aquisition-siderophore_synthesis", "iron_aquisition-siderophore_transport",
                    "iron_aquisition-siderophore_transport_potential", "iron_gene_regulation", "iron_oxidation",
                    "possible_iron_oxidation_and_possible_iron_reduction", "probable_iron_reduction",
                    "iron_reduction", "iron_storage", "magnetosome_formation"]

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/FeGenie-geneSummary-clusters.csv" % outDirectory, "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if not re.search(r'#', i):
                    if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome" and ls[1] != "genome/assembly":
                        if not re.match(r'#', i):
                            process = ls[0]
                            cell = ls[1]
                            orf = ls[2]
                            gene = ls[3]
                            Dict[cell][process].append(gene)

            normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in os.listdir(args.bin_dir):
                if lastItem(i.split(".")) == args.bin_ext and not re.findall(r'-proteins.faa', i):
                    if args.orfs:
                        file = open("%s/%s" % (args.bin_dir, i), "r")
                    else:
                        file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                    file = fasta(file)
                    normDict[i] = len(file.keys())

            outHeat = open("%s/FeGenie-heatmap-data.csv" % (outDirectory), "w")
            outHeat.write("X")
            for i in sorted(Dict.keys()):
                outHeat.write("," + i)
            outHeat.write("\n")

            for i in cats:
                outHeat.write(i + ",")
                for j in sorted(Dict.keys()):
                    if not re.match(r'#', j):
                        if args.norm:
                            outHeat.write("," + str((len(Dict[j][i]) / int(normDict[j])) * float(args.inflation)))
                        else:
                            outHeat.write("," + str(len(Dict[j][i])))
                outHeat.write("\n")

            outHeat.close()

        elif args.bams != "NA":

            # COVERAGE-BASED ABUNDANCE WITH MULTIPLE BAM FILES

            cats = ["iron_aquisition-iron_transport", "iron_aquisition-heme_transport", "iron_aquisition-heme_oxygenase",
                    "iron_aquisition-siderophore_synthesis", "iron_aquisition-siderophore_transport",
                    "iron_aquisition-siderophore_transport_potential", "iron_gene_regulation", "iron_oxidation",
                    "possible_iron_oxidation_and_possible_iron_reduction", "probable_iron_reduction",
                    "iron_reduction", "iron_storage", "magnetosome_formation"]

            depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            BAMmapDict = defaultdict(lambda: defaultdict(lambda: "EMPTY"))
            normDict = defaultdict(lambda: defaultdict(lambda: "EMPTY"))
            BAMmap = open(args.bams)
            for bamLine in BAMmap:
                string = ''
                ls = bamLine.rstrip().split("\t")
                cell = ls[0]
                for j in ls[1:]:
                    string += " "
                    string += j

                try:
                    print("processing... " + cell)
                    depth = open("%s/contigDepths/%s.depth" % (args.out, cell))
                    total = 0
                    for k in depth:
                        LS = k.rstrip().split("\t")
                        if LS[0] != "contigName":
                            if args.which_bams == "average":
                                depthDict[cell][LS[0]] = LS[2]
                            else:
                                column = int(args.which_bams)*2 + 1
                                depthDict[cell][LS[0]] = LS[column]
                            total += float(LS[2])
                    normDict[cell] = total

                except FileNotFoundError:
                    os.system("jgi_summarize_bam_contig_depths --outputDepth %s/%s.depth%s" % (outDirectory, cell, string))
                    print("processing... " + cell)
                    depth = open("%s/%s.depth" % (outDirectory, cell))
                    total = 0
                    for k in depth:
                        LS = k.rstrip().split("\t")
                        if LS[0] != "contigName":
                            depthDict[cell][LS[0]] = LS[2]
                            total += float(LS[2])
                    normDict[cell] = total

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/FeGenie-geneSummary-clusters.csv" % outDirectory, "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if not re.search(r'#', i):
                    if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome" and ls[1] != "genome/assembly":
                        if not re.match(r'#', i):
                            process = ls[0]
                            cell = ls[1]
                            orf = ls[2]
                            contig = allButTheLast(orf, "_")
                            gene = ls[3]
                            Dict[cell][process].append(float(depthDict[cell][contig]))

            outHeat = open("%s/FeGenie-%s-heatmap-data.csv" % (outDirectory, args.which_bams), "w")
            outHeat.write("X")
            for i in sorted(Dict.keys()):
                outHeat.write("," + i)
            outHeat.write("\n")

            for i in cats:
                outHeat.write(i)
                for j in sorted(Dict.keys()):
                    if not re.match(r'#', j):
                        outHeat.write("," + str(SUM(Dict[j][i])))
                outHeat.write("\n")

            outHeat.close()

        elif args.bam != "NA":

            # COVERAGE-BASED ABUNDANCE USING ONLY ONE BAM FILE

            cats = ["iron_aquisition-iron_transport", "iron_aquisition-heme_transport", "iron_aquisition-heme_oxygenase",
                    "iron_aquisition-siderophore_synthesis", "iron_aquisition-siderophore_transport",
                    "iron_aquisition-siderophore_transport_potential", "iron_gene_regulation", "iron_oxidation",
                    "possible_iron_oxidation_and_possible_iron_reduction", "probable_iron_reduction",
                    "iron_reduction", "iron_storage", "magnetosome_formation"]

            depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))

            try:
                total = 0
                depth = open("%s.depth" % (args.bam))
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[LS[0]] = LS[2]
                        total += float(LS[2])

            except FileNotFoundError:
                os.system("jgi_summarize_bam_contig_depths --outputDepth %s.depth %s" % (args.bam, args.bam))
                depth = open("%s.depth" % (args.bam))
                total = 0
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[LS[0]] = LS[2]
                        total += float(LS[2])
                os.system("mv %s.depth %s/" % (args.bam, outDirectory))

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/FeGenie-geneSummary-clusters.csv" % outDirectory, "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if not re.search(r'#', i):
                    if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome" and ls[1] != "genome/assembly":
                        if not re.match(r'#', i):
                            process = ls[0]
                            cell = ls[1]
                            orf = ls[2]
                            contig = allButTheLast(orf, "_")
                            gene = ls[3]
                            Dict[cell][process].append(float(depthDict[contig]))

            outHeat = open("%s/FeGenie-heatmap-data.csv" % (outDirectory), "w")
            outHeat.write("X")
            for i in sorted(Dict.keys()):
                outHeat.write("," + i)
            outHeat.write("\n")

            for i in cats:
                outHeat.write(i)
                for j in sorted(Dict.keys()):
                    if not re.match(r'#', j):
                        outHeat.write("," + str(SUM(Dict[j][i])))
                outHeat.write("\n")

        outHeat.close()
        time.sleep(5)

        # ******** RUNNING RSCRIPT TO GENERATE PLOTS **************
        if args.makeplots:
            print("Running Rscript to generate plots. Do not be alarmed if you see Warning or Error messages from Rscript. "
                  "This will not affect any of the output data that was already created. If you see plots generated, great! "
                  "If not, you can plot the data as you wish on your own, or start an issue on FeGenie's GitHub repository")

            if args.bams != "NA":
                if args.norm:
                    os.system("Rscript --vanilla %s/DotPlot.R %s/FeGenie-%s-heatmap-data.csv %s/" % (rscriptDir, outDirectory, args.which_bams, outDirectory))
                    os.system(
                        "Rscript --vanilla %s/dendro-heatmap.R %s/FeGenie-%s-heatmap-data.csv %s/" % (rscriptDir, outDirectory, args.which_bams, outDirectory))
                else:
                    os.system(
                        "Rscript --vanilla %s/DotPlot-nonorm.R %s/FeGenie-%s-heatmap-data.csv %s/" % (rscriptDir, outDirectory, args.which_bams, outDirectory))
                    os.system(
                        "Rscript --vanilla %s/dendro-heatmap.R %s/FeGenie-%s-heatmap-data.csv %s/" % (rscriptDir, outDirectory, args.which_bams, outDirectory))

                os.system("mv %s/Fegenie-dotplot.tiff %s/Fegenie-%s-dotplot.tiff" % (outDirectory, outDirectory, args.which_bams))

            else:
                if args.norm:
                    os.system("Rscript --vanilla %s/DotPlot.R %s/FeGenie-heatmap-data.csv %s/" % (rscriptDir, outDirectory, outDirectory))
                    os.system(
                        "Rscript --vanilla %s/dendro-heatmap.R %s/FeGenie-heatmap-data.csv %s/" % (rscriptDir, outDirectory, outDirectory))
                else:
                    os.system(
                        "Rscript --vanilla %s/DotPlot-nonorm.R %s/FeGenie-heatmap-data.csv %s/" % (rscriptDir, outDirectory, outDirectory))
                    os.system(
                        "Rscript --vanilla %s/dendro-heatmap.R %s/FeGenie-heatmap-data.csv %s/" % (rscriptDir, outDirectory, outDirectory))

            print("\n\n\n")
            print("...")

            # ******** CHECKING ON SUCCESS OF RSCRIPTS ***********
            DIR = os.listdir(outDirectory)
            count = 0
            for i in DIR:
                if not re.match(r'\.', i):
                    count += 1

            if count == 2:
                print(
                    "Looks like Rscript has not performed succesfully. This, unfortunately, is a very finicky part of the pipeline. "
                    "The CSV files have, nonetheless, been successfully created, so you can take that data and plot if manually as you wish. "
                    "Also, feel free to start an Issue on FeGenie's GitHub page, by posting the error that was printed during the Rscript command.")

            if count > 2 and count < 5:
                print(
                    "Looks like at least one plot was generated by Rscript, but there was likely an error with one of the scripts. "
                    "Don't be alarmed if you see some error or warning messages in the terminal window. "
                    "The main CSV output should be present, so that you can plot the data as you wish on your own. "
                    "Also, feel free to start an Issue on FeGenie's GitHub page by posting the error that was printed during the Rscript command.")

            if count == 5:
                print(
                    "Looks like Rscript ran succesfully! Congrats on this. Hopefully, the resulting plots are of use to you.")

        print("")
        print("Pipeline finished without crashing!!! Thanks for using :)")

    else:
        # ****************************** CREATING A HEATMAP-COMPATIBLE CSV FILE *************************************
        if args.heme:
            print("Looking for heme-binding motifs")
            out = open("%s/FeGenie-hemeProteins.csv" % outDirectory, "w")
            out.write("genome/assembly,orf,heme_c_binding_motifs,heme_b_binding_motifs,seq\n")
            for i in binDirLS:
                if lastItem(i.split(".")) == args.bin_ext:  # FILTERING OUT ANY NON-BIN-RELATED FILES
                    if args.orfs:
                        BIN = open("%s/%s" % (binDir, i))
                        BIN = fasta(BIN)
                        for j in BIN.keys():
                            seq = BIN[j]
                            hemec = len(re.findall(r'C(..)CH', seq)) + len(re.findall(r'C(...)CH', seq)) \
                                    + len(re.findall(r'C(....)CH', seq)) + len(re.findall(r'C(..............)CH', seq)) \
                                    + len(re.findall(r'C(...............)CH', seq))
                            hemeb = len(re.findall(r'G(.)[HR](.)C[PLAV]G', seq))
                            if hemec > 0 or hemeb > 0:
                                out.write(i + "," + j + "," + str(hemec) + "," + str(hemeb) + "," + seq + "\n")

                    else:
                        BIN = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i))
                        BIN = fasta(BIN)
                        for j in BIN.keys():
                            seq = BIN[j]
                            hemec = len(re.findall(r'C(..)CH', seq)) + len(re.findall(r'C(...)CH', seq)) \
                                    + len(re.findall(r'C(....)CH', seq)) + len(re.findall(r'C(..............)CH', seq)) \
                                    + len(re.findall(r'C(...............)CH', seq))
                            hemeb = len(re.findall(r'G(.)[HR](.)C[PLAV]G', seq))
                            if hemec > 0 or hemeb > 0:
                                out.write(i + "," + j + "," + str(hemec) + "," + str(hemeb) + "," + seq + "\n")
            out.close()

        if args.hematite:
            print("Looking for hematite-binding motifs")
            positive = ["R", "H", "K"]
            negative = ["D", "E"]
            polar = ["S", "T", "N", "Q"]
            special = ["C", "G", "P"]
            hydrophobic = ["A", "V", "I", "L", "M", "F", "Y", "W"]

            out = open("%s/FeGenie-hematiteProteins.csv" % outDirectory, "w")
            out.write("genome/assembly,orf,hematite_binding_motifs,seq\n")
            for i in binDirLS:
                if lastItem(i.split(".")) == args.bin_ext:  # FILTERING OUT ANY NON-BIN-RELATED FILES
                    if args.orfs:
                        BIN = open("%s/%s" % (binDir, i))
                        BIN = fasta(BIN)
                        for j in BIN.keys():
                            seq = BIN[j]
                            hbm = len(re.findall(r'[STC][AVILMFYWH][ST]P[ST]', seq))
                            if hbm > 0:
                                out.write(i + "," + j + "," + str(hbm) + "," + seq + "\n")

                    else:
                        BIN = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i))
                        BIN = fasta(BIN)
                        for j in BIN.keys():
                            seq = BIN[j]
                            hbm = len(re.findall(r'[STC][AVILMFYWH][ST]P[ST]', seq))
                            if hbm > 0:
                                out.write(i + "," + j + "," + str(hbm) + "," + seq + "\n")
            out.close()

        time.sleep(5)

        print("Writing heatmap-formatted output file: %s/FeGenie-heatmap-data.csv\n" % outDirectory)
        # GENE-COUNTS BASED ABUNDANCE
        if args.bam == "NA" and args.bams == "NA":
            cats = ["iron_aquisition-iron_transport", "iron_aquisition-heme_transport",
                    "iron_aquisition-heme_oxygenase",
                    "iron_aquisition-siderophore_synthesis", "iron_aquisition-siderophore_transport",
                    "iron_aquisition-siderophore_transport_potential", "iron_gene_regulation", "iron_oxidation",
                    "possible_iron_oxidation_and_possible_iron_reduction", "probable_iron_reduction",
                    "iron_reduction", "iron_storage", "magnetosome_formation"]

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/FeGenie-geneSummary-clusters.csv" % outDirectory, "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if not re.search(r'#', i):
                    if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome" and ls[1] != "genome/assembly":
                        if not re.match(r'#', i):
                            process = ls[0]
                            cell = ls[1]
                            orf = ls[2]
                            gene = ls[3]
                            Dict[cell][process].append(gene)

            normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            for i in os.listdir(args.bin_dir):
                if lastItem(i.split(".")) == args.bin_ext and not re.findall(r'-proteins.faa', i):
                    if args.orfs:
                        file = open("%s/%s" % (args.bin_dir, i), "r")
                    else:
                        file = open("%s/ORF_calls/%s-proteins.faa" % (outDirectory, i), "r")
                    file = fasta(file)
                    normDict[i] = len(file.keys())

            outHeat = open("%s/FeGenie-heatmap-data.csv" % (outDirectory), "w")
            outHeat.write("X")
            for i in sorted(Dict.keys()):
                outHeat.write("," + i)
            outHeat.write("\n")

            for i in cats:
                outHeat.write(i)
                for j in sorted(Dict.keys()):
                    if not re.match(r'#', j):
                        if args.norm:
                            outHeat.write("," + str((len(Dict[j][i]) / int(normDict[j])) * float(args.inflation)))
                        else:
                            outHeat.write("," + str(len(Dict[j][i])))
                outHeat.write("\n")

            outHeat.close()

        elif args.bams != "NA":

            # COVERAGE-BASED ABUNDANCE WITH MULTIPLE BAM FILES

            cats = ["iron_aquisition-iron_transport", "iron_aquisition-heme_transport",
                    "iron_aquisition-heme_oxygenase",
                    "iron_aquisition-siderophore_synthesis", "iron_aquisition-siderophore_transport",
                    "iron_aquisition-siderophore_transport_potential", "iron_gene_regulation", "iron_oxidation",
                    "possible_iron_oxidation_and_possible_iron_reduction", "probable_iron_reduction",
                    "iron_reduction", "iron_storage", "magnetosome_formation"]

            depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
            BAMmapDict = defaultdict(lambda: defaultdict(lambda: "EMPTY"))
            normDict = defaultdict(lambda: defaultdict(lambda: "EMPTY"))
            BAMmap = open(args.bams)
            for i in BAMmap:
                string = ''
                ls = i.rstrip().split("\t")
                cell = ls[0]
                for j in ls[1:]:
                    string += " "
                    string += j

                try:
                    print("processing... " + cell)
                    depth = open("%s/contigDepths/%s.depth" % (args.out, cell))
                    total = 0
                    for k in depth:
                        LS = k.rstrip().split("\t")
                        if LS[0] != "contigName":
                            depthDict[cell][LS[0]] = LS[2]
                            total += float(LS[2])
                    normDict[cell] = total

                except FileNotFoundError:
                    os.system(
                        "jgi_summarize_bam_contig_depths --outputDepth %s/%s.depth%s" % (outDirectory, cell, string))
                    print("processing... " + cell)
                    depth = open("%s/%s.depth" % (outDirectory, cell))
                    total = 0
                    for k in depth:
                        LS = k.rstrip().split("\t")
                        if LS[0] != "contigName":
                            if args.which_bams == "average":
                                depthDict[cell][LS[0]] = LS[2]
                            else:
                                column = int(args.which_bams)*2 + 1
                                depthDict[cell][LS[0]] = LS[column]
                            total += float(LS[2])
                    normDict[cell] = total

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/FeGenie-geneSummary-clusters.csv" % outDirectory, "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if not re.search(r'#', i):
                    if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome" and ls[1] != "genome/assembly":
                        if not re.match(r'#', i):
                            process = ls[0]
                            cell = ls[1]
                            orf = ls[2]
                            contig = allButTheLast(orf, "_")
                            gene = ls[3]
                            Dict[cell][process].append(float(depthDict[cell][contig]))

            outHeat = open("%s/FeGenie-%s-heatmap-data.csv" % (outDirectory, args.which_bams), "w")
            outHeat.write("X")
            for i in sorted(Dict.keys()):
                outHeat.write("," + i)
            outHeat.write("\n")

            for i in cats:
                outHeat.write(i)
                for j in sorted(Dict.keys()):
                    if not re.match(r'#', j):
                        outHeat.write("," + str(SUM(Dict[j][i])))
                outHeat.write("\n")

            outHeat.close()

        elif args.bam != "NA":

            # COVERAGE-BASED ABUNDANCE USING ONLY ONE BAM FILE

            cats = ["iron_aquisition-iron_transport", "iron_aquisition-heme_transport",
                    "iron_aquisition-heme_oxygenase",
                    "iron_aquisition-siderophore_synthesis", "iron_aquisition-siderophore_transport",
                    "iron_aquisition-siderophore_transport_potential", "iron_gene_regulation", "iron_oxidation",
                    "possible_iron_oxidation_and_possible_iron_reduction", "probable_iron_reduction",
                    "iron_reduction", "iron_storage", "magnetosome_formation"]

            depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))

            try:
                total = 0
                depth = open("%s.depth" % (args.bam))
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[LS[0]] = LS[2]
                        total += float(LS[2])

            except FileNotFoundError:
                os.system("jgi_summarize_bam_contig_depths --outputDepth %s.depth %s" % (args.bam, args.bam))
                depth = open("%s.depth" % (args.bam))
                total = 0
                for k in depth:
                    LS = k.rstrip().split("\t")
                    if LS[0] != "contigName":
                        depthDict[LS[0]] = LS[2]
                        total += float(LS[2])
                os.system("mv %s.depth %s/" % (args.bam, outDirectory))

            Dict = defaultdict(lambda: defaultdict(list))
            final = open("%s/FeGenie-geneSummary-clusters.csv" % outDirectory, "r")
            for i in final:
                ls = (i.rstrip().split(","))
                if not re.search(r'#', i):
                    if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome" and ls[1] != "genome/assembly":
                        if not re.match(r'#', i):
                            process = ls[0]
                            cell = ls[1]
                            orf = ls[2]
                            contig = allButTheLast(orf, "_")
                            gene = ls[3]
                            Dict[cell][process].append(float(depthDict[contig]))

            outHeat = open("%s/FeGenie-heatmap-data.csv" % (outDirectory), "w")
            outHeat.write("X")
            for i in sorted(Dict.keys()):
                outHeat.write("," + i)
            outHeat.write("\n")

            for i in cats:
                outHeat.write(i)
                for j in sorted(Dict.keys()):
                    if not re.match(r'#', j):
                        outHeat.write("," + str(SUM(Dict[j][i])))
                outHeat.write("\n")

        outHeat.close()
        time.sleep(5)

        # ******** RUNNING RSCRIPT TO GENERATE PLOTS **************
        if args.makeplots:
            print(
                "Running Rscript to generate plots. Do not be alarmed if you see Warning or Error messages from Rscript. "
                "This will not affect any of the output data that was already created. If you see plots generated, great! "
                "If not, you can plot the data as you wish on your own, or start an issue on FeGenie's GitHub repository")

            if args.bams != "NA":
                if args.norm:
                    os.system("Rscript --vanilla %s/DotPlot.R %s/FeGenie-%s-heatmap-data.csv %s/" % (rscriptDir, outDirectory, args.which_bams, outDirectory))
                    os.system(
                        "Rscript --vanilla %s/dendro-heatmap.R %s/FeGenie-%s-heatmap-data.csv %s/" % (rscriptDir, outDirectory, args.which_bams, outDirectory))
                else:
                    os.system(
                        "Rscript --vanilla %s/DotPlot-nonorm.R %s/FeGenie-%s-heatmap-data.csv %s/" % (rscriptDir, outDirectory, args.which_bams, outDirectory))
                    os.system(
                        "Rscript --vanilla %s/dendro-heatmap.R %s/FeGenie-%s-heatmap-data.csv %s/" % (rscriptDir, outDirectory, args.which_bams, outDirectory))

                os.system("mv %s/Fegenie-dotplot.tiff %s/Fegenie-%s-dotplot.tiff" % (outDirectory, outDirectory, args.which_bams))

            else:

                if args.norm:
                    os.system("Rscript --vanilla %s/DotPlot.R %s/FeGenie-heatmap-data.csv %s/" % (
                    rscriptDir, outDirectory, outDirectory))

                    os.system("Rscript --vanilla %s/dendro-heatmap.R %s/FeGenie-heatmap-data.csv %s/" % (
                        rscriptDir, outDirectory, outDirectory))
                else:
                    os.system("Rscript --vanilla %s/DotPlot-nonorm.R %s/FeGenie-heatmap-data.csv %s/" % (
                        rscriptDir, outDirectory, outDirectory))

                    os.system("Rscript --vanilla %s/dendro-heatmap.R %s/FeGenie-heatmap-data.csv %s/" % (
                        rscriptDir, outDirectory, outDirectory))

            print("\n\n\n")
            print("...")

            # ******** CHECKING ON SUCCESS OF RSCRIPTS ***********
            DIR = os.listdir(outDirectory)
            count = 0
            for i in DIR:
                if not re.match(r'\.', i):
                    count += 1

            if count == 2:
                print(
                    "Looks like Rscript has not performed succesfully. This, unfortunately, is a very finicky part of the pipeline. "
                    "The CSV files have, nonetheless, been successfully created, so you can take that data and plot if manually as you wish. "
                    "Also, feel free to start an Issue on FeGenie's GitHub page, by posting the error that was printed during the Rscript command.")

            if count > 2 and count < 5:
                print(
                    "Looks like at least one plot was generated by Rscript, but there was likely an error with one of the scripts. "
                    "Don't be alarmed if you see some error or warning messages in the terminal window. "
                    "The main CSV output should be present, so that you can plot the data as you wish on your own. "
                    "Also, feel free to start an Issue on FeGenie's GitHub page by posting the error that was printed during the Rscript command.")

            if count == 5:
                print(
                    "Looks like Rscript ran succesfully! Congrats on this. Hopefully, the resulting plots are of use to you.")

        print("")
        print("Pipeline finished without crashing!!! Thanks for using :)")


if __name__ == '__main__':
    main()



