import sys
import glob
import re
import itertools
import pandas as pd
import numpy as np

CHROMOSOMES = list(range(1, 8))
SUBGENOMES = ['A', 'B', 'D']

chroms = []
chromset = []
stringy = []
for chrom in CHROMOSOMES:
    alist = list(map(''.join, itertools.chain(itertools.product([str(chrom)], SUBGENOMES))))
    chromset.append(alist)
    for _ in chromset:
        stringy.append(" ".join(str(elm) for elm in _))
    for _ in alist:
        chroms.append(_)

comb = {}
for chrom in CHROMOSOMES:
    comb[chrom] = []
    chromss = list(map(''.join, itertools.chain(itertools.product([str(chrom)], SUBGENOMES))))
    # print(chromss)
    for x, y in itertools.permutations(chromss, 2):
        comb[chrom].append((x, y))

allcomb = [comb[x] for x in comb]

TRIPLETS = None

#alis, cpares = [], [], [], []
#for (x, y) in allcomb:
  #alis.append("{}_on_{}.gff3".format(x,y))
  # err_alis.append("{}_on_{}.err".format(x,y))
  # mikados.append("{}_on_{}.compare.log".format(x,y))
  #cpares.append("{}_on_{}.compare.stats".format(x,y))


rule all:
    input:
        chroms = chroms,
        chromset = chromset,
        stringy = stringy,
        CHROMOSOMES = CHROMOSOMES,
        SUBGENOMES = SUBGENOMES,
        fasta = expand("chr{nom}.fasta", nom=chroms),
        index = expand("Index/chr{nom}/chr{nom}.salcpchilddc", nom=chroms),
        #err = expand("{nom}.err", nom=chroms),
        #start = expand("{nom}.log", nom=chroms),
        #ref_gff3 = expand("{nom}.ref_gff3.log", nom=chroms),
        gff3 = expand("chr{nom}.reference.gff3", nom=chroms),
        cdna = expand("chr{nom}.cdnas.fasta", nom=chroms),
        cds = expand("chr{nom}.cds.fasta", nom=chroms),
        prot = expand("chr{nom}.proteins.fasta", nom=chroms),
        ali = expand("{chrom}{sub_query}_on_{chrom}{sub_target}.gff3", chrom=CHROMOSOMES, sub_query=SUBGENOMES, sub_target=SUBGENOMES),
        pare = expand("{chrom}{sub_query}_on_{chrom}{sub_target}.compare", chrom=CHROMOSOMES, sub_query=SUBGENOMES, sub_target=SUBGENOMES),
        triplets_hc = expand("triplets_HC_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        triplets_lc = expand("triplets_LC_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        triplets_all = expand("triplets_all_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        triplets_hc_eq = expand("triplets_HC_eq_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        triplets_lc_eq = expand("triplets_LC_eq_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        triplets_all_eq = expand("triplets_all_eq_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        crossc_hc = expand("crosscheck_HC_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        crossc_lc = expand("crosscheck_LC_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        crossc_all = expand("crosscheck_all_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        crossc_hc_eq = expand("crosscheck_HC_eq_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        crossc_all_eq = expand("crosscheck_all_eq_{chrom}A-{chrom}B-{chrom}D.tsv", chrom=CHROMOSOMES),
        tripin = expand("{chrom}A {chrom}b {chrom}D", chrom=CHROMOSOMES)
    output: touch("all.done")

##Retrieve & index the chromosomes and extract GFF3S
rule start:
    input: "."
    #    genome="../references/161010_Chinese_Spring_v1.0_pseudomolecules.fasta",
    #    annot="../references/iwgsc_refseqv1.0_UTR_2017May05.gff3"
    output: "."
    #    fasta="chr{nom}.fasta",
    #    index="Index/chr{nom}/chr{nom}.salcpchilddc"
    params:
        chrom=lambda wildcards: "chr%s" % wildcards.nom
    log: "{nom}.log"
    shell: """set +u && source samtools-1.3_dm && source gmap-20170508 && samtools faidx {input.genome} chr2A > {output.fasta} && mkdir -p Index && gmap_build -D $(pwd)/Index -d {params.chroms} {output.fasta} > {log} 2> {log} && set -u"""


rule ref_gff3:
    input:
        annot="../references/iwgsc_refseqv1.0_UTR_2017May05.gff3"
    output:
        gff3="chr{nom}.reference.gff3"
    log:
        "{nom}.ref_gff3.log"
    shell: """egrep \"^(#|{input.chroms})\" {input.annot} | uniq > {output.gff3}"""

##Retrieve FASTAs
rule fastas:
    input:
        rules.start.output.fasta,
        rules.ref_gff3.output.gff3
    output:
        cdna="chr{nom}.cdnas.fasta",
        cds="chr{nom}.cds.fasta",
        prot="chr{nom}.proteins.fasta"
    shell:
        "set +u && cufflinks-2.2.2_beta_20150904 && gffread -g {input.fasta} -w {output.cdna} -x {output.cds} -y {output.prot} {input.gff3} && set -u"


##Perform Alignments
rule gmap:
    input:
        cds=rules.fastas.output.cds
    output:
        ali="{chrom}{sub_query}_on_{chrom}{sub_target}.gff3"
    shell: "set +u && source gmap-20170508 && gmap -F -D $(pwd)/Index -d {input.chroms} --max-intronlength-middle=71000 --max-intronlength-ends=71000 --no-chimeras -t 10 -f gff3_gene -n 1 {input.cds} 2> {output.err_ali} > {output.ali} && set -u"


###Get initial pairwise comparison files (*.compare.*)
rule mikado:
    input:
        ali=rules.gmap.output.ali,
        gff3 = rules.fastas.output.gff3
    output:
        cpare="{chrom}{sub_query}_on_{chrom}{sub_target}.compare"
    log: "{}_on_{}.compare.log"
    shell: "set +u && source mikado-1.0.1 &&mikado compare -r {input.gff3} -p {input.ali} --log {log} -o {output.cpare} && set -u"

rule triplets:
    input:
        tripin="{chrom}A {chrom}b {chrom}D"
    output:
        triplets_hc = "triplets_HC_{chrom}A-{chrom}B-{chrom}D.tsv",
        triplets_lc = "triplets_LC_{chrom}A-{chrom}B-{chrom}D.tsv",
        triplets_all = "triplets_all_{chrom}A-{chrom}B-{chrom}D.tsv",
        triplets_hc_eq = "triplets_HC_eq_{chrom}A-{chrom}B-{chrom}D.tsv",
        triplets_lc_eq = "triplets_LC_eq_{chrom}A-{chrom}B-{chrom}D.tsv",
        triplets_all_eq = "triplets_all_eq_{chrom}A-{chrom}B-{chrom}D.tsv",
        crossc_hc = "crosscheck_HC_{chrom}A-{chrom}B-{chrom}D.tsv",
        crossc_lc = "crosscheck_LC_{chrom}A-{chrom}B-{chrom}D.tsv",
        crossc_all = "crosscheck_all_{chrom}A-{chrom}B-{chrom}D.tsv",
        crossc_hc_eq = "crosscheck_HC_eq_{chrom}A-{chrom}B-{chrom}D.tsv",
        crossc_all_eq = "crosscheck_all_eq_{chrom}A-{chrom}B-{chrom}D.tsv"
    shell:
        "python3 main.py {input.tripin}"
        ###| for i in {1..7}; do echo $i\A $i\B $i\D; done && set u-"