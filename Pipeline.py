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
    #print(chromss)
    for x, y in itertools.permutations(chromss, 2):
        comb[chrom].append((x, y))

allcomb = [comb[x] for x in comb]

TRIPLETS = None

rule all:
    input:
        #chroms = chroms,
        #chromset = chromset,
        #stringy = stringy,
        #CHROMOSOMES = CHROMOSOMES,
        #SUBGENOMES = SUBGENOMES,
        fasta = expand("chr{chrom}{sub_target}.fasta", chrom=CHROMOSOMES, sub_target=SUBGENOMES),
        index = expand("Index/chr{chrom}{sub_target}/chr{chrom}{sub_target}.salcpchilddc", chrom=CHROMOSOMES, sub_target=SUBGENOMES),
        #err = expand("{chrom}{sub_target}.err", chrom=CHROMOSOMES, sub_target=SUBGENOMES),
        #start = expand("{chrom}{sub_target}.log", chrom=CHROMOSOMES, sub_target=SUBGENOMES),
        #ref_gff3 = expand("{chrom}{sub_target}.ref_gff3.log", chrom=CHROMOSOMES, sub_target=SUBGENOMES),
        gff3 = expand("chr{chrom}{sub_target}.reference.gff3", chrom=CHROMOSOMES, sub_target=SUBGENOMES),
        cdna = expand("chr{chrom}{sub_target}.cdnas.fasta", chrom=CHROMOSOMES, sub_target=SUBGENOMES),
        cds = expand("chr{chrom}{sub_target}.cds.fasta", chrom=CHROMOSOMES, sub_target=SUBGENOMES),
        prot = expand("chr{chrom}{sub_target}.proteins.fasta", chrom=CHROMOSOMES, sub_target=SUBGENOMES),
        ali = expand("{chrom}{sub_query}_on_{chrom}{sub_target}.gff3", chrom=CHROMOSOMES, sub_query=SUBGENOMES, sub_target=SUBGENOMES),
        cpare = expand("{chrom}{sub_query}_on_{chrom}{sub_target}.compare", chrom=CHROMOSOMES, sub_query=SUBGENOMES, sub_target=SUBGENOMES),
        stats_cm=expand("{chrom}{sub_query}_on_{chrom}{sub_target}.stats.tsv",chrom=CHROMOSOMES, sub_query=SUBGENOMES, sub_target=SUBGENOMES),
        stats_cx=expand("{chrom}{sub_query}_on_{chrom}{sub_target}.stats.txt",chrom=CHROMOSOMES, sub_query=SUBGENOMES, sub_target=SUBGENOMES),
        stats_rf=expand("{chrom}{sub_target}.reference.stats.tsv",chrom=CHROMOSOMES, sub_target=SUBGENOMES),
        stats_rx=expand("{chrom}{sub_target}.reference.stats.txt",chrom=CHROMOSOMES, sub_target=SUBGENOMES),
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
        #tripin = expand("{chrom}A {chrom}B {chrom}D", chrom=CHROMOSOMES)
    output: touch("all.done")

##Retrieve & index the chromosomes and extract GFF3S
rule start:
    input: genome="../references/161010_Chinese_Spring_v1.0_pseudomolecules.fasta"
    output:
        fasta="chr{chrom}{sub_target}.fasta",
    params:
        chroms=lambda wildcards: "chr%s" % wildcards.chrom,
	subg=lambda wildcards: wildcards.sub_target
    log: "{chrom}{sub_target}.log"
    shell: """set +u && source samtools-1.3_dm && source gmap-20170508 && samtools faidx {input.genome} {params.chroms}{params.subg} > {output.fasta} && mkdir -p Index && gmap_build -D $(pwd)/Index -d {params.chroms}{params.subg} {output.fasta} > {log} 2> {log} && set -u"""

rule infas:
    input: fasta = rules.start.output.fasta
    output:
        index="Index/chr{chrom}{sub_target}/chr{chrom}{sub_target}.salcpchilddc"
    params:
        chroms=lambda wildcards: "chr%s" % wildcards.chrom,
        subg=lambda wildcards: wildcards.sub_target
    log: "{chrom}{sub_target}.log"
    shell: """set +u && source samtools-1.3_dm && source gmap-20170508 && mkdir -p Index && gmap_build -D $(pwd)/Index -d {params.chroms}{params.subg} {input.fasta} > {log} 2> {log} && set -u"""

rule ref_gff3:
    input:
        annot="../references/iwgsc_refseqv1.0_UTR_2017May05.gff3"
    output:
        gff3="chr{chrom}{sub_target}.reference.gff3"
    params:
        chroms=lambda wildcards: "chr%s" % wildcards.chrom,
	subg=lambda wildcards: wildcards.sub_target
    log:
        "{chrom}{sub_target}.ref_gff3.log"
    shell:"egrep \"^(#|{params.chroms}{params.subg})\" {input.annot} | uniq > {output.gff3}"

##Retrieve FASTAs
rule fastas:
    input:
        fasta = rules.start.output.fasta,
        gff3 = rules.ref_gff3.output.gff3
    output:
        cdna="chr{chrom}{sub_target}.cdnas.fa",
        cds="chr{chrom}{sub_target}.cds.fa",
        prot="chr{chrom}{sub_target}.proteins.fa"
    shell:
       "" "set +u && cufflinks-2.2.2_beta_20150904 && gffread -g {input.fasta} -w {output.cdna} -x {output.cds} -y {output.prot} {input.gff3} && set -u"""


##Perform Alignments
rule gmap:
    input:
        cds=rules.fastas.output.cds
    output:
        ali="{chrom}{sub_query}_on_{chrom}{sub_target}.gff3",
    params:
        chroms=lambda wildcards: "%s" % wildcards.chrom,
	subg=lambda wildcards: wildcards.sub_target
    shell: """set +u && source gmap-20170508 && gmap -F -D $(pwd)/Index -d {params.chroms}{params.subg} --max-intronlength-middle=71000 --max-intronlength-ends=71000 --no-chimeras -t 10 -f gff3_gene -n 1 {input.cds} 2> {output.ali} && set -u"""


###Get initial pairwise comparison files (*.compare.*)
rule mikado:
    input:
        ali=rules.gmap.output.ali,
        gff3 = rules.ref_gff3.output.gff3
    output:
        cpare="{chrom}{sub_query}_on_{chrom}{sub_target}.compare"
    log: "{}_on_{}.compare.log"
    shell: """set +u && source mikado-1.0.1 && mikado compare -r {input.gff3} -p {input.ali} --log {log} -o {output.cpare} && set -u"""

rule statsdo:
    input:
        ali = rules.gmap.output.ali
    output:
        stats_cm="{chrom}{sub_query}_on_{chrom}{sub_target}.stats.tsv",
        stats_cx="{chrom}{sub_query}_on_{chrom}{sub_target}.stats.txt"
    shell: "set +u && source mikado-1.0.1 && mikado util stats --tab-stats {output.stats_cm} {input.ali} {output.stats_cx} && set -u"

rule ref_stats:
    input:
       gff3 = rules.ref_gff3.output.gff3
    output:
       stats_rf="{chrom}{sub_target}.reference.stats.tsv",
       stats_rx="{chrom}{sub_target}.reference.stats.txt"
    shell: "set +u && source mikado-1.0.1 && mikado util stats --tab-stats {output.stats_rf} {input.gff3} {output.stats_rx} && set -u"


rule triplets:
    input:
        #tripin="{chrom}A {chrom}B {chrom}D"
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
    params: chroms=lambda wildcards: "%s" % wildcards.chrom
    shell:
        "python3 main.py {params.chroms}A {params.chroms}B {params.chroms}D"
