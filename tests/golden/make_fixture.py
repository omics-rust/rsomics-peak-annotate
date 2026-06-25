#!/usr/bin/env python3
"""Deterministic synthetic GTF + peaks exercising every ChIPseeker annotatePeak
category on both strands. CDS phases are 0 (valid 0-2) and CDS sit inside exons
so makeTxDbFromGFF keeps every transcript and derives 5'/3' UTRs from exon-vs-CDS."""
CHROM = "chr1"
gtf = []
def row(feat, s, e, strand, attrs, frame="."):
    gtf.append(f'{CHROM}\thand\t{feat}\t{s}\t{e}\t.\t{strand}\t{frame}\t{attrs}')

def gene(gid, tid, strand, gs, ge, exons, cds):
    row("gene", gs, ge, strand, f'gene_id "{gid}"; gene_name "{gid}";')
    row("transcript", gs, ge, strand, f'gene_id "{gid}"; transcript_id "{tid}";')
    for i,(s,e) in enumerate(exons,1):
        row("exon", s, e, strand, f'gene_id "{gid}"; transcript_id "{tid}"; exon_number "{i}";')
    for (s,e) in cds:
        row("CDS", s, e, strand, f'gene_id "{gid}"; transcript_id "{tid}";', frame="0")

# geneA: + strand, TSS 10000; 5'UTR 10000-10199, 3'UTR 19801-20000
gene("GA","TA","+",10000,20000,
     [(10000,10500),(12000,12500),(19000,20000)],
     [(10200,10500),(12000,12500),(19000,19800)])
# geneB: - strand, TSS 50000; 5'UTR 49801-50000, 3'UTR 40000-40199
gene("GB","TB","-",40000,50000,
     [(49500,50000),(45000,45500),(40000,41000)],
     [(49500,49800),(45000,45500),(40200,41000)])
# geneC: + strand single-exon, TSS 70000
gene("GC","TC","+",70000,75000,[(70000,75000)],[(70100,74900)])


# geneD: + strand, two transcripts with distinct TSS (multi-transcript probe)
row("gene",80000,90000,"+",'gene_id "GD"; gene_name "GD";')
row("transcript",80000,90000,"+",'gene_id "GD"; transcript_id "TD1";')
for i,(s,e) in enumerate([(80000,81000),(85000,86000),(89000,90000)],1):
    row("exon",s,e,"+",f'gene_id "GD"; transcript_id "TD1"; exon_number "{i}";')
for (s,e) in [(80200,81000),(85000,86000),(89000,89800)]:
    row("CDS",s,e,"+",'gene_id "GD"; transcript_id "TD1";',frame="0")
row("transcript",83000,90000,"+",'gene_id "GD"; transcript_id "TD2";')
for i,(s,e) in enumerate([(83000,84000),(89000,90000)],1):
    row("exon",s,e,"+",f'gene_id "GD"; transcript_id "TD2"; exon_number "{i}";')
for (s,e) in [(83200,84000),(89000,89800)]:
    row("CDS",s,e,"+",'gene_id "GD"; transcript_id "TD2";',frame="0")

# geneH: + strand, two isoforms sharing a far exon (320000-321000) with identical
# coordinates but distinct TSS — tests the category/gene decoupling at coord-equal exons
row("gene",300000,321000,"+",'gene_id "GH"; gene_name "GH";')
row("transcript",300000,321000,"+",'gene_id "GH"; transcript_id "TH1";')
for i,(s,e) in enumerate([(300000,301000),(320000,321000)],1):
    row("exon",s,e,"+",f'gene_id "GH"; transcript_id "TH1"; exon_number "{i}";')
for (s,e) in [(300200,301000),(320000,320800)]:
    row("CDS",s,e,"+",'gene_id "GH"; transcript_id "TH1";',frame="0")
row("transcript",310000,321000,"+",'gene_id "GH"; transcript_id "TH2";')
for i,(s,e) in enumerate([(310000,311000),(320000,321000)],1):
    row("exon",s,e,"+",f'gene_id "GH"; transcript_id "TH2"; exon_number "{i}";')
for (s,e) in [(310200,311000),(320000,320800)]:
    row("CDS",s,e,"+",'gene_id "GH"; transcript_id "TH2";',frame="0")
# geneI cluster: two close TSS, for a wide peak spanning both
gene("GI1","TI1","+",400000,405000,[(400000,405000)],[(400200,404800)])
gene("GI2","TI2","+",402000,407000,[(402000,407000)],[(402200,406800)])

open("anno.gtf","w").write("\n".join(gtf)+"\n")

peaks = [
    ("p1",  9400, 9600),   # GA promoter (TSS 10000)
    ("p2",  8200, 8400),   # GA promoter 1-2kb
    ("p3",  7100, 7300),   # GA promoter 2-3kb
    ("p4", 12100,12300),   # GA exon2
    ("p5", 13000,13200),   # GA intron
    ("p6", 19850,19950),   # GA 3'UTR
    ("p7", 21500,21700),   # GA downstream
    ("p8", 30000,30200),   # distal intergenic
    ("p9", 50300,50500),   # GB promoter (- strand TSS 50000)
    ("pA", 44000,44200),   # GB intron
    ("pB", 19200,19400),   # GA exon3 far from TSS (>3kb) -> Exon
    ("pC", 20050,20150),   # GA just downstream of end 20000 (50bp) -> Downstream?
    ("pD",  9960,10040),   # spans GA TSS 10000 -> distance 0
    ("pE", 45100,45300),   # GB exon2 far from TSS -> Exon (- strand rank)
    ("pF", 82500,82700),   # between GD TD1(80000)/TD2(83000) TSS, nearer TD2
    ("pG", 80100,80300),   # GD TD1 promoter
    ("pH", 320400,320600), # GH shared coord-identical exon, far from both TSS -> Exon
    ("pI", 401000,403000), # 2kb peak spanning GI1(400000)/GI2(402000) TSS -> Promoter dist 0
]
with open("peaks.bed","w") as f:
    for n,s,e in peaks:
        f.write(f"{CHROM}\t{s}\t{e}\t{n}\t0\t.\n")
print(f"wrote anno.gtf ({len(gtf)} lines), peaks.bed ({len(peaks)} peaks)")
