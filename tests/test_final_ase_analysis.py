"""Scientific invariants for the final background ASE analysis."""
import importlib.util
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pysam

SCRIPTS = Path(__file__).resolve().parents[1]/"pipeline/scripts"
sys.path.insert(0, str(SCRIPTS))
import calibrate_ase_shared_dispersion as cal
from final_ase_postprocess import prepare_shared


def test_dispersion_moments_recover_overdispersion_with_unequal_depths():
    rng = np.random.default_rng(23)
    true_rho = 0.06
    mu = rng.uniform(0.05, 0.95, (10000,1))
    n = rng.integers(10, 500, (10000,3))
    c = 1/true_rho-1
    y = rng.binomial(n, rng.beta(mu*c, (1-mu)*c, size=n.shape))
    moments = np.array([cal.dispersion_moments(a,b) for a,b in zip(y,n)])
    estimate = moments[:,0].sum()/moments[:,1].sum()
    assert abs(estimate-true_rho) < 0.005


def test_simulation_is_reproducible_nonzero_and_allele_symmetric():
    y, n, rho = np.array([90,88,92]), np.array([100,100,100]), np.repeat(0.02,3)
    result = cal.simulate_pvalue(y,n,rho,None,9999,42)
    assert result == cal.simulate_pvalue(y,n,rho,None,9999,42)
    assert result[0] >= 1/10000
    assert result[0] < 0.01
    assert np.isclose(cal.statistic(y,n,rho),cal.statistic(n-y,n,rho))
    balanced = cal.simulate_pvalue(np.array([50,49,51]),n,rho,None,9999,42)
    assert balanced[0] == 1 and balanced[-1]


def test_shared_model_null_error_is_not_inflated():
    rng = np.random.default_rng(58)
    n, rho = np.array([40,150,500]), np.repeat(0.05,3)
    c = 1/rho-1
    calls = 0
    for i in range(250):
        y = rng.binomial(n,rng.beta(0.5*c,0.5*c))
        calls += cal.simulate_pvalue(y,n,rho,None,1999,i)[0] < 0.05
    assert calls/250 <= 0.09


def test_differential_statistic_is_group_symmetric():
    y, n = np.array([80,79,81,20,21,19]), np.repeat(100,6)
    rho = np.array([0.03]*3+[0.05]*3)
    first = cal.statistic(y,n,rho,3)
    assert np.isclose(first,cal.statistic(np.r_[y[3:],y[:3]],n,np.r_[rho[3:],rho[:3]],3))
    assert cal.simulate_pvalue(y,n,rho,3,9999,91)[0] < 0.01


def test_ambiguous_exonic_snp_is_excluded_before_fragment_counting(tmp_path):
    fixture = Path(__file__).parent
    bam = tmp_path/"input.bam"
    with pysam.AlignmentFile(str(fixture/"gene_counter_fixture.sam"),"r") as source:
        with pysam.AlignmentFile(str(bam),"wb",template=source) as out:
            for read in source:
                out.write(read)
    pysam.index(str(bam))
    vcf = tmp_path/"snps.vcf.gz"
    pysam.tabix_compress(str(fixture/"gene_counter_fixture.vcf"),str(vcf))
    pysam.tabix_index(str(vcf),preset="vcf")
    gtf = tmp_path/"genes.gtf"
    gtf.write_text((fixture/"gene_counter_fixture.gtf").read_text()+
        'chr1\ttest\texon\t99\t101\t.\t+\t.\tgene_id "gene2"; transcript_id "tx2";\n')
    genes, sites = tmp_path/"genes.tsv",tmp_path/"sites.tsv"
    subprocess.run([sys.executable,str(SCRIPTS/"count_gene_alleles.py"),"--bam",str(bam),"--vcf",str(vcf),
        "--gtf",str(gtf),"--output",str(genes),"--site-output",str(sites),"--exclude-multigene-snps"],check=True)
    g,s = pd.read_csv(genes,sep="\t"),pd.read_csv(sites,sep="\t")
    assert s.position.tolist() == [110]
    assert g.gene_id.tolist() == ["gene1"]
    assert g.maternal_fragments.tolist() == [1]
    assert g.paternal_fragments.tolist() == [2]
    assert g.ambiguous_fragments.tolist() == [0]


def test_common_sites_require_identical_alleles_orientation_and_support(tmp_path):
    (tmp_path/"reference").mkdir()
    (tmp_path/"variants/informative").mkdir(parents=True)
    (tmp_path/"reference/ChineseLong_gene_v3.gtf").write_text(
        'chr1\ttest\texon\t1\t1000\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
    header = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'+
        ''.join(f'##INFO=<ID={key},Number=1,Type=Integer,Description="test">\n' for key in ["MAT_ALLELE","PAT_ALLELE","MAT_N","PAT_N"])+
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    for cross in ("CD","CS"):
        rows=[]
        for pos in (10,20,30,40):
            mat,pat = (1,0) if pos == 20 and cross == "CS" else (0,1)
            alt = "T" if pos == 30 and cross == "CS" else "G"
            support=3 if pos == 40 else 4
            rows.append(f'chr1\t{pos}\t.\tA\t{alt}\t100\tPASS\tMAT_ALLELE={mat};PAT_ALLELE={pat};MAT_N={support};PAT_N={support}\n')
        src=tmp_path/f"{cross}.vcf"
        src.write_text(header+''.join(rows))
        dest=tmp_path/f"variants/informative/{cross}.informative_snps.vcf.gz"
        pysam.tabix_compress(str(src),str(dest))
        pysam.tabix_index(str(dest),preset="vcf")
    prepare_shared(tmp_path)
    for support,positions in [(3,[10,40]),(4,[10])]:
        with pysam.VariantFile(tmp_path/f"common_support{support}/shared_informative_snps.vcf.gz") as vcf:
            assert [r.pos for r in vcf] == positions
