import os
import sys

KATUALI_HOME = os.path.split(sys.exec_prefix)[0]
configfile: os.path.join(KATUALI_HOME, "config.yaml")

SCRAPPIE_EXEC = os.path.join(os.path.expanduser(config["SCRAPPIE"]), "build", "scrappie")
SCRAPPIE_JSON_TO_TSV = os.path.join(os.path.expanduser(config["SCRAPPIE"]), "misc", "json_to_tsv.py")
NANOPOLISH_EXEC = os.path.join(os.path.expanduser(config["NP"]), "nanopolish")
IN_POMOXIS = os.path.expanduser(config["IN_POMOXIS"])
IN_MEDAKA = os.path.expanduser(config["IN_MEDAKA"])
IN_MIYAGI = os.path.expanduser(config["IN_MIYAGI"])

# NOTE on virtual environments
# Snakemake uses bash strict mode, virtualenv violates bash strict mode.
# https://snakemake.readthedocs.io/en/stable/project_info/faq.html#my-shell-command-fails-with-with-errors-about-an-unbound-variable-what-s-wrong
# so any activation commands must be wrapped as such:
#      set +u; {config[IN_MEDAKA]}; set -u; 

rule help:
    run:
        print(open(os.path.join(KATUALI_HOME, "README.md")).read())

def get_contig_opt(wildcards):
    if wildcards.contig == "all_contigs":
        contig_opt = ""
    else:
        contig_opt = "-r {}".format(wildcards.contig)
    return contig_opt

# support either processing all contigs or a single contig
ALL_CONTIGS = "all_contigs"
if "contig" not in config:
    config["contig"] = ALL_CONTIGS

CONTIGOPTION = "" if config["contig"] == ALL_CONTIGS else "-r {}".format(config["contig"])

rule basecall_scrappie:
    input:
        scrappie = SCRAPPIE_EXEC,
        fast5 = config["READS"]
    output:
        fasta = "analysis/basecall/scrappie{suffix,[^/]*}/basecalls.fasta",
    log:
        "analysis/basecall/scrappie{suffix,[^/]*}/scrappie.log",
    params:
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]) 
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        echo "{input.scrappie} {config[SCRAPPIE_OPTS]} "> {log}
        find -L {input.fast5} -name '*.fast5' | xargs {input.scrappie} {config[SCRAPPIE_OPTS]} --threads {threads} >{output.fasta} 2>> {log}
        """
    

# TODO rule basecall_guppy:

rule scrappie_summary:
    input:
        json_to_tsv = SCRAPPIE_JSON_TO_TSV,
        fasta = "{dir}/basecalls.fasta",
    output:
        summary = "{dir}/sequencing_summary.txt"
    params:
        sge = "m_mem_free=1G,gpu=0"
    shell:
    	"cat {input.fasta} | grep '^>' | cut -d ' ' -f 2- | python {input.json_to_tsv} > {output.summary} &&"
    	# most tools expect read_id, not uuid
	    "sed -i '1 s/uuid/read_id/' {output.summary}"

rule align_to_ref:
    input:
        venv = IN_POMOXIS,
        basecalls = "{bc_dir}/basecalls.fasta",
        ref = config["REFERENCE"]
    params:
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]) 
    output:
        bam = "{bc_dir}/align/calls2ref.bam"
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
    	mini_align -i {input.basecalls} -r {input.ref} -p {wildcards.bc_dir}/align/calls2ref -P -t {config[THREADS_PER_JOB]}
        """

rule align_to_draft:
    input:
        venv = IN_POMOXIS,
        basecalls = "{dir}/{subdir}/basecalls.fasta",
        draft = "{dir}/consensus.fasta"
    output:
        bam = "{dir}/{subdir}/calls2draft.bam"
    params:
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]) 
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
    	mini_align -i {input.basecalls} -r {input.draft} -p {wildcards.dir}/{wildcards.subdir}/calls2draft -P -t {config[THREADS_PER_JOB]}
        """

rule get_depth:
    input:
        venv = IN_POMOXIS,
        bam = "{dir}/calls2ref.bam"
    output:
        depth = directory("{dir}/depth")
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
    	mkdir -p {output.depth} && coverage_from_bam {input.bam} -s 1000 -p {output.depth}/ &>{output.depth}/depth.log
        """

rule subsample_bam:
    input:
        venv = IN_POMOXIS,
        bam = "{dir}/calls2ref.bam"
    output:
        fasta = "{dir}/{contig}/{depth}X/basecalls.fasta",
    log:
        "{dir}/{contig}/{depth}X/subsample.log"
    params:
        contig_opt = get_contig_opt,
        prefix = lambda w: "{dir}/{contig}/{depth}X/sub_sample".format(**dict(w)),
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]) 
    threads: config["THREADS_PER_JOB"]
    shell:
	    """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        subsample_bam {input.bam} {wildcards.depth} {params[contig_opt]} -o {params[prefix]} -t {threads} &>{log};
        sleep 5;
        for i in {params[prefix]}*.bam; do samtools fasta $i ; done > {output.fasta}
        """

rule ref_guided_racon:
    input:
        venv = IN_POMOXIS,
        basecalls = "{dir}/basecalls.fasta",
        ref = config["REFERENCE"]
    output:
        consensus = "{dir}/ref_guided_racon{suffix,[^/]*}/consensus.fasta",
        basecalls = "{dir}/ref_guided_racon{suffix,[^/]*}/basecalls.fasta"
    log:
        "{dir}/ref_guided_racon{suffix}.log"
    threads: config["THREADS_PER_JOB"]
    params:
        racon_dir = lambda w: "{dir}/ref_guided_racon{suffix}".format(**dict(w)),
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]) 
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        # snakemake will create the output dir, mini_assemble will fail if it exists..
        rm -r {params[racon_dir]} && 
        mini_assemble -i {input.basecalls} -r {input.ref} -o {params[racon_dir]} -t {threads} -p assm {config[MINI_ASSEMBLE_OPTS]} &&
        # rename output
        mv {params[racon_dir]}/assm_final.fa {output.consensus}
        # keep a link of basecalls with the consensus
        ln -s $PWD/{input.basecalls} $PWD/{params[racon_dir]}/basecalls.fasta &&
        # sync timestamps, without following basecalls link (otherwise consensus will be older than basecalls)
        touch --no-dereference $PWD/{params[racon_dir]}/*
        """

rule medaka_consensus:
    input:
        venv = IN_MEDAKA,
        draft = "{dir}/consensus.fasta",
        basecalls = "{dir}/basecalls.fasta"
    output:
        consensus = "{dir}/medaka/consensus.fasta"
    log:
        "{dir}/medaka.log"
    threads: config["THREADS_PER_JOB"]
    params:
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]) 
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        # snakemake will create the output dir, medaka_consensus will fail if it exists..
        rm -r {wildcards.dir}/medaka &&
        medaka_consensus -i {input.basecalls} -d {input.draft} -o {wildcards.dir}/medaka -t {threads} &> {log}
        # keep a link of basecalls with the consensus
        ln -s $PWD/{input.basecalls} $PWD/{wildcards.dir}/medaka/basecalls.fasta
        """

rule nanopolish_index:
    input:
        nanopolish = NANOPOLISH_EXEC,
        fast5 = config["READS"],
        summary = "analysis/basecall/{basecaller}/sequencing_summary.txt",
        basecalls = "analysis/basecall/{basecaller}/align/{subdir}/basecalls.fasta",
    output:
        "analysis/basecall/{basecaller}/align/{subdir}/nanopolish/basecalls.fasta",
        "analysis/basecall/{basecaller}/align/{subdir}/nanopolish/basecalls.fasta.index",
        "analysis/basecall/{basecaller}/align/{subdir}/nanopolish/basecalls.fasta.index.gzi",
        "analysis/basecall/{basecaller}/align/{subdir}/nanopolish/basecalls.fasta.index.fai",
        "analysis/basecall/{basecaller}/align/{subdir}/nanopolish/basecalls.fasta.index.readdb",
    log:
        "analysis/basecall/{basecaller}/align/{subdir}/nanopolish/nanopolish_index.log"
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
    	# create indices then synchronise time stamps
	    """
        ln -s $PWD/{input.basecalls} $PWD/{output[0]} && sleep 1 && 
        {input.nanopolish} index -d {input.fast5} -s {input.summary} {output[0]} &> {log} && sleep 5 &&
    	touch --no-dereference {output[0]}*
        """

rule fast5_list:
    input:
        fast5 = config["READS"],
    output:
        filelist = os.path.join(config["READS"], "reads.txt")
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        "find -L `readlink -f {config[READS]}` -name '*.fast5' > {output.filelist}"

rule miyagi_index:
    input:
        filelist = os.path.join(config["READS"], "reads.txt"),
        summary = "analysis/basecall/{basecaller}/sequencing_summary.txt",
        basecalls = "analysis/basecall/{basecaller}/align/{subdir}/basecalls.fasta",
    output:
        index = "analysis/basecall/{basecaller}/align/{subdir}/miyagi/miyagi.index",
        basecalls = "analysis/basecall/{basecaller}/align/{subdir}/miyagi/basecalls.fasta",
    params:
        sge = "m_mem_free=1G,gpu=0" 
    run:
        # create miyagi index linking read_id to file path
        import os
        import pandas as pd
        summary = pd.read_table(input.summary)
        paths = pd.read_table(input.filelist, names=["path"])
        paths["filename"] = paths["path"].map(os.path.basename)
        index = pd.merge(summary, paths, on="filename")
        cols = ["read_id", "path"]
        index[cols].to_csv(output.index, index=False, sep="\t")
        # create symlink to basecalls
        os.symlink(os.path.abspath(input.basecalls), os.path.abspath(output.basecalls))

rule nanopolish_vcf:
    input:
        nanopolish = NANOPOLISH_EXEC,
        draft = "analysis/basecall/{basecaller}/align/{subdir}/consensus.fasta",
        basecalls = "analysis/basecall/{basecaller}/align/{subdir}/nanopolish/basecalls.fasta",
        bam = "analysis/basecall/{basecaller}/align/{subdir}/nanopolish/calls2draft.bam",
        index = "analysis/basecall/{basecaller}/align/{subdir}/nanopolish/basecalls.fasta.index",
        index_gzi = "analysis/basecall/{basecaller}/align/{subdir}/nanopolish/basecalls.fasta.index.gzi",
        index_fai = "analysis/basecall/{basecaller}/align/{subdir}/nanopolish/basecalls.fasta.index.fai",
        index_readdb = "analysis/basecall/{basecaller}/align/{subdir}/nanopolish/basecalls.fasta.index.readdb",
    output:
        vcf = "analysis/basecall/{basecaller}/align/{subdir}/nanopolish/regions/{region}.vcf",
    log:
        "analysis/basecall/{basecaller}/align/{subdir}/nanopolish/regions/{region}.vcf.log"
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
	    """
        {input.nanopolish} variants --consensus -o {output.vcf} -w {wildcards.region} -r {input.basecalls} -b {input.bam} -g {input.draft} -t 1 {config[NP_OPTS]} &> {log}
        """

rule nanopolish_regions:
    input:
        venv = IN_POMOXIS, # use pomoxis python as this has all requirements of the script
        make_range = os.path.join(os.path.expanduser(config["NP"]), "scripts", "nanopolish_makerange.py"),
        draft = "{dir}/consensus.fasta",
    output:
        # Number of regions is unknown ahead of time, so use dynamic keyword to delay evaluation
        regions = dynamic("{dir}/nanopolish/regions/{region}.region"),
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        python {input.make_range} {input.draft} > {wildcards.dir}/nanopolish/regions.list &&
        for i in `cat {wildcards.dir}/nanopolish/regions.list`; do touch {wildcards.dir}/nanopolish/regions/$i.region; done
        """

rule nanopolish:
    input:
        nanopolish = NANOPOLISH_EXEC,
        draft = "{dir}/consensus.fasta",
        # Number of regions is unknown ahead of time, so use dynamic keyword to delay evaluation
        regions = dynamic("{dir}/nanopolish/regions/{region}.region"),
        vcfs = dynamic("{dir}/nanopolish/regions/{region}.vcf"),
    output:
        consensus = "{dir}/nanopolish/consensus.fasta",
    log:
        "{dir}/nanopolish/vcf2fasta.log"
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        "{input.nanopolish} vcf2fasta -g {input.draft} {input.vcfs} > {output.consensus} 2> {log}"

rule miyagi_hp_vcf:
    input:
        venv = IN_MIYAGI,
        draft = "{dir}/consensus.fasta",
    output:
        vcf = "{dir}/miyagi/regions/{region}.hp.vcf"
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        miyagi_create_homop_vcf {input.draft} {wildcards.region} {output.vcf} 
        """

rule miyagi_score_reads:
    input:
        venv = IN_MIYAGI,
        draft = "{dir}/consensus.fasta",
        #region = "{dir}/miyagi/regions/{region}.region",
        bam = "{dir}/miyagi/calls2draft.bam",
        index = "{dir}/miyagi/miyagi.index",
        vcf = "{dir}/miyagi/regions/{region}.hp.vcf"
    output:
        hdf = "{dir}/miyagi/regions/parts/{region}_part_{n}_of_{N}_scores_per_read.hdf",
    log:
        "{dir}/miyagi/regions/parts/{region}_part_{n}_of_{N}_scores_per_read.log"
    params:
        prefix = "{dir}/miyagi/regions/parts/{region}_part_{n}_of_{N}",  # miyagi will append _scores_per_read to this
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        export OPENBLAS_NUM_THREADS=1;
        miyagi_score_homop_reads {wildcards.region} {input.bam} {input.draft} {input.index} {input.vcf} {params.prefix} {config[MIYAGI_SCORE_READS_OPTS]} --workers 1 --batch {wildcards.n} {wildcards.N} &> {log}
        """

def wildcards_2_score_targets(wildcards):
    """Generate a list of part_{n}_of_{N}_scores_per_read.hdf targets for a single region"""
    region_bounds = [int(d) for d in wildcards.region.split(":")[-1].split("-")]
    region_len = region_bounds[1] - region_bounds[0]
    n_batches = max(1, round(region_len * (10 ** -6) * config["MIYAGI_BATCHES_PER_REF_MB"]))
    template = "{}/miyagi/regions/parts/{}_part_{}_of_{}_scores_per_read.hdf"
    return [template.format(wildcards.dir, wildcards.region, n, n_batches) for n in range(n_batches)]

def wildcards_2_sorted_score_targets(wildcards):
    """Generate a list of {region}_sorted_scores_per_read.hdf for all regions/contigs in {wildcards[dir]}/consensus.fasta"""
    import pysam
    draft = os.path.join(wildcards.dir, "consensus.fasta")
    dir = os.path.join(wildcards.dir, "miyagi", "regions")
    with pysam.FastaFile(draft) as fasta:
        template = "{}/{}:0-{}_sorted_scores_per_read.hdf"
        targets = [template.format(dir, ref, l) for ref, l in zip(fasta.references, fasta.lengths)]
    return targets

rule miyagi_concat_region:
    input:
        wildcards_2_score_targets,
        venv = IN_MIYAGI,
    output:
        hdf = "{dir}/miyagi/regions/{region}_sorted_scores_per_read.hdf",
    log:
        "{dir}/miyagi/regions/{region}_sorted_scores_per_read.log"
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        miyagi_concat_scores --output {output.hdf} {wildcards.dir}/miyagi/regions/parts/{wildcards.region}_part*.hdf --debug &> {log}
        """

rule miyagi_apply_model:
    input:
        venv = IN_MIYAGI,
        draft = "{dir}/consensus.fasta",
        hdfs = wildcards_2_sorted_score_targets,
    output:
        vcf = "{dir}/miyagi/corrections.vcf",
    log:
        "{dir}/miyagi/apply_model.log",
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        miyagi_apply_model {input.draft} {output.vcf} --score_input {input.hdfs} {config[MIYAGI_APPLY_MODEL_OPTS]} --debug &> {log}
        """

rule miyagi_apply_corrections:
    input:
        venv = IN_MIYAGI,
        draft = "{dir}/consensus.fasta",
        vcf = "{dir}/miyagi/corrections.vcf",
    output:
        consensus = "{dir}/miyagi/consensus.fasta",
    log:
        "{dir}/miyagi/apply_corrections.log",
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        miyagi_apply_corrections {input.draft} {input.vcf} {output.consensus} &> {log}
        """
