import os
import sys

KATUALI_HOME = os.path.split(sys.exec_prefix)[0]
configfile: os.path.join(KATUALI_HOME, "config.yaml")

GUPPY_EXEC = os.path.expanduser(config["GUPPY"])
SCRAPPIE_EXEC = os.path.join(os.path.expanduser(config["SCRAPPIE"]), "build", "scrappie")
SCRAPPIE_JSON_TO_TSV = os.path.join(os.path.expanduser(config["SCRAPPIE"]), "misc", "json_to_tsv.py")
NANOPOLISH_EXEC = os.path.join(os.path.expanduser(config["NP"]), "nanopolish")
NANOPOLISH_MAKE_RANGE = os.path.join(os.path.expanduser(config["NP"]), "scripts", "nanopolish_makerange.py")
IN_POMOXIS = os.path.expanduser(config["IN_POMOXIS"])
IN_MEDAKA = os.path.expanduser(config["IN_MEDAKA"])
IN_MIYAGI = os.path.expanduser(config["IN_MIYAGI"])
IN_RAY = os.path.expanduser(config["IN_RAY"])
CANU_EXEC = os.path.expanduser(config["CANU_EXEC"])

# NOTE on virtual environments
# Snakemake uses bash strict mode, virtualenv violates bash strict mode.
# https://snakemake.readthedocs.io/en/stable/project_info/faq.html#my-shell-command-fails-with-with-errors-about-an-unbound-variable-what-s-wrong
# so any activation commands must be wrapped as such:
#      set +u; {config[IN_MEDAKA]}; set -u; 

rule help:
    run:
        print(open(os.path.join(KATUALI_HOME, "README.md")).read())

rule fast_assm_polish:
    input:
        consensus = ancient("basecall/{BASECALLER}/miniasm_racon/medaka/consensus.fasta".format(**config))

rule standard_assm_polish:
    input:
        consensus = ancient("basecall/{BASECALLER}/canu/nanopolish/consensus.fasta".format(**config))

def get_contig_opt(wildcards):
    if wildcards.contig == "all_contigs":
        contig_opt = ""
    else:
        contig_opt = "-r {}".format(wildcards.contig)
    return contig_opt

rule basecall_scrappie:
    input:
        scrappie = ancient(SCRAPPIE_EXEC),
        fast5 = ancient(config["READS"]),
    output:
        fasta = "basecall/scrappie{suffix,[^/]*}/basecalls.fasta",
    log:
        "basecall/scrappie{suffix,[^/]*}/scrappie.log",
    params:
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]) 
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        echo "{input.scrappie} {config[SCRAPPIE_OPTS]} "> {log}
        find -L {input.fast5} -name '*.fast5' | xargs {input.scrappie} {config[SCRAPPIE_OPTS]} --threads {threads} >{output.fasta} 2>> {log}
        """
    
rule basecall_guppy:
    input:
        guppy = ancient(GUPPY_EXEC),
        fast5 = ancient(config["READS"]),
        venv = ancient(IN_POMOXIS),
    output:
        fasta = "basecall/guppy{suffix,[^/]*}/basecalls.fasta",
        summary = "basecall/guppy{suffix,[^/]*}/sequencing_summary.txt",
    log:
        "basecall/guppy{suffix,[^/]*}.log",
    params:
        sge = "m_mem_free=1G,gpu=1 -pe mt {}".format(config["GUPPY_SLOTS"]),
        output_dir = lambda w: "basecall/guppy{suffix}".format(**dict(w))
    shell:
        """
        # snakemake will create the output dir, guppy will fail if it exists..
        rm -r {params[output_dir]}
 
        echo "{input.guppy} {config[GUPPY_OPTS]} " > {log}
        echo "GPU status before" >> {log}
        gpustat >> {log}
         
        if [ -z ${{SGE_HGR_gpu+x}} ]; then 
          echo "SGE_HGR_gpu was not set" >> {log}
          # pick GPU with lowest memory usage
          SGE_HGR_gpu=`gpustat | awk '{{ if(NR==2){{gpu0=$11}}; if(NR==3){{gpu1=$11}} }}END{{if(gpu0 > gpu1){{print "1"}} else{{print "0"}} }}'`
          echo "set GPU to $SGE_HGR_gpu based on memory usage" >> {log}
        fi
        SGE_HGR_gpu="${{SGE_HGR_gpu#cuda}}"

        echo "Runnning on host $HOSTNAME GPU $SGE_HGR_gpu" >> {log}

        {input.guppy} -s {params.output_dir} -r -i {input.fast5} -x cuda:$SGE_HGR_gpu {config[GUPPY_OPTS]} &>> {log}

        echo "gpustat after" >> {log}
        gpustat >> {log}

        # convert fastq to fasta
        sleep 5
        echo "Combining the following fastq files into {output.fasta}" >> {log}
        ls {params[output_dir]}/*.fastq >> {log}
        set +u; {config[SOURCE]} {input.venv}; set -u;
        seqkit fq2fa {params.output_dir}/*.fastq > {output.fasta}
       
        # update time stamp of summary otherwise it will be older than basecalls
        if [ ! -f {output.summary} ]; then
          echo "{output.summary} not found!" >> {log}
        else
          echo "Updating the time stamp of the sequencing summary {output.summary}" >> {log}
          touch {output.summary}  
        fi
        """

rule scrappie_summary:
    input:
        json_to_tsv = ancient(SCRAPPIE_JSON_TO_TSV),
        fasta = ancient("scrappie/basecalls.fasta"),
    output:
        summary = "scrappie/sequencing_summary.txt"
    params:
        sge = "m_mem_free=1G,gpu=0"
    shell:
    	"cat {input.fasta} | grep '^>' | cut -d ' ' -f 2- | python {input.json_to_tsv} > {output.summary} &&"
    	# most tools expect read_id, not uuid
	    "sed -i '1 s/uuid/read_id/' {output.summary}"

rule align_to_ref:
    input:
        venv = ancient(IN_POMOXIS),
        basecalls = ancient("{bc_dir}/basecalls.fasta"),
        ref = ancient(config["REFERENCE"]),
    output:
        bam = "{bc_dir}/align/calls2ref.bam"
    params:
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]),
        prefix = lambda w, output: os.path.splitext(output.bam)[0],
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
    	mini_align -i {input.basecalls} -r {input.ref} -p {params.prefix} -P -t {threads}
        """

rule align_to_draft:
    input:
        venv = ancient(IN_POMOXIS),
        basecalls = ancient("{dir}/{subdir}/basecalls.fasta"),
        draft = ancient("{dir}/consensus.fasta"),
    output:
        bam = "{dir}/{subdir}/calls2draft.bam"
    params:
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]) 
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
    	mini_align -i {input.basecalls} -r {input.draft} -p {wildcards.dir}/{wildcards.subdir}/calls2draft -P -t {threads}
        """

rule assess_consensus:
    input:
        venv = ancient(IN_POMOXIS),
        consensus = ancient("{dir}/consensus.fasta"),
        truth = ancient(config["REFERENCE"]),
    output:
        summ = "{dir}/consensus_to_truth_summ.txt",
        bam = "{dir}/consensus_to_truth.bam",
        stats = "{dir}/consensus_to_truth_stats.txt",
    params:
        prefix = "{dir}/consensus_to_truth",
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]) 
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        assess_assembly -i {input.consensus} -r {input.truth} -p {params.prefix} -t {threads} {config[ASSESS_ASSM_OPTS]} 
        """

rule ray_catalogue:
    input:
        venv = ancient(IN_RAY),
        bam = ancient("{dir}/{prefix}.bam"),
        truth = ancient(config["REFERENCE"]),
    output:
        catalogue = "{dir}/{prefix}_ray_catalogue.txt"
    log:
        "{dir}/{prefix}_ray_catalogue.log"
    params:
        prefix = "{dir}/{prefix}_ray",
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]) 
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        ray call {input.bam} {input.truth} --threads {threads} --output_prefix {params.prefix} --catalogue &> {log}
        """

rule hp_acc_vs_length:
    input:
        catalogue = ancient("{dir}/{prefix}_ray_catalogue.txt"),
    output:
        hp_acc_sum = "{dir}/{prefix}_ray_summary.txt"
    run:
        from collections import defaultdict
        import pandas as pd

        def get_acc(df):
            correct = df['ref_hp_len'] == df['q_hp_len']
            return 100 * float(len(df[correct])) / len(df)

        df = pd.read_table(input.catalogue)
        accs = defaultdict(dict)
        for hp_len, df_l in df.groupby('ref_hp_len'):
               accs['acc_all_bases'][hp_len] = get_acc(df_l)
               accs['n_all_bases'][hp_len] = len(df_l)
               for base, df_b in df_l.groupby('base'):
                   accs['acc_{}'.format(base)][hp_len] = get_acc(df_b)
                   accs['n_{}'.format(base)][hp_len] = len(df_b)
        summ = pd.DataFrame(accs).reset_index().rename(columns={'index': 'hp_len'})
        pd.DataFrame(summ).to_csv(output.hp_acc_sum, sep=',', index=False)
            

rule get_depth:
    input:
        venv = ancient(IN_POMOXIS),
        bam = ancient("{dir}/calls2ref.bam"),
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
        venv = ancient(IN_POMOXIS),
        bam = ancient("{dir}/calls2ref.bam"),
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
        venv = ancient(IN_POMOXIS),
        basecalls = ancient("{dir}/basecalls.fasta"),
        ref = ancient(config["REFERENCE"])
    output:
        consensus = "{dir}/ref_guided_racon{suffix,[^/]*}/consensus.fasta",
        basecalls = "{dir}/ref_guided_racon{suffix,[^/]*}/basecalls.fasta"
    log:
        "{dir}/ref_guided_racon{suffix}.log"
    threads: config["THREADS_PER_JOB"]
    params:
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]), 
        output_dir = lambda w: "{dir}/ref_guided_racon{suffix}".format(**dict(w))
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        # snakemake will create the output dir, mini_assemble will fail if it exists..
        rm -r {params[output_dir]} && 
        mini_assemble -i {input.basecalls} -r {input.ref} -o {params[output_dir]} -t {threads} -p assm {config[MINI_ASSEMBLE_OPTS]} &&
        # rename output
        mv {params[output_dir]}/assm_final.fa {output.consensus}
        # keep a link of basecalls with the consensus
        ln -s $PWD/{input.basecalls} $PWD/{params[output_dir]}/basecalls.fasta &&
        # sync timestamps, without following basecalls link (otherwise consensus will be older than basecalls)
        touch --no-dereference $PWD/{params[output_dir]}/*
        """

rule miniasm_racon:
    input:
        venv = ancient(IN_POMOXIS),
        basecalls = ancient("{dir}/basecalls.fasta"),
    output:
        consensus = "{dir}/miniasm_racon{suffix,[^/]*}/consensus.fasta",
        basecalls = "{dir}/miniasm_racon{suffix,[^/]*}/basecalls.fasta"
    log:
        "{dir}/miniasm_racon{suffix}.log"
    params:
        output_dir = lambda w: "{dir}/miniasm_racon{suffix}".format(**dict(w)),
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]), 
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        # snakemake will create the output dir, mini_assemble will fail if it exists..
        rm -r {params[output_dir]} && 
        mini_assemble -i {input.basecalls} -o {params[output_dir]} -t {threads} -p assm {config[MINI_ASSEMBLE_OPTS]} &&
        # rename output
        mv {params[output_dir]}/assm_final.fa {output.consensus}
        # keep a link of basecalls with the consensus
        ln -s $PWD/{input.basecalls} $PWD/{params[output_dir]}/basecalls.fasta &&
        # sync timestamps, without following basecalls link (otherwise consensus will be older than basecalls)
        touch --no-dereference $PWD/{params[output_dir]}/*
        """

rule canu:
    input:
        canu = ancient(CANU_EXEC),
        basecalls = ancient("{dir}/basecalls.fasta"),
    output:
        consensus = "{dir}/canu{suffix,[^/]*}/consensus.fasta",
        basecalls = "{dir}/canu{suffix,[^/]*}/basecalls.fasta"
    log:
        "{dir}/canu{suffix}.log"
    params:
        output_dir = lambda w: "{dir}/canu{suffix}".format(**dict(w)),
        genome_sz = config["CANU_GENOME_SIZE"],
        exec_opts = config["CANU_EXEC_OPTS"],
        prefix = "canu",
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]), 
    shell:
        """
        # snakemake will create the output dir, canu will fail if it exists..
        #rm -r {params[output_dir]}
        {input.canu} -d {params.output_dir} -p {params.prefix} genomeSize={params.genome_sz} -nanopore-raw {input.basecalls} {params.exec_opts} &&
        mv {params.output_dir}/{params.prefix}.contigs.fasta {output.consensus} &&
        ln -s $PWD/{input.basecalls} $PWD/{params[output_dir]}/basecalls.fasta &&
        # sync timestamps, without following basecalls link (otherwise consensus will be older than basecalls)
        touch --no-dereference $PWD/{output.consensus} $PWD/{output.basecalls}
        """

rule medaka_consensus:
    input:
        venv = ancient(IN_MEDAKA),
        draft = ancient("{dir}/consensus.fasta"),
        basecalls = ancient("{dir}/basecalls.fasta"),
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

rule nanopolish_basecalls:
    # nanopolish index can't seem to cope with fasta headers
    input:
        ancient("basecall/{subdir}/basecalls.fasta"),
    output:
        "basecall/{subdir}/nanopolish/basecalls.fasta",
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        """
        cut -d' ' -f1 < {input} > {output}
        """

rule nanopolish_index:
    input:
        nanopolish = ancient(NANOPOLISH_EXEC),
        fast5 = ancient(config["READS"]),
        summary = ancient("basecall/{basecaller}/sequencing_summary.txt"),
        basecalls = ancient("basecall/{basecaller}/{subdir}/nanopolish/basecalls.fasta"),
    output:
        "basecall/{basecaller,[^/]*}/{subdir}/nanopolish/basecalls.fasta.index",
        "basecall/{basecaller,[^/]*}/{subdir}/nanopolish/basecalls.fasta.index.gzi",
        "basecall/{basecaller,[^/]*}/{subdir}/nanopolish/basecalls.fasta.index.fai",
        "basecall/{basecaller,[^/]*}/{subdir}/nanopolish/basecalls.fasta.index.readdb",
    log:
        "basecall/{basecaller,[^/]*}/{subdir}/nanopolish/nanopolish_index.log"
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
    	# create indices then synchronise time stamps
	    """
        {input.nanopolish} index -d {input.fast5} -s {input.summary} {input.basecalls} &> {log} && sleep 5 &&
    	touch --no-dereference {output[0]}*
        """

rule fast5_list:
    input:
        fast5 = ancient(config["READS"]),
    output:
        filelist = os.path.join(config["READS"], "reads.txt")
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        "find -L `readlink -f {config[READS]}` -name '*.fast5' > {output.filelist}"

rule miyagi_index:
    input:
        filelist = ancient(os.path.join(config["READS"], "reads.txt")),
        summary = ancient("basecall/{basecaller}/sequencing_summary.txt"),
        basecalls = ancient("basecall/{basecaller}/align/{subdir}/basecalls.fasta"),
    output:
        index = "basecall/{basecaller}/align/{subdir}/miyagi/miyagi.index",
        basecalls = "basecall/{basecaller}/align/{subdir}/miyagi/basecalls.fasta",
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
        nanopolish = ancient(NANOPOLISH_EXEC),
        draft = ancient("basecall/{basecaller}/{subdir}/consensus.fasta"),
        basecalls = ancient("basecall/{basecaller}/{subdir}/nanopolish/basecalls.fasta"),
        bam = ancient("basecall/{basecaller}/{subdir}/nanopolish/calls2draft.bam"),
        index = ancient("basecall/{basecaller}/{subdir}/nanopolish/basecalls.fasta.index"),
        index_gzi = ancient("basecall/{basecaller}/{subdir}/nanopolish/basecalls.fasta.index.gzi"),
        index_fai = ancient("basecall/{basecaller}/{subdir}/nanopolish/basecalls.fasta.index.fai"),
        index_readdb = ancient("basecall/{basecaller}/{subdir}/nanopolish/basecalls.fasta.index.readdb"),
    output:
        vcf = "basecall/{basecaller,[^/]*}/{subdir}/nanopolish/regions/{region}.vcf",
    log:
        "basecall/{basecaller,[^/]*}/{subdir}/nanopolish/regions/{region}.vcf.log"
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
	    """
        {input.nanopolish} variants --consensus -o {output.vcf} -w {wildcards.region} -r {input.basecalls} -b {input.bam} -g {input.draft} -t 1 {config[NP_OPTS]} &> {log}
        """

rule nanopolish_regions:
    input:
        # use pomoxis python as this has all requirements of the script
        venv = ancient(IN_POMOXIS), 
        make_range = ancient(NANOPOLISH_MAKE_RANGE),
        draft = ancient("{dir}/consensus.fasta"),
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
        nanopolish = ancient(NANOPOLISH_EXEC),
        draft = ancient("{dir}/consensus.fasta"),
        # Number of regions is unknown ahead of time, so use dynamic keyword to delay evaluation
        regions = ancient(dynamic("{dir}/nanopolish/regions/{region}.region")),
        vcfs = ancient(dynamic("{dir}/nanopolish/regions/{region}.vcf")),
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
        venv = ancient(IN_MIYAGI),
        draft = ancient("{dir}/consensus.fasta"),
    output:
        vcf = "{dir}/miyagi/regions/{region}.hp.vcf"
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        miyagi_create_homop_vcf {input.draft} {wildcards.region} {output.vcf} {config[MIYAGI_HOMOP_VCF_OPTS]} 
        """

rule miyagi_score_reads:
    input:
        venv = ancient(IN_MIYAGI),
        draft = ancient("{dir}/consensus.fasta"),
        #region = ancient("{dir}/miyagi/regions/{region}.region"),
        bam = ancient("{dir}/miyagi/calls2draft.bam"),
        index = ancient("{dir}/miyagi/miyagi.index"),
        vcf = ancient("{dir}/miyagi/regions/{region}.hp.vcf"),
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
        ancient(wildcards_2_score_targets),
        venv = ancient(IN_MIYAGI),
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
        venv = ancient(IN_MIYAGI),
        draft = ancient("{dir}/consensus.fasta"),
        hdfs = ancient(wildcards_2_sorted_score_targets),
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
        venv = ancient(IN_MIYAGI),
        draft = ancient("{dir}/consensus.fasta"),
        vcf = ancient("{dir}/miyagi/corrections.vcf"),
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
