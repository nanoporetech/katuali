from functools import partial
import itertools
import os
import sys

import gpustat
import pysam

KATUALI_VERSION = "0.1"

KATUALI_HOME = os.path.split(sys.exec_prefix)[0]

PICK_GPU = os.path.join(KATUALI_HOME, 'scripts', 'pick_gpu.py') 

GUPPY_EXEC = os.path.expanduser(config["GUPPY"])
SCRAPPIE_EXEC = os.path.join(os.path.expanduser(config["SCRAPPIE"]), "build", "scrappie")
FLAPPIE_EXEC = os.path.join(os.path.expanduser(config["FLAPPIE"]), "flappie")
SCRAPPIE_JSON_TO_TSV = os.path.join(os.path.expanduser(config["SCRAPPIE"]), "misc", "json_to_tsv.py")
NANOPOLISH_EXEC = os.path.join(os.path.expanduser(config["NANOPOLISH"]), "nanopolish")
NANOPOLISH_MAKE_RANGE = os.path.join(os.path.expanduser(config["NANOPOLISH"]), "scripts", "nanopolish_makerange.py")
IN_POMOXIS = os.path.expanduser(config["IN_POMOXIS"])
IN_MEDAKA = os.path.expanduser(config["IN_MEDAKA"])
IN_RAY = os.path.expanduser(config["IN_RAY"])
CANU_EXEC = os.path.expanduser(config["CANU_EXEC"])

# Always set the environment
LD_LIBRARY_PATH = config.get("LD_LIBRARY_PATH")
shell.prefix("export LD_LIBRARY_PATH=\"%s\"; " % (LD_LIBRARY_PATH))

config["THREADS_PER_JOB"] = int(config["THREADS_PER_JOB"])

# NOTE on virtual environments
# Snakemake uses bash strict mode, virtualenv violates bash strict mode.
# https://snakemake.readthedocs.io/en/stable/project_info/faq.html#my-shell-command-fails-with-with-errors-about-an-unbound-variable-what-s-wrong
# so any activation commands must be wrapped as such:
#      set +u; {config[IN_MEDAKA]}; set -u; 

rule version:
    run:
        print('katuali version {}'.format(KATUALI_VERSION))

rule help:
    run:
        print(open(os.path.join(KATUALI_HOME, "README.md")).read())

rule environment:
    # Just output the contents of PATH and LD_LIBRARY_PATH to log
    log:
        "katuali.sge.log",
    params:
        sge = "m_mem_free=1G,gpu=0", 
    shell:
        'echo $PATH; echo $LD_LIBRARY_PATH'


# shortcut targets to process all runs in config["DATA"] with a standard pipeline
rule fast_assm_polish:
    input:
        consensus = ancient(["{runid}/basecall/guppy/miniasm_racon/medaka/consensus.fasta".format(runid=r) for r in config["DATA"].keys()])

rule standard_assm_polish:
    input: 
        consensus = ancient(["{runid}/basecall/guppy_flipflop/canu_gsz_{gs}/racon/medaka_flipflop/consensus.fasta".format(runid=r, gs=d["GENOME_SIZE"]) for r,d in config["DATA"].items() if "GENOME_SIZE" in d] )

rule standard_assm_nanopolish:
    input: 
        consensus = ancient(["{runid}/basecall/{bc}/canu_gsz_{gs}/racon/nanopolish_hp/consensus.fasta".format(runid=r, gs=d["GENOME_SIZE"], bc=config["BASECALLER"]) for r,d in config["DATA"].items() if "GENOME_SIZE" in d])


def get_reference(wildcards, config):
    return config["DATA"][wildcards["runid"]]["REFERENCE"]


def get_contig_opt(wildcards, config):
    if "REGIONS" in config and "REGIONS" != "":
        contig_opt = "-r {}".format(config["REGIONS"])
    elif wildcards.contig == "all_contigs":
        contig_opt = ""
    else:
        contig_opt = "-r {}".format(wildcards.contig)
    logger.run_info("Setting region option to {}".format(contig_opt))
    return contig_opt


def get_opts(wildcards, config, config_key):
    default_key = ''

    if config_key not in config:
        raise KeyError('{} not in config'.format(config_key))

    suffix = wildcards["suffix"]

    if not isinstance(config[config_key], dict):
        opts = config[config_key]
        logger.run_info("{} parameters were not nested, using {}".format(config_key, opts))

    elif suffix in config[config_key] and suffix != default_key:
        opts = config[config_key][suffix]
        logger.run_info("Using {} parameters specified by suffix {}".format(config_key, suffix))
    else:
        opts = config[config_key][default_key]
        logger.run_info("Using default {} parameters".format(config_key))

    return opts
    
    
rule basecall_scrappie:
    input:
        scrappie = ancient(SCRAPPIE_EXEC),
        fast5 = ancient("{runid}/reads"),
    output:
        fasta = "{runid}/basecall/scrappie{suffix,[^/]*}/basecalls.fasta",
    log:
        "{runid}/basecall/scrappie{suffix,[^/]*}/scrappie.log",
    params:
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]), 
        opts = partial(get_opts, config=config, config_key="SCRAPPIE_OPTS")
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        echo "{input.scrappie} {config[SCRAPPIE_OPTS]} "> {log}
        find -L {input.fast5} -name '*.fast5' | xargs {input.scrappie} {params[opts]} --threads {threads} >{output.fasta} 2>> {log}
        """


rule basecall_flappie:
    input:
        flappie = ancient(FLAPPIE_EXEC),
        fast5 = ancient("{runid}/reads"),
        venv = ancient(IN_POMOXIS),
    output:
        fastq = "{runid}/basecall/flappie{suffix,[^/]*}/basecalls.fastq",
        fasta = "{runid}/basecall/flappie{suffix,[^/]*}/basecalls.fasta",
    log:
        "{runid}/basecall/flappie{suffix,[^/]*}/flappie.log",
    params:
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]), 
        opts = partial(get_opts, config=config, config_key="FLAPPIE_OPTS"),
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        export OPENBLAS_NUM_THREADS=1
        echo "{input.flappie} {params[opts]} "> {log}
        #find -L {input.fast5} -name '*.fast5' | parallel -P {threads} -X {input.flappie} {params[opts]} >{output.fastq} 2>> {log}
        find -L {input.fast5} -name '*.fast5' | parallel -j {threads} -n 10 -X {input.flappie} {params[opts]} >{output.fastq} 2>> {log}
        set +u; {config[SOURCE]} {input.venv}; set -u;
        sleep 5
        seqkit fq2fa {output.fastq} > {output.fasta}
        """


rule basecall_guppy:
    input:
        guppy = ancient(GUPPY_EXEC),
        PICK_GPU = ancient(PICK_GPU),
        fast5 = ancient("{runid}/reads"),
        venv = ancient(IN_POMOXIS),
    output:
        fasta = "{runid}/basecall/guppy{suffix,[^/]*}/basecalls.fasta",
        summary = "{runid}/basecall/guppy{suffix,[^/]*}/sequencing_summary.txt",
    log:
        "{runid}/basecall/guppy{suffix,[^/]*}.log",
    params:
        sge = "m_mem_free=1G,gpu=1 -pe mt {}".format(config["GUPPY_SLOTS"]),
        output_dir = lambda w: "{runid}/basecall/guppy{suffix}".format(**dict(w)),
        opts = partial(get_opts, config=config, config_key="GUPPY_OPTS"),
    shell:
        """
        # snakemake will create the output dir, guppy will fail if it exists..
        rm -r {params[output_dir]}
 
        echo "GPU status before" >> {log}
        gpustat >> {log}
    
        sleep $(((RANDOM % 30)  + 1 ))        

        GPU=$({input.PICK_GPU} 2>> {log})

        echo "Runnning on host $HOSTNAME GPU $GPU" >> {log}

        {input.guppy} -s {params.output_dir} -r -i {input.fast5} -x cuda:$GPU {params.opts} &>> {log}

        echo "gpustat after" >> {log}
        gpustat >> {log}

        # convert fastq to fasta
        sleep 5
        echo "Combining the following fastq files into {output.fasta}" >> {log}
        ls {params[output_dir]}/*.fastq >> {log}
        set +u; {config[SOURCE]} {input.venv}; set -u;
        seqkit fq2fa {params.output_dir}/*.fastq > {output.fasta}
        """


rule scrappie_summary:
    input:
        json_to_tsv = ancient(SCRAPPIE_JSON_TO_TSV),
        fasta = ancient("{runid}/basecall/scrappie/basecalls.fasta"),
    output:
        summary = "{runid}/basecall/scrappie/sequencing_summary.txt"
    params:
        sge = "m_mem_free=1G,gpu=0"
    shell:
    	"cat {input.fasta} | grep '^>' | cut -d ' ' -f 2- | python {input.json_to_tsv} > {output.summary} &&"
    	# most tools expect read_id, not uuid
	    "sed -i '1 s/uuid/read_id/' {output.summary}"


rule align_to_ref:
    input:
        venv = ancient(IN_POMOXIS),
        basecalls = ancient("{runid}/basecall/{bc_dir}/basecalls.fasta"),
        ref = ancient(partial(get_reference, config=config)),
    output:
        bam = "{runid}/basecall/{bc_dir}/align/calls2ref.bam"
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
        consensus = ancient("{runid,[^/]+}/{dir}/consensus.fasta"),
        truth = ancient(partial(get_reference, config=config)),
    output:
        summ = "{runid,[^/]+}/{dir}/consensus_to_truth_summ.txt",
        bam = "{runid,[^/]+}/{dir}/consensus_to_truth.bam",
        stats = "{runid,[^/]+}/{dir}/consensus_to_truth_stats.txt",
    params:
        prefix = "{runid}/{dir}/consensus_to_truth",
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]), 
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        assess_assembly -i {input.consensus} -r {input.truth} -p {params.prefix} -t {threads} {config[ASSESS_ASSM_OPTS]} 
        """


rule ray_catalogue:
    input:
        venv = ancient(IN_RAY),
        bam = ancient("{runid}/{dir}/{prefix}.bam"),
        truth = ancient(partial(get_reference, config=config)),
    output:
        catalogue = "{runid,[^/]+}/{dir}/{prefix}_ray_catalogue.txt"
    log:
        "{runid}/{dir}/{prefix}_ray_catalogue.log"
    params:
        prefix = "{runid}/{dir}/{prefix}_ray",
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]) 
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        ray call {input.bam} {input.truth} --threads {threads} --output_prefix {params.prefix} --catalogue &> {log}
        """

rule hp_acc_vs_length:
    input:
        catalogue = ancient("{runid}/{dir}/{prefix}_ray_catalogue.txt"),
    output:
        hp_acc_sum = "{runid}/{dir}/{prefix}_ray_summary.txt"
    params:
        sge = "m_mem_free=1G,gpu=0"
    run:
        from collections import defaultdict
        import pandas as pd

        def get_acc(df):
            correct = df['ref_hp_len'] == df['q_hp_len']
            return 100 * float(len(df[correct])) / len(df)

        def get_summ(df):
            accs = defaultdict(dict)
            for hp_len, df_l in df.groupby('ref_hp_len'):
                   accs['acc_all_bases'][hp_len] = get_acc(df_l)
                   accs['n_all_bases'][hp_len] = len(df_l)
                   for base, df_b in df_l.groupby('ref_base'):
                       accs['acc_{}'.format(base)][hp_len] = get_acc(df_b)
                       accs['n_{}'.format(base)][hp_len] = len(df_b)
            summ = pd.DataFrame(accs).reset_index().rename(columns={'index': 'hp_len'})
            return summ

        df = pd.read_table(input.catalogue)
        # create a summary over all refs, and one per reference
        get_summ(df).to_csv(output.hp_acc_sum, sep=',', index=False)
        for ref, d in df.groupby('record_name'):
            out = output.hp_acc_sum.replace('_ray_summary.txt', '_{}_ray_summary.txt'.format(ref))
            get_summ(d).to_csv(out, sep=',', index=False)
            

rule get_depth:
    input:
        venv = ancient(IN_POMOXIS),
        bam = ancient("{dir}/calls2ref.bam"),
    output:
        depth = directory("{dir}/depth")
    params:
        contig_opt = partial(get_contig_opt, config=config),
        sge = "m_mem_free=1G,gpu=0", 
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
    	mkdir -p {output.depth} && coverage_from_bam {input.bam} {params.contig_opt} -s 1000 -p {output.depth}/ &>{output.depth}/depth.log
        """


rule get_basecall_stats:
    input:
        venv = ancient(IN_POMOXIS),
        bam = ancient("{dir}/calls2ref.bam"),
    output:
        stats = "{dir}/calls2ref_stats.txt"
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
    	stats_from_bam {input.bam} -o {output.stats}
        """

rule subsample_bam:
    input:
        venv = ancient(IN_POMOXIS),
        bam = ancient("{dir}/calls2ref.bam"),
    output:
        fasta = "{dir}/{contig}/{depth,[0-9]+}X{suffix,[^/]*}/basecalls.fasta",
    log:
        "{dir}/{contig}/{depth}X{suffix}/subsample.log"
    params:
        contig_opt = partial(get_contig_opt, config=config),
        prefix = lambda w: "{dir}/{contig}/{depth}X{suffix}/sub_sample".format(**dict(w)),
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]),
        opts = partial(get_opts, config=config, config_key="SUBSAMPLE_BAM_OPTS"),
    threads: config["THREADS_PER_JOB"]
    shell:
	    """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        subsample_bam {input.bam} {wildcards.depth} {params[contig_opt]} -o {params[prefix]} {params[opts]} -t {threads} &>{log};
        sleep 5;
        for i in {params[prefix]}*.bam; do samtools fasta $i ; done > {output.fasta}
        """


rule racon:
    input:
        venv = ancient(IN_POMOXIS),
        draft = ancient("{dir}/consensus.fasta"),
        basecalls = ancient("{dir}/basecalls.fasta"),
    output:
        consensus = "{dir}/racon{suffix,[^/]*}/consensus.fasta",
        basecalls = "{dir}/racon{suffix,[^/]*}/basecalls.fasta"
    log:
        "{dir}/racon{suffix}.log"
    threads: config["THREADS_PER_JOB"]
    params:
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]), 
        output_dir = lambda w: "{dir}/racon{suffix}".format(**dict(w)),
        opts = partial(get_opts, config=config, config_key="MINI_ASSEMBLE_OPTS"),

    threads: config["THREADS_PER_JOB"]
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        # snakemake will create the output dir, mini_assemble will fail if it exists..
        rm -r {params[output_dir]} && 
        mini_assemble -i {input.basecalls} -r {input.draft} -o {params[output_dir]} -t {threads} -p assm {params[opts]} &> {log}
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
        opts = partial(get_opts, config=config, config_key="MINI_ASSEMBLE_OPTS"),
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        # snakemake will create the output dir, mini_assemble will fail if it exists..
        rm -r {params[output_dir]} && 
        mini_assemble -i {input.basecalls} -o {params[output_dir]} -t {threads} -p assm {params[opts]} &> {log}
        # rename output
        mv {params[output_dir]}/assm_final.fa {output.consensus}
        # keep a link of basecalls with the consensus
        ln -s $PWD/{input.basecalls} $PWD/{params[output_dir]}/basecalls.fasta &&
        # sync timestamps, without following basecalls link (otherwise consensus will be older than basecalls)
        touch --no-dereference $PWD/{params[output_dir]}/*
        """

rule canu:
    # TODO can we make use of github.com/thiesgehrmann/FungalAssemblerPipeline ?
    input:
        canu = ancient(CANU_EXEC),
        basecalls = ancient("{dir}/basecalls.fasta"),
    output:
        consensus = "{dir}/canu{suffix,[^/]*}_gsz_{genome_size,[^/]*}/consensus.fasta",
        basecalls = "{dir}/canu{suffix,[^/]*}_gsz_{genome_size,[^/]*}/basecalls.fasta"
    log:
        "{dir}/canu{suffix}gsz_{genome_size}.log"
    params:
        gsz = lambda w: "genomeSize={genome_size}".format(**dict(w)),
        output_dir = lambda w: "{dir}/canu{suffix}gsz_{genome_size}".format(**dict(w)),
        exec_opts = config["CANU_EXEC_OPTS"],
        opts = partial(get_opts, config=config, config_key="CANU_OPTS"),
        prefix = "canu",
        mt="maxThreads={}".format(config["THREADS_PER_JOB"]),
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]),
    threads: config["THREADS_PER_JOB"]
    shell:
        """
        # snakemake will create the output dir, canu will fail if it exists..
        #rm -r {params[output_dir]}
        {input.canu} -d {params.output_dir} -p {params.prefix} {params.gsz} -nanopore-raw {input.basecalls} {params.exec_opts} {params[opts]} {params.mt} &> {log}
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
        consensus = "{dir}/medaka{suffix,[^/]*}/consensus.fasta",
        basecalls = "{dir}/medaka{suffix,[^/]*}/basecalls.fasta"
    log:
        "{dir}/medaka{suffix}.log"
    threads: 8
    params:
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]),
        opts = partial(get_opts, config=config, config_key="MEDAKA_OPTS"),
        output_dir = lambda w: "{dir}/medaka{suffix}".format(**dict(w))
    shell:
        """
        set +u; {config[SOURCE]} {input.venv}; set -u;
        # snakemake will create the output dir if it does not exist, remove it it exists. 
        rm -r {params[output_dir]} 

        medaka_consensus -i {input.basecalls} -d {input.draft} -o {params[output_dir]} -t {threads} {params[opts]} &> {log}

        # keep a link of basecalls with the consensus
        ln -s $PWD/{input.basecalls} $PWD/{output.basecalls}
        """

rule nanopolish_basecalls:
    # nanopolish index can't seem to cope with fasta headers
    input:
        ancient("{runid}/basecall/{subdir}/basecalls.fasta"),
    output:
        "{runid}/basecall/{subdir}/nanopolish/basecalls.fasta",
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        """
        cut -d' ' -f1 < {input} > {output}
        """

rule nanopolish_index:
    input:
        nanopolish = ancient(NANOPOLISH_EXEC),
        fast5 = ancient("{runid}/reads"),
        summary = ancient("{runid}/basecall/{basecaller}/sequencing_summary.txt"),
        basecalls = ancient("{runid}/basecall/{basecaller}/{subdir}/nanopolish/basecalls.fasta"),
    output:
        "{runid}/basecall/{basecaller,[^/]+}/{subdir}/nanopolish/basecalls.fasta.index",
        "{runid}/basecall/{basecaller,[^/]+}/{subdir}/nanopolish/basecalls.fasta.index.gzi",
        "{runid}/basecall/{basecaller,[^/]+}/{subdir}/nanopolish/basecalls.fasta.index.fai",
        "{runid}/basecall/{basecaller,[^/]+}/{subdir}/nanopolish/basecalls.fasta.index.readdb",
    log:
        "{runid}/basecall/{basecaller}/{subdir}/nanopolish/nanopolish_index.log"
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
	    """
        {input.nanopolish} index -d {input.fast5} -s {input.summary} {input.basecalls} &> {log}
        """

rule fast5_list:
    input:
        fast5 = ancient("{runid}/reads")
    output:
        filelist = "{runid}/reads/reads.txt"
    params:
        sge = "m_mem_free=1G,gpu=0" 
    shell:
        """
        find -L `readlink -f {input.fast5}` -name '*.fast5' > {output.filelist}
        """

rule nanopolish_vcf:
    input:
        nanopolish = ancient(NANOPOLISH_EXEC),
        draft = ancient("{dir}/consensus.fasta"),
        basecalls = ancient("{dir}/nanopolish/basecalls.fasta"),
        bam = ancient("{dir}/nanopolish/calls2draft.bam"),
        index = ancient("{dir}/nanopolish/basecalls.fasta.index"),
        index_gzi = ancient("{dir}/nanopolish/basecalls.fasta.index.gzi"),
        index_fai = ancient("{dir}/nanopolish/basecalls.fasta.index.fai"),
        index_readdb = ancient("{dir}/nanopolish/basecalls.fasta.index.readdb"),
    output:
        vcf = "{dir}/nanopolish/regions/{region}.vcf",
    log:
        "{dir}/nanopolish/regions/{region}.vcf.log"
    params:
        sge = "m_mem_free=1G,gpu=0", 
        # wildcards in dynamic files cannot be constrained => we can't safely extract a
        # suffix from dynamic nanopolish targets to use to use nested config
        opts = config["NANOPOLISH_OPTS"],
    shell:
	    """
        {input.nanopolish} variants --consensus -o {output.vcf} -w {wildcards.region} -r {input.basecalls} -b {input.bam} -g {input.draft} -t 1 {params.opts} &> {log}
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
        sge = "m_mem_free=1G,gpu=0", 
    shell:
        "{input.nanopolish} vcf2fasta -g {input.draft} {input.vcfs} > {output.consensus} 2> {log}"


rule medaka_train_features:
    input:
        in_medaka = ancient(IN_MEDAKA),
        in_pomoxis = ancient(IN_POMOXIS),
        draft = ancient("{runid}/{dir}/consensus.fasta"),
        basecalls = ancient("{runid}/{dir}/basecalls.fasta"),
        truth = ancient(partial(get_reference, config=config)),
    output:
        features = "{runid,[^/]+}/{dir}/medaka_train{suffix,[^/]*}/medaka_train.hdf",
        rc_features = "{runid,[^/]+}/{dir}/medaka_train{suffix,[^/]*}/medaka_train_rc.hdf",
    log:
        "{runid}/{dir}/medaka_train{suffix}.log"
    threads: config["THREADS_PER_JOB"]
    params:
        output_dir = lambda w, output: os.path.dirname(output.features),
        bam = lambda w, output: os.path.join(os.path.dirname(output.features), "calls2ref"),
        truth_bam = lambda w, output: os.path.join(os.path.dirname(output.features), "truth2ref"),
        rc_bam = lambda w, output: os.path.join(os.path.dirname(output.features), "callsrc2ref"),
        rc_truth_bam = lambda w, output: os.path.join(os.path.dirname(output.features), "truth2refrc"),
        rc_draft = lambda w, output: os.path.join(os.path.dirname(output.features), "draftrc.fasta"),
        rc_truth = lambda w, output: os.path.join(os.path.dirname(output.features), "truthrc.fasta"),
        sge = "m_mem_free=1G,gpu=0 -pe mt {}".format(config["THREADS_PER_JOB"]),
        opts = partial(get_opts, config=config, config_key="MEDAKA_TRAIN_FEAT_OPTS"),
        truth_chunk = config["MEDAKA_TRAIN_TRUTH_CHUNK"],
    shell:
        """
        set +u; {config[SOURCE]} {input.in_pomoxis}; set -u;

        # keep a link of basecalls with the consensus
        
        ln -sf $PWD/{input.basecalls} $PWD/{params.output_dir}/basecalls.fasta
        sleep 1

        echo "aligning basecalls to draft" >{log}
        mini_align -i {input.basecalls} -r {input.draft} -p {params.bam} -t {threads} -P -m &>> {log}
        sleep 5

        echo "aligning truth to draft" >>{log}
        mini_align -i {input.truth} -r {input.draft} -p {params.truth_bam} -t {threads} -P -m -c {params.truth_chunk} &>> {log}
        sleep 5

        echo "reverse complement the draft and align reads" >> {log}
        seqkit seq --complement --reverse {input.draft} -o {params.rc_draft} &>> {log}
        sleep 5

        echo "aligning basecalls to rc draft" >>{log}
        mini_align -i {input.basecalls} -r {params.rc_draft} -p {params.rc_bam} -t {threads} -P -m &>> {log}
        sleep 5

        echo "aligning truth to rc draft" >> {log}
        mini_align -i {input.truth} -r {params.rc_draft} -p {params.rc_truth_bam} -t {threads} -P -m -c {params.truth_chunk} &>> {log}
        sleep 5

        echo "creating features" >> {log}
        set +u; {config[SOURCE]} {input.in_medaka} set -u;
        medaka features {params.bam}.bam {output.features} --truth {params.truth_bam}.bam {params[opts]} --threads {threads} &>> {log}
        echo "creating rc features"
        medaka features {params.rc_bam}.bam {output.rc_features} --truth {params.rc_truth_bam}.bam {params[opts]} --threads {threads} &>> {log}

        """


def get_region_len(region, ref_fasta):
    """Get region length from samtools region string, using start and end if present, else obtaining contig length from reference fasta. 
    """
    if ':' in region:
        if not '-' in region:
            raise ValueError('Regions must be specified just as contig or contig:start-end')
        start, end = map(int, region.split(':')[1].split('-'))
        region_len = end -start
        if region_len < 1:
            raise ValueError('Region length < 1, check your region specification')
    else:
        # we have a full contig
        with pysam.FastaFile(ref_fasta) as fa:
            lengths = dict(zip(fa.references, fa.lengths))
        if region not in fa.references:
            raise KeyError('Region {} is not a contig in the reference'.format(region))
        region_len = lengths[region]
        
    return region_len


def format_genome_len(l, fmt='{:.1f}'):
    """Format a genome length supplied as int as string with either k, m or g suffix"""
    ranges = [('k', 10**3), ('m', 10**6), ('g', 10**9),]
    for suffix, scale in ranges:
        scaled = l / scale
        if scaled < 100 or scale == ranges[-1][1]:
            return fmt.format(scaled) + suffix


def get_medaka_features_targets(config):
    template = "{runid}/basecall/{basecaller}/align/{region}/{depth}X_prop/canu_gsz_{genome_size}/racon/medaka_train/{feature_file}"
    targets = []
    try:
        for runid, d in config["DATA"].items():
            for region in d["MEDAKA_TRAIN_REGIONS"]:
                region_size = get_region_len(region, d["REFERENCE"])  # an int
                region_size = format_genome_len(region_size)  # as str
                targets.extend([template.format(runid=runid,
                                                basecaller=config["BASECALLER"],
                                                region=region,
                                                depth=depth,
                                                genome_size=region_size,
                                                feature_file=f)
                                for depth in config["DEPTHS"] for f in ["medaka_train.hdf", "medaka_train_rc.hdf"]
                ])
    except:
        pass

    if config["USE_ONLY_EXISTING_MEDAKA_FEAT"]:
        logger.run_info("WARNING: USE_ONLY_EXISTING_MEDAKA_FEAT set to true, only using existing medaka training features.")
        n_all_feat = len(targets)
        targets = [f for f in targets if os.path.isfile(f)]
        if len(targets) == 0:
            raise ValueError("USE_ONLY_EXISTING_MEDAKA_FEAT set to true in config, but no features were found")
        logger.run_info("Found {}/{} of the medaka training feature files.".format(len(targets), n_all_feat))

    return targets


rule all_medaka_train_features:
    input:
        venv = ancient(IN_MEDAKA),
        features = ancient(get_medaka_features_targets(config)),
    log:
        "medaka_train_features.log"


rule train_medaka:
    input:
        PICK_GPU = ancient(PICK_GPU),
        venv = ancient(IN_MEDAKA),
        features = ancient(get_medaka_features_targets(config)),
    output:
        train_dir = directory("medaka_train_{suffix,[^/]*}")
    log:
        "medaka_train_{suffix}.log",
    params:
        sge = "m_mem_free=1G,gpu=1",
        opts = partial(get_opts, config=config, config_key="MEDAKA_TRAIN_OPTS"),
    shell:
        """
        echo "GPU status before" >> {log}
        gpustat >> {log}

        GPU=$({input.PICK_GPU} 2>> {log})

        echo "Runnning on host $HOSTNAME GPU $GPU" >> {log}

        if [ "{config[SCRATCH]}" != "" ];then
            if [ "{config[TMPSCRATCH]}" != "" ];then
                tmpscr={config[TMPSCRATCH]}
                echo "Using provided TMPSCRATCH {config[TMPSCRATCH]}"
                mkdir -p $tmpscr
            else
                userscr={config[SCRATCH]}/$USER
                mkdir -p $userscr
                tmpscr=$(mktemp -d -p $userscr)
            fi
            t=$(date +"%T")
            echo "$t: Copying feature files to specified scratch directory: $tmpscr."
            for f in {input.features}; do
                d=$tmpscr/$f;
                t0=$(date +"%T");
                if [[ ! -e $d ]]; then 
                    echo "$t0 Copying $d";
                    mkdir -p $(dirname $d) && sleep 1 && cp $f $d; 
                    t1=$(date +"%T");
                    echo "$t1 Copied $d";
                else
                    echo "Found existing file: $d"
                fi
            done
            features=$(for f in {input.features}; do echo $tmpscr/$f; done)
        else
            features={input.features}
        fi

        set +u; {config[SOURCE]} {input.venv} set -u;
        medaka train $features --train_name {output.train_dir} {params.opts} --device $GPU &>> {log}
        """

rule medaka_train_replicates:
    input:
        ancient(expand("medaka_train_{replicate}", replicate=config["MEDAKA_TRAIN_REPLICATES"]))
    log:
        "medaka_train_replicates.log",


def canu_racon_medaka_twice_targets(config):
    template = "{runid}/basecall/{basecaller}/align/{region}/{depth}X/canu_gsz_{genome_size}/racon/medaka{suff}/medaka{suff}/consensus_to_truth_summ.txt"
    targets = []
    try:
        for runid, d in config["DATA"].items():
            for region in d["MEDAKA_EVAL_REGIONS"]:
                region_size = get_region_len(region, d["REFERENCE"])  # an int
                region_size = format_genome_len(region_size)  # as str
                targets.extend([template.format(runid=runid,
                                                basecaller=config["BASECALLER"],
                                                region=region,
                                                depth=depth,
                                                genome_size=region_size,
                                                suff=suffix)
                                for depth in config["DEPTHS"] for suffix in config["MEDAKA_EVAL_SUFFIXES"]
                ])
    except:
        pass

    return targets


rule all_medaka_eval:
    input:
        targets = ancient(canu_racon_medaka_twice_targets(config)),
    log:
        "medaka_eval.log"

