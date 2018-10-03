.PHONY: install

venv: venv/bin/activate
IN_VENV=. ./venv/bin/activate
SNAKEFILE=${PWD}/Snakefile
JOBS=$(shell nproc)

DATADIR=test/data
TESTDIR=${DATADIR}/basecall

TEST=${IN_VENV} && cd test/data && snakemake -s ${SNAKEFILE} --printshellcmds -j ${JOBS}

OPT='--config SCRAPPIE_OPTS="raw -H mean --model rgrgr_r94 --local 10.0 --uuid --temperature1 0.65 --temperature2 1.7"'


venv/bin/activate:
	test -d venv || virtualenv venv --prompt '(katuali) ' --python=python3
	${IN_VENV} && pip install pip --upgrade
	${IN_VENV} && pip install -r requirements.txt

install: venv

reads:
	cd test/data && tar -xvf reads.tgz

test: install test_basecall_suffix test_align test_subsample test_racon test_racon_suffix test_medaka test_nanopolish test_canu test_miniasm_racon test_nanopolish_from_scratch 

test_basecall: clean reads
	${TEST} basecall/scrappie/basecalls.fasta

test_basecall_suffix: clean reads
	${TEST} basecall/scrappie_alt/basecalls.fasta

test_align: reads 
	${TEST} basecall/scrappie/align/calls2ref.bam

test_subsample: reads
	${TEST} basecall/scrappie/align/all_contigs/25X/basecalls.fasta

test_racon: reads
	${TEST} basecall/scrappie/align/all_contigs/25X/ref_guided_racon/consensus.fasta

test_racon_one_contig: reads
	${TEST} basecall/scrappie/align/ecoli_SCS110_plasmid2/25X/ref_guided_racon/consensus.fasta

test_racon_suffix: reads
	${TEST} basecall/scrappie/align/all_contigs/25X/ref_guided_racon_alt/consensus.fasta

test_medaka: reads
	${TEST} basecall/scrappie/align/all_contigs/25X/ref_guided_racon/medaka/consensus.fasta

test_nanopolish: reads
	${TEST} basecall/scrappie/align/all_contigs/25X/ref_guided_racon/nanopolish/consensus.fasta

test_canu: reads
	${TEST} basecall/scrappie/canu/consensus.fasta

test_miniasm_racon: reads
	${TEST} basecall/scrappie/miniasm_racon/consensus.fasta

test_nanopolish_from_scratch: clean test_nanopolish

clean:
	rm -rf ${TESTDIR} ${DATADIR}/reads ${DATADIR}/ref.fasta.*