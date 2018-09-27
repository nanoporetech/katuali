.PHONY: install

venv: venv/bin/activate
IN_VENV=. ./venv/bin/activate

DATADIR=test/data
TESTDIR=${DATADIR}/analysis

TEST=${IN_VENV} && cd test/data && snakemake -s ~/git/katuali/Snakefile --printshellcmds -j 8

OPT='--config SCRAPPIE_OPTS="raw -H mean --model rgrgr_r94 --local 10.0 --uuid --temperature1 0.65 --temperature2 1.7"'


venv/bin/activate:
	test -d venv || virtualenv venv --prompt '(katuali) ' --python=python3
	${IN_VENV} && pip install pip --upgrade
	${IN_VENV} && pip install -r requirements.txt

install: venv

reads:
	cd test/data && tar -xvf reads.tgz

test: install test_basecall_suffix test_align test_subsample test_racon test_racon_suffix test_medaka test_nanopolish test_nanopolish_from_scratch 

test_basecall: clean reads
	${TEST} analysis/basecall/scrappie/basecalls.fasta

test_basecall_suffix: clean reads
	${TEST} analysis/basecall/scrappie_alt/basecalls.fasta

test_align: reads 
	${TEST} analysis/basecall/scrappie/align/calls2ref.bam

test_subsample: reads
	${TEST} analysis/basecall/scrappie/align/all_contigs/25X/basecalls.fasta

test_racon: reads
	${TEST} analysis/basecall/scrappie/align/all_contigs/25X/ref_guided_racon/consensus.fasta

test_racon_one_contig: reads
	${TEST} analysis/basecall/scrappie/align/ecoli_SCS110_plasmid2/25X/ref_guided_racon/consensus.fasta

test_racon_suffix: reads
	${TEST} analysis/basecall/scrappie/align/all_contigs/25X/ref_guided_racon_alt/consensus.fasta

test_medaka: reads
	${TEST} analysis/basecall/scrappie/align/all_contigs/25X/ref_guided_racon/medaka/consensus.fasta

test_nanopolish: reads
	${TEST} analysis/basecall/scrappie/align/all_contigs/25X/ref_guided_racon/nanopolish/consensus.fasta

test_nanopolish_from_scratch: clean test_nanopolish

clean:
	rm -rf ${TESTDIR} ${DATADIR}/reads ${DATADIR}/ref.fasta.*
