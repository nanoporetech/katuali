.PHONY: install install_scripts docs

venv: venv/bin/activate
IN_VENV=. ./venv/bin/activate
SNAKEFILE=${PWD}/Snakefile
JOBS=$(shell nproc)

DATADIR=test/data
TESTDIR=${DATADIR}/run
CONFIGDIR=${DATADIR}/configs
LOGDIR=${DATADIR}/logs

TEST=${IN_VENV} && cd test/data && katuali -s ${SNAKEFILE} --configfile config.yaml --printshellcmds -j ${JOBS} --config THREADS_PER_JOB=${JOBS}

OPT='--config SCRAPPIE_OPTS="raw -H mean --model rgrgr_r94 --local 10.0 --uuid --temperature1 0.65 --temperature2 1.7"'


venv/bin/activate:
	test -d venv || virtualenv venv --prompt '(katuali) ' --python=python3
	${IN_VENV} && pip install pip --upgrade
	${IN_VENV} && pip install -r requirements.txt


install: venv install_scripts

install_scripts:
	cp scripts/* venv/bin/

reads:
	cd test/data && mkdir -p run && cd run && tar -xf ../reads.tgz

test: install test_basecall test_align test_subsample test_canu test_racon test_medaka test_nanopolish test_miniasm_racon test_nanopolish_from_scratch check

check:
	grep '100%' ${LOGDIR}/*.log
	grep 'Exception' ${LOGDIR}/*.log || echo "No exceptions found"

test_basecall: clean reads
	${TEST} run/basecall/scrappie/basecalls.fasta

test_align: reads 
	${TEST} run/basecall/scrappie/align/calls2ref.bam

test_subsample: reads
	${TEST} run/basecall/scrappie/align/all_contigs/25X/basecalls.fasta

test_canu: reads
	${TEST} run/basecall/scrappie/canu_gsz_50k/consensus.fasta

test_racon: reads
	${TEST} run/basecall/scrappie/canu_gsz_50k/racon/consensus.fasta

test_medaka: reads
	${TEST} run/basecall/scrappie/canu_gsz_50k/racon/medaka/consensus.fasta

test_nanopolish: reads
	${TEST} run/basecall/scrappie/canu_gsz_50k/racon/nanopolish/consensus.fasta

test_miniasm_racon: reads
	${TEST} run/basecall/scrappie/miniasm_racon/consensus.fasta

test_nanopolish_from_scratch: clean test_nanopolish

clean:
	rm -rf ${TESTDIR} ${LOGDIR} ${CONFIGDIR} ${DATADIR}/reads ${DATADIR}/ref.fasta.*

update:
	git stash
	git checkout master
	git pull origin master
	${IN_VENV} && pip install -r requirements.txt


# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =
BUILDDIR      = _build

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

DOCSRC = docs

docs: venv
	${IN_VENV} && pip install sphinx sphinx_rtd_theme sphinx-argparse
	${IN_VENV} && cd $(DOCSRC) && $(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	rm -rf docs/modules.rst docs/medaka.rst  
	@echo
	@echo "Build finished. The HTML pages are in $(DOCSRC)/$(BUILDDIR)/html."
	touch $(DOCSRC)/$(BUILDDIR)/html/.nojekyll

