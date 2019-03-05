.PHONY: install update reads clean_test test_basecall test_align test_subsample test_canu test_racon test_medaka test_nanopolish test_miniasm_racon test_nanopolish_from_scratch check

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
	JOBS=$(shell sysctl -n hw.physicalcpu)
endif
ifeq ($(UNAME), Linux)
	JOBS=$(shell nproc)
endif

IN_VENV=. ./venv/bin/activate
TEST=${IN_VENV} && cd test && katuali -s Snakefile --configfile config.yaml --printshellcmds -j ${JOBS} --config THREADS_PER_JOB=${JOBS}
OPT='--config SCRAPPIE_OPTS="raw -H mean --model rgrgr_r94 --local 10.0 --uuid --temperature1 0.65 --temperature2 1.7"'

venv/bin/activate:
	test -d venv || virtualenv venv --prompt '(katuali) ' --python=python3
	${IN_VENV} && pip install pip --upgrade
	${IN_VENV} && pip install -r requirements.txt

install: venv/bin/activate
	${IN_VENV} && python setup.py install

update: venv/bin/activate
	git stash
	git checkout master
	git pull origin master
	${IN_VENV} && pip install -r requirements.txt


test: install clean_test test_basecall test_align test_subsample test_canu test_racon test_medaka test_nanopolish test_miniasm_racon test_nanopolish_from_scratch check

test/config.yaml:
	mkdir -p test
	${IN_VENV} && katuali_config $@ 

test/Snakefile: test/config.yaml
	${IN_VENV} && katuali_datafile Snakefile | xargs -I {} cp {} $@ 

reads:
	mkdir -p test/run
	${IN_VENV} && cd test/run && katuali_datafile test/reads.tgz | xargs -I {} tar -xf {}


# The following targets step through a pipeline

test_basecall: reads test/config.yaml test/Snakefile
	${TEST} run/basecall/scrappie/basecalls.fasta

test_align: reads test/config.yaml test/Snakefile
	${TEST} run/basecall/scrappie/align/calls2ref.bam

test_subsample: reads test/config.yaml test/Snakefile
	${TEST} run/basecall/scrappie/align/all_contigs/25X/basecalls.fasta

test_canu: reads test/config.yaml test/Snakefile
	${TEST} run/basecall/scrappie/canu_gsz_50k/consensus.fasta

test_racon: reads test/config.yaml test/Snakefile
	${TEST} run/basecall/scrappie/canu_gsz_50k/racon/consensus.fasta

test_medaka: reads test/config.yaml test/Snakefile
	${TEST} run/basecall/scrappie/canu_gsz_50k/racon/medaka/consensus.fasta

test_nanopolish: reads test/config.yaml test/Snakefile
	${TEST} run/basecall/scrappie/canu_gsz_50k/racon/nanopolish/consensus.fasta

test_miniasm_racon: reads test/config.yaml test/Snakefile
	${TEST} run/basecall/scrappie/miniasm_racon/consensus.fasta

test_nanopolish_from_scratch: clean_test test_nanopolish test/config.yaml test/Snakefile

clean_test:
	rm -rf test

check:
	grep '100%' test/logs/*.log
	grep 'Exception' test/logs/*.log || echo "No exceptions found"


# Build docs
SPHINXOPTS      =
SPHINXBUILD     = sphinx-build
PAPER           =
BUILDDIR        = _build
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .
DOCSRC          = docs

docs: venv/bin/activate 
	${IN_VENV} && pip install sphinx sphinx_rtd_theme sphinx-argparse
	${IN_VENV} && cd $(DOCSRC) && $(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	rm -rf docs/modules.rst docs/medaka.rst  
	@echo
	@echo "Build finished. The HTML pages are in $(DOCSRC)/$(BUILDDIR)/html."
	touch $(DOCSRC)/$(BUILDDIR)/html/.nojekyll

