all:  externalRepos/pbgreedyphase/partitionByPhasedSNVs 

SHELL := /bin/bash

externalRepos/pbgreedyphase/partitionByPhasedSNVs:
	cd externalRepos/pbgreedyphase && make





.PHONY: testcase
testcase:
	mkdir TestCases/testcase1/
	mkdir TestCases/testcase1/config
	cd TestCases/testcase1/
	wget https://eichlerlab.gs.washington.edu/help/mvollger/SDA/testcase1/ref.fasta .
	wget https://eichlerlab.gs.washington.edu/help/mvollger/SDA/testcase1/reads.tar.gz .
	tar -zxvf reads.tar.gz
	wget https://eichlerlab.gs.washington.edu/help/mvollger/SDA/testcase1/denovo.setup.config.json config/.


