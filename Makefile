all: scripts/readToSNVList \
	SDAtest \
	ymls/python2.done \
	ymls/python3.done \
	TestCases/SDAtest/ref.fasta

SHELL := /bin/bash


# 
# Get readToSNVList
#
externalRepos/pbgreedyphase/readToSNVList:
	cd externalRepos/pbgreedyphase && make readToSNVList

scripts/readToSNVList: externalRepos/pbgreedyphase/readToSNVList
	cp externalRepos/pbgreedyphase/readToSNVList scripts/readToSNVList


#
# make conda envs
#
ymls/python2.done: 
	conda env create -f ymls/python2.yml && touch ymls/python2.done

ymls/python3.done:
	conda env create -f ymls/python3.yml && touch ymls/python3.done 


#
# Make the test case
#
TestCases/SDAtest/ref.fasta:
	rm -rf TestCases/SDAtest* 
	cd TestCases/ && \
	wget https://eichlerlab.gs.washington.edu/help/mvollger/SDA/SDAtest.tar.gz  && \
	tar -zxvf SDAtest.tar.gz


.PHONY: testcase, SDAtest
testcase:
	mkdir -p TestCases/testcase1/
	rm -rf TestCases/testcase1/* 
	mkdir -p TestCases/testcase1/config
	cd TestCases/testcase1/ && \
		wget https://eichlerlab.gs.washington.edu/help/mvollger/SDA/testcase1/ref.fasta  && \
		wget https://eichlerlab.gs.washington.edu/help/mvollger/SDA/testcase1/reads.fofn  && \
		wget https://eichlerlab.gs.washington.edu/help/mvollger/SDA/testcase1/reads.tar.gz  && \
		tar -zxvf reads.tar.gz && \
		cd config && \
		wget https://eichlerlab.gs.washington.edu/help/mvollger/SDA/testcase1/denovo.setup.config.json config/


