#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )

if [ -z ${EDGE_HOME+x} ]; then
	EDGE_HOME="$rootdir/../../"
fi

test_result(){
	Test=$rootdir/TestOutput/SNP_Phylogeny/basesUsed.txt
	Expect=$rootdir/basesUsed.txt
	testName="EDGE Phylogenetic Analysis test";
	if cmp -s "$Test" "$Expect"
	then
		echo "$testName passed!"
		touch "$rootdir/TestOutput/test.success"
	else
		echo "$testName failed!"
		touch "$rootdir/TestOutput/test.fail"
	fi
}

cd $rootdir
echo "Working Dir: $rootdir";
echo "EDGE HOME Dir: $EDGE_HOME";

perl $EDGE_HOME/runPipeline -c $rootdir/config.txt -o $rootdir/TestOutput -cpu 4 -noColorLog  -p $rootdir/ebola_R1.fastq $rootdir/ebola_R2.fastq || true
rm -rf $rootdir/TestOutput/QcReads

test_result;