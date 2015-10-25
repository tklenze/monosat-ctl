#!/bin/bash
for filename in inputs/sat/*; do
    bash -c './Debug/monosat $0 >/tmp/monosat-regression-test 2>/tmp/monosat-regression-test-errors' $filename >/tmp/monosat-regression-test-errors-bash 2>&1
    #./Debug/monosat "$filename" 2> /tmp/monosat-regression-test | grep UNSAT
    cat /tmp/monosat-regression-test* | grep UNSAT
    if [[ $? -ne 1 ]]; then
	echo "REGRESSION PROBLEM"
	echo "File $filename is supposed to be SAT, but is UNSAT"
    fi
    cat /tmp/monosat-regression-test-error* | grep -v "FPU to use double precision" # yeah, yeah, we know...
    if [[ $? -ne 1 ]]; then
	echo "REGRESSION PROBLEM"
	echo "File $filename has a problem"
    fi
    # Now, let's check the solution with NuSMV
    bash -c 'NuSMV inputConvertedToNuSMVInput.txt >/tmp/monosat-regression-test-nusmv 2>/tmp/monosat-regression-test-nusmv-errors' $filename >/tmp/monosat-regression-test-nusmv-errors-bash 2>&1
    cat /tmp/monosat-regression-test-nusmv* | grep Counterexample
    if [[ $? -ne 1 ]]; then
	echo "REGRESSION PROBLEM detected in NuSMV"
	echo "File $filename is supposed to be SAT, but is UNSAT in NuSMV"
    fi
    cat /tmp/monosat-regression-test-nusmv-error* | grep -v FOOBARBLUB
    if [[ $? -ne 1 ]]; then
	echo "REGRESSION PROBLEM detected in NuSMV"
	echo "File $filename has a problem in NuSMV"
    fi
done
for filename in inputs/unsat/*; do
    bash -c './Debug/monosat $0 >/tmp/monosat-regression-test 2>/tmp/monosat-regression-test-errors' $filename >/tmp/monosat-regression-test-errors-bash 2>&1
    cat /tmp/monosat-regression-test* | grep SATISFIABLE | grep -v UNSAT
    if [[ $? -ne 1 ]]; then
	echo "REGRESSION PROBLEM"
	echo "File $filename is supposed to be UNSAT, but is SAT"
    fi
    cat /tmp/monosat-regression-test-error* | grep -v "FPU to use double precision"
    if [[ $? -ne 1 ]]; then
	echo "REGRESSION PROBLEM"
	echo "File $filename has a problem"
    fi
done
