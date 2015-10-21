#!/bin/bash
for filename in inputs/sat/*; do
    ./Debug/monosat "$filename" 2> /tmp/monosat-regression-test | grep UNSAT
    if [[ $? -ne 1 ]]; then
	echo "REGRESSION PROBLEM"
	echo "File $filename is supposed to be SAT, but is UNSAT"
    fi
    cat /tmp/monosat-regression-test | grep -v "FPU to use double precision" # yeah, yeah, we know...
    if [[ $? -ne 1 ]]; then
	echo "REGRESSION PROBLEM"
	echo "File $filename violates an Assertion"
    fi
done
for filename in inputs/unsat/*; do
    ./Debug/monosat "$filename" 2> /tmp/monosat-regression-test | grep SATISFIABLE | grep -v UNSAT
    if [[ $? -ne 1 ]]; then
	echo "REGRESSION PROBLEM"
	echo "File $filename is supposed to be UNSAT, but is SAT"
    fi
    cat /tmp/monosat-regression-test | grep -v "FPU to use double precision"
    if [[ $? -ne 1 ]]; then
	echo "REGRESSION PROBLEM"
	echo "File $filename violates an Assertion"
    fi
done
