#!/bin/bash
for filename in inputs/sat/*; do
    ./Debug/monosat "$filename" 2> /dev/null | grep UNSAT
    if [[ $? -ne 1 ]]; then
	echo "REGRESSION PROBLEM"
	echo "File $filename is supposed to be SAT, but is UNSAT"
    fi
done
for filename in inputs/unsat/*; do
    ./Debug/monosat "$filename" 2> /dev/null | grep SATISFIABLE | grep -v UNSAT
    if [[ $? -ne 1 ]]; then
	echo "REGRESSION PROBLEM"
	echo "File $filename is supposed to be UNSAT, but is SAT"
    fi
done
