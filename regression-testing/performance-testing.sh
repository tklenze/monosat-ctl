#!/bin/bash
mode=0
echo $mode
for filename in sat-inputs/*; do
    echo $filename
# running these 20 times to make more significant time results
    /usr/bin/time -f "%E" bash -c '../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; ../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; exit 0' $filename $mode
done
for filename in unsat-inputs/*; do
    echo $filename
    /usr/bin/time -f "%E" bash -c '../Debug/monosat -use-symmetry-reduction=$1 $0 >/dev/null 2>&1; exit 0' $filename $mode
done
