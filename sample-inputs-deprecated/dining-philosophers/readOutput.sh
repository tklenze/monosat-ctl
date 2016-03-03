# Make solution human readable, by replacing the outbut bitvector
# T stands for: this fork is on the table
# L: This is taken by the left person
# R: This for is taken by the right person
# Note that a philosopher sitting in between two forks is eating when the left for has the state R and the right fork has the state L (modulo 5).
cat solution | sed 's/{100/L\ /g' | sed 's/{010/T\ /g' | sed 's/{001/R\ /g' | sed 's/\ 100/\ L\ /g' | sed 's/\ 010/\ T\ /g' | sed 's/\ 001/\ R\ /g' | sed 's/\ 100/\ L\ /g' | sed 's/\ 010/\ T\ /g' | sed 's/\ 001/\ R\ /g' | sed 's/\ 100/\ L\ /g' | sed 's/\ 010/\ T\ /g' | sed 's/\ 001/\ R\ /g' | sed 's/\ 100/\ L\ /g' | sed 's/\ 010/\ T\ /g' | sed 's/\ 001/\ R\ /g'
