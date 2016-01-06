# this is extremely hackish and requires some on-hand work to make the brackets match up afterwards

#just the formula, without and
cat MESIv1-abstract-spec | grep "(" | sed 's/M1/0/g' | sed 's/E1/1/g' | sed 's/S1/2/g' | sed 's/I1/3/g' | sed 's/M2/4/g' | sed 's/E2/5/g' | sed 's/S2/6/g' | sed 's/I2/7/g'

echo ""

#AND it together. Will require adding more brackets
cat MESIv1-abstract-spec | grep "(" | sed 's/M1/0/g' | sed 's/E1/1/g' | sed 's/S1/2/g' | sed 's/I1/3/g' | sed 's/M2/4/g' | sed 's/E2/5/g' | sed 's/S2/6/g' | sed 's/I2/7/g' | sed 's/^/AND\ (\ /g'
