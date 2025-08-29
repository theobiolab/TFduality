n=1000000

declare -a names=("constraints_binding" "constraints_chr")
declare -a seeds=(1 20)
len=${#names[@]}


i=0
seed=${seeds[$i]}
name=${names[$i]}
seed2=$((seed+1))
seed3=$((seed+2))
seed4=$((seed+3))
seed5=$((seed+4))
seed6=$((seed+5))
seed7=$((seed+6))
seed8=$((seed+7))
seed9=$((seed+8))
seed10=$((seed+9))
echo "index $i, name $name, seed $seed, seed2 $seed2, seed3 $seed3, seed4 $seed4 seed5 $seed5"
../cfO2/enumConvexPolytopeVertices.o ${name}.txt > ${name}.vertices
../cfO2/sampleFromConvexPolytopeRMC.o ${name}.vertices $n $seed> x1x2x3x4_inc_1.points
echo 2
../cfO2/sampleFromConvexPolytopeRMC.o ${name}.vertices $n $seed2> x1x2x3x4_inc_2.points
echo 3
../cfO2/sampleFromConvexPolytopeRMC.o ${name}.vertices $n $seed3> x1x2x3x4_inc_3.points
echo 4
../cfO2/sampleFromConvexPolytopeRMC.o ${name}.vertices $n $seed4> x1x2x3x4_inc_4.points
echo 5
../cfO2/sampleFromConvexPolytopeRMC.o ${name}.vertices $n $seed5> x1x2x3x4_inc_5.points
echo 6
../cfO2/sampleFromConvexPolytopeRMC.o ${name}.vertices $n $seed6> x1x2x3x4_inc_6.points
echo 7
../cfO2/sampleFromConvexPolytopeRMC.o ${name}.vertices $n $seed7> x1x2x3x4_inc_7.points


i=1
seed=${seeds[$i]}
name=${names[$i]}
seed2=$((seed+1))
seed3=$((seed+2))
seed4=$((seed+3))
seed5=$((seed+4))
echo "index $i, name $name, seed $seed, seed2 $seed2, seed3 $seed3" 
../cfO2/enumConvexPolytopeVertices.o ${name}.txt > ${name}.vertices
../cfO2/sampleFromConvexPolytopeRMC.o ${name}.vertices $n $seed> x1x2x3x4_relative_1.points
echo 2
../cfO2/sampleFromConvexPolytopeRMC.o ${name}.vertices $n $seed2> x1x2x3x4_relative_2.points
echo 3
../cfO2/sampleFromConvexPolytopeRMC.o ${name}.vertices $n $seed3> x1x2x3x4_relative_3.points
echo 4
../cfO2/sampleFromConvexPolytopeRMC.o ${name}.vertices $n $seed4> x1x2x3x4_relative_4.points

