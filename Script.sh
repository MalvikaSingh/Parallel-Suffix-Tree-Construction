#! /bin/bash
#! /usr/bin/gnuplot
clear
gcc -fopenmp Serial.c -lm -o output_s			#Serial object file
gcc -fopenmp Parallel.c -lm -o output_p			#Parallel object file

for ((i=3;i<=9;i++))
do
	j=$((10**$i))					
	for((k=1;k<=3;k++))			
	do						#Loop for number of cores
		ts="$(/home/201401098/output_s $j)"	#Specify the path of output files		
		tp="$(/home/201401098/output_p $j $k)"
		echo $k\	$(echo "$ts/$tp"|bc -l|awk '{printf "%f", $0}')>>$(echo file.$i.txt)	#For different cores vs speedup
	done
	ts="$(/home/201401098/output_s $j)"		#For 4 cores, serial and parallel time and speedup
	echo $i\	$ts>>Serial_sizevt.txt
	tp="$(/home/201401098/output_p $j 4)"
	echo $i\	$tp>>Parallel_sizevt.txt
	echo $i\	$(echo "$ts/$tp"|bc -l|awk '{printf "%f", $0}')>>Speedup.txt	#Calculating speedup
	echo 4\	$(echo "$ts/$tp"|bc -l|awk '{printf "%f", $0}')>>$(echo file.$i.txt)
	for((k=5;k<=12;k++))		
	do						#Loop for 5 to 12 cores
		ts="$(/home/201401098/output_s $j)"	#Specify the path		
		tp="$(/home/201401098/output_p $j $k)"
		echo $k\	$(echo "$ts/$tp"|bc -l|awk '{printf "%f", $0}')>>$(echo file.$i.txt)	#Files containing cores v/s speedup
	done
done
exit 0
