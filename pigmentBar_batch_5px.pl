#!/usr/bin/perl


##############################
#
#For bar pigment phenotypes in parental and hybrid offspring of Aulonocara koningsii and Metriaclima mbenjii
#
#inputs a directory with .csv files outputted from ImageJ 
#
#for each of the plot profiles for the bar, calculate
#
#average darkness of the bar
#average darkness of non bar regions 
#overall % bar (len sum of bars/total Len)
#number of bars
#average width of bars
#
#to count as a bar, needs to be at least 5 pixels in width
#if there's a break of 5 pixels that are lighter, counts this as a non-bar
#if <5 pixels in a row that are dark, counts as non-bar
##############################


if (@ARGV != 2) {
    die "\nUsage: perl pigmentBar_update.pl DataDirectory BarOutputResults.txt\n\n";
}

$directory = shift @ARGV;		
$out1 = shift @ARGV;				



opendir(DIR, $directory) or die ("Couldn't open directory $directory\n");
@files = readdir(DIR);
closedir(DIR);


open (OUT1, ">$out1");
print OUT1 "ID\tNumberBars\tAveWidthBar\tAveIntensityBar\tAveWidthNonbar\tAveIntensityNonbar\tPctBars\n";

foreach $filename (@files) {
	@all = ();

	#############################
	#read in file from directory
	#data from cdv file goes into 2 by n array @all with format [pixel#][intensity]
	#############################

	if($filename =~ /.csv$/){
		
		@name = split /_/, $filename;
		#format: species_ID_cold_CB_layer#_barname.jpg

		print OUT1 "$name[1]\t";

		$count=0;
		$totalintensity=$average=0;
		open (IN, "<$directory/$filename") or die ("Couldn't open file $filename\n");
		$line=<IN>;


		while ($line=<IN>){
			@temp=split (',', $line);
			$count++;
			$all[$count][0]=$count; #each pixel is 0.002857 inches
			$all[$count][1]=$temp[1];
			$totalintensity=$totalintensity+$temp[1];
		}


		#calculate average intensity to use as a cutoff
		$average=$totalintensity/$count;


		#############################
		#go through entire array, combining series of bar (pigment <=$average) or non-bar (pigment >$average) into one element
		#put into 3 by n array @combo with format [pixel_length][sum_intensity][1=bar,0=nonbar]
		#############################

		$bpcount=$nbpcount=$barrunintensity=$nonbarrunintensity=$combolines=0;
		@combo=();

		for ($m=1; $m<=$count; $m++){ 

			if($all[$m][1]<=$average){

				$bpcount++;
				$barrunintensity =$barrunintensity +$all[$m][1];
		
				#print non-bar info to @combo array
				if($nbpcount>0){
					$combo[$combolines][0]=$nbpcount;
					$combo[$combolines][1]=$nonbarrunintensity;
					$combo[$combolines][2]=0;
					$nbpcount=$nonbarrunintensity =0;
					$combolines++;
				}	
			}


			if($all[$m][1]>$average){
				$nbpcount ++;
				$nonbarrunintensity =$nonbarrunintensity +$all[$m][1];
	
				#print bar info to @combo array
				if($bpcount>0){
					$combo[$combolines][0]=$bpcount;
					$combo[$combolines][1]=$barrunintensity;
					$combo[$combolines][2]=1;
					$bpcount=$barrunintensity =0;
					$combolines++;
				}	
			}
		}

		#last line needs to be printed

		if($nbpcount>0){
			$combo[$combolines][0]=$nbpcount;
			$combo[$combolines][1]=$nonbarrunintensity;
			$combo[$combolines][2]=0;
			$combolines++;
		}

		if($bpcount>0){
			$combo[$combolines][0]=$bpcount;
			$combo[$combolines][1]=$barrunintensity;
			$combo[$combolines][2]=1;
			$combolines++;
		}	


		#############################
		#go through combo array
		#bar must be >=5 pixels to count as a bar
		#if bar is interrupted by >=5 pixels, then this is a non-bar
		#any dark pixels <=5 wide count as non-bar		
		#put into 3 by n array @final with format [pixel_length][sum_intensity][1=bar,0=nonbar]
		#############################

		$pixelcount2=$intensity2=$finallines=0;
		@final = ();
		$flag=-1; #don't know if it's bar (1) or non bar (0) to start

		for ($n=0; $n<$combolines; $n++){ 
			#determine if we start with a bar or non-bar
			if ($flag==-1){
				$flag=$combo[$n][2];
			}	

			if($flag==$combo[$n][2]){
				$pixelcount2=$pixelcount2+$combo[$n][0];
				$intensity2=$intensity2+ $combo[$n][1];

			}

			#if new row doesn't match, check to see if it's short (5 or less pixels) and should be combined
			if ($flag!=$combo[$n][2]){
				if ($combo[$n][0]<6){
					$pixelcount2=$pixelcount2+$combo[$n][0];
					$intensity2=$intensity2+ $combo[$n][1];
				}
			
				#if it's long enough, then 
				else{
					$final[$finallines][0]=$pixelcount2;
					$final[$finallines][1]=$intensity2;
					$final[$finallines][2]=$flag;
					$finallines++;
					$pixelcount2=$intensity2=0;
					$pixelcount2=$pixelcount2+$combo[$n][0];
					$intensity2=$intensity2+ $combo[$n][1];
					$flag=$combo[$n][2];
				}
			}
		}

		#print last line
		if ($pixelcount2>0){
			$final[$finallines][0]=$pixelcount2;
			$final[$finallines][1]=$intensity2;
			$final[$finallines][2]=$flag;
			$finallines++;
		}


		#############################
		#go through final array @final with format [pixel_length][sum_intensity][1=bar,0=nonbar]
		#calculate:
		#1. Total number of bars
		#2. Ave width of bars
		#3. Ave intensity of bars
		#4. Ave width non bars
		#5. Ave intensity of non bars
		#6. % of length that is barred
		#############################

		$numberbars=$numbernonbars=$avewidthbar=$avewidthnonbar=$aveintensitybar=$aveintensitynonbar=$pctbars=0;
		$barintensitysum=$nonbarintensitysum=$barpixels=$nonbarpixels=0;

		for ($t=0; $t<$finallines; $t++){ 
			#if bar, add to bar running tallies

			if($final[$t][2]==1){
				$numberbars++;
				$barpixels=$barpixels+$final[$t][0];
				$barintensitysum=$barintensitysum+$final[$t][1];

			}
			else{
				$numbernonbars++;
				$nonbarpixels=$nonbarpixels+$final[$t][0];
				$nonbarintensitysum=$nonbarintensitysum+$final[$t][1];

			}
		}

		if ($barpixels>0){
			$avewidthbar=$barpixels/$numberbars;
			$aveintensitybar=$barintensitysum/$barpixels;
			$totalpixels=$barpixels+$nonbarpixels;
			$pctbars=$barpixels/$totalpixels;
		}
		
		if ($nonbarpixels>0){
			$avewidthnonbar=$nonbarpixels/$numbernonbars;
			$aveintensitynonbar=$nonbarintensitysum/$nonbarpixels;
		}


		print OUT1 "$numberbars\t$avewidthbar\t$aveintensitybar\t$avewidthnonbar\t$aveintensitynonbar\t$pctbars\n";

	close IN;
	} #if .csv
} #foreach file

close OUT1;
exit;