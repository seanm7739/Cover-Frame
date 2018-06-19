


require "./aaindexcz544.pl";
#require "./aaindexhash.pl";


my @bp;
my @sp;
my @fp;
my @ans;


#print "$para[0]{A}\n";
#exit;
#sub practice 

#input: $l=ligand $r=receptor output:pdbqt GPF DPF DLG 

#ex perl msub.pl 2LY4 2Z7X 1 0 0 0 0 0 0

#by seanm 2014-1-16



my $indir = shift;
my $l = shift;  #ligand
my $r = shift;  #receptor
my $spacing = shift; #spacing
my $mode = shift;      #0:normal 1:lh_finder 2:random center try 3:segment try
my $lows = shift;      #low start 0~1
#my $lowb = shift;      #low bound 0~1
my $lowb = $lows;      #low bound 0~1
my $highs = shift;     #high start 0~1
#my $highb = shift;     #high bound 0~1
my $highb = $highs;     #high bound 0~1
my $Boxnumber = shift; #boxnumber
my $al = shift;
my @gridxyz;
my $gnes = shift;	#ga_num_evals
my $rad = 5;
my $factor = shift; #factor
my $factor2 = shift; #factor
my $sw = shift;  #0 = lr ver; 1 = SD ver; 2 = dr ver
my $sw2 = shift; #0 = pr of max degree; 1 = ranking
my $sw3 = shift; #0 = %blueatom; 1~n = choose how much bule atom
my $dpr = shift; #r = pr of rec length (0~0.5) 
my $mp = shift;
my $dr = shift; #dr = the multi factor of dpr
my $pwrl = shift; #print wrl
my $candi = shift; #candiate box upper by ligand lr
my $bw = shift; #candiate weight
#print "$sw3  \n";
#system("cp Ligand/pdbqt/$l.pdbqt ../jp");
#print "$dr";
#system("cp Receptor/pdbqt/$r.pdbqt ../jp");

#print "$indir $l $r  $spacing $mode  $lows  $lowb $highs $Boxnumber $al $gnes  $factor  $factor2 $sw $sw2 $sw3  $dpr  $mp   $dr  $pwrl  $candi  $bw  \n";

#system("./pdb2pdbqt.sh $l $r");



#-----------debug message---------

#print "min x:$gridxyz[0][0] y:$gridxyz[1][0] z:$gridxyz[2][0]\n";

#print "max x:$gridxyz[0][1] y:$gridxyz[1][1] z:$gridxyz[2][1]\n";

#print "len x:$gridxyz[3][0] y:$gridxyz[3][1] z:$gridxyz[3][2]\n";

#print "center x:$gridxyz[4][0] y:$gridxyz[4][1] z:$gridxyz[4][2]\n";

#------------------------

#pass gridsize,gridcenter and spacing into sub Rallin

if($mode==0){

	@gridxyz = g_finder($l);  #call g_finder to ger a matrix about min max len and center
	my @temp=($gridxyz[3][0]*2, $gridxyz[3][1]*2, $gridxyz[3][2]*2);
	@gridxyz = g_finder($r);
	$gridxyz[3][0]=$gridxyz[3][0]+$temp[0]; $gridxyz[3][1]=$gridxyz[3][1]+$temp[1]; $gridxyz[3][2]=$gridxyz[3][2]+$temp[2];
	Rallin($indir,$gridxyz[3][0],$gridxyz[3][1],$gridxyz[3][2],$gridxyz[4][0],$gridxyz[4][1],$gridxyz[4][2],$spacing,$gnes);
	#system("mv svm $indir");
}

elsif($mode==1){
	#@gridxyz = g_finder($r);
	#$dr = ($gridxyz[3][0]+$gridxyz[3][1]+$gridxyz[3][2])/3;
	#$dr = $dr*$dpr;
	@gridxyz = g_finder($l);  #call g_finder to ger a matrix about min max len and center
	#$gridxyz[3][0]=$gridxyz[3][0]*$factor; $gridxyz[3][1]=$gridxyz[3][1]*$factor; $gridxyz[3][2]=$gridxyz[3][2]*$factor;
	lh_finder($indir,$l,$r,$gridxyz[3][0],$gridxyz[3][1],$gridxyz[3][2],$spacing,$lows,$lowb,$highs,$highb,$Boxnumber,$candi,$bw);
	#lh_finder($indir,$l,$r,$gridxyz[3][0],$gridxyz[3][1],$gridxyz[3][2],$spacing,$lows,$lowb,$highs,$highb,$Boxnumber);
	#system("mv svm $indir");
}

elsif($mode==2){

	#$gnes=int($gnes/$Boxnumber);

	@gridxyz = g_finder($l);  #call g_finder to ger a matrix about min max len and center

	$gridxyz[3][0]=$gridxyz[3][0]*$factor; $gridxyz[3][1]=$gridxyz[3][1]*$factor; $gridxyz[3][2]=$gridxyz[3][2]*$factor;

	r_try($l,$r,$gridxyz[3][0],$gridxyz[3][1],$gridxyz[3][2],$spacing,$Boxnumber);

}

elsif($mode==3){

	#$gnes=int($gnes/$Boxnumber);

	@gridxyz = g_finder($l);  #call g_finder to ger a matrix about min max len and center

	$gridxyz[3][0]=$gridxyz[3][0]*$factor; $gridxyz[3][1]=$gridxyz[3][1]*$factor; $gridxyz[3][2]=$gridxyz[3][2]*$factor;

	seg_try($l,$r,$gridxyz[3][0],$gridxyz[3][1],$gridxyz[3][2],$spacing,$Boxnumber);

}



sub Rallin{ #recepter in a big grid

	my $id = shift;
	my $gsx = shift;  #grid size x length
	my $gsy = shift;  #grid size y length
	my $gsz = shift;  #grid size z length
	my $gcx = shift;  #grid center x
	my $gcy = shift;  #grid center y
	my $gcz = shift;  #grid center z
	my $spacing = shift;        #angstrom
	my $gnes = shift; #prode num

	my @gridpts;            #number of xyz intervals
	my $factor = 0;

	#compute x y z intervals and use gauss ceil

	$gridpts[0]=int($gsx/$spacing+0.99);  

	$gridpts[1]=int($gsy/$spacing+0.99);

	$gridpts[2]=int($gsz/$spacing+0.99);

	#------------------------------------------------

	#----debug message------------

	#print "$gsx $gsy $gsz\n";

	#print "$gcx $gcy $gcz\n";

	#print "$spacing\n";	

	#print "$gridpts[0] $gridpts[1] $gridpts[2]\n";

	#----------------------------------

	#prepare GPF and DPF file

	#my $dir="$l"."_"."$r"."_Rallin";
	#my $dir="result";
	my $dir = $id;
	mkdir "$dir";
	mkdir "$dir/wrl";
	mkdir "$dir/wrl/1";
	mkdir "$dir/wrl/2";

	#system("cp rungrids_docking.sh $dir/");	

	#print"./prepare_gpf_dpf.sh $l $r $gridpts[0],$gridpts[1],$gridpts[2] $spacing $gcx,$gcy,$gcz all $spacing $dir $gnes \n";
	system("./prepare_gpf_dpf.sh $l $r $gridpts[0],$gridpts[1],$gridpts[2] $spacing $gcx,$gcy,$gcz all $spacing $dir $gnes");

	#system("./rungrids_docking.sh $l $r");

	#open FILE, ">>$dir/$l"."_"."$r"."_"."all.ctf"; #hope this could split on hadoop
	open FILE, ">>$dir/goodbox.ctf"; #hope this could split on hadoop

	print FILE "PUDOCK $l $r all $spacing\n";

	close FILE;

	#system("chmod 755 $l"."_"."$r"."_"."all.sh"); #hope this could split on hadoop

	system("cp $l.pdbqt $r.pdbqt $dir");

}



sub Rprin{ #receptor in several grid boxs and based on Rallin

#pass gridsize,gridcenter and spacing into sub Rallin



}



#sub getGoodBoxes{ #bacon's work: to find good probe method 

#Caller format: @arr = getGoodBoxes( $receptor, $low, $high, $Gsx, $Gsy, $Gsz,$Boxnumber)

#ex @arr = getGoodBoxes( 2Z7X, 0.3, 0.7, 10, 11, 12, 1);

# @arr contains $Boxnumber points (centers of good boxes)

#	my $receptor = shift;  #pdbqt file

#	my $low = shift;  #prl 0~1

#	my $high = shift; #prh 0~1

#	my $gsx = shift;  #grid size x

#	my $gsy = shift;  #grid size y

#	my $gsz = shift;  #grid size z

#	#Pseudo Call, bacon, please change the following codes:

#	return [(1, 2, 3)];

#}



sub g_finder{  #find a grid based on receptor and return its detail

	#find max grid size by receptor based

	my $pdbqt_name = shift;

	#my $mode = shift;  

	my $pdbqt;

	my @gridxyz; 

	#$gridxyz[axis][min_max]  axis = x=0 y=1 z=2 ;min_max = min=0 max=1

	#$gridxyz[3][0~2]=x y z length

	#$gridxyz[4][0~2]=x y z center point

	my $i=0;

	my $j=0;

	open FILE, "$pdbqt_name.pdbqt";

	@pdbqt=<FILE>;

	close FILE;



	#print @atom;

	#$grid[][];

	#print $#atom;



	#----initial the gridxyz box min_max----

	foreach $i(0..$#pdbqt){

		if(substr($pdbqt[$i],30,8)=~/\d.\d{3}$/){

			#----initial the gridxyz box min_max----

			$gridxyz[0][0]=substr($pdbqt[$i],30,8);

			$gridxyz[0][1]=$gridxyz[0][0];

			$gridxyz[1][0]=substr($pdbqt[$i],38,8);

			$gridxyz[1][1]=$gridxyz[1][0];

			$gridxyz[2][0]=substr($pdbqt[$i],46,8);

			$gridxyz[2][1]=$gridxyz[2][0];

			#-------------------------------------------------

			$j=$i+1;

			last;

		}

	}

	#-------------------------------------------------

	#print "$j\n";





	#----find the gridxyz box min_max----

	foreach $i($j..$#pdbqt){

		if(substr($pdbqt[$i],30,8)=~/\d.\d{3}$/){

			if(substr($pdbqt[$i],30,8)<$gridxyz[0][0]){

				$gridxyz[0][0]=substr($pdbqt[$i],30,8);

			}

			elsif(substr($pdbqt[$i],30,8)>$gridxyz[0][1]){

				$gridxyz[0][1]=substr($pdbqt[$i],30,8);

			}

			if(substr($pdbqt[$i],38,8)<$gridxyz[1][0]){

				$gridxyz[1][0]=substr($pdbqt[$i],38,8);

			}

			elsif(substr($pdbqt[$i],38,8)>$gridxyz[1][1]){

				$gridxyz[1][1]=substr($pdbqt[$i],38,8);

			}

			if(substr($pdbqt[$i],46,8)<$gridxyz[2][0]){

				$gridxyz[2][0]=substr($pdbqt[$i],46,8);

			}

			elsif(substr($pdbqt[$i],46,8)>$gridxyz[2][1]){

				$gridxyz[2][1]=substr($pdbqt[$i],46,8);

			}

		}

	}

	#-----------------------------------------------

	#----find the box's length of xyz-------------------

	$gridxyz[3][0]=$gridxyz[0][1]-$gridxyz[0][0];

	$gridxyz[3][1]=$gridxyz[1][1]-$gridxyz[1][0];

	$gridxyz[3][2]=$gridxyz[2][1]-$gridxyz[2][0];

	#-----------------------------------------

	#----find the box's center of xyz-------------------

	$gridxyz[4][0]=($gridxyz[0][1]+$gridxyz[0][0])/2;

	$gridxyz[4][1]=($gridxyz[1][1]+$gridxyz[1][0])/2;

	$gridxyz[4][2]=($gridxyz[2][1]+$gridxyz[2][0])/2;

	#--------debug----------------------------

	#print "min x:$gridxyz[0][0] y:$gridxyz[1][0] z:$gridxyz[2][0]\n";

	#print "max x:$gridxyz[0][1] y:$gridxyz[1][1] z:$gridxyz[2][1]\n";

	#print "len x:$gridxyz[3][0] y:$gridxyz[3][1] z:$gridxyz[3][2]\n";

	#print "center x:$gridxyz[4][0] y:$gridxyz[4][1] z:$gridxyz[4][2]\n";

	#---------------------------------------

	return @gridxyz;

}



sub lh_finder{  #test different lh to get better lh pair bacon method
	my $id = shift;        #in dir
	my $l = shift;         #ligand
	my $r = shift;         #receptor
	my $gsx = shift;       #grid size x
	my $gsy = shift;       #grid size y
	my $gsz = shift;       #grid size z
	my $spacing = shift;   #spacing
	my $lows = shift;      #low start
	my $lowb = shift;      #low bound
	my $highs = shift;     #high start
	my $highb = shift;     #high bound
	my $Boxnumber = shift;     #boxnumber
	my $candi = shift;	#candi box
	my $bw =shift;		#candiate weigt
	my @goodcenter;        #goodcenter
	my @gridpts;           #number of xyz intervals
	#my $dir="$l"."_"."$r"."_lh_finder";
	#my $dir = "result";
	my $dir = $id;
	mkdir "$dir";
        mkdir "$dir/wrl";
        mkdir "$dir/wrl/1";
        mkdir "$dir/wrl/2";
#print "hello \n";
	#system("cp rungrids_docking.sh $dir/");
	#compute x y z intervals and use gauss ceil
#print "$spacing $gsx $gsy $gsz\n";
        my $maxgs;
        if (($gsx<$gsy)&&($gsz<$gsy)){
                $maxgs=$gsy;
        }
        elsif (($gsy<$gsx)&&($gsz<$gsx)){
                $maxgs=$gsx;
        }
        elsif (($gsy<$gsz)&&($gsx<$gsz)){
                $maxgs=$gsz;
        }
	my $dpr = $maxgs;
	$dpr = $dpr*$dr;
        $gsx = $maxgs; $gsy = $maxgs; $gsz = $maxgs;

	my $lr=($gsx+$gsy+$gsz)/3;

		
	$gsx1=($gsx*$factor2+7.5)*$mp;
	$gsy1=($gsy*$factor2+7.5)*$mp;
	$gsz1=($gsz*$factor2+7.5)*$mp;
	$gsx=$gsx*$factor+7.5;
	$gsy=$gsy*$factor+7.5;
	$gsz=$gsz*$factor+7.5;
	$gridpts[0]=int($gsx1/$spacing+0.99);  
	$gridpts[1]=int($gsy1/$spacing+0.99);
	$gridpts[2]=int($gsz1/$spacing+0.99);
#print "$spacing $gridpts[0] $gridpts[1] $gridpts[2]\n";
	#system("./pdb2pdbqt.sh $l $r"); #prepare pdbqt file
	$lows = $lows*100;
	$lowb = $lowb*100;
	$highs = $highs*100;
	$highb = $highb*100;
	#open FI, ">>$dir/$l"."_"."$r"."_"."$lows"."_"."$lowb"."_"."$highs"."_"."$highb"."_finder.ctf";
	open FI, ">>$dir/goodbox.ctf";
my $wrl = "$r"."_"."$l";
if($pwrl == 1){
	system("perl pdb2vrml.pl $r.pdbqt > result/wrl/1/$wrl.wrl");
	system("perl pdb2vrml.pl $r.pdbqt > result/wrl/2/$wrl.wrl");
}
	foreach $low($lows...$lowb){
		foreach $high($highs...$highb){
			$prl = $low/100;
			$prh = $high/100;
#print "go into \n";
#print "$r,$prl,$prh,$gsx,$gsy,$gsz,$Boxnumber,$lr,$al,$factor,$sw,$sw2,$dpr,$sw3,$candi,$bw\n";
			@goodcenter = connectboxes($r,$prl,$prh,$gsx,$gsy,$gsz,$Boxnumber,$lr,$al,$factor,$sw,$sw2,$dpr,$sw3,$candi,$bw);
#print "go out\n";
			#print "$r,$prl,$prh,$gsx,$gsy,$gsz,$Boxnumber,$lr,$al,$factor\n";
			foreach $i(0...$#goodcenter){
				#system("sh prepare_gpf_dpf.sh $l $r $gridpts[0],$gridpts[1],$gridpts[2] $spacing $goodcenter[$i][0],$goodcenter[$i][1],$goodcenter[$i][2] $low $spacing"."_$high"."_$i"."_$factor"."_$factor2"."_$al"."_$gnes $dir $gnes");
				system("sh prepare_gpf_dpf.sh $l $r $gridpts[0],$gridpts[1],$gridpts[2] $spacing $goodcenter[$i][0],$goodcenter[$i][1],$goodcenter[$i][2] $low $spacing"."_$high"."_$i"."_$factor"."_$factor2"."_$dr"."_$gnes $dir $gnes"); #speed up for plot 2016/06/20 by seanm
				#print"sh prepare_gpf_dpf.sh $l $r $gridpts[0],$gridpts[1],$gridpts[2] $spacing $goodcenter[$i][0],$goodcenter[$i][1],$goodcenter[$i][2] $low $spacing"."_$high"."_$i"."_$factor"."_$factor2"."_$dr"."_$gnes $dir $gnes \n";
if($pwrl == 1){
				system("perl addbox.pl $goodcenter[$i][0] $goodcenter[$i][1] $goodcenter[$i][2] $gsx >> result/wrl/1/$wrl.wrl");
				system("perl addbox.pl $goodcenter[$i][0] $goodcenter[$i][1] $goodcenter[$i][2] $gsx1 >> result/wrl/2/$wrl.wrl");
}
				#open FILE, ">>$l_$r_$lows_$lowb_$highs_$highb.sh";
				#print FI "PUDOCK $l $r $low $spacing"."_$high"."_$i"."_$factor"."_$factor2"."_$al"."_$gnes\n";
				print FI "PUDOCK $l $r $low $spacing"."_$high"."_$i"."_$factor"."_$factor2"."_$dr"."_$gnes\n";
			}
		}
	}
	close FI;
	#system("chmod 755 $l"."_"."$r"."_"."$lows"."_"."$lowb"."_"."$highs"."_"."$highb.sh"); #hope this could split on hadoop
	system("cp $l.pdbqt $r.pdbqt $dir");
}


sub seg_try{  #segment center try

	my $l = shift;         #ligand

	my $r = shift;         #receptor

	my $gsx = shift;       #grid size x

	my $gsy = shift;       #grid size y

	my $gsz = shift;       #grid size z

	my $spacing = shift;   #spacing

	my $Boxnumber = shift; #Boxnumber

	my @center;            #gridcenter xyz

	my $cut = $Boxnumber+1;

	my @gridpts;           #number of xyz intervals

	my $dir="$l"."_"."$r"."_seg_try";

	mkdir "$dir";

	#system("cp rungrids_docking@gridxyz = g_finder($r);.sh $dir/");

	#compute x y z intervals and use gauss ceil

	$gridpts[0]=int($gsx/$spacing+0.99);

	$gridpts[1]=int($gsy/$spacing+0.99);

	$gridpts[2]=int($gsz/$spacing+0.99);

	open FILE, "$r.pdbqt";

	@tmp=<FILE>;

	close FILE;

	#my $rand = int (rand $#tmp-2);

	my $peace = $#tmp/$cut;  #place range

	my $point = 0;

	open FILE, ">$dir/$l"."_"."$r"."_"."seg.ctf"; #hope this could split on hadoop

	foreach $i(0...$Boxnumber-1){

		$point=$point+$peace;

		if(substr($tmp[$point],30,8)=~/\d.\d{3}$/){

			$center[$i][0]=int(substr($tmp[$point],30,8));

			$center[$i][1]=int(substr($tmp[$point],38,8));

			$center[$i][2]=int(substr($tmp[$point],46,8));

		}

		else{

			$center[$i][0]=int(substr($tmp[$point+1],30,8));

			$center[$i][1]=int(substr($tmp[$point+1],38,8));

			$center[$i][2]=int(substr($tmp[$point+1],46,8));

		}

		system("./prepare_gpf_dpf.sh $l $r $gridpts[0],$gridpts[1],$gridpts[2] $spacing $center[$i][0],$center[$i][1],$center[$i][2] seg $i $dir $gnes");

		print FILE "PUDOCK $l $r seg $i\n";

	}

	close FILE;

	#system("chmod 755 $l"."_"."$r"."_"."rand.sh"); #hope this could split on hadoop

	#system("./$l"."_"."$r"."_"."rand.sh"); #hope this could split on hadoop

	#system("cp $l.pdbqt $r.pdbqt $dir");

}



sub r_try{  #random center try

	my $l = shift;         #ligand

	my $r = shift;         #receptor

	my $gsx = shift;       #grid size x

	my $gsy = shift;       #grid size y

	my $gsz = shift;       #grid size z

	my $spacing = shift;   #spacing

	my $Boxnumber = shift; #Boxnumber

	my @center;            #gridcenter xyz

	#my $rtimes = $Boxnumber+1;

	my @gridpts;           #number of xyz intervals

	my $dir="$l"."_"."$r"."_r_try";

	mkdir "$dir";

	#system("cp rungrids_docking.sh $dir/");

	#compute x y z intervals and use gauss ceil

	$gridpts[0]=int($gsx/$spacing+0.99);

	$gridpts[1]=int($gsy/$spacing+0.99);

	$gridpts[2]=int($gsz/$spacing+0.99);

	open FILE, "$r.pdbqt";

	@tmp=<FILE>;

	close FILE;

	my $rand = int (rand $#tmp-2);

	#my $peace = $#tmp/$cut;  #place range

	#my $point = 0;

	open FILE, ">$dir/$l"."_"."$r"."_"."rand.ctf"; #hope this could split on hadoop

	foreach $i(0...$Boxnumber-1){

		#$point=$point+$peace;

		my $rand = int (rand $#tmp-2);

		if(substr($tmp[$point],30,8)=~/\d.\d{3}$/){

			$center[$i][0]=int(substr($tmp[$rand],30,8));

			$center[$i][1]=int(substr($tmp[$rand],38,8));

			$center[$i][2]=int(substr($tmp[$rand],46,8));

		}

		else{

			$center[$i][0]=int(substr($tmp[$rand+1],30,8));

			$center[$i][1]=int(substr($tmp[$rand+1],38,8));

			$center[$i][2]=int(substr($tmp[$rand+1],46,8));

		}

		system("./prepare_gpf_dpf.sh $l $r $gridpts[0],$gridpts[1],$gridpts[2] $spacing $center[$i][0],$center[$i][1],$center[$i][2] rand $i $dir $gnes");

		print FILE "PUDOCK $l $r rand $i\n";

	}

	close FILE;

	#system("chmod 755 $l"."_"."$r"."_"."rand.sh"); #hope this could split on hadoop

	#system("./$l"."_"."$r"."_"."rand.sh"); #hope this could split on hadoop

	#system("cp $l.pdbqt $r.pdbqt $dir");

}
