

sub connectboxes{ #The main procedure: to get [$boxnum] good boxes above receptors [$fln]
	my $filename = shift; 	#pdbqt file name
	my $prl = shift;	#neighbor percentile low
	my $prh = shift;	#neighbor percentile hight
	my $gx = shift;
	my $gy = shift;
	my $gz = shift;
	my $boxnum = shift;
	my $lr = shift;		#ligand radius
	my $alpha = shift;	
	my $beta = shift;
	my $rmodule = shift;	# 0=lr mod, 1=BestSd mod, 2=manual radius
	my $gbmod = shift;	# get prh nebor atom. 0=prh of full ball, 1=prh of all nebor of all atom
	my $mr = shift;		# manual radius
	my $pushmod = shift;	#pushmod 0=gamma% blueatmom, 1=number of puhmod blueatom 
	my @a = `java ConnectBoxes $filename $prl $prh $gx $gy $gz $boxnum $lr $beta $alpha $rmodule $gbmod $mr $pushmod`;
	chomp @a;
	my @x;
	my @y;
	my @z;
	my $g=$beta+1;
	$x[0]=$gx/2; $x[1]=$gx/2*-1;
	$y[0]=$gy/2; $y[1]=$gy/2*-1;
	$z[0]=$gz/2; $z[1]=$gz/2*-1;
	#my @bp; #blue point the atom inside the mol
	#my @sp; #silver point the atom on the surface
	#my @bkp; #black point the atom for candidate
	#my @fp; #first point

		$t = shift @a;
		#print "$t\n";
	my $bound = $#a+1;
	for(my $i=0 ; $i<$bound ; $i++){
		$t = shift @a;
		@ta = split(" ",$t);
		if($ta[0]=~/Surface/){
			#print "$ta[0]\n";
			last;
		}
		#if($t[0] == "null"){
		#	last;
		#}
		$bp[$i][0] = $ta[0];
		$bp[$i][1] = $ta[1];
		$bp[$i][2] = $ta[2];
	}
	$bound = $#a+1;
        for(my $i=0 ; $i<$bound ; $i++){
                $t = shift @a;
                #@ta = split(" ",$t);
                if($t=~/First/){
			#print "$ta[0]\n";
                        last;
                }
                #if($t[0] == "null"){
                #       last;
                #}
                $sp[$i][0] = $t;
                $sp[$i][1] = 0;  # 0=red 1=black 2=yellow 
        }
	$t = shift @a;
	@fp = split(" ",$t);   #get the first box
#print "$bp[0][0] $bp[0][1] $bp[0][2]\n";
#print "$sp[0][0] $sp[0][1] $sp[0][2]\n";
#print "$fp[0] $fp[1] $fp[2]\n";

############## initial	
	bomb($t,$g,$x[0],$x[1],$y[0],$y[1],$z[0],$z[1]);
	push(@ans,$t);
####################
	my $hasb=1;
	while($hasb==1){
		$hasb = neighbor_select($g,$x[0],$x[1],$y[0],$y[1],$z[0],$z[1]);
		if($hasb==0){
			$hasb = redcheck($g,$x[0],$x[1],$y[0],$y[1],$z[0],$z[1]);
		}
	}
	#foreach $str(@ans){
	#	print "$str \n";	
	#}
	#print "$#ans \n";
	#$a = shift ;
	print "@ans \n";
	pushans($pushmod, $lr, $alpha);
	
	print "@ans \n";
	my @result;
	foreach $i(0...$#ans){
		my @t=split(" ",$ans[$i]);
		$result[$i][0]=$t[0];  $result[$i][1]=$t[1];  $result[$i][2]=$t[2]; 
	}
#print "$#result";
#print "$result[0][0] $result[0][1] $result[0][2]";
	return @result;
}

sub bomb{
	my $bc = shift;
	my @x,@y,@z;
	my $dd = shift;
	$x[0]=shift; $x[1]=shift;
        $y[0]=shift; $y[1]=shift;
        $z[0]=shift; $z[1]=shift;

	my @ctr = split(" ",$bc);
	foreach $i(0...$#sp){
		@pt=split(" ",$sp[$i][0]);
		if(($sp[$i][1]!=2)&&(($pt[0]-$ctr[0]<$x[0])&&($pt[0]-$ctr[0]>$x[1]))&&(($pt[1]-$ctr[1]<$y[0])&&($pt[1]-$ctr[1]>$y[1]))&&(($pt[2]-$ctr[2]<$z[0])&&($pt[2]-$ctr[2]>$z[1])))
		{
			$sp[$i][1]=2;  #bombed 2
		}
	}
        $x[0]=$x[0]*$dd; $x[1]=$x[1]*$dd;
        $y[0]=$y[0]*$dd; $y[1]=$y[1]*$dd;
        $z[0]=$z[0]*$dd; $z[1]=$z[1]*$dd;

        foreach $i(0...$#sp){
                @pt=split(" ",$sp[$i][0]);
		if(($sp[$i][1]==0)&&(($pt[0]-$ctr[0]<$x[0])&&($pt[0]-$ctr[0]>$x[1]))&&(($pt[1]-$ctr[1]<$y[0])&&($pt[1]-$ctr[1]>$y[1]))&&(($pt[2]-$ctr[2]<$z[0])&&($pt[2]-$ctr[2]>$z[1])))
                {
                        $sp[$i][1]=1;  #new candidate 1
                }
        }
}




sub neighbor_select{

        my @x,@y,@z;
        my $dd = shift;
        $x[0]=shift; $x[1]=shift;
        $y[0]=shift; $y[1]=shift;
        $z[0]=shift; $z[1]=shift;
	my $max=-1;
	my $maxnum=-1;
	my @ctr;	
	my $haveb=0;

	foreach $i(0...$#sp){
		if($sp[$i][1]==1){
			my $count=0;
			@ctr = split(" ",$sp[$i][0]);
			foreach $j(0...$#sp){
				@pt=split(" ",$sp[$j][0]);
				if(($sp[$i][1]!=2)&&(($pt[0]-$ctr[0]<$x[0])&&($pt[0]-$ctr[0]>$x[1]))&&(($pt[1]-$ctr[1]<$y[0])&&($pt[1]-$ctr[1]>$y[1]))&&(($pt[2]-$ctr[2]<$z[0])&&($pt[2]-$ctr[2]>$z[1])))
               			{
					$count++;
		                }
			}
			if($count>$max){
				$max=$count;
				$maxnum=$i;
			}
		}
	}
	if($maxnum == -1){
		return 0;
	}
	else{
		push(@ans,$sp[$maxnum][0]);
		bomb($sp[$maxnum][0],$dd,$x[0],$x[1],$y[0],$y[1],$z[0],$z[1]);
		return 1;
	}


}

sub redcheck{

        my @x,@y,@z;
        my $dd = shift;
        $x[0]=shift; $x[1]=shift;
        $y[0]=shift; $y[1]=shift;
        $z[0]=shift; $z[1]=shift;
        my $max=-1;
        my $maxnum=-1;
        my @ctr;

        foreach $i(0...$#sp){
                if($sp[$i][1]==0){
                        my $count=0;
                        @ctr = split(" ",$sp[$i][0]);
                        foreach $j(0...$#sp){
                                @pt=split(" ",$sp[$j][0]);
                                if(($sp[$i][1]==0)&&(($pt[0]-$ctr[0]<$x[0])&&($pt[0]-$ctr[0]>$x[1]))&&(($pt[1]-$ctr[1]<$y[0])&&($pt[1]-$ctr[1]>$y[1]))&&(($pt[2]-$ctr[2]<$z[0])&&($pt[2]-$ctr[2]>$z[1])))
                                {
                                        $count++;
                                }
                        }
                        if($count>$max){
                                $max=$count;
                                $maxnum=$i;
                        }
                }
        }
        if($maxnum == -1){
                return 0;
        }
        else{
                push(@ans,$sp[$maxnum][0]);
                bomb($sp[$maxnum][0],$dd,$x[0],$x[1],$y[0],$y[1],$z[0],$z[1]);
                return 1;
        }






}


sub pushans{

        my $bn = shift;
        my $lr = shift;
        my $alpha = shift;
        my $l = $lr*$alpha;

        foreach $pt(0...$#ans){
                my @t = split(" ",$ans[$pt]);
                my $db = 0;
                my $min = 9999;
                my @sum;
				my $uv;
				$sum[0][0] = 0; $sum[0][1] = 0; $sum[0][2] = 0;
                for(my $i=0;$i<$bn;$i++){
						$min = 9999;
						foreach(@bp){
							my @tt = split(" ",$_);
							my $dt = (($t[0]-$tt[0])**2+($t[1]-$tt[1])**2+($t[2]-$tt[2])**2)**0.5;
							if( $dt<$min && $dt>$db ){
								$min = $dt;
								$sum[1][0] = $tt[0]; $sum[1][1] = $tt[1]; $sum[1][2] = $tt[2];
							}
						}
						$db = $min;
						$sum[0][0] = $sum[0][0] + $sum[1][0]; 
						$sum[0][1] = $sum[0][1] + $sum[1][1]; 
						$sum[0][2] = $sum[0][2] + $sum[1][2]; 
                }
				$sum[0][0] = $sum[0][0]/$bn; $sum[0][1] = $sum[0][1]/$bn; $sum[0][2] = $sum[0][2]/$bn;
				$sum[0][0] = $t[0]-$sum[0][0]; 
				$sum[0][1] = $t[1]-$sum[0][1]; 
				$sum[0][2] = $t[2]-$sum[0][2];
				$uv = ($sum[0][0]**2+$sum[0][1]**2+$sum[0][2]**2)**0.5;
				$sum[0][0] = $t[0]+($sum[0][0]/$uv)*$l; 
				$sum[0][1] = $t[1]+($sum[0][1]/$uv)*$l; 
				$sum[0][2] = $t[2]+($sum[0][2]/$uv)*$l;
				$ans[$pt] = "$sum[0][0] $sum[0][1] $sum[0][2]";
        }

}


















1


