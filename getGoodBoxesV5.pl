

#($filename,$prl,$prh,$gx,$gy,$gz,$rintval,$box) = ('TEST.pdbqt', 0.1, 0.9, 6, 6, 6, 1, 3);

#print ($filename,$prl,$prh,$gx,$gy,$gz,$rintval,$box);

#my @a = &getGoodBoxes($filename,$prl,$prh,$gx,$gy,$gz,$rintval,$box);
#foreach(@a){
#	print "$_\n";
#}


sub getGoodBoxes{
	my ($filename,$prl,$prh,$gx,$gy,$gz,$rintval,$box)=@_;
	($gx,$gy,$gz) = ($gx/2,$gy/2,$gz/2);
	my $maxdistan = 0;
	my $mindistan = 10000000000000000;
	my $continu = 5;
	my @neborNum;
	my @array_1;

	#Open PDB file
	open(FILE, $filename) || die "Can't open file $filename : $!\n";
	my @pdb = <FILE>;
	close(FILE);

	#Get Atom array
	my @Atom;
	foreach (@pdb){
		if($_ =~ /^ATOM/){
			if($_=~ /CA/){
			my ($num,$x,$y,$z)=(substr($_,6,5),substr($_,30,8),substr($_,38,8),substr($_,46,8));
			push @Atom, "$num $x $y $z";}
		}
	}

	#Generate Matrix by Atom array
	for(my $i=0;$i<=$#Atom;$i++){
		my ($iAtom,$ix,$iy,$iz)=split ' ',$Atom[$i];
		for(my $j=$i;$j<=$#Atom;$j++){
			my ($jAtom,$jx,$jy,$jz)=split ' ',$Atom[$j];
			my $d = getDistan ($ix, $iy, $iz, $jx, $jy, $jz);
			$array_1[$i][$j]=$d;
			$array_1[$j][$i]=$d;
			if($d<$mindistan){$mindistan=$d;}
			elsif($d>$maxdistan){$maxdistan=$d;}
		}
	}

	#Sort Matrix
	for(my $i=0;$i<=$#Atom;$i++){
		arraysort($array_1[$i]);
	}

	($mindistan<1)?($mindistan=1):($mindistan=int($mindistan));
	$maxdistan=int($maxdistan+0.99);

	my $max_sd;
	my $good_r = $mindistan;
	my $cont = $continu;
	for(my $r=$mindistan;$r<=$maxdistan;$r+=$rintval){
		my @newNeborNum;
		my @nb_num=();
		for(my $i=0;$i<=$#Atom;$i++){
			my ($iAtom,$ix,$iy,$iz)=split ' ',$Atom[$i];
			my $AtomNbNum = getaround($array_1[$i],$r);
			push @newNeborNum,"$iAtom $AtomNbNum $ix $iy $iz";
			push @nb_num,$AtomNbNum;
		}
		my ($avg,$sum1,$sum2)=(0,0,0);
		foreach my $temp(@newNeborNum){
			if($temp =~ /^\d+\s+(\d+)/){
				$sum1+=$1;
				$sum2+=($1**2);
			}else{
				die "Can't find nb_num in newNeborNum, Error line :$temp\n";
			}
		}

		$avg = $sum1 / ($#newNeborNum+1);
		my $now_sd = ((($sum2/($#newNeborNum+1))-($avg**2))**0.5);
	
		if($now_sd > $max_sd){
			$max_sd = $now_sd;
			$good_r = $r;
			@neborNum = ();
			@neborNum = @newNeborNum[sort {$nb_num[$a] <=> $nb_num[$b]} 0...$#newNeborNum];
			$cont = $continu; # should reset            by ken
		}
		else{
			$cont--;
			if($cont<=0){
				last;
			}
		}
	}
	my @tempNeborNum = @neborNum;
	my $lindex = int($prl*($#Atom+1)+0.99);
	my $hindex = int($prh*($#Atom+1));
	my @ans;
	for(my $i=0; $i<$box ; $i++){
		my $sent = getGoodSentinel ($lindex, $lindex, $hindex, $gx, $gy, $gz);
		($reN,$reG,$reX,$reY,$reZ) = split /\s+/,$neborNum[$sent];
		#push @ans, "$reX $reY $reZ";
		$ans[$i][0] = $reX; $ans[$i][1] = $reY; $ans[$i][2] = $reZ;
	}

	
	return @ans;
	
	sub arraysort{
		my $array = shift;
		@$array = sort {$a<=>$b} @$array;
	}

	sub getaround{
		my ($array,$num) = @_;
		if($$array[$#$array]<=$num){return $#$array;}
		my @t_array = sort {$a<=>$b}(@$array,$num);
		my ($low,$mid,$high)=(0,int(($#t_array+1)/2),$#t_array);
		while($t_array[$mid] != $num){
			($t_array[$mid]>$num)?(($low,$mid,$high)=($low,int(($low+$mid)/2),$mid)):(($low,$mid,$high)=($mid,int(($mid+$high)/2),$high));}
		while($t_array[$mid]==$t_array[$mid+1]){$mid+=1;}
		return $mid-1; # -1 for drop push $num forward
	}

	sub getGoodSentinel {
		
		my ($i,$lindex,$hindex,$gx,$gy,$gz) = @_;
		my ($iAtom,$in,$ix,$iy,$iz) = split ' ',$tempNeborNum[$i];
		my @notAtom = ();
		for(my $j=$lindex ; $j<=$hindex ; $j++){
			my ($jAtom,$jn,$jx,$jy,$jz) = split ' ',$tempNeborNum[$j];
			my ($x,$y,$z)=(abs($ix-$jx),abs($iy-$jy),abs($iz-$jz));
			if($x<$gx && $y<$gy && $z<$gz){
				if($in<$jn){
					return getGoodSentinel ($j, $lindex, $hindex, $gx, $gy, $gz);
				}
				else{
					push @notAtom, $j;
				}
			}
		}
		$tempNeborNum;
		deleteNotGoodSentinel (\@notAtom);
		return $i;
	}

	sub deleteNotGoodSentinel{
		my $nAtom = shift;
		my @nAtom = @$nAtom;
		
		foreach(@nAtom){
		
			my @line = split(' ',$tempNeborNum[$_]);
			splice @line,1,1,0;
			$tempNeborNum[$_]=join ' ',@line;
		}
	}

	sub getDistan {
		my ($ix,$iy,$iz,$jx,$jy,$jz)=@_;
		return (($ix-$jx)**2+($iy-$jy)**2+($iz-$jz)**2)**0.5;              
	}
}


1
