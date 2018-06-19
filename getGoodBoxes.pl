
sub getGoodBoxes{ #The main procedure: to get [$boxnum] good boxes above receptors [$fln]
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
	my @a = `java GetGoodBoxes $filename $prl $prh $gx $gy $gz $boxnum $lr $beta $alpha $rmodule $gbmod $mr $pushmod`;
	chomp @a;
	my @result;
	
	for(my $i=0 ; $i<$boxnum ; $i++){
		for(my $j=0 ; $j<3 ; $j++){
			$t = shift @a;
			if($t == "null"){
				last;
			}
			$result[$i][$j] = $t;
		}
	}
	$a = shift ;
	return @result;
}

1


