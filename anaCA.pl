


my $filename = shift;
my @p;
my $max=0;
my $num=0;
my @ca;

        open(QT, "$filename".".pdbqt") || die "Can't open file $filename : $!\n";
        my @pdb = <QT>;
        close(QT);


        foreach my $j(30...600){
		my $count=0;
		foreach $i(0...$#pdb*0.5){
                if($pdb[$i] =~ /^ATOM/){
			my ($n,$x,$y,$z)=(substr($pdb[$i],23,3),substr($pdb[$i],30,8),substr($pdb[$i],38,8),substr($pdb[$i],46,8));
                        if($pdb[$i]=~ /CA/ && $n == $j){
				$ca[0]=$x; $ca[1]=$y; $ca[2]=$z;
                        #print "$sp[$i][0]\n";
				#print "$n $x $y $z\n";
                        }
			elsif($n==$j){
				$p[$count][0]=$x; $p[$count][1]=$y; $p[$count][2]=$z;
				$count++;
			}
			elsif($n>$j){
				
				foreach $k(0...$#p){
					$temp=(($ca[0]-$p[$k][0])**2 + ($ca[1]-$p[$k][1])**2 + ($ca[2]-$p[$k][2])**2)**0.5;
					if($temp>$max){
						#print "ca $ca[0] $ca[1] $ca[2]  \n";
						#print "p $p[$k][0] $p[$k][1] $p[$k][2]\n";
						$max=$temp;
						print "$n $temp \n";
					}
				}
				@p=();
				$last;
			}
                }
		}
        }

