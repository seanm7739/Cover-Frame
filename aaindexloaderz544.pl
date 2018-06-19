
#my %para; #parameters
my @para;
my $pc=0; #parameter count
#my @aa=["SEC","PYL","GLY","ALA","VAL","LEU","ILE","SER","THR","CYS","MET","PRO","LYS","ARG","PHE","TYR","TRP","HIS","ASP","GLU","ASN","GLN"];
#print "$aa[0][1] \n";
my %aacode=("1"=>"A","2"=>"R","3"=>"N","4"=>"D","5"=>"C","6"=>"Q","7"=>"E","8"=>"G","9"=>"H","10"=>"I","11"=>"L","12"=>"K","13"=>"M","14"=>"F","15"=>"P","16"=>"S","17"=>"T","18"=>"W","19"=>"Y","20"=>"V");


open FILE, "aaindex544";
@temp=<FILE>;
close FILE;
#print @temp;
my $f=0;
print 'my @para=(';
#print "\n";

foreach $i(0...$#temp){
	if($temp[$i]=~/^I    /){
		if($f==0){
			$f=1;
		}
        	else{
                	print ", ";
        	}
		@tmp1=split('\s+',$temp[$i+1]);
		@tmp2=split('\s+',$temp[$i+2]);
		#%para{aaa}[$pc] = "";
		my $sum = 0;
		for ( @tmp1 ) {
			$sum += $_;
		}
                for ( @tmp2 ) {
                        $sum += $_;
                }
		my $avg=$sum/20;
		my $std=0;
		for (1..$#tmp1){
			$std += ($tmp1[$_]-$avg)**2;
			$tt=($tmp1[$_]-$avg)**2;
			#print "\n $_ $avg  $tt\n";
		}
		for (1..$#tmp2){
			$std += ($tmp2[$_]-$avg)**2;
			$tt=($tmp2[$_]-$avg)**2;
                        #print "\n$tt\n";
		}
		$std = (($std/20))**0.5;
		for (1..$#tmp1){
                        $tmp1[$_] = ($tmp1[$_]-$avg)/$std;
                }
                for (1..$#tmp2){
                        $tmp2[$_] = ($tmp2[$_]-$avg)/$std;
                }

		#print "\n\n\n $avg $std \n\n\n";
		
		print "{";
		foreach $j(0..9){
			#$para[$pc][$j]="$tmp1[$j+1]";
			#$para[$pc][$j+10]="$tmp2[$j+1]";
			#print "\"$aacode{$j+1}\"=>\"$tmp1[$j+1]\", \"$aacode{$j+11}\"=>\"$tmp2[$j+1]\", ";
			if($j==9){
				#print "$tmp1[$j+1], $tmp2[$j+1]";
				print "\"$aacode{$j+1}\"=>\"$tmp1[$j+1]\", \"$aacode{$j+11}\"=>\"$tmp2[$j+1]\"";
			}
			elsif($j!=9)
			{
				print "\"$aacode{$j+1}\"=>\"$tmp1[$j+1]\", \"$aacode{$j+11}\"=>\"$tmp2[$j+1]\", ";
			}
		}
		#if($i!=$#temp){
		#	print "});\n";
		#}
		#else{
			print "}";
		#}
		$pc++;
		
				
	}




}

print ");\n";
print "print \"\$para[0]{A} \\n\";";

print "\n  $pc  \n";














