
#my %para; #parameters
my @para;
my $pc=0; #parameter count
#my @aa=["SEC","PYL","GLY","ALA","VAL","LEU","ILE","SER","THR","CYS","MET","PRO","LYS","ARG","PHE","TYR","TRP","HIS","ASP","GLU","ASN","GLN"];
#print "$aa[0][1] \n";
my %aacode=("1"=>"A","2"=>"R","3"=>"N","4"=>"D","5"=>"C","6"=>"Q","7"=>"E","8"=>"G","9"=>"H","10"=>"I","11"=>"L","12"=>"K","13"=>"M","14"=>"F","15"=>"P","16"=>"S","17"=>"T","18"=>"W","19"=>"Y","20"=>"V");


open FILE, "aaindexfav";
@temp=<FILE>;
close FILE;
#print @temp;
my $f=0;
print 'my @para=(';
#print "\n";

foreach $i(0...$#temp){
	if($temp[$i]=~/%%/){
		if($f==0){
			$f=1;
		}
        	else{
                	print ", ";
        	}
		@tmp1=split('\s+',$temp[$i+2]);
		@tmp2=split('\s+',$temp[$i+3]);
		#%para{aaa}[$pc] = "";
		print "{";
		foreach $j(0...9){
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














