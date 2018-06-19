
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
#print "\n";

foreach $i(0...$#temp){
	if($temp[$i]=~/^H /){
		@tmp1=split('\s+',$temp[$i]);
		print "$tmp1[1] ";
	}
}
















