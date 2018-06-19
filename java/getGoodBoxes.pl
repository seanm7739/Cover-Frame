sub getGoodBoxes{
	my $filename = shift;
	my $prl = shift;
	my $prh = shift;
	my $gx = shift;
	my $gy = shift;
	my $gz = shift;
	my $boxnum = shift;
	my @a = `java GetGoodBoxes $filename $prl $prh $gx $gy $gz $boxnum`;

	return @a;
}

1