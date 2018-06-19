#c1 = centerx; c2 = centery; c3 = centerz; bs=boxsize 
#2014-9-11
my $c1 = shift;
my $c2 = shift;
my $c3 = shift;
my $bs = shift;
my $bp = $bs/2;
my $bn = $bp*-1;
print <<EOT;
Transform {
EOT

print "   translation $c1 $c2 $c3 \n";

print <<EOT;
   children [
Shape {
	appearance Appearance {
		material Material {
			transparency 0.01
		}
	}
	geometry IndexedLineSet {
		coord Coordinate {
			point [
EOT

print " 		 		$bn $bn $bp,	\n";
print "				$bp $bn $bp,	\n";
print "				$bp $bn $bn,	\n";
print "				$bn $bn $bn,	\n";			
print "				$bn $bp $bp,	\n";
print "				$bp $bp $bp,	\n";
print "				$bp $bp $bn,	\n";
print "				$bn $bp $bn	\n";

print <<EOT;
            ]
        }
		colorPerVertex FALSE 
		color Color {
			color [
                    1 0 0,
                    0 1 0,
                    0 0 1,
            ]
        }
		colorIndex [ 1 1 0 0 2 2 ] 
		coordIndex [
			 0, 3, 2, 1, -1,
			 4, 7, 6, 5, -1,
			 2, 3, 7, 6, -1,
			 0, 1, 5, 4, -1,
			 3, 7, 4, 0, -1,
			 1, 2, 6, 5, -1, 
		]	
	}
}
     ]
}
EOT
