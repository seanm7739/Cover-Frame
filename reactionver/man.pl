
require "./connectboxes.pl";
require "./cfver0722.pl";
require "./cfsb.pl";

my @bp;
my @sp;
my @fp;
my @ans;
my @degree;

#connectboxes("2Z7X",0,0.4,13.52,27.58,13.20,9999,10,0.4,3.8,1,1,0,3);
#cf0722("2Z7X",0,0.4,13.52,27.58,13.20,9999,10,0.4,3.8,1,1,0,3);
cfsb("2Z7X",0,0.4,13.52,27.58,13.20,9999,10,0.4,3.8,1,1,0,3);
#cfsb("2Z7X",0,0.4,50,50.58,50.20,9999,10,0.4,3.8,1,1,0,3);

#print "@ans";
