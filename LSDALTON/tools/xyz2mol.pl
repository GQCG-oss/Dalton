#!/usr/bin/perl -w
# A script for converting xyz input (Gaussian-type) to dalton mol input.
# First two lines are treated as a comment, all the rest as the atom data.
# Pawel Salek, 2006.
# Extend to convert Spartan "input" as well.
#
my %charge = ( 'H' => 1, 'C' => 6, 'N' => '7', 'O' => 8,
               'S' => 16 );

my $i;

die "Usage: xyz2mol.pl XYZ-FILE ...\n" unless $#ARGV >=0;
foreach $i (@ARGV) {
    unless(open FL, $i) {
	warn "cannot open $i for reading\n";
	next;
    }
    my $l1 = <FL>; chomp $l1;
    my $l2 = <FL>; chomp $l2;
    my $l3 = <FL>; chomp $l3;
    my ($totcharge, $spin) = split ' ', $l3;
    my %k = (); my %c = ();
    my $atoms = 0;
    while(<FL>){
        last if /^ENDCART/;
	s/^\s+//g;
	s/\s+$//g;
	($a) = ($_ =~ /^([a-zA-Z0-9]+)/);
        $k{$a} = '' unless defined $k{$a};
        $_ =~ s/^([a-zA-Z0-9]+)/$1  /g; # add spaces after atom name.
	$k{$a} = $k{$a} . "\n" . $_;
	$c{$a}++;
        $atoms++;
    }
    close FL;
    $ne = $i; $ne =~ s/\.xyz$/.mol/g;
    unless(open FL, ">$ne") {
	warn "cannot open $ne for writing.\n";
	next;
    }
    $at = scalar keys %c;
    print FL "ATOMBASIS\nLine1: $l1 Line2: $l2 Line3: $l3\n===========\n",
    "Atomtypes=$at Charge=$totcharge Angstrom\n";
    foreach $i (sort keys %k) {
        if($i>0) {
            $ch = 0+$i;
        } else {
            print "Atom '$i' unknown, file=$ne\n" unless  $charge{$i}>0;
            $ch = $charge{$i};
        }
	print FL "Charge=$ch Atoms=".$c{$i},
	" Basis=6-31G Aux=AhlrichsDenFit", $k{$i},"\n";
    }
    close FL;
    print "File: $i converted to $ne, $atoms atoms in total.\n";
}
