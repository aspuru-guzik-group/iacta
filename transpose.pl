#!/usr/bin/perl
use strict;
use warnings;

sub distribute {
    my $fh = shift;
    my $concat = shift;    
    # Get the number of atoms
    my $mol = shift;
    while (1) {
	my $mol_index = sprintf("%4.4i", $mol);
	open(my $outfh, $concat, "mol$mol_index.xyz");
	my $line = readline($fh);
	my $natoms = int($line);
	print $outfh $line;

	for (my $i=0; $i<($natoms+1); $i++) {
	    my $line = readline($fh);
	    print $outfh $line;
	}

	close $outfh;
	$mol++;
	last if eof $fh;
    }
    return $mol;
}


my $mols = 0;
foreach my $mtd_index (@ARGV) {

    # Load the reactants file
    my $index_str = sprintf("%2.2i", $mtd_index);
    {
	my $react = "reactants_$index_str.xyz";
	open(my $fh, "<", $react) or die "can't open file $react, $!";
	distribute($fh, ">", $mols);
	close($fh);
    }

    my $popi = 0;
    while (1) {
	my $popindx = sprintf("%3.3i", $popi);
	my $fn = "prop$index_str" ."_$popindx" . ".xyz";
	open(my $fh, "<", $fn) or last;
	distribute($fh, ">>", $mols);
	close($fh);
	print $fn . "\n";
	$popi++;
    };

    {
	my $prod = "products_$index_str.xyz";
	open(my $fh, "<", $prod) or die "can't open file $prod, $!";
	$mols = distribute($fh, ">>", $mols);
	close($fh);
    }
}

`obabel reactants_* -o smiles -xn > reactants.smi`;
`obabel products_* -o smiles -xn > products.smi`;
