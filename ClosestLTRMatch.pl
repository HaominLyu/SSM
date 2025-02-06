#!perl -w
use strict;

# perl this_code.pl *.genome.fa.out.LAI.LTR.ava.out > *.genome.fa.out.LAI.LTR.ava.out.closest

my %done;
while (<>) {
	my ($ltr1, $ltr2) = split /\t/;
	next if $ltr1 eq $ltr2;
	next if defined $done{$ltr1};
	print and $done{$ltr1} = $ltr2;
}


