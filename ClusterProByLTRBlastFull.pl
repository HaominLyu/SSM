#!perl -w
use strict;


open CTL,"<","$ARGV[0]" or die;
my (@uplimit, @downlimit, @chrline);
while (<CTL>){
	chomp;
	next if /^#/;
	my @line = split /:\s+/;
	@chrline = split /\s*,\s*/,$line[1] if $line[0] eq 'Chromosome';
	@uplimit = split /\s*,\s*/,$line[1] if $line[0] eq 'UpLimit';
	@downlimit = split /\s*,\s*/,$line[1] if $line[0] eq 'DownLimit';
}
close CTL;


die "Check your ctl file: Please use correct format and input more than two chromosomes!" if scalar @chrline <= 1;
die "Check your ctl file: Please use correct format and input more than one up limit!" if scalar @uplimit <= 0;
die "Check your ctl file: Please use correct format and input more than one down limit!" if scalar @downlimit <= 0;
die "Check your ctl file: Please input the same numbers of up and down limits!" if (scalar @uplimit) != (scalar @downlimit);


my %chr;
for (@chrline) { $chr{$_} = 1; }
my %limit;
for (1..(scalar @uplimit)) { ${$limit{$_}}{'up'} = $uplimit[$_-1]; }
for (1..(scalar @downlimit)) { ${$limit{$_}}{'down'} = $downlimit[$_-1]; }


open LAI,"<","$ARGV[1]" or die;
my %count;
my %total;
while (<LAI>) {
	my ($name1, $name2, $ident) = (split /\t/)[0..2];
	$name1 =~ s/:.*//;
	$name2 =~ s/:.*//;
	next unless defined $chr{$name1} and defined $chr{$name2};
	my @name = sort ($name1, $name2);
	my $name = join(":",@name);
	$total{$name} ++;

	for my $no (keys %limit) {
		if ($ident <= ${$limit{$no}}{'up'} and $ident >= ${$limit{$no}}{'down'}) {
			my $label = "${$limit{$no}}{'down'}".".To."."${$limit{$no}}{'up'}";
			${$count{$label}}{$name} ++;
		}
	}
}
close LAI;


my @allchr = sort keys %chr;
my $outhead = join("\t",@allchr);


for my $outident (keys %count) {
	open OUT,">>","SSM.$outident.proportion" or die;
	print OUT "\t$outhead\n";
	for my $outchr1 (@allchr) {
		print OUT "$outchr1";
		for my $outchr2 (@allchr) {
			my @outchr = sort ($outchr1, $outchr2);
			my $outchr = join(":", @outchr);
			my $num = 0;
			$num = ${$count{$outident}}{$outchr} if defined ${$count{$outident}}{$outchr};
			my $outpro = $num / $total{$outchr};
			print OUT "\t$outpro";
		}
		print OUT "\n";
	}
	close OUT;
}


