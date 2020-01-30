#!/usr/bin/perl

use strict;
use warnings;

print STDERR "Enter name of bed file:\n";
my $bedfile = <STDIN>;
chomp $bedfile;
open (BED, "$bedfile") || die "Could not open .bed file: $bedfile\n";

print STDERR "Enter name of matching mutation sequence file:\n";
my $nucfile = <STDIN>;
chomp $nucfile;
open (NUC, "$nucfile") || die "Could not open .fa file: $nucfile\n";

my $head = "";
my $nucseq = "";
my $match_flag = 1;
while( my $line = <BED> )
{
	chomp $line;
	my @field = split /\t/, $line;
	if ( $match_flag )
	{	
		$head = <NUC>;
		chomp $head;
		$nucseq = <NUC>;
		chomp $nucseq;
		$head =~ s/^>//;
	}
	if ( $head eq $field[3] )
	{
		$match_flag = 1;
		my $strand = "";
		my $mid = substr $nucseq, 1, 1;
		my $refbase = substr $field[4], 0, 1;
		if ( $mid ne $refbase )
		{
			die "Mismatched ref base!\n";
		}

		if ( $mid ne "C" && $mid ne "T" )
		{
			die "Nucseq is $nucseq for misformatted line: $line\n";
		}

		print "$field[0]\t$field[1]\t$field[2]\t$nucseq\t$field[4]\t$field[5]\n";
	}
	else
	{
		print STDERR "Error: sequence read names didn't match: $field[3] and $head\n";
		$match_flag = 0;
	}
}

