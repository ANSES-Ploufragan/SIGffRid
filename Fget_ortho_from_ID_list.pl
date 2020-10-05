#!/usr/bin/perl -w

# fonctionne
# retoune la liste d' orthologues correspondant
use strict;

# if( ($#ARGV != 3)or(!-e $ARGV[0])or(!-e $ARGV[1])or($ARGV[2] !~ /^[12]$/)or(! -d $ARGV[3]) ){
if( ($#ARGV != 2)or(!-e $ARGV[0])or(!-e $ARGV[1])or(! -d $ARGV[2]) ){
    die "To use: $0\n<file with ID list we have to found in the file of orthologues>\n<file of orthologues-MBGD_file>\n<result directory>";
#     die "To use: $0\n<file with ID list we have to found in the file of orthologues>\n<file of orthologues>\n<column number>\n<result directory>";
}
my ($file,$ortho_file,$result_dir) = @ARGV;
# my ($file,$ortho_file,$No_col,$result_dir) = @ARGV;

$file =~ /(.*?)([^\/]+)$/;

# for user, col 1 is the first (not col 0)
# $No_col++;
my $dir = '';
my $short_file = '';
if(defined $2)
{
    $dir = $1;
    $short_file = $2;
    $short_file =~ s/\.txt$//;
}

$ortho_file =~ /.*?([^\/]+)$/;

my $short_ortho_file = $1;
$short_ortho_file =~ s/\.txt$//;

$result_dir =~ s/^(.+?[^\/])$/$1\//;
my $name_fout = $result_dir."ortho_pairs_for_".$short_file."_in_".$short_ortho_file.".txt";
# die "result_dir $result_dir\n";
# my $regular_exp_searched = '';
my $regular_exp_searched = '([^\s]+)';


# if($No_col == 1){ $regular_exp_searched = '^([^\s]+)'; }
# else   { $regular_exp_searched = '^[^\s]+\s+([^\s]+)'; }

open(INF,"< $file")or die "Impossible to open $file: $!";

my @fout = (); # output written in file (we store it to sort and eliminate doubles)
#open(INF2,"< $ortho_file")or die "Impossible to open $ortho_file: $!";
#$fout[0] = <INF2>;
#my $prem_l = <INF2>;
#die @fout;
#close INF2;

my $prev_fout = -1;

open(ORTH,"< $ortho_file")or die "$!\n";
my @orth = (<ORTH>);

my $mes = shift @orth;
chomp $mes;
print "$mes line removed (first line)\n";
close ORTH;

while(my $line = <INF>)
{
    if($line =~ /^>/)
    {
	print "ligne intitule $line\n";
	next; 
    }
    elsif($line =~ /$regular_exp_searched/)
    {
	my $expression = $1.'(?:\s|\,|$)';
	# open(GREP,"grep '$1' $ortho_file | ");

	# my @G = (<GREP>);
	# grep $1, @orth;
	push @fout, grep /$expression/, @orth;

	# push @fout, "@G";
	# if($#G == -1)
	if($#fout == $prev_fout)
	{
	    print "on ne trouve PAS $1\n";
	}
	else
	{
	    print "on trouve $1 dans $fout[$#fout]";
	}
	# close GREP;
	$prev_fout = $#fout;
    }
    else
    {
	print "cas $line not treated\n$line does not look like $regular_exp_searched\n";
    }
}

print "$name_fout created\n";
close INF;

@fout = sort @fout;
open(FOUT,"> $name_fout")or die "Impossible to create $name_fout: $!";
#print FOUT $prem_l;
print FOUT "$mes\n";
for my $i(1..$#fout)
{
    ($fout[$i] ne $fout[$i-1])and print FOUT $fout[$i];
}
close FOUT;
