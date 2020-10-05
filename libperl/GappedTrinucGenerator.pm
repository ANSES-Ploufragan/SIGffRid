#!/usr/bin/perl -w


# FILL WORD TABLE WITH EVERY COMBINATIONS OF LETTERS (JOKER ALLOWED)
# FILL all structures needed by rech_tri...49.pl
# SORT REGULAR EXPRESSIONS BEFORE TREATMENT ACCORDING TO THEIR LENGTHS 
use strict; 
# use diagnostics;
#die "@INC\n";
# $|++;

# package Gapped_trinuc_generator;
my $joker_l = 'n';
my $fixed_l = 'w';
my $not_verbose = 1;
my @alphabet = ('a','c','g','t');
my $size_alphabet = $#alphabet + 1;

# partie permettant (par appel) de preciser que le tri se fait sur la longueur des chaines
sub longueur{ length($a) <=> length($b) }

# r_regular_exp : ref to a 'regular expression' only with n (gap) or w (each of four letters)
sub proc_gapped_trinuc_generator(\@\@\@\$)
{
    my ($r_regular_exp, $r_trinuc_table, $r_trinuc_lengths, $r_str_trinuc_give_VAR) = @_;

    # hashage storing every trinuc lengths met (to store it in r_trinuc_lengths table at the end)

    @$r_regular_exp = sort longueur @$r_regular_exp;

    my %H_trinuc_lengths = ();

    foreach my $ind_re(0..$#$r_regular_exp)
    {
	# we store in anonymous table pointed by hash table the indice of the regular expression involved
	print "we put for length ".length( $r_regular_exp->[$ind_re] )." ind_re $ind_re in H for $r_regular_exp->[$ind_re]\n";
	push @{$H_trinuc_lengths{length( $r_regular_exp->[$ind_re] )}}, $ind_re;
	my @letters = split //, $r_regular_exp->[$ind_re];
	
	($not_verbose)or print "regexp  $r_regular_exp->[$ind_re], letters @letters\n";
	
	my @trinuc_by_regexp = ();

	for my $gp_l(@letters){  &letter_treatment($gp_l, \@trinuc_by_regexp);   }
	push @$r_trinuc_table, @trinuc_by_regexp;
	$$r_str_trinuc_give_VAR .= ($ind_re)x @trinuc_by_regexp;
    }
    for my $c(sort keys %H_trinuc_lengths )
    {
	# $c = length, $v = ref to table of regexp indices
	push @{$r_trinuc_lengths->[0]}, $c;
	# at the end, we only keep first and last indices of regexp related to this length (as they are sorted, we do not need indices of the others)
	
	push @{$r_trinuc_lengths->[1]}, [ $H_trinuc_lengths{$c}[ 0 ], $H_trinuc_lengths{$c}[ $#{ $H_trinuc_lengths{$c} } ] ];

	if(! $not_verbose)
	{
	    print 'length '.$r_trinuc_lengths->[0][ $#{$r_trinuc_lengths->[0]} ]."\n";
	    print "indice deb $r_trinuc_lengths->[1][ $#{$r_trinuc_lengths->[1]} ][0] var_trinuc $r_regular_exp->[ $r_trinuc_lengths->[1][ $#{$r_trinuc_lengths->[1]} ][0] ]\n";
	    print "indice fin $r_trinuc_lengths->[1][ $#{$r_trinuc_lengths->[1]} ][1] var_trinuc $r_regular_exp->[ $r_trinuc_lengths->[1][ $#{$r_trinuc_lengths->[1]} ][1] ]\n";
	}
	# print 'length ',$c," regexp indices @{$H_trinuc_lengths{$c}[1]} in H\n";
    }
    
    ($not_verbose)or &verif($r_trinuc_table,$r_trinuc_lengths);
}

sub letter_treatment(\$\@)
{
    my ($ll, $r_trinuc_by_regexp) = @_;

    # gap treatment
    # we add 'n' (joker letter) at the end of each word (we create a word 'n' if we do not have already word
    if($ll eq $joker_l)
    {
	if($#$r_trinuc_by_regexp != -1)
	{
	    foreach(@$r_trinuc_by_regexp){ $_ .= $joker_l; }
	}
	else
	{
	    @$r_trinuc_by_regexp = $joker_l;
	}
    }
    # every letter (4 possibilities)
    # each word is copied to have 4 corresponding words (size of the alphabet in fact)
    else
    {
	($ll eq $fixed_l)or print "caution: $ll treated as $fixed_l (not gap, choice 4 motifs, one foreach DNA letter\n";

	# ajout simple de l'alphabet si on a une lettre quelqconque en debut d'expression reguliere
	if($#$r_trinuc_by_regexp == -1)
	{
	    push @$r_trinuc_by_regexp, @alphabet;
	}
	else
	{
	    # augmente le tableau pour avoir les variations de lettres ajoutees dans la boucle qui suit
	    # ($not_verbose)or print "r_trinuc_by_regexp @$r_trinuc_by_regexp AV\n";
	    my @new_trinuc = (); 
	    for(@$r_trinuc_by_regexp) { push @new_trinuc, ($_)x $size_alphabet; }
	    @$r_trinuc_by_regexp = @new_trinuc;
	    # ($not_verbose)or print "r_trinuc_by_regexp @$r_trinuc_by_regexp AP\n";
	    
	    my $cpt = 0;

	    # ajoute 'a' a une ligne du tableau, 'c' a la suivante, 'g' a la suivante, 't' a la suivante et recommence jusqu'a la fin du tableau
	    foreach(@$r_trinuc_by_regexp){    $_ .= $alphabet[$cpt++ % $size_alphabet ];  }
	}
    }
}

sub verif(\@\@)
{
    my ($r_trinuc_table,$r_trinuc_lengths) = @_;

    my $cpt = 0;
    print "TRINUC\n";
    foreach(@$r_trinuc_table)
    {
	print $_,' ';
	(++$cpt % 4 == 0)and 	print "\n";
    }
    ($cpt % 4 == 0)or 	print "\n";
    print 'nb words '.($#{$_[0]}+1)."\n";

    $cpt = 0;
    print "LENGTHS\n";
    for my $i(0..$#{$r_trinuc_lengths->[0]})
    {
	print 'length ',$r_trinuc_lengths->[0][$i]," regexp indices @{$r_trinuc_lengths->[1][$i]}\n";
	$cpt++;
    }
    print "nb LENGTHS $cpt\n";
}


# test MAIN **********************************

if(! $not_verbose)
{
    ($#ARGV != -1) or die "$!\nUsage $0 < regexp [nw]+ >\n";
    my @re_tab = @ARGV;
    
    my @tab = ();
    my @l = ();
    
    &proc_gapped_trinuc_generator(\@re_tab,\@tab,\@l);
}

1;
