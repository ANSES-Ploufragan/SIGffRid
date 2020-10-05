#!/usr/local/bin/perl -w

use strict;

# ******************************************************************************
# ***************** DATE
# ******************************************************************************
sub date ($)
{
    my ($separateur) = @_;
    (defined $separateur)or $separateur = '';

#   0      1     2        3      4       5       6       7       8
# my($sec, $min, $heure, $mjour, $mois, $annee, $sjour, $ajour, $est_dst) = localtime;
# die localtime( time());
    my ($ce_jour_lettre, $ce_jour_num, $ce_mois, $ce_mois_num, $cette_annee, $date);
    $ce_jour_lettre = ( qw(sunday monday tuesday wednesday thursday friday saturday))[(localtime)[6]];
    $ce_jour_num    = sprintf "%02u",(localtime)[3];
    $ce_mois_num    = sprintf "%02u",(localtime)[4]+1;
    # $ce_mois        = ( qw(januar februar march april may june july august september october november december))[(localtime)[4]];
    $cette_annee    = 1900 + (localtime)[5];
    $date = $ce_jour_num.$separateur.$ce_mois_num.$separateur.substr($cette_annee,2,2);
    return $date;
    #"DATE: $ce_jour_lettre, the $ce_jour_num of $ce_mois $cette_annee\n";
}    
# ******************************************************************************

# ******************************************************************************
# TRI RAPIDE CHIFFRE                ********************************************
# ******************************************************************************

{
    no strict;
    
    sub tri_rapide($$\@@) # ($$\@\@\@\@)
    {
	local($min, $max, $r_arrayTri, @arrayVal) = @_;
	local($pos);
	# print "min $min, max $max; ";
	# print "r_arrayTri $r_arrayTri->[0]\n";
	# die;
	# print "arrayval $#arrayVal tri_rapide devrait etre ".($#_-2)."\n";
	# Est-ce que la position actuelle du tri est plus petite que la taille du tableau à trier ?
	# Oui alors on commence/continue le tri
	if ( $min < $max )
	{
	    # we get the position where last sort stop
	    $pos = &partitionner($min, $max, $r_arrayTri, @arrayVal);
	    
	    # we sort the left section of the remaining section in the table
	    &tri_rapide($min, $pos, $r_arrayTri, @arrayVal);
	    # we sort the right section of the remaining section in the table
	    &tri_rapide($pos+1, $max, $r_arrayTri, @arrayVal);
	}
	# Non, on peut donc s'arreter, le tri est terminé !
	else
	{
	    return( 1 );
	}
    }

sub partitionner($$\@\@)
{ 
    local($min, $max, $r_arrayTri, @arrayVal) = @_;
    local($i, $j, $centre, $tmp);
    
    # we record as temp. value for the center this of the left value
    $centre = $r_arrayTri->[$min];
    $i = $min - 1;
    $j = $max + 1;
    
    
    while ( 1 )
    {
	$j--;
	$i++;
	# print "arrayTri de $j $arrayTri[$j]\n";

	until (( $r_arrayTri->[$j] <= $centre )or( $j == 0 )) { $j--; }
	until (( $r_arrayTri->[$i] >= $centre )or( $i == $#$r_arrayTri )) { $i++; }
	if(! defined $r_arrayTri->[$j]){ die ; }
	if(! defined $r_arrayTri->[$i]){ die ; }
	# die;
	if ( $i < $j )
	{
	    $tmp = $r_arrayTri->[$j];
	    $r_arrayTri->[$j] = $r_arrayTri->[$i];
	    $r_arrayTri->[$i] = $tmp;
	    # print "arrayval $#arrayVal partitionner devrait etre ".($#_-2)."\n";
	    # print "arrayval $#arrayVal concerne les données:\n";
	    # die;
	    for(my $i_param=0; $i_param <= $#arrayVal; $i_param++)
	    {
		# print "arrayval $arrayVal[$i_param][$j]\n";
		$tmp = $arrayVal[$i_param][$j];
		$arrayVal[$i_param][$j] = $arrayVal[$i_param][$i];
		$arrayVal[$i_param][$i] = $tmp;
	    }
	    # die;
	}
	else
	{
	    return( $j );
	}
    }
}


# **********************************END TRI RAPIDE CHIFFRE**********************
# ******************************************************************************

# ******************************************************************************
# **********************************DEB TRI RAPIDE LETTRE***********************

sub tri_rapide_lettre($$\@\@)
{
    local($min, $max, $r_arrayTri, @arrayVal) = @_;
    local($pos);
    # print "min $min, max $max; ";
    # print "r_arrayTri $r_arrayTri->[0]\n";
    # die;
    # print "arrayval $#arrayVal tri_rapide devrait etre ".($#_-2)."\n";
    # Est-ce que la position actuelle du tri est plus petite que la taille du tableau à trier ?
    # Oui alors on commence/continue le tri
    if ( $min < $max )
    {
	# we get the position where last sort stop
	$pos = &partitionner_lettre($min, $max, $r_arrayTri, @arrayVal);
	
	# we sort the left section of the remaining section in the table
	&tri_rapide_lettre($min, $pos, $r_arrayTri, @arrayVal);
	# we sort the right section of the remaining section in the table
	&tri_rapide_lettre($pos+1, $max, $r_arrayTri, @arrayVal);
    }
    # Non, on peut donc s'arreter, le tri est terminé !
    else
    {
	return( 1 );
    }
}

sub partitionner_lettre($$\@\@)
{ 
    local($min, $max, $r_arrayTri, @arrayVal) = @_;
    local($i, $j, $centre, $tmp);
    
    # we record as temp. value for the center this of the left value
    $centre = $r_arrayTri->[$min];
    $i = $min - 1;
    $j = $max + 1;
    
    
    while ( 1 )
    {
	$j--;
	$i++;
	# print "arrayTri de $j $arrayTri[$j]\n";

	until (( $r_arrayTri->[$j] le $centre )or( $j == 0 )) { $j--; }
	until (( $r_arrayTri->[$i] ge $centre )or( $i == $#$r_arrayTri )) { $i++; }
	if(! defined $r_arrayTri->[$j]){ die ; }
	if(! defined $r_arrayTri->[$i]){ die ; }
	# die;
	if ( $i < $j )
	{
	    $tmp = $r_arrayTri->[$j];
	    $r_arrayTri->[$j] = $r_arrayTri->[$i];
	    $r_arrayTri->[$i] = $tmp;
	    # print "arrayval $#arrayVal partitionner devrait etre ".($#_-2)."\n";
	    # print "arrayval $#arrayVal concerne les données:\n";
	    # die;
	    for(my $i_param=0; $i_param <= $#arrayVal; $i_param++)
	    {
		# print "arrayval $arrayVal[$i_param][$j]\n";
		$tmp = $arrayVal[$i_param][$j];
		$arrayVal[$i_param][$j] = $arrayVal[$i_param][$i];
		$arrayVal[$i_param][$i] = $tmp;
	    }
	    # die;
	}
	else
	{
	    return( $j );
	}
    }
}

# end no strict

}


# ******************************************************************************
# clone an array or H (multidimension possible) and return a ref to the duplicated struct
# ******************************************************************************

sub clone(\$){
    my ($r_original_array) = @_;

    my $ref = ref($r_original_array);

    if(not $ref){
	# print "oui, ref $ref r_original_array $r_original_array\n";
	return $r_original_array;
    }
    else{
	if($ref eq "ARRAY"){
	    my @anon_array = ();
	    for(@$r_original_array){ 
	      push @anon_array, &clone($_);
	    }
	    # print "anon_array: @anon_array\n";
	    return \@anon_array;
	}
	elsif($ref eq "HASH"){
	  my %anon_H = ();
	  for(keys %$r_original_array){
	    # print "c $_, $r_original_array->{$_}\n";
	    $anon_H{$_} = &clone($r_original_array->{$_});
	  }
	  #  print "anon_H: \n";
	  # 	  for( my ($c,$v) = each %anon_H){
	  # 	    print "c: $c, v: $v\n";
	  # 	  }
	  # 	  print "\n";
	  return \%anon_H;
	}
	elsif($ref eq "REF"){
	  my $new_ref = '';
	  # print "appel clone $$r_original_array\n";
	  $new_ref = \&clone($$r_original_array);
	  return \$new_ref;
	}
	else{
	    die "kind of data not treated ($_) in clone sub line ".__LINE__."\n";
	}
    }
    return $r_original_array;
}

# test

# array

# my @toto;
# push @{ $toto[0]}, (1,2);
# push @{ $toto[1]}, (6,7);
# my @toto2 = ();
# my $rtoto2;
# $rtoto2 = &clone(\@toto);
# @toto2 = @$rtoto2;
# $toto[0][0] = 3;
# my $totodisplay = '';
# my $toto2display = '';
# for my $i(0..1){
#   for my $j(0..1){
#     $totodisplay .= "toto $toto[$i][$j],";
#     $toto2display .= "toto $toto2[$i][$j],";
#   }
#   $totodisplay .= "\n";
#   $toto2display .= "\n";
# }
# print "toto @toto:\n";
# print $totodisplay;
# print "toto2 @toto2:\n";
# print $toto2display;
# exit;



# hash

# my %toto;
# %{ $toto{'a'}} = (
# 		c => 12,
# 		g => 4
# 		);
# %{ $toto{'b'} } = (
# 		 u => 6,
# 		 t => 7
# 		 );
# my %toto2 = ();
# print "toto %toto\n";
# my $rtoto2 = &clone(\%toto);
# %toto2 = %$rtoto2;
# $toto{'a'}{'c'} = 3;
# my $totodisplay = '';
# my $toto2display = '';
# for my $i('a','b'){
#     for my $j(sort keys %{ $toto{$i} } ){
# 	$totodisplay .= "toto $toto{$i}{$j},";
# 	$toto2display .= "toto $toto2{$i}{$j},";
#     }
#     $totodisplay .= "\n";
#     $toto2display .= "\n";
# }
# print "toto %toto:\n";
# print $totodisplay;
# print "toto2 %toto2:\n";
# print $toto2display;
# exit;

# ******************************************************************************



1;
