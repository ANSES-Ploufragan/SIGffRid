#!/usr/bin/perl -w

use strict;
use warnings;
use Env qw(SIGFFRID_DIR); #  ANSES_MEMORY_LIMIT);
use Getopt::Long;
use lib $SIGFFRID_DIR."libperl";
use SysUtilities;          # qw(tri_rapide);

$|++; # ensure display of print/errors in the order of appearing

# **********************************************************************

# **********************************************************************

=head1  NAME

Fsort_by_khi2score.pl

=head1 DESCRIPTION

Get motifs found by SIGffRid program by parsing grep results.


=head1  USAGE

=over

=item -grep_f <s>

file corresponding to a grep of result file on 'MOTIF' for 1 bacteria.
Redondant motifs in a same function directory are removed.

=back

example:
perl Fsort_by_khi2score.pl -grep_f grep_results_in_MOTIF.txt

=cut

# **********************************************************************

my $verbose                        = 0;
my $b_force                        = 0;
my $b_test                         = 0;

my $grep_f = undef;

# ($#ARGV == 0)or die "to use this programme\n$0 <file corresponding to a grep of result file on 'MOTIF' for 1 bacteria>\nredondant motifs in a same function directory are removed\n";

# **********************************************************************
# CHECK OPTIONS
# **********************************************************************
my $nbargmini = 1;
if(scalar(@ARGV) < $nbargmini){ 
  print "Bad number of arguments: ".scalar(@ARGV)." found, at least $nbargmini wanted\n";
  foreach(0..$#ARGV){
    print "$_: $ARGV[$_]\n";
  }
  die `perldoc $0`;
}
GetOptions(
    "grep_f=s"                 => \$grep_f
);
# **********************************************************************

open(F,'<', $grep_f)or die "$grep_f file cannot be read:$!, line ".__LINE__."\n";

$grep_f =~ /([^\/]+)\.txt$/;

my $input_file_name = $1;
my $path = $`;
# print "input_file_name $input_file_name, path s$path\n";
my @motifs = ();
my @R      = ();
my @LRT    = ();
my @ligne  = <F>;
my @dir    = ();
my @NBR    = ();

my $bool_dir = 0;

close F;

for my $ligne(@ligne){
#     $_ =~ /\/?(\d+\/)?\d+.*?MOTIF\s+(.+?)\s*, R: ([^\s]+).*?,\s+(?:LRT|Tr) is ([^\s]+)[^\)]+\)(.+)$/; 
    if($ligne =~ /\/?(\d+\/)?.*?MOTIF\s+([^\s]+)\s*, R: ([^\s]+).*?,\s+(?:LRT|Tr) is ([^\s]+)[^\)]+\)(.+)$/){
		$verbose and print "ligne $ligne\ndol: "; 
		
		if(defined $1){ push @dir, $1; 
				chop($dir[$#dir]); 
				$verbose and print "dir $1, ";
				$bool_dir = 1;
		}
		else{
		    push @dir, 0;
		}
		push @motifs, $2; 
		push @R ,     $3;
		push @LRT,    $4;
		push @NBR,    $5;
		if($verbose){
		    print "motif $2, R $3, LRT $4\n";
		    if($LRT[$#LRT] !~ /[0-9]\.+[0-9]/){
				die "LRT non num:".$LRT[$#LRT]."FFFF\n";
		    }
		}
		# }
		# else{ push @dir, 0; 
		#  push @motifs, $1; 
		#  push @R , $2;
		#  push @LRT, $3;
		#  print "ligne $_\ndol: motif $1, R $2, LRT $3\n";
		# }
		
		# die;
		# $_ = $&."\n";
		$_ = $5;
	# exit;
    }
    else{
		$verbose and print "ligne ne correspond pas à l'expression régulière:$ligne\n";
    }
}

# tri according to motif (NOT LRT, as IF 2 motifs have same LRT????)
# after LRT, we have to sort by motif to remove only same motifs for SAME dir!!
my $deb_tri = 0;
my $fin_tri = 0;
my $bool_diff = 0;
# sort only if we have dir of course :o)

my $f_sor_tr = $path.'res_sort_LRT_for_'.$input_file_name.'.txt';
my $f_sor_r  = $path.'res_sort_R_for_'.$input_file_name.'.txt';
open(F_SOR_TR,'>', $f_sor_tr)or die "$f_sor_tr cannot be opened:$!, line ".__LINE__."\n";
open(F_SOR_R,'>', $f_sor_r)or die "$f_sor_r cannot be opened:$!, line ".__LINE__."\n";

if($bool_dir) {
    for my $interv_tri(0..$#LRT-1){
	
		if( $LRT[$interv_tri] == $LRT[$interv_tri+1] ){
		    
		    if($bool_diff){ 
				$deb_tri = $interv_tri + 1;
				$bool_diff = 0;
		    }
		    $verbose and print "deb_tri, fin_tri: $deb_tri, $fin_tri, ligne ".__LINE__."\n";
		    ($deb_tri == $fin_tri)or &tri_rapide(  $deb_tri, $fin_tri, \@dir, \@NBR, \@motifs, \@LRT, \@R);
		    
		}
		else{
		    
		    $fin_tri = $interv_tri;
		    $verbose and print "deb_tri, fin_tri: $deb_tri, $fin_tri, ligne ".__LINE__."\n";
		    ($fin_tri != $#LRT)or($deb_tri == $fin_tri)or &tri_rapide(  $deb_tri, $fin_tri, \@dir, \@NBR, \@motifs, \@LRT, \@R);
		    $bool_diff = 1;
		}
    }

       
# after LRT, and motifs we have to sort by dir to remove only same motifs for SAME dir!!
    $deb_tri = $fin_tri = 0;

# sort only if we have dir of course :o)
    for my $interv_tri(0..$#dir-1){
		die "2EME BOUCLE\n";
		if( $motifs[$interv_tri] ne $motifs[$interv_tri+1] ){
		    
		    $deb_tri = $interv_tri + 1;
		    ($deb_tri == $fin_tri)or &TRI_R_LRT_SS_DIR($deb_tri, $fin_tri); # , \@dir, \@R, \@ligne, \@motifs, \@LRT);
		}
		else{
		    
		    $fin_tri = $interv_tri + 1;
		    ($fin_tri != $#dir)or($deb_tri == $fin_tri)or &TRI_R_LRT_SS_DIR($deb_tri, $fin_tri); # , \@dir, \@R, \@ligne, \@motifs, \@LRT); 
		}
    }
}
else{
    &TRI_R_LRT_SS_DIR(0, $#LRT); # , \@dir, \@R, \@ligne, \@motifs, \@LRT);
}


sub TRI_NOM_FINAL($$){
    my ($deb_subtri, $fin_subtri) = @_;
    my $prevRsub = -1;
    my $boolEqPrevR = 0;
    my $deb_subsubtri = 0;
    my $fin_subsubtri = 0;
    for my $j($deb_subtri..$fin_subtri){
		if($R[$j] == $prevRsub){
		    if($boolEqPrevR){
			$fin_subsubtri = $j;
			($j == $fin_subtri) and &tri_rapide_lettre($deb_subsubtri, $fin_subsubtri, \@motifs, \@dir, \@LRT, \@R, \@NBR);
		    }
		    else{
			$boolEqPrevR = 1;
			$deb_subsubtri = $j-1;
			$fin_subsubtri = $j;
		    }
		}
		else{
		    if($boolEqPrevR){
			&tri_rapide_lettre($deb_subsubtri, $fin_subsubtri, \@motifs, \@dir, \@LRT, \@R, \@NBR);
			$boolEqPrevR = 0;
		    }
		}
		$prevRsub = $R[$j];
		
    }
}


sub TRI_R_LRT_SS_DIR($$){
    my ($deb_tri, $fin_tri) = @_;

    $verbose and print "deb_tri, fin_tri: $deb_tri, $fin_tri ligne ".__LINE__."\n";
#    print "LRT @LRT\n";
# tri according to LRT, but IF 2 motifs have same LRT????... need other sorts
    &tri_rapide($deb_tri, $fin_tri, \@LRT, \@NBR, \@motifs, \@R, \@dir);
#     die "LRT: @LRT\n";
    # die "ON EST DEDANS, taille LRT $#LRT\n";
    my $prevLRT = -1;
    my $bool_to_sort = 0;
    my $deb_subtri = 0;
    my $fin_subtri = 0;
    for my $i ($deb_tri..$fin_tri){
		if($LRT[$i] == $prevLRT){
		    if($bool_to_sort){
			$fin_subtri = $i;
		    }
		    else{
			$bool_to_sort = 1;
			$deb_subtri = $i-1;
			$fin_subtri = $i;
		    }
	
		    if($i == $fin_tri){
			&tri_rapide($deb_subtri, $fin_subtri, \@R, \@NBR, \@motifs, \@dir, \@LRT);
			&TRI_NOM_FINAL($deb_subtri, $fin_subtri);
		    }
		}
		else{
		    if($bool_to_sort){
			&tri_rapide($deb_subtri, $fin_subtri, \@R, \@NBR, \@motifs, \@dir, \@LRT);
			&TRI_NOM_FINAL($deb_subtri, $fin_subtri);
			$bool_to_sort = 0;
		    }
		    $prevLRT = $LRT[$i];
		}
    }



    my $prev_motif = 'zzzzzzzzzz';
   #  my $prev_dir   = -1;
    my $cptsplice = 0;
    for my $i (reverse $deb_tri..$fin_tri){
	
		if($motifs[$i] eq $prev_motif){
		    # if(($dir[$i] !~ /(?:^|,)$prev_dir(?:$|,)/)and($#dir == -1)){
			# $dir[$i-1] .= ','.$dir[$i];
		    # }
		    splice @LRT   , $i, 1;
		    splice @R     , $i, 1;
		    splice @motifs, $i, 1;
		    splice @dir, $i, 1;
		    splice @NBR, $i, 1;
		    $cptsplice++;
		    # redo;
		}
		else{
		    # print F_SOR_TR 'DIR $dir[$i], MOTIF '.$motifs[$i].', R: $R[$i], LRT is LRT[$i]\n"; not reqd
		    print F_SOR_TR "DIR $dir[$i]; MOTIF $motifs[$i]; R: $R[$i]; LRT is $LRT[$i]; $NBR[$i]\n";
		    $prev_motif = $motifs[$i];
		    # ($#dir == -1)or $prev_dir   = $dir[$i];
		}
    }
    $fin_tri -= $cptsplice;
# print "$path".'res_sort_LRT_for_'.$input_file_name.".txt file created\n";
    
# tri suivant R
    $verbose and print "R $R[0], LRT $LRT[0], NBR $NBR[0], motifs $motifs[0], dir $dir[0]\n";
    &tri_rapide($deb_tri, $fin_tri, \@R, \@LRT, \@NBR, \@motifs, \@dir);
    
    my $prevR = -1;
    $bool_to_sort = 0;
    $deb_subtri = 0;
    $fin_subtri = 0;
    for my $i ($deb_tri..$fin_tri){
	$verbose and print "R $R[$i], LRT $LRT[$i], NBR $NBR[$i], motifs $motifs[$i], dir $dir[$i]\n";
	if(not defined $R[$i]){die;}
	if($R[$i] == $prevR){
	    if($bool_to_sort){
			$fin_subtri = $i;
	    }
	    else{
			$bool_to_sort = 1;
			$deb_subtri = $i-1;
			$fin_subtri = $i;
	    }

	    if($i == $fin_tri){
			&tri_rapide($deb_subtri, $fin_subtri, \@LRT, \@NBR, \@motifs, \@dir, \@R);
	    }
	}
	else{
	    if($bool_to_sort){
			&tri_rapide($deb_subtri, $fin_subtri, \@LRT, \@NBR ,\@motifs, \@dir, \@R);
			$bool_to_sort = 0;
	    }
	    $prevR = $R[$i];
	}
    }

    for my $i (reverse 0..$#R){
	
		print F_SOR_R "DIR $dir[$i]; MOTIF $motifs[$i]; R: $R[$i]; LRT is $LRT[$i]; $NBR[$i]\n";
		#  print F_SOR_R $ligne[$i];
    }
}

close F_SOR_TR;
close F_SOR_R;
print "$f_sor_tr file created\n";
print "$f_sor_r file created\n";
