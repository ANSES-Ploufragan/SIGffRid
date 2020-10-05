#!/usr/bin/perl -w

#genere les fichier de modele de Markov utilises par MEME pour le "bruit de fond" a partir de un ou plusieurs fichiers embl et/ou genbank et/ou (multi)fasta

# if we have embl or gbk file, DIRECT AND REVERSE sens are used
# if we have fasta file, only 
##Charu has made some changes

use strict;
use PerlIO::gzip;

sub sortir{
  die "This programme must be used as following:\n$0\n< Markov model order ( < 8 ) >\n(< EMBL or GenBank or fasta file where extract genome sequence or small individual sequences>)n\n< pathway to reach directory where you want results (/ for the current one) >\n[ direct sens only (1) (default for fasta files) or direct and reverse sens (0) (default for embl or gbk files) ]\n";
}
## Charu's change
# here we do nto have embl or gbk file but instead we areusing the fasta file for the entire sequence..so int his case the above die function states wrong. In the program we are choosing boolean as zero 

my $prog_tag = '[FmarkovBackgroundSeveralBacteria4.pl]';
my $ordre = 0;
my $bool_not_reverse = 1;

($#ARGV > 0) or &sortir;

( ($ARGV[0] =~ /^\d$/)and($ARGV[0] < 8) )or &sortir;
$ordre = $ARGV[0];

my $chemin_dest;
if($ARGV[$#ARGV] =~ /^[01]$/){ $chemin_dest = $ARGV[$#ARGV-1]; }
else{ $chemin_dest = $ARGV[$#ARGV]; }
if($chemin_dest !~ /\/$/){
    $chemin_dest .= '/';
}
elsif($chemin_dest eq '/'){
    $chemin_dest = '';
}


my %tab_ordre0 = ();
my %tab_ordre1 = ();
my %tab_ordre2 = ();
my %tab_ordre3 = ();
my %tab_ordre4 = ();
my %tab_ordre5 = ();
my %tab_ordre6 = ();
my %tab_ordre7 = ();
# my $total1 = 0;
# my $total2 = 0;
# my $total3 = 0;
# my $total4 = 0;
# my $total5 = 0;
# my $total6 = 0;
# my $total7 = 0;

my $partie_nom_fic_sor = '';
my $b_gz = 0;
my $der_b;
if($ARGV[$#ARGV] =~ /^[01]$/){ $der_b = $#ARGV-2; }
else{ $der_b = $#ARGV-1; }
 
for my $num_fic(1..$der_b){
  if(! -e $ARGV[$num_fic]){
    print "$prog_tag Unknown file $ARGV[$num_fic]\. Give an existing one!\n";
    &sortir;
  }
  $ARGV[$num_fic] =~ /([^\/]+)\.(?:embl|gbk|gb|txt|fasta|fa|GBK|fsa|fa\.gz)$/;
  if(defined $1)
  {
      $partie_nom_fic_sor .= "_$1";
	print "$prog_tag Output file name is: $partie_nom_fic_sor\n";
	
  }
  if($ARGV[$num_fic] =~ /\.gz/)
  {
      $b_gz = 1;
  }
  #die "dol1 $1\n";
}
print "$prog_tag der_b:$der_b, line ".__LINE__."\n";
print "$prog_tag partie_nom_fic_sor:$partie_nom_fic_sor, line ".__LINE__."\n";
if($ARGV[1] =~ /\.(?:embl|gbk|gb)$/){$bool_not_reverse = 0; }

if(-e $chemin_dest."markov_mod_order".$ordre.$partie_nom_fic_sor.".txt"){
  # print "markov_mod_order".$ordre.$partie_nom_fic_sor.".txt already exist, do you want overwrite(y for yes, n for no) ?\n";

  # my $reponse="y";
  # if($reponse =~ /^n$/ ){ 
  #   die "Program stopped not to overwrite existing file\n"
  # }
  # else{
    system("rm -f markov_mod_order$ordre$partie_nom_fic_sor.txt");
  # }
}
if($ARGV[$#ARGV] =~ /^([01])$/)
{
    $bool_not_reverse = $1;
    pop @ARGV;
}

for my $num_fic(1..$#ARGV-1){
  if(! -e $ARGV[$num_fic]){
    print "Unknown file $ARGV[$num_fic]\. Give an existing one!\n";
    &sortir;
  }
  &traitement_un_fic($ARGV[$num_fic]);
}

#getCosmidSequence [f] returns the complete cosmid nucleotid sequence from EMBL format file [f]. Resulting sequence is one big lowercase nucleotid string.
sub getCosmidSequence ($) {
  my $fileName = shift @_;
  my $result = ""; #Result buffer.
  my ($currentLine,$sequenceBlock);
  
  open(EMBLFILE, '<', $fileName) or die "$prog_tag [Error] Can't open EMBL GBK file \"$fileName\": $!, line ".__LINE__."\n"; #EMBL GBK format file.
  if($fileName =~ /\.embl$/){
      while($currentLine = <EMBLFILE>)
      {
	  if(($sequenceBlock) = ($currentLine =~ /^\s+([atgc\s]+)\d+$/)) #Current line is a DNA sequence line.
	  {
	      $sequenceBlock =~ s/\s+//g; #Eat up whitespace.
	      
	      #a chaque occurrence trouvee pour un \s, on la remplace par "pas de caract".On enleve les espaces. 	    
	      $result .= $sequenceBlock ;#Add the DNA block to final sequence.	    
	      #les lignes de la sequence entiere du cosmide sont mises bout a bout.
	  }
      }
  }
  elsif($fileName =~ /\.(?:gbk|gb)$/){
      while($currentLine = <EMBLFILE>)
      {
	  if(($sequenceBlock) = ($currentLine =~ /^\s+\d+([\satgc]+)$/)) #Current line is a DNA sequence line.
	  {
	      #a chaque occurrence trouvee pour un \s, on la remplace par "pas de caract".On enleve les espaces.
	      $sequenceBlock =~ s/\s+//g; #Eat up whitespace.
	      
	      
	      $result .= $sequenceBlock ;#Add the DNA block to final sequence.	      
	      #les lignes de la sequence entiere du cosmide sont mises bout a bout.
	  }
      }
      
  }
  else{
      die "You must give a filename with extensions .embl or .gbk !\n";
  }
  close EMBLFILE;
  return $result;
}

#l'argument correspond a la sequence a traiter
sub traitement_seq($){
  my $sequence = shift @_;
 	$sequence=~tr/ACGT/acgt/;
     	
  my $longueur = length($sequence);
  my $rev_seq;
  #if($bool_not_reverse){ print "DIRECT sens ONLY\n"; } 
  #else{                  print "DIRECT AND REVERSE sens\n"; } 
  if(! $bool_not_reverse)
  {   $rev_seq = reverse "$sequence";
      $rev_seq =~ tr/acgtACGT/tgcaTGCA/;
	 my $rev_longueur = length($rev_seq);
	
      #open(FS1,"> gbk0.txt");
      #print FS1 $sequence;
      #close FS1;
      #open(FS2,"> gbk1.txt");
      #print FS2 $rev_seq;
      #close FS2;
      #die;
  }

  #print "Treatment of sequence who begins with\n";
  #print substr($sequence,0,20).' and end with '.substr($sequence, $longueur-20, 20)."\n";
  #print "and the reverse one who begins with\n";
  #print substr($rev_seq,0,20).' and end with '.substr($rev_seq, $longueur-20, 20)."\n";
  
 BOUCLE:for(my $i=0; $i < $longueur; $i++){
    
    $tab_ordre0{ substr($sequence,$i,1) }++;

    ($bool_not_reverse)or $tab_ordre0{ substr($rev_seq,$i,1) }++;

    if($ordre > 0){
      ($i < $longueur-1)or next BOUCLE;
      $tab_ordre1{ substr($sequence,$i,2) }++;
      ($bool_not_reverse)or $tab_ordre1{ substr($rev_seq,$i,2) }++;
     # $total1 = $total1 + 2;
    }
    if($ordre > 1){
      ($i < $longueur-2)or next BOUCLE;
      $tab_ordre2{ substr($sequence,$i,3) }++;
      ($bool_not_reverse)or $tab_ordre2{ substr($rev_seq,$i,3) }++;
     # $total2 = $total2 + 2;
    }
    if($ordre > 2){
      ($i < $longueur-3)or next BOUCLE;
      $tab_ordre3{ substr($sequence,$i,4) }++;
      ($bool_not_reverse)or $tab_ordre3{ substr($rev_seq,$i,4) }++;
     # $total3 = $total3 + 2;
    }
    if($ordre > 3){
      ($i < $longueur-4)or next BOUCLE;
      $tab_ordre4{ substr($sequence,$i,5) }++;
      ($bool_not_reverse)or $tab_ordre4{ substr($rev_seq,$i,5) }++;
     # $total4 = $total4 + 2;
    }
    if($ordre > 4){
      ($i < $longueur-5)or next BOUCLE;
      $tab_ordre5{ substr($sequence,$i,6) }++;
      ($bool_not_reverse)or $tab_ordre5{ substr($rev_seq,$i,6) }++;
     # $total5 = $total5 + 2;
    }
    if($ordre > 5){
      ($i < $longueur-6)or next BOUCLE;
      $tab_ordre6{ substr($sequence,$i,7) }++;
      ($bool_not_reverse)or $tab_ordre6{ substr($rev_seq,$i,7) }++;
     # $total6 = $total6 + 2;
    }
    if($ordre > 6){
      ($i < $longueur-7)or next BOUCLE;
      $tab_ordre7{ substr($sequence,$i,8) }++;
      ($bool_not_reverse)or $tab_ordre7{ substr($rev_seq,$i,8) }++;
     # $total7 = $total7 + 2;
    }
  }
}

sub traitement_un_fic ($){
    
  my $fic = $_[0];

  print "Treatment of $fic file\n";
  #traitement s'il s'agit d'un fichier embl ou genbank
  if($fic =~ /\.(?:embl|gb|gbk)$/){
    #recuperation de la sequence ADN complete du fichier
    &traitement_seq(&getCosmidSequence($fic));
  }
  else{
      #traitement s'il s'agit d'un fichier fasta
      if($b_gz){  open(FIC,'<:gzip', $fic)or die "$prog_tag [Error] Impossible to open $fic file:$!, line ".__LINE__."\n"; }
      else     {  open(FIC,'<', $fic)or die "$prog_tag [Error] Impossible to open $fic file:$!, line ".__LINE__."\n";      }
    my $seq = '';
    while(my $ligne = <FIC>){
      if($ligne =~ /^>/ ){
	($seq eq '')or &traitement_seq($seq);
	$seq = '';
      }
      elsif($ligne =~ /^\w+/ ){
	chomp($seq .= $ligne);
        }
	
     }#end of while
   
   ($seq eq '')or &traitement_seq($seq);
    $seq = '';

  }
   
  print "$fic file treated\n";
}


open(FSOR,"> $chemin_dest"."markov_mod_order".$ordre.$partie_nom_fic_sor.".txt")or die "$prog_tag [Error] Impossible to create and open $chemin_dest"."markov_mod_order".$ordre.$partie_nom_fic_sor.".txt line ".__LINE__."\n";

print "Printing result file\n";

my @tab_n = ('a','c','g','t');

my $total = $tab_ordre0{ a }+$tab_ordre0{ c }+$tab_ordre0{ g }+$tab_ordre0{ t };
print FSOR "#ordre 0 $total\n";
printf FSOR "a %4.8f\n",($tab_ordre0{ a }/$total);
printf FSOR "c %4.8f\n",($tab_ordre0{ c }/$total);
printf FSOR "g %4.8f\n",($tab_ordre0{ g }/$total);
printf FSOR "t %4.8f\n",($tab_ordre0{ t }/$total);


my $totalMM = 0;  
my $prev_w = '';

if($ordre > 0){
    
    print FSOR "#ordre 1 $total\n";
    for my $tab(@tab_n){

	$prev_w = $tab;

	&TOTALMM(\$totalMM, \$prev_w, \%tab_ordre1, \@tab_n);

	for my $tab2(@tab_n){
	    if(exists $tab_ordre1{ $prev_w.$tab2 }){
		printf FSOR $prev_w.$tab2." %4.8f\n",($tab_ordre1{ $prev_w.$tab2 }/$totalMM);
		#Il peut arriver qu'une combinaison de lettre ne soit pas presente dans un genome, dans ce cas la case correspondante ne sera pas defini (mais c'est bien 0 qui est affcihe dans le fichier). Les lignes qui suivent servent simplement a indiquer a l'utilisateur ces "combinaisons nucleotidiques" absentes
	    }
	    else{
		printf FSOR  $prev_w.$tab2." 0.00000000\n";
		print "$prev_w$tab2  missing!! (without consequences on the result)\n";
	    }
	}
    }
}

if($ordre > 1){
    
    print FSOR "#ordre 2 $total\n";
    for my $tab(@tab_n){
	for my $tab2(@tab_n){

	    $prev_w = $tab.$tab2;

	    &TOTALMM(\$totalMM, \$prev_w, \%tab_ordre2, \@tab_n);

	    for my $tab3(@tab_n){
		if(exists $tab_ordre2{ $prev_w.$tab3 }){
		    printf FSOR $prev_w.$tab3." %4.8f\n",($tab_ordre2{ $prev_w.$tab3 }/$totalMM);
		}
		else{
		    printf FSOR $prev_w.$tab3." 0.00000000\n";
		    print "$prev_w$tab3  missing!! (without consequences on the result)\n";
		}
	    }
	}
    }
}


if($ordre > 2){
    print FSOR "#ordre 3 $total\n";
    for my $tab(@tab_n){
	for my $tab2(@tab_n){
	    for my $tab3(@tab_n){

		$prev_w = $tab.$tab2.$tab3;
		
		&TOTALMM(\$totalMM, \$prev_w, \%tab_ordre3, \@tab_n);

		for my $tab4(@tab_n){
		    if(exists $tab_ordre3{ $prev_w.$tab4 }){
			printf FSOR $prev_w.$tab4." %4.8f\n",($tab_ordre3{ $prev_w.$tab4 }/$totalMM);
		    }
		    else{
			printf FSOR  $prev_w.$tab4." 0.00000000\n";
			print "$prev_w$tab4  missing!! (without consequences on the result)\n";
		    }
		}
	    }
	}
    } 
}

if($ordre > 3){
    print FSOR "#ordre 4 $total\n";
    for my $tab(@tab_n){
	for my $tab2(@tab_n){
	    for my $tab3(@tab_n){
		for my $tab4(@tab_n){

		    $prev_w = $tab.$tab2.$tab3.$tab4;
		    
		    &TOTALMM(\$totalMM, \$prev_w, \%tab_ordre4, \@tab_n);

		    for my $tab5(@tab_n){
			if(exists $tab_ordre4{ $prev_w.$tab5 }){
			    printf FSOR $prev_w.$tab5." %4.8f\n",($tab_ordre4{ $prev_w.$tab5 }/$totalMM);
	      }
			else{
			    printf FSOR  $prev_w.$tab5." 0.00000000\n";
			    print "$prev_w$tab5  missing!! (without consequences on the result)\n";
			}
		    }
		}
	    }
	}
    }
}

if($ordre > 4){
    print FSOR "#ordre 5 $total\n";
    for my $tab(@tab_n){
	  for my $tab2(@tab_n){
	      for my $tab3(@tab_n){
		  for my $tab4(@tab_n){
		      for my $tab5(@tab_n){

			  $prev_w = $tab.$tab2.$tab3.$tab4.$tab5;
			  
			  &TOTALMM(\$totalMM, \$prev_w, \%tab_ordre5, \@tab_n);

			  for my $tab6(@tab_n){
			      if(exists $tab_ordre5{ $prev_w.$tab6 }){
				  printf FSOR $prev_w.$tab6." %4.8f\n",($tab_ordre5{ $prev_w.$tab6 }/$totalMM);
			      }
			      else{
				  printf FSOR  $prev_w.$tab6." 0.00000000\n";
				  print "$prev_w$tab6 missing!! (without consequences on the result)\n";
			      }
			      
			  }
		      }
		  }
	      }
	  }
      }
}

if($ordre > 5){
    print FSOR "#ordre 6 $total\n";
    for my $tab(@tab_n){
	for my $tab2(@tab_n){
	    for my $tab3(@tab_n){
		for my $tab4(@tab_n){
		    for my $tab5(@tab_n){
			for my $tab6(@tab_n){
			    
			    $prev_w = $tab.$tab2.$tab3.$tab4.$tab5;
			    
			    &TOTALMM(\$totalMM, \$prev_w, \%tab_ordre6, \@tab_n);
			    
			    for my $tab7(@tab_n){
				if(exists $tab_ordre6{ $prev_w.$tab7 }){
				    printf FSOR $prev_w.$tab7." %4.8f\n",($tab_ordre6{ $prev_w.$tab7 }/$totalMM);
				}
				else{
				    printf FSOR  $prev_w.$tab7." 0.00000000\n";
				    print "$prev_w$tab7 missing!! (without consequences on the result)\n";
				}
			    }
			}
		    }
		}
	    }
	}
    }
}

if($ordre > 6){
    print FSOR "#ordre 7 $total\n";
    for my $tab(@tab_n){
	for my $tab2(@tab_n){
	    for my $tab3(@tab_n){
		for my $tab4(@tab_n){
		    for my $tab5(@tab_n){
			for my $tab6(@tab_n){
			    for my $tab7(@tab_n){

				$prev_w = $tab.$tab2.$tab3.$tab4.$tab5;
				
				&TOTALMM(\$totalMM, \$prev_w, \%tab_ordre7, \@tab_n);

				for my $tab8(@tab_n){
				      if(exists $tab_ordre7{ $prev_w.$tab8 }){
					  printf FSOR $prev_w.$tab8." %4.8f\n",($tab_ordre7{ $prev_w.$tab8 }/$totalMM);
				      }
				      else{
					  printf FSOR  $prev_w.$tab8." 0.00000000\n";
					  print "$prev_w$tab8 missing!! (without consequences on the result)\n";
				      }
				  }
			    }
			}
		    }
		}
	    }
	}
    }
}

close FSOR;
print $chemin_dest."markov_mod_order".$ordre.$partie_nom_fic_sor.".txt was created\n";


# compute total of words for MM compute (to have something else than Bernouilli)
sub TOTALMM(\$\$\%\@)
{
    my ($r_totalMM, $r_prev_w, $r_tab_ordre, $r_tab_n) = @_;
    # we count subtotal of occurrences for previous letter(s)
    $$r_totalMM = 0; 
    for my $tab(@$r_tab_n){
	if(exists $r_tab_ordre->{ $$r_prev_w.$tab }){
	    $$r_totalMM += $r_tab_ordre->{ $$r_prev_w.$tab };
	}
    }
}
