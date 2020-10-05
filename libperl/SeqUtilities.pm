#!/usr/local/bin/perl -w

use strict;

my $debug = 0;

# ************************************************************************************************    
#getCosmidSequence [f] returns the complete cosmid nucleotid sequence from EMBL or GBK format file [f]. Resulting sequence is one big lowercase nucleotid string.
sub get_cosmid_seq ($$) {
  my ($fileName, $bool_extended_alphabet) = @_;
  my $result = ""; #Result buffer.
  my ($currentLine,$sequenceBlock, $lg);

  my $normal_alphabet = 'acgt';
  my $extended_alphabet = 'abcgkmnrstwy';
  my $alphabet;
  if($bool_extended_alphabet)
  {
      $alphabet = $extended_alphabet;
  }
  else
  {
      $alphabet = $normal_alphabet;
  }

  open EMBLFILE, "< $fileName" or die "Can't open EMBL GBK file \"$fileName\": $!\n"; #EMBL GBK format file.
  #my $bool_ori = 0;
  # my %lettre = ();
  if($fileName =~ /\.(gbk|gb)$/)
  {
      
	while($currentLine = <EMBLFILE>)
	{
	    if($currentLine =~ /^ORIGIN/){ last;}
	 # if($currentLine =~ /^\s+(\d+)?([\satgc]+)(\d+)?$/) #Current line is a DNA sequence line.
	} 
	if(eof(EMBLFILE))
	{
	    die("ORIGIN dans fichier non trouvé pour début de la sequence\n");
	}
	while($currentLine = <EMBLFILE>)
	{
	    if($currentLine =~ /^\s+(\d+)([\s$alphabet]+)$/) #Current line is a DNA sequence line.
	    {
		$sequenceBlock = $2;
		$lg  = $1;
		#a chaque occurrence trouvee pour un \s, on la remplace par "pas de caract".On enleve les espaces. 	
		$sequenceBlock =~ s/\s+//g; #Eat up whitespace.
		
		#les lignes de la sequence entiere du cosmide sont mises bout a bout.     
		$result .= $sequenceBlock ;#Add the DNA block to final sequence.
	    }
	}
	if( ($lg + length($sequenceBlock)-1) > length($result))
	{
	    die (" Longueur de la sequence extraite(".length($result).") inferieur de la longueur attendue(".($lg + length($sequenceBlock)-1).")\n Utilisez l'alphabet etendu\n");
	}
	elsif( ($lg + length($sequenceBlock)-1) < length($result))
	{
	    die (" Longueur de la sequence extraite(".length($result).") superieur de la longueur attendue(".($lg + length($sequenceBlock)-1).")\n");
	}
	  
	  
	   # elsif(($bool_ori == 1 )and($currentLine =~ /^\s+(\d+)([\w\s]+)/)){
	#	$currentLine =~ s/[acgt\d ]//g;
	#	for my $l(split //, $currentLine){
	#	    if(not exists $lettre{$l}){
	#		$lettre{$l} = 1;
	#		print "lettre $l\n";
	#	    }
	#	}
	#    }
	

  }
  else{
      die "You must give a filename with extensions .embl or .gbk ! (you gave $fileName)\n";
  }
  close EMBLFILE;
  # print $result;
  # die "\n";
  return $result;
}
# ************************************************************************************************    

# ************************************************************************************************    
# "complementSequence [seq]" returns the complement of nucleotid sequence [seq]. Nucleotids may be uppercase or lowercase.
# DOES NOT CHECK SEQUENCE FOR ILLEGAL NUCLEOTIDS.
# CAUTION, THERE IS A MAX SIZE for string passed as parameter!!!!!!!!!!!!!!!!!!!!!
sub complement_seq ($) {
    my $dnaSequence = shift @_;
    my $reverseSeq;
  
    $reverseSeq = reverse "$dnaSequence"; #Sequence is reversed
    $reverseSeq =~ tr/atgcATGC/tacgTACG/; #Individual nucleotids are complemented
    return $reverseSeq;
}
# sub complementSequence(\$) 
#		       {
#			   my $r_dnaSequence = shift @_;
#			   $$r_dnaSequence = reverse "$$r_dnaSequence"; # Sequence is reversed
#			   $$r_dnaSequence =~ tr/atgcATGC/tacgTACG/; # Individual nucleotids are complemented
			   # return $reverseSeq;
#		       }

# ************************************************************************************************    
#permet l'affichage de sequences par series de 50 nucleotides (retour à la ligne tous les 50 nucleotides). Cela evite de devoir passer les fichiers par READSEQ pour les transformer en fichiers au format FASTA
# param 1: booleen pour savoir si l'on veut un affichage au format fasta ou non
# param 2: reference pointant sur la sequence traitee
my $largeur = 50;

sub print50 ($$)
{
    my ($bool_fasta, $seq_print) = @_;

    ($bool_fasta)or return lc($seq_print)."\n";
    
    (! defined $seq_print)and die "[print50] arg seq_print non def pour print50, line ".__LINE__."\n";
    
    my $res;
    my $longueur_seq = length($seq_print);

    if($longueur_seq >= $largeur)
    {
	for my $affichage(0..($longueur_seq / $largeur-1))
	{
	    $res .= substr($seq_print,($largeur*$affichage),$largeur)."\n";
	}
    }
    $res .= substr($seq_print,$longueur_seq-($longueur_seq % $largeur),$longueur_seq % $largeur);
    ($longueur_seq % $largeur == 0)or $res .= "\n";
     return lc($res);
}
# ************************************************************************************************    

# ******************************************************************************
# ******* sub to create an index file (in another file) for fast search in FASTA file
sub create_index_file($$)
{
    my ($fasta_f, $index_f) = @_;
# current pos in file (byte)

    open(FASTA_F, "< $fasta_f")or die "Can not open $fasta_f in $0, line ".__LINE__."\n";
    open(INDEX_F, "> $index_f")or die "Can not open $index_f in $0, line ".__LINE__."\n";
    my $prev_tell = tell FASTA_F;

    # my $file = 'toto';
    # (-e $file.'.txt')and $file .= 1;
    # $file .= '.txt';
    # open(TOTO, "> $file")or die "Can not open $file line ".__LINE__."\n";

    # my $cptdebug = 1;
    while(<FASTA_F>){
	/^>.*?\_(.+?)\s/ and do { # print "previous byte for $1: $prev_tell\n"; 
	    print INDEX_F $prev_tell,' ';
	    # print TOTO "line $cptdebug, octet $prev_tell, ID $_\n";
	    $prev_tell = tell FASTA_F;
	    # print "new_tell $prev_tell\n";
	    # $cptdebug++;
	};
	$prev_tell = tell FASTA_F;
	# print "tell line without > $prev_tell\n";

    }
    # we do not need to record $prev_tell because last line is a real annotation (so, already done)
    # the last record corresponds to the end byte
    print INDEX_F $prev_tell;
    # print "line $cptdebug, octet $prev_tell\n";


    close FASTA_F;
    close INDEX_F;

    # tell handle 
    
    # seek handle, decallage, depart (0 deb);
    # read handle, scalaire, longueur (byte) [, decallage];

    # while (read DEPUIS, $tampon, 16384) {
	# print VERS $tampon;
    # }
    # 0 => fin de fichier
    # undef => erreur
}

# ******************************************************************************
# ******************************************************************************

sub load_index_f($\@)
{
    my ($index_f, $r_index_table) = @_;

    open(INDEX_F, "< $index_f")or die "Can not open $index_f in $0, line ".__LINE__."\n";
    # push @$r_index_table, 0;
    while(<INDEX_F>){
	push @$r_index_table, split ' ', $_;
    }
    # shift @$r_index_table;
    close INDEX_F;
}


if($debug){
    my $fasta_f = "../../donnees_brutes/010205/annot_promot/annot_SCO.txt";
    my $index_f = "test.txt";
    my @index_table = ();
    
    &create_index_file($fasta_f, $index_f);
    &load_index_f($index_f, \@index_table);
    &verif_annot_byte_position(\@index_table);
}
# for my $l(0..$#index_table){
    #   print "line $l, octet $index_table[$l]\n";
# }
# my $tampon = '';
# open(FANNOT, $fasta_f)or die __LINE__;

sub verif_annot_byte_position(\@){
    my ($r_index_table) = @_;

    print "verif annot\n\n";
    for my $l (0..$#$r_index_table){
	print "l vaut $l, octets table $r_index_table->[$l]\n";
	#  seek(FANNOT, $index_table[$l], 0);
	
	#  read FANNOT, $tampon, $index_table[$l+1] - $index_table[$l];
	#  print "line $l, $tampon\n";
    }
    print "END verif annot\n";
}

1;

